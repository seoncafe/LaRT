module h2_mod
!===============================================================================
!  Molecular-hydrogen (H2) effect on Lyman-alpha radiative transfer.
!
!  Some H2 Lyman (B-X) and Werner (C-X) band lines lie within a few hundred km/s
!  of H Ly-alpha.  A Ly-alpha photon at the right local-frame frequency can be
!  absorbed into an excited electronic state of H2 which then fluoresces to other
!  wavelengths (destruction) or decays back to the pumping level (resonance
!  scattering near Ly-alpha).  This module supplies the H2 line opacity and the
!  branching data; the scattering kinematics live in scattering_car.f90.
!
!  Phase 1: Neufeld (1990, ApJ 350, 216) two-line treatment -- the two B-X lines
!  closest to Ly-alpha, R(6) (+14.1 km/s) and P(5) (+99.2 km/s), both v=1 <- v=2.
!  Uniform molecular fraction (par%f_H2) and a single H2 temperature
!  (par%h2_temperature) that sets BOTH the LTE level populations and the H2
!  Doppler width.  The opacity is a constant multiple of the local H I opacity
!  grid%rhokap (Neufeld's dimensionless strength s_i), so no cell-by-cell array.
!
!  The design mirrors the existing H+D (line_type = 7) second-species path: a
!  co-moving absorber with its own Doppler width and Voigt-a, added on top of the
!  H I Voigt profile.
!===============================================================================
  use define
  use voigt_mod, only : voigt
  implicit none
  private

  !--- one H2 absorption line near Ly-alpha
  type :: h2_line_t
     real(kind=wp) :: dv_kms      ! velocity offset from Ly-alpha (redward = +)
     real(kind=wp) :: dnu_Hz      ! nu_line - nu_Lya = -(dv/c) nu_Lya  [Hz]
     integer       :: vl, Jl, Ju  ! lower vib., lower/upper rot. quantum numbers
     real(kind=wp) :: lambda_A    ! rest wavelength [Angstrom]
     real(kind=wp) :: A_ul        ! Einstein A of this transition [s^-1]
     real(kind=wp) :: A_tot_up    ! total radiative decay rate of the upper level [s^-1]
     real(kind=wp) :: f_osc       ! absorption oscillator strength (lower -> upper)
     real(kind=wp) :: pop         ! LTE fractional population of the lower level
     real(kind=wp) :: strength    ! s_i = f_H2 * pop * f_osc/f_osc_Lya (rel. to H I)
     real(kind=wp) :: a_damp      ! Voigt a-parameter in H2 Doppler units
     real(kind=wp) :: p_scat      ! return-to-Ly-alpha probability = A_ul / A_tot_up
  end type h2_line_t

  integer,                       public :: n_h2_lines = 0
  type(h2_line_t), allocatable,  public :: h2l(:)
  real(kind=wp),                 public :: h2_Dfreq_Hz = 0.0_wp   ! H2 Doppler width [Hz]
  logical,                       public :: h2_on = .false.
  !--- pumping-weight accumulator for each line (kept here, not in par, so the params
  !--- namelist read is not broken by an allocatable derived-type field).
  real(kind=wp), allocatable,    public :: W_H2pump(:)

  real(kind=wp), parameter :: f_osc_Lya   = 0.4162_wp
  real(kind=wp), parameter :: hc_over_k   = 1.4387769_wp          ! [cm K] for E in cm^-1
  real(kind=wp), parameter :: f_osc_const = 1.4992e-16_wp         ! f = const*(gu/gl)*lam_A^2*A_ul

  public :: h2_init, h2_kappa, h2_select_line, h2_kappa_D, h2_select_line_D

contains

  !=============================================================================
  function h2_kappa_D(xfreq, Dfreq) result(kap)
  !--- Core H2 opacity multiplier (of the H I opacity) at frequency xfreq in a
  !--- cell whose H I Doppler width is Dfreq [Hz].  Grid-agnostic so both the
  !--- Cartesian (grid%Dfreq(i,j,k)) and AMR (amr%Dfreq(il)) paths can call it.
  real(kind=wp), intent(in) :: xfreq, Dfreq
  real(kind=wp) :: kap, ratio, dx, x_h2
  integer :: il
  kap = 0.0_wp
  if (.not. h2_on) return
  if (par%h2_hi_width) then
     ratio = 1.0_wp
  else
     ratio = Dfreq / h2_Dfreq_Hz
  end if
  do il = 1, n_h2_lines
     dx   = h2l(il)%dnu_Hz / Dfreq
     x_h2 = (xfreq - dx) * ratio
     kap  = kap + h2l(il)%strength * ratio * voigt(x_h2, h2l(il)%a_damp)
  end do
  end function h2_kappa_D

  !=============================================================================
  function h2_select_line_D(xfreq, Dfreq) result(il_sel)
  !--- Core line selector (opacity-weighted) for a cell of H I Doppler width Dfreq.
  use random, only : rand_number
  real(kind=wp), intent(in) :: xfreq, Dfreq
  integer :: il_sel, il
  real(kind=wp) :: ratio, dx, x_h2, w(64), wtot, r, acc
  if (par%h2_hi_width) then
     ratio = 1.0_wp
  else
     ratio = Dfreq / h2_Dfreq_Hz
  end if
  wtot = 0.0_wp
  do il = 1, n_h2_lines
     dx    = h2l(il)%dnu_Hz / Dfreq
     x_h2  = (xfreq - dx) * ratio
     w(il) = h2l(il)%strength * ratio * voigt(x_h2, h2l(il)%a_damp)
     wtot  = wtot + w(il)
  end do
  il_sel = 1
  if (wtot <= 0.0_wp) return
  r   = rand_number() * wtot
  acc = 0.0_wp
  do il = 1, n_h2_lines
     acc = acc + w(il)
     if (r <= acc) then
        il_sel = il
        return
     end if
  end do
  il_sel = n_h2_lines
  end function h2_select_line_D

  !=============================================================================
  subroutine h2_init(data_dir)
  !--- Build the Phase-1 line table, read X-state energies for the LTE partition
  !--- function, and precompute strengths / damping / branching at
  !--- par%h2_temperature and par%f_H2.  Called once from setup (all ranks).
  character(len=*), intent(in) :: data_dir
  real(kind=wp) :: nu_Lya_Hz, vth1_H2, vth_H2_kms, T, bturb2
  real(kind=wp) :: Zpart, Ecm, gns
  real(kind=wp), allocatable :: Ev(:), Jv(:), Eng(:)
  integer :: il, i, nlev
  character(len=256) :: fname

  h2_on = .false.
  if (trim(par%h2_model) == 'none') return

  !--- Ly-alpha line-center frequency [Hz] from the H rest wavelength [um].
  nu_Lya_Hz = (speedc*1.0e5_wp) / (line%wavelength0 * um2m * 1.0e2_wp)  ! c[cm/s]/lambda[cm]

  !--- H2 thermal velocity at 1 K [km/s]: scale the H value by sqrt(m_H/m_H2).
  !--- Turbulence (par%bturb) applies equally to all species; the default
  !--- sentinel par%bturb <= 0 means "no turbulence".
  vth1_H2    = line%vtherm1 * sqrt(line%mass_amu / (2.0_wp*line%mass_amu))
  T          = par%h2_temperature
  bturb2     = 0.0_wp
  if (par%bturb > 0.0_wp) bturb2 = par%bturb**2
  vth_H2_kms = sqrt( (vth1_H2*sqrt(T))**2 + bturb2 )
  h2_Dfreq_Hz = nu_Lya_Hz * vth_H2_kms / speedc                        ! [Hz]

  !--- Phase-1 line table: the two Neufeld B-X lines closest to Ly-alpha
  !--- (B v=1 <- X v=2). Data: dv from CLOUDY energies; A from Abgrall+00;
  !--- A_tot_up = total radiative decay of the upper B(v=1,J) level.
  n_h2_lines = 2
  if (allocated(h2l)) deallocate(h2l)
  allocate(h2l(n_h2_lines))
  h2l(1)%dv_kms = 14.140_wp                    ! R(6): B(1,7) <- X(2,6)
  h2l(1)%vl = 2; h2l(1)%Jl = 6; h2l(1)%Ju = 7
  h2l(1)%lambda_A = 1215.72534_wp; h2l(1)%A_ul = 1.36e8_wp; h2l(1)%A_tot_up = 1.6825e9_wp
  h2l(2)%dv_kms = 99.229_wp                    ! P(5): B(1,4) <- X(2,5)
  h2l(2)%vl = 2; h2l(2)%Jl = 5; h2l(2)%Ju = 4
  h2l(2)%lambda_A = 1216.07038_wp; h2l(2)%A_ul = 1.59e8_wp; h2l(2)%A_tot_up = 1.7199e9_wp

  !--- Read X-state rovibrational energies for the LTE partition function.
  fname = trim(data_dir)//'/energy_X.dat'
  call read_energy_X(trim(fname), Ev, Jv, Eng, nlev)

  !--- LTE partition function Z = sum g_ns(J) (2J+1) exp(-hc E / kT), E in cm^-1.
  !--- Nuclear-spin weight: para-H2 (J even) = 1, ortho-H2 (J odd) = 3.
  Zpart = 0.0_wp
  do i = 1, nlev
     gns   = 1.0_wp
     if (mod(nint(Jv(i)), 2) /= 0) gns = 3.0_wp
     Zpart = Zpart + gns * (2.0_wp*Jv(i)+1.0_wp) * exp(-hc_over_k*Eng(i)/T)
  end do

  !--- Per-line derived quantities.
  do il = 1, n_h2_lines
     h2l(il)%dnu_Hz  = -(h2l(il)%dv_kms/speedc) * nu_Lya_Hz
     !--- absorption oscillator strength from the tabulated A_ul
     h2l(il)%f_osc   = f_osc_const * (2.0_wp*h2l(il)%Ju+1.0_wp)/(2.0_wp*h2l(il)%Jl+1.0_wp) &
                       * h2l(il)%lambda_A**2 * h2l(il)%A_ul
     !--- LTE population of the lower level (vl, Jl)
     Ecm = level_energy_X(Ev, Jv, Eng, nlev, h2l(il)%vl, h2l(il)%Jl)
     gns = 1.0_wp
     if (mod(h2l(il)%Jl, 2) /= 0) gns = 3.0_wp
     h2l(il)%pop     = gns * (2.0_wp*h2l(il)%Jl+1.0_wp) * exp(-hc_over_k*Ecm/T) / Zpart
     !--- Neufeld dimensionless strength relative to the H I line-center opacity
     h2l(il)%strength = par%f_H2 * h2l(il)%pop * (h2l(il)%f_osc / f_osc_Lya)
     !--- Voigt a-parameter in H2 Doppler units
     h2l(il)%a_damp  = h2l(il)%A_tot_up / (4.0_wp*pi*h2_Dfreq_Hz)
     !--- return-to-Ly-alpha (resonance-scatter) probability
     if (par%h2_pure_absorption) then
        h2l(il)%p_scat = 0.0_wp
     else
        h2l(il)%p_scat = h2l(il)%A_ul / h2l(il)%A_tot_up
     end if
  end do

  if (allocated(W_H2pump)) deallocate(W_H2pump)
  allocate(W_H2pump(n_h2_lines))
  W_H2pump  = 0.0_wp
  par%W_H2abs  = 0.0_wp
  par%W_H2scat = 0.0_wp
  h2_on = .true.

  if (mpar%p_rank == 0) then
     write(*,'(a)')      '--- H2 effect on Ly-alpha enabled (Phase-1 Neufeld two-line) ---'
     write(*,'(a,f9.1,a,es11.3)') '    T_H2 =', par%h2_temperature, ' K   f_H2 =', par%f_H2
     write(*,'(a,es11.3,a)')      '    H2 Doppler width =', h2_Dfreq_Hz, ' Hz'
     do il = 1, n_h2_lines
        write(*,'(a,i1,a,f8.3,a,f8.5,a,es10.3,a,f7.4)') &
           '    line ', il, ': dv=', h2l(il)%dv_kms, ' km/s  f_osc=', h2l(il)%f_osc, &
           '  pop=', h2l(il)%pop, '  p_scat=', h2l(il)%p_scat
     end do
  end if
  end subroutine h2_init

  !=============================================================================
  function h2_kappa(grid, xfreq, icell, jcell, kcell) result(kap)
  !--- H2 line opacity in a cell, expressed as a multiplier of grid%rhokap:
  !---   opacity_H2(x) = grid%rhokap * h2_kappa(x)
  !--- Mirrors calc_voigt_HD: each line is a Voigt in H2 Doppler units, added
  !--- with the Doppler-width ratio Dfreq_H(cell)/Dfreq_H2.
  type(grid_type), intent(in) :: grid
  real(kind=wp),   intent(in) :: xfreq
  integer,         intent(in) :: icell, jcell, kcell
  real(kind=wp) :: kap
  kap = h2_kappa_D(xfreq, grid%Dfreq(icell,jcell,kcell))
  end function h2_kappa

  !=============================================================================
  function h2_select_line(grid, xfreq, icell, jcell, kcell) result(il_sel)
  !--- At an H2 interaction, pick which line absorbed the photon, with
  !--- probability proportional to each line's local opacity contribution.
  use random, only : rand_number
  type(grid_type), intent(in) :: grid
  real(kind=wp),   intent(in) :: xfreq
  integer,         intent(in) :: icell, jcell, kcell
  integer :: il_sel
  il_sel = h2_select_line_D(xfreq, grid%Dfreq(icell,jcell,kcell))
  end function h2_select_line

  !=============================================================================
  subroutine read_energy_X(fname, Ev, Jv, Eng, nlev)
  !--- Read the CLOUDY energy_X.dat file: a magic-number first line, '#' comment
  !--- lines, then rows "V  J  Energy[cm^-1]".
  character(len=*), intent(in) :: fname
  real(kind=wp), allocatable, intent(out) :: Ev(:), Jv(:), Eng(:)
  integer, intent(out) :: nlev
  integer :: unit, ios, n
  real(kind=wp) :: vv, jj, ee
  character(len=256) :: linebuf
  logical :: exists

  inquire(file=fname, exist=exists)
  if (.not. exists) then
     write(*,'(2a)') 'ERROR (h2_mod): cannot find ', trim(fname)
     stop
  end if

  !--- first pass: count data rows
  open(newunit=unit, file=fname, status='old', action='read')
  n = 0
  do
     read(unit,'(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     if (index(linebuf, '//') > 0) cycle          ! magic-number line
     read(linebuf, *, iostat=ios) vv, jj, ee
     if (ios /= 0) cycle
     n = n + 1
  end do
  rewind(unit)

  nlev = n
  allocate(Ev(n), Jv(n), Eng(n))
  n = 0
  do
     read(unit,'(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     if (index(linebuf, '//') > 0) cycle
     read(linebuf, *, iostat=ios) vv, jj, ee
     if (ios /= 0) cycle
     n = n + 1
     Ev(n) = vv; Jv(n) = jj; Eng(n) = ee
  end do
  close(unit)
  end subroutine read_energy_X

  !=============================================================================
  function level_energy_X(Ev, Jv, Eng, nlev, v, J) result(Ecm)
  !--- Look up the X-state energy [cm^-1] of level (v, J).
  real(kind=wp), intent(in) :: Ev(:), Jv(:), Eng(:)
  integer,       intent(in) :: nlev, v, J
  real(kind=wp) :: Ecm
  integer :: i
  Ecm = -1.0_wp
  do i = 1, nlev
     if (nint(Ev(i)) == v .and. nint(Jv(i)) == J) then
        Ecm = Eng(i)
        return
     end if
  end do
  if (Ecm < 0.0_wp) then
     write(*,'(a,i0,a,i0,a)') 'ERROR (h2_mod): X-level (v=',v,', J=',J,') not found in energy_X.dat'
     stop
  end if
  end function level_energy_X

end module h2_mod
