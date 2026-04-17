module setup_mod
contains
  !+++++++++++++++++++++++++++++++++++++++++++
  subroutine read_input
  use define
  use read_mod
  use utility
  use mpi
  implicit none

! local variables
  character(len=128) :: model_infile, arg
  integer :: unit
  integer :: ierr
  !real(kind=wp) :: dx,dy,dz
  integer       :: nx,ny,nz
  real(kind=wp) :: xmax,ymax,zmax
  real(kind=wp) :: vtherm0
  integer :: status

  !--- Read in parameters from params.par using namelist command
  namelist /parameters/ par

  !--- read parameters
  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     write(*,*) 'Usage: ',trim(arg),' input_file.'
     stop
  endif
  !--- newunit specifier is instrodueced in Fortran 2008.
  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,parameters)
  close(unit)

  if (par%no_photons > 1) par%nphotons = par%no_photons
  if (par%no_print   > 1) par%nprint   = par%no_print
  if (par%nprint >= par%nphotons) par%nprint = par%nphotons/10
  if (par%no_photons <= 10) par%nprint = 1

  !--- change the input strings into lower cases (2020.08.30).
  par%geometry        = strlowcase(par%geometry)
  par%source_geometry = strlowcase(par%source_geometry)
  par%velocity_type   = strlowcase(par%velocity_type)
  par%distance_unit   = strlowcase(par%distance_unit)
  par%spectral_type   = strlowcase(par%spectral_type)

  !-- setup line data --------
  call setup_resonance_line()
  !---------------------------

  par%nscatt_dust = 0.0_wp
  par%nscatt_gas  = 0.0_wp
  par%nscatt_tot  = 0.0_wp

  !--- for exoplanet atmosphere (2021.04.22)
  if (trim(par%geometry) == 'plane_atmosphere') then
     par%xy_periodic  = .true.
     par%xyz_symmetry = .false.
     par%xmax         = par%zmax
     par%ymax         = par%zmax
     par%nx           = 1
     par%ny           = 1
  else if (trim(par%geometry) == 'spherical_atmosphere') then
     !--- for now, xyz_symmetry is not allowed. (2021.05.10)
     par%xy_periodic  = .false.
     par%xyz_symmetry = .false.
     if (par%xy_symmetry) then
        par%nx           = maxval([par%nx, par%ny])
        par%ny           = par%nx
        par%xmax         = maxval([par%xmax, par%ymax])
        par%ymax         = par%xmax
     else
        par%nx           = maxval([par%nx, par%ny, par%nz])
        par%ny           = par%nx
        par%nz           = par%nx
        par%xmax         = maxval([par%xmax, par%ymax, par%zmax])
        par%ymax         = par%xmax
        par%zmax         = par%xmax
     endif
  endif

  !--- temperature0 is a temperature of the photon source, which is independent of the cell temperature.
  !--- comment added, 2020.09.02.
  if (par%temperature0 <= 0.0_wp) par%temperature0 = par%temperature
  vtherm0       = line%vtherm1*sqrt(par%temperature0)
  par%Dfreq0    = vtherm0/(line%lambda0*um2km)
  par%voigt_a0  = (line%A21/fourpi)/par%Dfreq0

  !--- continuum parameters
  cpar%slope     = par%continuum_slope
  cpar%wgt0      = par%continuum_wgt0
  cpar%abs_depth = cpar%wgt0 * par%continuum_abs_depth
  cpar%abs_width = par%continuum_abs_width

  !--- bug-fixed, 2020.11.28 (MPI-related variables should be defined before they are used.)
  !--- MPI-related parameters (added SAME_HRANK_COMM, 2020-11-01).
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc,  ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  !The third argument, key, determines the ordering (rank) within each new communicator.
  !call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpar%p_rank, MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_RANK(mpar%hostcomm, mpar%h_rank, ierr)
  call MPI_COMM_SIZE(mpar%hostcomm, mpar%h_nproc, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%p_rank, mpar%SAME_HRANK_COMM, ierr)
  call MPI_COMM_SIZE(mpar%SAME_HRANK_COMM, mpar%SAME_HRANK_NPROC, ierr)

  if (par%nx == 1 .or. par%ny == 1 .or. par%nz == 1) par%xyz_symmetry = .false.
#ifdef PEELINGOFF
  !--- let's do not consider xyz_symmetry condition. (2017-08-19) (2021-05-10).
  !--- later, the xyz_symmetry condition will be implemented.
  !--- for a goemetry with xyz_symmetry, the peeling-off technique can be used only for special cases.
  if (par%xyz_symmetry) then
     par%xyz_symmetry = .false.
     if (mpar%p_rank == 0) then
        write(*,'(a)') '--->'
        write(*,'(a)') '---> xyz_symmetry is not allowed when peelingoff is performed'
        write(*,'(a)') '---> setting xyz_symmetry = .false. and doubling nx, ny, nz.'
        write(*,'(a)') '--->'
     endif
     if ((par%nx/2)*2 == par%nx) then
        par%nx = par%nx * 2
     else
        par%nx = (par%nx-1)*2 + 1
     endif
     if ((par%ny/2)*2 == par%ny) then
        par%ny = par%ny * 2
     else
        par%ny = (par%ny-1)*2 + 1
     endif
     if ((par%nz/2)*2 == par%nz) then
        par%nz = par%nz * 2
     else
        par%nz = (par%nz-1)*2 + 1
     endif
  endif
#endif

  !--- external density, temperature, and velocity field files for TIGRESS model.
  par%input_field  = trim(par%input_field)
  if (len_trim(par%input_field) > 0) then
     par%dens_file = trim(par%input_field)//'.dens.fits.gz'
     par%temp_file = trim(par%input_field)//'.temp.fits.gz'
     par%velo_file = trim(par%input_field)//'.velo.fits.gz'
     call get_dimension(trim(par%dens_file),par%nx,par%ny,par%nz,par%xmax,par%ymax,par%zmax,status,reduce_factor=par%reduce_factor)
     if (status /= 0 .and.  mpar%p_rank == 0) then
        write(*,*) 'Density File ',trim(par%dens_file),' does not have Dimensional Information!'
        write(*,*) 'You need to specify xmax, ymax, zmax, nx, ny, nz in input parameter file!'
     endif
  endif

  !--- DEPRECATION: par%metalZ will be removed later.
  !--- use par%DGR instead.
  if (par%metalZ > 0.0_wp) par%DGR           = par%metalZ

  !--- par%tau0 and par%N_HI are reagared to be par%taumax and par%N_HImax, respectively.
  if (par%tau0 > 0.0_wp .and. par%taumax  < 0.0_wp) par%taumax  = par%tau0
  if (par%N_HI > 0.0_wp .and. par%N_HImax < 0.0_wp) par%N_HImax = par%N_HI
  if (par%N_HImax  > 0.0_wp .and. par%N_gasmax  < 0.0_wp) par%N_gasmax  = par%N_HImax
  if (par%N_HIhomo > 0.0_wp .and. par%N_gashomo < 0.0_wp) par%N_gashomo = par%N_HIhomo

  if (par%cext_dust <= 0.0_wp) par%DGR       = 0.0_wp
  if (par%DGR == 0.0_wp)       par%save_Jabs = .false.
  if (par%core_skip_global)    par%core_skip = .true.

  ! if rmax > 0, then the system is regarded as a sphere or cylinder.
  if (par%nr > 1) then
     par%nx = par%nr
     par%ny = par%nr
     if (trim(par%geometry) /= 'cylinder') par%nz = par%nr
  endif
  if (trim(par%geometry) /= 'cylinder' .and. par%rmax > 0.0_wp) then
     par%xmax = par%rmax
     par%ymax = par%rmax
     par%zmax = par%rmax
     if (.not. par%xy_symmetry) then
        par%nx   = maxval([par%nx,par%ny,par%nz])
        par%ny   = par%nx
        par%nz   = par%nx
     endif
     if (par%source_rmax < 0.0_wp) par%source_rmax = par%rmax
  endif

#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  !--- for Jx(x,y,z) and Pa(x,y,z)
  if (par%geometry_JPa < -1 .or. par%geometry_JPa > 3) then
     par%geometry_JPa  = 1
     if (par%save_all) par%geometry_JPa = 3
     if (par%geometry_JPa /= 3) then
        if (par%xy_periodic) then
           par%geometry_JPa  = -1
        else if (par%xyz_symmetry) then
           par%geometry_JPa  = 1
        else if (trim(par%geometry) == 'cylinder') then
           par%geometry_JPa  = 2
        else if (trim(par%geometry) == 'spherical_atmosphere') then
           par%geometry_JPa  = 2
        endif
     endif
  endif
#endif

  !--- galaxy model in Song, Seon, & Hwang (2020).
  if (trim(par%source_geometry) == 'ssh') then
     par%velocity_type = 'ssh'
     par%sersic_m      = 1.0_wp
     par%Reff          = 1.67834607093866_wp * par%source_rscale
  endif

  if (par%distance2cm < 0.0_wp) then
     select case(trim(par%distance_unit))
        case ('kpc')
           par%distance2cm = kpc2cm
        case ('pc')
           par%distance2cm = pc2cm
        case ('au')
           par%distance2cm = au2cm
        ! 2017-06-27
        case ('')
           par%distance2cm = 1.0_wp
        case default
           par%distance2cm = kpc2cm
     end select
  endif

  if (len_trim(par%scatt_mat_file) > 0 .and. par%DGR > 0.0_wp) then
     call setup_scattering_matrix(trim(par%scatt_mat_file))
  else
     if (par%use_stokes) par%DGR = 0.0_wp
  endif

  if (trim(par%spectral_type) == 'line_prof_file') then
     if (len_trim(par%line_prof_file) > 0) then
        call setup_line_profile(trim(par%line_prof_file))
     else
        par%spectral_type = 'monochromatic'
     endif
  endif

  if (len_trim(par%out_file) == 0) then
     par%base_name = trim(get_base_input_name(model_infile))
     par%out_file  = trim(par%base_name)//'.fits.gz'
  else
     par%base_name = trim(get_base_name(trim(par%out_file)))
  endif

  if (mpar%p_rank == 0) then
     write(*,'(a)')        ''
     write(*,'(3a)')       '+++++ ',trim(model_infile),' +++++'
     write(*,'(a,L1)')     'Use_Master_Slave          : ', par%use_master_slave
     write(*,'(a,a)')      'Line ID                   : ', trim(par%line_id)
     write(*,'(a,2f7.4)')  'Dust Parameters(a, g)     : ', par%albedo, par%hgg
     write(*,'(a,es12.3)') 'Dust Extinction per H     : ', par%cext_dust
     write(*,'(a,es12.3)') 'Dust-to-Gas Ratio         : ', par%DGR
     write(*,'(a,es12.3)') 'Total number of photons   : ', par%no_photons
  endif

  return
  end subroutine read_input
  !---------------------------------------
  subroutine setup_resonance_line()
  use define
  implicit none
  real(kind=wp), parameter :: speedc_cm = 2.99792458e10_wp
  !-- line id is defined using the shorter wavelength of the doublet.
  !-- the central wavelength is the harmonic mean of the two wavelengths.
  !-- the data are obtained from physics.nist.gov.
  !-- cross0  = cross-section at line center in units of cm^2
  !-- sigma_0 = pi x e^2/(m_e c) = pi x c x R_e = 0.02656 (cm^2 Hz)
  select case (trim(par%line_id))
  case ('CIV_1548')
     !-- vacuum wavelengths 1548.187 (short, K, S1/2-P3/2), 1550.772 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1548.202, 1550.774
     !-- Aki = 2.65e8 (short), 2.64e8 (long), fik = 1.90e-1 (short), 9.52e-2 (long)
     line%ion_id    = 'C IV'
     line%lambda0   = 0.1548187_wp
     line%f12       = 0.2852_wp
     line%A21       = 2.647e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(12.0_wp)
     line%DnuHK_Hz  = speedc_cm * (64591.7_wp - 64484.0_wp)
     line%g_recoil0 = (h_planck/massH/12.0_wp)/(line%lambda0*um2m)**2
  case ('NV_1239')
     !-- vacuum wavelengths 1238.821 (short, K, S1/2-P3/2), 1242.804 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1238.921, 1242.804
     !-- Aki = 3.40e8 (short), 3.37e8 (long), fik = 1.56e-1 (short), 7.80e-2 (long)
     line%ion_id    = 'N V'
     line%lambda0   = 0.1238821_wp
     line%f12       = 0.2340_wp
     line%A21       = 3.390e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(14.0_wp)
     line%DnuHK_Hz  = speedc_cm * (80721.9_wp - 80463.2_wp)
     line%g_recoil0 = (h_planck/massH/14.0_wp)/(line%lambda0*um2m)**2
  case ('OVI_1032')
     !-- vacuum wavelengths 1031.912 (short, K, S1/2-P3/2), 1037.613 (long, H, S1/2-P1/2)
     !-- observed wavelengths ??? (not shown in NIST)
     !-- Aki = 4.16e8 (short), 4.09e8 (long), fik = 1.33e-1 (short), 6.60e-2 (long)
     line%ion_id    = 'O VI'
     line%lambda0   = 0.1031912_wp
     line%f12       = 0.1990_wp
     line%A21       = 4.137e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(16.0_wp)
     line%DnuHK_Hz  = speedc_cm * (96907.5_wp - 96375.0_wp)
     line%g_recoil0 = (h_planck/massH/16.0_wp)/(line%lambda0*um2m)**2
  case ('NaI_D')
     !-- vacuum wavelengths 5891.583253 (short, K, S1/2-P3/2), 5897.558147 (long, H, S1/2-P1/2)
     !-- observed wavelengths 5891.583264 (short), 5897.558147 (long)
     !-- Aki = 6.16e7 (short), 6.14e7 (long), fik = 6.41e-1 (short), 3.20e-1 (long)
     line%ion_id    = 'Na I'
     line%lambda0   = 0.5891583253_wp
     line%f12       = 0.9610_wp
     line%A21       = 6.153e7_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(22.0_wp)
     line%DnuHK_Hz  = speedc_cm * (16973.36619_wp - 16956.17025_wp)
     line%g_recoil0 = (h_planck/massH/22.0_wp)/(line%lambda0*um2m)**2
  case ('MgII_2796')
     !-- vacuum wavelengths 2796.352 (short, K, S1/2-P3/2), 2803.531 (long, H, S1/2-P1/2)
     !-- observed wavelengths 2796.352 (short), 2803.530 (long)
     !-- Aki = 2.60e8 (short), 2.57e8 (long), fik = 6.08e-1 (short), 3.03e-1 (long)
     line%ion_id    = 'Mg II'
     line%lambda0   = 0.2796352_wp
     line%f12       = 0.9110_wp
     line%A21       = 2.590e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(24.0_wp)
     line%DnuHK_Hz  = speedc_cm * (35760.88_wp - 35669.31_wp)
     line%g_recoil0 = (h_planck/massH/24.0_wp)/(line%lambda0*um2m)**2
  case ('SiIV_1394')
     !-- vacuum wavelengths 1393.755 (short, K, S1/2-P3/2), 1402.770 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1393.76 (short), 1402.77 (long)
     !-- Aki = 8.80e8 (short), 8.63e8 (long), fik = 5.13e-1 (short), 2.55e-1 (long)
     line%ion_id    = 'Si IV'
     line%lambda0   = 0.1393755_wp
     line%f12       = 0.768_wp
     line%A21       = 8.743e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(28.0_wp)
     line%DnuHK_Hz  = speedc_cm * (71748.64_wp - 71287.54_wp)
     line%g_recoil0 = (h_planck/massH/28.0_wp)/(line%lambda0*um2m)**2
  case ('AlII_1671')
     !-- vacuum wavelengths 1670.7874 (singlet, S0-P1)
     !-- observed wavelengths 1670.7867
     !-- Aki = 1.41e9, fik = 1.77
     line%ion_id    = 'Al II'
     line%lambda0   = 0.16707874_wp
     line%f12       = 1.77_wp
     line%A21       = 1.41e9_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp/sqrt(32.0_wp)
     line%DnuHK_Hz  = 0.0_wp
     line%g_recoil0 = (h_planck/massH/32.0_wp)/(line%lambda0*um2m)**2
  case default
     !-- Lyman-alpha
     !-- line%lambda0 and lambda is in units of micron.
     !-- vacuum wavelengths 1215.668237310 (short, K, S1/2-P3/2), 1215.673644608 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1215.6699 (short), 1215.6699 (long)
     !-- Aki = 6.2648e8 (short), 6.2649e8 (long), fik = 2.7760e-1 (short), 1.3881e-1 (long)
     line%ion_id    = 'H I'
     line%lambda0   = 0.1215668237310_wp
     line%f12       = 0.41641_wp
     line%A21       = 6.26483e8_wp
     line%cross0    = 0.02656_wp/sqrt(pi)*line%f12
     line%vtherm1   = 0.12843374_wp
     line%DnuHK_Hz  = speedc_cm * (82259.2850014_wp - 82258.9191133_wp)
     line%g_recoil0 = (h_planck/massH)/(line%lambda0*um2m)**2
  endselect
  end subroutine setup_resonance_line
  !---------------------------------------
  subroutine setup_scattering_matrix(scatt_mat_file)
  use define
  use mathlib
  use memory_mod
  use mpi
  implicit none
  character(len=*), intent(in) :: scatt_mat_file

  character(len=128) :: string_tmp
  real(kind=wp), allocatable :: phase_CDFi(:)
  real(kind=wp) :: lambda,cext,albedo,hgg
  real(kind=wp) :: S11_norm
  integer       :: i, unit, ierr

  !--- setup Mueller Matrix
  !--- In this code, we will assume spherical dust grains.
  !--- For non-spherical grains, we also need to consider componets (i.e., S22, S21, S13, etc)
  if (len_trim(scatt_mat_file) > 0) then
     if (mpar%h_rank == 0) then
        open(newunit=unit,file=trim(scatt_mat_file),status='old')
        read(unit,*) string_tmp
        read(unit,*) lambda,cext,albedo,hgg,scatt_mat%nPDF
        !read(unit,*) lambda,albedo,hgg,scatt_mat%nPDF
        par%albedo  = albedo
        par%hgg     = hgg
        ! line%lambda0 and lambda is in units of micron.
        !line%lambda0 = lambda
        read(unit,*) string_tmp
     endif
     call MPI_BCAST(par%albedo,     1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(par%hgg,        1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     !call MPI_BCAST(line%lambda0,   1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
     call MPI_BCAST(scatt_mat%nPDF, 1, MPI_INTEGER,          0, mpar%hostcomm, ierr)

     call create_shared_mem(scatt_mat%coss,      [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S11,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S12,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S33,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%S34,       [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%cosp,      [scatt_mat%nPDF])
     call create_shared_mem(scatt_mat%phase_CDF, [scatt_mat%nPDF])

     if (mpar%h_rank == 0) then
        do i=1,scatt_mat%nPDF
           read(unit,*) scatt_mat%coss(i),scatt_mat%S11(i),scatt_mat%S12(i),scatt_mat%S33(i),scatt_mat%S34(i)
        enddo

        ! Normalize the Mueller matrixes
        S11_norm      = calc_Integral(scatt_mat%coss,scatt_mat%S11)
        scatt_mat%S11 = scatt_mat%S11/S11_norm
        scatt_mat%S12 = scatt_mat%S12/S11_norm
        scatt_mat%S33 = scatt_mat%S33/S11_norm
        scatt_mat%S34 = scatt_mat%S34/S11_norm

        ! To make generating a random direction fast, we need to calculate intervals of uniformly-distributed cos(theta).
        if (.not. allocated(phase_CDFi)) allocate(phase_CDFi(scatt_mat%nPDF))
        call calc_CDF(scatt_mat%coss,scatt_mat%S11,phase_CDFi)
        do i=1,scatt_mat%nPDF
           scatt_mat%phase_CDF(i) = real(i-1,kind=wp)/real(scatt_mat%nPDF-1,kind=wp)
           call interp(phase_CDFi,scatt_mat%coss,scatt_mat%phase_CDF(i),scatt_mat%cosp(i))
        enddo
        if (allocated(phase_CDFi)) deallocate(phase_CDFi)
        close(unit)
     endif
  else
     scatt_mat%nPDF = 0
  endif
  end subroutine setup_scattering_matrix
  !---------------------------------------
  subroutine setup_line_profile(line_prof_file)
  use define
  use mathlib
  use memory_mod
  use mpi
  implicit none
  character(len=*), intent(in) :: line_prof_file

  real(kind=wp), allocatable :: xfreq(:), PDF(:), CDF(:)
  real(kind=wp) :: freq,prof,PDF_norm,Dfreq_ref
  integer       :: i, unit, ierr

  if (len_trim(line_prof_file) > 0) then
     if (mpar%h_rank == 0) then
        open(newunit=unit,file=trim(line_prof_file),status='old')
        line_prof%nPDF = 0
        do while(.true.)
           !--- please update this part, depending on the data file format.
           read(unit,*,iostat=ierr) freq,prof
           if (ierr /= 0) exit
           line_prof%nPDF = line_prof%nPDF + 1
        enddo
        close(unit)
     endif

     call MPI_BCAST(line_prof%nPDF, 1, MPI_INTEGER, 0, mpar%hostcomm, ierr)
     call create_shared_mem(line_prof%xfreq, [line_prof%nPDF])
     call create_shared_mem(line_prof%CDF,   [line_prof%nPDF])

     !--- reference temperature and Doppler frequencey.
     Dfreq_ref = line%vtherm1*sqrt(par%temperature)/(line%lambda0*um2km)

     if (mpar%h_rank == 0) then
        if (.not.allocated(xfreq)) allocate(xfreq(line_prof%nPDF))
        if (.not.allocated(PDF))   allocate(PDF(line_prof%nPDF))
        if (.not.allocated(CDF))   allocate(CDF(line_prof%nPDF))

        open(newunit=unit,file=trim(line_prof_file),status='old')
        do i=1,line_prof%nPDF
           !--- please update this part, depending on the data file format.
           !--- speedc is given in km and lambda_LyaH is given in units of micron.
           read(unit,*) freq,prof
           if (par%line_prof_file_type == 1) then
              xfreq(i) = freq
           else
              xfreq(i) = (freq - speedc/(lambda_LyaH*um2km))/Dfreq_ref
           endif
           PDF(i)   = prof
        enddo
        close(unit)

        !--- Normalize the probability distribution function
        PDF_norm = calc_Integral(xfreq,PDF)
        PDF(:)   = PDF(:)/PDF_norm

        !--- To make the random-number generation fast, we need to calculate intervals of uniformly-distributed xfreq.
        call calc_CDF(xfreq,PDF,CDF)
        do i=1,line_prof%nPDF
           line_prof%CDF(i) = real(i-1,kind=wp)/real(line_prof%nPDF-1,kind=wp)
           call interp(CDF,xfreq,line_prof%CDF(i),line_prof%xfreq(i))
        enddo
        if (allocated(xfreq)) deallocate(xfreq)
        if (allocated(PDF))   deallocate(PDF)
        if (allocated(CDF))   deallocate(CDF)
     endif
  else
     line_prof%nPDF = 0
  endif
  end subroutine setup_line_profile
  !+++++++++++++++++++++++++++++++++++++++++++
  subroutine setup_procedure
  use define
  use raytrace
  use scatter_mod
  use run_simulation_mod
  use write_mod
  implicit none

  !--- Initialize Random Number Generator
  call init_random_seed(par%iseed)

  !--- procedure pointer for raytrace routine
  if (par%xyz_symmetry) then
     raytrace_to_tau => raytrace_to_tau_car_xyzsym
  else if (par%xy_symmetry) then
     raytrace_to_tau => raytrace_to_tau_car_xysym
  else if (par%xy_periodic) then
     if (par%nx == 1 .and. par%ny == 1) then
        if (trim(par%geometry) == 'plane_atmosphere') then
           raytrace_to_tau => raytrace_to_tau_car_zonly_atmosphere
        else
           raytrace_to_tau => raytrace_to_tau_car_zonly
        endif
     else
        if (par%Omega /= 0.0_wp) then
           ! 2017-07-14, raytracing routine for the case of shearing boundary.
           raytrace_to_tau => raytrace_to_tau_car_xyper_shear
        else
           raytrace_to_tau => raytrace_to_tau_car_xyper
        endif
     endif
  else
     if (trim(par%geometry) == 'spherical_atmosphere') then
        if (par%xy_symmetry) then
           raytrace_to_tau => raytrace_to_tau_car_xysym_atmosphere
        else
           raytrace_to_tau  => raytrace_to_tau_car_atmosphere
        endif
     else
        raytrace_to_tau  => raytrace_to_tau_car
     endif
  endif

  if (trim(par%geometry) == 'spherical_atmosphere') then
     raytrace_to_edge         => raytrace_to_edge_car_atmosphere
     raytrace_to_edge_tau_gas => raytrace_to_edge_car_tau_gas_atmosphere
     raytrace_to_edge_column  => raytrace_to_edge_car_column_atmosphere
  else
     raytrace_to_edge         => raytrace_to_edge_car
     raytrace_to_edge_tau_gas => raytrace_to_edge_car_tau_gas
     raytrace_to_edge_column  => raytrace_to_edge_car_column
  endif

  !--- procedure pointer for scattering routine
  if (par%use_stokes) then
     scatter_dust      => scatter_dust_stokes
     scatter_resonance => scatter_resonance_stokes
  else
     scatter_dust      => scatter_dust_nostokes
     scatter_resonance => scatter_resonance_nostokes
  endif

  !--- procedure pointer for simulation run
  if (par%use_master_slave) then
     run_simulation => run_master_slave
  else
     run_simulation => run_equal_number
  endif

  !--- procedure pointer to write output
  write_output => write_output_car

  end subroutine setup_procedure
  !+++++++++++++++++++++++++++++++++++++++++++
end module setup_mod
