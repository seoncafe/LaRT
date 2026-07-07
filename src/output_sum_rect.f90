module output_sum_rect
  use define
  use mpi
  use memory_mod
contains
  subroutine output_reduce_outside(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variables
  integer :: k, ierr

  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_dust, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_gas,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nrejected,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (trim(par%source_geometry) == 'stellar_illumination' .or. trim(par%source_geometry) == 'point_illumination') then
     call MPI_ALLREDUCE(MPI_IN_PLACE,par%flux_factor, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     par%flux_factor = par%flux_factor/(par%no_photons + par%nrejected)
  endif
  par%acceptance_rate = par%no_photons/(par%no_photons + par%nrejected)

  call reduce_mem(grid%Jout)
  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     call reduce_mem(grid%Jabs)
  endif
  !--- Ly-beta fluorescence (line_type = 8): band-2 spectra + weight bookkeeping.
  if (line%line_type == 8) then
     if (associated(grid%Jout_Ha)) call reduce_mem(grid%Jout_Ha)
     if (associated(grid%Jabs_Ha)) call reduce_mem(grid%Jabs_Ha)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_conv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_esc1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_abs1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_esc2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_abs2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     call reduce_mem(grid%Jabs2)
  endif
  if (par%save_Jin) then
     call reduce_mem(grid%Jin)
  endif
  if (par%save_Jmu .and. associated(grid%Jmu)) then
     call reduce_mem(grid%Jmu)
  endif

#ifdef CALCJ
  select case (grid%geometry_JPa)
  case (3)
     call reduce_mem(grid%J)
  case (2)
     call reduce_mem(grid%J2)
  case default
     call reduce_mem(grid%J1)
  end select
#endif
#ifdef CALCP
  select case (grid%geometry_JPa)
  case (3)
     call reduce_mem(grid%Pa)
  case (2)
     call reduce_mem(grid%P2)
  case default
     call reduce_mem(grid%P1)
  end select
  !--- Ly-beta (line_type = 8) conversion-rate maps.
  if (line%line_type == 8) then
     if (associated(grid%Pc))  call reduce_mem(grid%Pc)
     if (associated(grid%Pc2)) call reduce_mem(grid%Pc2)
     if (associated(grid%Pc1)) call reduce_mem(grid%Pc1)
  endif
#endif
#ifdef CALCPnew
  select case (grid%geometry_JPa)
  case (3)
     call reduce_mem(grid%Pa_new)
  case (2)
     call reduce_mem(grid%P2_new)
  case default
     call reduce_mem(grid%P1_new)
  end select
#endif

  !--- 2D image
  if (par%save_peeloff_2D) then
     do k=1,par%nobs
        call reduce_mem(observer(k)%scatt_2D)
        call reduce_mem(observer(k)%direc_2D)
        if (par%save_direc0) then
           call reduce_mem(observer(k)%direc0_2D)
        endif
        if (par%use_stokes) then
           call reduce_mem(observer(k)%I_2D)
           call reduce_mem(observer(k)%Q_2D)
           call reduce_mem(observer(k)%U_2D)
           call reduce_mem(observer(k)%V_2D)
        endif
     enddo
  endif

  !--- 3D spectral image
  if (par%save_peeloff_3D) then
     do k=1,par%nobs
        call reduce_mem(observer(k)%scatt)
        call reduce_mem(observer(k)%direc)
        if (associated(observer(k)%peel_Ha)) then
           call reduce_mem(observer(k)%peel_Ha)
        endif
        if (par%save_direc0) then
           call reduce_mem(observer(k)%direc0)
        endif
        if (par%save_dust_scattered) then
           call reduce_mem(observer(k)%scatt_dust)
        endif
        if (par%use_stokes) then
           call reduce_mem(observer(k)%I)
           call reduce_mem(observer(k)%Q)
           call reduce_mem(observer(k)%U)
           call reduce_mem(observer(k)%V)
        endif
     enddo
  endif

  if (par%save_all_photons) then
     !--- Note that allph arrays are declared to be shared memory
     !--- so that "reduce" should be performed only in shared_memory = .true. (2020-11-02).
     call reduce_mem(allph%rp,          shared_memory = .true.)
     call reduce_mem(allph%xfreq1,      shared_memory = .true.)
     call reduce_mem(allph%xfreq2,      shared_memory = .true.)
     call reduce_mem(allph%nscatt_gas,  shared_memory = .true.)
     call reduce_mem(allph%nscatt_dust, shared_memory = .true.)
     if (trim(par%source_geometry) /= 'point') then
        call reduce_mem(allph%rp0, shared_memory = .true.)
     endif
     if (par%use_stokes) then
        call reduce_mem(allph%I, shared_memory = .true.)
        call reduce_mem(allph%Q, shared_memory = .true.)
        call reduce_mem(allph%U, shared_memory = .true.)
        call reduce_mem(allph%V, shared_memory = .true.)
     endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end subroutine output_reduce_outside
!-----------------------------------------------------------
  subroutine output_normalize_outside(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variable
  real(kind=wp) :: dVol, area
  real(kind=wp) :: scale_factor
  real(kind=wp) :: intensity_bin_unit
  real(kind=wp) :: intensity_bin_unit_Ha, nph
  integer       :: i,j,k,ierr

  !--- slightly modified, 2020.09.16
  par%nscatt_dust = par%nscatt_dust/par%nphotons
  par%nscatt_gas  = par%nscatt_gas /par%nphotons
  par%nscatt_tot  = par%nscatt_gas + par%nscatt_dust

  !--- intensity unit.
  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dwave
  else
     intensity_bin_unit = grid%dxfreq
  endif

  if (par%xy_periodic) then
     ! slab geometry
     !-- bug-fixed, 2020.10.23, (distance2cm)**2 was missing.
     area            = grid%xrange*grid%yrange * (par%distance2cm)**2
     dVol            = grid%dx*grid%dy*grid%dz * (par%distance2cm)**2
     ! solid angle = 2pi : all photons are propagating outward only. no photons get scattered back through the surface.
     ! 2.0 is for two (top and bottom) boundaries.
     ! luminosity is assumed to be 1 photons/cm^2.
     ! J(nu,x,y,z) is essentially a sum of path lengths, so that distance2cm should be muliplied to convert to cgs unit.
     !-grid%Jout(:)    = grid%Jout(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
     !-grid%Jin(:)     = grid%Jin(:) /(par%nphotons*grid%dxfreq*twopi*2.0_wp)
     grid%Jout(:)    = grid%Jout(:)/(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
     if (associated(grid%Jin)) &
        grid%Jin(:)  = grid%Jin(:) /(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
     !--- Jmu normalized so that each mu bin = Jout for homogeneous isotropic case
     if (associated(grid%Jmu)) &
        grid%Jmu(:,:) = grid%Jmu(:,:) * par%nmu / (par%nphotons*intensity_bin_unit*twopi*2.0_wp)
  else
     ! sphere or box geometry
     ! luminosity is assumed to be 1 photons/whole volume.
     if (trim(par%geometry) == 'sphere') then
        area = fourpi * grid%rmax**2 * (par%distance2cm)**2
     else
        area = (grid%xmax*grid%ymax + grid%ymax*grid%zmax + grid%zmax*grid%xmax)*8.0_wp * (par%distance2cm)**2
     endif
     dVol            = grid%dx*grid%dy*grid%dz * (par%distance2cm)**2
     !-grid%Jout(:)    = grid%Jout(:)/(par%nphotons*grid%dxfreq*twopi*area)
     !-grid%Jin(:)     = grid%Jin(:) /(par%nphotons*grid%dxfreq*twopi*area)
     grid%Jout(:)    = grid%Jout(:)/(par%nphotons*intensity_bin_unit*twopi*area)
     if (associated(grid%Jin)) &
        grid%Jin(:)  = grid%Jin(:) /(par%nphotons*intensity_bin_unit*twopi*area)
     !--- Jmu normalized so that each mu bin = Jout for homogeneous isotropic case
     if (associated(grid%Jmu)) &
        grid%Jmu(:,:) = grid%Jmu(:,:) * par%nmu / (par%nphotons*intensity_bin_unit*twopi*area)
  endif

  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     if (par%xy_periodic) then
        ! slab geometry
        !-grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
        grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
     else
        ! sphere or box geometry
        !-grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*grid%dxfreq*twopi*area)
        grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*intensity_bin_unit*twopi*area)
     endif
  endif
  !--- Ly-beta fluorescence (line_type = 8): band-2 (H-alpha) spectra.
  !--- Same nphotons/area conventions as band 1, but with the band-2 bin unit
  !--- (dxfreq_Ha / dwave_Ha), so band ratios are luminosity ratios per unit
  !--- bandwidth. xy_periodic is rejected for ly_beta in setup.f90, so the
  !--- sphere/box `area` computed above always applies here.
  intensity_bin_unit_Ha = intensity_bin_unit
  if (line%line_type == 8) then
     if (par%intensity_unit == 1) then
        intensity_bin_unit_Ha = grid%dwave_Ha
     else
        intensity_bin_unit_Ha = grid%dxfreq_Ha
     endif
     if (associated(grid%Jout_Ha)) &
        grid%Jout_Ha(:) = grid%Jout_Ha(:)/(par%nphotons*intensity_bin_unit_Ha*twopi*area)
     if (associated(grid%Jabs_Ha)) &
        grid%Jabs_Ha(:) = grid%Jabs_Ha(:)/(par%nphotons*intensity_bin_unit_Ha*twopi*area)
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     if (par%xy_periodic) then
        ! slab geometry
        !-grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
     else
        ! sphere or box geometry
        !-grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*area)
        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*intensity_bin_unit*twopi*area)
     endif
  endif
  if ((trim(par%spectral_type) == 'continuum' .or. &
       trim(par%spectral_type) == 'continuum+gaussian') .and. par%continuum_normalize) then
     if (.not. associated(grid%Jin)) then
        if (mpar%p_rank == 0) write(*,*) &
           'ERROR: continuum_normalize=T requires save_Jin=T'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     endif
     !--- For continuum+gaussian, the continuum level is (1 - f_line) of the
     !--- average Jin.  Dividing by this gives I/I_continuum = 1 for the
     !--- flat continuum part.
     if (par%f_line > 0.0_wp .and. par%f_line < 1.0_wp) then
        scale_factor = sum(grid%Jin)/size(grid%Jin) * (1.0_wp - par%f_line)
     else
        scale_factor = sum(grid%Jin)/size(grid%Jin)
     endif
     grid%Jout(:) = grid%Jout(:)/scale_factor
     grid%Jin(:)  = grid%Jin(:) /scale_factor
     if (associated(grid%Jabs))  grid%Jabs(:)  = grid%Jabs(:) /scale_factor
     if (associated(grid%Jabs2)) grid%Jabs2(:) = grid%Jabs2(:)/scale_factor
     if (associated(grid%Jmu))   grid%Jmu(:,:) = grid%Jmu(:,:)/scale_factor
     !--- ly_beta band 2: divide by the same band-1 continuum level so the
     !--- H-alpha-to-continuum ratio stays meaningful.
     if (associated(grid%Jout_Ha)) grid%Jout_Ha(:) = grid%Jout_Ha(:)/scale_factor
     if (associated(grid%Jabs_Ha)) grid%Jabs_Ha(:) = grid%Jabs_Ha(:)/scale_factor
  endif

#ifdef CALCJ
  select case (grid%geometry_JPa)
  case (3)
     if (par%xy_periodic) then
        !-grid%J(:,:,:,:) = grid%J(:,:,:,:) * (area/(fourpi*dVol*par%nphotons*grid%dxfreq))
        grid%J(:,:,:,:) = grid%J(:,:,:,:) * (area/(fourpi*dVol*par%nphotons*intensity_bin_unit))
     else
        !-grid%J(:,:,:,:) = grid%J(:,:,:,:) / (fourpi*dVol * par%nphotons*grid%dxfreq)
        grid%J(:,:,:,:) = grid%J(:,:,:,:) / (fourpi*dVol * par%nphotons*intensity_bin_unit)
     endif
     if (par%xyz_symmetry) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           grid%J(:,i,j,k) = grid%J(:,i,j,k)/grid%ncount3D(i,j,k)
        enddo
        enddo
        enddo
     endif
  case (2)
     do k=1,grid%nz
     do j=1,grid%nr
        !-if (grid%ncount_cyl(j,k) > 0) grid%J2(:,j,k) = grid%J2(:,j,k)/grid%ncount_cyl(j,k) / (fourpi*dVol*par%nphotons*grid%dxfreq)
        if (grid%ncount_cyl(j,k) > 0) grid%J2(:,j,k) = grid%J2(:,j,k)/grid%ncount_cyl(j,k) / (fourpi*dVol*par%nphotons*intensity_bin_unit)
     enddo
     enddo
  case (1)
     !--- Note that the factors 2 and 8 in the case of xyz_symmetry are considered in add_to_J, add_to_Pnew and add_to_Pa.
     do k=1,grid%nr
        !-if (grid%ncount_sph(k) > 0) grid%J1(:,k) = grid%J1(:,k)/grid%ncount_sph(k) / (fourpi*dVol*par%nphotons*grid%dxfreq)
        if (grid%ncount_sph(k) > 0) grid%J1(:,k) = grid%J1(:,k)/grid%ncount_sph(k) / (fourpi*dVol*par%nphotons*intensity_bin_unit)
     enddo
  case (-1)
     !grid%J1 = grid%J1/(grid%nx*grid%ny) * (area/(fourpi*dVol*par%nphotons*grid%dxfreq))
     !grid%J1 = grid%J1/(grid%nx*grid%ny) * (area/(fourpi*dVol*par%nphotons*intensity_bin_unit))
     do k=1,grid%nz
        !-if (grid%ncount_plane(k) > 0) grid%J1(:,k) = grid%J1(:,k)/grid%ncount_plane(k)*(area/(fourpi*dVol*par%nphotons*grid%dxfreq))
        if (grid%ncount_plane(k) > 0) grid%J1(:,k) = grid%J1(:,k)/grid%ncount_plane(k)*(area/(fourpi*dVol*par%nphotons*intensity_bin_unit))
     enddo
  end select
#endif

#ifdef CALCP
  select case (grid%geometry_JPa)
  case (3)
     if (par%xy_periodic) then
        grid%Pa(:,:,:) = grid%Pa(:,:,:)*(area/(dVol*par%nphotons))
     else
        grid%Pa(:,:,:) = grid%Pa(:,:,:)/(dVol*par%nphotons)
     endif
     if (par%xyz_symmetry) then
        grid%Pa(:,:,:)  = grid%Pa(:,:,:) / grid%ncount3D(:,:,:)
     endif
  case (2)
     do k=1,grid%nz
     do j=1,grid%nr
        if (grid%ncount_cyl(j,k) > 0) grid%P2(j,k) = grid%P2(j,k)/grid%ncount_cyl(j,k) / (dVol * par%nphotons)
     enddo
     enddo
  case (1)
     !--- Note that the factors 2 and 8 in the case of xyz_symmetry are considered in add_to_J, add_to_Pnew and add_to_Pa.
     do k=1,grid%nr
        if (grid%ncount_sph(k) > 0) grid%P1(k) = grid%P1(k)/ grid%ncount_sph(k) / (dVol * par%nphotons)
     enddo
  case (-1)
     !grid%P1 = grid%P1/(grid%nx*grid%ny) * (area/(dVol*par%nphotons))
     do k=1,grid%nz
        if (grid%ncount_plane(k) > 0) grid%P1(k) = grid%P1(k)/grid%ncount_plane(k) * (area/(dVol*par%nphotons))
     enddo
  end select
  !--- Ly-beta (line_type = 8) conversion-rate maps: identical normalization
  !--- to Pa/P1/P2 (per-atom rate). Expectation: P_conv/Pa -> 0.11834.
  if (line%line_type == 8) then
     select case (grid%geometry_JPa)
     case (3)
        if (par%xy_periodic) then
           grid%Pc(:,:,:) = grid%Pc(:,:,:)*(area/(dVol*par%nphotons))
        else
           grid%Pc(:,:,:) = grid%Pc(:,:,:)/(dVol*par%nphotons)
        endif
        if (par%xyz_symmetry) then
           grid%Pc(:,:,:)  = grid%Pc(:,:,:) / grid%ncount3D(:,:,:)
        endif
     case (2)
        do k=1,grid%nz
        do j=1,grid%nr
           if (grid%ncount_cyl(j,k) > 0) grid%Pc2(j,k) = grid%Pc2(j,k)/grid%ncount_cyl(j,k) / (dVol * par%nphotons)
        enddo
        enddo
     case (1)
        do k=1,grid%nr
           if (grid%ncount_sph(k) > 0) grid%Pc1(k) = grid%Pc1(k)/ grid%ncount_sph(k) / (dVol * par%nphotons)
        enddo
     case (-1)
        do k=1,grid%nz
           if (grid%ncount_plane(k) > 0) grid%Pc1(k) = grid%Pc1(k)/grid%ncount_plane(k) * (area/(dVol*par%nphotons))
        enddo
     end select
  endif
#endif

#ifdef CALCPnew
  select case (grid%geometry_JPa)
  case (3)
     if (par%xy_periodic) then
        grid%Pa_new(:,:,:) = grid%Pa_new(:,:,:)*(area/(dVol*par%nphotons))
     else
        grid%Pa_new(:,:,:) = grid%Pa_new(:,:,:)/(dVol*par%nphotons)
     endif
     if (par%xyz_symmetry) then
        grid%Pa_new(:,:,:)  = grid%Pa_new(:,:,:) / grid%ncount3D(:,:,:)
     endif
  case (2)
     do k=1,grid%nz
     do j=1,grid%nr
        if (grid%ncount_cyl(j,k) > 0) grid%P2_new(j,k) = grid%P2_new(j,k)/grid%ncount_cyl(j,k) / (dVol * par%nphotons)
     enddo
     enddo
  case (1)
     !--- Note that the factors 2 and 8 in the case of xyz_symmetry are considered in add_to_J, add_to_Pnew and add_to_Pa.
     do k=1,grid%nr
        if (grid%ncount_sph(k) > 0) grid%P1_new(k) = grid%P1_new(k)/ grid%ncount_sph(k) / (dVol * par%nphotons)
     enddo
  case (-1)
     !grid%P1_new = grid%P1_new/(grid%nx*grid%ny) * (area/(dVol*par%nphotons))
     do k=1,grid%nz
        if (grid%ncount_plane(k) > 0) grid%P1_new(k) = grid%P1_new(k)/grid%ncount_plane(k) * (area/(dVol*par%nphotons))
     enddo
  end select
#endif

  !--- 2D image
  if (par%save_peeloff_2D) then
     do k=1,par%nobs
        scale_factor = par%no_photons*observer(k)%steradian_pix * par%distance2cm**2
        observer(k)%scatt_2D(:,:)  = observer(k)%scatt_2D(:,:) / scale_factor
        observer(k)%direc_2D(:,:)  = observer(k)%direc_2D(:,:) / scale_factor

        if (par%save_direc0) then
           observer(k)%direc0_2D(:,:) = observer(k)%direc0_2D(:,:)/ scale_factor
        endif
        if (par%use_stokes) then
           observer(k)%I_2D(:,:) = observer(k)%I_2D(:,:) / scale_factor
           observer(k)%Q_2D(:,:) = observer(k)%Q_2D(:,:) / scale_factor
           observer(k)%U_2D(:,:) = observer(k)%U_2D(:,:) / scale_factor
           observer(k)%V_2D(:,:) = observer(k)%V_2D(:,:) / scale_factor
        endif
     enddo
  endif
  !--- 3D spectral image
  if (par%save_peeloff_3D) then
     do k=1,par%nobs
        !-scale_factor = par%no_photons*observer(k)%steradian_pix*grid%dxfreq * par%distance2cm**2
        scale_factor = par%no_photons*observer(k)%steradian_pix*intensity_bin_unit * par%distance2cm**2
        observer(k)%scatt(:,:,:)  = observer(k)%scatt(:,:,:) / scale_factor
        observer(k)%direc(:,:,:)  = observer(k)%direc(:,:,:) / scale_factor

        !--- ly_beta band-2 (H-alpha) peel cube: same convention, band-2 bin unit.
        if (associated(observer(k)%peel_Ha)) then
           observer(k)%peel_Ha(:,:,:) = observer(k)%peel_Ha(:,:,:) / &
              (par%no_photons*observer(k)%steradian_pix*intensity_bin_unit_Ha * par%distance2cm**2)
        endif
        if (par%save_direc0) then
           observer(k)%direc0(:,:,:) = observer(k)%direc0(:,:,:)/ scale_factor
        endif
        if (par%save_dust_scattered) then
           observer(k)%scatt_dust(:,:,:) = observer(k)%scatt_dust(:,:,:)/ scale_factor
        endif
        if (par%use_stokes) then
           observer(k)%I(:,:,:) = observer(k)%I(:,:,:) / scale_factor
           observer(k)%Q(:,:,:) = observer(k)%Q(:,:,:) / scale_factor
           observer(k)%U(:,:,:) = observer(k)%U(:,:,:) / scale_factor
           observer(k)%V(:,:,:) = observer(k)%V(:,:,:) / scale_factor
        endif
     enddo
  endif

  !--- Ly-beta fluorescence: per-band weight bookkeeping (per source photon).
  !--- Conservation: band1_esc + band1_abs + conv = 1 ; band2_esc + band2_abs = conv.
  if (line%line_type == 8 .and. mpar%p_rank == 0) then
     nph = dble(par%nphotons)
     write(*,'(a)')        '--- ly_beta bookkeeping (weights per source photon) ---'
     write(*,'(a,es14.6)') '  band-1 (Ly-beta) escaped        : ', par%W_esc1/nph
     write(*,'(a,es14.6)') '  band-1 (Ly-beta) dust-absorbed  : ', par%W_abs1/nph
     write(*,'(a,es14.6)') '  conversions (3p->2s)            : ', par%W_conv/nph
     write(*,'(a,es14.6)') '  band-2 (H-alpha) escaped        : ', par%W_esc2/nph
     write(*,'(a,es14.6)') '  band-2 (H-alpha) dust-absorbed  : ', par%W_abs2/nph
     write(*,'(a,es14.6)') '  band1_esc + band1_abs + conv (=1)     : ', &
                           (par%W_esc1 + par%W_abs1 + par%W_conv)/nph
     write(*,'(a,es14.6)') '  band2_esc + band2_abs (= conv)        : ', &
                           (par%W_esc2 + par%W_abs2)/nph
     write(*,'(a)')        '--------------------------------------------------------'
  endif

  end subroutine output_normalize_outside
!--------------------------------------------------------------
  subroutine make_radial_intensity(grid,obs,use_2D_data)
  implicit none
  type(grid_type),     intent(in)    :: grid
  type(observer_type), intent(inout) :: obs
  logical, optional,   intent(in)    :: use_2D_data

  logical              :: using_2D_data
  integer              :: i,j,ir,nr
  real(kind=wp)        :: xcen,ycen
  real(kind=wp)        :: xx,yy,rr
  real(kind=wp)        :: intensity_bin_unit
  integer, allocatable :: ncount(:)

  using_2D_data = .false.
  if (present(use_2D_data)) using_2D_data = use_2D_data

  !--- intensity unit.
  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dwave
  else
     intensity_bin_unit = grid%dxfreq
  endif

  nr = (maxval([obs%nxim,obs%nyim])+1)/2
  if (.not. associated(obs%radial_I)) allocate(obs%radial_I(nr))
  if (.not. allocated(ncount))        allocate(ncount(nr))
  obs%radial_I(:) = 0.0_wp
  ncount(:)       = 0
  if (.not. associated(obs%radial_r)) allocate(obs%radial_r(nr))
  if ((nr/2)*2 == nr) then
     obs%radial_r(:) = [((i-0.5_wp)/nr,          i=1,nr)]
  else
     obs%radial_r(:) = [((i-1.0_wp)/(nr-0.5_wp), i=1,nr)]
  endif

  xcen = (obs%nxim+1.0_wp)/2.0_wp
  ycen = (obs%nxim+1.0_wp)/2.0_wp
  do j=1,obs%nyim
     yy = j - ycen
     do i=1,obs%nxim
        xx = i - xcen
        rr = sqrt(xx**2 + yy**2)
        if ((nr/2)*2 == nr) then
           ir = floor(rr) + 1
        else
           ir = floor(rr+0.5_wp) + 1
        endif
        if (ir >= 1 .and. ir <= nr) then
           ncount(ir) = ncount(ir) + 1
           if (using_2D_data) then
              obs%radial_I(ir) = obs%radial_I(ir) + obs%scatt_2D(i,j) + obs%direc_2D(i,j)
           else
              !--- integral over the frequency range (2020.09.30).
              !-obs%radial_I(ir) = obs%radial_I(ir) + (sum(obs%scatt(:,i,j)) + sum(obs%direc(:,i,j))) * grid%dxfreq
              obs%radial_I(ir) = obs%radial_I(ir) + (sum(obs%scatt(:,i,j)) + sum(obs%direc(:,i,j))) * intensity_bin_unit
           endif
        endif
     enddo
  enddo

  do ir=1,nr
     if (ncount(ir) > 0) obs%radial_I(ir) = obs%radial_I(ir)/ncount(ir)
  enddo
  if (allocated(ncount)) deallocate(ncount)
  end subroutine make_radial_intensity
!--------------------------------------------------------------
  subroutine make_radial_stokes(grid,obs,use_2D_data)
  implicit none
  type(grid_type),     intent(in)    :: grid
  type(observer_type), intent(inout) :: obs
  logical, optional,   intent(in)    :: use_2D_data

  logical              :: using_2D_data
  integer              :: i,j,ir,nr
  real(kind=wp)        :: xcen,ycen
  real(kind=wp)        :: xx,yy,rr
  real(kind=wp)        :: cosp,sinp,cos2p,sin2p
  real(kind=wp)        :: intensity_bin_unit
  integer, allocatable :: ncount(:)

  using_2D_data = .false.
  if (present(use_2D_data)) using_2D_data = use_2D_data

  !--- intensity unit.
  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dwave
  else
     intensity_bin_unit = grid%dxfreq
  endif

  nr = (maxval([obs%nxim,obs%nyim])+1)/2
  if (.not. associated(obs%radial_pol)) allocate(obs%radial_pol(nr))
  if (.not. associated(obs%radial_r))   allocate(obs%radial_r(nr))
  if (.not. associated(obs%radial_I))   allocate(obs%radial_I(nr))
  if (.not. associated(obs%radial_Q))   allocate(obs%radial_Q(nr))
  if (.not. associated(obs%radial_U))   allocate(obs%radial_U(nr))
  if (.not. associated(obs%radial_V))   allocate(obs%radial_V(nr))
  if (.not. allocated(ncount))          allocate(ncount(nr))
  obs%radial_pol(:) = 0.0_wp
  obs%radial_r(:)   = 0.0_wp
  obs%radial_I(:)   = 0.0_wp
  obs%radial_Q(:)   = 0.0_wp
  obs%radial_U(:)   = 0.0_wp
  obs%radial_V(:)   = 0.0_wp
  ncount(:)         = 0

  if ((nr/2)*2 == nr) then
     obs%radial_r(:) = [((i-0.5_wp)/nr,          i=1,nr)]
  else
     obs%radial_r(:) = [((i-1.0_wp)/(nr-0.5_wp), i=1,nr)]
  endif

  xcen = (obs%nxim+1.0_wp)/2.0_wp
  ycen = (obs%nxim+1.0_wp)/2.0_wp
  do j=1,obs%nyim
     yy = j - ycen
     do i=1,obs%nxim
        xx = i - xcen
        rr = sqrt(xx**2 + yy**2)
        if ((nr/2)*2 == nr) then
           ir = floor(rr) + 1
        else
           ir = floor(rr+0.5_wp) + 1
        endif

        if (ir >= 1 .and. ir <= nr) then
           if (rr == 0.0_wp) then
              cosp = 1.0_wp
              sinp = 0.0_wp
           else
              cosp  =  yy/rr
              sinp  = -xx/rr
           endif
           cos2p = 2.0_wp * cosp**2 - 1.0_wp
           sin2p = 2.0_wp * cosp * sinp
           ncount(ir)       = ncount(ir) + 1
           if (using_2D_data) then
              obs%radial_I(ir) = obs%radial_I(ir) + obs%I_2D(i,j)
              obs%radial_Q(ir) = obs%radial_Q(ir) + obs%Q_2D(i,j)*cos2p + obs%U_2D(i,j)*sin2p
              obs%radial_U(ir) = obs%radial_U(ir) - obs%Q_2D(i,j)*sin2p + obs%U_2D(i,j)*cos2p
              obs%radial_V(ir) = obs%radial_V(ir) + obs%V_2D(i,j)
           else
              !--- integral over the frequency range (2020.09.30).
              !-obs%radial_I(ir) = obs%radial_I(ir) + sum(obs%I(:,i,j)) * grid%dxfreq
              !-obs%radial_Q(ir) = obs%radial_Q(ir) + sum( obs%Q(:,i,j)*cos2p + obs%U(:,i,j)*sin2p) * grid%dxfreq
              !-obs%radial_U(ir) = obs%radial_U(ir) + sum(-obs%Q(:,i,j)*sin2p + obs%U(:,i,j)*cos2p) * grid%dxfreq
              !-obs%radial_V(ir) = obs%radial_V(ir) + sum(obs%V(:,i,j)) * grid%dxfreq
              obs%radial_I(ir) = obs%radial_I(ir) + sum(obs%I(:,i,j)) * intensity_bin_unit
              obs%radial_Q(ir) = obs%radial_Q(ir) + sum( obs%Q(:,i,j)*cos2p + obs%U(:,i,j)*sin2p) * intensity_bin_unit
              obs%radial_U(ir) = obs%radial_U(ir) + sum(-obs%Q(:,i,j)*sin2p + obs%U(:,i,j)*cos2p) * intensity_bin_unit
              obs%radial_V(ir) = obs%radial_V(ir) + sum(obs%V(:,i,j)) * intensity_bin_unit
           endif
        endif
     enddo
  enddo

  do ir=1,nr
     if (ncount(ir) > 0) then
        obs%radial_I(ir)   = obs%radial_I(ir) / ncount(ir)
        obs%radial_Q(ir)   = obs%radial_Q(ir) / ncount(ir)
        obs%radial_U(ir)   = obs%radial_U(ir) / ncount(ir)
        obs%radial_V(ir)   = obs%radial_V(ir) / ncount(ir)
        !--- ncount is a number of pixels and hence radial_I = 0 can happens (2020.11.08).
        if (obs%radial_I(ir) > 0.0_wp) then
           obs%radial_pol(ir) = sqrt(obs%radial_Q(ir)**2 + obs%radial_U(ir)**2)/obs%radial_I(ir)
        endif
     endif
  enddo
  if (allocated(ncount)) deallocate(ncount)
  end subroutine make_radial_stokes
!--------------------------------------------------------------
  subroutine output_reduce_amr(grid)
  use octree_mod, only: amr_grid
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: ierr, k

  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nscatt_dust,      1,               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nscatt_gas,       1,               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nrejected,        1,               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jout(1),     amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jin)) &
      call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jin(1),  amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jabs)) &
      call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jabs(1), amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jmu)) &
      call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jmu(1,1), amr_grid%nxfreq*size(amr_grid%Jmu,2), &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  !--- Ly-beta fluorescence (line_type = 8): band-2 spectra + weight bookkeeping.
  if (line%line_type == 8) then
     if (associated(amr_grid%Jout_Ha)) &
         call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jout_Ha(1), amr_grid%nxfreq_Ha, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     if (associated(amr_grid%Jabs_Ha)) &
         call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jabs_Ha(1), amr_grid%nxfreq_Ha, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_conv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_esc1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_abs1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_esc2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, par%W_abs2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  endif

#ifdef CALCJ
  select case (amr_grid%geometry_JPa)
     case (3);     call reduce_mem(amr_grid%J)
     case (2);     call reduce_mem(amr_grid%J2)
     case default; call reduce_mem(amr_grid%J1)
  end select
#endif
#ifdef CALCP
  select case (amr_grid%geometry_JPa)
     case (3);     call reduce_mem(amr_grid%Pa)
     case (2);     call reduce_mem(amr_grid%P2)
     case default; call reduce_mem(amr_grid%P1)
  end select
  !--- Ly-beta (line_type = 8) conversion-rate maps.
  if (line%line_type == 8) then
     if (associated(amr_grid%Pc))  call reduce_mem(amr_grid%Pc)
     if (associated(amr_grid%Pc2)) call reduce_mem(amr_grid%Pc2)
     if (associated(amr_grid%Pc1)) call reduce_mem(amr_grid%Pc1)
  endif
#endif
#ifdef CALCPnew
  select case (amr_grid%geometry_JPa)
     case (3);     call reduce_mem(amr_grid%Pa_new)
     case (2);     call reduce_mem(amr_grid%P2_new)
     case default; call reduce_mem(amr_grid%P1_new)
  end select
#endif

  ! Copy reduced amr_grid spectra into grid pointers for write_output
  grid%Jout = amr_grid%Jout
  if (allocated(amr_grid%Jin))  grid%Jin  = amr_grid%Jin
  if (allocated(amr_grid%Jabs)) grid%Jabs = amr_grid%Jabs
  if (allocated(amr_grid%Jmu))  grid%Jmu  = amr_grid%Jmu
  !--- Ly-beta (line_type = 8): point grid band-2 pointers at reduced amr_grid data.
  if (line%line_type == 8) then
     if (associated(amr_grid%Jout_Ha)) grid%Jout_Ha => amr_grid%Jout_Ha
     if (associated(amr_grid%Jabs_Ha)) grid%Jabs_Ha => amr_grid%Jabs_Ha
  endif

  ! Peel-off observer array reductions (same as output_reduce_outside)
  if (par%save_peeloff_2D) then
     do k = 1, par%nobs
        call reduce_mem(observer(k)%scatt_2D)
        call reduce_mem(observer(k)%direc_2D)
        if (par%save_direc0)  call reduce_mem(observer(k)%direc0_2D)
        if (par%use_stokes) then
           call reduce_mem(observer(k)%I_2D)
           call reduce_mem(observer(k)%Q_2D)
           call reduce_mem(observer(k)%U_2D)
           call reduce_mem(observer(k)%V_2D)
        end if
     end do
  end if
  if (par%save_peeloff_3D) then
     do k = 1, par%nobs
        call reduce_mem(observer(k)%scatt)
        call reduce_mem(observer(k)%direc)
        if (associated(observer(k)%peel_Ha)) call reduce_mem(observer(k)%peel_Ha)
        if (par%save_direc0)         call reduce_mem(observer(k)%direc0)
        if (par%save_dust_scattered) call reduce_mem(observer(k)%scatt_dust)
        if (par%use_stokes) then
           call reduce_mem(observer(k)%I)
           call reduce_mem(observer(k)%Q)
           call reduce_mem(observer(k)%U)
           call reduce_mem(observer(k)%V)
        end if
     end do
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine output_reduce_amr
!--------------------------------------------------------------
  subroutine output_normalize_amr(grid)
  use octree_mod, only: amr_grid
  implicit none
  type(grid_type), intent(inout) :: grid
  real(wp) :: area, norm_out, scale_factor, intensity_bin_unit
  real(wp) :: intensity_bin_unit_Ha, norm_out_Ha
  integer  :: k, ierr
#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  real(wp) :: d2
  integer  :: il, ir, iz
#endif

  par%nscatt_dust = par%nscatt_dust / par%nphotons
  par%nscatt_gas  = par%nscatt_gas  / par%nphotons
  par%nscatt_tot  = par%nscatt_gas  + par%nscatt_dust

  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dwave
  else
     intensity_bin_unit = grid%dxfreq
  end if

  ! Geometry-aware surface-area / solid-angle normalization, mirroring the
  ! Cartesian output_normalize_outside convention so that AMR and Cartesian
  ! runs with the SAME par%geometry agree:
  !   par%xy_periodic        -> slab: 2 faces (top+bottom), no area factor
  !   trim(par%geometry)=='sphere' -> area = 4*pi*rmax^2  (rmax = L_box/2)
  !   otherwise ('rectangle'/'box') -> 6-face box area
  ! NOTE: setup.f90 maps the par%geometry='' default to 'sphere', so default
  ! AMR runs keep the sphere convention (unchanged from before this fix); the
  ! box branch is reached only for an explicit par%geometry='box'/'rectangle'.
  ! The box dimensions come from the octree grid (amr_grid), not the grid_type
  ! argument.  A single divisor `norm_out` is valid for Jout/Jin/Jabs/Jmu.
  if (par%xy_periodic) then
     ! slab geometry: 2.0 = two (top+bottom) boundaries, solid angle 2*pi
     norm_out = par%nphotons * intensity_bin_unit * twopi * 2.0_wp
  else
     if (trim(par%geometry) == 'sphere') then
        area = fourpi * par%rmax**2 * par%distance2cm**2
     else
        ! 6-face area of the AMR box; half-widths are xrange/2, ... so that
        ! (hx*hy + hy*hz + hz*hx)*8 == 2*(xrange*yrange + yrange*zrange + zrange*xrange)
        area = 2.0_wp*(amr_grid%xrange*amr_grid%yrange + amr_grid%yrange*amr_grid%zrange + &
                       amr_grid%zrange*amr_grid%xrange) * par%distance2cm**2
     endif
     norm_out = par%nphotons * intensity_bin_unit * twopi * area
  endif

  grid%Jout(:) = grid%Jout(:) / norm_out
  if (associated(grid%Jin)) &
      grid%Jin(:)  = grid%Jin(:)  / norm_out
  if (associated(grid%Jabs)) &
      grid%Jabs(:) = grid%Jabs(:) / norm_out
  if (associated(grid%Jmu)) &
      grid%Jmu(:,:) = grid%Jmu(:,:) * par%nmu / norm_out

  !--- Ly-beta (line_type = 8): band-2 (H-alpha) spectra. Same nphotons/area
  !--- convention as band 1, but with the band-2 bin unit so band ratios are
  !--- luminosity ratios per unit bandwidth (mirrors output_normalize_outside).
  intensity_bin_unit_Ha = intensity_bin_unit
  if (line%line_type == 8) then
     if (par%intensity_unit == 1) then
        intensity_bin_unit_Ha = grid%dwave_Ha
     else
        intensity_bin_unit_Ha = grid%dxfreq_Ha
     end if
     norm_out_Ha = norm_out * (intensity_bin_unit_Ha / intensity_bin_unit)
     if (associated(grid%Jout_Ha)) grid%Jout_Ha(:) = grid%Jout_Ha(:) / norm_out_Ha
     if (associated(grid%Jabs_Ha)) grid%Jabs_Ha(:) = grid%Jabs_Ha(:) / norm_out_Ha
  end if

  ! Continuum normalization (same logic as output_normalize_outside)
  if ((trim(par%spectral_type) == 'continuum' .or. &
       trim(par%spectral_type) == 'continuum+gaussian') .and. par%continuum_normalize) then
     if (.not. associated(grid%Jin)) then
        if (mpar%p_rank == 0) write(*,*) &
           'ERROR: continuum_normalize=T requires save_Jin=T'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     endif
     if (par%f_line > 0.0_wp .and. par%f_line < 1.0_wp) then
        scale_factor = sum(grid%Jin)/size(grid%Jin) * (1.0_wp - par%f_line)
     else
        scale_factor = sum(grid%Jin)/size(grid%Jin)
     endif
     grid%Jout(:) = grid%Jout(:)/scale_factor
     grid%Jin(:)  = grid%Jin(:) /scale_factor
     if (associated(grid%Jabs))  grid%Jabs(:)  = grid%Jabs(:) /scale_factor
     if (associated(grid%Jabs2)) grid%Jabs2(:) = grid%Jabs2(:)/scale_factor
     if (associated(grid%Jmu))   grid%Jmu(:,:) = grid%Jmu(:,:)/scale_factor
     !--- ly_beta band 2: divide by the same band-1 continuum level.
     if (associated(grid%Jout_Ha)) grid%Jout_Ha(:) = grid%Jout_Ha(:)/scale_factor
     if (associated(grid%Jabs_Ha)) grid%Jabs_Ha(:) = grid%Jabs_Ha(:)/scale_factor
  endif

  !--- CALCJ / CALCP / CALCPnew : volume-weighted normalization on AMR leaves.
  !    Operates on the module-global amr_grid (this routine takes a grid_type
  !    argument only for the spectrum copy).  Mirrors the Cartesian convention:
  !    the per-bin sum of leaf volumes (vol_*) replaces ncount_*, and d2 plays
  !    the role of the Cartesian dVol's distance factor.  nadd is pre-folded.
#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  d2 = par%distance2cm**2
#endif
#ifdef CALCJ
  select case (amr_grid%geometry_JPa)
  case (3)
     !--- xy_periodic (slab): luminosity is 1 photon/cm^2, so multiply by the
     !--- slab area (Cartesian: J*(area/(fourpi*dVol*...)); d2 cancels here).
     do il = 1, amr_grid%nleaf
        if (amr_grid%vol_leaf(il) > 0.0_wp) then
           if (par%xy_periodic) then
              amr_grid%J(:,il) = amr_grid%J(:,il) * (amr_grid%xrange*amr_grid%yrange) / &
                 (fourpi * amr_grid%vol_leaf(il) * par%nphotons * intensity_bin_unit)
           else
              amr_grid%J(:,il) = amr_grid%J(:,il) / &
                 (fourpi * amr_grid%vol_leaf(il) * d2 * par%nphotons * intensity_bin_unit)
           end if
        end if
     end do
  case (2)
     do iz = 1, amr_grid%nz_JPa
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_cyl(ir,iz) > 0.0_wp) &
           amr_grid%J2(:,ir,iz) = amr_grid%J2(:,ir,iz) / &
              (fourpi * amr_grid%vol_cyl(ir,iz) * d2 * par%nphotons * intensity_bin_unit)
     end do
     end do
  case (-1)
     do iz = 1, amr_grid%nz_JPa
        if (amr_grid%vol_plane(iz) > 0.0_wp) &
           amr_grid%J1(:,iz) = amr_grid%J1(:,iz) * (amr_grid%xrange*amr_grid%yrange) / &
              (amr_grid%vol_plane(iz) * fourpi * par%nphotons * intensity_bin_unit)
     end do
  case default   ! geometry_JPa == 1 (radial / sphere)
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_sph(ir) > 0.0_wp) &
           amr_grid%J1(:,ir) = amr_grid%J1(:,ir) / &
              (fourpi * amr_grid%vol_sph(ir) * d2 * par%nphotons * intensity_bin_unit)
     end do
  end select
#endif
#ifdef CALCP
  select case (amr_grid%geometry_JPa)
  case (3)
     do il = 1, amr_grid%nleaf
        if (amr_grid%vol_leaf(il) > 0.0_wp) then
           if (par%xy_periodic) then
              amr_grid%Pa(il) = amr_grid%Pa(il) * (amr_grid%xrange*amr_grid%yrange) / &
                 (amr_grid%vol_leaf(il) * par%nphotons)
           else
              amr_grid%Pa(il) = amr_grid%Pa(il) / (amr_grid%vol_leaf(il) * d2 * par%nphotons)
           end if
        end if
     end do
  case (2)
     do iz = 1, amr_grid%nz_JPa
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_cyl(ir,iz) > 0.0_wp) &
           amr_grid%P2(ir,iz) = amr_grid%P2(ir,iz) / (amr_grid%vol_cyl(ir,iz) * d2 * par%nphotons)
     end do
     end do
  case (-1)
     do iz = 1, amr_grid%nz_JPa
        if (amr_grid%vol_plane(iz) > 0.0_wp) &
           amr_grid%P1(iz) = amr_grid%P1(iz) * (amr_grid%xrange*amr_grid%yrange) / &
              (amr_grid%vol_plane(iz) * par%nphotons)
     end do
  case default   ! geometry_JPa == 1 (radial / sphere)
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_sph(ir) > 0.0_wp) &
           amr_grid%P1(ir) = amr_grid%P1(ir) / (amr_grid%vol_sph(ir) * d2 * par%nphotons)
     end do
  end select
  !--- Ly-beta (line_type = 8) conversion-rate maps: identical per-atom-rate
  !--- volume-weighted normalization as Pa (expectation P_conv/Pa -> 0.11834).
  if (line%line_type == 8) then
     select case (amr_grid%geometry_JPa)
     case (3)
        do il = 1, amr_grid%nleaf
           if (amr_grid%vol_leaf(il) > 0.0_wp) then
              if (par%xy_periodic) then
                 amr_grid%Pc(il) = amr_grid%Pc(il) * (amr_grid%xrange*amr_grid%yrange) / &
                    (amr_grid%vol_leaf(il) * par%nphotons)
              else
                 amr_grid%Pc(il) = amr_grid%Pc(il) / (amr_grid%vol_leaf(il) * d2 * par%nphotons)
              end if
           end if
        end do
     case (2)
        do iz = 1, amr_grid%nz_JPa
        do ir = 1, amr_grid%nr_JPa
           if (amr_grid%vol_cyl(ir,iz) > 0.0_wp) &
              amr_grid%Pc2(ir,iz) = amr_grid%Pc2(ir,iz) / (amr_grid%vol_cyl(ir,iz) * d2 * par%nphotons)
        end do
        end do
     case (-1)
        do iz = 1, amr_grid%nz_JPa
           if (amr_grid%vol_plane(iz) > 0.0_wp) &
              amr_grid%Pc1(iz) = amr_grid%Pc1(iz) * (amr_grid%xrange*amr_grid%yrange) / &
                 (amr_grid%vol_plane(iz) * par%nphotons)
        end do
     case default   ! geometry_JPa == 1 (radial / sphere)
        do ir = 1, amr_grid%nr_JPa
           if (amr_grid%vol_sph(ir) > 0.0_wp) &
              amr_grid%Pc1(ir) = amr_grid%Pc1(ir) / (amr_grid%vol_sph(ir) * d2 * par%nphotons)
        end do
     end select
  endif
#endif
#ifdef CALCPnew
  select case (amr_grid%geometry_JPa)
  case (3)
     do il = 1, amr_grid%nleaf
        if (amr_grid%vol_leaf(il) > 0.0_wp) then
           if (par%xy_periodic) then
              amr_grid%Pa_new(il) = amr_grid%Pa_new(il) * (amr_grid%xrange*amr_grid%yrange) / &
                 (amr_grid%vol_leaf(il) * par%nphotons)
           else
              amr_grid%Pa_new(il) = amr_grid%Pa_new(il) / (amr_grid%vol_leaf(il) * d2 * par%nphotons)
           end if
        end if
     end do
  case (2)
     do iz = 1, amr_grid%nz_JPa
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_cyl(ir,iz) > 0.0_wp) &
           amr_grid%P2_new(ir,iz) = amr_grid%P2_new(ir,iz) / (amr_grid%vol_cyl(ir,iz) * d2 * par%nphotons)
     end do
     end do
  case (-1)
     do iz = 1, amr_grid%nz_JPa
        if (amr_grid%vol_plane(iz) > 0.0_wp) &
           amr_grid%P1_new(iz) = amr_grid%P1_new(iz) * (amr_grid%xrange*amr_grid%yrange) / &
              (amr_grid%vol_plane(iz) * par%nphotons)
     end do
  case default   ! geometry_JPa == 1 (radial / sphere)
     do ir = 1, amr_grid%nr_JPa
        if (amr_grid%vol_sph(ir) > 0.0_wp) &
           amr_grid%P1_new(ir) = amr_grid%P1_new(ir) / (amr_grid%vol_sph(ir) * d2 * par%nphotons)
     end do
  end select
#endif

  ! Peel-off observer array normalizations
  if (par%save_peeloff_2D) then
     do k = 1, par%nobs
        scale_factor = par%no_photons * observer(k)%steradian_pix * par%distance2cm**2
        observer(k)%scatt_2D(:,:) = observer(k)%scatt_2D(:,:) / scale_factor
        observer(k)%direc_2D(:,:) = observer(k)%direc_2D(:,:) / scale_factor
        if (par%save_direc0) &
            observer(k)%direc0_2D(:,:) = observer(k)%direc0_2D(:,:) / scale_factor
        if (par%use_stokes) then
           observer(k)%I_2D(:,:) = observer(k)%I_2D(:,:) / scale_factor
           observer(k)%Q_2D(:,:) = observer(k)%Q_2D(:,:) / scale_factor
           observer(k)%U_2D(:,:) = observer(k)%U_2D(:,:) / scale_factor
           observer(k)%V_2D(:,:) = observer(k)%V_2D(:,:) / scale_factor
        end if
     end do
  end if
  if (par%save_peeloff_3D) then
     do k = 1, par%nobs
        scale_factor = par%no_photons * observer(k)%steradian_pix * intensity_bin_unit * par%distance2cm**2
        observer(k)%scatt(:,:,:) = observer(k)%scatt(:,:,:) / scale_factor
        observer(k)%direc(:,:,:) = observer(k)%direc(:,:,:) / scale_factor
        !--- ly_beta band-2 (H-alpha) peel cube: same convention, band-2 bin unit.
        if (associated(observer(k)%peel_Ha)) then
           observer(k)%peel_Ha(:,:,:) = observer(k)%peel_Ha(:,:,:) / &
              (par%no_photons * observer(k)%steradian_pix * intensity_bin_unit_Ha * par%distance2cm**2)
        endif
        if (par%save_direc0) &
            observer(k)%direc0(:,:,:) = observer(k)%direc0(:,:,:) / scale_factor
        if (par%save_dust_scattered) &
            observer(k)%scatt_dust(:,:,:) = observer(k)%scatt_dust(:,:,:) / scale_factor
        if (par%use_stokes) then
           observer(k)%I(:,:,:) = observer(k)%I(:,:,:) / scale_factor
           observer(k)%Q(:,:,:) = observer(k)%Q(:,:,:) / scale_factor
           observer(k)%U(:,:,:) = observer(k)%U(:,:,:) / scale_factor
           observer(k)%V(:,:,:) = observer(k)%V(:,:,:) / scale_factor
        end if
     end do
  end if

  !--- Ly-beta fluorescence: per-band weight bookkeeping (per source photon).
  !--- Conservation: band1_esc + band1_abs + conv = 1 ; band2_esc + band2_abs = conv.
  if (line%line_type == 8 .and. mpar%p_rank == 0) then
     write(*,'(a)')        '--- ly_beta bookkeeping (weights per source photon) ---'
     write(*,'(a,es14.6)') '  band-1 (Ly-beta) escaped        : ', par%W_esc1/dble(par%nphotons)
     write(*,'(a,es14.6)') '  band-1 (Ly-beta) dust-absorbed  : ', par%W_abs1/dble(par%nphotons)
     write(*,'(a,es14.6)') '  conversions (3p->2s)            : ', par%W_conv/dble(par%nphotons)
     write(*,'(a,es14.6)') '  band-2 (H-alpha) escaped        : ', par%W_esc2/dble(par%nphotons)
     write(*,'(a,es14.6)') '  band-2 (H-alpha) dust-absorbed  : ', par%W_abs2/dble(par%nphotons)
     write(*,'(a,es14.6)') '  band1_esc + band1_abs + conv (=1)     : ', &
                           (par%W_esc1 + par%W_abs1 + par%W_conv)/dble(par%nphotons)
     write(*,'(a,es14.6)') '  band2_esc + band2_abs (= conv)        : ', &
                           (par%W_esc2 + par%W_abs2)/dble(par%nphotons)
     write(*,'(a)')        '--------------------------------------------------------'
  endif
  end subroutine output_normalize_amr
!--------------------------------------------------------------
end module output_sum_rect
