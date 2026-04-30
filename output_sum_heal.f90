module output_sum_heal
  use define
  use mpi
  use memory_mod
contains
  subroutine output_reduce_inside(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variables
  integer :: k, ierr

  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_dust, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_gas,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nrejected,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  call reduce_mem(grid%Jout)
  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     call reduce_mem(grid%Jabs)
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
        call reduce_mem(observer(k)%scatt_heal_2D)
        call reduce_mem(observer(k)%direc_heal_2D)
        if (par%save_direc0) then
           call reduce_mem(observer(k)%direc0_heal_2D)
        endif
        !if (par%use_stokes) then
        !   call reduce_mem(observer(k)%I_2D)
        !   call reduce_mem(observer(k)%Q_2D)
        !   call reduce_mem(observer(k)%U_2D)
        !   call reduce_mem(observer(k)%V_2D)
        !endif
     enddo
  endif

  !--- 3D spectral image
  if (par%save_peeloff_3D) then
     do k=1,par%nobs
        call reduce_mem(observer(k)%scatt_heal)
        call reduce_mem(observer(k)%direc_heal)
        if (par%save_direc0) then
           call reduce_mem(observer(k)%direc0_heal)
        endif
        !if (par%use_stokes) then
        !   call reduce_mem(observer(k)%I)
        !   call reduce_mem(observer(k)%Q)
        !   call reduce_mem(observer(k)%U)
        !   call reduce_mem(observer(k)%V)
        !endif
     enddo
  endif

  !if (par%save_all_photons) then
  !   !--- Note that allph arrays are declared to be shared memory
  !   !--- so that "reduce" should be performed only in shared_memory = .true. (2020-11-02).
  !   call reduce_mem(allph%rp,          shared_memory = .true.)
  !   call reduce_mem(allph%xfreq1,      shared_memory = .true.)
  !   call reduce_mem(allph%xfreq2,      shared_memory = .true.)
  !   call reduce_mem(allph%nscatt_gas,  shared_memory = .true.)
  !   call reduce_mem(allph%nscatt_dust, shared_memory = .true.)
  !   if (trim(par%source_geometry) /= 'point') then
  !      call reduce_mem(allph%rp0, shared_memory = .true.)
  !   endif
  !   if (par%use_stokes) then
  !      call reduce_mem(allph%I, shared_memory = .true.)
  !      call reduce_mem(allph%Q, shared_memory = .true.)
  !      call reduce_mem(allph%U, shared_memory = .true.)
  !      call reduce_mem(allph%V, shared_memory = .true.)
  !   endif
  !endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end subroutine output_reduce_inside
!-----------------------------------------------------------
  subroutine output_normalize_inside(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variable
  real(kind=wp) :: dVol, area
  real(kind=wp) :: scale_factor
  real(kind=wp) :: intensity_bin_unit
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
     if (associated(grid%Jmu)) &
        grid%Jmu(:,:) = grid%Jmu(:,:) * par%nmu / (par%nphotons*intensity_bin_unit*twopi*2.0_wp)
  else
     ! sphere or box geometry
     ! luminosity is assumed to be 1 photons/whole volume.
     if (par%rmax > 0.0_wp .and. (grid%nx==grid%ny .and. grid%nx==grid%nz)) then
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
!  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
!     if (par%xy_periodic) then
!        ! slab geometry
!        !-grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
!        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
!     else
!        ! sphere or box geometry
!        !-grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*area)
!        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*intensity_bin_unit*twopi*area)
!     endif
!  endif
  if (trim(par%spectral_type) == 'continuum' .and. par%continuum_normalize) then
     if (.not. associated(grid%Jin)) then
        if (mpar%p_rank == 0) write(*,*) &
           'ERROR: continuum_normalize=T requires save_Jin=T'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     endif
     scale_factor = sum(grid%Jin)/size(grid%Jin)
     grid%Jout(:) = grid%Jout(:)/scale_factor
     grid%Jin(:)  = grid%Jin(:) /scale_factor
     if (associated(grid%Jabs))  grid%Jabs(:)  = grid%Jabs(:) /scale_factor
     if (associated(grid%Jabs2)) grid%Jabs2(:) = grid%Jabs2(:)/scale_factor
     if (associated(grid%Jmu))   grid%Jmu(:,:) = grid%Jmu(:,:)/scale_factor
  endif

#ifdef CALCJ
  select case (grid%geometry_JPa)
  case (3)
     if (par%xy_periodic) then
        grid%J(:,:,:,:) = grid%J(:,:,:,:) * (area/(fourpi*dVol*par%nphotons*intensity_bin_unit))
     else
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
        if (grid%ncount_cyl(j,k) > 0) grid%J2(:,j,k) = grid%J2(:,j,k)/grid%ncount_cyl(j,k) / (fourpi*dVol*par%nphotons*intensity_bin_unit)
     enddo
     enddo
  case (1)
     !--- Note that the factors 2 and 8 in the case of xyz_symmetry are considered in add_to_J, add_to_Pnew and add_to_Pa.
     do k=1,grid%nr
        if (grid%ncount_sph(k) > 0) grid%J1(:,k) = grid%J1(:,k)/grid%ncount_sph(k) / (fourpi*dVol*par%nphotons*intensity_bin_unit)
     enddo
  case (-1)
     !grid%J1 = grid%J1/(grid%nx*grid%ny) * (area/(fourpi*dVol*par%nphotons*intensity_bin_unit))
     do k=1,grid%nz
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
        observer(k)%scatt_heal_2D(:)  = observer(k)%scatt_heal_2D(:) / scale_factor
        observer(k)%direc_heal_2D(:)  = observer(k)%direc_heal_2D(:) / scale_factor

        if (par%save_direc0) then
           observer(k)%direc0_heal_2D(:) = observer(k)%direc0_heal_2D(:)/ scale_factor
        endif
        !if (par%use_stokes) then
        !   observer(k)%I_2D(:) = observer(k)%I_2D(:) / scale_factor
        !   observer(k)%Q_2D(:) = observer(k)%Q_2D(:) / scale_factor
        !   observer(k)%U_2D(:) = observer(k)%U_2D(:) / scale_factor
        !   observer(k)%V_2D(:) = observer(k)%V_2D(:) / scale_factor
        !endif
     enddo
  endif
  !--- 3D spectral image
  if (par%save_peeloff_3D) then
     do k=1,par%nobs
        !-scale_factor = par%no_photons*observer(k)%steradian_pix*grid%dxfreq * par%distance2cm**2
        scale_factor = par%no_photons*observer(k)%steradian_pix*intensity_bin_unit * par%distance2cm**2
        observer(k)%scatt_heal(:,:)  = observer(k)%scatt_heal(:,:) / scale_factor
        observer(k)%direc_heal(:,:)  = observer(k)%direc_heal(:,:) / scale_factor

        if (par%save_direc0) then
           observer(k)%direc0_heal(:,:) = observer(k)%direc0_heal(:,:)/ scale_factor
        endif
        !if (par%use_stokes) then
        !   observer(k)%I(:,:) = observer(k)%I(:,:) / scale_factor
        !   observer(k)%Q(:,:) = observer(k)%Q(:,:) / scale_factor
        !   observer(k)%U(:,:) = observer(k)%U(:,:) / scale_factor
        !   observer(k)%V(:,:) = observer(k)%V(:,:) / scale_factor
        !endif
     enddo
  endif

  end subroutine output_normalize_inside
!--------------------------------------------------------------
  subroutine output_reduce_inside_amr(grid)
  !-- AMR + inside-observer reduction.
  !-- Reduces amr_grid%J* (then copies to grid%J*) plus HEALPix observer arrays.
  use octree_mod, only: amr_grid
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: k, ierr

  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nscatt_dust, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nscatt_gas,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, par%nrejected,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jout(1), amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jin)) &
     call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jin(1),  amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jabs)) &
     call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jabs(1), amr_grid%nxfreq, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (allocated(amr_grid%Jmu)) &
     call MPI_ALLREDUCE(MPI_IN_PLACE, amr_grid%Jmu(1,1), amr_grid%nxfreq*size(amr_grid%Jmu,2), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  grid%Jout = amr_grid%Jout
  if (allocated(amr_grid%Jin))  grid%Jin  = amr_grid%Jin
  if (allocated(amr_grid%Jabs)) grid%Jabs = amr_grid%Jabs
  if (allocated(amr_grid%Jmu))  grid%Jmu  = amr_grid%Jmu

  if (par%save_peeloff_2D) then
     do k = 1, par%nobs
        call reduce_mem(observer(k)%scatt_heal_2D)
        call reduce_mem(observer(k)%direc_heal_2D)
        if (par%save_direc0) call reduce_mem(observer(k)%direc0_heal_2D)
     enddo
  endif
  if (par%save_peeloff_3D) then
     do k = 1, par%nobs
        call reduce_mem(observer(k)%scatt_heal)
        call reduce_mem(observer(k)%direc_heal)
        if (par%save_direc0) call reduce_mem(observer(k)%direc0_heal)
     enddo
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine output_reduce_inside_amr
!--------------------------------------------------------------
  subroutine output_normalize_inside_amr(grid)
  !-- AMR + inside-observer normalization.
  !-- Uses the AMR-convention area = fourpi * distance2cm**2 (matches output_normalize_amr),
  !-- so Jout/Jin are independent of L_box. HEALPix observer arrays are
  !-- normalized identically to the Cartesian _inside path.
  implicit none
  type(grid_type), intent(inout) :: grid
  real(wp) :: area, scale_factor, intensity_bin_unit
  integer  :: k

  par%nscatt_dust = par%nscatt_dust / par%nphotons
  par%nscatt_gas  = par%nscatt_gas  / par%nphotons
  par%nscatt_tot  = par%nscatt_gas  + par%nscatt_dust

  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dwave
  else
     intensity_bin_unit = grid%dxfreq
  endif

  area = fourpi * par%distance2cm**2

  grid%Jout(:) = grid%Jout(:) / (par%nphotons * intensity_bin_unit * twopi * area)
  if (associated(grid%Jin)) &
     grid%Jin(:)  = grid%Jin(:)  / (par%nphotons * intensity_bin_unit * twopi * area)
  if (associated(grid%Jabs)) &
     grid%Jabs(:) = grid%Jabs(:) / (par%nphotons * intensity_bin_unit * twopi * area)
  if (associated(grid%Jmu)) &
     grid%Jmu(:,:) = grid%Jmu(:,:) * par%nmu / (par%nphotons * intensity_bin_unit * twopi * area)

  if (trim(par%spectral_type) == 'continuum' .and. par%continuum_normalize) then
     if (.not. associated(grid%Jin)) then
        if (mpar%p_rank == 0) write(*,*) 'ERROR: continuum_normalize=T requires save_Jin=T'
        call MPI_ABORT(MPI_COMM_WORLD, 1, k)
     endif
     scale_factor = sum(grid%Jin)/size(grid%Jin)
     grid%Jout(:) = grid%Jout(:)/scale_factor
     grid%Jin(:)  = grid%Jin(:) /scale_factor
     if (associated(grid%Jabs)) grid%Jabs(:) = grid%Jabs(:)/scale_factor
     if (associated(grid%Jmu))  grid%Jmu(:,:) = grid%Jmu(:,:)/scale_factor
  endif

  if (par%save_peeloff_2D) then
     do k = 1, par%nobs
        scale_factor = par%no_photons*observer(k)%steradian_pix * par%distance2cm**2
        observer(k)%scatt_heal_2D(:) = observer(k)%scatt_heal_2D(:) / scale_factor
        observer(k)%direc_heal_2D(:) = observer(k)%direc_heal_2D(:) / scale_factor
        if (par%save_direc0) &
           observer(k)%direc0_heal_2D(:) = observer(k)%direc0_heal_2D(:) / scale_factor
     enddo
  endif
  if (par%save_peeloff_3D) then
     do k = 1, par%nobs
        scale_factor = par%no_photons*observer(k)%steradian_pix * intensity_bin_unit * par%distance2cm**2
        observer(k)%scatt_heal(:,:) = observer(k)%scatt_heal(:,:) / scale_factor
        observer(k)%direc_heal(:,:) = observer(k)%direc_heal(:,:) / scale_factor
        if (par%save_direc0) &
           observer(k)%direc0_heal(:,:) = observer(k)%direc0_heal(:,:) / scale_factor
     enddo
  endif
  end subroutine output_normalize_inside_amr
!--------------------------------------------------------------
end module output_sum_heal
