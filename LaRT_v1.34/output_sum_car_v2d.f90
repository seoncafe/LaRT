module output_sum
  use define
  use mpi
  use memory_mod
contains
  subroutine output_reduce(grid)
  implicit none
  type(grid_type), intent(inout) :: grid
  !--- local variables
  integer :: k, ierr

  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_dust, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_HI,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nrejected,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call MPI_ALLREDUCE(MPI_IN_PLACE,par%flux_factor, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     par%flux_factor = par%flux_factor/(par%no_photons + par%nrejected)
  endif
  par%acceptance_rate = par%no_photons/(par%no_photons + par%nrejected)

  call reduce_mem(grid%Jout)
  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     call reduce_mem(grid%Jabs)
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     call reduce_mem(grid%Jabs2)
  endif
  if (par%save_Jin) then
     call reduce_mem(grid%Jin)
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
        if (par%save_direc0) then
           call reduce_mem(observer(k)%direc0)
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
     call reduce_mem(allph%nscatt_HI,   shared_memory = .true.)
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

  end subroutine output_reduce
!-----------------------------------------------------------
  subroutine output_normalize(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variable
  real(kind=wp) :: dVol, area
  real(kind=wp) :: scale_factor
  real(kind=wp) :: intensity_bin_unit
  integer       :: i,j,k

  !--- slightly modified, 2020.09.16
  par%nscatt_dust = par%nscatt_dust/par%nphotons
  par%nscatt_HI   = par%nscatt_HI  /par%nphotons
  par%nscatt_tot  = par%nscatt_HI + par%nscatt_dust

  !--- intensity unit.
  if (par%intensity_unit == 1) then
     intensity_bin_unit = grid%dlambda
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
     grid%Jin(:)     = grid%Jin(:) /(par%nphotons*intensity_bin_unit*twopi*2.0_wp)
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
     grid%Jin(:)     = grid%Jin(:) /(par%nphotons*intensity_bin_unit*twopi*area)
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
  if (trim(par%spectral_type) == 'continuum' .and. par%continuum_normalize) then
     scale_factor = sum(grid%Jin)/size(grid%Jin)
     grid%Jout(:) = grid%Jout(:)/scale_factor
     grid%Jin(:)  = grid%Jin(:) /scale_factor
     if (associated(grid%Jabs))  grid%Jabs(:)  = grid%Jabs(:) /scale_factor
     if (associated(grid%Jabs2)) grid%Jabs2(:) = grid%Jabs2(:)/scale_factor
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

        if (par%save_direc0) then
           observer(k)%direc0(:,:,:) = observer(k)%direc0(:,:,:)/ scale_factor
        endif
        if (par%use_stokes) then
           observer(k)%I(:,:,:) = observer(k)%I(:,:,:) / scale_factor
           observer(k)%Q(:,:,:) = observer(k)%Q(:,:,:) / scale_factor
           observer(k)%U(:,:,:) = observer(k)%U(:,:,:) / scale_factor
           observer(k)%V(:,:,:) = observer(k)%V(:,:,:) / scale_factor
        endif
     enddo
  endif

  end subroutine output_normalize
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
     intensity_bin_unit = grid%dlambda
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
     intensity_bin_unit = grid%dlambda
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
end module output_sum
