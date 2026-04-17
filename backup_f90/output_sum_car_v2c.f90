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
  integer :: nobs_reduce

  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_dust, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,par%nscatt_HI,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

#ifdef CALCJ
  call reduce_mem(grid%J)
#endif
#ifdef CALCP
  call reduce_mem(grid%Pa)
#endif
#ifdef CALCPnew
  call reduce_mem(grid%Pa_new)
#endif
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
#ifdef PEELINGOFF
  if (par%nobs > 1 .and. par%peel_average) then
     do k=2,par%nobs
        observer(1)%scatt(:,:,:)  = observer(1)%scatt(:,:,:)  + observer(k)%scatt(:,:,:)
        observer(1)%direc(:,:,:)  = observer(1)%direc(:,:,:)  + observer(k)%direc(:,:,:)
        if (par%save_direc0) then
           observer(1)%direc0(:,:,:) = observer(1)%direc0(:,:,:) + observer(k)%direc0(:,:,:)
        endif
        if (par%use_stokes) then
           observer(1)%I(:,:,:) = observer(1)%I(:,:,:) + observer(k)%I(:,:,:)
           observer(1)%Q(:,:,:) = observer(1)%Q(:,:,:) + observer(k)%Q(:,:,:)
           observer(1)%U(:,:,:) = observer(1)%U(:,:,:) + observer(k)%U(:,:,:)
           observer(1)%V(:,:,:) = observer(1)%V(:,:,:) + observer(k)%V(:,:,:)
        endif
        if (par%save_sightline_tau) then
           observer(1)%tau_HI(:,:,:) = observer(1)%tau_HI(:,:,:) + observer(k)%tau_HI(:,:,:)
           observer(1)%N_HI(:,:)     = observer(1)%N_HI(:,:)     + observer(k)%N_HI(:,:)
           if (par%DGR > 0.0_wp) observer(1)%tau_dust(:,:) = observer(1)%tau_dust(:,:) + observer(k)%tau_dust(:,:)
        endif
     enddo
     observer(1)%scatt(:,:,:)  = observer(1)%scatt(:,:,:)  / par%nobs
     observer(1)%direc(:,:,:)  = observer(1)%direc(:,:,:)  / par%nobs
     if (par%save_direc0) then
        observer(1)%direc0(:,:,:) = observer(1)%direc0(:,:,:) / par%nobs
     endif
     if (par%use_stokes) then
        observer(1)%I(:,:,:)   = observer(1)%I(:,:,:) / par%nobs
        observer(1)%Q(:,:,:)   = observer(1)%Q(:,:,:) / par%nobs
        observer(1)%U(:,:,:)   = observer(1)%U(:,:,:) / par%nobs
        observer(1)%V(:,:,:)   = observer(1)%V(:,:,:) / par%nobs
        if (observer(1)%nxim == observer(1)%nyim) call make_radial_stokes(grid,observer(1))
     endif
     if (observer(1)%nxim == observer(1)%nyim) call make_radial_spectrum(observer(1))
     if (par%save_sightline_tau) then
        observer(1)%tau_HI(:,:,:) = observer(1)%tau_HI(:,:,:) / par%nobs
        observer(1)%N_HI(:,:)     = observer(1)%N_HI(:,:)     / par%nobs
        if (par%DGR > 0.0_wp) observer(1)%tau_dust(:,:) = observer(1)%tau_dust(:,:) / par%nobs
     endif
     nobs_reduce = 1
  else
     nobs_reduce = par%nobs
  endif
  do k=1,nobs_reduce
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
#endif
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
  integer       :: k, nobs_reduce

  !--- slightly modified, 2020.09.16
  par%nscatt_dust = par%nscatt_dust/par%nphotons
  par%nscatt_HI   = par%nscatt_HI  /par%nphotons
  par%nscatt_tot  = par%nscatt_HI + par%nscatt_dust

  ! Convert rhokap into density unit * distance2cm. (number/cm^3) * distance2cm.
  ! 2017-06-27
  do k=1,grid%nz
     grid%rhokap(:,:,k) = grid%rhokap(:,:,k)*grid%Dfreq(:,:,k)/par%cross0
  enddo

  if (par%xy_periodic) then
     ! slab geometry
     !-- bug-fixed, 2020.10.23, (distance2cm)**2 was missing.
     area            = grid%xrange*grid%yrange * (par%distance2cm)**2
     dVol            = grid%dx*grid%dy*grid%dz * (par%distance2cm)**2
     ! solid angle = 2pi : all photons are propagating outward only. no photons get scattered back through the surface.
     ! 2.0 is for two (top and bottom) boundaries.
     ! luminosity is assumed to be 1 photons/cm^2.
     ! J(nu,x,y,z) is essentially a sum of path lengths, so that distance2cm should be muliplied to convert to cgs unit.
     grid%Jout(:)    = grid%Jout(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
     grid%Jin(:)     = grid%Jin(:) /(par%nphotons*grid%dxfreq*twopi*2.0_wp)
#ifdef CALCJ
     grid%J(:,:,:,:) = grid%J(:,:,:,:) * (area/(fourpi*dVol*par%nphotons*grid%dxfreq))
#endif
#ifdef CALCP
     ! note: macbook pro 2012 (24 cores) caused segmentation error when whole 3D array is used in where statement (2017-06-13).
     do k=1,grid%nz
        where (grid%rhokap(:,:,k) > 0.0_wp)
           grid%Pa(:,:,k) = (grid%Pa(:,:,k)/grid%rhokap(:,:,k))*(area/(dVol*par%nphotons))
        endwhere
     enddo
#endif
#ifdef CALCPnew
     ! note: macbook pro 2012 (24 cores) caused segmentation error when whole 3D array is used in where statement (2017-06-13).
     do k=1,grid%nz
        where (grid%rhokap(:,:,k) > 0.0_wp)
           grid%Pa_new(:,:,k) = (grid%Pa_new(:,:,k)/grid%rhokap(:,:,k))*(area/(dVol*par%nphotons))
        endwhere
     enddo
#endif
  else
     ! sphere or box geometry
     ! luminosity is assumed to be 1 photons/whole volume.
     if (par%rmax > 0.0_wp .and. (grid%nx==grid%ny .and. grid%nx==grid%nz)) then
        area = fourpi * grid%rmax**2 * (par%distance2cm)**2
     else
        area = (grid%xmax*grid%ymax + grid%ymax*grid%zmax + grid%zmax*grid%xmax)*8.0_wp * (par%distance2cm)**2
     endif
     dVol            = grid%dx*grid%dy*grid%dz * (par%distance2cm)**2
     grid%Jout(:)    = grid%Jout(:)/(par%nphotons*grid%dxfreq*twopi*area)
     grid%Jin(:)     = grid%Jin(:) /(par%nphotons*grid%dxfreq*twopi*area)
#ifdef CALCJ
     grid%J(:,:,:,:) = grid%J(:,:,:,:) / (fourpi*dVol * par%nphotons*grid%dxfreq)
#endif
#ifdef CALCP
     do k=1,grid%nz
        where (grid%rhokap(:,:,k) > 0.0_wp)
           grid%Pa(:,:,k) = (grid%Pa(:,:,k)/grid%rhokap(:,:,k))/(dVol*par%nphotons)
        endwhere
     enddo
#endif
#ifdef CALCPnew
     do k=1,grid%nz
        where (grid%rhokap(:,:,k) > 0.0_wp)
           grid%Pa_new(:,:,k) = (grid%Pa_new(:,:,k)/grid%rhokap(:,:,k))/(dVol*par%nphotons)
        endwhere
     enddo
#endif
  endif

  if (par%xyz_symmetry) then
#ifdef CALCP
     grid%Pa(:,:,:)  = grid%Pa(:,:,:) /8.0_wp
     if ((par%nx/2)*2 /= par%nx) then
        grid%Pa(1,:,:)  = grid%Pa(1,:,:) *2.0_wp
     endif
     if ((par%ny/2)*2 /= par%ny) then
        grid%Pa(:,1,:)  = grid%Pa(:,1,:) *2.0_wp
     endif
     if ((par%nz/2)*2 /= par%nz) then
        grid%Pa(:,:,1)  = grid%Pa(:,:,1) *2.0_wp
     endif
#endif
#ifdef CALCPnew
  if (par%xyz_symmetry) then
     grid%Pa_new(:,:,:)  = grid%Pa_new(:,:,:) /8.0_wp
     if ((par%nx/2)*2 /= par%nx) then
        grid%Pa_new(1,:,:)  = grid%Pa_new(1,:,:) *2.0_wp
     endif
     if ((par%ny/2)*2 /= par%ny) then
        grid%Pa_new(:,1,:)  = grid%Pa_new(:,1,:) *2.0_wp
     endif
     if ((par%nz/2)*2 /= par%nz) then
        grid%Pa_new(:,:,1)  = grid%Pa_new(:,:,1) *2.0_wp
     endif
  endif
#endif
#ifdef CALCJ
     grid%J(:,:,:,:) = grid%J(:,:,:,:)/8.0_wp
     if ((par%nx/2)*2 /= par%nx) then
        grid%J(:,1,:,:) = grid%J(:,1,:,:)*2.0_wp
     endif
     if ((par%ny/2)*2 /= par%ny) then
        grid%J(:,:,1,:) = grid%J(:,:,1,:)*2.0_wp
     endif
     if ((par%nz/2)*2 /= par%nz) then
        grid%J(:,:,:,1) = grid%J(:,:,:,1)*2.0_wp
     endif
#endif
  endif

  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     if (par%xy_periodic) then
        ! slab geometry
        grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
     else
        ! sphere or box geometry
        grid%Jabs(:) = grid%Jabs(:)/(par%nphotons*grid%dxfreq*twopi*area)
     endif
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     if (par%xy_periodic) then
        ! slab geometry
        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*2.0_wp)
     else
        ! sphere or box geometry
        grid%Jabs2(:) = grid%Jabs2(:)/(par%nphotons*grid%dxfreq*twopi*area)
     endif
  endif

  if (grid%nx == grid%ny .and. grid%nx == grid%nz) call make_radial_profile(grid)
  if (par%xy_periodic) call make_z_profile(grid)

#ifdef PEELINGOFF
  if (par%peel_average) then
     nobs_reduce = 1
  else
     nobs_reduce = par%nobs
  endif
  do k=1,nobs_reduce
     scale_factor = par%no_photons*observer(k)%steradian_pix*grid%dxfreq * par%distance2cm**2
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
        if (observer(k)%nxim == observer(k)%nyim) call make_radial_stokes(grid,observer(k))
     endif
     if (observer(k)%nxim == observer(k)%nyim) call make_radial_spectrum(observer(k))
  enddo
#endif

  end subroutine output_normalize
!--------------------------------------------------------------
#ifdef PEELINGOFF
  subroutine make_radial_spectrum(obs)
  implicit none
  type(observer_type), intent(inout) :: obs
  integer       :: i,j,ir,nr
  real(kind=wp) :: xcen,ycen
  real(kind=wp) :: xx,yy,rr
  integer, allocatable :: ncount(:)

  nr = (maxval([obs%nxim,obs%nyim])+1)/2
  if (.not. associated(obs%radial_spec)) allocate(obs%radial_spec(obs%nxfreq,nr))
  if (.not. allocated(ncount))           allocate(ncount(nr))
  obs%radial_spec(:,:) = 0.0_wp
  ncount(:)            = 0

  xcen = (obs%nxim+1.0_wp)/2.0_wp
  ycen = (obs%nxim+1.0_wp)/2.0_wp
  do j=1,obs%nyim
     yy = j - ycen
     do i=1,obs%nxim
        xx = i - xcen
        rr = sqrt(xx**2 + yy**2)
        !ir = floor(rr) + 1
        if ((nr/2)*2 == nr) then
           ir = floor(rr) + 1
        else
           ir = floor(rr+0.5_wp) + 1
        endif
        if (ir >= 1 .and. ir <= nr) then
           ncount(ir) = ncount(ir) + 1
           obs%radial_spec(:,ir) = obs%radial_spec(:,ir) + obs%scatt(:,i,j) + obs%direc(:,i,j)
        endif
     enddo
  enddo

  do ir=1,nr
     if (ncount(ir) > 0) obs%radial_spec(:,ir) = obs%radial_spec(:,ir)/ncount(ir)
  enddo
  if (allocated(ncount)) deallocate(ncount)
  end subroutine make_radial_spectrum
!--------------------------------------------------------------
  subroutine make_radial_stokes(grid,obs)
  implicit none
  type(grid_type),     intent(in)    :: grid
  type(observer_type), intent(inout) :: obs
  integer              :: i,j,ir,nr
  real(kind=wp)        :: xcen,ycen
  real(kind=wp)        :: xx,yy,rr
  !real(kind=wp)        :: sumI
  real(kind=wp)        :: cosp,sinp,cos2p,sin2p
  integer, allocatable :: ncount(:)

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
           obs%radial_I(ir) = obs%radial_I(ir) + sum(obs%I(:,i,j))
           obs%radial_Q(ir) = obs%radial_Q(ir) + sum( obs%Q(:,i,j)*cos2p + obs%U(:,i,j)*sin2p)
           obs%radial_U(ir) = obs%radial_U(ir) + sum(-obs%Q(:,i,j)*sin2p + obs%U(:,i,j)*cos2p)
           obs%radial_V(ir) = obs%radial_V(ir) + sum(obs%V(:,i,j))
        endif
     enddo
  enddo

  do ir=1,nr
     if (ncount(ir) > 0) then
        !--- integral over the frequency range (2020.09.30).
        obs%radial_I(ir)   = obs%radial_I(ir) / ncount(ir) * grid%dxfreq
        obs%radial_Q(ir)   = obs%radial_Q(ir) / ncount(ir) * grid%dxfreq
        obs%radial_U(ir)   = obs%radial_U(ir) / ncount(ir) * grid%dxfreq
        obs%radial_V(ir)   = obs%radial_V(ir) / ncount(ir) * grid%dxfreq
        !--- ncount is a number of pixels and hence radial_I = 0 can happens (2020.11.08).
        if (obs%radial_I(ir) > 0.0_wp) then
           obs%radial_pol(ir) = sqrt(obs%radial_Q(ir)**2 + obs%radial_U(ir)**2)/obs%radial_I(ir)
        endif
     endif
  enddo
  if (allocated(ncount)) deallocate(ncount)
  end subroutine make_radial_stokes
#endif
!--------------------------------------------------------------
  subroutine make_radial_profile(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  integer       :: i,j,k,ir,nr
  real(kind=wp) :: xx,yy,zz,rr,dr,rmin,rmax
  integer, allocatable :: ncount(:)

  nr   = maxval([grid%nx,grid%ny,grid%nz])
  rmax = minval([grid%xmax,grid%ymax,grid%zmax])
  if ((nr/2)*2 == nr) then
     dr   = rmax/nr
     rmin = 0.0_wp
  else
     dr   = rmax/(nr-0.5_wp)
     rmin = -dr/2.0_wp
  endif

#ifdef CALCJ
  if (.not. associated(grid%J1)) allocate(grid%J1(grid%nxfreq,nr))
  grid%J1(:,:) = 0.0_wp
#endif
#ifdef CALCP
  if (.not. associated(grid%P1)) allocate(grid%P1(nr))
  grid%P1(:)   = 0.0_wp
#endif
#ifdef CALCPnew
  if (.not. associated(grid%P1_new)) allocate(grid%P1_new(nr))
  grid%P1_new(:)   = 0.0_wp
#endif
  if (.not. allocated(ncount))   allocate(ncount(nr))
  ncount(:)    = 0

  do k=1,grid%nz
    zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
    do j=1,grid%ny
      yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
      do i=1,grid%nx
        xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
        rr = sqrt(xx**2 + yy**2 + zz**2)
        ir = floor((rr-rmin)/dr) + 1
        if (grid%rhokap(i,j,k) > 0.0_wp .and. ir >= 1 .and. ir <= nr) then
           ncount(ir)    = ncount(ir)    + 1
#ifdef CALCP
           grid%P1(ir)   = grid%P1(ir)   + grid%Pa(i,j,k)
#endif
#ifdef CALCPnew
           grid%P1_new(ir)   = grid%P1_new(ir)   + grid%Pa_new(i,j,k)
#endif
#ifdef CALCJ
           grid%J1(:,ir) = grid%J1(:,ir) + grid%J(:,i,j,k)
#endif
        endif
      enddo
    enddo
  enddo

#ifdef CALCP
  where(ncount > 0)
     grid%P1(:) = grid%P1(:)/ncount(:)
  endwhere
#endif
#ifdef CALCPnew
  where(ncount > 0)
     grid%P1_new(:) = grid%P1_new(:)/ncount(:)
  endwhere
#endif
#ifdef CALCJ
  do k=1,nr
     if (ncount(k) > 0) grid%J1(:,k) = grid%J1(:,k)/ncount(k)
  enddo
#endif
  if (allocated(ncount)) deallocate(ncount)

  end subroutine make_radial_profile
!--------------------------------------------------------------
  subroutine make_z_profile(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  integer :: k

#ifdef CALCP
  if (.not. associated(grid%P1)) allocate(grid%P1(grid%nz))
  do k=1,grid%nz
     grid%P1(k)   = sum(grid%Pa(:,:,k))
  enddo
  grid%P1 = grid%P1/(grid%nx*grid%ny)
#endif
#ifdef CALCPnew
  if (.not. associated(grid%P1_new)) allocate(grid%P1_new(grid%nz))
  do k=1,grid%nz
     grid%P1_new(k)   = sum(grid%Pa_new(:,:,k))
  enddo
  grid%P1_new = grid%P1_new/(grid%nx*grid%ny)
#endif

#ifdef CALCJ
  if (.not. associated(grid%J1)) allocate(grid%J1(grid%nxfreq,grid%nz))
  do k=1,grid%nz
     grid%J1(:,k) = sum(sum(grid%J(:,:,:,k),3),2)
  enddo
  grid%J1 = grid%J1/(grid%nx*grid%ny)
#endif

  end subroutine make_z_profile
!--------------------------------------------------------------
end module output_sum
