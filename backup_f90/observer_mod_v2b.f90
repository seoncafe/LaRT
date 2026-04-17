module observer_mod
  use define
  use utility
  use random
  use memory_mod
  implicit none
contains
!---------------------------------------------------------------------------------
#ifdef PEELINGOFF
  !-----------------
  subroutine observer_create()
  use mpi
  implicit none
  real(kind=wp), allocatable :: cosa(:),cosb(:),cosg(:)
  real(kind=wp), allocatable :: sina(:),sinb(:),sing(:)
  real(kind=wp) :: dist_scale
  real(kind=wp) :: cost, sint, phi
  integer       :: ierr, i
  !-- 8 vertices of a cube
  real(kind=wp), parameter :: vertex_x(8) = [1.0,  1.0,  1.0, -1.0, -1.0, -1.0,  1.0,  -1.0]
  real(kind=wp), parameter :: vertex_y(8) = [1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  -1.0]
  real(kind=wp), parameter :: vertex_z(8) = [1.0, -1.0,  1.0,  1.0,  1.0, -1.0, -1.0,  -1.0]
  real(kind=wp) :: px, py, pz, kx, ky, kz, ang_x, ang_y, max_ang_x, max_ang_y
  integer       :: iv

  !--- setup observer and image plane
  ! alpha = -phase angle
  ! beta  = -inclination angle
  ! gamma = -posision angle
  !         Position angle is measured counter-clock wise from x-axis of the detector plane.
  ! Rmaxtrix is the matrix for a rotation starting from the grid frame (initial fixed frame) to observer frame.
  ! rotation sequence is (1) alpha about z-axis,
  !                      (2) beta  about new y-axis, and
  !                      (3) gamma about new z-axis.
  !---
  if (any(is_finite(par%phase_angle)))       par%alpha = -par%phase_angle
  if (any(is_finite(par%inclination_angle))) par%beta  = -par%inclination_angle
  if (any(is_finite(par%position_angle)))    par%gamma = -par%position_angle

  !--- Default setting: when nothing is given, only a single observer located at a coordinate (0,0,100)*box size.
  if (.not.  (is_finite(par%alpha(1)) .and. is_finite(par%beta(1)) .and. is_finite(par%gamma(1))) .and. &
      .not.  (is_finite(par%obsx(1))  .and. is_finite(par%obsy(1)) .and. is_finite(par%obsz(1)))   ) then
     if (.not. is_finite(par%distance)) then
        par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
     endif
     par%obsx(1)  = 0.0
     par%obsy(1)  = 0.0
     par%obsz(1)  = par%distance

     par%alpha(1) = 0.0
     par%beta(1)  = 0.0
     par%gamma(1) = 0.0
  endif

  if (is_finite(par%obsx(1)) .and. is_finite(par%obsy(1)) .and. is_finite(par%obsz(1))) then
     par%nobs = count(is_finite(par%obsx) .and. is_finite(par%obsy) .and. is_finite(par%obsz))
     if (.not.allocated(cosa))     allocate(cosa(par%nobs))
     if (.not.allocated(cosb))     allocate(cosb(par%nobs))
     if (.not.allocated(cosg))     allocate(cosg(par%nobs))
     if (.not.allocated(sina))     allocate(sina(par%nobs))
     if (.not.allocated(sinb))     allocate(sinb(par%nobs))
     if (.not.allocated(sing))     allocate(sing(par%nobs))
     if (.not.allocated(observer)) allocate(observer(par%nobs))


     if (.not. is_finite(par%distance)) then
        par%distance = sqrt(par%obsx(1)**2 + par%obsy(1)**2 + par%obsz(1)**2)
        if (par%distance < 10.0_wp*maxval([par%xmax, par%ymax, par%zmax])) then
           par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
        endif
     endif

     do i=1, par%nobs
        if (.not. is_finite(par%gamma(i))) par%gamma(i) = 0.0_wp

        ! observer's coordinates represented in the grid system.
        dist_scale    = par%distance/sqrt(par%obsx(i)**2 + par%obsy(i)**2 + par%obsz(i)**2)
        observer(i)%x = par%obsx(i) * dist_scale
        observer(i)%y = par%obsy(i) * dist_scale
        observer(i)%z = par%obsz(i) * dist_scale

        cosb(i) = observer(i)%z/par%distance
        sinb(i) = sqrt(1.0_wp - cosb(i)**2)
        if (sinb(i) == 0.0_wp) then
           cosa(i) = 1.0_wp
           sina(i) = 0.0_wp
        else
           cosa(i) = observer(i)%x/(par%distance*sinb(i))
           sina(i) = observer(i)%y/(par%distance*sinb(i))
        endif
        cosg(i) = cos(par%gamma(i)*deg2rad)
        sing(i) = sin(par%gamma(i)*deg2rad)
     enddo
  else
     par%nobs = count(is_finite(par%alpha) .and. is_finite(par%beta) .and. is_finite(par%gamma))
     if (.not.allocated(cosa))     allocate(cosa(par%nobs))
     if (.not.allocated(cosb))     allocate(cosb(par%nobs))
     if (.not.allocated(cosg))     allocate(cosg(par%nobs))
     if (.not.allocated(sina))     allocate(sina(par%nobs))
     if (.not.allocated(sinb))     allocate(sinb(par%nobs))
     if (.not.allocated(sing))     allocate(sing(par%nobs))
     if (.not.allocated(observer)) allocate(observer(par%nobs))

     if (.not. is_finite(par%distance)) then
        par%distance = maxval([par%xmax, par%ymax, par%zmax]) * 100.0_wp
     endif

     !--- Note that the previous convention was to rotate the coordinate system starting from galaxy to observer.
     !--- To following the previous convention, we need to change the sign of the three angles.
     !--- comment added on 2020.08.15
     do i=1, par%nobs
        cosa(i) = cos(par%alpha(i)*deg2rad)
        sina(i) = sin(par%alpha(i)*deg2rad)
        cosb(i) = cos(par%beta(i) *deg2rad)
        sinb(i) = sin(par%beta(i) *deg2rad)
        cosg(i) = cos(par%gamma(i)*deg2rad)
        sing(i) = sin(par%gamma(i)*deg2rad)

        ! observer's coordinates represented in the grid system.
        observer(i)%x = par%distance*cosa(i)*sinb(i)
        observer(i)%y = par%distance*sina(i)*sinb(i)
        observer(i)%z = par%distance*cosb(i)
     enddo
  endif

  !--- rotation matrix from Grid system to Observer system.
  do i=1, par%nobs
     observer(i)%rmatrix(1,1) =  cosa(i)*cosb(i)*cosg(i) - sina(i)*sing(i)
     observer(i)%rmatrix(1,2) =  sina(i)*cosb(i)*cosg(i) + cosa(i)*sing(i)
     observer(i)%rmatrix(1,3) = -sinb(i)*cosg(i)
     observer(i)%rmatrix(2,1) = -cosa(i)*cosb(i)*sing(i) - sina(i)*cosg(i)
     observer(i)%rmatrix(2,2) = -sina(i)*cosb(i)*sing(i) + cosa(i)*cosg(i)
     observer(i)%rmatrix(2,3) =  sinb(i)*sing(i)
     observer(i)%rmatrix(3,1) =  cosa(i)*sinb(i)
     observer(i)%rmatrix(3,2) =  sina(i)*sinb(i)
     observer(i)%rmatrix(3,3) =  cosb(i)
  enddo

  !--- parameters for observer's image plane
  !--- 2020.10.20, determine the image size that convers whole cloud's grid system.
  if (.not. (is_finite(par%dxim) .and. is_finite(par%dyim))) then
     !--- For a spherical geometry, the maximum angular size is always the radius, regardless of the observer direction.
     if (par%rmax > 0.0_wp) then
        par%dxim = atan2(par%rmax,par%distance)/(par%nxim/2.0_wp) * rad2deg
        par%dyim = atan2(par%rmax,par%distance)/(par%nyim/2.0_wp) * rad2deg
     else
        max_ang_x = -999.0
        max_ang_y = -999.0
        do i = 1, par%nobs
        do iv = 1,8
           px = observer(i)%x - vertex_x(iv) * par%xmax
           py = observer(i)%y - vertex_y(iv) * par%ymax
           pz = observer(i)%z - vertex_z(iv) * par%zmax
           kx = observer(i)%rmatrix(1,1)*px + observer(i)%rmatrix(1,2)*py + observer(i)%rmatrix(1,3)*pz
           ky = observer(i)%rmatrix(2,1)*px + observer(i)%rmatrix(2,2)*py + observer(i)%rmatrix(2,3)*pz
           kz = observer(i)%rmatrix(3,1)*px + observer(i)%rmatrix(3,2)*py + observer(i)%rmatrix(3,3)*pz
           ang_x = abs(atan2(-kx, kz))
           ang_y = abs(atan2(-ky, kz))
           if (max_ang_x < ang_x) max_ang_x = ang_x
           if (max_ang_y < ang_y) max_ang_y = ang_y
        enddo
        enddo
        if (par%nxim == par%nyim) then
           par%dxim = maxval([max_ang_x, max_ang_y]) / (par%nxim/2.0_wp) * rad2deg
           par%dyim = maxval([max_ang_x, max_ang_y]) / (par%nyim/2.0_wp) * rad2deg
        else
           if (.not. is_finite(par%dxim)) par%dxim = max_ang_x / (par%nxim/2.0_wp) * rad2deg
           if (.not. is_finite(par%dyim)) par%dyim = max_ang_y / (par%nyim/2.0_wp) * rad2deg
        endif
     endif
     !if (.not. is_finite(par%dxim)) par%dxim = atan2(par%xmax,par%distance)/(par%nxim/2.0_wp) * rad2deg
     !if (.not. is_finite(par%dyim)) par%dyim = atan2(par%ymax,par%distance)/(par%nyim/2.0_wp) * rad2deg
  endif

  do i=1, par%nobs
     observer(i)%dxim = par%dxim
     observer(i)%dyim = par%dyim
     observer(i)%nxim = par%nxim
     observer(i)%nyim = par%nyim

     observer(i)%steradian_pix = par%dxim * par%dyim * (deg2rad**2)
     observer(i)%nxfreq        = par%nxfreq

     ! memory allocations of output images
     call create_mem(observer(i)%scatt, [par%nxfreq,par%nxim,par%nyim])
     call create_mem(observer(i)%direc, [par%nxfreq,par%nxim,par%nyim])
     if (par%save_direc0) then
        call create_mem(observer(i)%direc0,[par%nxfreq,par%nxim,par%nyim])
     endif
     if (par%use_stokes) then
        call create_mem(observer(i)%I, [par%nxfreq,par%nxim,par%nyim])
        call create_mem(observer(i)%Q, [par%nxfreq,par%nxim,par%nyim])
        call create_mem(observer(i)%U, [par%nxfreq,par%nxim,par%nyim])
        call create_mem(observer(i)%V, [par%nxfreq,par%nxim,par%nyim])
     endif
     if (par%save_sightline_tau) then
        call create_mem(observer(i)%tau_HI, [par%nxfreq,par%nxim,par%nyim])
        call create_mem(observer(i)%N_HI,   [par%nxim,par%nyim])
        if (par%DGR > 0.0_wp) then
           call create_mem(observer(i)%tau_dust, [par%nxim,par%nyim])
        endif
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !--- deallocate cosines and sines of alpha,beta,gamma.
  if (allocated(cosa)) deallocate(cosa)
  if (allocated(cosb)) deallocate(cosb)
  if (allocated(cosg)) deallocate(cosg)
  if (allocated(sina)) deallocate(sina)
  if (allocated(sinb)) deallocate(sinb)
  if (allocated(sing)) deallocate(sing)

  end subroutine observer_create
  !-----------------
  subroutine observer_destroy()
  use define
  use mpi
  implicit none
  integer :: ierr, i

  call destroy_shared_mem_all()

  do i=1, par%nobs
     if (associated(observer(i)%tau_HI)) deallocate(observer(i)%tau_HI)
     if (associated(observer(i)%N_HI))   deallocate(observer(i)%N_HI)

     if (associated(observer(i)%scatt))  deallocate(observer(i)%scatt)
     if (associated(observer(i)%direc))  deallocate(observer(i)%direc)
     if (associated(observer(i)%direc0)) deallocate(observer(i)%direc0)

     if (associated(observer(i)%I)) deallocate(observer(i)%I)
     if (associated(observer(i)%Q)) deallocate(observer(i)%Q)
     if (associated(observer(i)%U)) deallocate(observer(i)%U)
     if (associated(observer(i)%V)) deallocate(observer(i)%V)

     if (associated(observer(i)%radial_spec)) deallocate(observer(i)%radial_spec)
     if (associated(observer(i)%radial_pol))  deallocate(observer(i)%radial_pol)
     if (associated(observer(i)%radial_r))    deallocate(observer(i)%radial_r)
     if (associated(observer(i)%radial_I))    deallocate(observer(i)%radial_I)
     if (associated(observer(i)%radial_Q))    deallocate(observer(i)%radial_Q)
     if (associated(observer(i)%radial_U))    deallocate(observer(i)%radial_U)
     if (associated(observer(i)%radial_V))    deallocate(observer(i)%radial_V)
  enddo
  end subroutine observer_destroy
  !----------------------------------------
#endif
end module observer_mod
