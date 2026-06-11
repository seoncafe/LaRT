module observer_heal
  use define
  use utility
  use random
  use memory_mod
  implicit none
contains
  !---------------------------------------------------------------------------------
  !-----------------
  subroutine observer_create_inside()
  use mpi
  use healpix
  implicit none
  real(kind=wp) :: cost, sint, phi
  integer       :: ierr, i

  if (par%observer_located_inside) then
     if (.not.(is_finite(par%obsx(1)) .and. is_finite(par%obsy(1)) .and. is_finite(par%obsz(1)))) then
         par%obsx(1) = 0.0
         par%obsy(1) = 0.0
         par%obsz(1) = 0.0
     endif
  endif

  if (is_finite(par%obsx(1)) .and. is_finite(par%obsy(1)) .and. is_finite(par%obsz(1))) then
     par%nobs = count(is_finite(par%obsx) .and. is_finite(par%obsy) .and. is_finite(par%obsz))
     if (.not.allocated(observer)) allocate(observer(par%nobs))

     par%npix = nside2npix(par%nside)
     do i=1, par%nobs
        !-- (par%obsx, par%obsy, par%obsz) are the actual coordinates of observer.
        observer(i)%x     = par%obsx(i)
        observer(i)%y     = par%obsy(i)
        observer(i)%z     = par%obsz(i)
        observer(i)%nside = par%nside
        observer(i)%npix  = par%npix
     enddo
  else
     par%nobs            = 0
     par%save_peeloff    = .false.
     par%save_peeloff_2D = .false.
     par%save_peeloff_3D = .false.
  endif

  do i=1, par%nobs
     observer(i)%steradian_pix = 4.0d0*pi/observer(i)%npix
     observer(i)%nxfreq        = par%nxfreq

     ! memory allocations of output 2D images
     if (par%save_peeloff_2D) then
        call create_mem(observer(i)%scatt_heal_2D, [par%npix])
        call create_mem(observer(i)%direc_heal_2D, [par%npix])
        if (par%save_direc0) then
           call create_mem(observer(i)%direc0_heal_2D,[par%npix])
        endif
        !--- To be done. We need to add this part for Stokes!
        !if (par%use_stokes) then
        !endif
     endif

     ! memory allocations of output 3D spectral images
     if (par%save_peeloff_3D) then
        call create_mem(observer(i)%scatt_heal, [par%nxfreq,par%npix])
        call create_mem(observer(i)%direc_heal, [par%nxfreq,par%npix])
        if (par%save_direc0) then
           call create_mem(observer(i)%direc0_heal,[par%nxfreq,par%npix])
        endif
        !if (par%use_stokes) then
        !endif
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine observer_create_inside
  !-----------------
  subroutine observer_destroy_inside()
  use define
  use mpi
  implicit none
  integer :: ierr, i

  call destroy_shared_mem_all()

  do i=1, par%nobs
     if (associated(observer(i)%scatt_heal))  deallocate(observer(i)%scatt_heal)
     if (associated(observer(i)%direc_heal))  deallocate(observer(i)%direc_heal)
     if (associated(observer(i)%direc0_heal)) deallocate(observer(i)%direc0_heal)

     !if (associated(observer(i)%I)) deallocate(observer(i)%I)
     !if (associated(observer(i)%Q)) deallocate(observer(i)%Q)
     !if (associated(observer(i)%U)) deallocate(observer(i)%U)
     !if (associated(observer(i)%V)) deallocate(observer(i)%V)

     if (associated(observer(i)%scatt_heal_2D))  deallocate(observer(i)%scatt_heal_2D)
     if (associated(observer(i)%direc_heal_2D))  deallocate(observer(i)%direc_heal_2D)
     if (associated(observer(i)%direc0_heal_2D)) deallocate(observer(i)%direc0_heal_2D)

     !if (associated(observer(i)%I_2D)) deallocate(observer(i)%I_2D)
     !if (associated(observer(i)%Q_2D)) deallocate(observer(i)%Q_2D)
     !if (associated(observer(i)%U_2D)) deallocate(observer(i)%U_2D)
     !if (associated(observer(i)%V_2D)) deallocate(observer(i)%V_2D)
  enddo
  end subroutine observer_destroy_inside
  !----------------------------------------
end module observer_heal
