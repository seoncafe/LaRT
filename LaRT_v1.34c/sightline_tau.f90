!-- Modification History
!   2021-06-11, now the routines for sightline_tau are gathered in a separate module.
!--
module sightline_tau_mod
  use define
  use utility
  use fitsio_mod
  use memory_mod
contains
  !--------------------------------------------------
  subroutine make_sightline_tau(grid)
  !--- calculate the optical depth along sight lines that are projected to a detector-plane pixel (2020/09/20).
  use mpi
  use write_mod
  implicit none
  type(grid_type),  intent(in) :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist
  real(kind=wp) :: tau_HI,N_HI,tau_dust
  real(kind=wp) :: kx,ky,kz,kr,u1
  integer       :: ix,iy
  integer       :: loop,loop1,loop2,nsize
  integer       :: ierr
  integer       :: i,jj,kk,j0
  character(len=4) :: filename_end

  !--- memories allocation
  do i=1,par%nobs
     call create_mem(observer(i)%tau_HI, [par%nxfreq,par%nxim,par%nyim])
     call create_mem(observer(i)%N_HI,   [par%nxim,par%nyim])
     if (par%DGR > 0.0_wp) then
        call create_mem(observer(i)%tau_dust, [par%nxim,par%nyim])
     endif
  enddo

  !--- calculate tau, N(HI) maps.
  do i=1,par%nobs
    nsize = observer(i)%nxim * observer(i)%nyim
    call loop_divide(nsize,mpar%nproc,mpar%p_rank,loop1,loop2)
    do loop=loop1,loop2
       call array_2D_indices(observer(i)%nxim,observer(i)%nyim,loop,ix,iy)
       kx = tan((ix - (observer(i)%nxim+1.0_wp)/2.0_wp) * observer(i)%dxim/rad2deg)
       ky = tan((iy - (observer(i)%nyim+1.0_wp)/2.0_wp) * observer(i)%dyim/rad2deg)
       kz = -1.0_wp
       kr = sqrt(kx*kx + ky*ky + kz*kz)
       kx = kx/kr
       ky = ky/kr
       kz = kz/kr
       pobs%kx = observer(i)%rmatrix(1,1)*kx + observer(i)%rmatrix(2,1)*ky + observer(i)%rmatrix(3,1)*kz
       pobs%ky = observer(i)%rmatrix(1,2)*kx + observer(i)%rmatrix(2,2)*ky + observer(i)%rmatrix(3,2)*kz
       pobs%kz = observer(i)%rmatrix(1,3)*kx + observer(i)%rmatrix(2,3)*ky + observer(i)%rmatrix(3,3)*kz
       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest
          delt(2) = hugest
       else
          delt(1) = (grid%xmax-observer(i)%x)/pobs%kx
          delt(2) = (grid%xmin-observer(i)%x)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest
          delt(4) = hugest
       else
          delt(3) = (grid%ymax-observer(i)%y)/pobs%ky
          delt(4) = (grid%ymin-observer(i)%y)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest
          delt(6) = hugest
       else
          delt(5) = (grid%zmax-observer(i)%z)/pobs%kz
          delt(6) = (grid%zmin-observer(i)%z)/pobs%kz
       endif

       !-- Find the farthest boundary where the ray touches the grid system.
       !-- We measure optical depth from the distant universe toward Earth. (2020.10.20)
       dist = -999.9_wp
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             pobs%x     = observer(i)%x + pobs%kx * delt(jj)
             pobs%y     = observer(i)%y + pobs%ky * delt(jj)
             pobs%z     = observer(i)%z + pobs%kz * delt(jj)
             pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
             pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
             pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
             !--- to reduce numerical errors (2021.08.05)
             if (jj == 1) pobs%icell = grid%nx+1
             if (jj == 2) pobs%icell = 1
             if (jj == 3) pobs%jcell = grid%ny+1
             if (jj == 4) pobs%jcell = 1
             if (jj == 5) pobs%kcell = grid%nz+1
             if (jj == 6) pobs%kcell = 1
             if (pobs%icell >=1 .and. pobs%icell <= grid%nx+1 .and. &
                 pobs%jcell >=1 .and. pobs%jcell <= grid%ny+1 .and. &
                 pobs%kcell >=1 .and. pobs%kcell <= grid%nz+1) then
                if (delt(jj) > dist) then
                   dist = delt(jj)
                   j0   = jj
                endif
             endif
          endif
       enddo

       if (dist > 0.0_wp .and. dist < hugest) then
          pobs%x     = observer(i)%x + pobs%kx * dist
          pobs%y     = observer(i)%y + pobs%ky * dist
          pobs%z     = observer(i)%z + pobs%kz * dist
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1

          !--- After finding the starting position of the ray, the direction vector should be revsersed. (2020.10.20)
          pobs%kx = -pobs%kx
          pobs%ky = -pobs%ky
          pobs%kz = -pobs%kz

          !--- to reduce numerical errors (2021.08.05)
          if (j0 == 1) then
             pobs%icell = grid%nx+1
             pobs%x     = grid%xface(grid%nx+1)
          else if (j0 == 2) then
             pobs%icell = 1
             pobs%x     = grid%xface(1)
          else if (j0 == 3) then
             pobs%jcell = grid%ny+1
             pobs%y     = grid%yface(grid%ny+1)
          else if (j0 == 4) then
             pobs%jcell = 1
             pobs%y     = grid%yface(1)
          else if (j0 == 5) then
             pobs%kcell = grid%nz+1
             pobs%z     = grid%zface(grid%nz+1)
          else if (j0 == 6) then
             pobs%kcell = 1
             pobs%z     = grid%zface(1)
          endif

          if (pobs%icell == grid%nx+1 .and. pobs%kx < 0.0_wp) pobs%icell = grid%nx
          if (pobs%jcell == grid%ny+1 .and. pobs%ky < 0.0_wp) pobs%jcell = grid%ny
          if (pobs%kcell == grid%nz+1 .and. pobs%kz < 0.0_wp) pobs%kcell = grid%nz

          if (pobs%icell >=1 .and. pobs%icell <= grid%nx .and. &
              pobs%jcell >=1 .and. pobs%jcell <= grid%ny .and. &
              pobs%kcell >=1 .and. pobs%kcell <= grid%nz) then
             !--- Calculate the optical depth for the whole range of frequency (2020.10.19).
             do kk=1, observer(i)%nxfreq
                pobs%xfreq = grid%xfreq(kk)
                u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
                     grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
                     grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
                !--- xfreq should be expressed in the local cell frame (2021.08.05).
                pobs%xfreq = pobs%xfreq * grid%Dfreq_ref/grid%Dfreq(pobs%icell,pobs%jcell,pobs%kcell) - u1

                !--- note that we are using raytrace_to_edge_tau_HI
                call raytrace_to_edge_tau_HI(pobs,grid,tau_HI)
                observer(i)%tau_HI(kk,ix,iy) = tau_HI
             enddo
             call raytrace_to_edge_column(pobs,grid,N_HI,tau_dust)
             observer(i)%N_HI(ix,iy)  = N_HI
             if (par%DGR > 0.0_wp) observer(i)%tau_dust(ix,iy) = tau_dust
          endif
       endif
    enddo

    call reduce_mem(observer(i)%tau_HI)
    call reduce_mem(observer(i)%N_HI)
    if (par%DGR > 0.0_wp) call reduce_mem(observer(i)%tau_dust)
  enddo

  if (mpar%p_rank == 0) then
     !--- write tau and N(HI) maps
     do i = 1, par%nobs
        if (par%nobs == 1) then
           filename_end = ''
        else
           write(filename_end,'(a,i3.3)') '_',i
        endif
        call write_sightline_tau(trim(par%out_file),grid,observer(i), suffix=trim(filename_end))
     enddo
  endif

  !--- memories deallocation
  do i=1, par%nobs
     if (associated(observer(i)%tau_HI))   deallocate(observer(i)%tau_HI)
     if (associated(observer(i)%tau_dust)) deallocate(observer(i)%tau_dust)
     if (associated(observer(i)%N_HI))     deallocate(observer(i)%N_HI)
  enddo
  end subroutine make_sightline_tau
  !-------------------------------------------------------
  subroutine write_sightline_tau(filename,grid,obs,suffix)
  implicit none
  character(len=*),    intent(in) :: filename
  type(grid_type),     intent(in) :: grid
  type(observer_type), intent(in) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  integer            :: unit,status=0
  character(len=128) :: filename1, filename_end
  real(real64)       :: cd1_1, cd1_2, cd2_1, cd2_2
  real(real64)       :: crpix1, crpix2, crval1, crval2
  integer            :: equinox = 2000

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_tau'//trim(filename_end)//'.fits.gz'

  !--- open FITS file.
  call fits_open_new(unit,trim(filename1),status)

  !--- write Images for optical depths
  call fits_append_image(unit,obs%tau_HI,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','TAU_HI','HI optical depth',status)

  !--- write keywords
  cd1_1  = par%dxim
  cd1_2  = 0.0_wp
  cd2_1  = 0.0_wp
  cd2_2  = par%dyim
  crpix1 = (par%nxim+1)/2.0_wp
  crpix2 = (par%nyim+1)/2.0_wp
  crval1 = 0.0_wp
  crval2 = 0.0_wp
  call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
  call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
  call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
  call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
  call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
  call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
  call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
  call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)

  !--- write Images for optical depths
  call fits_append_image(unit,obs%N_HI,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','N_HI','HI column density',status)

  !--- write Images for optical depths
  if (par%DGR > 0.0_wp) then
     call fits_append_image(unit,obs%tau_dust,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','TAU_dust','dust column density',status)
  endif

  !--- close fits file
  call fits_close(unit,status)
  end subroutine write_sightline_tau
  !-------------------------------------------------------
end module sightline_tau_mod
