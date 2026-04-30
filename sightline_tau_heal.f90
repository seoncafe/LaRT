!-- Modification History
!   2021-06-11, now the routines for sightline_tau are gathered in a separate module.
!--
module sightline_tau_heal
  use define
  use utility
  use fitsio_mod
  use memory_mod
  use octree_mod, only: amr_grid, amr_find_leaf
contains
  !--------------------------------------------------
  subroutine make_sightline_tau_inside(grid)
  !--- calculate the optical depth along sight lines that are projected to a detector-plane pixel (2020/09/20).
  !use mpi
  use healpix
  !--- now, write subrotines are defined in define.f90 (2023.01.16).
  !use write_mod
  implicit none
  type(grid_type),  intent(in) :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist
  real(kind=wp) :: tau_gas,N_gas,tau_dust
  real(kind=wp) :: kx,ky,kz,kr,u1
  integer       :: ipix
  integer       :: loop,loop1,loop2
  integer       :: ierr
  integer       :: i,jj,kk,j0
  character(len=4) :: filename_end

  !--- memories allocation
  do i=1,par%nobs
     call create_shared_mem(observer(i)%tau_gas_heal, [par%nxfreq,par%npix])
     call create_shared_mem(observer(i)%N_gas_heal,   [par%npix])
     if (par%DGR > 0.0_wp) then
        call create_shared_mem(observer(i)%tau_dust_heal, [par%npix])
     endif
  enddo

  !--- calculate tau, N(gas) maps.
  do i=1,par%nobs
    call loop_divide(par%npix,mpar%nproc,mpar%p_rank,loop1,loop2)
    do ipix=loop1,loop2
       call pix2vec(observer(i)%nside, ipix, pobs%kx,pobs%ky,pobs%kz)
       !--- (pobs%kx, pobs%ky, pobs%kz) is the direction from the observer to the boundary.
       if (pobs%kx == 0.0_wp) then
          delt(1) = -999.0
          delt(2) = -999.0
       else
          delt(1) = (grid%xmax-observer(i)%x)/pobs%kx
          delt(2) = (grid%xmin-observer(i)%x)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = -999.0
          delt(4) = -999.0
       else
          delt(3) = (grid%ymax-observer(i)%y)/pobs%ky
          delt(4) = (grid%ymin-observer(i)%y)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = -999.0
          delt(6) = -999.0
       else
          delt(5) = (grid%zmax-observer(i)%z)/pobs%kz
          delt(6) = (grid%zmin-observer(i)%z)/pobs%kz
       endif

       !-- Find the closest boundary where the ray touches the grid system.
       !-- We measure optical depth from the distant universe toward Earth. (2020.10.20)
       dist = hugest
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             if (delt(jj) < dist) then
                dist = delt(jj)
                j0   = jj
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

          if (pobs%icell >=1 .and. pobs%icell <= grid%nx+1 .and. &
              pobs%jcell >=1 .and. pobs%jcell <= grid%ny+1 .and. &
              pobs%kcell >=1 .and. pobs%kcell <= grid%nz+1) then
             !--- Calculate the optical depth for the whole range of frequency (2020.10.19).
             do kk=1, observer(i)%nxfreq
                pobs%xfreq = grid%xfreq(kk)
                u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
                     grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
                     grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
                !--- xfreq should be expressed in the local cell frame (2021.08.05).
                pobs%xfreq = pobs%xfreq * grid%Dfreq_ref/grid%Dfreq(pobs%icell,pobs%jcell,pobs%kcell) - u1

                !--- note that we are using raytrace_to_dist_tau_gas
                call raytrace_to_dist_tau_gas(pobs,grid,dist,tau_gas)
                observer(i)%tau_gas_heal(kk,ipix) = tau_gas
             enddo
             call raytrace_to_dist_column(pobs,grid,dist,N_gas,tau_dust)
             observer(i)%N_gas_heal(ipix)  = N_gas
             if (par%DGR > 0.0_wp) observer(i)%tau_dust_heal(ipix) = tau_dust
          endif
       endif
       !print*,ipix,N_gas
    enddo

    call reduce_mem(observer(i)%tau_gas_heal, shared_memory=.true.)
    call reduce_mem(observer(i)%N_gas_heal,   shared_memory=.true.)
    if (par%DGR > 0.0_wp) call reduce_mem(observer(i)%tau_dust_heal, shared_memory=.true.)
  enddo

  if (mpar%p_rank == 0) then
     !--- write tau and N(gas) maps
     do i = 1, par%nobs
        if (par%nobs == 1) then
           filename_end = ''
        else
           write(filename_end,'(a,i3.3)') '_',i
        endif
        call write_sightline_tau_inside(trim(par%out_file),grid,observer(i), suffix=trim(filename_end))
     enddo
  endif

  !--- memories deallocation
  do i=1, par%nobs
     call destroy_mem(observer(i)%tau_gas_heal)
     call destroy_mem(observer(i)%N_gas_heal)
     if (par%DGR > 0.0_wp) call destroy_mem(observer(i)%tau_dust_heal)
  enddo
  end subroutine make_sightline_tau_inside
  !-------------------------------------------------------
  subroutine make_sightline_tau_inside_amr(grid)
  !--- AMR equivalent of make_sightline_tau_inside.
  !--- Cell lookup uses amr_find_leaf; raytrace_to_dist_* are procedure
  !--- pointers already routed to AMR variants in setup_procedure().
  use healpix
  implicit none
  type(grid_type),  intent(in) :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist
  real(kind=wp) :: tau_gas,N_gas,tau_dust
  real(kind=wp) :: u1
  integer       :: ipix
  integer       :: loop1,loop2
  integer       :: ierr
  integer       :: i,jj,kk,j0,il
  character(len=4) :: filename_end

  do i=1,par%nobs
     call create_shared_mem(observer(i)%tau_gas_heal, [par%nxfreq,par%npix])
     call create_shared_mem(observer(i)%N_gas_heal,   [par%npix])
     if (par%DGR > 0.0_wp) then
        call create_shared_mem(observer(i)%tau_dust_heal, [par%npix])
     endif
  enddo

  do i=1,par%nobs
    call loop_divide(par%npix,mpar%nproc,mpar%p_rank,loop1,loop2)
    do ipix=loop1,loop2
       call pix2vec(observer(i)%nside, ipix, pobs%kx,pobs%ky,pobs%kz)
       if (pobs%kx == 0.0_wp) then
          delt(1) = -999.0; delt(2) = -999.0
       else
          delt(1) = (grid%xmax-observer(i)%x)/pobs%kx
          delt(2) = (grid%xmin-observer(i)%x)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = -999.0; delt(4) = -999.0
       else
          delt(3) = (grid%ymax-observer(i)%y)/pobs%ky
          delt(4) = (grid%ymin-observer(i)%y)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = -999.0; delt(6) = -999.0
       else
          delt(5) = (grid%zmax-observer(i)%z)/pobs%kz
          delt(6) = (grid%zmin-observer(i)%z)/pobs%kz
       endif

       dist = hugest
       j0   = 0
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             if (delt(jj) < dist) then
                dist = delt(jj)
                j0   = jj
             endif
          endif
       enddo

       if (dist > 0.0_wp .and. dist < hugest) then
          pobs%x = observer(i)%x + pobs%kx * dist
          pobs%y = observer(i)%y + pobs%ky * dist
          pobs%z = observer(i)%z + pobs%kz * dist

          !--- reverse direction so the ray points from box surface toward observer.
          pobs%kx = -pobs%kx
          pobs%ky = -pobs%ky
          pobs%kz = -pobs%kz

          !--- nudge slightly inside the box face we just hit, then look up the leaf.
          if      (j0 == 1) then
             pobs%x = grid%xmax - tiny(1.0_wp)*max(1.0_wp,abs(grid%xmax))
          else if (j0 == 2) then
             pobs%x = grid%xmin + tiny(1.0_wp)*max(1.0_wp,abs(grid%xmin))
          else if (j0 == 3) then
             pobs%y = grid%ymax - tiny(1.0_wp)*max(1.0_wp,abs(grid%ymax))
          else if (j0 == 4) then
             pobs%y = grid%ymin + tiny(1.0_wp)*max(1.0_wp,abs(grid%ymin))
          else if (j0 == 5) then
             pobs%z = grid%zmax - tiny(1.0_wp)*max(1.0_wp,abs(grid%zmax))
          else if (j0 == 6) then
             pobs%z = grid%zmin + tiny(1.0_wp)*max(1.0_wp,abs(grid%zmin))
          endif
          il = amr_find_leaf(pobs%x, pobs%y, pobs%z)
          pobs%icell_amr = il
          pobs%icell = 1; pobs%jcell = 1; pobs%kcell = 1

          if (il > 0) then
             do kk=1, observer(i)%nxfreq
                pobs%xfreq = grid%xfreq(kk)
                u1 = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
                pobs%xfreq = pobs%xfreq * grid%Dfreq_ref/amr_grid%Dfreq(il) - u1

                call raytrace_to_dist_tau_gas(pobs,grid,dist,tau_gas)
                observer(i)%tau_gas_heal(kk,ipix) = tau_gas
             enddo
             call raytrace_to_dist_column(pobs,grid,dist,N_gas,tau_dust)
             observer(i)%N_gas_heal(ipix) = N_gas
             if (par%DGR > 0.0_wp) observer(i)%tau_dust_heal(ipix) = tau_dust
          endif
       endif
    enddo

    call reduce_mem(observer(i)%tau_gas_heal, shared_memory=.true.)
    call reduce_mem(observer(i)%N_gas_heal,   shared_memory=.true.)
    if (par%DGR > 0.0_wp) call reduce_mem(observer(i)%tau_dust_heal, shared_memory=.true.)
  enddo

  if (mpar%p_rank == 0) then
     do i = 1, par%nobs
        if (par%nobs == 1) then
           filename_end = ''
        else
           write(filename_end,'(a,i3.3)') '_',i
        endif
        call write_sightline_tau_inside(trim(par%out_file),grid,observer(i), suffix=trim(filename_end))
     enddo
  endif

  do i=1, par%nobs
     call destroy_mem(observer(i)%tau_gas_heal)
     call destroy_mem(observer(i)%N_gas_heal)
     if (par%DGR > 0.0_wp) call destroy_mem(observer(i)%tau_dust_heal)
  enddo
  end subroutine make_sightline_tau_inside_amr
  !-------------------------------------------------------
  subroutine write_sightline_tau_inside(filename,grid,obs,suffix)
  implicit none
  character(len=*),    intent(in) :: filename
  type(grid_type),     intent(in) :: grid
  type(observer_type), intent(in) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  integer            :: unit,status=0
  character(len=128) :: filename1, filename_end
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

  !--- write Images for gas optical depth
  call fits_append_image(unit,obs%tau_gas_heal,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','TAU_gas','gas optical depth',status)

  !--- write keywords
  call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
  call fits_put_keyword(unit,'NSIDE',  obs%nside, 'HEALPIX nside',status)
  call fits_put_keyword(unit,'NPIX',   obs%npix,  'HEALPIX npix',status)
  call fits_put_keyword(unit,'OBSX',   obs%x,     'observer x',status)
  call fits_put_keyword(unit,'OBSY',   obs%y,     'observer y',status)
  call fits_put_keyword(unit,'OBSZ',   obs%z,     'observer z',status)
  call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)

  !--- write Images for gas column density
  call fits_append_image(unit,obs%N_gas_heal,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','N_gas','gas column density',status)

  !--- write Images for dust optical depth
  if (par%DGR > 0.0_wp) then
     call fits_append_image(unit,obs%tau_dust_heal,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','TAU_dust','dust optical depth',status)
  endif

  !--- close fits file
  call fits_close(unit,status)
  end subroutine write_sightline_tau_inside
  !-------------------------------------------------------
end module sightline_tau_heal
