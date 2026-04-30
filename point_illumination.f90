module point_illumination_mod
  !--
  !-- 2022.06.09: Written by K.I. Seon.
  !--
  use define
  use octree_mod, only: amr_grid, amr_find_leaf
  implicit none
  public random_point_illumination
  public peeling_direct_point_illumination
  public peeling_direct_point_illumination_amr
  private
contains
  !================================================
  subroutine random_point_illumination(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cost, sint, phi, cosp, sinp
  real(kind=wp) :: alpha, beta, dist
  logical,       save :: parameters_initialized = .false.
  real(kind=wp), save :: flux_fac1, costm, dist_wall
  !$OMP THREADPRIVATE(parameters_initialized, flux_fac1, costm, dist_wall)

  !--- This routine works only for the case where the point source is located on the z-axis.
  !--- This must be updated for more general cases (comments added on 2023.02.01).
  if (.not. parameters_initialized) then
     !-- In this method, flux_fac1 = solid angle of the plate subtended by the point source.
     dist_wall = abs(par%zs_point) - grid%zmax
     !-- bug-fixed, 2023.02.01 (see https://vixra.org/pdf/2001.0603v1.pdf)
     alpha     = grid%xmax/dist_wall
     beta      = grid%ymax/dist_wall
     !-- flux_fac1 = solid angle/(4 pi) (2023.02.01)
     !solid angle = 4.0_wp * atan(alpha*beta/sqrt(1.0_wp + alpha**2 + beta**2))
     flux_fac1 = atan(alpha*beta/sqrt(1.0_wp + alpha**2 + beta**2))/pi
     costm     = dist_wall/sqrt(dist_wall**2 + grid%xmax**2 + grid%ymax**2)
     parameters_initialized = .true.
  endif

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp
  do while(.true.)
     !-- photon direction vector.
     cost = rand_number() * (1.0_wp - costm) + costm
     sint = sqrt(1.0_wp - cost**2)
     phi  = twopi * rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     !-- Find the location where the ray meets the medium (grid system).
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     dist      = dist_wall/cost
     photon%x  = dist * photon%kx
     photon%y  = dist * photon%ky

     if (photon%x >= grid%xmin .and. photon%x <= grid%xmax .and. &
         photon%y >= grid%ymin .and. photon%y <= grid%ymax) then
        if (par%zs_point < 0.0_wp) then
           photon%z  = grid%zmin
           photon%kz = cost
        else
           photon%z  = grid%zmax
           photon%kz = -cost
        endif
        if (par%use_amr_grid) then
           photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
           photon%icell = 1; photon%jcell = 1; photon%kcell = 1
        else
           photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
           photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
           if (par%zs_point < 0.0_wp) then
              photon%kcell = 1
           else
              photon%kcell = grid%nz
           endif
        endif
        photon%flux_factor = flux_fac1 * photon%wgt
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo

  if (par%use_stokes) then
     cost = photon%kz
     sint = sqrt(1.0_wp - cost**2)
     if (sint > 0.0_wp) then
        cosp = photon%kx / sint
        sinp = photon%ky / sint
     else
        cosp = 1.0_wp
        sinp = 0.0_wp
     endif

     !--- Set the polarization basis vectors (m) and (n) perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
  endif
  end subroutine random_point_illumination
  !================================================
  subroutine peeling_direct_point_illumination(photon,grid)
  use define
  use random
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: rr,dist,delt(6)
  real(kind=wp) :: r2,r,wgt,tau
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  integer :: ix,iy,ixf
  integer :: i,jj,j0

  !-- take "frequency" from the input photon.
  !-- But, note that the input frequency is expressed in the local cell of the medium.
  !-- Thus, transform the frequency back to the lab frame.
  pobs      = photon
  u1        = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
              grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
              grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
  xfreq_ref = (photon%xfreq + u1) * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)

  !-- reset photon%wgt.
  pobs%wgt = 1.0_wp

  !-- Location in spectral bins. (note that external radiation is emitted in a non-comoving frame).
  ixf      = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

  do i=1,par%nobs
    pobs%kx = observer(i)%x
    pobs%ky = observer(i)%y
    pobs%kz = observer(i)%z - par%zs_point
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !-- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !-- Calculate direc0 (spectral) image.
       if (par%save_direc0) then
          !--- 2D image
          if (par%save_peeloff_2D) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt
          endif

          !--- 3D spectral image
          if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt
          endif
       endif

       !pobs%kx = observer(i)%rmatrix(1,1)*kx + observer(i)%rmatrix(2,1)*ky + observer(i)%rmatrix(3,1)*kz
       !pobs%ky = observer(i)%rmatrix(1,2)*kx + observer(i)%rmatrix(2,2)*ky + observer(i)%rmatrix(3,2)*kz
       !pobs%kz = observer(i)%rmatrix(1,3)*kx + observer(i)%rmatrix(2,3)*ky + observer(i)%rmatrix(3,3)*kz
       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest
          delt(2) = hugest
       else
          delt(1) = grid%xmax/pobs%kx
          delt(2) = grid%xmin/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest
          delt(4) = hugest
       else
          delt(3) = grid%ymax/pobs%ky
          delt(4) = grid%ymin/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest
          delt(6) = hugest
       else
          delt(5) = (grid%zmax-par%zs_point)/pobs%kz
          delt(6) = (grid%zmin-par%zs_point)/pobs%kz
       endif

       !-- Find the closest boundary where the ray touches the grid system.
       dist = hugest
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             pobs%x     = pobs%kx * delt(jj)
             pobs%y     = pobs%ky * delt(jj)
             pobs%z     = pobs%kz * delt(jj) + par%zs_point
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
                if (delt(jj) < dist) then
                   dist = delt(jj)
                   j0   = jj
                endif
             endif
          endif
       enddo

       if (dist < hugest) then
          pobs%x     = pobs%kx * dist
          pobs%y     = pobs%ky * dist
          pobs%z     = pobs%kz * dist + par%zs_point
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1

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

          !-- xfreq_ref is expressed in the lab frame. Here, we need to transform the frequency to the fluid rest frame.
          !-- note that pobs%icell /= photon%icell, etc.
          u2 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
          pobs%xfreq = xfreq_ref * grid%Dfreq_ref/grid%Dfreq(pobs%icell,pobs%jcell,pobs%kcell) - u2

          call raytrace_to_edge(pobs,grid,tau)
          wgt = exp(-tau)/(fourpi*r2) * pobs%wgt
       else
          wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
       endif

       !-- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I_2D(ix,iy) = observer(i)%I_2D(ix,iy) + wgt
          endif
       endif

       !-- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
          endif
       endif
    endif
  enddo
  end subroutine peeling_direct_point_illumination
  !================================================
  subroutine peeling_direct_point_illumination_amr(photon,grid)
  !-- AMR equivalent of peeling_direct_point_illumination.
  !-- Source: external point at (0,0,par%zs_point); ray peels off through the AMR box.
  use define
  use random
  use raytrace_amr_mod, only: raytrace_to_edge_amr
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  type(photon_type) :: pobs
  real(kind=wp) :: dist, delt(6)
  real(kind=wp) :: r2,r,wgt,tau
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  integer :: ix,iy,ixf
  integer :: i,jj,j0,il,il_e

  il = photon%icell_amr

  pobs      = photon
  u1        = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
  xfreq_ref = (photon%xfreq + u1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
  pobs%wgt  = 1.0_wp

  ixf       = floor((xfreq_ref - amr_grid%xfreq_min)/amr_grid%dxfreq)+1

  do i=1,par%nobs
    pobs%kx = observer(i)%x
    pobs%ky = observer(i)%y
    pobs%kz = observer(i)%z - par%zs_point
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       if (par%save_direc0) then
          if (par%save_peeloff_2D) then
             wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt
          endif
          if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= amr_grid%nxfreq) then
             wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt
          endif
       endif

       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest;  delt(2) = hugest
       else
          delt(1) = amr_grid%xmax/pobs%kx
          delt(2) = amr_grid%xmin/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest;  delt(4) = hugest
       else
          delt(3) = amr_grid%ymax/pobs%ky
          delt(4) = amr_grid%ymin/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest;  delt(6) = hugest
       else
          delt(5) = (amr_grid%zmax-par%zs_point)/pobs%kz
          delt(6) = (amr_grid%zmin-par%zs_point)/pobs%kz
       endif

       dist = hugest
       j0   = 0
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             pobs%x = pobs%kx * delt(jj)
             pobs%y = pobs%ky * delt(jj)
             pobs%z = pobs%kz * delt(jj) + par%zs_point
             if (pobs%x >= amr_grid%xmin .and. pobs%x <= amr_grid%xmax .and. &
                 pobs%y >= amr_grid%ymin .and. pobs%y <= amr_grid%ymax .and. &
                 pobs%z >= amr_grid%zmin .and. pobs%z <= amr_grid%zmax) then
                if (delt(jj) < dist) then
                   dist = delt(jj)
                   j0   = jj
                endif
             endif
          endif
       enddo

       if (dist < hugest) then
          pobs%x = pobs%kx * dist
          pobs%y = pobs%ky * dist
          pobs%z = pobs%kz * dist + par%zs_point
          !-- snap to box face exactly to reduce numerical errors.
          if (j0 == 1) pobs%x = amr_grid%xmax
          if (j0 == 2) pobs%x = amr_grid%xmin
          if (j0 == 3) pobs%y = amr_grid%ymax
          if (j0 == 4) pobs%y = amr_grid%ymin
          if (j0 == 5) pobs%z = amr_grid%zmax
          if (j0 == 6) pobs%z = amr_grid%zmin

          il_e = amr_find_leaf(pobs%x, pobs%y, pobs%z)
          if (il_e <= 0) then
             wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
          else
             pobs%icell_amr = il_e
             pobs%icell = 1; pobs%jcell = 1; pobs%kcell = 1
             u2 = amr_grid%vfx(il_e)*pobs%kx + amr_grid%vfy(il_e)*pobs%ky + amr_grid%vfz(il_e)*pobs%kz
             pobs%xfreq = xfreq_ref * amr_grid%Dfreq_ref/amr_grid%Dfreq(il_e) - u2
             call raytrace_to_edge_amr(pobs, grid, tau)
             wgt = exp(-tau)/(fourpi*r2) * pobs%wgt
          endif
       else
          wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
       endif

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I_2D(ix,iy) = observer(i)%I_2D(ix,iy) + wgt
          endif
       endif
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= amr_grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
          endif
       endif
    endif
  enddo
  end subroutine peeling_direct_point_illumination_amr
  !================================================
end module point_illumination_mod
