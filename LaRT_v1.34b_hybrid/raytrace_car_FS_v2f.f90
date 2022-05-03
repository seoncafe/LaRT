module raytrace
  use voigt_mod
contains

  subroutine raytrace_to_edge_car(photon0,grid,tau)
!--- Find the cumulative optical depth to the edge of grid.
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: tau

! local variables
  !-- note that exp(-tau_huge) = 0.0
  !-- at this large optical depth peeling-off contribution is practically zero.
  !-- Hence, no need to integrate further.
  real(kind=wp), parameter :: tau_huge = 745.2_wp
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap, DnuHK, u1,u2,xfreq
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau        = 0.0_wp
  d          = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz
  xfreq = photon0%xfreq

  !--- integrate through grid
  do while(.true.)
     DnuHK  = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokap = grid%rhokap(icell,jcell,kcell)*(voigt(xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                              voigt(xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) rhokap = rhokap + grid%rhokapD(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau   = tau + (tx - d) * rhokap
        d     = tx
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        tau   = tau + (ty - d) * rhokap
        d     = ty
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        tau   = tau + (tz - d) * rhokap
        d     = tz
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif
     if (tau >= tau_huge) exit

     u2    = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     xfreq = (xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
  enddo
  return
  end subroutine raytrace_to_edge_car

  subroutine raytrace_to_edge_car_tau_HI(photon0,grid,tau_HI)
!--- Find the cumulative optical depth due to HI and dust, column density of HI to the edge of grid. (2020.10.21)
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: tau_HI

! local variables
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap, DnuHK, u1,u2,xfreq
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau_HI  = 0.0_wp
  d       = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz
  xfreq = photon0%xfreq

  !--- integrate through grid
  do while(.true.)
     DnuHK  = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokap = grid%rhokap(icell,jcell,kcell)*(voigt(xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                              voigt(xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau_HI = tau_HI + (tx - d) * rhokap
        d      = tx
        icell  = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        tau_HI = tau_HI + (ty - d) * rhokap
        d      = ty
        jcell  = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        tau_HI = tau_HI + (tz - d) * rhokap
        d      = tz
        kcell  = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif

     u2    = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     xfreq = (xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
  enddo
  return
  end subroutine raytrace_to_edge_car_tau_HI

  subroutine raytrace_to_edge_car_column(photon0,grid,N_HI,tau_dust)
!--- Find the cumulative optical depth due to HI and dust, column density of HI to the edge of grid. (2020.10.21)
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: N_HI, tau_dust

! local variables
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: u1,u2,xfreq
  real(kind=wp) :: rho,rhokapD,del
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  N_HI     = 0.0_wp
  tau_dust = 0.0_wp
  d        = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell

  !--- integrate through grid
  do while(.true.)
     rho    = grid%rhokap(icell,jcell,kcell)*grid%Dfreq(icell,jcell,kcell) / par%cross0
     if (par%DGR > 0.0_wp) rhokapD = grid%rhokapD(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del    = tx - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = tx
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = ty
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        del    = tz - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = tz
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif

     iold = icell
     jold = jcell
     kold = kcell
  enddo
  return
  end subroutine raytrace_to_edge_car_column

  subroutine raytrace_to_tau_car(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           photon%inside = .false.
           return
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

!--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del   = tx - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx = tx + delx
     else if (idx_min == 2) then
        del   = ty - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty = ty + dely
     else
        del   = tz - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d   - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
     !dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car

  subroutine raytrace_to_tau_car_xyzsym(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- 2017-05-15, Reflected at xy-, xz- and yz-planes.

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     istep = -1
     delx  = -grid%dx/kx
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%i0
           istep = 1
           kx    = abs(kx)
        endif
     endif
     if (istep == 1) then
        tx = (grid%xface(icell+1)-xp)/kx
     else
        tx = (grid%xface(icell)-xp)/kx
     endif
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     jstep = -1
     dely  = -grid%dy/ky
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%j0
           jstep = 1
           ky    = abs(ky)
        endif
     endif
     if (jstep == 1) then
        ty = (grid%yface(jcell+1)-yp)/ky
     else
        ty = (grid%yface(jcell)-yp)/ky
     endif
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     kstep = -1
     delz  = -grid%dz/kz
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           kcell = grid%k0
           kstep = 1
           kz    = abs(kz)
        endif
     endif
     if (kstep == 1) then
        tz = (grid%zface(kcell+1)-zp)/kz
     else
        tz = (grid%zface(kcell)-zp)/kz
     endif
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  !tauold = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

!--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del    = tx - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        icell = icell + istep
        if (icell < 1) then
           icell = grid%i0
           istep = 1
           kx    = -kx
        else if (icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx     = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1) then
           jcell = grid%j0
           jstep = 1
           ky    = -ky
        else if (jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty     = ty + dely
     else
        del    = tz - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1) then
           !kcell = 1
           kcell = grid%k0
           kstep = 1
           kz    = -kz
        else if (kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz    = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold   = icell
     jold   = jcell
     kold   = kcell
     u1     = u2
     !dold   = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     !u1 = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! here, photon%kx, etc should be used instead of kx,ky,kz
  photon%x  = photon%x + d*photon%kx
  photon%y  = photon%y + d*photon%ky
  photon%z  = photon%z + d*photon%kz
  if (photon%x < grid%xmin) photon%x = -photon%x
  if (photon%y < grid%ymin) photon%y = -photon%y
  if (photon%z < grid%zmin) photon%z = -photon%z
  photon%kx = kx
  photon%ky = ky
  photon%kz = kz
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xyzsym

  subroutine raytrace_to_tau_car_xysym(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- 2017-05-15, Reflected at xy-, xz- and yz-planes.

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     istep = -1
     delx  = -grid%dx/kx
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%i0
           istep = 1
           kx    = abs(kx)
        endif
     endif
     if (istep == 1) then
        tx = (grid%xface(icell+1)-xp)/kx
     else
        tx = (grid%xface(icell)-xp)/kx
     endif
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     jstep = -1
     dely  = -grid%dy/ky
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%j0
           jstep = 1
           ky    = abs(ky)
        endif
     endif
     if (jstep == 1) then
        ty = (grid%yface(jcell+1)-yp)/ky
     else
        ty = (grid%yface(jcell)-yp)/ky
     endif
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     kstep = -1
     delz  = -grid%dz/kz
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           kcell = grid%k0
           kstep = 1
           kz    = abs(kz)
        endif
     endif
     if (kstep == 1) then
        tz = (grid%zface(kcell+1)-zp)/kz
     else
        tz = (grid%zface(kcell)-zp)/kz
     endif
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  !tauold = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

!--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del    = tx - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        icell = icell + istep
        if (icell < 1) then
           icell = grid%i0
           istep = 1
           kx    = -kx
        else if (icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx     = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1) then
           jcell = grid%j0
           jstep = 1
           ky    = -ky
        else if (jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty     = ty + dely
     else
        del    = tz - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz    = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold   = icell
     jold   = jcell
     kold   = kcell
     u1     = u2
     !dold   = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     !u1 = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! here, photon%kx, etc should be used instead of kx,ky,kz
  photon%x  = photon%x + d*photon%kx
  photon%y  = photon%y + d*photon%ky
  photon%z  = photon%z + d*photon%kz
  if (photon%x < grid%xmin) photon%x = -photon%x
  if (photon%y < grid%ymin) photon%y = -photon%y
  photon%kx = kx
  photon%ky = ky
  photon%kz = kz
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xysym

  subroutine raytrace_to_tau_car_xyper(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Period Boundary Codition around x, y axis for plane-paralle geometry (2016-12-07)

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        icell = 1
        xp    = grid%xface(icell)
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%nx
           xp    = grid%xface(icell+1)
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        jcell = 1
        yp    = grid%yface(jcell)
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%ny
           yp    = grid%yface(jcell+1)
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  d     = 0.0_wp
  tau   = 0.0_wp
  dold  = 0.0_wp
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold, jold, kold) *photon%kx + grid%vfy(iold, jold, kold) *photon%ky + grid%vfz(iold, jold, kold) *photon%kz

  !--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx, ty, tz], dim=1)
     if (idx_min == 1) then
        del    = tx - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1)       icell = grid%nx
        if (icell > grid%nx) icell = 1
        tx = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1)       jcell = grid%ny
        if (jcell > grid%ny) jcell = 1
        ty = ty + dely
     else
        del    = tz - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     ! transform to the rest-frame of new cell
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
     dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  !--- (2016-12-07)
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! fold back into the system (2016-12-11)
  photon%x  = xp - floor((xp-grid%xmin)/grid%xrange)*grid%xrange
  photon%y  = yp - floor((yp-grid%ymin)/grid%yrange)*grid%yrange
  photon%z  = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xyper

  subroutine raytrace_to_tau_car_zonly(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tz,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

!--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     del    = tz - d
     tau    = tau + del * rhokap
     dtauH  = del * rhokapH
     d      = tz
     if (tau >= tau_in) then
        if (rhokap > 0.0) then
           d_overshoot = (tau - tau_in)/rhokap
           d     = d - d_overshoot
           del   = del - d_overshoot
           dtauH = del * rhokapH
        endif
        xp = xp + d * kx
        yp = yp + d * ky
        zp = zp + d * kz
        exit
     endif
     kcell = kcell + kstep
     if (kcell < 1 .or. kcell > grid%nz) then
        photon%inside = .false.
        exit
     endif
     tz = tz + delz

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     kold = kcell
     u1   = u2
     !dold = d
  enddo

  if (.not. photon%inside) then
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_zonly

  subroutine raytrace_to_tau_car_xyper_shear(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Period Boundary Codition around x, y axis for plane-paralle geometry (2016-12-07)
!--- Shear Effect is taken into account (2017-07-14).

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        icell = 1
        xp    = grid%xface(icell)
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%nx
           xp    = grid%xface(icell+1)
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        jcell = 1
        yp    = grid%yface(jcell)
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%ny
           yp    = grid%yface(jcell+1)
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  d     = 0.0_wp
  tau   = 0.0_wp
  dold  = 0.0_wp
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold,jold,kold)*photon%kx + (grid%vfy(iold,jold,kold)+photon%vfy_shear)*photon%ky + &
          grid%vfz(iold,jold,kold)*photon%kz

  !--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx, ty, tz], dim=1)
     if (idx_min == 1) then
        del    = tx - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1) then
           icell            = grid%nx
           ! shear effect
           photon%vfy_shear = photon%vfy_shear - par%Omega
        endif
        if (icell > grid%nx) then
           icell            = 1
           ! shear effect
           photon%vfy_shear = photon%vfy_shear + par%Omega
        endif
        tx = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1)       jcell = grid%ny
        if (jcell > grid%ny) jcell = 1
        ty = ty + dely
     else
        del    = tz - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*photon%kx + (grid%vfy(icell,jcell,kcell)+photon%vfy_shear)*photon%ky + &
          grid%vfz(icell,jcell,kcell)*photon%kz
     ! transform to the rest-frame of new cell
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
     dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  !--- (2017-07-14)
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*photon%kx + (grid%vfy(icell,jcell,kcell)+photon%vfy_shear)*photon%ky + &
          grid%vfz(icell,jcell,kcell)*photon%kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! fold back into the system (2016-12-11)
  photon%x  = xp - floor((xp-grid%xmin)/grid%xrange)*grid%xrange
  photon%y  = yp - floor((yp-grid%ymin)/grid%yrange)*grid%yrange
  photon%z  = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xyper_shear

  subroutine raytrace_to_tau_car_zonly_atmosphere(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tz,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

!--- integrate through grid
  do while(photon%inside)

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     del    = tz - d
     tau    = tau + del * rhokap
     dtauH  = del * rhokapH
     d      = tz
     if (tau >= tau_in) then
        if (rhokap > 0.0) then
           d_overshoot = (tau - tau_in)/rhokap
           d     = d - d_overshoot
           del   = del - d_overshoot
           dtauH = del * rhokapH
        endif
        xp = xp + d * kx
        yp = yp + d * ky
        zp = zp + d * kz
        exit
     endif
     kcell = kcell + kstep
     if (kcell < 1 .or. kcell > grid%nz) then
        photon%inside = .false.
        exit
     endif
     tz = tz + delz

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     kold = kcell
     u1   = u2
     !dold = d
  enddo

  if (.not. photon%inside) then
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        if (kcell > 1) then
           !$OMP ATOMIC UPDATE
           grid%Jout(ix)  = grid%Jout(ix)  + photon%wgt
        else
           !$OMP ATOMIC UPDATE
           grid%Jabs2(ix) = grid%Jabs2(ix) + photon%wgt
        endif
     endif
  endif

!--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_zonly_atmosphere

  subroutine raytrace_to_tau_car_atmosphere(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
  !--- for photon destory
  logical       :: destroyed
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

  !--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx   = (grid%xface(icell+1)-xp)/kx
     delx = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           photon%inside = .false.
           return
        endif
     endif
     istep = -1
     tx   = (grid%xface(icell)-xp)/kx
     delx = -grid%dx/kx
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty   = (grid%yface(jcell+1)-yp)/ky
     dely = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     jstep = -1
     ty   = (grid%yface(jcell)-yp)/ky
     dely = -grid%dy/ky
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           photon%inside = .false.
           return
        endif
     endif
     kstep = -1
     tz   = (grid%zface(kcell)-zp)/kz
     delz = -grid%dz/kz
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

  !--- destroyed --> Is the photon destroyed or not?
  destroyed = .false.

  !--- integrate through grid
  do while(photon%inside)

     !--- The Lya photon is assumed to be destroyed by planet's dense molecular layer.
     if (grid%mask(icell,jcell,kcell) == -1_int8) then
        destroyed = .true.
        exit
     endif

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del   = tx - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx = tx + delx
     else if (idx_min == 2) then
        del   = ty - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty = ty + dely
     else
        del   = tz - d
        tau   = tau + del * rhokap
        dtauH = del * rhokapH
        d     = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d   - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           xp = xp + d * kx
           yp = yp + d * ky
           zp = zp + d * kz
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
     !dold = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

  !--- if photon has been destroyed, add to Jabs2
  if (destroyed) then
     u1               = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq     = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix               = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jabs2(ix) = grid%Jabs2(ix) + photon%wgt
     endif
     photon%inside = .false.
  endif

  !--- xp, yp, zp = the current photon coordinates.
  photon%x = xp
  photon%y = yp
  photon%z = zp
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_atmosphere

  subroutine raytrace_to_tau_car_xysym_atmosphere(photon,grid,tau_in)
!--- Find the coordinates and cell indices corresponding to an input optical depth tau_in.
!---
!--- The algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012-09-24
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- 2017-05-15, Reflected at xy-, xz- and yz-planes.

  use define
  implicit none

  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

! Local variables
  integer :: icell,jcell,kcell
  integer :: iold, jold, kold
  integer :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: tau,d_overshoot,d,del
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap,u1,u2,DnuHK
  real(kind=wp) :: rhokapH
  ! for mean intensity calculation
  integer       :: ix
  real(kind=wp) :: dold,dtauH
  integer       :: idx_min
  !--- for photon destory
  logical       :: destroyed
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- xp, yp, zp = the current photon coordinates.
  xp = photon%x
  yp = photon%y
  zp = photon%z
  kx = photon%kx
  ky = photon%ky
  kz = photon%kz
  icell = photon%icell
  jcell = photon%jcell
  kcell = photon%kcell

  if (kx > 0.0_wp) then
     if (icell > grid%nx .and. grid%xface(icell) == xp) then
        ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
        ! Ignore the case of going out of the grid system.
        photon%inside = .false.
        return
     endif
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     istep = -1
     delx  = -grid%dx/kx
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           icell = grid%i0
           istep = 1
           kx    = abs(kx)
        endif
     endif
     if (istep == 1) then
        tx = (grid%xface(icell+1)-xp)/kx
     else
        tx = (grid%xface(icell)-xp)/kx
     endif
  else
     istep = 0
     tx   = hugest
     delx = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) == yp) then
        photon%inside = .false.
        return
     endif
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     jstep = -1
     dely  = -grid%dy/ky
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           jcell = grid%j0
           jstep = 1
           ky    = abs(ky)
        endif
     endif
     if (jstep == 1) then
        ty = (grid%yface(jcell+1)-yp)/ky
     else
        ty = (grid%yface(jcell)-yp)/ky
     endif
  else
     jstep = 0
     ty   = hugest
     dely = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) == zp) then
        photon%inside = .false.
        return
     endif
     kstep = 1
     tz   = (grid%zface(kcell+1)-zp)/kz
     delz = grid%dz/kz
  else if (kz < 0.0_wp) then
     kstep = -1
     delz  = -grid%dz/kz
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           kcell = grid%k0
           kstep = 1
           kz    = abs(kz)
        endif
     endif
     if (kstep == 1) then
        tz = (grid%zface(kcell+1)-zp)/kz
     else
        tz = (grid%zface(kcell)-zp)/kz
     endif
  else
     kstep = 0
     tz   = hugest
     delz = hugest
  endif

  !--- d   = cumulative path length
  !--- tau = cumulative optical depth
  tau    = 0.0_wp
  d      = 0.0_wp
  dold   = 0.0_wp
  !tauold = 0.0_wp
  iold   = icell
  jold   = jcell
  kold   = kcell
  u1     = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz

  !--- destroyed --> Is the photon destroyed or not?
  destroyed = .false.

!--- integrate through grid
  do while(photon%inside)

     !--- The Lya photon is assumed to be destroyed by planet's dense molecular layer.
     if (grid%mask(icell,jcell,kcell) == -1_int8) then
        destroyed = .true.
        exit
     endif

     DnuHK   = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokapH = grid%rhokap(icell,jcell,kcell)*(voigt(photon%xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                               voigt(photon%xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) then
        rhokap = rhokapH + grid%rhokapD(icell,jcell,kcell)
     else
        rhokap = rhokapH
     endif

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del    = tx - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tx
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        icell = icell + istep
        if (icell < 1) then
           icell = grid%i0
           istep = 1
           kx    = -kx
        else if (icell > grid%nx) then
           photon%inside = .false.
           exit
        endif
        tx     = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = ty
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        jcell = jcell + jstep
        if (jcell < 1) then
           jcell = grid%j0
           jstep = 1
           ky    = -ky
        else if (jcell > grid%ny) then
           photon%inside = .false.
           exit
        endif
        ty     = ty + dely
     else
        del    = tz - d
        tau    = tau + del * rhokap
        dtauH  = del * rhokapH
        d      = tz
        if (tau >= tau_in) then
           if (rhokap > 0.0) then
              d_overshoot = (tau - tau_in)/rhokap
              d     = d - d_overshoot
              del   = del - d_overshoot
              dtauH = del * rhokapH
           endif
           exit
        endif
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) then
           photon%inside = .false.
           exit
        endif
        tz    = tz + delz
     endif

#ifdef CALCJ
     call add_to_J(photon,grid,iold,jold,kold,del)
#endif
#ifdef CALCPnew
     call add_to_Pnew(photon,grid,iold,jold,kold,dtauH)
#endif

     u2 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = (photon%xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold   = icell
     jold   = jcell
     kold   = kcell
     u1     = u2
     !dold   = d
  enddo

  if (.not. photon%inside) then
     icell = iold
     jcell = jold
     kcell = kold
  endif

#ifdef CALCJ
  call add_to_J(photon,grid,icell,jcell,kcell,del)
#endif
#ifdef CALCPnew
  call add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
#endif

  ! if photon escapes from the system, then transform the frequency from the fluid rest frame to the lab frame.
  ! comment added on 2017-04-28.
  if (.not. photon%inside) then
     !u1 = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     u1 = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix) = grid%Jout(ix) + photon%wgt
     endif
  endif

  !--- if photon has been destroyed, add to Jabs2
  if (destroyed) then
     u1               = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     photon%xfreq     = photon%xfreq + u1
     photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix               = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jabs2(ix) = grid%Jabs2(ix) + photon%wgt
     endif
     photon%inside = .false.
  endif

!--- xp, yp, zp = the current photon coordinates.
  ! here, photon%kx, etc should be used instead of kx,ky,kz
  photon%x  = photon%x + d*photon%kx
  photon%y  = photon%y + d*photon%ky
  photon%z  = photon%z + d*photon%kz
  if (photon%x < grid%xmin) photon%x = -photon%x
  if (photon%y < grid%ymin) photon%y = -photon%y
  photon%kx = kx
  photon%ky = ky
  photon%kz = kz
  photon%icell = icell
  photon%jcell = jcell
  photon%kcell = kcell

  return
  end subroutine raytrace_to_tau_car_xysym_atmosphere

  subroutine raytrace_to_edge_car_atmosphere(photon0,grid,tau)
!--- Find the cumulative optical depth to the edge of grid.
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  use, intrinsic :: ieee_arithmetic
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: tau

! local variables
  !-- note that exp(-tau_huge) = 0.0
  !-- at this large optical depth peeling-off contribution is practically zero.
  !-- Hence, no need to integrate further.
  real(kind=wp), parameter :: tau_huge = 745.2_wp
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap, DnuHK, u1,u2,xfreq
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau        = 0.0_wp
  d          = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz
  xfreq = photon0%xfreq

  !--- integrate through grid
  do while(.true.)

     !--- The Lya photon is assumed to be destroyed by planet's dense molecular layer.
     if (grid%mask(icell,jcell,kcell) == -1_int8) then
        tau = ieee_value(0.0_wp, ieee_positive_inf)
        exit
     endif

     DnuHK  = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokap = grid%rhokap(icell,jcell,kcell)*(voigt(xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                              voigt(xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))
     if (par%DGR > 0.0_wp) rhokap = rhokap + grid%rhokapD(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau   = tau + (tx - d) * rhokap
        d     = tx
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        tau   = tau + (ty - d) * rhokap
        d     = ty
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        tau   = tau + (tz - d) * rhokap
        d     = tz
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif
     if (tau >= tau_huge) exit

     u2    = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     xfreq = (xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
  enddo
  return
  end subroutine raytrace_to_edge_car_atmosphere

  subroutine raytrace_to_edge_car_tau_HI_atmosphere(photon0,grid,tau_HI)
!--- Find the cumulative optical depth due to HI and dust, column density of HI to the edge of grid. (2020.10.21)
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  use, intrinsic :: ieee_arithmetic
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: tau_HI

! local variables
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: rhokap, DnuHK, u1,u2,xfreq
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  tau_HI  = 0.0_wp
  d       = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell
  u1    = grid%vfx(iold,jold,kold)*kx + grid%vfy(iold,jold,kold)*ky + grid%vfz(iold,jold,kold)*kz
  xfreq = photon0%xfreq

  !--- integrate through grid
  do while(.true.)

     !--- The Lya photon is assumed to be destroyed by planet's dense molecular layer.
     if (grid%mask(icell,jcell,kcell) == -1_int8) then
        tau_HI = ieee_value(0.0_wp, ieee_positive_inf)
        exit
     endif

     DnuHK  = grid%DnuHK_ref_half/grid%Dfreq(icell,jcell,kcell)
     rhokap = grid%rhokap(icell,jcell,kcell)*(voigt(xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))/3d0 + &
                                              voigt(xfreq-DnuHK, grid%voigt_a(icell,jcell,kcell))*(2d0/3d0))

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        tau_HI = tau_HI + (tx - d) * rhokap
        d      = tx
        icell  = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        tau_HI = tau_HI + (ty - d) * rhokap
        d      = ty
        jcell  = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        tau_HI = tau_HI + (tz - d) * rhokap
        d      = tz
        kcell  = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif

     u2    = grid%vfx(icell,jcell,kcell)*kx + grid%vfy(icell,jcell,kcell)*ky + grid%vfz(icell,jcell,kcell)*kz
     xfreq = (xfreq + u1) * grid%Dfreq(iold,jold,kold)/grid%Dfreq(icell,jcell,kcell) - u2

     iold = icell
     jold = jcell
     kold = kcell
     u1   = u2
  enddo
  return
  end subroutine raytrace_to_edge_car_tau_HI_atmosphere

  subroutine raytrace_to_edge_car_column_atmosphere(photon0,grid,N_HI,tau_dust)
!--- Find the cumulative optical depth due to HI and dust, column density of HI to the edge of grid. (2020.10.21)
!--- Note that photon0.(x0,y0,z0) and photon0.(icell0,jcell0,kcell0) do not change.
!---
!--- Basic algorithm is based on the fast Voxel traversal algorithm
!--- by Amanatides & Woo (1987), Proc. of Eurographics, 87, 3
!---
!--- Author: Kwang-Il Seon, 2012
!--- Fixed minor bugs at grid boundaries: 2013-05-24
!--- Frequency shift due to the spatial variation of temperature and fluid velocity: 2018-01-10

  use define
  use, intrinsic :: ieee_arithmetic
  implicit none

  type(photon_type), intent(in) :: photon0
  type(grid_type),   intent(in) :: grid
  real(kind=wp),    intent(out) :: N_HI, tau_dust

! local variables
  integer       :: icell,jcell,kcell
  integer       :: iold,jold,kold
  integer       :: istep,jstep,kstep
  real(kind=wp) :: xp,yp,zp,kx,ky,kz
  real(kind=wp) :: d
  real(kind=wp) :: tx,ty,tz,delx,dely,delz
  real(kind=wp) :: u1,u2,xfreq
  real(kind=wp) :: rho,rhokapD,del
  integer       :: idx_min
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

!--- (xp,yp,zp) and (icell,jcell,kcell) = the photon coordinates and cell indices
  xp = photon0%x
  yp = photon0%y
  zp = photon0%z
  kx = photon0%kx
  ky = photon0%ky
  kz = photon0%kz
  icell = photon0%icell
  jcell = photon0%jcell
  kcell = photon0%kcell

!--- tau = cumulative optical depth
!--- d   = cumulative path length
  N_HI     = 0.0_wp
  tau_dust = 0.0_wp
  d        = 0.0_wp

  if (kx > 0.0_wp) then
     ! icell = nx+1 & xface(nx+1) = xp & kx > 0.0
     ! Ignore the case of going out of the grid system.
     if (icell > grid%nx .and. grid%xface(icell) <= xp) return
     istep = 1
     tx    = (grid%xface(icell+1)-xp)/kx
     delx  = grid%dx/kx
  else if (kx < 0.0_wp) then
     if (grid%xface(icell) == xp) then
        if (icell > 1) then
           ! icell = 2, 3, or, ... & xface(icell) = xp & kx < 0.0
           ! Skip the case of tau = 0
           icell = icell - 1
        else
           ! icell = 1 & xface(1) = xp & kx < 0.0
           ! Ignore the case of going out of the grid system.
           return
        endif
     endif
     istep = -1
     tx    = (grid%xface(icell)-xp)/kx
     delx  = -grid%dx/kx
  else
     istep = 0
     tx    = hugest
     delx  = hugest
  endif

  if (ky > 0.0_wp) then
     if (jcell > grid%ny .and. grid%yface(jcell) <= yp) return
     jstep = 1
     ty    = (grid%yface(jcell+1)-yp)/ky
     dely  = grid%dy/ky
  else if (ky < 0.0_wp) then
     if (grid%yface(jcell) == yp) then
        if (jcell > 1) then
           jcell = jcell - 1
        else
           return
        endif
     endif
     jstep = -1
     ty    = (grid%yface(jcell)-yp)/ky
     dely  = -grid%dy/ky
  else
     jstep = 0
     ty    = hugest
     dely  = hugest
  endif

  if (kz > 0.0_wp) then
     if (kcell > grid%nz .and. grid%zface(kcell) <= zp) return
     kstep = 1
     tz    = (grid%zface(kcell+1)-zp)/kz
     delz  = grid%dz/kz
  else if (kz < 0.0_wp) then
     if (grid%zface(kcell) == zp) then
        if (kcell > 1) then
           kcell = kcell - 1
        else
           return
        endif
     endif
     kstep = -1
     tz    = (grid%zface(kcell)-zp)/kz
     delz  = -grid%dz/kz
  else
     kstep = 0
     tz    = hugest
     delz  = hugest
  endif

  !--- 2018-01-10, frequency shift due to spatial variation of velocity and temperature had not been took into account.
  iold  = icell
  jold  = jcell
  kold  = kcell

  !--- integrate through grid
  do while(.true.)

     !--- The Lya photon is assumed to be destroyed by planet's dense molecular layer.
     if (grid%mask(icell,jcell,kcell) == -1_int8) then
        N_HI = ieee_value(0.0d0, ieee_positive_inf)
        if (par%DGR > 0.0_wp) tau_dust = ieee_value(0.0d0, ieee_positive_inf)
        exit
     endif

     rho    = grid%rhokap(icell,jcell,kcell)*grid%Dfreq(icell,jcell,kcell) / par%cross0
     if (par%DGR > 0.0_wp) rhokapD = grid%rhokapD(icell,jcell,kcell)

     idx_min = minloc([tx,ty,tz], dim=1)

     if (idx_min == 1) then
        del    = tx - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = tx
        icell = icell + istep
        if (icell < 1 .or. icell > grid%nx) exit
        tx    = tx + delx
     else if (idx_min == 2) then
        del    = ty - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = ty
        jcell = jcell + jstep
        if (jcell < 1 .or. jcell > grid%ny) exit
        ty = ty + dely
     else
        del    = tz - d
        N_HI   = N_HI   + del * rho
        if (par%DGR > 0.0_wp) tau_dust = tau_dust + del * rhokapD
        d     = tz
        kcell = kcell + kstep
        if (kcell < 1 .or. kcell > grid%nz) exit
        tz    = tz + delz
     endif

     iold = icell
     jold = jcell
     kold = kcell
  enddo
  return
  end subroutine raytrace_to_edge_car_column_atmosphere

#ifdef CALCJ
  subroutine add_to_J(photon,grid,icell,jcell,kcell,del)
  use define
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  integer,           intent(in)    :: icell,jcell,kcell
  real(kind=wp),     intent(in)    :: del
  integer                          :: ix

  photon%xfreq_ref = photon%xfreq * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
  ix               = floor((photon%xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
  if (ix > 0 .and. ix <= grid%nxfreq .and. &
      icell > 0 .and. icell <= grid%nx .and. jcell > 0 .and. jcell <= grid%ny .and. kcell > 0 .and. kcell <= grid%nz) then
     if (grid%rhokap(icell,jcell,kcell) > 0.0_wp) then
        select case (grid%geometry_JPa)
        case (3)
           !$OMP ATOMIC UPDATE
           grid%J(ix,icell,jcell,kcell)                = grid%J(ix,icell,jcell,kcell) + del * photon%wgt
        case (2)
           !$OMP ATOMIC UPDATE
           grid%J2(ix,grid%ind_cyl(icell,jcell),kcell) = grid%J2(ix,grid%grid%ind_cyl(icell,jcell),kcell) + del * photon%wgt
        case (1)
           !$OMP ATOMIC UPDATE
           grid%J1(ix,grid%ind_sph(icell,jcell,kcell)) = grid%J1(ix,grid%ind_sph(icell,jcell,kcell)) + del * photon%wgt
        case (-1)
           !$OMP ATOMIC UPDATE
           grid%J1(ix,kcell)                           = grid%J1(ix,kcell) + del * photon%wgt
        end select
     endif
  endif
  end subroutine add_to_J
#endif

#ifdef CALCPnew
  subroutine add_to_Pnew(photon,grid,icell,jcell,kcell,dtauH)
  use define
  type(photon_type), intent(in)    :: photon
  type(grid_type),   intent(inout) :: grid
  integer,           intent(in)    :: icell,jcell,kcell
  real(kind=wp),     intent(in)    :: dtauH
  real(kind=wp)                    :: rhokap
  if (icell > 0 .and. icell <= grid%nx .and. jcell > 0 .and. jcell <= grid%ny .and. kcell > 0 .and. kcell <= grid%nz) then
     if (grid%rhokap(icell,jcell,kcell) > 0.0_wp) then
        !--- Convert rhokap into density unit * distance2cm. (number/cm^3) * distance2cm.
        rhokap = grid%rhokap(icell,jcell,kcell) * grid%Dfreq(icell,jcell,kcell)/par%cross0
        select case (grid%geometry_JPa)
        case (3)
           !$OMP ATOMIC UPDATE
           grid%Pa_new(icell,jcell,kcell)               = grid%Pa_new(icell,jcell,kcell) + dtauH * photon%wgt / rhokap
        case (2)
           !$OMP ATOMIC UPDATE
           grid%P2_new(grid%ind_cyl(icell,jcell),kcell) = grid%P2_new(grid%ind_cyl(icell,jcell),kcell) + dtauH * photon%wgt / rhokap
        case (1)
           !$OMP ATOMIC UPDATE
           grid%P1_new(grid%ind_sph(icell,jcell,kcell)) = grid%P1_new(grid%ind_sph(icell,jcell,kcell)) + dtauH * photon%wgt / rhokap
        case (-1)
           !$OMP ATOMIC UPDATE
           grid%P1_new(kcell)                           = grid%P1_new(kcell) + dtauH * photon%wgt / rhokap
        end select
     endif
  endif
  end subroutine add_to_Pnew
#endif

end module raytrace
