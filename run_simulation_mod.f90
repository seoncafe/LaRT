module run_simulation_mod
  use define
  use grid_mod
  use generate_photon_mod
  use random
  use scatter_mod
  use utility
  use mpi
contains
  !============================================
  subroutine run_master_slave(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, ii
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp) :: dtime
  logical :: first

  !--- for MPI
  integer :: ierr
  integer :: master, worker, irank, ans, tag
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=int64) :: numsent, numreceived, numdone
  !integer(kind=int64) :: numsent, numreceived

  !--- Photon loop
  master   = 0
  if (mpar%p_rank == master) then
     numsent = 0_int64
     !--- note 1,2,...,nproc-1
     do irank = 1, min(par%nphotons,int(mpar%nproc-1,int64))
        tag     = 1
        numsent = numsent + par%num_send_at_once
        call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
     enddo
     numdone = 0_int64
     do ip = 1, par%nphotons, par%num_send_at_once
        call MPI_RECV(ans,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        worker  = status(MPI_SOURCE)
        !--- cap printed counter at nphotons (the final chunk may be partial)
        numdone = min(numdone + par%num_send_at_once, par%nphotons)
        if (mod(numdone,int(par%nprint,int64)) == 0) then
           call time_stamp(dtime)
           write(6,'(es14.5,a,f8.3,a)') dble(numdone),' photons calculated: ',dtime/60.0_wp,' mins'
        endif
        if (numsent < par%nphotons) then
           !--- modified to send numsent (2019-08-19)
           tag     = 1
           numsent = numsent + par%num_send_at_once
           call MPI_SEND(numsent,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        else
           tag = 0
           call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        endif
     enddo
  else
     do while(.true.)
        call MPI_RECV(numreceived,1,MPI_INTEGER8,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        if (status(MPI_TAG) == 0) exit
        do ii = 1, par%num_send_at_once
           ip = numreceived - par%num_send_at_once + ii
           !--- skip indices that exceed the requested total (final chunk may
           !    be partial when nphotons is not a multiple of num_send_at_once)
           if (ip > par%nphotons) exit
           !++++++++++++++++++++++++++++++++
           !--- Main Part of the simulation
           !--- Release photon
           photon%id = ip
           call generate_photon(grid,photon)
           if (par%save_all_photons) call make_all_initial_photons(photon)

           first = .true.
           do while(photon%inside)
              !--- Find scattering location of tau
              !if (photon%nscatt_gas + photon%nscatt_dust == 0.0_wp) then
              if (first) then
                 !--- Force photon to scatter at optical depth tau before edge of grid
                 call raytrace_to_edge(photon,grid,tau0)
                 call add_escaped_fraction_to_Jout(photon,grid,tau0)
                 wgt1       = 1.0_wp - exp(-tau0)
                 photon%wgt = photon%wgt * wgt1
                 !--- tau0 = 0 & wgt1 = 0 occurs when photons are generated at locations where rho = 0. (2025.10.11)
                 !tau        = -log(1.0_wp - rand_number()*wgt1)
                 if (tau0 > 0.0_wp) then
                    tau = -log(1.0_wp - rand_number()*wgt1)
                 else
                    tau = hugest
                 endif
                 first = .false.
              else
                 tau = -log(rand_number())
              endif
              call raytrace_to_tau(photon,grid,tau)

              if (photon%inside) then
                 call scattering(photon,grid)
              endif
           enddo
           !--- End of the main part of the simulation
           !++++++++++++++++++++++++++++++++
           par%nscatt_gas  = par%nscatt_gas  + photon%nscatt_gas
           par%nscatt_dust = par%nscatt_dust + photon%nscatt_dust
           par%nscatt_tot  = par%nscatt_gas  + par%nscatt_dust
           par%nrejected   = par%nrejected   + photon%nrejected
           if (trim(par%source_geometry) == 'stellar_illumination' .or. trim(par%source_geometry) == 'point_illumination') then
              par%flux_factor = par%flux_factor + photon%flux_factor
           endif
           if (par%save_all_photons) call make_all_photons(photon)
        enddo
        !if (mod(ip,int(par%nprint,int64)) == 0) then
        !   call time_stamp(dtime)
        !   write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
        !endif
        ans = 1
        tag = 1
        call MPI_SEND(ans,1,MPI_INTEGER,master,tag,MPI_COMM_WORLD,ierr)
     enddo
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine run_master_slave

  !============================================
  subroutine run_equal_number(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp) :: dtime
  logical :: first

  !--- for MPI
  integer :: ierr

  !--- Photon loop
  !do ip=mpar%p_rank+1,par%nphotons,mpar%nproc
  do ip=mpar%p_rank+1,par%nphotons+mpar%nproc,mpar%nproc
     if (ip > par%nphotons) exit
     !++++++++++++++++++++++++++++++++
     !--- Main Part of the simulation
     !--- Release photon
     photon%id = ip
     call generate_photon(grid,photon)
     if (par%save_all_photons) call make_all_initial_photons(photon)

     first = .true.
     do while(photon%inside)
        !--- Find scattering location of tau
        !if (photon%nscatt_gas + photon%nscatt_dust == 0) then
        if (first) then
           !--- Force photon to scatter at optical depth tau before edge of grid
           call raytrace_to_edge(photon,grid,tau0)
           call add_escaped_fraction_to_Jout(photon,grid,tau0)
           wgt1       = 1.0_wp - exp(-tau0)
           photon%wgt = photon%wgt * wgt1
           !tau        = -log(1.0_wp - rand_number()*wgt1)
           if (tau0 > 0.0_wp) then
              tau = -log(1.0_wp - rand_number()*wgt1)
           else
              tau = hugest
           endif
           first = .false.
        else
           tau = -log(rand_number())
        endif
        call raytrace_to_tau(photon,grid,tau)

        if (photon%inside) then
           call scattering(photon,grid)
        endif
     enddo
     !--- End of the main part of the simulation
     !++++++++++++++++++++++++++++++++
     par%nscatt_gas  = par%nscatt_gas  + photon%nscatt_gas
     par%nscatt_dust = par%nscatt_dust + photon%nscatt_dust
     par%nscatt_tot  = par%nscatt_gas  + par%nscatt_dust
     par%nrejected   = par%nrejected   + photon%nrejected
     if (trim(par%source_geometry) == 'stellar_illumination' .or. trim(par%source_geometry) == 'point_illumination') then
        par%flux_factor = par%flux_factor + photon%flux_factor
     endif
     if (par%save_all_photons) call make_all_photons(photon)

     if (mod(ip,int(par%nprint,int64)) == 0 .or. ip == par%nphotons) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine run_equal_number
  !============================================
  !--- add escaped photon fraction to Jout.
  !--- This is done when the forced first scattering technique is used (2021.08.10).
  subroutine add_escaped_fraction_to_Jout(photon,grid,tau0)
  use octree_mod, only: amr_grid
  implicit none
  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(in)    :: photon
  real(wp),          intent(in)    :: tau0
  real(wp) :: u1, xfreq_ref, wgt_esc, mu_val
  integer  :: icell, jcell, kcell, il, ix, imu
  wgt_esc = photon%wgt * exp(-tau0)
  if (par%save_Jmu) then
     mu_val = photon%kz
     if (par%xyz_symmetry) mu_val = abs(mu_val)
     imu = floor((mu_val - par%mu_min)/par%dmu) + 1
     if (imu < 1)       imu = 1
     if (imu > par%nmu) imu = par%nmu
  endif
  if (par%use_amr_grid) then
     il = photon%icell_amr
     u1 = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
     xfreq_ref = (photon%xfreq + u1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
     ix = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        amr_grid%Jout(ix) = amr_grid%Jout(ix) + wgt_esc
        if (par%save_Jmu .and. allocated(amr_grid%Jmu)) &
           amr_grid%Jmu(ix, imu) = amr_grid%Jmu(ix, imu) + wgt_esc
     end if
  else
     icell     = photon%icell
     jcell     = photon%jcell
     kcell     = photon%kcell
     u1        = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
     xfreq_ref = (photon%xfreq + u1) * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
     ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        grid%Jout(ix) = grid%Jout(ix) + wgt_esc
        if (par%save_Jmu .and. associated(grid%Jmu)) &
           grid%Jmu(ix, imu) = grid%Jmu(ix, imu) + wgt_esc
     endif
  end if
  end subroutine add_escaped_fraction_to_Jout
  !--------------------------------------------------
  subroutine make_all_initial_photons(photon)
  use utility
  implicit none
  type(photon_type), intent(inout) :: photon
  real(wp) :: dist, xp, yp, zp
  real(wp) :: rr, rk, det
  real(wp) :: mx, my, mz, mm

  allph%xfreq1(photon%id) = photon%xfreq

  if (trim(par%source_geometry) /= 'point') then
     if (par%rmax > 0.0_wp) then
        rr  = photon%x**2 + photon%y**2 + photon%z**2
        if (rr > par%rmax**2) then
           rk  = photon%x * photon%kx + photon%y * photon%ky + photon%z * photon%kz
           det = rk**2 - (rr - par%rmax**2)
           if (det < 0.0_wp) then
              !--- 2020.09.20
              !--- this condition can occur because we are implementing a sphere in a rectangular coordinate system.
              dist = 0.0_wp
           else
              dist = -rk + sqrt(max(0.0_wp, det))
           endif
        else
           dist = 0.0_wp
        endif
        xp = photon%x + dist * photon%kx
        yp = photon%y + dist * photon%ky
        zp = photon%z + dist * photon%kz
     else
        xp = photon%x
        yp = photon%y
        zp = photon%z
     endif

     rk = xp*photon%kx + yp*photon%ky + zp*photon%kz
     mx = xp - rk * photon%kx
     my = yp - rk * photon%ky
     mz = zp - rk * photon%kz
     mm = sqrt(mx**2 + my**2 + mz**2)

     allph%rp0(photon%id) = mm
  endif
  end subroutine make_all_initial_photons
  !--------------------------------------------------
  subroutine make_all_photons(photon)
  use utility
  implicit none
  type(photon_type), intent(in) :: photon
  real(wp) :: dist, xp, yp, zp, cosp, sinp, cos2p, sin2p
  real(wp) :: rr, rk, det
  real(wp) :: mx, my, mz, mm

  if (par%rmax > 0.0_wp) then
     rr  = photon%x**2 + photon%y**2 + photon%z**2
     if (rr > par%rmax**2) then
        rk  = photon%x * photon%kx + photon%y * photon%ky + photon%z * photon%kz
        det = rk**2 - (rr - par%rmax**2)
        if (det < 0.0_wp) then
           !--- 2020.09.20
           !--- this condition can occur because we are implementing a sphere in a rectangular coordinate system.
           dist = 0.0_wp
        else
           dist = -rk + sqrt(max(0.0_wp, det))
        endif
     else
        dist = 0.0_wp
     endif
     xp = photon%x + dist * photon%kx
     yp = photon%y + dist * photon%ky
     zp = photon%z + dist * photon%kz
  else
     xp = photon%x
     yp = photon%y
     zp = photon%z
  endif

  rk = xp*photon%kx + yp*photon%ky + zp*photon%kz
  mx = xp - rk * photon%kx
  my = yp - rk * photon%ky
  mz = zp - rk * photon%kz
  mm = sqrt(mx**2 + my**2 + mz**2)

  allph%rp(photon%id)          = mm
  allph%xfreq2(photon%id)      = photon%xfreq_ref
  allph%nscatt_gas(photon%id)  = photon%nscatt_gas
  allph%nscatt_dust(photon%id) = photon%nscatt_dust

  if (par%use_stokes) then
     if (mm > 0.0_wp) then
        mx = mx/mm
        my = my/mm
        mz = mz/mm
        cosp  = mx*photon%mx + my*photon%my + mz*photon%mz
        sinp  = mx*photon%nx + my*photon%ny + mz*photon%nz
        cos2p = 2.0_wp*cosp*cosp - 1.0_wp
        sin2p = 2.0_wp*sinp*cosp
     else
        cosp  = 1.0_wp
        sinp  = 0.0_wp
        cos2p = 1.0_wp
        sin2p = 0.0_wp
     endif

     allph%I(photon%id)        = photon%wgt
     allph%Q(photon%id)        = ( cos2p*photon%Q + sin2p*photon%U)*photon%wgt
     allph%U(photon%id)        = (-sin2p*photon%Q + cos2p*photon%U)*photon%wgt
     allph%V(photon%id)        = photon%V*photon%wgt
  endif
  end subroutine make_all_photons
!--------------------------------------------------

end module run_simulation_mod
