module run_simulation_mod
  use define
  use grid_mod
  use photon_mod
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
     do irank = 1, min(par%nphotons,mpar%nproc-1)
        tag     = 1
        numsent = numsent + par%num_send_at_once
        call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
     enddo
     numdone = 0_int64
     do ip = 1, par%nphotons, par%num_send_at_once
        call MPI_RECV(ans,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
        worker  = status(MPI_SOURCE)
        numdone = numdone + par%num_send_at_once
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
           !++++++++++++++++++++++++++++++++
           !--- Main Part of the simulation
           !--- Release photon
           photon%id = ip
           call gen_photon(grid,photon)
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
     call gen_photon(grid,photon)
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

     if (mod(ip,par%nprint) == 0 .or. ip == par%nphotons) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine run_equal_number

  !============================================
  ! This method seems to be slower than the others. (2023.08.06)
  ! Why?
  subroutine run_method3(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, ii
  real(kind=wp) :: tau, tau0, wgt1
  real(kind=wp) :: dtime
  logical :: first

  !--- for MPI
  integer :: win,ierr
  integer :: master
  integer(kind=int64) :: ncount
  integer(kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofint

  CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lowerbound,sizeofint, ierr)
  CALL MPI_WIN_CREATE(ncount,sizeofint,int(sizeofint),MPI_INFO_NULL,MPI_COMM_WORLD,win,ierr)

  master = 0
  ncount = 0
  do while(.true.)
     call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,master,0,win,ierr)
     call MPI_Fetch_and_Op(par%num_send_at_once,ncount,MPI_INTEGER8,master,0_MPI_ADDRESS_KIND, MPI_SUM,win,ierr)
     call MPI_WIN_UNLOCK(master,win,ierr)
     if (ncount >= par%nphotons) exit

     do ii = 1, par%num_send_at_once
        ip = ncount + ii
        !++++++++++++++++++++++++++++++++
        !--- Main Part of the simulation
        !--- Release photon
        photon%id = ip
        call gen_photon(grid,photon)
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
     enddo

     if (mod(ip,int(par%nprint,int64)) == 0) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine run_method3

  !============================================
  !--- add escaped photon fraction to Jout.
  !--- This is done when the forced first scattering technique is used (2021.08.10).
  subroutine add_escaped_fraction_to_Jout(photon,grid,tau0)
  implicit none
  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(in)    :: photon
  real(wp),          intent(in)    :: tau0
  real(wp) :: u1, xfreq_ref
  integer  :: icell, jcell, kcell, ix
  icell     = photon%icell
  jcell     = photon%jcell
  kcell     = photon%kcell
  u1        = grid%vfx(icell,jcell,kcell)*photon%kx + grid%vfy(icell,jcell,kcell)*photon%ky + grid%vfz(icell,jcell,kcell)*photon%kz
  xfreq_ref = (photon%xfreq + u1) * (grid%Dfreq(icell,jcell,kcell) / grid%Dfreq_ref)
  ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
  if (ix >= 1 .and. ix <= grid%nxfreq) then
     grid%Jout(ix) = grid%Jout(ix) + photon%wgt * exp(-tau0)
  endif
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
              dist = -rk + sqrt(det)
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
           dist = -rk + sqrt(det)
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
  !============================================
end module run_simulation_mod
