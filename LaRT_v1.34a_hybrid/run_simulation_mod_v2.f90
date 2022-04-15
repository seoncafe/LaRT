module run_simulation_mod
  use define
  use grid_mod
  use photon_mod
  use random
  use scatter_mod
  use utility
  use mpi
  use omp_lib
contains
  !============================================
  !--- Now, this routine works fine (2021.05.05).
  !--- 2021.05.05 (bug-fixed) prevent a thread from receiving messages sent by itself.
  !--- 2020.12.17 (bug-fixed) "ans" should be int64 because several threads in a single node do "send" and "receive" simultaneously.
  !--- But, why int32 ans worked in the combination of "intel 2018 + ubuntu 16.04"?
  subroutine run_master_slave(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, ii
  real(kind=wp) :: tau, nscatt_HI, nscatt_dust, nrejected, flux_factor
  real(kind=wp) :: dtime

  !--- for MPI & OpenMP
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr
  integer :: my_threadid, master, irank, ithread, ithread0
  integer :: tag, worker
  integer(kind=int64) :: ans
  integer(kind=int64) :: numsent, numreceived, numdone

  !--- Initialize
  my_threadid = omp_get_thread_num()
  master      = 0
  nscatt_HI   = 0.0
  nscatt_dust = 0.0
  nrejected   = 0.0
  flux_factor = 0.0

  !--- master part
  if (mpar%p_rank == master .and. my_threadid == 0) then
     numsent = 0_int64
     !--- note 0,1,2,...,num_nodes-1
     do irank = 0, mpar%num_nodes-1
        if (irank == 0) then
           ithread0 = 1
        else
           ithread0 = 0
        endif
        do ithread = ithread0, mpar%num_threads(irank+1)-1
           tag     = ithread
           numsent = numsent + par%num_send_at_once
           if (numsent > par%nphotons) numsent = -1_int64
           call MPI_SEND(numsent,1,MPI_INTEGER8,irank,tag,MPI_COMM_WORLD,ierr)
        enddo
     enddo
     numdone = 0_int64
     do ip = 1, par%nphotons, par%num_send_at_once
        do while (.true.)
           call MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
           worker = status(MPI_SOURCE)
           tag    = status(MPI_TAG)
           !--- do not receive a message sent by the master itself (p_rank = 0, threadid = 0).
           !--- receive messages only from other threads (tag >= tau_offset) (2021.05.05).
           if (tag >= mpar%max_num_threads) then
              call MPI_RECV(ans,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,status,ierr)
              !--- will send the threadid. (Note a negative tag is not allowed in MPI. 2021.05.05)
              tag = tag - mpar%max_num_threads
              exit
           endif
        enddo
        numdone = numdone + par%num_send_at_once
        if (mod(numdone,int(par%nprint,int64)) == 0) then
           call time_stamp(dtime)
           write(6,'(es14.5,a,f8.3,a)') dble(numdone),' photons calculated: ',dtime/60.0_wp,' mins'
        endif
        if (numsent < par%nphotons) then
           numsent = numsent + par%num_send_at_once
           call MPI_SEND(numsent,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        else
           call MPI_SEND(-1_int64,1,MPI_INTEGER8,worker,tag,MPI_COMM_WORLD,ierr)
        endif
     enddo
  !--- worker (slave) part
  else
     do while(.true.)
        !--- here, will receive messages with tag = my_threadid (2021.05.05)
        call MPI_RECV(numreceived,1,MPI_INTEGER8,master,my_threadid,MPI_COMM_WORLD,status,ierr)
        if (numreceived == -1_int64) exit
        do ii = 1, par%num_send_at_once
           ip = numreceived - par%num_send_at_once + ii
           !--- Condition for the case of mod(par%nphotons, par%num_send_at_once) /= 0. (2020.12.05)
           if (ip > par%nphotons) exit
           !++++++++++++++++++++++++++++++++
           !--- Main Part of the simulation
           !--- Release photon
           photon%id = ip
           call gen_photon(grid,photon)
           if (par%save_all_photons) call make_all_initial_photons(photon)

           do while(photon%inside)
              !--- Find scattering location of tau
              tau = -log(rand_number())
              call raytrace_to_tau(photon,grid,tau)

              if (photon%inside) then
                 call scattering(photon,grid)
              endif
           enddo
           !--- End of the main part of the simulation
           !++++++++++++++++++++++++++++++++
           nscatt_HI   = nscatt_HI   + photon%nscatt_HI
           nscatt_dust = nscatt_dust + photon%nscatt_dust
           nrejected   = nrejected   + photon%nrejected
           if (trim(par%source_geometry) == 'stellar_illumination') flux_factor = flux_factor + photon%flux_factor
           if (par%save_all_photons) call make_all_photons(photon)
           !++++++++++++++++++++++++++++++++
        enddo
        !--- do not send a tag that can be received by itself which sended the message.
        !--- will send the threadid + tag_offset. (Note a negative tag is not allowed. 2021.05.05)
        ans = 1_int64
        tag = my_threadid + mpar%max_num_threads
        call MPI_SEND(ans,1,MPI_INTEGER8,master,tag,MPI_COMM_WORLD,ierr)
     enddo
  endif
  !$OMP BARRIER
  !$OMP ATOMIC
  par%nscatt_HI   = par%nscatt_HI   + nscatt_HI
  !$OMP ATOMIC
  par%nscatt_dust = par%nscatt_dust + nscatt_dust
  !$OMP ATOMIC
  par%nrejected   = par%nrejected   + nrejected
  !$OMP ATOMIC
  par%flux_factor = par%flux_factor + flux_factor
  !$OMP BARRIER
  par%nscatt_tot  = par%nscatt_HI + par%nscatt_dust
  end subroutine run_master_slave
  !============================================
  subroutine run_equal_number(grid)
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  type(photon_type)   :: photon
  integer(kind=int64) :: ip
  real(kind=wp) :: tau
#ifdef __GFORTRAN__
  !--- reduction variable should be private in outer context? what about intel compiler?
  real(kind=wp), save :: nscatt_HI, nscatt_dust
#else
  real(kind=wp) :: nscatt_HI, nscatt_dust, nrejected, flux_factor
#endif
  real(kind=wp) :: dtime

  !--- Photon loop
  nscatt_HI   = 0.0
  nscatt_dust = 0.0
  nrejected   = 0.0
  flux_factor = 0.0
  !$OMP DO &
  !$OMP private(tau, photon, dtime, ip) &
  !$OMP reduction(+:nscatt_HI, nscatt_dust, nrejected, flux_factor) &
  !$OMP schedule(dynamic,1)
  do ip = mpar%p_rank + 1, par%nphotons, mpar%num_nodes
     !++++++++++++++++++++++++++++++++
     !--- Main Part of the simulation
     !--- Release photon
     photon%id = ip
     call gen_photon(grid,photon)
     if (par%save_all_photons) call make_all_initial_photons(photon)

     do while(photon%inside)
        !--- Find scattering location of tau
        tau = -log(rand_number())
        call raytrace_to_tau(photon,grid,tau)

        if (photon%inside) then
           call scattering(photon,grid)
        endif
     enddo
     !--- End of the main part of the simulation
     !++++++++++++++++++++++++++++++++
     nscatt_HI   = nscatt_HI   + photon%nscatt_HI
     nscatt_dust = nscatt_dust + photon%nscatt_dust
     nrejected   = nrejected   + photon%nrejected
     if (trim(par%source_geometry) == 'stellar_illumination') flux_factor = flux_factor + photon%flux_factor
     if (par%save_all_photons) call make_all_photons(photon)

     if (mod(ip, par%nprint) == 0 .or. ip == par%nphotons) then
        call time_stamp(dtime)
        write(6,'(es14.5,a,f8.3,a)') dble(ip),' photons calculated: ',dtime/60.0_wp,' mins'
     endif
  enddo
  !$OMP END DO
  !--- Note that reduced value is stored in MASTER thread.
  !$OMP MASTER
  par%nscatt_HI   = nscatt_HI
  par%nscatt_dust = nscatt_dust
  par%nscatt_tot  = par%nscatt_HI + par%nscatt_dust
  par%nrejected   = nrejected
  par%flux_factor = flux_factor
  !$OMP END MASTER
  end subroutine run_equal_number
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
  allph%nscatt_HI(photon%id)   = photon%nscatt_HI
  allph%nscatt_dust(photon%id) = photon%nscatt_dust

  if (par%use_stokes) then
     if (mm > 0.0_wp) then
        mx    = mx/mm
        my    = my/mm
        mz    = mz/mm
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

     allph%I(photon%id) = photon%wgt
     allph%Q(photon%id) = ( cos2p*photon%Q + sin2p*photon%U)*photon%wgt
     allph%U(photon%id) = (-sin2p*photon%Q + cos2p*photon%U)*photon%wgt
     allph%V(photon%id) = photon%V*photon%wgt
  endif
  end subroutine make_all_photons
!--------------------------------------------------

end module run_simulation_mod
