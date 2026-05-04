module sightline_tau_clump_mod
  !-----------------------------------------------------------------------
  ! Clump-medium equivalent of make_sightline_tau_outside (sightline_tau_rect).
  !
  ! For each observer pixel, a virtual photon is placed at the far face of
  ! the bounding sphere of radius `sphere_R` (clump_mod) and traced inward
  ! through the CSR clump grid to accumulate:
  !   tau_gas(nu, ix, iy)  -- gas optical depth per frequency bin
  !   N_gas(ix, iy)        -- HI column density
  !   tau_dust(ix, iy)     -- dust optical depth (if DGR > 0; currently 0
  !                          for clump mode since the inter-clump medium
  !                          is vacuum and the clump-internal dust path is
  !                          not yet implemented)
  !
  ! The pixel-loop body and FITS output mirror sightline_tau_rect; only
  ! the geometry (sphere, not box) and the underlying raytrace
  ! (raytrace_to_edge_*_clump) differ.
  !
  ! Output is written by write_sightline_tau_outside (sightline_tau_rect),
  ! which reads only par% and obs% fields (the grid argument is unused).
  !-----------------------------------------------------------------------
  use define
  use utility
  use fitsio_mod
  use memory_mod
  use clump_mod, only: sphere_R
  implicit none
  private
  public :: make_sightline_tau_clump

contains

  subroutine make_sightline_tau_clump(g)
    use mpi
    use raytrace_clump_mod, only: raytrace_to_edge_tau_gas_clump, &
                                  raytrace_to_edge_column_clump
    use sightline_tau_rect, only: write_sightline_tau_outside
    implicit none
    type(grid_type), intent(in) :: g    ! passed to write routine; unused there

    type(photon_type) :: pobs
    real(wp) :: tau_gas, N_gas, tau_dust
    real(wp) :: kx, ky, kz, kr
    real(wp) :: ox, oy, oz, b, c, disc, t_far
    real(wp) :: eps_pos
    integer  :: ix, iy, kk, i
    integer  :: loop, loop1, loop2, nsize
    character(len=4) :: filename_end

    eps_pos = 1.0e-9_wp * sphere_R

    !--- Allocate per-observer sight-line arrays (MPI shared memory).
    do i = 1, par%nobs
      call create_shared_mem(observer(i)%tau_gas, [par%nxfreq, par%nxim, par%nyim])
      call create_shared_mem(observer(i)%N_gas,   [par%nxim, par%nyim])
      if (par%DGR > 0.0_wp) &
        call create_shared_mem(observer(i)%tau_dust, [par%nxim, par%nyim])
    end do

    do i = 1, par%nobs
      nsize = observer(i)%nxim * observer(i)%nyim
      call loop_divide(nsize, mpar%nproc, mpar%p_rank, loop1, loop2)

      do loop = loop1, loop2
        call array_2D_indices(observer(i)%nxim, observer(i)%nyim, loop, ix, iy)

        !--- Pixel direction in the observer frame.
        kx = tan((ix - (observer(i)%nxim + 1.0_wp)/2.0_wp) * observer(i)%dxim / rad2deg)
        ky = tan((iy - (observer(i)%nyim + 1.0_wp)/2.0_wp) * observer(i)%dyim / rad2deg)
        kz = -1.0_wp
        kr = sqrt(kx*kx + ky*ky + kz*kz)
        kx = kx/kr;  ky = ky/kr;  kz = kz/kr

        !--- Rotate to the global (clump) frame.
        pobs%kx = observer(i)%rmatrix(1,1)*kx + observer(i)%rmatrix(2,1)*ky + observer(i)%rmatrix(3,1)*kz
        pobs%ky = observer(i)%rmatrix(1,2)*kx + observer(i)%rmatrix(2,2)*ky + observer(i)%rmatrix(3,2)*kz
        pobs%kz = observer(i)%rmatrix(1,3)*kx + observer(i)%rmatrix(2,3)*ky + observer(i)%rmatrix(3,3)*kz

        !--- Ray-sphere intersection: |obs + t*k|^2 = sphere_R^2
        !    Solve  t^2 + 2 b t + c = 0    where
        !       b = obs . k
        !       c = |obs|^2 - sphere_R^2
        !    The two real roots (when disc >= 0) bracket the chord
        !    through the sphere; the larger root is the far intersection.
        ox = observer(i)%x
        oy = observer(i)%y
        oz = observer(i)%z
        b  = ox*pobs%kx + oy*pobs%ky + oz*pobs%kz
        c  = ox*ox + oy*oy + oz*oz - sphere_R*sphere_R
        disc = b*b - c
        if (disc < 0.0_wp) cycle               ! sight-line misses the sphere
        t_far = -b + sqrt(disc)                ! far-side entry (looking back toward observer)
        if (t_far <= 0.0_wp) cycle              ! observer is past the sphere on this ray

        !--- Place the virtual photon at the far surface point.
        pobs%x = ox + pobs%kx * t_far
        pobs%y = oy + pobs%ky * t_far
        pobs%z = oz + pobs%kz * t_far

        !--- Reverse the direction: now points inward, toward the observer.
        pobs%kx = -pobs%kx
        pobs%ky = -pobs%ky
        pobs%kz = -pobs%kz

        !--- Nudge slightly inward so the entry is unambiguously inside the
        !    bounding sphere even with float round-off.
        pobs%x = pobs%x + eps_pos * pobs%kx
        pobs%y = pobs%y + eps_pos * pobs%ky
        pobs%z = pobs%z + eps_pos * pobs%kz

        !--- Inter-clump medium is vacuum.  The far-side surface is
        !    guaranteed to be in vacuum (clumps are interior to the
        !    sphere), so the photon starts outside any clump.
        pobs%icell_clump = 0

        !--- Gas optical depth per frequency bin.  raytrace_to_edge_tau_gas_clump
        !    treats pobs%xfreq as the reference-frame x; internally it shifts
        !    by the per-clump bulk velocity at every clump entry/exit.
        do kk = 1, observer(i)%nxfreq
          pobs%xfreq = g%xfreq(kk)
          call raytrace_to_edge_tau_gas_clump(pobs, g, tau_gas)
          observer(i)%tau_gas(kk, ix, iy) = tau_gas
        end do

        !--- Frequency-independent column density and dust tau.
        call raytrace_to_edge_column_clump(pobs, g, N_gas, tau_dust)
        observer(i)%N_gas(ix, iy) = N_gas
        if (par%DGR > 0.0_wp) observer(i)%tau_dust(ix, iy) = tau_dust

      end do  ! pixels

      call reduce_mem(observer(i)%tau_gas, shared_memory=.true.)
      call reduce_mem(observer(i)%N_gas,   shared_memory=.true.)
      if (par%DGR > 0.0_wp) call reduce_mem(observer(i)%tau_dust, shared_memory=.true.)
    end do  ! observers

    !--- Write output (rank 0 only).
    if (mpar%p_rank == 0) then
      do i = 1, par%nobs
        if (par%nobs == 1) then
          filename_end = ''
        else
          write(filename_end,'(a,i3.3)') '_', i
        end if
        call write_sightline_tau_outside(trim(par%out_file), g, observer(i), &
                                          suffix=trim(filename_end))
      end do
    end if

    !--- Deallocate.
    do i = 1, par%nobs
      call destroy_mem(observer(i)%tau_gas)
      call destroy_mem(observer(i)%N_gas)
      if (par%DGR > 0.0_wp) call destroy_mem(observer(i)%tau_dust)
    end do

  end subroutine make_sightline_tau_clump

end module sightline_tau_clump_mod
