module read_text_data
  use define
  use mathlib
contains
  !=================================================================
  !--- for par%geometry == 'plane_atmosphere'
  subroutine read_plane_data(fname,array,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,xnew,ynew,x2
  integer :: unit, info
  integer :: n, i

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  !$OMP parallel do &
  !$OMP default(shared) &
  !$OMP private(i,xnew,ynew)
  do i=1, grid%nz
     xnew = (grid%zface(i)+grid%zface(i+1))/2.0
     if (xnew >= x1 .and. xnew <= x2) then
        call interp(xarr,yarr,xnew,ynew)
     endif
     array(:,:,i) = ynew
  enddo
  !$OMP end parallel do
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_plane_data
  !=================================================================
  !--- for par%geometry == 'spherical_atmosphere'
  subroutine read_spherical_data(fname,array,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,x,y,z,r,ynew,x2
  integer :: unit, info
  integer :: n, i, j, k

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  !$OMP parallel do &
  !$OMP default(shared) &
  !$OMP private(i,j,k,x,y,z,r,ynew)
  do k=1, grid%nz
  do j=1, grid%ny
  do i=1, grid%nx
     x = (grid%xface(i)+grid%xface(i+1))/2.0
     y = (grid%yface(j)+grid%yface(j+1))/2.0
     z = (grid%zface(k)+grid%zface(k+1))/2.0
     r = sqrt(x**2 + y**2 + z**2)
     if (r >= x1 .and. r <= x2) then
        call interp(xarr,yarr,r,ynew)
        array(i,j,k) = ynew
     endif
  enddo
  enddo
  enddo
  !$OMP end parallel do
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_spherical_data
  !=================================================================
  !--- for par%geometry == 'spherical_atmosphere'
  subroutine read_spherical_velocity(fname,vx,vy,vz,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,x,y,z,r,ynew,x2
  integer :: unit, info
  integer :: n, i, j, k

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  !$OMP parallel do &
  !$OMP default(shared) &
  !$OMP private(i,j,k,x,y,z,r,ynew)
  do k=1, grid%nz
  do j=1, grid%ny
  do i=1, grid%nx
     x = (grid%xface(i)+grid%xface(i+1))/2.0
     y = (grid%yface(j)+grid%yface(j+1))/2.0
     z = (grid%zface(k)+grid%zface(k+1))/2.0
     r = sqrt(x**2 + y**2 + z**2)
     call interp(xarr,yarr,r,ynew)
     if (r > 0.0_wp .and. r >= x1 .and. r <= x2) then
        vx(i,j,k) = x/r * ynew
        vy(i,j,k) = y/r * ynew
        vz(i,j,k) = z/r * ynew
     endif
  enddo
  enddo
  enddo
  !$OMP end parallel do
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_spherical_velocity
  !=================================================================
  subroutine setup_plane_emissivity(fname,emiss,grid,sampling_method,f_composite)
  use random
  use memory_mod
  implicit none
  character(len=*),             intent(in)    :: fname
  type(emiss_1D_profile_type),  intent(inout) :: emiss
  type(grid_type),              intent(in)    :: grid
  integer,       optional, intent(in) :: sampling_method
  real(kind=wp), optional, intent(in) :: f_composite

  !-- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1, y1, xmax, f_comp, f_comp1
  integer       :: unit, info, i, sample_method
  real(kind=wp) :: Psum, Pcomp

  xmax = minval([grid%xmax, grid%ymax, grid%zmax])

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  emiss%nPDF = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     if (y1 < 0.0_wp) y1 = 0.0_wp
     emiss%nPDF = emiss%nPDF + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
     if (x1 >= xmax) exit
  enddo
  close(unit)

  call create_mem(emiss%axis,       [emiss%nPDF])
  call create_mem(emiss%prob,       [emiss%nPDF])
  call create_mem(emiss%prob_alias, [emiss%nPDF-1])
  call create_mem(emiss%alias,      [emiss%nPDF-1])

  do i=1,emiss%nPDF
     emiss%axis(i) = xarr(i)
     emiss%prob(i) = yarr(i)
  enddo
  if (xarr(emiss%nPDF) > xmax) then
     emiss%axis(emiss%nPDF) = xmax
     call interp(xarr,yarr,xmax,emiss%prob(emiss%nPDF))
  endif
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)

  do i=1,emiss%nPDF-1
     !-- we need to deal with the case where x-axis values are not equally spaced.
     emiss%prob_alias(i) = (emiss%prob(i) + emiss%prob(i+1))/2.0_wp * (emiss%axis(i+1) - emiss%axis(i))
  enddo
  Psum             = sum(emiss%prob_alias)
  emiss%prob_alias = emiss%prob_alias / Psum
  emiss%prob       = emiss%prob       / Psum

  !-- sampling method = 0: use the original PDF.
  !--                 = 1: use the composite method
  sample_method = 0
  f_comp        = 0.5_wp
  if (present(sampling_method)) then
     sample_method = sampling_method
     if (present(f_composite)) f_comp = f_composite
  endif

  !-- composite method
  if (sample_method > 0) then
     call create_mem(emiss%wgt, [emiss%nPDF])
     !-- we need to deal with the case where x-axis values are not equally spaced.
     !-- prob is the probability density at a given single point.
     !-- prob_alias is the probability for each bin.
     f_comp1 = 1.0_wp - f_comp
     Psum    = 0.0_wp
     do i=1,emiss%nPDF-1
        if (emiss%prob_alias(i) > 0.0_wp) then
           Psum = Psum + emiss%axis(i+1) - emiss%axis(i)
        endif
     enddo
     do i=1,emiss%nPDF-1
        if (emiss%prob_alias(i) > 0.0_wp) then
           Pcomp               = (emiss%axis(i+1) - emiss%axis(i))/Psum
           emiss%prob_alias(i) = emiss%prob_alias(i) * f_comp1 + f_comp * Pcomp
        endif
     enddo
     emiss%wgt(:)  = emiss%prob(:) / (emiss%prob(:) * f_comp1 + f_comp / Psum)
     emiss%prob(:) = emiss%prob(:) * f_comp1 + f_comp / Psum
  endif

  call random_alias_setup(emiss%prob_alias, emiss%alias)
  end subroutine setup_plane_emissivity
  !=================================================================
  subroutine setup_spherical_emissivity(fname,emiss,grid,sampling_method,f_composite)
  use random
  use memory_mod
  implicit none
  character(len=*),             intent(in)    :: fname
  type(emiss_1D_profile_type),  intent(inout) :: emiss
  type(grid_type),              intent(in)    :: grid
  integer,       optional, intent(in) :: sampling_method
  real(kind=wp), optional, intent(in) :: f_composite

  !-- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1, y1, xmax, f_comp, f_comp1
  integer       :: unit, info, i, sample_method
  real(kind=wp) :: Psum, Pcomp

  xmax = minval([grid%xmax, grid%ymax, grid%zmax])

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  emiss%nPDF = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     if (y1 < 0.0_wp) y1 = 0.0_wp
     emiss%nPDF = emiss%nPDF + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
     if (x1 >= xmax) exit
  enddo
  close(unit)

  call create_mem(emiss%axis,       [emiss%nPDF])
  call create_mem(emiss%prob,       [emiss%nPDF])
  call create_mem(emiss%prob_alias, [emiss%nPDF-1])
  call create_mem(emiss%alias,      [emiss%nPDF-1])

  yarr(:) = yarr(:) * xarr(:)**2
  do i=1,emiss%nPDF
     emiss%axis(i) = xarr(i)
     emiss%prob(i) = yarr(i)
  enddo
  if (xarr(emiss%nPDF) > xmax) then
     emiss%axis(emiss%nPDF) = xmax
     call interp(xarr,yarr,xmax,emiss%prob(emiss%nPDF))
  endif
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)

  do i=1,emiss%nPDF-1
     !-- we need to deal with the case where x-axis values are not equally spaced.
     emiss%prob_alias(i) = (emiss%prob(i) + emiss%prob(i+1))/2.0_wp * (emiss%axis(i+1) - emiss%axis(i))
  enddo
  Psum             = sum(emiss%prob_alias)
  emiss%prob_alias = emiss%prob_alias / Psum
  emiss%prob       = emiss%prob       / Psum

  !-- sampling method = 0: use the original PDF.
  !--                 = 1: use the composite method
  sample_method = 0
  f_comp        = 0.5_wp
  if (present(sampling_method)) then
     sample_method = sampling_method
     if (present(f_composite)) f_comp = f_composite
  endif

  !-- composite method
  if (sample_method > 0) then
     call create_mem(emiss%wgt, [emiss%nPDF])
     !-- we need to deal with the case where x-axis values are not equally spaced.
     !-- prob is the probability density at a given single point.
     !-- prob_alias is the probability for each bin.
     f_comp1 = 1.0_wp - f_comp
     Psum    = 0.0_wp
     do i=1,emiss%nPDF-1
        if (emiss%prob_alias(i) > 0.0_wp) then
           Psum = Psum + emiss%axis(i+1) - emiss%axis(i)
        endif
     enddo
     do i=1,emiss%nPDF-1
        if (emiss%prob_alias(i) > 0.0_wp) then
           Pcomp               = (emiss%axis(i+1) - emiss%axis(i))/Psum
           emiss%prob_alias(i) = emiss%prob_alias(i) * f_comp1 + f_comp * Pcomp
        endif
     enddo
     emiss%wgt(:)  = emiss%prob(:) / (emiss%prob(:) * f_comp1 + f_comp / Psum)
     emiss%prob(:) = emiss%prob(:) * f_comp1 + f_comp / Psum
  endif

  call random_alias_setup(emiss%prob_alias, emiss%alias)
  end subroutine setup_spherical_emissivity
  !=================================================================
  subroutine read_stars(fname, star, sampling_method, f_composite)
  use mpi
  use random
  use memory_mod
  implicit none
  character(*),    intent(in)  :: fname
  type(star_type), intent(out) :: star
  integer,       optional, intent(in) :: sampling_method
  real(kind=wp), optional, intent(in) :: f_composite

  !-- local variables
  real(kind=wp) :: x1, y1, z1, lum1, f_comp
  integer       :: unit, info, ierr, i, sample_method, ncount

  open(newunit=unit,file=trim(fname),status='old')
  star%nstars = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1,z1,lum1
     if (info /= 0) exit
     star%nstars = star%nstars + 1
  enddo
  close(unit)

  call create_mem(star%x,     [star%nstars])
  call create_mem(star%y,     [star%nstars])
  call create_mem(star%z,     [star%nstars])
  call create_mem(star%lum,   [star%nstars])
  call create_mem(star%prob,  [star%nstars])
  call create_mem(star%alias, [star%nstars])

  open(newunit=unit,file=trim(fname),status='old')
  do i=1, star%nstars
     read(unit,*,iostat=info) x1,y1,z1,lum1
     if (lum1 < 0.0_wp) lum1 = 0.0_wp
     star%x(i)   = x1
     star%y(i)   = y1
     star%z(i)   = z1
     star%lum(i) = lum1
  enddo
  close(unit)
  star%prob(:) = star%lum(:) / sum(star%lum)

  !-- sampling method = 0: use the original PDF.
  !--                 = 1: use the composite method
  sample_method = 0
  f_comp        = 0.5_wp
  if (present(sampling_method)) then
     sample_method = sampling_method
     if (present(f_composite)) f_comp = f_composite
  endif

  if (sample_method > 0) then
     call create_mem(star%wgt,  [star%nstars])
     ncount = count(star%prob > 0.0_wp)
     do i=1, star%nstars
        if (star%prob(i) > 0.0_wp) then
           star%wgt(i)  = star%prob(i) / (star%prob(i) * (1.0_wp - f_comp) + f_comp / ncount)
           star%prob(i) = star%prob(i) * (1.0_wp - f_comp) + f_comp / ncount
        endif
     enddo
     call random_alias_setup(star%prob, star%alias)
  endif
  end subroutine read_stars
  !=================================================================
end module read_text_data
