module utility
   use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
   implicit none
   private
   public :: strupcase, strlowcase
   public :: array_2D_indices, array_3D_indices
   public :: loop_divide
   public :: get_base_name, get_base_input_name, get_extension
   public :: name_for_backup
   public :: is_finite, finite
   public :: copy_file
   public :: time_stamp
   public :: get_date_time
   character(*), private, parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
   character(*), private, parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

   interface finite
      module procedure is_finite32,is_finite64
   end interface finite
   interface is_finite
      module procedure is_finite32,is_finite64
   end interface is_finite
   interface loop_divide
      module procedure loop_divide_int32, loop_divide_int64, &
                       loop_divide4_int32, loop_divide4_int64
   end interface loop_divide
   interface arctan
      module procedure arctan32, arctan64
   end interface arctan
contains

   elemental function arctan32(y,x) result(angle)
      real(kind=real32), intent(in) :: x,y
      real(kind=real32) :: angle
      real(kind=real32), parameter :: twopi   = 6.283185307179586476925286766559005768394_real32
      angle = atan2(y,x)
      if (angle < 0.0_real32) angle = angle + twopi
   end function arctan32
   elemental function arctan64(y,x) result(angle)
      real(kind=real64), intent(in) :: x,y
      real(kind=real64) :: angle
      real(kind=real64), parameter :: twopi   = 6.283185307179586476925286766559005768394_real64
      angle = atan2(y,x)
      if (angle < 0.0_real64) angle = angle + twopi
   end function arctan64
   elemental logical function is_finite32(x)
      use, intrinsic :: ieee_arithmetic
      real(real32),intent(in) :: x
      !is_finite32 = .not. (isnan(x) .or. abs(x) > huge(x))
      !--- in polaris machine, it is not allowed to compare NaN with huge(x).
      !is_finite32 = .not. isnan(x)
      !--- in PGI fortran, isnan is not defined. (2020.11.05)
      !is_finite32 = (x == x)
      !if (is_finite32) is_finite32 = .not. (abs(x) > huge(x))
      !--- the following is the easiest way to do (2020.11.05).
      !--- But, this does not work when gfortran is used with -Ofast option (2022.08.20).
      !is_finite32 = ieee_is_finite(x)
      !--- All the above does not work when gfortran is used with -Ofast option (2022.08.20).
      is_finite32 = (ieee_class(x) /= ieee_quiet_nan) .and.&
                    (ieee_class(x) /= ieee_negative_inf) .and. (ieee_class(x) /= ieee_positive_inf)
   end function is_finite32
   elemental logical function is_finite64(x)
      use, intrinsic :: ieee_arithmetic
      real(real64),intent(in) :: x
      !is_finite64 = .not. (isnan(x) .or. abs(x) > huge(x))
      !--- in polaris machine, it is not allowed to compare NaN with huge(x).
      !is_finite64 = .not. isnan(x)
      !--- in PGI fortran, isnan is not defined. (2020.11.05)
      !is_finite64 = (x == x)
      !if (is_finite64) is_finite64 = .not. (abs(x) > huge(x))
      !--- the following is the easiest way to do (2020.11.05).
      !is_finite64 = ieee_is_finite(x)
      !--- All the above does not work when gfortran is used with -Ofast option (2022.08.20).
      is_finite64 = (ieee_class(x) /= ieee_quiet_nan) .and.&
                    (ieee_class(x) /= ieee_negative_inf) .and. (ieee_class(x) /= ieee_positive_inf)
   end function is_finite64
 
   function strupcase(input_string) result(output_string)
      !--- Argument and result
      character(*), intent(in)     :: input_string
      character(len(input_string)) :: output_string
      !--- Local variables
      integer :: i, n
      !--- Copy input string
      output_string = input_string
      !--- Loop over string elements
      do i = 1, len(output_string)
         !--- Find location of letter in lower case constant string
         n = INDEX(LOWER_CASE, output_string(i:i))
         !--- If current substring is a lower case letter, make it upper case
         if (n /= 0) output_string(i:i) = UPPER_CASE(n:n)
      enddo
   end function strupcase
   function strlowcase(input_string) result(output_string)
      !--- Argument and result
      character(*), intent(in)     :: input_string
      character(len(input_string)) :: output_string
      !--- Local variables
      integer :: i, n
      !--- Copy input string
      output_string = input_string
      !--- Loop over string elements
      do i = 1, len(output_string)
         !--- Find location of letter in upper case constant string
         n = index(UPPER_CASE, output_string(i:i))
         !--- If current substring is an upper case letter, make it lower case
         if (n /= 0) output_string(i:i) = LOWER_CASE(n:n)
      enddo
   end function strlowcase
   !------------------------------------------------
   function get_extension(str) result(res)
   character(len=*)  :: str
   character(len=10) :: ext
   character(len=:), allocatable :: str_tmp, res
   integer :: ii, jj
   str_tmp = trim(str)
   jj  = len_trim(str_tmp)
   ii  = index(str_tmp,'.', back=.true.)+1
   ext = str_tmp(ii:jj)
   if (ext /= 'fits' .and. ext /= 'txt') then
      jj = ii-2
      ii = index(str_tmp(:jj),'.',back=.true.)+1
      if (ii > 1) ext = str_tmp(ii:jj)
   endif
   res = trim(ext)
   end function get_extension
   !------------------------------------------------
   subroutine array_3D_indices(nx,ny,nz,loc,i,j,k)
   !  This converts one-dimensional subscript (loc) of an array into corresponding 3D subscripts (i,j,k).
   implicit none
   integer, intent(in)  :: nx,ny,nz
   integer, intent(in)  :: loc
   integer, intent(out) :: i,j,k
   integer :: nxy
   nxy = nx*ny
   k   = (loc-1)/nxy+1
   j   = (loc-(k-1)*nxy-1)/nx+1
   i   = loc-(j-1)*nx-(k-1)*nxy
   return
   end subroutine array_3D_indices
   !------------------------------------------------
   subroutine array_2D_indices(nx,ny,loc,i,j)
   !  This converts one-dimensional subscript (loc) of an array into corresponding 2D subscripts (i,j).
   implicit none
   integer, intent(in)  :: nx,ny
   integer, intent(in)  :: loc
   integer, intent(out) :: i,j
   j   = (loc-1)/nx+1
   i   = loc-(j-1)*nx
   return
   end subroutine array_2D_indices
   !------------------------------------------------
   subroutine loop_divide4_int64(nloop,nnodes,myid,my_nloop)
   implicit none
   integer(kind=int64), intent(in)  :: nloop
   integer,             intent(in)  :: nnodes, myid
   integer(kind=int64), intent(out) :: my_nloop
   integer(int64) :: remainder

   remainder = mod(nloop, nnodes)
   if (remainder > myid) then
      my_nloop = (nloop/nnodes) - myid + 1
   else
      my_nloop = (nloop/nnodes) - remainder
   endif
   end subroutine loop_divide4_int64
   !------------------------------------------------
   subroutine loop_divide4_int32(nloop,nnodes,myid,my_nloop)
   implicit none
   integer, intent(in)  :: nloop
   integer, intent(in)  :: nnodes, myid
   integer, intent(out) :: my_nloop
   integer :: remainder

   remainder = mod(nloop, nnodes)
   if (remainder > myid) then
      my_nloop = (nloop/nnodes) - myid + 1
   else
      my_nloop = (nloop/nnodes) - remainder
   endif
   end subroutine loop_divide4_int32
   !------------------------------------------------
   subroutine loop_divide_int64(nloop,nthreads,myid,mystart,myend)
   implicit none
   integer(kind=int64), intent(in)  :: nloop
   integer,             intent(in)  :: nthreads, myid
   integer(kind=int64), intent(out) :: mystart, myend
   integer :: remainder

   remainder = mod(nloop, nthreads)
   mystart   = myid * (nloop/nthreads) + 1
   if (remainder > myid) then
      mystart = mystart + myid
      myend   = mystart + (nloop/nthreads)
   else
      mystart = mystart + remainder
      myend   = mystart + (nloop/nthreads) - 1
   endif
   end subroutine loop_divide_int64
   !------------------------------------------------
   subroutine loop_divide_int32(nloop,nthreads,myid,mystart,myend)
   implicit none
   integer, intent(in)  :: nloop
   integer, intent(in)  :: nthreads, myid
   integer, intent(out) :: mystart, myend
   integer :: remainder

   remainder = mod(nloop, nthreads)
   mystart   = myid * (nloop/nthreads) + 1
   if (remainder > myid) then
      mystart = mystart + myid
      myend   = mystart + (nloop/nthreads)
   else
      mystart = mystart + remainder
      myend   = mystart + (nloop/nthreads) - 1
   endif
   end subroutine loop_divide_int32
   !------------------------------------------------
   function get_base_name(filename) result(base_name)
   implicit none
   character(len=*), intent(in) :: filename
   character(len(filename))     :: base_name
   integer :: i1,i2
   i1        = index(filename,'/', back=.true.)+1
   i2        = index(filename,'.fits')-1
   base_name = trim(filename(i1:i2))
   return
   end function get_base_name
   !------------------------------------------------
   function get_base_input_name(filename) result(base_name)
   character(len=*), intent(in) :: filename
   character(len=128) :: base_name
   integer :: i1,i2
   i1        = index(filename,'/', back=.true.)+1
   i2        = index(filename,'.in')-1
   base_name = trim(filename(i1:i2))
   return
   end function get_base_input_name
   !------------------------------------------------
   function name_for_backup(fname) result(fname_next)
   !--- this routine gives a name for backup file for a given base name of a fits file.
   !--- for example, fname = 'DL19' for 'DL19.fits.gz'
   implicit none
   character(len=*), intent(in) :: fname
   character(len=128)           :: filename, fname_next
   logical :: file_exists
   integer :: i

   do i = 0,99
      if (i == 0) then
         fname_next = trim(fname)
      else
         !write(fname_next,'(a,i3.3)') trim(fname)//'_',i
         write(fname_next,'(a,i2.2)') trim(fname)//'_',i
      endif
      filename = trim(fname_next)//'.fits.gz'
      inquire(file=trim(filename), exist=file_exists)
      if (.not.file_exists) exit
   enddo
   end function name_for_backup
   !------------------------------------------------
   subroutine copy_file(infile, outfile, ierr)
   use, intrinsic :: iso_c_binding, only : C_CHAR, C_INT
   implicit none
   character(len=*),  intent(in)  :: infile
   character(len=*),  intent(in)  :: outfile
   integer, optional, intent(out) :: ierr
   integer :: c_err

   interface
      function c_copy_file(c_infile, c_outfile) bind(c, name="copy_file") result(c_err)
      import C_CHAR, C_INT
      character(kind=C_CHAR, len=1), intent(in) :: c_infile(*), c_outfile(*)
      integer(C_INT) :: c_err
      end function c_copy_file
   end interface

   c_err = c_copy_file(str2arr(trim(infile)), str2arr(trim(outfile)))
   if (present(ierr)) ierr = c_err
   end subroutine copy_file
   !------------------------------------------------
   pure function str2arr(string) result (array)
   use, intrinsic :: iso_c_binding, only : C_CHAR, C_NULL_CHAR
   implicit none
   character(len=*),intent(in)  :: string
   character(len=1,kind=C_CHAR) :: array(len(string)+1)
   integer                      :: i

   do i = 1,len_trim(string)
      array(i) = string(i:i)
   enddo
   array(i:i) = C_NULL_CHAR
   end function str2arr
   !------------------------------------------------
   subroutine time_stamp(dtime, reset)
#ifdef MPI
   use mpi
#elif _OPENMP
   use omp_lib
#endif
   implicit none
   real(kind=real64), intent(out) :: dtime
   logical, optional, intent(in)  :: reset
   !--- local variables
   real(kind=real64)  :: time2
   logical            :: time_stamp_reset
   real(real64), save :: time1
   logical,      save :: first_time_stamp__ = .true.

   if (present(reset)) then
      time_stamp_reset = reset
   else
      time_stamp_reset = .false.
   endif

   if (first_time_stamp__ .or. time_stamp_reset) then
#ifdef MPI
      time1 = MPI_WTIME()
#elif _OPENMP
      time1 = omp_get_wtime()
#else
      call cpu_time(time1)
#endif
      first_time_stamp__ = .false.
      dtime              = 0.0_real64
      return
   endif

#ifdef MPI
   time2 = MPI_WTIME()
#elif _OPENMP
   time2 = omp_get_wtime()
#else
   call cpu_time(time2)
#endif
   dtime = time2 - time1
   end subroutine time_stamp
   !------------------------------------------------
   function get_date_time() result(date_time)
   implicit none
   integer,dimension(8) :: values
   character(len=19)    :: date_time
   call date_and_time(VALUES=values)
   write(date_time,'(i4,a,i2.2,a,i2.2,x,i2.2,a,i2.2,a,i2.2)') &
               values(1),'/',values(2),'/',values(3), values(5),':',values(6),':',values(7)
   return
   end function get_date_time
   !------------------------------------------------
end module utility
