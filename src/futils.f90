
module futils
  use iso_fortran_env, only: dp => real64
  use futils_fastmath, only: rebin
  use futils_brent, only: brent_class
  implicit none
  
  private
  
  public :: dp
  ! misc
  public :: Timer, printf, is_close
  ! interpolation
  public :: addpnt, inter2, rebin, interp
  ! strings
  public :: replaceStr
  ! sortting
  public :: argsort, sort
  ! root finding
  public :: brent_class
  
  type Timer
    integer :: cr, cm, c1, c2
    real(dp) :: time
  contains
    procedure :: start => Timer_start
    procedure :: finish => Timer_finish
  end type
  
  interface argsort
    module procedure iargsort, rargsort
  end interface
  
  interface sort
    module procedure sortNums, sortINums
  end interface
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!
  !!! misc utilities !!!
  !!!!!!!!!!!!!!!!!!!!!!
  
  subroutine Timer_start(self)
    class(Timer), intent(inout) :: self
    call system_clock(count = self%c1, count_rate = self%cr, count_max = self%cm)
  end subroutine
  
  subroutine Timer_finish(self, msg, niters)
    class(Timer), intent(inout) :: self
    character(len=*), optional, intent(in) :: msg
    integer, optional, intent(in) :: niters
    call system_clock(count = self%c2)
    
    if (present(niters)) then
      self%time = ((self%c2-self%c1)/real(self%cr))/real(niters)
    else
      self%time = (self%c2-self%c1)/real(self%cr)
    endif
    if (present(msg)) then
      print"(A,1x,es10.4)",trim(msg),self%time
    endif
  end subroutine
  
  subroutine printf(str)
    ! print a string C-style. "\n" indicates a new line.
    use iso_fortran_env, only: output_unit
    character(*), intent(in) :: str
    integer :: i, ii
    ii = 1
    do i = 1,len(str)
      if (str(i:i+1) == '\n') then
        write(output_unit,'(a)') str(ii:i-1)
        ii = i + 2
      endif
    enddo
    if (ii <= len(str)) then
      write(output_unit,'(a)',advance='no') str(ii:)
    endif
  end subroutine
  
  pure elemental function is_close(val1, val2, tol) result(res)
    real(dp), intent(in) :: val1, val2
    real(dp), optional, intent(in) :: tol
    
    logical :: res
    
    real(dp) :: tol_

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = 1.0e-5_dp
    endif
    
    if (val1 < val2 + (val2*tol_) .and. val1 > val2 - (val2*tol_)) then
      res = .true.
    else
      res = .false.
    endif

  end function
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Interpolation and binning !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine addpnt(x, y, ld, n, xnew, ynew, ierr)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!=         program                                                           =*
!=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!=         N < LD.  On exit, N is incremented by 1.                          =*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*
    ! calling parameters
    real(dp), intent(inout) :: x(ld), y(ld)
    integer, intent(in) :: ld
    integer, intent(inout) :: n
    real(dp), intent(in) :: xnew, ynew
    integer, intent(out) :: ierr
    ! local variables
    integer :: insert
    integer :: i
    !-----------------------------------------------------------------------
    ! initialize error flag
    ierr = 0
    ! check n<ld to make sure x will hold another point
    if ( n>=ld ) then
      write (0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
      write (0,*) '                        All elements used.'
      ierr = 1
    endif

    insert = 1
    i = 2
    do
      ! check, whether x is already sorted.
      ! also, use this loop to find the point at which xnew needs to be inserted
      ! into vector x, if x is sorted.
      if ( i<n ) then
        if ( x(i)<x(i-1) ) then
          write (0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//&
                      'in ascending order!' , i , x(i) , x(i-1)
          ierr = 2
        else
          if ( xnew>x(i) ) insert = i + 1
        endif
        i = i + 1
        cycle
      endif
      ! if <xnew,ynew> needs to be appended at the end, just do so,
      ! otherwise, insert <xnew,ynew> at position INSERT
      if ( xnew>x(n) ) then

        x(n+1) = xnew
        y(n+1) = ynew

      else
        ! shift all existing points one index up
        do i = n , insert , -1
          x(i+1) = x(i)
          y(i+1) = y(i)
        enddo
        ! insert new point
        x(insert) = xnew
        y(insert) = ynew
      endif
      ! increase total number of elements in x, y
      n = n + 1
      exit
    enddo

  end subroutine

  subroutine inter2(ng, xg, yg, n, x, y, ierr)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*
    ! input:
    integer, intent(in) :: ng , n
    real(dp), intent(in) :: x(n) , y(n) , xg(ng)
    integer, intent(out) :: ierr
    ! output:
    real(dp), intent(out) :: yg(ng)

    ! local:
    real(dp) :: area , xgl , xgu
    real(dp) :: darea , slope
    real(dp) :: a1 , a2 , b1 , b2
    integer :: ngintv
    integer :: i , k , jstart

    ierr = 0
    !_______________________________________________________________________
    !  test for correct ordering of data, by increasing value of x
    do i = 2 , n
      if ( x(i)<=x(i-1) ) then
        ierr = 1
        write (*,*) 'data not sorted' , i , x(i) , x(i-1)
        return
      endif
    enddo
  
    do i = 2 , ng
      if ( xg(i)<=xg(i-1) ) then
        ierr = 2
        write (0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
        return
      endif
    enddo
    ! check for xg-values outside the x-range

    if ( (x(1)>xg(1)) .or. (x(n)<xg(ng)) ) then
      write (0,*) '>>> ERROR (inter2) <<<  Data do not span '//       &
   &              'grid.  '
      write (0,*) '                        Use ADDPNT to '//          &
   &              'expand data and re-run.'
      return
    endif

    !  find the integral of each grid interval and use this to
    !  calculate the average y value for the interval
    !  xgl and xgu are the lower and upper limits of the grid interval

    jstart = 1
    ngintv = ng - 1
    do i = 1 , ngintv

      ! initalize:
      area = 0.0
      xgl = xg(i)
      xgu = xg(i+1)

      !  discard data before the first grid interval and after the
      !  last grid interval
      !  for internal grid intervals, start calculating area by interpolating
      !  between the last point which lies in the previous interval and the
      !  first point inside the current interval

      k = jstart

      if ( k<=n-1 ) then

        !  if both points are before the first grid, go to the next point
        do while ( x(k+1)<=xgl )
          jstart = k - 1
          k = k + 1
          if ( k>n-1 ) exit
        enddo

        !  if the last point is beyond the end of the grid, complete and go to the next
        !  grid
        do while ( (k<=n-1) .and. (x(k)<xgu) )

          jstart = k - 1

          ! compute x-coordinates of increment

          a1 = max(x(k),xgl)
          a2 = min(x(k+1),xgu)
          !  if points coincide, contribution is zero
          if ( x(k+1)==x(k) ) then
            darea = 0.e0_dp
          else
            slope = (y(k+1)-y(k))/(x(k+1)-x(k))
            b1 = y(k) + slope*(a1-x(k))
            b2 = y(k) + slope*(a2-x(k))
            darea = (a2-a1)*(b2+b1)/2.
          endif
          !  find the area under the trapezoid from a1 to a2
          area = area + darea
          ! go to next point
          k = k + 1
        enddo
      endif
      !  calculate the average y after summing the areas in the interval
      yg(i) = area/(xgu-xgl)
    enddo
    !_______________________________________________________________________

  end subroutine

  ! 1D linear interpolation with constant extrapolation.
  subroutine interp(ng, n, xg, x, y, yg, ierr)
    implicit none
    integer, intent(in) :: ng ! length of new grid (we interpolate to this new grid)
    integer, intent(in) :: n ! length of old grid
    real(dp), intent(in) :: xg(ng) ! new grid
    real(dp), intent(in) :: x(n), y(n) ! old data
    
    real(dp), intent(out) :: yg(ng) ! new data 
    integer, optional, intent(out) :: ierr 
    
    real(dp) :: slope
    integer :: i, j, nn
    
    if (present(ierr)) then
      ierr = 0
      do i = 1,n-1
        if (x(i+1) <= x(i)) then
          ierr = -1
          return
        endif
      enddo
      do i = 1,ng-1
        if (xg(i+1) <= xg(i)) then
          ierr = -2
          return
        endif
      enddo
    endif
    
    nn = 1
    do i = 1,ng
      if (xg(i) <= x(1)) then
        ! contant extrapolation below
        yg(i) = y(1)
      elseif ((xg(i) > x(1)) .and. (xg(i) < x(n))) then
        ! if within the data points, then linear interpolation
        do j = nn,n
          if ((xg(i) >= x(j)) .and. (xg(i) <= x(j+1))) then
            slope = (y(j+1)-y(j))/(x(j+1)-x(j))
            yg(i) = y(j) + slope*(xg(i)-x(j))
            nn = j
            exit
          endif
        enddo
      elseif (xg(i) >= x(n)) then
        ! contant extrapolation above
        yg(i) = y(n)
      endif
    enddo
    
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!! String Utilities !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  pure recursive function replaceStr(string,search,substitute) result(modifiedString)
    character(len=*), intent(in)  :: string, search, substitute
    character(len=:), allocatable :: modifiedString
    integer                       :: i, stringLen, searchLen
    stringLen = len(string)
    searchLen = len(search)
    if (stringLen==0 .or. searchLen==0) then
      modifiedString = ""
      return
    elseif (stringLen<searchLen) then
      modifiedString = string
      return
    end if
    i = 1
    do
      if (string(i:i+searchLen-1)==search) then
        modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
        exit
      end if
      if (i+searchLen>stringLen) then
        modifiedString = string
        exit
      end if
      i = i + 1
      cycle
    end do
  end function replaceStr
  
  !!!!!!!!!!!!!!!
  !!! Sorting !!!
  !!!!!!!!!!!!!!!

  subroutine sortNums(nums)
    ! sorts array of numbers, nums, from smallest to largest
    real(dp), intent(inout):: nums(:)   ! array of numbers
    nums = nums(argsort(nums))
  end subroutine

  subroutine sortINums(nums)
    ! sorts array of inegers, nums, from smallest to largest
    integer, intent(inout):: nums(:)    ! array of numbers
    nums = nums(argsort(nums))
  end subroutine

  function iargsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    integer, intent(in):: a(:)    ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! iargsort([10, 9, 8, 7, 6])   ! Returns [5, 4, 3, 2, 1]

    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp                        ! temporary
    integer :: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1

        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
            temp = b(i); b(i) = b(imin); b(imin) = temp
        end if
    end do
  end function

  function rargsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    real(dp), intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    real(dp) :: temp2
    real(dp) :: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1
        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
        end if
    end do
  end function
  
end module