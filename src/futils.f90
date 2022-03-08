
module futils
  use iso_fortran_env, only: dp => real64
  implicit none
  
  type Timer
    integer :: cr, cm, c1, c2
    real(dp) :: time
  contains
    procedure :: start => Timer_start
    procedure :: finish => Timer_finish
  end type
  
  ! Sorting
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
  
  subroutine Timer_finish(self, msg)
    class(Timer), intent(inout) :: self
    character(len=*), optional, intent(in) :: msg
    call system_clock(count = self%c2)
    self%time = (self%c2-self%c1)/real(self%cr)
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Interpolation and binning !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! binning.f90 is same as binning.f 
  ! adding to module means explicit interface to function
  ! which allows for more compiler errors to be caught
  include "binning.f90"

  ! 1D linear interpolation with constant extrapolation.
  subroutine interp(ng, n, xg, x, y, yg, err)
    implicit none
    integer, intent(in) :: ng ! length of new grid (we interpolate to this new grid)
    integer, intent(in) :: n ! length of old grid
    real(dp), intent(in) :: xg(ng) ! new grid
    real(dp), intent(in) :: x(n), y(n) ! old data
    
    real(dp), intent(out) :: yg(ng) ! new data 
    character(:), allocatable, intent(out) :: err 
    
    real(dp) :: slope
    integer :: i, j, nn
    
    do i = 1,n-1
      if (x(i+1) <= x(i)) then
        err = 'x must be sorted.'
        return
      endif
    enddo
    do i = 1,ng-1
      if (xg(i+1) <= xg(i)) then
        err = 'xg must be sorted.'
        return
      endif
    enddo
    
    nn = 1
    do i = 1,ng
      if (xg(i) <= x(1)) then
        yg(i) = y(1)
      elseif ((xg(i) > x(1)) .and. (xg(i) < x(n))) then
        do j = nn,n
          if ((xg(i) >= x(j)) .and. (xg(i) <= x(j+1))) then
            slope = (y(j+1)-y(j))/(x(j+1)-x(j))
            yg(i) = y(j) + slope*(xg(i)-x(j))
            nn = j
            exit
          endif
        enddo
      elseif (xg(i) >= x(n)) then
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!! Sorting !!!
  !!!!!!!!!!!!!!!!!!!!!!!!

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