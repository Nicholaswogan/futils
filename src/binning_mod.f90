subroutine addpnt(x,y,ld,n,xnew,ynew,ierr)
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
  implicit none
  ! calling parameters
  integer, intent(in) :: ld
  integer, intent(inout) :: n
  real(dp), intent(inout) :: x(ld), y(ld)
  real(dp), intent(in) :: xnew, ynew
  integer, intent(out) :: ierr
  ! local variables
  integer insert
  integer i
  !-----------------------------------------------------------------------
  ! initialize error flag
  ierr = 0
  ! check n<ld to make sure x will hold another point
  if (n >= ld) then
    ! write(*,*) '>>> error (addpnt) <<<  cannot expand array '
    ! write(*,*) '                        all elements used.'
    ierr = 1
    return
  endif
 
  insert = 1
 
  ! check, whether x is already sorted.
  ! also, use this loop to find the point at which xnew needs to be inserted
  ! into vector x, if x is sorted.
  do i = 2, n - 1
    if ( x(i).lt.x(i-1) ) then
      ! write(*,*) '>>> error (addpnt) <<<  x-data must be '//     &
      !           &'in ascending order!' , i , x(i) , x(i-1)
      ierr = 2
      return
    else
      if ( xnew.gt.x(i) ) insert = i + 1
    endif
  enddo

  ! if <xnew,ynew> needs to be appended at the end, just do so,
  ! otherwise, insert <xnew,ynew> at position INSERT
 
  if ( xnew.gt.x(n) ) then
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
 
end subroutine


subroutine inter2_1(ng,xg,yg,n,x,y,ierr)

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
  implicit none
  ! input:
  integer, intent(in) :: ng, n
  real(dp), intent(in) :: x(n), y(n), xg(ng)
  ! output:
  real(dp), intent(out) :: yg(ng)
  integer, intent(out) :: ierr
  ! local:
  real(dp) :: area, xgl, xgu
  real(dp) :: darea, slope
  real(dp) :: a1, a2, b1, b2
  integer :: ngintv
  integer :: i, k, jstart

  ierr = 0
  !_______________________________________________________________________
  !  test for correct ordering of data, by increasing value of x
  do i = 2 , n
    if ( x(i).le.x(i-1) ) then
      ierr = 1
      ! write (*,*) 'data not sorted' , i , x(i) , x(i-1)
      return
    endif
  enddo

  do i = 2 , ng
    if ( xg(i).le.xg(i-1) ) then
      ierr = 2
      ! write (0,*) '>>> error (inter2) <<<  xg-grid not sorted!'
      return
    endif
  enddo

  ! check for xg-values outside the x-range
  if ( (x(1).gt.xg(1)) .or. (x(n).lt.xg(ng)) ) then
    ! write (0,*) '>>> error (inter2) <<<  data do not span '//      &
    !           &'grid.  '
    ! write (0,*) '                        use addpnt to '//         &
    !           &'expand data and re-run.'
    ierr = 3
    return
  endif

  !  find the integral of each grid interval and use this to
  !  calculate the average y value for the interval
  !  xgl and xgu are the lower and upper limits of the grid interval
  jstart = 1
  ngintv = Ng - 1
  do i = 1 , ngintv
    ! initalize:
    area = 0.0
    xgl = Xg(i)
    xgu = Xg(i+1)
    ! discard data before the first grid interval and after the
    ! last grid interval
    ! for internal grid intervals, start calculating area by interpolating
    ! between the last point which lies in the previous interval and the
    ! first point inside the current interval
    k = jstart
    if ( k <= n-1 ) then

      ! if both points are before the first grid, go to the next point

  
20         IF ( X(k+1).LE.xgl ) THEN
         jstart = k - 1
         k = k + 1
         IF ( k.LE.N-1 ) GOTO 20
      ENDIF

!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
40         IF ( (k.LE.N-1) .AND. (X(k).LT.xgu) ) THEN

         jstart = k - 1

! compute x-coordinates of increment

         a1 = MAX(X(k),xgl)
         a2 = MIN(X(k+1),xgu)

!  if points coincide, contribution is zero

         IF ( X(k+1).EQ.X(k) ) THEN
            darea = 0.E0
         ELSE
            slope = (Y(k+1)-Y(k))/(X(k+1)-X(k))
            b1 = Y(k) + slope*(a1-X(k))
            b2 = Y(k) + slope*(a2-X(k))
            darea = (a2-a1)*(b2+b1)/2.
!                       print *,a2,a1,k,y(k),slope,b2,b1,darea
         ENDIF


!  find the area under the trapezoid from a1 to a2

         area = area + darea

! go to next point

         k = k + 1
         GOTO 40

      ENDIF

   ENDIF

!  calculate the average y after summing the areas in the interval
   Yg(i) = area/(xgu-xgl)


ENDDO
!_______________________________________________________________________

END
