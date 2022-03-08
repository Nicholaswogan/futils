
program test_futils
  use futils
  implicit none
    
  call test_addpnt()
contains
  
  subroutine test_addpnt()
    integer, parameter :: ld = 5
    real(dp) :: x(ld), y(ld)
    real(dp) :: x_new, y_new
    integer :: n, ierr
    
    n = 4
    x = [1.0, 2.0, 3.0, 4.0, 0.0]
    y = [100.0, 90.0, 110.0, 200.0, 0.0]
    
    x_new = 2.5_dp
    y_new = 9999.0_dp
    
    call addpnt(x, y, ld, n, x_new, y_new, ierr)
    
    if (x(3) /= x_new) error stop "test_addpnt: x(3) /= x_new"
    if (y(3) /= y_new) error stop "test_addpnt: y(3) /= y_new"
    if (n /= ld) error stop "test_addpnt: n /= 5"
  
  end subroutine
  
  
end program