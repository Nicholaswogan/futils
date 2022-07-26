module futils_misc
  use iso_fortran_env, only: dp => real64
  implicit none
  public

contains

  !> coppied from fortran stdlib v0.2.0
  elemental function is_close(a, b, tol, abs_tol, equal_nan) result(close)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol, abs_tol
    logical, intent(in), optional :: equal_nan
    logical :: close

    real(dp) :: rel_tol_, abs_tol_
    logical :: equal_nan_

    if (present(tol)) then
      rel_tol_ = tol
    else
      rel_tol_ = 1.0e-5_dp
    endif

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    else
      abs_tol_ = 0.0_dp
    endif

    if (present(equal_nan)) then
      equal_nan_ = equal_nan
    else
      equal_nan_ = .false.
    endif

    if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
        close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
    else
        close = abs(a - b) <= max(abs(rel_tol_*max(abs(a), abs(b))), abs(abs_tol_))
    end if     

  end function

end module