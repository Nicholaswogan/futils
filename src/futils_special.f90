module futils_special
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  ! coppied from fortran stdlib v0.2.0

  public :: legendre, dlegendre
  public :: gauss_legendre

  real(dp), parameter :: pi = acos(-1._dp)
  real(dp), parameter :: tolerance = 4._dp*epsilon(1._dp)
  integer, parameter :: newton_iters = 100

contains

  ! same as scipy.special.roots_legendre
  pure module subroutine gauss_legendre(x, w, interval)
    real(dp), intent(out) :: x(:), w(:)
    real(dp), intent(in), optional :: interval(2)

    associate (n => size(x)-1 )
    select case (n)
      case (0)
        x = 0
        w = 2
      case (1)
        x(1) = -sqrt(1._dp/3._dp)
        x(2) = -x(1)
        w = 1
      case default
        block
        integer :: i,j
        real(dp) :: leg, dleg, delta

        do i = 0, (n+1)/2 - 1
          ! use Gauss-Chebyshev points as an initial guess
          x(i+1) = -cos((2*i+1)/(2._dp*n+2._dp) * pi)
          do j = 1, newton_iters
            leg  = legendre(n+1,x(i+1))
            dleg = dlegendre(n+1,x(i+1))
            delta = -leg/dleg
            x(i+1) = x(i+1) + delta
            if ( abs(delta) <= tolerance * abs(x(i+1)) )  exit
          end do
          x(n-i+1) = -x(i+1)

          dleg = dlegendre(n+1,x(i+1))
          w(i+1)   = 2._dp/((1-x(i+1)**2)*dleg**2) 
          w(n-i+1) = w(i+1)
        end do

        if (mod(n,2) == 0) then
          x(n/2+1) = 0

          dleg = dlegendre(n+1, 0.0_dp)
          w(n/2+1) = 2._dp/(dleg**2) 
        end if
        end block
    end select
    end associate

    if (present(interval)) then
      associate ( a => interval(1) , b => interval(2) )
        x = 0.5_dp*(b-a)*x+0.5_dp*(b+a)
        x(1)       = interval(1)
        x(size(x)) = interval(2)
        w = 0.5_dp*(b-a)*w
      end associate
    end if
  end subroutine

  ! derivatives of legegendre polynomials
  ! unspecified behaviour if n is negative
  pure elemental function dlegendre(n,x) result(dleg)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: dleg

    select case(n)
      case(0)
        dleg = 0
      case(1)
        dleg = 1
      case default
        block
          real(dp) :: leg_down1, leg_down2, leg
          real(dp) :: dleg_down1, dleg_down2
          integer :: i 

          leg_down1  = x
          dleg_down1 = 1

          leg_down2  = 1
          dleg_down2 = 0

          do i = 2, n
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
            dleg = dleg_down2 + (2*i-1)*leg_down1
            leg_down2 = leg_down1
            leg_down1 = leg
            dleg_down2 = dleg_down1
            dleg_down1 = dleg
          end do
        end block
    end select
  end function

  ! legegendre polynomials
  ! unspecified behaviour if n is negative
  pure elemental function legendre(n,x) result(leg)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: leg
    select case(n)
      case(0)
        leg  = 1
      case(1)
        leg  = x
      case default
        block
          real(dp) :: leg_down1, leg_down2
          integer :: i 

          leg_down1  = x
          leg_down2  = 1

          do i = 2, n
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
            leg_down2 = leg_down1
            leg_down1 = leg
          end do
        end block
    end select
  end function

end module
