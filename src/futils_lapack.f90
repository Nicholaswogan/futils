module futils_lapack
  use iso_fortran_env, only: dp => real64
  implicit none
  public

interface
    subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      import :: dp
      integer :: INFO, LDA, LDB, N, NRHS
      integer :: IPIV(*)
      real(dp) :: A(LDA,*), B(LDB,*)
    end subroutine DGESV

    subroutine DGBSV(N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
      import :: dp
      integer :: INFO, KL, KU, LDAB, LDB, N, NRHS
      integer :: IPIV(*)
      real(dp) :: AB(LDAB,*), B(LDB,*)
    end subroutine DGBSV

end interface

end module futils_lapack