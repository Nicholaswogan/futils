
module futils_fastmath
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: rebin
  
  !! Routines that will be compiled with -ffast-math
  !! for an extra speed boost.
  
contains

  subroutine rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    real(dp), intent(in) :: old_bins(:)
    real(dp), intent(in) :: old_vals(:)
    real(dp), intent(in) :: new_bins(:)
    real(dp), intent(out) :: new_vals(:)
    integer, optional, intent(out) :: ierr
    
    integer :: i, j, l, n_old, n_new
    real(dp) :: b1, b2
    
    n_old = size(old_vals)
    n_new = size(new_vals)
    
    ! option to check inputs.
    if (present(ierr)) then
      ierr = 0
      
      ! check shape
      if (n_old+1 /= size(old_bins)) then
        ierr = -1
        return
      endif
      
      if (n_new+1 /= size(new_bins)) then
        ierr = -2
        return
      endif
      
      ! check that new bins overlap with the old bins
      if (new_bins(1) < old_bins(1) .and. new_bins(2) < old_bins(1)) then
        ierr = -3
        return
      endif
      
      if (new_bins(n_new) > old_bins(n_old+1) .and. new_bins(n_new+1) > old_bins(n_old+1)) then
        ierr = -4
        return
      endif
      
      ! check that bins are all increasing
      do i = 1,n_old
        if (old_bins(i+1) <= old_bins(i)) then
          ierr = -5
          return
        endif
      enddo
      
      do i = 1,n_new
        if (new_bins(i+1) <= new_bins(i)) then
          ierr = -6
          return
        endif
      enddo
      
    endif
    
    new_vals = 0.0_dp
    l = 1
    
    do i = 1,n_new
      b1 = new_bins(i+1) - new_bins(i)
      do j = l,n_old
        ! Several different cases
        
        !    ____    (old bin)
        !  ________  (new bin)
        if (old_bins(j) > new_bins(i) .and. old_bins(j+1) < new_bins(i+1)) then
          
          b2 = old_bins(j+1) - old_bins(j)
          new_vals(i) = new_vals(i) + (b2/b1)*old_vals(j)
        
        ! ______       (old bin)
        !     ________ (new bin)
        elseif (old_bins(j) <= new_bins(i) .and. &
            old_bins(j+1) > new_bins(i) .and. old_bins(j+1) < new_bins(i+1)) then
          
          b2 = old_bins(j+1) - new_bins(i)
          new_vals(i) = new_vals(i) + (b2/b1)*old_vals(j)
          
        !        ____    (old bin)
        !  ________      (new bin)
        elseif (old_bins(j) >= new_bins(i) .and. &
                old_bins(j) < new_bins(i+1) .and. old_bins(j+1) > new_bins(i+1)) then
        
          b2 = new_bins(i+1) - old_bins(j)
          new_vals(i) = new_vals(i) + (b2/b1)*old_vals(j)
        
        ! ________    (old bin)
        !   ____      (new bin) 
        elseif (old_bins(j) <= new_bins(i) .and. old_bins(j+1) >= new_bins(i+1)) then
          
          new_vals(i) = new_vals(i) + old_vals(j)
          
        !         _____   (old bin)
        !   ____          (new bin) 
        ! elseif (old_bins(j) > new_bins(i+1)) then
        else
          ! time to move onto the next new_bin
          l = j - 1
          exit
          
          ! This is a bit unsafe. If we have the case
          !  ____           (old bin)
          !        _____    (new bin) 
          ! on j == 1, then l = 0, and we will read
          ! outside of the old bins. This is only possible if
          ! ierr is NOT an argument. Error checking determines
          ! if there is some overlap between the bins.
          
        endif
      enddo
    enddo
    
  end subroutine
  
end module