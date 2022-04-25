
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
    real(dp) :: b2, b1_inv, v_new
    
    real(dp) :: b_old0, b_old1
    real(dp) :: b_new0, b_new1
    real(dp) :: v_old
    
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
    
    l = 1
    
    do i = 1,n_new
      b_new0 = new_bins(i)
      b_new1 = new_bins(i+1)
            
      b1_inv = 1.0_dp/(b_new1 - b_new0)
      v_new = 0.0_dp
      do j = l,n_old

        b_old0 = old_bins(j)
        b_old1 = old_bins(j+1)
        v_old = old_vals(j)
        ! Several different cases
        
        !    ____    (old bin)
        !  ________  (new bin)
        if (b_old0 > b_new0 .and. b_old1 < b_new1) then
          
          b2 = b_old1 - b_old0
          v_new = v_new + (b2*b1_inv)*v_old
        
        ! ______       (old bin)
        !     ________ (new bin)
        elseif (b_old0 <= b_new0 .and. &
            b_old1 > b_new0 .and. b_old1 < b_new1) then
          
          b2 = b_old1 - b_new0
          v_new = v_new + (b2*b1_inv)*v_old
          
        !        ____    (old bin)
        !  ________      (new bin)
        elseif (b_old0 >= b_new0 .and. &
                b_old0 < b_new1 .and. b_old1 > b_new1) then
        
          b2 = b_new1 - b_old0
          v_new = v_new + (b2*b1_inv)*v_old
          
          ! need to move to the next new_bin
          l = j
          exit
        
        ! ________    (old bin)
        !   ____      (new bin) 
        elseif (b_old0 <= b_new0 .and. b_old1 >= b_new1) then
          
          v_new = v_new + v_old
          
          ! need to move to the next new_bin
          l = j
          exit
            
        !  There are two remaining cases:
        !         _____   (old bin)
        !   ____          (new bin) 
        ! 
        ! and
        !  ____           (old bin)
        !        _____    (new bin) 
        ! 
        ! These cases can occur on either end of re-binning.
        ! However this is OK. The part of old_bin that does not
        ! overlap with new_bin on the ends, will just be ignored.
          
        endif
      enddo
      new_vals(i) = v_new
    enddo
    
  end subroutine
  
end module