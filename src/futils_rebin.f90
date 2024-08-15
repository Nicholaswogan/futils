
module futils_rebin
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: rebin, rebin_with_errors, conserving_rebin

contains

  !> Rebins `old_vals` defined on `old_bins` to `new_bins`. An example is
  !> rebinning a high resolution spectra of infrared emission of Earth 
  !> to a lower resolution. I have optimized the routine for downbinning
  !> data. Upbinning seems like an unlikely application.
  subroutine rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    real(dp), intent(in) :: old_bins(:) !! Edges of bins for which old_vals are defined
    real(dp), intent(in) :: old_vals(:) !! Values defined on old_bins.
    real(dp), intent(in) :: new_bins(:) !! Edges of target bin that you want to rebin to.
    real(dp), intent(out) :: new_vals(:) !! Values defined on new_bins (output).
    integer, optional, intent(out) :: ierr !! Inputs will be checked if ierr
                                           !! is passes as an argument. if ierr < 0
                                           !! on return, then there is an issue with
                                           !! the inputs.
    
    integer :: i, j, l, n_old, n_new
    real(dp) :: b2, b1_inv, v_new
    
    real(dp) :: b_old0, b_old1
    real(dp) :: b_new0, b_new1
    real(dp) :: v_old
    
    n_old = size(old_vals)
    n_new = size(new_vals)
    
    ! option to check inputs.
    if (present(ierr)) then
      call check_rebin_inputs(old_bins, old_vals, new_bins, new_vals, ierr)
      if (ierr /= 0) return
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
            b_old1 > b_new0 .and. b_old1 <= b_new1) then
          
          b2 = b_old1 - b_new0
          v_new = v_new + (b2*b1_inv)*v_old
          
        !        ____    (old bin)
        !  ________      (new bin)
        elseif (b_old0 >= b_new0 .and. &
                b_old0 < b_new1 .and. b_old1 >= b_new1) then
        
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

  !> Rebins `old_vals` and `old_errs` defined on `old_bins` to `new_bins`. This function has
  !> the same behavior as `rebin`, except it also rebins errors.
  subroutine rebin_with_errors(old_bins, old_vals, old_errs, new_bins, new_vals, new_errs, ierr)
    real(dp), intent(in) :: old_bins(:) !! Edges of bins for which old_vals are defined
    real(dp), intent(in) :: old_vals(:) !! Values defined on old_bins.
    real(dp), intent(in) :: old_errs(:) !! Standard deviations defined on old_bins.
    real(dp), intent(in) :: new_bins(:) !! Edges of target bin that you want to rebin to.
    real(dp), intent(out) :: new_vals(:) !! Values defined on new_bins (output).
    real(dp), intent(out) :: new_errs(:) !! Standard deviations defined on new_bins (output).
    integer, optional, intent(out) :: ierr !! Inputs will be checked if ierr
                                           !! is passes as an argument. if ierr < 0
                                           !! on return, then there is an issue with
                                           !! the inputs.
  
    integer :: i, j, l, n_old, n_new
    real(dp) :: b2, b1_inv, v_new, var_new
    
    real(dp) :: b_old0, b_old1
    real(dp) :: b_new0, b_new1
    real(dp) :: v_old, var_old
    
    n_old = size(old_vals)
    n_new = size(new_vals)

    ! option to check inputs.
    if (present(ierr)) then
      call check_rebin_inputs(old_bins, old_vals, new_bins, new_vals, ierr)
      if (ierr /= 0) return

      ! Additionally, check shape of errors
      if (n_old /= size(old_errs)) then
        ierr = -5
        return
      endif

      if (n_new /= size(new_errs)) then
        ierr = -5
        return
      endif

      ! Check to make sure errors are positive
      if (any(old_errs < 0.0_dp)) then
        ierr = -6
        return
      endif

      ! Check to make sure that new bins don't extend beyond old bins
      if (new_bins(1) < old_bins(1) .or. new_bins(n_new+1) > old_bins(n_old+1)) then
        ierr = -7
        return
      endif

    endif
    
    l = 1
    
    do i = 1,n_new
      b_new0 = new_bins(i)
      b_new1 = new_bins(i+1)
            
      b1_inv = 1.0_dp/(b_new1 - b_new0)
      v_new = 0.0_dp
      var_new = 0.0_dp
      do j = l,n_old

        b_old0 = old_bins(j)
        b_old1 = old_bins(j+1)
        v_old = old_vals(j)
        var_old = old_errs(j)**2.0_dp
        ! Several different cases
        
        !    ____    (old bin)
        !  ________  (new bin)
        if (b_old0 > b_new0 .and. b_old1 < b_new1) then
          
          b2 = b_old1 - b_old0
          v_new = v_new + (b2*b1_inv)*v_old
          var_new = var_new + (b2*b1_inv)**2.0_dp*var_old
        
        ! ______       (old bin)
        !     ________ (new bin)
        elseif (b_old0 <= b_new0 .and. &
            b_old1 > b_new0 .and. b_old1 <= b_new1) then
          
          b2 = b_old1 - b_new0
          v_new = v_new + (b2*b1_inv)*v_old
          var_new = var_new + (b2*b1_inv)**2.0_dp*var_old
          
        !        ____    (old bin)
        !  ________      (new bin)
        elseif (b_old0 >= b_new0 .and. &
                b_old0 < b_new1 .and. b_old1 >= b_new1) then
        
          b2 = b_new1 - b_old0
          v_new = v_new + (b2*b1_inv)*v_old
          var_new = var_new + (b2*b1_inv)**2.0_dp*var_old
          
          ! need to move to the next new_bin
          l = j
          exit
        
        ! ________    (old bin)
        !   ____      (new bin) 
        elseif (b_old0 <= b_new0 .and. b_old1 >= b_new1) then
          
          v_new = v_new + v_old
          var_new = var_new + var_old
          
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
      new_errs(i) = sqrt(var_new)
    enddo
    
  end subroutine

  subroutine conserving_rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    real(dp), intent(in) :: old_bins(:) !! Edges of bins for which old_vals are defined
    real(dp), intent(in) :: old_vals(:) !! Values defined on old_bins.
    real(dp), intent(in) :: new_bins(:) !! Edges of target bin that you want to rebin to.
    real(dp), intent(out) :: new_vals(:) !! Values defined on new_bins (output).
    integer, optional, intent(out) :: ierr !! Inputs will be checked if ierr
                                           !! is passes as an argument. if ierr < 0
                                           !! on return, then there is an issue with
                                           !! the inputs.
    integer :: i, l, n_old, n_new
    real(dp) :: val
    
    n_old = size(old_vals)
    n_new = size(new_vals)

    call rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    if (present(ierr)) then
      if (ierr /= 0) return
    endif

    !!!!!!!!!!!!!!!!!
    !!! Left edge !!!
    !!!!!!!!!!!!!!!!!
    !  ______>    (old bins)
    !  ______>    (new bins)
    if (new_bins(1) == old_bins(1)) then
      ! Do nothing. mass is conserved.

    !    ____>    (old bins)
    !  ______>    (new bins)
    elseif (new_bins(1) < old_bins(1)) then
      ! mass is conserved, but lets distribute lowest bin in new_bin
      ! among all the lower bins, to avoid zero values

      do i = 1,n_new
        if (new_bins(i) > old_bins(1)) then
          l = i - 1
          exit
        endif
      enddo
      new_vals(1:l) = new_vals(l)*((new_bins(l+1) - new_bins(l))/(new_bins(l+1) - new_bins(1)))

    !  ______>    (old bins)
    !    ____>    (new bins)
    elseif (new_bins(1) > old_bins(1)) then
      ! mass has been lost. Find the lost mass, and add it to the
      ! first grid cell

      val = 0.0_dp
      do i = 1,n_old
        if (old_bins(i+1) > new_bins(1)) then
          val = val + old_vals(i)*(new_bins(1) - old_bins(i))
          exit
        else
          val = val + old_vals(i)*(old_bins(i+1)-old_bins(i))
        endif
      enddo
      new_vals(1) = new_vals(1) + val/(new_bins(2) - new_bins(1))

    endif

    !!!!!!!!!!!!!!!!!!
    !!! Right edge !!!
    !!!!!!!!!!!!!!!!!!
    !  <______    (old bins)
    !  <______    (new bins)
    if (new_bins(n_new+1) == old_bins(n_old+1)) then
      ! Do nothing

    !  <____      (old bins)
    !  <______    (new bins)
    elseif (new_bins(n_new+1) > old_bins(n_old+1)) then
      ! Mass is conserved, but lets re-distribute mass
      ! so that there are no zero new_bins

      do i = n_new+1,1,-1
        if (new_bins(i) < old_bins(n_old+1)) then
          l = i
          exit
        endif
      enddo
      new_vals(l:) = new_vals(l)*((new_bins(l+1) - new_bins(l))/(new_bins(n_new+1) - new_bins(l)))
    
    !  <______    (old bins)
    !  <____      (new bins)
    elseif (new_bins(n_new+1) < old_bins(n_old+1)) then
      ! Mass is not conserved. Let's add lost mass to furthest right cell in 
      ! new_bins

      val = 0.0_dp
      do i = n_old,1,-1
        if (old_bins(i) < new_bins(n_new+1)) then
          val = val + old_vals(i)*(old_bins(i+1) - new_bins(n_new+1))
          exit
        else
          val = val + old_vals(i)*(old_bins(i+1)-old_bins(i))
        endif
      enddo
      new_vals(n_new) = new_vals(n_new) + val/(new_bins(n_new+1) - new_bins(n_new))
    endif

    if (present(ierr)) then; block
      use futils_misc, only: is_close
      real(dp) :: old_mass, new_mass
      ! check for mass conservation
      old_mass = sum(old_vals*(old_bins(2:)-old_bins(1:size(old_bins)-1)))
      new_mass = sum(new_vals*(new_bins(2:)-new_bins(1:size(new_bins)-1)))

      if (.not. is_close(old_mass, new_mass, tol = 1.0e-10_dp)) then
        ierr = -5
        return
      endif
    end block; endif

  end subroutine

  !> Checks inputs for problems with inputs
  subroutine check_rebin_inputs(old_bins, old_vals, new_bins, new_vals, ierr)
    real(dp), intent(in) :: old_bins(:) !! Edges of bins for which old_vals are defined
    real(dp), intent(in) :: old_vals(:) !! Values defined on old_bins.
    real(dp), intent(in) :: new_bins(:) !! Edges of target bin that you want to rebin to.
    real(dp), intent(in) :: new_vals(:) !! Values defined on new_bins.
    integer, intent(out) :: ierr !! if ierr < 0 on return, then there is an issue with
                                 !! the inputs.

    integer :: i, n_old, n_new

    n_old = size(old_vals)
    n_new = size(new_vals)
    
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
    
    ! check that bins are all increasing
    do i = 1,n_old
      if (old_bins(i+1) <= old_bins(i)) then
        ierr = -3
        return
      endif
    enddo
    
    do i = 1,n_new
      if (new_bins(i+1) <= new_bins(i)) then
        ierr = -4
        return
      endif
    enddo

  end subroutine
  
end module