
module futils_rebin
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: rebin, rebin_with_errors, conserving_rebin
  public :: make_bins, grid_at_resolution, rebin_error_message

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
        ierr = -6
        return
      endif

      if (n_new /= size(new_errs)) then
        ierr = -7
        return
      endif

      ! Check to make sure errors are positive
      if (any(old_errs < 0.0_dp)) then
        ierr = -8
        return
      endif

      ! Check to make sure that new bins don't extend beyond old bins
      if (new_bins(1) < old_bins(1) .or. new_bins(n_new+1) > old_bins(n_old+1)) then
        ierr = -9
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

  !> Similar to `rebin`, except if the `new_bins` have extents that do not match
  !> `old_bins`, then extra effort is made so that no mass is lost.
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

  !> Checks inputs
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

  !> Given a series of wavelength points, find the corresponding bin edges
  subroutine make_bins(wv, wavl, ierr)
    real(dp), intent(in) :: wv(:)
    real(dp), intent(out) :: wavl(:)
    integer, optional, intent(out) :: ierr

    integer :: i
    
    if (present(ierr)) then
      ierr = 0

      ! Check size
      if (size(wv)+1 /= size(wavl)) then
        ierr = -12
        return
      endif

      ! Must be increasing
      do i = 1,size(wv)-1
        if (wv(i+1) <= wv(i)) then
          ierr = -13
          return
        endif
      enddo

      if (size(wv) <= 1) then
        ierr = -14
        return
      endif

    endif

    wavl(1) = wv(1) - (wv(2) - wv(1))/2.0_dp
    wavl(size(wavl)) = wv(size(wv)) + (wv(size(wv)) - wv(size(wv)-1))/2.0_dp
    wavl(2:size(wavl)-1) = (wv(2:size(wv)) + wv(1:size(wv)-1))/2.0_dp

  end subroutine

  !> Computes a grid of bins at a given resolution `R`
  subroutine grid_at_resolution(wv_min, wv_max, R, wavl, ierr)
    real(dp), intent(in) :: wv_min !! Minimum bin extent
    real(dp), intent(in) :: wv_max !! Maximum bin extent
    real(dp), intent(in) :: R !! Bin resolution (dlam = lam/R)
    real(dp), allocatable, intent(out) :: wavl(:) !! Ouput grid
    integer, optional, intent(out) :: ierr !! If /= 0, then an error has occurred

    integer :: i, n
    real(dp) :: dlam, lam

    if (present(ierr)) then
      ierr = 0

      if (any([wv_min, wv_max, R] <= 0.0_dp)) then
        ierr = -10
        return
      endif

      if (wv_min >= wv_max) then
        ierr = -11
        return
      endif
    endif

    n = 1
    lam = wv_min
    do
      dlam = lam/R
      lam = lam + dlam
      n = n + 1
      if (lam >= wv_max) exit
    enddo

    allocate(wavl(n))
    wavl(1) = wv_min
    lam = wv_min
    do i = 2,n-1
      dlam = lam/R
      lam = lam + dlam
      wavl(i) = lam
    enddo
    wavl(n) = wv_max

  end subroutine

  function rebin_error_message(ierr) result(err)
    integer, intent(in) :: ierr
    character(:), allocatable :: err

    if (ierr >= 0) then
      err = ''
    elseif (ierr == -1) then
      err = '`old_bins` must have a length of `size(old_vals) + 1`'
    elseif (ierr == -2) then
      err = '`new_bins` must have a length of `size(new_vals) + 1`'
    elseif (ierr == -3) then
      err = '`old_bins` must be strictly increasing'
    elseif (ierr == -4) then
      err = '`new_bins` must be strictly increasing'
    elseif (ierr == -5) then
      err = '`conserving_rebin` failed to conserve'
    elseif (ierr == -6) then
      err = '`old_bins` must have a length of `size(old_errs) + 1`'
    elseif (ierr == -7) then
      err = '`new_bins` must have a length of `size(new_errs) + 1`'
    elseif (ierr == -8) then
      err = '`old_errs` can not contain negative values'
    elseif (ierr == -9) then
      err = '`new_bins` must have minimum and maximum extents that fall within `old_bins`'
    elseif (ierr == -10) then
      err = '`wv_min`, `wv_max` and `R` must all be positive'
    elseif (ierr == -11) then
      err = '`wv_max` must be >= `wv_min`'
    elseif (ierr == -12) then
      err = '`size(wv)+1` must equal `size(wavl)`'
    elseif (ierr == -13) then
      err = '`wv` must be strictly increasing'
    elseif (ierr == -14) then
      err = '`wv` must have a size >= 1'
    else
      err = 'Unknown error'
    endif

  end function

end module