
program test_futils
  use futils
  implicit none

  type, extends(brent_class) :: myfunc_type
    integer :: i = 0 !! function counter
  end type

  call test_sort()
  call test_searchsorted()
  call test_interp()
  call test_conserving_rebin()
  call test_expi()
  call test_gauss_legendre()
  call test_addpnt()
  call test_rebin_with_errors()
  call test_rebin()
  call test_grid_at_resolution()
  call test_make_bins()
  call test_brent()

contains

  subroutine test_sort()
    real(dp), allocatable :: arr(:), arr_sorted(:)
    integer, allocatable :: iarr(:), iarr_sorted(:)
    
    arr = [1.0_dp, 0.0_dp, 6.0_dp, 3.0_dp, 4.0_dp]
    arr_sorted = [0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp, 6.0_dp]
    call sort(arr)
    if (.not. all(is_close(arr, arr_sorted))) then
      error stop "test_sort test failed"
    endif

    iarr = [1, 0, 6, 3, 4, 10, 2]
    iarr_sorted = [0, 1, 2, 3, 4, 6, 10]
    call sort(iarr)
    if (.not. all(iarr == iarr_sorted)) then
      error stop "test_sort test failed"
    endif

  end subroutine

  subroutine test_searchsorted()
    real(dp), allocatable :: arr(:)
    real(dp) :: val
    integer :: ind

    arr = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    val = 5.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 6) then
      error stop "test_searchsorted error"
    endif

    val = 3.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 4) then
      error stop "test_searchsorted error"
    endif

    val = 3.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 3) then
      error stop "test_searchsorted error"
    endif

    val = 6.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 6) then
      error stop "test_searchsorted error"
    endif

    val = 7.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 7) then
      error stop "test_searchsorted error"
    endif

    val = 1.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 0.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 1.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 2) then
      error stop "test_searchsorted error"
    endif

    arr = [1.0_dp]
    val = 1.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 0.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 2.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 2) then
      error stop "test_searchsorted error"
    endif

    arr = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp]
    val = 2.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 2) then
      error stop "test_searchsorted error"
    endif

    val = 1.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 2) then
      error stop "test_searchsorted error"
    endif

    val = 1.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 0.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 1) then
      error stop "test_searchsorted error"
    endif

    val = 7.1_dp
    ind = searchsorted(arr, val)
    if (ind /= 8) then
      error stop "test_searchsorted error"
    endif

    val = 7.0_dp
    ind = searchsorted(arr, val)
    if (ind /= 7) then
      error stop "test_searchsorted error"
    endif

    val = 6.9_dp
    ind = searchsorted(arr, val)
    if (ind /= 7) then
      error stop "test_searchsorted error"
    endif

  end subroutine

  subroutine test_interp()
    real(dp), allocatable :: xg(:), yg(:), x(:), y(:), yg_scipy(:)
    integer :: ierr

    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    y = [2.0_dp, 3.0_dp, 2.8_dp, 3.2_dp, 4.6_dp]
    xg = [0.0_dp, 4.0_dp, 6.0_dp]
    allocate(yg(size(xg)))
    yg_scipy = [1.0_dp , 3.2_dp, 6.0_dp]

    call interp(xg, x, y, yg, linear_extrap=.true., ierr=ierr)
    if (ierr /= 0) then
      print*,ierr
      error stop "interp: returned error"
    endif
    if (.not. all(is_close(yg, yg_scipy))) then
      error stop "interp: returned incorrect values"
    endif

  end subroutine

  subroutine test_expi()
    real(dp), parameter :: scipy_result = 1665628.0229506120_dp
    if (.not.is_close(scipy_result, expi(17.1_dp),tol=1.0e-12_dp)) then
      error stop "expi does not match scipy"
    endif
  end subroutine

  subroutine test_gauss_legendre()
    real(dp) :: x(4), w(4)
    ! result from using scipy.special.roots_legendre
    real(dp), parameter :: x_scipy(4) = [-0.8611363115940526_dp, -0.3399810435848563_dp, &
                                          0.3399810435848563_dp, 0.8611363115940526_dp]
    real(dp), parameter :: w_scipy(4) = [0.3478548451374536_dp, 0.6521451548625464_dp, &
                                         0.6521451548625464_dp, 0.3478548451374536_dp]

    call gauss_legendre(x, w)

    if (.not. all(is_close(x, x_scipy))) then
      print*,x/x_scipy
      error stop 'gauss_legendre does not match scipy'
    endif
    if (.not. all(is_close(w, w_scipy))) then
      print*,w/w_scipy
      error stop 'gauss_legendre does not match scipy'
    endif

  end subroutine
  
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
        
    if (.not. is_close(x(3), x_new)) error stop "test_addpnt: x(3) /= x_new"
    if (.not. is_close(y(3), y_new)) error stop "test_addpnt: y(3) /= y_new"
    if (n /= ld) error stop "test_addpnt: n /= 5"
  
  end subroutine

  subroutine test_rebin_with_errors()

    integer :: ierr
    real(dp), allocatable :: flux_new(:), err_new(:)

    real(dp), parameter :: wavl(*) = &
    [ &
      1.0000000000000000e+00_dp, 1.2857142857142856e+00_dp, 1.5714285714285714e+00_dp, 1.8571428571428572e+00_dp, &
      2.1428571428571428e+00_dp, 2.4285714285714284e+00_dp, 2.7142857142857144e+00_dp, 3.0000000000000000e+00_dp, &
      3.2857142857142856e+00_dp, 3.5714285714285712e+00_dp, 3.8571428571428568e+00_dp, 4.1428571428571423e+00_dp, &
      4.4285714285714288e+00_dp, 4.7142857142857135e+00_dp, 5.0000000000000000e+00_dp &
    ]
    real(dp), parameter :: flux(*) = &
    [ &
      9.0982291294112394e-01_dp, 9.8990307637212382e-01_dp, 9.8972304885982143e-01_dp, 9.0929742682568171e-01_dp, &
      7.5514702623165808e-01_dp, 5.3977018240064412e-01_dp, 2.8062939951435684e-01_dp, -1.2644889303772900e-03_dp, &
      -2.8305585408225559e-01_dp, -5.4189737955940476e-01_dp, -7.5680249530792798e-01_dp, -9.1034694431078278e-01_dp, &
      -9.9008152109583547e-01_dp, -9.8954143884678003e-01_dp &
    ]
    real(dp), parameter :: err(*) = &
    [ &
      7.7440675196366238e-01_dp, 8.5759468318620979e-01_dp, 8.0138168803582199e-01_dp, 7.7244159149844838e-01_dp, &
      7.1182739966945241e-01_dp, 8.2294705653332811e-01_dp, 7.1879360563134620e-01_dp, 9.4588650039103994e-01_dp, &
      9.8183138025051464e-01_dp, 6.9172075941288891e-01_dp, 8.9586251904133229e-01_dp, 7.6444745987645224e-01_dp, &
      7.8402228054696610e-01_dp, 9.6279831914633052e-01_dp &
    ]
    real(dp), parameter :: wavl_new(*) = &
    [ &
      1.1000000000000001e+00_dp, 2.0000000000000000e+00_dp, 2.8999999999999999e+00_dp, 3.7999999999999998e+00_dp, &
      4.6999999999999993e+00_dp, 5.0000000000000000e+00_dp &
    ]
    ! Correct values were computed with https://github.com/ACCarnall/spectres
    real(dp), parameter :: flux_new_correct(*) = &
    [ &
      9.6052689906556088e-01_dp, 6.1332540689824611e-01_dp, -1.9670411327940693e-01_dp, -8.6225598748305898e-01_dp, &
      -9.8956715704911591e-01_dp &
    ]
    real(dp), parameter :: err_new_correct(*) = &
    [ &
      4.2357147603948003e-01_dp, 3.9541301062417333e-01_dp, 4.7387917075490427e-01_dp, 4.4454049374724319e-01_dp, &
      9.1771051546596749e-01_dp &
    ]

    allocate(flux_new(size(wavl_new)-1), err_new(size(wavl_new)-1))

    call rebin_with_errors(wavl, flux, err, wavl_new, flux_new, err_new, ierr)
    if (ierr /= 0) error stop "rebin: ierr not zero."

    if (.not. all(is_close(flux_new, flux_new_correct, tol = 1.0e-10_dp))) then
      error stop "rebin_with_errors: returned incorrect values"
    endif

    if (.not. all(is_close(err_new, err_new_correct, tol = 1.0e-10_dp))) then
      error stop "rebin_with_errors: returned incorrect errors"
    endif

  end subroutine
  
  subroutine test_rebin()
      
    real(dp), allocatable :: old_bins(:)
    real(dp), allocatable :: old_vals(:)
    real(dp), allocatable :: new_bins(:)
    real(dp), allocatable :: new_vals(:)
    
    real(dp), allocatable :: correct_vals(:)
    
    integer :: i, ierr
      
    old_bins = &
      [0.0000000000000000e+00_dp, 2.7301300317049026e-02_dp, 7.8484747558832169e-02_dp, &
       1.2966819480061531e-01_dp, 2.2562499716877937e-01_dp, 2.7680844441056252e-01_dp, &
       3.7276524677872658e-01_dp, 4.2394869402050972e-01_dp, 5.1990549638867378e-01_dp, &
       6.1586229875683784e-01_dp, 6.4316359907388687e-01_dp, 6.9434704631567001e-01_dp, &
       7.4553049355745316e-01_dp, 7.4696740403305739e-01_dp, 7.4966126971412450e-01_dp, &
       7.5235513539519161e-01_dp, 7.5504900107625872e-01_dp, 7.6009935897309333e-01_dp, &
       7.6514971686992794e-01_dp, 7.9245101718697697e-01_dp, 8.4363446442876011e-01_dp, &
       8.9481791167054325e-01_dp, 9.2211921198759228e-01_dp, 9.2355612246319652e-01_dp, &
       9.2624998814426363e-01_dp, 9.2894385382533073e-01_dp, 9.3399421172216535e-01_dp, &
       9.3904456961899996e-01_dp, 9.4173843530006707e-01_dp, 9.4317534577567130e-01_dp, &
       9.4586921145673841e-01_dp, 9.4856307713780552e-01_dp, 9.4999998761340976e-01_dp, &
       9.5007561447710032e-01_dp, 9.5021739687217632e-01_dp, 9.5035917926725233e-01_dp, &
       9.5305304494831944e-01_dp, 9.5810340284515405e-01_dp, 9.6315376074198866e-01_dp, &
       9.6584762642305577e-01_dp, 9.6598940881813178e-01_dp, 9.6625521712849149e-01_dp, &
       9.6652102543885121e-01_dp, 9.6921489111991832e-01_dp, 9.7426524901675293e-01_dp, &
       9.7931560691358754e-01_dp, 9.8200947259465465e-01_dp, 9.8215125498973066e-01_dp, &
       9.8241706330009038e-01_dp, 9.8268287161045009e-01_dp, 9.8411978208605433e-01_dp, &
       9.8681364776712144e-01_dp, 9.8950751344818855e-01_dp, 9.9094442392379278e-01_dp, &
       9.9102005078748334e-01_dp, 9.9116183318255935e-01_dp, 9.9130361557763536e-01_dp, &
       9.9274052605323959e-01_dp, 9.9543439173430670e-01_dp, 9.9812825741537381e-01_dp, &
       9.9956516789097805e-01_dp, 9.9964079475466860e-01_dp, 9.9978257714974461e-01_dp, &
       9.9992435954482062e-01_dp, 9.9999998640851118e-01_dp]
    
    old_vals = &
      [1.9533203169308873e-29_dp, 7.9347236391691476e-29_dp, 1.6193312297255380e-28_dp, &
       2.2174715619493640e-28_dp, 9.0332704119052523e-28_dp, 9.6314107441290783e-28_dp, &
       3.5483693504475345e-27_dp, 3.6907692702507794e-27_dp, 4.4321631884687508e-27_dp, &
       1.0321111180327492e-26_dp, 1.0380925213549875e-26_dp, 1.3849947327605718e-26_dp, &
       5.4412717310991178e-26_dp, 5.4472531344213560e-26_dp, 5.7941553458269403e-26_dp, &
       1.1200774782401671e-25_dp, 1.1206756185723909e-25_dp, 1.1553658397129493e-25_dp, &
       2.3876659651949951e-25_dp, 2.3890899643930275e-25_dp, 2.3965039035752072e-25_dp, &
       2.4906817449665769e-25_dp, 2.9315978062732138e-25_dp, 3.5075481114034691e-25_dp, &
       5.1177162738496365e-25_dp, 5.1183144141818603e-25_dp, 5.1530046353224188e-25_dp, &
       7.5051869070129385e-25_dp, 1.1964283773220593e-24_dp, 1.1965707772418626e-24_dp, &
       1.1973121711600805e-24_dp, 1.2067299552992175e-24_dp, 1.2508215614298812e-24_dp, &
       1.3084165919429067e-24_dp, 1.7081804715038537e-24_dp, 2.2177659996660244e-24_dp, &
       2.2179083995858276e-24_dp, 2.2186497935040456e-24_dp, 2.2280675776431826e-24_dp, &
       2.2721591837738463e-24_dp, 2.3297542142868718e-24_dp, 2.7295180938478187e-24_dp, &
       9.2521339464760921e-24_dp, 9.2522763463958953e-24_dp, 9.2530177403141133e-24_dp, &
       9.2624355244532503e-24_dp, 9.3065271305839140e-24_dp, 9.3641221610969395e-24_dp, &
       9.7638860406578864e-24_dp, 1.5311232927385396e-23_dp, 1.5311292741418618e-23_dp, &
       1.5314761763532674e-23_dp, 1.5549979990701726e-23_dp, 1.6507641771504286e-23_dp, &
       1.7528979393848251e-23_dp, 2.4563347340658319e-23_dp, 9.2925534444472239e-22_dp, &
       9.2925548684464219e-22_dp, 9.2925622823856041e-22_dp, 9.2926564602269955e-22_dp, &
       9.2930973762883021e-22_dp, 9.2936733265934324e-22_dp, 9.2976709653890418e-22_dp, &
       9.4456655783890462e-22_dp]
    
    new_bins = &
      [0.0000000000000000e+00_dp, 3.2258064516129031e-02_dp, 6.4516129032258063e-02_dp, &
       9.6774193548387094e-02_dp, 1.2903225806451613e-01_dp, 1.6129032258064516e-01_dp, &
       1.9354838709677419e-01_dp, 2.2580645161290322e-01_dp, 2.5806451612903225e-01_dp, &
       2.9032258064516125e-01_dp, 3.2258064516129031e-01_dp, 3.5483870967741937e-01_dp, &
       3.8709677419354838e-01_dp, 4.1935483870967738e-01_dp, 4.5161290322580644e-01_dp, &
       4.8387096774193550e-01_dp, 5.1612903225806450e-01_dp, 5.4838709677419351e-01_dp, &
       5.8064516129032251e-01_dp, 6.1290322580645162e-01_dp, 6.4516129032258063e-01_dp, &
       6.7741935483870963e-01_dp, 7.0967741935483875e-01_dp, 7.4193548387096775e-01_dp, &
       7.7419354838709675e-01_dp, 8.0645161290322576e-01_dp, 8.3870967741935476e-01_dp, &
       8.7096774193548387e-01_dp, 9.0322580645161288e-01_dp, 9.3548387096774188e-01_dp, &
       9.6774193548387100e-01_dp, 1.0000000000000000e+00_dp]
       
    allocate(new_vals(size(new_bins)-1))
    
    call rebin(old_bins, old_vals, new_bins, new_vals, ierr=ierr)
    if (ierr /= 0) error stop "rebin: ierr not zero."
    
    ! check results
    correct_vals = &
      [2.8724208982166802e-29_dp, 7.9347236391691476e-29_dp, 1.2617118986741737e-28_dp, &
       1.6193312297255380e-28_dp, 2.2056798002206489e-28_dp, 2.2174715619493640e-28_dp, &
       2.2558110286883035e-28_dp, 9.0332704119052523e-28_dp, 9.2838542599547624e-28_dp, &
       9.6314107441290783e-28_dp, 9.6314107441290783e-28_dp, 2.1116994416713538e-27_dp, &
       3.5483693504475345e-27_dp, 3.6704901667874219e-27_dp, 3.6907692702507794e-27_dp, &
       3.6907692702507794e-27_dp, 4.3453679147658879e-27_dp, 4.4321631884687508e-27_dp, &
       4.4321631884687508e-27_dp, 9.7846147414425785e-27_dp, 1.0380925213549875e-26_dp, &
       1.2029548709327527e-26_dp, 1.3849947327605718e-26_dp, 1.2528297955917184e-25_dp, &
       2.3882840071442248e-25_dp, 2.3890899643930275e-25_dp, 2.3953720293626467e-25_dp, &
       2.4210508623901413e-25_dp, 3.3489531942627177e-25_dp, 1.9104215635672634e-24_dp, &
       2.5900260554588492e-22_dp]
    
    do i = 1,size(new_vals)
      if (.not. is_close(new_vals(i), correct_vals(i), tol = 1.0e-10_dp)) then
        error stop "rebin: returned incorrect values"
      endif
    enddo
    
  end subroutine

  subroutine test_conserving_rebin()
      
    real(dp), allocatable :: old_bins(:)
    real(dp), allocatable :: old_vals(:)
    real(dp), allocatable :: new_bins(:)
    real(dp), allocatable :: new_vals(:)
    integer :: ierr
    real(dp), parameter :: tol = 1.0e-10_dp


    ! Case 1: perfect overlap
    old_bins = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    old_vals = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]

    new_bins = [0.0_dp, 1.1_dp, 3.1_dp, 3.2_dp, 4.0_dp]
    allocate(new_vals(size(new_bins)-1))

    call conserving_rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    if (ierr /= 0) error stop "conserving_rebin: ierr not zero."

    deallocate(old_bins, old_vals, new_bins, new_vals)

    ! Case 2: new bins shifted left
    old_bins = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    old_vals = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]

    new_bins = [-1.0_dp, -0.1_dp, 1.1_dp, 3.1_dp, 3.2_dp, 3.4_dp]
    allocate(new_vals(size(new_bins)-1))

    call conserving_rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    if (ierr /= 0) error stop "conserving_rebin: ierr not zero."

    deallocate(old_bins, old_vals, new_bins, new_vals)

    ! Case 3: new bins shifted right
    old_bins = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    old_vals = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]

    new_bins = [2.9_dp, 2.99_dp, 3.01_dp, 3.1_dp, 3.2_dp, 10.0_dp]
    allocate(new_vals(size(new_bins)-1))

    call conserving_rebin(old_bins, old_vals, new_bins, new_vals, ierr)
    if (ierr /= 0) error stop "conserving_rebin: ierr not zero."

    deallocate(old_bins, old_vals, new_bins, new_vals)

  end subroutine

  subroutine test_grid_at_resolution()

    real(dp), allocatable :: wavl(:), Rs(:)
    real(dp), parameter :: wv_min = 1.0_dp
    real(dp), parameter :: wv_max = 5.0_dp
    real(dp), parameter :: R = 10.0_dp
    integer :: ierr

    call grid_at_resolution(wv_min, wv_max, R, wavl, ierr)
    if (ierr /= 0) error stop "grid_at_resolution: ierr not zero."

    Rs = wavl(1:size(wavl)-1)/(wavl(2:) - wavl(1:size(wavl)-1))

    ! Gridding should be equal to R, except the last point.
    if (.not. all(is_close(Rs(1:size(Rs)-1), R))) then
      error stop "grid_at_resolution: returned incorrect values"
    endif

  end subroutine

  subroutine test_make_bins()
    real(dp), allocatable :: wv(:), wavl(:)
    integer :: ierr

    ! Correct values were computed with https://github.com/ACCarnall/spectres
    real(dp), parameter :: wavl_correct(*) = &
    [ &
    7.7777777777777779e-01_dp, 1.2222222222222223e+00_dp, 1.6666666666666665e+00_dp, 2.1111111111111107e+00_dp, &
    2.5555555555555554e+00_dp, 3.0000000000000000e+00_dp, 3.4444444444444446e+00_dp, 3.8888888888888884e+00_dp, &
    4.3333333333333330e+00_dp, 4.7777777777777777e+00_dp, 5.2222222222222223e+00_dp &
    ]

    allocate(wv(10))
    allocate(wavl(size(wv)+1))
    call linspace(1.0_dp, 5.0_dp, wv)
    call make_bins(wv, wavl, ierr)
    if (ierr /= 0) error stop "make_bins: ierr not zero."

    if (.not. all(is_close(wavl, wavl_correct, tol = 1.0e-10_dp))) then
      error stop "make_bins: returned incorrect errors"
    endif

  end subroutine

  subroutine test_brent()

    real(dp) :: r, fzero
    integer :: iflag

    real(dp), parameter :: ax = 0.0_dp
    real(dp), parameter :: pi = 3.141592653589793238462643383279_dp
    real(dp), parameter :: bx = 2.0_dp*pi
    real(dp), parameter :: tol = 1.0e-8_dp

    type(myfunc_type) :: myfunc

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' brent_test'
    write(*,*) '---------------'
    write(*,*) ''

    call myfunc%set_function(sin_func) ! set the function

    ! call fmin:
    ! [the minimum is at 270 deg]
    myfunc%i = 0
    r = myfunc%minimize(ax, bx, tol)
    write(*,*) 'minimum of sin(x) at: ', r*180.0_dp/pi,' deg'
    write(*,*) 'number of function calls: ', myfunc%i

    ! call zeroin:
    ! [the root is at pi]
    myfunc%i = 0
    call myfunc%find_zero(ax+0.0001_dp, bx/2.0_dp+0.0002, tol, r, fzero, iflag)
    if (iflag < 0) then
      error stop "find_zero failed to find root"
    endif
    write(*,*) 'root of sin(x) at: ', r*180.0_dp/pi,' deg'
    write(*,*) 'number of function calls: ', myfunc%i

  end subroutine

  function sin_func(me,x) result(f)
    !! Example function to minimize: sin(x)
    class(brent_class),intent(inout) :: me
    real(dp),intent(in) :: x
    real(dp) :: f

    f = sin(x)

    select type (me)
    class is (myfunc_type)
        me%i = me%i + 1 ! number of function calls
    end select

  end function
  
end program