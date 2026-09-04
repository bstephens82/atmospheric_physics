module clubb_mf

! =============================================================================== !
! Mass-flux module for use with CLUBB                                             !
! Together (CLUBB+MF) they comprise a eddy-diffusivity mass-flux approach (EDMF)  !
! =============================================================================== !
  use ccpp_kinds,    only: kind_phys
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun
  use time_manager,  only: is_first_step, get_nstep
  use spmd_utils,    only: iam
  use physconst,     only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                           cpwv, cpliq, rh2o, zvir, pi

  implicit none
  private
  save

  public :: integrate_mf, &
            clubb_mf_readnl, &
            do_clubb_mf, &
            do_clubb_mf_diag, &
            clubb_mf_nup, &
            do_clubb_mf_rad, &
            clubb_mf_Lopt, &
            clubb_mf_ddalph, &
            clubb_mf_up_ndt, &
            clubb_mf_cp_ndt, &
            do_clubb_mf_cmt, &
            do_clubb_mf_addtke, &
            clubb_mf_cldfrac_fac

  !
  ! Lopt 0 = fixed L0
  !      1 = tke_clubb L0
  !      2 = wpthlp_clubb L0
  !      3 = test plume L0
  !      4 = lel
  !      5 = cape
  !      6 = ztopm1
  !      7 = rel.hum. at 500 hPa
  !      8 = column int. rel.hum.
  integer  :: clubb_mf_Lopt    = 6
  real(kind_phys) :: clubb_mf_a0      = 0.15_kind_phys
  real(kind_phys) :: clubb_mf_b0      = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_L0      = 200._kind_phys
  real(kind_phys) :: clubb_mf_ent0    = 0.2_kind_phys
  real(kind_phys) :: clubb_mf_alphturb= 3.0_kind_phys
  real(kind_phys) :: clubb_mf_max_L0  = 1.e3_kind_phys
  real(kind_phys) :: clubb_mf_fdd     = 0.5_kind_phys
  real(kind_phys) :: clubb_mf_ddalph  = 12.0_kind_phys
  real(kind_phys) :: clubb_mf_ddbeta  = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_pwfac   = 1.0_kind_phys
  real(kind_phys) :: clubb_mf_ddexp   = 3.0_kind_phys
  real(kind_phys) :: clubb_mf_pwmin   = 1.5_kind_phys
  real(kind_phys) :: clubb_mf_pwmax   = 3.0_kind_phys
  real(kind_phys) :: clubb_mf_cldfrac_fac = 10._kind_phys
  integer  :: clubb_mf_up_ndt  = 1
  integer  :: clubb_mf_cp_ndt  = 1
  integer  :: clubb_mf_kseed = 1
  integer, protected :: clubb_mf_nup     = 0
  logical, protected :: do_clubb_mf = .false.
  logical, protected :: do_clubb_mf_diag = .false.
  logical, protected :: do_clubb_mf_rad = .true.
  logical, protected :: do_clubb_mf_addtke = .true.
  logical, protected :: do_clubb_mf_coldpool = .true.
  logical, protected :: do_clubb_mf_ustar = .false.
  logical, protected :: do_clubb_mf_mixd = .false.
  logical, protected :: do_clubb_mf_precip = .true.
  logical, protected :: do_clubb_mf_rhtke = .true.
  logical, protected :: do_clubb_mf_cmt = .false.
  logical, protected :: do_clubb_mf_aloft = .true.
  logical, protected :: do_clubb_mf_coldpool_init = .false.
  logical, protected :: do_clubb_mf_coldpool_perplume = .true.
  logical, protected :: do_clubb_mf_lscale_perplume = .true.
  logical :: tht_tweaks = .true.
  integer :: mf_num_cin = 5

  contains

  subroutine clubb_mf_readnl(nlfile)

  ! =============================================================================== !
  ! MF namelists                                                                    !
  ! =============================================================================== !

    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_integer, mpi_logical

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    character(len=*), parameter :: sub = 'clubb_mf_readnl'

    integer :: iunit, read_status, ierr


    namelist /clubb_mf_nl/ clubb_mf_Lopt, clubb_mf_a0, clubb_mf_b0, clubb_mf_L0, clubb_mf_ent0, clubb_mf_alphturb, &
                           clubb_mf_nup, clubb_mf_max_L0, do_clubb_mf, do_clubb_mf_diag, do_clubb_mf_precip, do_clubb_mf_rad, &
                           clubb_mf_fdd, do_clubb_mf_coldpool, clubb_mf_ddalph, clubb_mf_ddbeta, clubb_mf_pwfac, do_clubb_mf_ustar, &
                           clubb_mf_ddexp, do_clubb_mf_mixd, clubb_mf_up_ndt, clubb_mf_cp_ndt, do_clubb_mf_rhtke, do_clubb_mf_cmt, &
                           do_clubb_mf_coldpool_init, do_clubb_mf_coldpool_perplume, do_clubb_mf_lscale_perplume, clubb_mf_kseed, &
                           do_clubb_mf_addtke, do_clubb_mf_aloft, clubb_mf_pwmin, clubb_mf_pwmax, clubb_mf_cldfrac_fac

    if (masterproc) then
      open( newunit=iunit, file=trim(nlfile), status='old' )
      call find_group_name(iunit, 'clubb_mf_nl', status=read_status)
      if (read_status == 0) then
         read(iunit, clubb_mf_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun('clubb_mf_readnl: ERROR reading namelist')
         end if
      end if
      close(iunit)
    end if

    call mpi_bcast(clubb_mf_Lopt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_Lopt")
    call mpi_bcast(clubb_mf_a0,   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_a0")
    call mpi_bcast(clubb_mf_b0,   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_b0")
    call mpi_bcast(clubb_mf_L0,   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_L0")
    call mpi_bcast(clubb_mf_ent0, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ent0")
    call mpi_bcast(clubb_mf_alphturb,1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_alphturb")
    call mpi_bcast(clubb_mf_nup,  1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_nup")
    call mpi_bcast(clubb_mf_max_L0,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_max_L0")
    call mpi_bcast(do_clubb_mf,      1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf")
    call mpi_bcast(do_clubb_mf_diag, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_diag")
    call mpi_bcast(do_clubb_mf_precip, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_precip")
    call mpi_bcast(do_clubb_mf_rad, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_rad")
    call mpi_bcast(clubb_mf_fdd,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_fdd")
    call mpi_bcast(do_clubb_mf_coldpool, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_coldpool")
    call mpi_bcast(clubb_mf_ddalph,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddalph")
    call mpi_bcast(clubb_mf_ddbeta,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddbeta")
    call mpi_bcast(clubb_mf_pwfac,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_pwfac")
    call mpi_bcast(clubb_mf_up_ndt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_up_ndt")
    call mpi_bcast(clubb_mf_cp_ndt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_cp_ndt")
    call mpi_bcast(do_clubb_mf_ustar, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_ustar")
    call mpi_bcast(clubb_mf_ddexp,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddexp")
    call mpi_bcast(do_clubb_mf_mixd, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_mixd")
    call mpi_bcast(do_clubb_mf_rhtke, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_rhtke")
    call mpi_bcast(do_clubb_mf_cmt, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_cmt")
    call mpi_bcast(clubb_mf_kseed, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_kseed")
    call mpi_bcast(do_clubb_mf_coldpool_init, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_coldpool_init")
    call mpi_bcast(do_clubb_mf_coldpool_perplume, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_coldpool_perplume")
    call mpi_bcast(do_clubb_mf_lscale_perplume, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_lscale_perplume")
    call mpi_bcast(do_clubb_mf_addtke, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_addtke")
    call mpi_bcast(do_clubb_mf_aloft, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_aloft")
    call mpi_bcast(clubb_mf_pwmin,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_pwmin")
    call mpi_bcast(clubb_mf_pwmax,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_pwmax")
    call mpi_bcast(clubb_mf_cldfrac_fac,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_cldfrac_fac")

    if ((.not. do_clubb_mf) .and. do_clubb_mf_diag ) then
       call endrun('clubb_mf_readnl: Error - cannot turn on do_clubb_mf_diag without also turning on do_clubb_mf')
    end if

  end subroutine clubb_mf_readnl


  subroutine integrate_mf( nzm,     nzt,     dtime,                                 & ! input
                           rho_zm,           zm,      p_zm,      iexner_zm,         & ! input
                           rho_zt,  dzt,     zt,      p_zt,      iexner_zt,         & ! input
                           u,       v,       thl,     qt,        thv,               & ! input
                                    w,       th,      qv,        qc,                & ! input
                                             thl_zm,  qt_zm,     thv_zm,            & ! input
                                             th_zm,   qv_zm,     qc_zm,             & ! input
                           ustar,   ths,     wthl_sfc,wqt_sfc,   pblh,              & ! input
                           tke,     tpert,   rhinv,                                 & ! input
                           wpthlp_env,       wpthvp_env,         wpqtp_env,         & ! input
                           ztopm1,           ddcp,               cbm1,              & ! inout
                           mcape,                                                   & ! output
                           upa,     dna,                                            & ! output
                           upw,     dnw,                                            & ! output
                           upmf,                                                    & ! output
                           upqt,    dnqt,                                           & ! output
                           upthl,   dnthl,                                          & ! output
                           upthv,   dnthv,                                          & ! output
                           upth,    dnth,                                           & ! output
                           upqc,    dnqc,                                           & ! output
                           upbuoy,                                                  & ! output
                           upent,                                                   & ! output
                           updet,                                                   & ! output
                           dry_a,   moist_a,                                        & ! output
                           dry_w,   moist_w,                                        & ! output
                           dry_qt,  moist_qt,                                       & ! output
                           dry_thl, moist_thl,                                      & ! output
                           dry_u,   moist_u,                                        & ! output
                           dry_v,   moist_v,                                        & ! output
                                    moist_qc,                                       & ! output
                           ae,                                                      &
                           ac,      aup,     adn,                                   &
                           aw,      awup,    awdn,                                  &
                           aww,     awwup,   awwdn,                                 &
                           awthlup, awqtup,  awuup, awvup,                          & ! output
                           awthldn, awqtdn,  awudn, awvdn,                          & ! output
                           awthl,   awqt,                                           & ! output
                           awu,     awv,                                            & ! output
                           thlflxup,qtflxup, uflxup, vflxup,                        & ! output
                           thlflxdn,qtflxdn, uflxdn, vflxdn,                        & ! output
                           thlflx,  qtflx,   uflx,   vflx,                          & ! output - variables needed for solver
                           thvflx,                                                  & ! output
                           sqtup,   sthlup,                                         & ! output
                           sqtdn,   sthldn,                                         & ! output
                           sqt,     sthl,                                           & ! output - variables needed for solver
                           precc,                                                   & ! output
                           ztop,    dynamic_L0 )

  ! ================================================================================= !
  ! Mass-flux algorithm                                                               !
  !                                                                                   !
  ! Provides rtm and thl fluxes due to mass flux ensemble,                            !
  ! which are fed into the mixed explicit/implicit clubb solver as explicit terms     !
  !                                                                                   !
  ! Mass flux variables are computed on edges (i.e. momentum grid):                   !
  ! upa,upw,upqt,...                                                                  !
  ! dry_a,moist_a,dry_w,moist_w, ...                                                  !
  !                                                                                   !
  ! In CLUBB (unlike CAM) nlevs of momentum grid = nlevs of thermodynamic grid,       !
  ! due to a subsurface thermodynamic layer. To avoid confusion, below the variables  !
  ! are grouped by the grid they are on.                                              !
  !                                                                                   !
  ! *note that state on the lowest thermo level is equal to state on the lowest       !
  ! momentum level due to state_zt(1) = state_zt(2), and lowest momentum level        !
  ! is a weighted combination of the lowest two thermodynamic levels.                 !
  !                                                                                   !
  ! ---------------------------------Authors----------------------------------------  !
  ! Marcin Kurowski, JPL                                                              !
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM       !
  ! Additional modifications by Adam Herrington, NCAR                                 !
  ! ================================================================================= !

     use wv_saturation,      only : qsat

     integer,  intent(in)                 :: nzm, nzt
     real(kind_phys), intent(in)                 :: dtime                   ! host-model physics timestep [s]
     real(kind_phys), dimension(nzt), intent(in) :: u,      v,            & ! thermodynamic grid
                                             thl,    thv,          & ! thermodynamic grid
                                             th,     qv,           & ! thermodynamic grid
                                             qt,     qc,           & ! thermodynamic grid
                                             p_zt,   iexner_zt,    & ! thermodynamic grid
                                             dzt,    rho_zt,       & ! thermodynamic grid
                                             zt

     real(kind_phys), dimension(nzm), intent(in) :: thl_zm, thv_zm,       & ! momentum grid
                                             w,                    &
                                             th_zm,  qv_zm,        &
                                             qt_zm,  qc_zm,        & ! momentum grid
                                             p_zm,   iexner_zm,    & ! momentum grid
                                                     rho_zm,       & ! momentum grid
                                             zm,                   & ! momentum grid
                                             tke,    wpthlp_env,   & ! momentum grid
                                             wpthvp_env, wpqtp_env

     real(kind_phys), intent(in)                :: wthl_sfc,wqt_sfc
     real(kind_phys), intent(in)                :: pblh,tpert
     real(kind_phys), intent(in)                :: rhinv
     real(kind_phys), intent(in)                :: ths,ustar
     real(kind_phys), intent(inout)             :: cbm1

     real(kind_phys),dimension(clubb_mf_nup), intent(inout)  :: ztopm1,ddcp

     real(kind_phys),dimension(nzm,clubb_mf_nup), intent(out) :: upa,     & ! momentum grid
                                                          upw,     & ! momentum grid
                                                          upmf,    & ! momentum grid
                                                          upqt,    & ! momentum grid
                                                          upthl,   & ! momentum grid
                                                          upthv,   & ! momentum grid
                                                          upth,    & ! momentum grid
                                                          upqc,    & ! momentum grid
                                                          upbuoy,  & ! momentum grid
                                                          upent,   & ! momentum grid
                                                          updet

     real(kind_phys),dimension(nzm,clubb_mf_nup), intent(out) :: dna,     & ! momentum grid
                                                          dnw,     & ! momentum grid
                                                          dnqt,    & ! momentum grid
                                                          dnthl,   & ! momentum grid
                                                          dnthv,   & ! momentum grid
                                                          dnth,    & ! momentum grid
                                                          dnqc

     real(kind_phys),dimension(nzm), intent(out) :: dry_a,   moist_a,     & ! momentum grid
                                             dry_w,   moist_w,     & ! momentum grid
                                             dry_qt,  moist_qt,    & ! momentum grid
                                             dry_thl, moist_thl,   & ! momentum grid
                                             dry_u,   moist_u,     & ! momentum grid
                                             dry_v,   moist_v,     & ! momentum grid
                                                      moist_qc       ! momentum grid

     real(kind_phys),dimension(nzm), intent(out) :: ae,                                &
                                             ac,      aup,     adn,              &
                                             aw,      awup,    awdn,             &
                                             aww,     awwup,  awwdn,             &
                                             awthlup, awqtup, awuup, awvup,      & ! momentum grid
                                             awthldn, awqtdn, awudn, awvdn,      & ! momentum grid
                                             awthl,   awqt,                      & ! momentum grid
                                             awu,     awv,                       & ! momentum grid
                                             thlflxup,qtflxup, uflxup, vflxup,   & ! momentum grid
                                             thlflxdn,qtflxdn, uflxdn, vflxdn,   & ! momentum grid
                                             thlflx,  qtflx,   uflx,   vflx,     & ! momentum grid
                                             thvflx,  precc

     real(kind_phys),dimension(nzt), intent(out) :: sqtup,   sthlup,                   & ! thermodynamic grid
                                             sqtdn,   sthldn,                    & ! thermodynamic grid
                                             sqt,     sthl                         ! thermodynamic grid

     real(kind_phys),dimension(clubb_mf_nup), intent(out) :: ztop, dynamic_L0, mcape

     ! =============================================================================== !
     ! INTERNAL VARIABLES
     !

     ! =============================================================================== !
     ! GRID ORIENTATION GENERALIZATION VARIABLES
     ! ------------------------------------------------------------------------------- !
     ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
     ! duplicating code, these variables abstract the vertical loop bounds and slices.
     !
     ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
     ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
     ! kdir          : Directional step (+1 for moving up, -1 for moving down).
     !
     ! STAGGERED GRID INDEXING
     ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
     ! index of the cell center (zt) to the interface (zm) flips depending on whether
     ! memory is loaded top-down or bottom-up. These variables dynamically map them:
     !
     ! kt    : The active thermodynamic cell center associated with the current step.
     !         Upward Sweep:   kt = k - (1-kdir)/2
     !         Downward Sweep: kt = k - (1+kdir)/2
     !
     ! kn    : The NEXT momentum interface in the direction of the current sweep.
     !         Upward Sweep:   kn = k + kdir
     !         Downward Sweep: kn = k - kdir
     !
     ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
     !         kt_up = k - (1-kdir)/2
     !
     ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
     !         kt_dn = k - (1+kdir)/2
     ! =============================================================================== !
     integer :: ksfcm, ktopm, ksfct, ktopt, kdir, kt, kn, kt_up, kt_dn

     ! sums over all plumes
     real(kind_phys), dimension(nzm)              :: moist_th,   dry_th,      &
                                              thl_env,    qt_env,      &
                                              thv_env,                 &
                                              thvflxup,   thvflxdn,    &
                                              awthvup,    awthvdn
     ! updraft properties
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: upqv,     upqs,          & ! momentum grid
                                              upql,     upqi,          & ! momentum grid
                                              upu,      upv,           & ! momentum grid
                                              uplmix                     ! momentum grid
     ! downdraft properties
     real(kind_phys), dimension(nzm,clubb_mf_nup) ::           dnqs,          & ! momentum grid
                                              dnql,     dnqi,          & ! momentum grid
                                              dnu,      dnv,           & ! momentum grid
                                              dnlmix                     ! momentum grid
     ! microphyiscs terms
     real(kind_phys), dimension(nzt,clubb_mf_nup) :: supqt,    supthl,        & ! thermodynamic grid
                                              sdnqt,    sdnthl,        & ! thermodynamic grid
                                              upauto                     ! thermodynamic grid
     ! precipitation rates
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: uprr,     dnrr             ! momentum grid
     !
     ! grid-mean (area-weighted) rain reservoirs. uprr/dnrr above are
     ! per-unit-plume-area and drive the local evaporation physics; uprg/dnrg carry
     ! the same budgets in grid-mean units (weighted by the plume areas used in the
     ! sqtup/sqtdn accumulations) so that the area-weighted evaporation can never
     ! exceed the area-weighted rain supply. This guarantees the column integral of
     ! sqt (and hence surface precc / PRECC) is non-negative.
     real(kind_phys), dimension(nzm,clubb_mf_nup) :: uprg,     dnrg             ! momentum grid
     real(kind_phys)                              :: drytot,   lamk,   sdnmax
     real(kind_phys), parameter                   :: mf_dry_frac_max = 0.9_kind_phys   ! max fractional total water removable per timestep
     !        
     ! entrainment profiles
     real(kind_phys), dimension(nzt,clubb_mf_nup) :: entf,     mix              ! thermodynamic grid
     integer,  dimension(nzt,clubb_mf_nup) :: enti                       ! thermodynamic grid
     !
     ! other variables
     integer                              :: k,i,kstart,ddtopm,kcb,kpbl,kmid,nbot
     integer,  dimension(clubb_mf_nup)    :: ddbotm,kcbarr
     real(kind_phys), dimension(clubb_mf_nup)    :: zcb,cpfac
     real(kind_phys)                             :: zcb_unset,                &
                                             wthv_sfc, wthv,   wqt,    &
                                             ddint,   iddcp,   &
                                             wstar,  qstar,   thvstar, &
                                             sigmaw, sigmaqt, sigmathv,&
                                             convh,  wmin,    wmax,    &
                                             wlv,    wtv,     wp,      &
                                             B,                        & ! thermodynamic grid
                                             entexp, entexpu, entw,    & ! thermodynamic grid
                                             Mn,                       & ! momentum grid
                                             eturb,  det,     lmixt,   & ! thermodynamic grid
                                             qtovqs, sevap,   taum1,   & ! thermodynamic grid
                                             sqtint, sthlint, alphint, &
                                             qtmp,   betathl, betaqt,  & ! thermodynamic grid
                                             thln,   thvn,    thn,     & ! momentum grid
                                             qtn,    qsn,              & ! momentum grid
                                             qcn,    qln,     qin,     & ! momentum grid
                                             un,     vn,      wn2,     & ! momentum grid
                                             wn,                       & ! momentum grid
                                             lmixn,   srfarea,         & ! momentum grid
                                             srfwqtu, srfwthvu,        &
                                             facqtu,  facthvu,         &
                                             zsub,    wcb,    rh_L0,   &
                                             dzext

     ! limit convective area
     logical                                :: limarea = .false.
     real(kind_phys),parameter                     :: amax = 0.6_kind_phys

     ! buoyancy sorting variables
     logical                                :: bsort = .false.
     real(kind_phys),parameter                     :: rle = 0.1_kind_phys
     integer                                :: niter_xc = 1
     integer                                :: kk,      status,  iter_xc
     real(kind_phys)                               :: tlm,     excessm, qsm,     &
                                               tln,     excessn, es,      &
                                               xc,      xsat,    x_en,    &
                                               x_cu,    xs1,     xs2,     &
                                               aquad,   bquad,   cquad,   &
                                               thlxsat, thvxsat, qtxsat,  &
                                               thv_x0,  thv_x1,  cridis,  &
                                               thln0,   qtn0,    wn0,     &
                                               entn,    detn,    mfn,     &
                                               ee2,     ud2
     ! aloft trigger flag
     logical                                :: aloft = .false.
     ! rng seed
     real(kind_phys), dimension(4)                 :: u_seed
     
     ! alpha relates star qunataties to stddev after Suselj etal 2019
     real(kind_phys),parameter                     :: alphw   = 0.572_kind_phys,        &
                                               alphqt  = 2.890_kind_phys,        &
                                               alphthv = 2.890_kind_phys
     ! w' covariance after Suselj etal 2019
     real(kind_phys),parameter                     :: cwqt  = 0.32_kind_phys,           &
                                               cwthv = 0.58_kind_phys
     ! virtual mass coefficients for w-eqn after Suselj etal 2019
     real(kind_phys),parameter                     :: wa = 1.0_kind_phys,               &
                                               wb = 1.5_kind_phys
     ! min values to avoid singularities
     real(kind_phys),parameter                     :: wstarmin = 1.e-3_kind_phys,       &
                                               pblhmin  = 100._kind_phys
     ! evaporation efficiency after Suselj etal 2019
     real(kind_phys),parameter                     :: ke = 2.5e-4_kind_phys
     !
     ! height here downdrafts feel the surface
     real(kind_phys),parameter                     :: z00dn = 1.e3_kind_phys, &
                                               tinynum = 1.e-7_kind_phys
     ! to fix entrainmnet rate
     logical                                :: fixent = .false.
     !
     ! fixed entrainment rate
     real(kind_phys),parameter                     :: fixent_ent = 2.e-4_kind_phys
     !
     ! Arakawa and Schubert detrainment limiter
     logical                                :: do_aspd = .false.
     !
     ! Lower limit on entrainment length scale
     real(kind_phys),parameter                     :: min_L0 = 0.5_kind_phys
     !
     ! limiter for tke enahnced fractional entrainment
     ! (only used when do_aspd = .true.)
     real(kind_phys),parameter                     :: max_eturb = 10._kind_phys
     !
     ! to condensate or not to condensate
     logical                                :: do_condensation = .true.
     !
     ! use implicit method for plume updraft velocity
     logical                                :: do_implicit = .false.
     !
     ! to scale surface fluxes
     logical                                :: scalesrf = .false.
     !
     ! minimum downdraft speed
     real(kind_phys),parameter                     :: mindnw = 1.E-2_kind_phys
     !
     ! limiter on cold pool effects
     real(kind_phys),parameter                     :: max_cpfac = 5._kind_phys
     !
     ! max limiter on cold pool init effects
     real(kind_phys),parameter                     :: max_cpinit = 0.5_kind_phys

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! DETERMINE GRID ORIENTATION
     if (zm(1) < zm(nzm)) then
        ! Bottom-up
        ksfcm = 1
        ktopm = nzm
        kdir  = 1
     else
        ! Top-down
        ksfcm = nzm
        ktopm = 1
        kdir  = -1
     end if

     if (zt(1) < zt(nzt)) then
        ksfct = 1
        ktopt = nzt
     else
        ksfct = nzt
        ktopt = 1
     end if

     ! INITIALIZE OUTPUT VARIABLES
     dry_a     = 0._kind_phys
     moist_a   = 0._kind_phys
     dry_w     = 0._kind_phys
     moist_w   = 0._kind_phys
     dry_qt    = 0._kind_phys
     moist_qt  = 0._kind_phys
     dry_thl   = 0._kind_phys
     moist_thl = 0._kind_phys
     dry_u     = 0._kind_phys
     moist_u   = 0._kind_phys
     dry_v     = 0._kind_phys
     moist_v   = 0._kind_phys
     moist_qc  = 0._kind_phys

     ! this is the environmental area - by default 1.
     ae        = 1._kind_phys
     ac        = 0._kind_phys
     aup       = 0._kind_phys
     adn       = 0._kind_phys
     aw        = 0._kind_phys
     awup      = 0._kind_phys
     awdn      = 0._kind_phys
     aww       = 0._kind_phys
     awwup     = 0._kind_phys
     awwdn     = 0._kind_phys
     awuup     = 0._kind_phys
     awvup     = 0._kind_phys
     awudn     = 0._kind_phys
     awvdn     = 0._kind_phys
     awthvup   = 0._kind_phys
     awthvdn   = 0._kind_phys
     awthlup   = 0._kind_phys
     awqtup    = 0._kind_phys
     awthldn   = 0._kind_phys
     awqtdn    = 0._kind_phys
     awthl     = 0._kind_phys
     awqt      = 0._kind_phys
     awu       = 0._kind_phys
     awv       = 0._kind_phys
     thvflxup  = 0._kind_phys
     thlflxup  = 0._kind_phys
     qtflxup   = 0._kind_phys
     thvflxdn  = 0._kind_phys
     thlflxdn  = 0._kind_phys
     qtflxdn   = 0._kind_phys
     thvflx    = 0._kind_phys
     thlflx    = 0._kind_phys
     uflxup    = 0._kind_phys
     vflxup    = 0._kind_phys
     uflxdn    = 0._kind_phys
     vflxdn    = 0._kind_phys
     uflx      = 0._kind_phys
     vflx      = 0._kind_phys
     qtflx     = 0._kind_phys
     sqtup     = 0._kind_phys
     sthlup    = 0._kind_phys
     sqtdn     = 0._kind_phys
     sthldn    = 0._kind_phys
     sqt       = 0._kind_phys
     sthl      = 0._kind_phys
     precc     = 0._kind_phys

     mix       = 0._kind_phys
     entf      = 0._kind_phys
     enti      = 0
     det       = 0._kind_phys

     ! START MAIN COMPUTATION
     upw   = 0._kind_phys
     upth  = 0._kind_phys
     upthl = 0._kind_phys
     upthv = 0._kind_phys
     upqt  = 0._kind_phys
     upa   = 0._kind_phys
     upmf  = 0._kind_phys
     upu   = 0._kind_phys
     upv   = 0._kind_phys
     upqc  = 0._kind_phys
     upth  = 0._kind_phys
     upql  = 0._kind_phys
     upqi  = 0._kind_phys
     upqv  = 0._kind_phys
     upqs  = 0._kind_phys
     upbuoy= 0._kind_phys
     uplmix= 0._kind_phys
     uprr  = 0._kind_phys
     uprg  = 0._kind_phys
     supqt = 0._kind_phys
     supthl= 0._kind_phys
     upent = 0._kind_phys
     updet = 0._kind_phys
     upauto= 0._kind_phys

     dnw   = 0._kind_phys
     dna   = 0._kind_phys
     dnu   = 0._kind_phys
     dnv   = 0._kind_phys
     dnqt  = 0._kind_phys
     dnthl = 0._kind_phys
     dnthv = 0._kind_phys
     dnrr  = 0._kind_phys
     dnrg  = 0._kind_phys
     dnth  = 0._kind_phys
     dnqc  = 0._kind_phys
     dnql  = 0._kind_phys
     dnqi  = 0._kind_phys
     dnqs  = 0._kind_phys
     dnlmix= 0._kind_phys
     sdnqt = 0._kind_phys
     sdnthl= 0._kind_phys

     dynamic_L0 = 0._kind_phys
     ztop = 0._kind_phys
     mcape = 0._kind_phys
     ddbotm = 0

     if (bsort) then
       niter_xc = 3
       limarea = .true.
     end if

     ! unique identifier
     zcb_unset = 9999999._kind_phys
     zcb       = zcb_unset

     ! surface buoyancy flux
     !wthv = wthl+zvir*ths*wqt
     wthv_sfc = wthl_sfc+zvir*ths*wqt_sfc

     if (do_clubb_mf_aloft .and. wthv_sfc < 0.01_kind_phys) then
       aloft = .true.
       kpbl = ksfcm
       do while (zm(kpbl) < pblh)
         kpbl = kpbl + kdir
       end do

       kmid = ksfcm
       ! Use a pressure-based criterion to locate the mid-level within the
       ! troposphere. A threshold of ~500 hPa is preferred over the previous
       ! fixed height (9 km) because it better represents the tropopause
       ! location across different atmospheric conditions.
       do while (p_zm(kmid) > 500.E2_kind_phys)
         kmid = kmid + kdir
       end do

       ! Search absolute bounds by converting relative slice index
       kstart = maxloc(wpthvp_env(kpbl : kmid : kdir), DIM=1)
       kstart = kpbl + (kstart - 1) * kdir

       wthv = wpthvp_env(kstart)
       wqt = wpqtp_env(kstart)

       if (kstart == ktopm) then
         wthv = 0._kind_phys
         wqt  = 0._kind_phys
       end if

     else
       aloft = .false.
       kstart = ksfcm
       wthv = wthv_sfc
       wqt  = wqt_sfc
     end if

     ! if surface buoyancy is positive then do mass-flux
     !if ( wthv > 0._kind_phys ) then
     if ( wthv > 0._kind_phys .and. wqt > 0._kind_phys) then

       if (do_clubb_mf_mixd) then
         convh = max(cbm1,pblhmin)
       else
         convh = max(pblh,pblhmin)
       end if

       ! --------------------------------------------------------- !
       ! Initialize using Deardorff convective velocity scale      !
       ! --------------------------------------------------------- !
       wstar = max( wstarmin, (gravit/thv(kstart - (1-kdir)/2)*wthv*convh)**(1._kind_phys/3._kind_phys) )

       ! --------------------------------------------------------- !
       ! Compute cold pool feedback parameter                      !
       ! --------------------------------------------------------- !

       cpfac(:) = 1._kind_phys
       if (do_clubb_mf_coldpool) then
         do i=1,clubb_mf_nup
           cpfac(i) = min( (max(ddcp(i)/wstar,1._kind_phys))**clubb_mf_ddbeta, max_cpfac )
         end do
       end if

       ! --------------------------------------------------------- !
       ! Construct tri-variate PDF at the surface from wstar       !
       ! and initialize plume thv, qt, w                           !
       ! --------------------------------------------------------- !

       if (do_clubb_mf_ustar) then
         qstar   = wqt / max(wstarmin,ustar)
         thvstar = wthv / max(wstarmin,ustar)
       else
         qstar   = wqt / wstar
         thvstar = wthv / wstar
       end if

       do i=1,clubb_mf_nup

         if (do_clubb_mf_coldpool_init) then
           sigmaw   = alphw * wstar * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
           sigmaqt  = alphqt * abs(qstar) * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
           sigmathv = alphthv * abs(thvstar) * (1._kind_phys + max_cpinit*cpfac(i)/max_cpfac)
         else
           sigmaw   = alphw * wstar
           sigmaqt  = alphqt * abs(qstar)
           sigmathv = alphthv * abs(thvstar)
         end if

         wmin = sigmaw * clubb_mf_pwmin
         wmax = sigmaw * clubb_mf_pwmax

         wlv = wmin + (wmax-wmin) / (real(clubb_mf_nup,kind_phys)) * (real(i-1, kind_phys))
         wtv = wmin + (wmax-wmin) / (real(clubb_mf_nup,kind_phys)) * real(i,kind_phys)

         upw(kstart,i) = 0.5_kind_phys * (wlv+wtv)
         upa(kstart,i) = 0.5_kind_phys * erf( wtv/(sqrt(2._kind_phys)*sigmaw) ) &
                    - 0.5_kind_phys * erf( wlv/(sqrt(2._kind_phys)*sigmaw) )

         upu(kstart,i) = u(kstart - (1-kdir)/2)
         upv(kstart,i) = v(kstart - (1-kdir)/2)

         upqt(kstart,i)  = cwqt * upw(kstart,i) * sigmaqt/sigmaw
         upthv(kstart,i) = cwthv * upw(kstart,i) * sigmathv/sigmaw
       enddo

       facqtu=1._kind_phys
       facthvu=1._kind_phys
       if (scalesrf) then
         ! scale surface fluxes
         ! (req'd for conservation if not running with zero-flux b.c.'s)
         srfwqtu = 0._kind_phys
         srfwthvu = 0._kind_phys
         srfarea = 0._kind_phys
         do i=1,clubb_mf_nup
             srfwqtu=srfwqtu+upqt(kstart,i)*upw(kstart,i)*upa(kstart,i)
             srfwthvu=srfwthvu+upthv(kstart,i)*upw(kstart,i)*upa(kstart,i)
             srfarea = srfarea+upa(kstart,i)
         end do
         facqtu=srfarea*wqt/srfwqtu
         facthvu=srfarea*wthv/srfwthvu
       end if

       do i=1,clubb_mf_nup

         kt = kstart - (1-kdir)/2
         kn = kstart + kdir

         betaqt = (qt(kt+2*kdir)-qt(kt))/(0.5_kind_phys*(dzt(kt+2*kdir)+2._kind_phys*dzt(kt+kdir)+dzt(kt)))
         betathl = (thv(kt+2*kdir)-thv(kt))/(0.5_kind_phys*(dzt(kt+2*kdir)+2._kind_phys*dzt(kt+kdir)+dzt(kt)))

         if (.not.aloft) then
!jt The following 2 lines I think are correct since the surface is .5 grid down from the surface dzt
!jt           upqt(kstart,i)= qt(kt)-betaqt*0.5_kind_phys*dzt(kt)+facqtu*upqt(kstart,i)
!jt           upthv(kstart,i)= thv(kt)-betathl*0.5_kind_phys*dzt(kt)+facthvu*upthv(kstart,i)
!jt The following 2 lines are to restore bfb accuracy
           upqt(kstart,i)= qt(kt)-betaqt*dzt(kt)+facqtu*upqt(kstart,i)
           upthv(kstart,i)= thv(kt)-betathl*dzt(kt)+facthvu*upthv(kstart,i)
         else
           upqt(kstart,i)= qt(kt)+upqt(kstart,i)
           upthv(kstart,i)= thv(kt)+upthv(kstart,i)
           if (w(kstart) > 0._kind_phys) upw(kstart,i)= w(kstart)+upw(kstart,i)
         end if

         upthl(kstart,i) = upthv(kstart,i) / (1._kind_phys+zvir*upqt(kstart,i))
         upth(kstart,i)  = upthl(kstart,i)
         upmf(kstart,i) = rho_zm(kstart)*upa(kstart,i)*upw(kstart,i)

         ! get cloud, lowest momentum level
         if (do_condensation) then
           call condensation_mf(upqt(kstart,i), upthl(kstart,i), p_zm(kstart), iexner_zm(kstart), &
                                thvn, qcn, thn, qln, qin, qsn, lmixn)
           upthv(kstart,i) = thvn
           upqc(kstart,i)  = qcn
           upql(kstart,i)  = qln
           upqi(kstart,i)  = qin
           upqs(kstart,i)  = qsn
           upth(kstart,i)  = thn
           if (qcn > 0._kind_phys) zcb(i) = zm(kstart)
         else
           ! assume no cldliq
           upqc(kstart,i)  = 0._kind_phys
         end if
       end do

       ! if aloft extend the mass flux plume below kstart nbot levels
       if (aloft) then
         zsub = zm(kstart)
         dzext = 1000._kind_phys
         !if greater than dzext above the surface
         if (zsub > dzext) then
             ! find nbot levs below kstart
             nbot = 0
             do k = kstart-kdir, ksfcm, -kdir
               if ((zm(kstart)-zm(k)) < dzext) then
                 nbot = nbot + 1
               end if
             end do
         else
           !else set dxext to height above the suface
           dzext = zm(kstart) - zm(ksfcm)
           nbot = abs(kstart-ksfcm)
         end if

         zsub = zm(kstart)
         do i=1,clubb_mf_nup
           wcb  = upw(kstart,i)
           do k = kstart-kdir, kstart-nbot*kdir, -kdir
             upw(k,i) = wcb - (wcb/(dzext**clubb_mf_ddexp))*(zsub - zm(k))**clubb_mf_ddexp
             upa(k,i) = upa(kstart,i)
             upmf(k,i) = rho_zm(k)*upa(k,i)*upw(k,i)

             upu(k,i) = upu(kstart,i)
             upv(k,i) = upv(kstart,i)
             upqt(k,i) = upqt(kstart,i)
             upthv(k,i) = upthv(kstart,i)
             upthl(k,i) = upthl(kstart,i)
             upth(k,i)  = upth(kstart,i)
             upqc(k,i)  = upqc(kstart,i)
             upql(k,i)  = upql(kstart,i)
             upqi(k,i)  = upqi(kstart,i)
             upqs(k,i)  = upqs(kstart,i)
           end do
         end do
       end if

       do i=1,clubb_mf_nup
         ! --------------------------------------------------------- !
         ! Calculate ztop and dynamic_L based on value of namelist   !
         ! --------------------------------------------------------- !
         call get_Lscale (nzt, nzm, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                          wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                          do_condensation, qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1(i), dynamic_L0(i), ztop(i), mcape(i))

         ! cold pool feedback on the entrainmnet length scale
         dynamic_L0(i) = dynamic_L0(i) * cpfac(i)

         ! limit max/min
         dynamic_L0(i) = max(min_L0,dynamic_L0(i))
         dynamic_L0(i) = min(clubb_mf_max_L0,dynamic_L0(i))

         ! --------------------------------------------------------- !
         ! Stochastic entrainmnet calculation                        !
         ! From Suselj et al 2019, after Romps and Kuang 2010        !
         ! (ideally we wouldn't fill the entire arrray w/ the RNG,   !
         ! but the RNG doesn't work properly when it operates on     !
         ! the entire array. I'm not sure why this is happening.)    !
         ! --------------------------------------------------------- !
         do k = ksfct, ktopt-kdir, kdir
           ! get entrainment coefficient, dz/L0
           entf(k,i) = dzt(k) / dynamic_L0(i)
         end do
       end do

       ! get poisson, P(dz/L0)
       ! Grab the u-wind from the exact physical near-surface layers to preserve the PRNG seed
       u_seed(1) = u(ksfct)
       u_seed(2) = u(ksfct + (clubb_mf_kseed + 0)* kdir)
       u_seed(3) = u(ksfct + (clubb_mf_kseed + 1) * kdir)
       u_seed(4) = u(ksfct + (clubb_mf_kseed + 2) * kdir)

       ! get poisson, P(dz/L0)
       call poisson( nzt, clubb_mf_nup, ksfct, ktopt, kdir, entf, enti, u_seed )
       ! call poisson( nzt, clubb_mf_nup, entf, enti, u(clubb_mf_kseed+1:clubb_mf_kseed+4))

       ! --------------------------------------------------------- !
       ! Main upward sweep to compute updraft properties           !
       !                                                           !
       ! --------------------------------------------------------- !
       do i=1,clubb_mf_nup
         do k = kstart, ktopm-kdir, kdir

           kt = k - (1-kdir)/2
           kn = k + kdir

           ! get microphysics, autoconversion
           if (do_clubb_mf_precip .and. upqc(k,i) > 0._kind_phys) then
             call precip_mf(upqs(k,i),upqt(k,i),upw(k,i),dzt(kt),abs(zm(kn)-zcb(i)),supqt(kt,i))
             supthl(kt,i) = -1._kind_phys*lmixn*supqt(kt,i)*iexner_zt(kt)/cpair
           else
             supqt(kt,i)  = 0._kind_phys
             supthl(kt,i) = 0._kind_phys
           end if

           ! compute mixing rate
           if (fixent) then
             mix(kt,i) = fixent_ent
           else
             ! get entrainment, ent=ent0/dz*P(dz/L0)
             mix(kt,i) = real( enti(kt,i))*clubb_mf_ent0/dzt(kt)
           end if

           do iter_xc = 1, niter_xc

             if (bsort) then
               if (iter_xc==1) then
                 qtn  = upqt(k,i)
                 thln = upthl(k,i)
                 wn   = upw(k,i)
               else
                 qtn  = 0.5_kind_phys*(qtn + qtn0)
                 thln = 0.5_kind_phys*(thln + thln0)
                 wn = 0.5_kind_phys*(wn + wn0)
               end if

               ! save this iteration
               qtn0  = qtn
               thln0 = thln
               wn0 = wn

               ! --------------------------------------------------------- !
               ! Compute excess water to derive neutral mixing fraction    !
               ! after Bretherton et al 2014                               !
               ! --------------------------------------------------------- !

               ! qexcess of the envrionment
               tlm = thl_zm(kn)/iexner_zm(kn)
               call qsat(tlm,p_zm(kn),es,qsm)
               excessm = qt_zm(kn) - qsm

               ! qexcess in plume
               tln  = thln/iexner_zm(kn)
               call qsat(tln,p_zm(kn),es,qsn)
               excessn = qtn - qsn

               call condensation_mf(qtn, thln, p_zm(kn), iexner_zm(kn), &
                                    thvn, qcn, thn, qln, qin, qsn, lmixn)

               ! critical stopping distance
               cridis = rle*ztopm1(i)

               ! ----------------------------------------------------------------- !
               ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
               ! ----------------------------------------------------------------- !
               if (excessm*excessn > 0._kind_phys) then
                 xc = min(1._kind_phys,max(0._kind_phys,1._kind_phys-2._kind_phys*wa*gravit*cridis/wn**2._kind_phys*(1._kind_phys-thvn/thv_zm(kn))))
                 aquad = 0._kind_phys
                 bquad = 0._kind_phys
                 cquad = 0._kind_phys
               else
                 ! -------------------------------------------------- !
                 ! Case 2 : When either cumulus or env. is saturated. !
                 ! -------------------------------------------------- !
                 xsat    = excessn / ( excessn - excessm );
                 thlxsat = thln + xsat * ( thl_zm(kn) - thln );
                 qtxsat  = qtn  + xsat * ( qt_zm(kn) - qtn );
                 call condensation_mf(qtxsat, thlxsat, p_zm(kn), iexner_zm(kn), &
                                      thvxsat, qcn, thn, qln, qin, qsn, lmixn)
                 ! -------------------------------------------------- !
                 ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
                 ! -------------------------------------------------- !
                 do kk = 1, 2
                   if( kk .eq. 1 ) then
                     thv_x0 = thvn
                     thv_x1 = ( 1._kind_phys - 1._kind_phys/xsat ) * thvn + ( 1._kind_phys/xsat ) * thvxsat
                   else
                     thv_x1 = thv_zm(kn)
                     thv_x0 = ( xsat / ( xsat - 1._kind_phys ) ) * thv_zm(kn) + ( 1._kind_phys/( 1._kind_phys - xsat ) ) * thvxsat
                   endif
                   aquad =  wn**2
                   bquad =  2._kind_phys*wa*gravit*cridis*(thv_x1 - thv_x0)/thv_zm(kn) - 2._kind_phys*wn**2
                   cquad =  2._kind_phys*wa*gravit*cridis*(thv_x0 - thv_zm(kn))/thv_zm(kn)  + wn**2
                   if( kk .eq. 1 ) then
                     if( ( bquad**2-4._kind_phys*aquad*cquad ) .ge. 0._kind_phys ) then
                       call roots(aquad,bquad,cquad,xs1,xs2,status)
                       x_cu = min(1._kind_phys,max(0._kind_phys,min(xsat,min(xs1,xs2))))
                     else
                       x_cu = xsat
                     endif
                   else
                     if( ( bquad**2-4._kind_phys*aquad*cquad) .ge. 0._kind_phys ) then
                       call roots(aquad,bquad,cquad,xs1,xs2,status)
                       x_en = min(1._kind_phys,max(0._kind_phys,max(xsat,min(xs1,xs2))))
                     else
                       x_en = 1._kind_phys
                     endif
                   endif
                 enddo
                 if( x_cu .eq. xsat ) then
                   xc = max(x_cu, x_en)
                 else
                   xc = x_cu
                 endif
               endif

               ee2 = xc**2
               ud2 = 1._kind_phys - 2._kind_phys*xc + xc**2

               ! detrainment rate
               detn  = mix(kt,i) * ud2

             else !no bsort
               ee2 = 1._kind_phys
               ud2 = 1._kind_phys
             end if

             ! entrainment rate
             entn = mix(kt,i) * ee2

             ! --------------------------------------------------------- !
             ! TKE enhanced entrainment                                  !
             ! switches off when dynamic_L0 > max_L0                     !
             ! --------------------------------------------------------- !
             eturb = (1._kind_phys + clubb_mf_alphturb*sqrt(tke(k))/upw(k,i))
             if (do_clubb_mf_rhtke) then
               rh_L0 = 50._kind_phys*(rhinv**3._kind_phys)
               if (rh_L0 >= 733.34_kind_phys) eturb = 1._kind_phys
             else
               if (dynamic_L0(i) >= clubb_mf_max_L0) eturb = 1._kind_phys
             end if
             entn = entn * eturb

             ! integrate updraft
             entexp  = exp(-entn*eturb*dzt(kt))
             entexpu = exp(-entn*dzt(kt)/3._kind_phys)

             qtn  = qt(kt) *(1._kind_phys-entexp ) + upqt (k,i)*entexp + supqt(kt,i)
             thln = thl(kt)*(1._kind_phys-entexp ) + upthl(k,i)*entexp + supthl(kt,i)
             un   = u(kt)  *(1._kind_phys-entexpu) + upu  (k,i)*entexpu
             vn   = v(kt)  *(1._kind_phys-entexpu) + upv  (k,i)*entexpu

             ! convert source terms to a tendency (convert from S*dz/w to S)
             supqt(kt,i) = supqt(kt,i)*upw(k,i)/dzt(kt)
             upauto(kt,i) = supqt(kt,i)
             supthl(kt,i) = supthl(kt,i)*upw(k,i)/dzt(kt)

             ! get cloud, momentum levels
             if (do_condensation) then
               call condensation_mf(qtn, thln, p_zm(kn), iexner_zm(kn), &
                                    thvn, qcn, thn, qln, qin, qsn, lmixn)
               if (zcb(i).eq.zcb_unset .and. qcn > 0._kind_phys) zcb(i) = zm(kn)
             else
               thvn = thln*(1._kind_phys+zvir*qtn)
             end if

             ! get buoyancy
             B=gravit*(0.5_kind_phys*(thvn + upthv(k,i))/thv(kt)-1._kind_phys)

             if (do_implicit) then
               wp = clubb_mf_alphturb*wb*entn*sqrt(0.5_kind_phys*(tke(kn)+tke(k)))*dzt(kt)
               wn = (-wp + sqrt(wp**2._kind_phys + (1._kind_phys + 2._kind_phys*wb*entn*dzt(kt))* &
                     (upw(k,i)**2._kind_phys + 2._kind_phys*wa*B*dzt(kt))) )/(1._kind_phys + 2._kind_phys*wb*entn*dzt(kt))
             else
               ! get wn2
               wp = wb*entn*eturb
               if (wp==0._kind_phys) then
                 wn2 = upw(k,i)**2._kind_phys+2._kind_phys*wa*B*dzt(kt)
               else
                 entw = exp(-2._kind_phys*wp*dzt(kt))
                 wn2 = entw*upw(k,i)**2._kind_phys+(1._kind_phys-entw)*wa*B/wp
               end if
               wn = sqrt(max(wn2, 0._kind_phys))
             end if

           end do !iter_xc

           if (wn>0._kind_phys) then

             upthv(kn,i) = thvn
             upthl(kn,i) = thln
             upqt(kn,i)  = qtn
             upqc(kn,i)  = qcn
             upqs(kn,i)  = qsn
             upu(kn,i)   = un
             upv(kn,i)   = vn
             upql(kn,i)  = qln
             upqi(kn,i)  = qin
             upqv(kn,i)  = qtn - qcn
             uplmix(kn,i)= lmixn
             upth(kn,i)  = thn

             if (bsort) then
               mfn = upmf(k,i)*exp( dzt(kt)*( entn - detn ))
               upa(kn,i) = mfn/(wn*rho_zm(kn))
             else
               upa(kn,i) = upa(k,i)
               mfn = rho_zm(kn)*upa(kn,i)*wn
               detn = entn - (mfn - rho_zm(k)*upa(k,i)*upw(k,i)) &
                             /(rho_zm(k)*upa(k,i)*upw(k,i)*dzt(kt))
             end if

             upbuoy(kn,i)= B
             upw(kn,i)   = wn
             upmf(kn,i)  = mfn
             upent(kn,i) = entn
             updet(kn,i) = detn

           else
             ! zero out plumes that terminate at k<3
             if (abs(k-kstart+kdir)<4) then
               supqt(:,i) = 0._kind_phys
               upauto(:,i)= 0._kind_phys
               supthl(:,i)= 0._kind_phys
               upa(:,i)   = 0._kind_phys
               upbuoy(:,i)= 0._kind_phys
               upw(:,i)   = 0._kind_phys
               upmf(:,i)  = 0._kind_phys
               upent(:,i) = 0._kind_phys
               updet(:,i) = 0._kind_phys
               upthv(:,i) = 0._kind_phys
               upthl(:,i) = 0._kind_phys
               upqt(:,i)  = 0._kind_phys
               upqc(:,i)  = 0._kind_phys
               upqs(:,i)  = 0._kind_phys
               upu(:,i)   = 0._kind_phys
               upv(:,i)   = 0._kind_phys
               upql(:,i)  = 0._kind_phys
               upqi(:,i)  = 0._kind_phys
               upqv(:,i)  = 0._kind_phys
               uplmix(:,i)= 0._kind_phys
               upth(:,i)  = 0._kind_phys
             end if
             ! exit updraft integration
             exit
           end if
         enddo
       enddo

       ! ---------------------------------------------------------- !
       ! In-scheme production limiter                               !
       ! Cap the area-weighted plume drying at each level so the    !
       ! grid-mean microphysical sink cannot remove more than       !
       ! mf_dry_frac_max of the grid-cell total water in a timestep.!
       ! Because the cap is applied BEFORE the rain/evaporation     !
       ! sweeps and the downdraft section, all downstream rain      !
       ! accounting is consistent with the limited production. A    !
       ! common per-level factor lamk preserves the relative        !
       ! weighting of the plumes in the ensemble. The limiter       !
       ! establishes the invariant -sqt*dtime <= mf_dry_frac_max*qt !
       ! (asserted in clubb_intr.F90); the margin below 1.0 leaves  !
       ! room for other same-step moisture sinks and puts the       !
       ! assertion safely beyond roundoff.                          !
       ! ---------------------------------------------------------- !
       if (do_clubb_mf_precip) then
         do k = ksfcm+kdir, ktopm, kdir
           kt = k - (1+kdir)/2
           kn = k - kdir
           drytot = 0._kind_phys
           do i=1,clubb_mf_nup
             drytot = drytot + upa(kn,i)*max(-1._kind_phys*supqt(kt,i), 0._kind_phys)
           end do
           if (drytot*dtime > mf_dry_frac_max*max(qt(kt),0._kind_phys)) then
             lamk = mf_dry_frac_max*max(qt(kt),0._kind_phys)/(drytot*dtime)
             do i=1,clubb_mf_nup
               if (supqt(kt,i) < 0._kind_phys) then
                 supqt(kt,i)  = lamk*supqt(kt,i)
                 supthl(kt,i) = lamk*supthl(kt,i)
                 upauto(kt,i) = lamk*upauto(kt,i)
               end if
             end do
           end if
         end do
       end if

       ! --------------------------------------------------------- !
       ! downward sweep for rain evaporation, snow melting         !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip) then
         do i=1,clubb_mf_nup
           do k = ktopm, ksfcm+kdir, -kdir
             kt = k - (1+kdir)/2
             kn = k - kdir

             ! get rain evaporation
             if ((upqs(k,i) + upqs(kn,i)).le.0._kind_phys) then
               qtovqs = 0._kind_phys
             else
               qtovqs = (upqt(k,i) + upqt(kn,i))/(upqs(k,i) + upqs(kn,i))
             end if
             qtovqs = min(1._kind_phys,qtovqs)
             sevap = ke*(1._kind_phys - qtovqs)*sqrt(max(uprr(k,i),0._kind_phys))

             ! limit evaporation to available precip
             sevap = min(sevap,( uprr(k,i)/(rho_zt(kt)*dzt(kt)) - supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) ))

             ! further limit the area-weighted evaporation to the
             ! area-weighted (grid-mean) rain supply. Weight is upa
             ! at the interface below the cell, matching the sqtup accumulation
             if (upa(kn,i) > 0._kind_phys) then
               sevap = min(sevap, uprg(k,i)/(rho_zt(kt)*dzt(kt)*upa(kn,i)) &
                                  - supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) )
             end if
             sevap = max(sevap, 0._kind_phys)

             ! get rain rate
             uprr(kn,i) = uprr(k,i) &
                         - rho_zt(kt)*dzt(kt)*( supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) + sevap )

             ! grid-mean rain reservoir (max() only guards roundoff;
             ! the cap above keeps this non-negative by construction)
             uprg(kn,i) = max( uprg(k,i) &
                         - rho_zt(kt)*dzt(kt)*upa(kn,i)*( supqt(kt,i)*(1._kind_phys-clubb_mf_fdd) + sevap ), 0._kind_phys )

             ! update source terms
             lmixt = 0.5_kind_phys*(uplmix(k,i)+uplmix(kn,i))
             supqt(kt,i) = supqt(kt,i) + sevap
             supthl(kt,i) = supthl(kt,i) - lmixt*sevap*iexner_zt(kt)/cpair
           end do
         end do
       end if

       ! --------------------------------------------------------- !
       ! begin computing downdrafts                                !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip .and. clubb_mf_fdd > 0._kind_phys) then

         do i=1,clubb_mf_nup

           ! find cloud base
           do k = ksfcm, ktopm, kdir
             if (upqc(k,i) > 0._kind_phys) then
               ddbotm(i) = k
               exit
             end if
           end do

           ! find cloud top
           ddtopm = 0
           do k = ksfcm, ktopm, kdir
             if (uprr(k,i) > 0._kind_phys) ddtopm = k
           end do

           if (ddtopm /= 0) then
             ! initilaize downdrafts

             ! Kay initializes using negative of the updraft velocity
             ! this causes anomalouly large downdrafts at the initializaiton level
             ! I am intializing with zero velocity as that is more physically defensible
             dnw(ddtopm,i)   = -1._kind_phys*mindnw
             dna(ddtopm,i)   = upa(ddtopm,i)
             dnu(ddtopm,i)   = 0.5_kind_phys*(u(ddtopm - (1-kdir)/2)+u(ddtopm - (1+kdir)/2))
             dnv(ddtopm,i)   = 0.5_kind_phys*(v(ddtopm - (1-kdir)/2)+v(ddtopm - (1+kdir)/2))
             dnqt(ddtopm,i)  = qt_zm(ddtopm)

             ! no cloud in downdrafts, set to cloud free thl
             dnthl(ddtopm,i) = thl_zm(ddtopm)
             dnthv(ddtopm,i) = thv_zm(ddtopm) ! includes condensate loading (!)

             ! get rain generated in the updraft, appropriate it to the downdraft
             ! (we are using an upward sweep notation to be consistent with how uprr(ddtop) was computed)
             dnrr(ddtopm,i)  = -1._kind_phys*dzt(ddtopm - (1-kdir)/2)*rho_zt(ddtopm - (1-kdir)/2)*upauto(ddtopm - (1-kdir)/2,i)*clubb_mf_fdd

             ! grid-mean downdraft rain reservoir seed. area weight is
             ! upa at the interface below the generating cell, which is ddtopm
             ! itself (matching the sqtup accumulation weighting of the
             ! production that this rain came from).
             dnrg(ddtopm,i)  = -1._kind_phys*dzt(ddtopm - (1-kdir)/2)*rho_zt(ddtopm - (1-kdir)/2) &
                               *upa(ddtopm,i)*upauto(ddtopm - (1-kdir)/2,i)*clubb_mf_fdd

             if (fixent) then
               entn = fixent_ent
             else
               ! use deterministic mean entrainment
               entn = clubb_mf_ent0/dynamic_L0(i)
             end if

             ! downdraft qsat
             call qsat(dnthl(ddtopm,i)/iexner_zm(ddtopm),p_zm(ddtopm),es,dnqs(ddtopm,i))

             do k = ddtopm, ksfcm+kdir, -kdir

               kt = k - (1+kdir)/2
               kn = k - kdir

               ! assume fixed area with height
               dna(kn,i) = dna(k,i)

               ! get rain evaporation in integrated form
               taum1 = ke*sqrt(dnrr(k,i))/dnqs(k,i)
               alphint = exp(dzt(kt)*taum1/dnw(k,i))
               sqtint = max( (dnqs(k,i) - dnqt(k,i))*(1._kind_phys - alphint) ,0._kind_phys)

               ! limit to available rain
               sqtint = min( sqtint, -1._kind_phys*dnrr(k,i) / (rho_zt(kt)*dzt(kt)*dnw(k,i)) )
               sthlint = -1._kind_phys*latvap*sqtint*iexner_zt(kt)/cpair

               ! get rain evaporation in tendency form
               sdnqt(kt,i) = max( (dnqs(k,i) - dnqt(k,i))*taum1, 0._kind_phys )

               ! cap the tendency-form downdraft evaporation by the
               ! area-weighted rain actually available to the downdraft.
               ! weights: dna at the interface above the cell
               ! (matching the sqtdn accumulation), upa at the interface below
               ! (matching sqtup).
               if (dna(k,i) > 0._kind_phys) then
                 sdnmax = ( dnrg(k,i)/(rho_zt(kt)*dzt(kt)) &
                            - upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ) / dna(k,i)
                 sdnqt(kt,i) = min( sdnqt(kt,i), max(sdnmax, 0._kind_phys) )
               end if
               sdnthl(kt,i) = -1._kind_phys*latvap*sdnqt(kt,i)*iexner_zt(kt)/cpair

               ! compute rain rate (rain above - evaporation + appropriate updraft rain)
               dnrr(kn,i) = max( dnrr(k,i) &
                                - rho_zt(kt)*dzt(kt)*(sdnqt(kt,i) + upauto(kt,i)*clubb_mf_fdd) , 0._kind_phys )

               ! include eturb?
               entexp  = exp(-1._kind_phys*entn*eturb*dzt(kt))
               entexpu = exp(-1._kind_phys*entn*dzt(kt)/3._kind_phys)

               ! integrate downward
               dnu(kn,i)   = u(kt)  *(1._kind_phys-entexpu) + dnu  (k,i)*entexpu
               dnv(kn,i)   = v(kt)  *(1._kind_phys-entexpu) + dnv  (k,i)*entexpu
               dnqt(kn,i)  = qt(kt) *(1._kind_phys-entexp ) + dnqt (k,i)*entexp + sqtint
               dnthl(kn,i) = thl(kt)*(1._kind_phys-entexp ) + dnthl(k,i)*entexp + sthlint

               ! get qsat
               call qsat(dnthl(kn,i)/iexner_zm(kn),p_zm(kn),es,dnqs(kn,i))

               ! no supersaturation in downdrafts
               if (dnqt(kn,i) > dnqs(kn,i)) then
                 ! set qt to saturation vapor pressure
                 dnqt(kn,i) = dnqs(kn,i)

                 ! find evaporation that gives saturation vapor pressure
                 sqtint = dnqt(kn,i) - (qt(kt) *(1._kind_phys-entexp ) + dnqt (k,i)*entexp)

                 ! limit to available rain
                 sqtint = min( sqtint, -1._kind_phys*dnrr(k,i) / (rho_zt(kt)*dzt(kt)*dnw(k,i)) )
                 sthlint = -1._kind_phys*latvap*sqtint*iexner_zt(kt)/cpair

                 ! find new evap tendency
                 if ((alphint - 1._kind_phys) /= 0._kind_phys) then
                   qtmp = dnqs(k,i) + sqtint/(alphint - 1._kind_phys)
                   sdnqt(kt,i) = max( (dnqs(k,i) - qtmp)*taum1, 0._kind_phys )
                 else
                   sdnqt(kt,i) = 0._kind_phys
                 end if

                 ! apply the same area-weighted available-rain cap to
                 ! the recomputed (saturation-adjusted) evaporation tendency
                 if (dna(k,i) > 0._kind_phys) then
                   sdnmax = ( dnrg(k,i)/(rho_zt(kt)*dzt(kt)) &
                              - upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ) / dna(k,i)
                   sdnqt(kt,i) = min( sdnqt(kt,i), max(sdnmax, 0._kind_phys) )
                 end if
                 sdnthl(kt,i) = -1._kind_phys*latvap*sdnqt(kt,i)*iexner_zt(kt)/cpair

                 ! re-compute thl with new evaporation rate
                 dnthl(kn,i) = thl(kt)*(1._kind_phys-entexp ) + dnthl(k,i)*entexp + sthlint

                 ! adjust rain
                 dnrr(kn,i) = max( dnrr(k,i) &
                                  - rho_zt(kt)*dzt(kt)*(sdnqt(kt,i) + upauto(kt,i)*clubb_mf_fdd) , 0._kind_phys )
               end if

               ! grid-mean downdraft rain reservoir, using the final
               ! (capped) evaporation tendency. max() only guards roundoff; the
               ! sdnqt cap keeps this non-negative by construction.
               dnrg(kn,i) = max( dnrg(k,i) &
                    - rho_zt(kt)*dzt(kt)*( dna(k,i)*sdnqt(kt,i) &
                                           + upa(kn,i)*upauto(kt,i)*clubb_mf_fdd ), 0._kind_phys )

               ! get virtual temperature
               dnthv(kn,i) = dnthl(kn,i)*(1._kind_phys+zvir*dnqt(kn,i))

               if ((kn - ddbotm(i))*kdir > 0) then
                 ! get virtual temperature
                 dnthv(kn,i) = dnthl(kn,i)*(1._kind_phys+zvir*dnqt(kn,i))

                 ! get buoyancy
                 ! (midpoint k is surrounded by interface k and k-1,
                 ! and therefore we can't compute B at the midpoint properly)
                 B = gravit*(dnthv(kn,i)/thv(kt)-1._kind_phys)

                 ! get wn2
                 wp = wb*entn*eturb  &
                      + clubb_mf_pwfac/( 2._kind_phys*zm(kn + kdir)+tinynum ) * max( 1._kind_phys - exp( zm(kn + kdir)/z00dn-1._kind_phys), 0._kind_phys )
                 if (wp==0._kind_phys) then
                   wn2 = dnw(k,i)**2._kind_phys-2._kind_phys*wa*B*dzt(kt)
                 else
                   entw = exp(-2._kind_phys*wp*dzt(kt))
                   wn2 = entw*dnw(k,i)**2._kind_phys-(1._kind_phys-entw)*wa*B/wp
                 end if
                 wn2 = max(wn2,mindnw**2._kind_phys)
                 dnw(kn,i) = -1._kind_phys*sqrt(wn2)

               else
                 zsub = zm(ddbotm(i)+kdir)
                 wcb  = dnw(ddbotm(i)+kdir,i)
                 dnw(kn,i) = wcb - (wcb/(zsub**clubb_mf_ddexp))*(zsub - zm(kn))**clubb_mf_ddexp
                 dnw(kn,i) = min(dnw(kn,i),-1._kind_phys*mindnw)
               end if

             end do!k

           end if

         end do!i

!+++ARH this should be changed to only zero out above the downdraft (dnw<-mindw)
!+++ARH also this should zero out dna as well
         ! zero out downdraft fluxes for dnw == -mindnw
         do i=1,clubb_mf_nup
           do k=ksfcm, ktopm, kdir
             if ( dnw(k,i) == -1._kind_phys*mindnw ) then
               dnw(k,i) = 0._kind_phys
               dna(k,i) = 0._kind_phys
             end if
           end do
         end do

       end if
       ! end computing downdrafts

       ! --------------------------------------------------------- !
       ! AS.pd limiter                                             !
       ! --------------------------------------------------------- !
       if (do_aspd) then
         do k = ksfcm, ktopm-kdir, kdir
           kt = k - (1-kdir)/2
           kn = k + kdir
           do i=1,clubb_mf_nup
             if (upw(kn,i)>0._kind_phys) then
               ! diagnose detrainment
               Mn = rho_zm(k)*upa(k,i)*upw(k,i)
               det = upent(kn,i)*eturb - (rho_zm(kn)*upa(kn,i)*upw(kn,i) - Mn) &
                             /(Mn*dzt(kt))
               if (det < 0._kind_phys) then
                 ! diagnose area to eliminate detrainment and conserve mass
                 Mn = rho_zm(k)*upa(k,i)*upw(k,i)*exp(upent(kn,i)*eturb*dzt(kt))
                 upa(kn,i) = Mn/(rho_zm(kn)*upw(kn,i))
               end if
             end if
           end do
         end do
       end if

       ! --------------------------------------------------------- !
       ! integrate for total convective area                       !
       ! --------------------------------------------------------- !
       do k=ksfcm, ktopm-kdir, kdir
         do i=1,clubb_mf_nup
           aup(k) = aup(k) + upa(k,i)
           adn(k) = adn(k) + dna(k,i)
         end do
         ac(k) = aup(k) + adn(k)
         if (limarea .and. ac(k) > amax) then
           upa(k,:) = upa(k,:)*amax/ac(k)
           ac(k) = amax
         end if
         ae(k) = ae(k) - ac(k)
       end do

       ! --------------------------------------------------------- !
       ! updraft properties for output                             !
       ! --------------------------------------------------------- !
       do k=ksfcm, ktopm, kdir

         ! first sum over all i-updrafts
         do i=1,clubb_mf_nup
           if (upqc(k,i)>0._kind_phys) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0._kind_phys ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0._kind_phys
           dry_qt(k)  = 0._kind_phys
           dry_thl(k) = 0._kind_phys
           dry_u(k)   = 0._kind_phys
           dry_v(k)   = 0._kind_phys
         endif

         if ( moist_a(k) > 0._kind_phys ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0._kind_phys
           moist_qt(k)  = 0._kind_phys
           moist_thl(k) = 0._kind_phys
           moist_u(k)   = 0._kind_phys
           moist_v(k)   = 0._kind_phys
           moist_qc(k)  = 0._kind_phys
         endif

       enddo
       ! --------------------------------------------------------- !
       ! get ensemble mean                                         !
       ! --------------------------------------------------------- !

       ! 1. Momentum Grid Accumulations (Interfaces)
       ! Iterates over all interfaces (1 to 59)
       do k=ksfcm, ktopm, kdir
         do i=1,clubb_mf_nup
           awup(k) = awup(k) + upa(k,i)*upw(k,i)
           awdn(k) = awdn(k) + dna(k,i)*dnw(k,i)

           awwup(k) = awwup(k) + upa(k,i)*upw(k,i)*upw(k,i)
           awwdn(k) = awwdn(k) + dna(k,i)*dnw(k,i)*dnw(k,i)

           awuup(k) = awuup(k) + upa(k,i)*upw(k,i)*upu(k,i)
           awudn(k) = awudn(k) + dna(k,i)*dnw(k,i)*dnu(k,i)
           awvup(k) = awvup(k) + upa(k,i)*upw(k,i)*upv(k,i)
           awvdn(k) = awvdn(k) + dna(k,i)*dnw(k,i)*dnv(k,i)

           awthvdn(k)= awthvdn(k)+ dna(k,i)*dnw(k,i)*dnthv(k,i)
           awthldn(k)= awthldn(k)+ dna(k,i)*dnw(k,i)*dnthl(k,i)
           awqtdn(k) = awqtdn(k) + dna(k,i)*dnw(k,i)*dnqt(k,i)

           awthvup(k)= awthvup(k)+ upa(k,i)*upw(k,i)*upthv(k,i)
           awthlup(k)= awthlup(k)+ upa(k,i)*upw(k,i)*upthl(k,i)
           awqtup(k) = awqtup(k) + upa(k,i)*upw(k,i)*upqt(k,i)
         enddo

         aw (k) = awup(k)+ awdn(k)
         aww(k) = awwup(k)+ awwdn(k)
         awu(k) = awuup(k)+ awudn(k)
         awv(k) = awvup(k)+ awvdn(k)
       enddo

       ! 2. Thermodynamic Grid Accumulations (Cells)
       ! Iterates over active atmospheric cells only
       ! Top-Down: 58 down to 1. Bottom-Up: 2 up to nzm.
       do k = ksfcm+kdir, ktopm, kdir
         ! The interface directly BELOW the current cell
         ! In top-down (-1): cell 58 -> interface 59
         ! In bottom-up (1): cell 2 -> interface 1
         kn = k - kdir
         kt_dn = k - (1+kdir)/2

         do i=1,clubb_mf_nup
           sqtup(kt_dn)  = sqtup(kt_dn)  + upa(kn,i)*supqt(kt_dn,i)
           sthlup(kt_dn) = sthlup(kt_dn) + upa(kn,i)*supthl(kt_dn,i)

           sqtdn(kt_dn)  = sqtdn(kt_dn)  + dna(k,i)*sdnqt(kt_dn,i)
           sthldn(kt_dn) = sthldn(kt_dn) + dna(k,i)*sdnthl(kt_dn,i)
         end do

         sqt(kt_dn)  = sqtup(kt_dn)  + sqtdn(kt_dn)
         sthl(kt_dn) = sthlup(kt_dn) + sthldn(kt_dn)
       enddo
       ! --------------------------------------------------------- !
       ! ztopm1 calculation                                        !
       ! --------------------------------------------------------- !
       do i=1,clubb_mf_nup
         do k=kstart, ktopm, kdir
           ! return if no convection at first level above surface
           if (k == (ksfcm+kdir) .and. ac(k) == 0._kind_phys .and. .not.aloft) then
             sqt(k) = 0._kind_phys
             sthl(k) = 0._kind_phys
             ztopm1(:) = zm(ksfcm)
             ddcp(:) = 0._kind_phys
             return
           end if
           ! height of the plume ensemble
           if (do_clubb_mf_lscale_perplume) then
             if ((upa(k,i)+dna(k,i)) > 0._kind_phys) ztopm1(i) = zm(k)
           else
             if (ac(k) > 0._kind_phys) ztopm1(:) = zm(k)
           end if
         end do
       end do

       !subtract init level from ztop for aloft plumes
       if (aloft) then
         ztopm1 = ztopm1 - zm(kstart-nbot*kdir)
         do i=1,clubb_mf_nup
           if (ztopm1(i) < zm(ksfcm)) ztopm1(i) = zm(ksfcm)
         end do
       end if

       ! --------------------------------------------------------- !
       ! cloud base / mixing depth calculation                     !
       ! --------------------------------------------------------- !
       cbm1 = 0._kind_phys
       do i=1,clubb_mf_nup
         kcbarr(i) = 0
         do k=kstart, ktopm, kdir
           if (upqc(k,i) > 0._kind_phys) then
             kcbarr(i) = k
             exit
           end if
         end do

         ! find height of dry plumes
         if (kcbarr(i) == 0) then
           do k=kstart, ktopm, kdir
             if (upw(k,i) <= 0._kind_phys) then
               kcbarr(i) = k
               exit
             end if
           end do
         end if

         cbm1 = cbm1 + zm(kcbarr(i))

       end do
       cbm1 = cbm1/REAL(clubb_mf_nup)

       ! --------------------------------------------------------- !
       ! bulk downdraft velocity for coldpool parameterization     !
       ! --------------------------------------------------------- !
!+++ARH
!       ! reset ddcp
!       ddcp = 0._kind_phys
!       do i=1,clubb_mf_nup
!         ! find cloud base
!         kcb = 0
!         do k=1,nz
!           if (upqc(k,i) > 0._kind_phys) then
!             kcb = k
!             exit
!           end if
!         end do
!
!         ! reset iddcp
!         iddcp = 0._kind_phys
!         if (kcb == 0) then
!           continue
!         else if (kcb == 1) then
!           iddcp = iddcp + dna(k,i)*dnw(k,i)
!           continue
!         else
!           ddint = 0._kind_phys
!           do k=1,kcb-1
!             ddint = ddint + dna(k,i)*dnw(k,i)*dzt(k+1)
!           end do
!           iddcp = iddcp + -1._kind_phys*ddint/zm(kcb)
!         end if
!         ddcp = ddcp + iddcp
!         !
!       end do
!

       ddcp(:) = 0._kind_phys
       if (do_clubb_mf_coldpool .and. clubb_mf_fdd > 0._kind_phys) then
         ! use single level for cold pool param.
         ! reset ddcp
         do i=1,clubb_mf_nup
           if (ddbotm(i) == 0) then
             continue
           else
             if (do_clubb_mf_coldpool_perplume) then
               ddcp(i) = -1._kind_phys*dnw(ddbotm(i)+kdir,i)
             else
               ddcp(:) = -1._kind_phys*dna(ddbotm(i)+kdir,i)*dnw(ddbotm(i)+kdir,i) + ddcp(:)
             end if
           end if
         end do
       end if

       ! --------------------------------------------------------- !
       ! downward sweep to get ensemble mean precip                !
       ! --------------------------------------------------------- !
       do k = ktopm, ksfcm+kdir, -kdir
         kt_dn = k - (1+kdir)/2
         precc(k-kdir) = precc(k) - rho_zt(kt_dn)*dzt(kt_dn)*sqt(kt_dn)
       end do

       ! this clamp only removes roundoff (upstream limiters conserve mass)
       precc(ksfcm) = max( precc(ksfcm), 0._kind_phys )

       ! --------------------------------------------------------- !
       ! get turbulent fluxes                                      !
       ! --------------------------------------------------------- !
       kstart = ksfcm + kdir
       if (scalesrf) then
         kstart = ksfcm
       end if

       do k=kstart, ktopm-kdir, kdir

         ! Secure boundary cells to prevent array out-of-bounds on zero-flux boundaries
         kt_up = max(1, min(nzt, k - (1-kdir)/2))
         kt_dn = max(1, min(nzt, k - (1+kdir)/2))

         thvflxup(k)= awthvup(k) - awup(k)*thv(kt_up)
         thlflxup(k)= awthlup(k) - awup(k)*thl(kt_up)
         qtflxup (k)= awqtup (k) - awup(k)*qt(kt_up)

         uflxup  (k)= awuup(k) - awup(k)*u(kt_up)
         vflxup  (k)= awvup(k) - awup(k)*v(kt_up)

         thvflxdn(k)= awthvdn(k) - awdn(k)*thv(kt_dn)
         thlflxdn(k)= awthldn(k) - awdn(k)*thl(kt_dn)
         qtflxdn (k)= awqtdn (k) - awdn(k)*qt(kt_dn)

         uflxdn  (k)= awudn(k) - awdn(k)*u(kt_dn)
         vflxdn  (k)= awvdn(k) - awdn(k)*v(kt_dn)

         thvflx(k)  = thvflxup(k) + thvflxdn(k)
         thlflx(k)  = thlflxup(k) + thlflxdn(k)
         qtflx (k)  = qtflxup (k) + qtflxdn (k)

         uflx(k)    = uflxup(k) + uflxdn(k)
         vflx(k)    = vflxup(k) + vflxdn(k)
       enddo
     else
       ddcp(:) = 0._kind_phys
       ztopm1(:) = zm(ksfcm)
     end if

  end subroutine integrate_mf


  subroutine get_Lscale(nzt, nzm, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                        wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                        do_condensation, qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1, dynamic_L0, ztop, mcape)
  ! --------------------------------------------------------- !
  ! Calculate ztop and dynamic_L based on value of namelist   !
  ! --------------------------------------------------------- !
     integer,  intent(in)                 :: nzt, nzm
     real(kind_phys), dimension(nzt), intent(in) :: thl,    thv,          &
                                             th,                   &
                                             qt,     qv,           &
                                             p_zt,   iexner_zt,    &
                                             dzt,    zt

     real(kind_phys), dimension(nzm), intent(in) :: p_zm,   iexner_zm,    &
                                             zm,     tke,    wpthlp_env

     real(kind_phys), intent(in) ::                wmax,   wmin,       tpert, &
                                            sigmaw, sigmaqt, sigmathv, &
                                            cwqt,   cwthv,  zcb_unset, &
                                            wa,     wb,        ztopm1, &
                                            pblh,   convh,     rhinv

     logical, intent(in) ::                 do_condensation

     real(kind_phys), intent(out) ::               dynamic_L0, ztop, mcape

     ! local variables
     ! =============================================================================== !
     ! GRID ORIENTATION GENERALIZATION VARIABLES
     ! ------------------------------------------------------------------------------- !
     ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
     ! duplicating code, these variables abstract the vertical loop bounds and slices.
     !
     ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
     ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
     ! kdir          : Directional step (+1 for moving up, -1 for moving down).
     !
     ! STAGGERED GRID INDEXING
     ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
     ! index of the cell center (zt) to the interface (zm) flips depending on whether
     ! memory is loaded top-down or bottom-up. These variables dynamically map them:
     !
     ! kt    : The active thermodynamic cell center associated with the current step.
     !         Upward Sweep:   kt = k - (1-kdir)/2
     !         Downward Sweep: kt = k - (1+kdir)/2
     !
     ! kn    : The NEXT momentum interface in the direction of the current sweep.
     !         Upward Sweep:   kn = k + kdir
     !         Downward Sweep: kn = k - kdir
     !
     ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
     !         kt_up = k - (1-kdir)/2
     !
     ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
     !         kt_dn = k - (1+kdir)/2
     ! =============================================================================== !
     integer :: ksfc, ktop, kdir


     real(kind_phys), dimension(nzt)               :: t_zt
     real(kind_phys), dimension(nzt)               :: tp,       qstp
     real(kind_phys), dimension(nzt,1)             :: dmpdz
     real(kind_phys), dimension(1)                 :: tl,                     &
                                               cape,     cin
     integer,  dimension(1)                 :: lcl,      lel
     real(kind_phys)                               :: landfrac
     integer                                :: kpbl,     msg,          &
                                               lon,      mx,           &
                                               k

     ! intialize local variables
     cape      = 0._kind_phys
     mcape     = 0._kind_phys
     ztop      = 0._kind_phys
     dmpdz     = 0._kind_phys

     if (zt(1) < zt(nzt)) then
        ksfc = 1
        ktop = nzt
        kdir = 1
     else
        ksfc = nzt
        ktop = 1
        kdir = -1
     end if

     if (clubb_mf_Lopt==0) then
       !Constant L0
       dynamic_L0 = clubb_mf_L0
       ztop = clubb_mf_L0
     else if (clubb_mf_Lopt==1) then
       !TKE
       do k = ktop-2*kdir, ksfc+kdir, -kdir
         if (zm(k) < 20000 .and. tke(k) - tke(k+kdir) > 1e-5) then
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==2) then
       !Heat flux
       do k = ktop-2*kdir, ksfc+kdir, -kdir
         if (zm(k) < 20000 .and. abs(abs(wpthlp_env(k))-abs(wpthlp_env(k-kdir))) > 1e-4) then
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==3) then
       !Test plume
       call oneplume( nzm, nzt, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                      wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                      wa, wb, tke, do_condensation, do_clubb_mf_precip, ztop )

       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
       !ztop = ztop - 1600._kind_phys
       !if (ztop < 1._kind_phys) then
       !  dynamic_L0 = clubb_mf_a0
       !else
       !  dynamic_L0 = min(35._kind_phys,clubb_mf_a0*(ztop**clubb_mf_b0))
       !end if
     else if (clubb_mf_Lopt==4 .or. clubb_mf_Lopt==5) then
       !dilute cape calculation
       !dmpdz = -1._kind_phys*ent_zt(2:nz,:)
       dmpdz(:,:) = -1.E-3_kind_phys
       t_zt = th/iexner_zt
       landfrac = 1._kind_phys

       do k = ksfc+kdir, ktop, kdir
         if (zt(k-kdir) <= pblh) then
           kpbl = k
         end if
       end do

       do k = ksfc, ktop, kdir
         if (p_zt(k) > 40.e2_kind_phys) then
           msg = k
         end if
       end do

       call buoyan_dilute(nzt, nzm, 1          ,dmpdz , &
                          qv         ,t_zt       ,p_zt*0.01_kind_phys       ,zt       ,p_zm*0.01_kind_phys , &
                          tp         ,qstp       ,tl         ,cape         ,cin  , &
                          kpbl-kdir  ,lcl        ,lel        ,lon      ,mx   , &
                          msg-kdir   ,tpert      ,landfrac )

       !do i=1,clubb_mf_nup
       !  mcape = mcape + cape(i)
       !end do
       !mcape = mcape/REAL(clubb_mf_nup)
       mcape = max(cape(1),25._kind_phys)

       if (clubb_mf_Lopt==4) then
         ztop = max(zt(lel(1)+kdir),convh)
       else if (clubb_mf_Lopt==5) then
         ztop = mcape
       end if
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==6) then
       ! grab ztop from max height of ensemble in prior time-step(s)
       ztop = ztopm1
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     else if (clubb_mf_Lopt==7 .or. clubb_mf_Lopt==8) then
       ztop = rhinv
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     end if

  end subroutine get_Lscale

  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi, qs, lmix )
  ! =============================================================================== !
  ! zero or one condensation for edmf: calculates thv and qc                        !
  ! =============================================================================== !
     use physconst,          only: cpair, zvir, h2otrip
     use wv_saturation,      only : qsat

     real(kind_phys),intent(in) :: qt,thl,p,iex
     real(kind_phys),intent(out):: thv,qc,th,ql,qi,qs,lmix

     !local variables
     integer  :: niter,i
     real(kind_phys) :: diff,t,qstmp,qcold,es,wf
     logical  :: noice = .true.
     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5_kind_phys

     qc=0._kind_phys
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter

       if (noice) then
         wf = 1._kind_phys
       else
         wf = get_watf(t)
       end if
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qstmp)
       qcold = qc
       qc = max(0.5_kind_phys*qc+0.5_kind_phys*(qt-qstmp),0._kind_phys)
       if (abs(qc-qcold)<diff) exit
     enddo

     if (noice) then
       wf = 1._kind_phys
     else
       wf = get_watf(t)
     end if
     t = thl/iex+get_alhl(wf)/cpair*qc

     call qsat(t,p,es,qs)
     qc = max(qt-qs,0._kind_phys)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1._kind_phys+zvir*(qt-qc)-qc)
     lmix = get_alhl(wf)
     th = t*iex
     qi = qc*(1._kind_phys-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(kind_phys)            :: t,get_watf,tc
       real(kind_phys), parameter :: &
                              tmax=-10._kind_phys, &
                              tmin=-40._kind_phys

       tc=t-h2otrip

       if (tc>tmax) then
         get_watf=1._kind_phys
       else if (tc<tmin) then
         get_watf=0._kind_phys
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice
       real(kind_phys) :: get_alhl,wf

       get_alhl = wf*latvap+(1._kind_phys-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf


  subroutine precip_mf(qs,qt,w,dz,dzcld,Supqt)
  !**********************************************************************
  ! Precipitation microphysics
  ! By Adam Herrington, after Kay Suselj
  !**********************************************************************

       real(kind_phys),intent(in)  :: qs,qt,w,dz,dzcld
       real(kind_phys),intent(out) :: Supqt
       ! ! local vars
       real(kind_phys)            :: tauwgt, tau,       & ! time-scale vars
                              qstar                ! excess cloud liquid

       real(kind_phys),parameter  :: tau0  = 15._kind_phys,    & ! base time-scale
                              zmin  = 300._kind_phys,   & ! small cloud thick
                              zmax  = 3000._kind_phys,  & ! large cloud thick
                              qcmin = 0.00125_kind_phys   ! supersat threshold

       qstar = qs+qcmin

       if (qt > qstar) then
         ! get precip efficiency
         tauwgt = (dzcld-zmin)/(zmax-zmin)
         tauwgt = min(max(tauwgt,0._kind_phys),1._kind_phys)
         tau    = tauwgt/tau0

         ! get source for updraft
         Supqt = (qstar-qt)*(1._kind_phys - exp(-1._kind_phys*tau*dz/w))
       else
         Supqt = 0._kind_phys
       end if

  end subroutine precip_mf

  subroutine poisson(nz, nup, ksfc, ktop, kdir, lambda, poi, state)
  !**********************************************************************
  ! Set a unique (but reproduceble) seed for the kiss RNG
  ! Call Poisson deviate
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       integer,                     intent(in)  :: nz, nup, ksfc, ktop, kdir
       real(kind_phys), dimension(4),      intent(in)  :: state
       real(kind_phys), dimension(nz,nup), intent(in)  :: lambda
       integer,  dimension(nz,nup), intent(out) :: poi
       integer,  dimension(1,4)                 :: tmpseed
       integer                                  :: i,j
       type(ShrKissRandGen)                     :: kiss_gen

       ! Compute seed
       tmpseed(1,1) = int((state(1) - int(state(1))) * 1000000000._kind_phys)
       tmpseed(1,2) = int((state(2) - int(state(2))) * 1000000000._kind_phys)
       tmpseed(1,3) = int((state(3) - int(state(3))) * 1000000000._kind_phys)
       tmpseed(1,4) = int((state(4) - int(state(4))) * 1000000000._kind_phys)

       ! Set seed
       kiss_gen = ShrKissRandGen(tmpseed)

       ! Loop from the SURFACE to the TOP to preserve the exact PRNG sequence
       do i = ksfc, ktop, kdir
         do j = 1, nup
           call hybridRNG(kiss_gen, lambda(i,j), poi(i,j))
         enddo
       enddo

  end subroutine poisson


  subroutine hybridRNG(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Interface for the two poisson rng subroutines
  ! chooses the appropriate subroutine based on the value of lambda
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(kind_phys),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       if (lambda < 10._kind_phys) then
          call knuth(kiss_gen,lambda,kout)
       else
          call hormann(kiss_gen,lambda,kout)
       end if

  end subroutine hybridRNG


  subroutine knuth(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(kind_phys),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       ! Local variables
       real(kind_phys), dimension(1,1) :: tmpuni
       real(kind_phys)                 :: puni, explam
       integer                  :: k

       k = 0
       explam = exp(-1._kind_phys*lambda)
       puni = 1._kind_phys
       do while (puni > explam)
         k = k + 1
         call kiss_gen%random(tmpuni)
         puni = puni*tmpuni(1,1)
       end do
       kout = k - 1

  end subroutine knuth


  subroutine hormann(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson
  ! Implements Poisson Transformed Rejection with Squeeze (PTRS)
  ! from W. Hormann Insurance: Mathematics and Economics 12, 39-45 (1993)
  ! By Jake Reschke
  !**********************************************************************
  use shr_RandNum_mod, only: ShrKissRandGen

      type(ShrKissRandGen), intent(inout) :: kiss_gen
      real(kind_phys),             intent(in)    :: lambda
      integer,              intent(out)   :: kout

      ! Local variables
      real(kind_phys), dimension(1,1) :: U,V
      real(kind_phys)                 :: a,b,vr,alphinv,us,loggam
      integer                  :: k,i

      b = 0.931_kind_phys + 2.53_kind_phys*sqrt(lambda)
      a = -0.059_kind_phys + 0.02483_kind_phys*b
      vr = 0.9277_kind_phys - 3.6224_kind_phys/(b - 2._kind_phys)
      alphinv = 1.1239_kind_phys + 1.1328_kind_phys/(b - 3.4_kind_phys)

      do
         call kiss_gen%random(U)
         call kiss_gen%random(V)
         U(1,1) = U(1,1) - 0.5_kind_phys
         us = 0.5_kind_phys - abs(U(1,1))
         k = floor( (2._kind_phys*a/us + b)*U(1,1) + lambda + 0.43_kind_phys )
         if (us >= 0.07_kind_phys .and.  V(1,1) <= vr) then
            kout = k
            exit
         end if
         if (k <= 0 .or. (us < 0.013_kind_phys .and. V(1,1) > us)) then
            cycle
         end if
         ! compute log(k!). If k >=10 use stirling's approximation
         if (k < 10) then
            loggam = 0._kind_phys
            do i = 1, k
               loggam = loggam + log(1._kind_phys*i)
            end do
         else
            loggam = log(sqrt(2._kind_phys*pi)) + (k + 0.5_kind_phys)*log(1._kind_phys*k) - k + (1._kind_phys/12._kind_phys - 1._kind_phys/(360._kind_phys*k*k))/k
         end if
         if (log( V(1,1)*alphinv/(a/(us*us) + b) ) <= -1._kind_phys*lambda + k*log(lambda) - loggam) then
            kout = k
            exit
         end if
      end do

  end subroutine hormann

  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! after uwshcu.F90                                          !
  ! --------------------------------------------------------- !
    real(kind_phys), intent(in)  :: a
    real(kind_phys), intent(in)  :: b
    real(kind_phys), intent(in)  :: c
    real(kind_phys), intent(out) :: r1
    real(kind_phys), intent(out) :: r2
    integer , intent(out) :: status
    real(kind_phys)              :: q

    status = 0

    if( a .eq. 0._kind_phys ) then                     ! Form b*x + c = 0
        if( b .eq. 0._kind_phys ) then                        ! Failure: c = 0
            status = 1
        else                                           ! b*x + c = 0
            r1 = -c/b
        endif
        r2 = r1
    else
        if( b .eq. 0._kind_phys ) then                        ! Form a*x**2 + c = 0
            if( a*c .gt. 0._kind_phys ) then                  ! Failure: x**2 = -c/a < 0
                status = 2
            else                                       ! x**2 = -c/a
                r1 = sqrt(-c/a)
            endif
            r2 = -r1
       else                                            ! Form a*x**2 + b*x + c = 0
            if( (b**2 - 4._kind_phys*a*c) .lt. 0._kind_phys ) then   ! Failure, no real roots
                 status = 3
            else
                 q  = -0.5_kind_phys*(b + sign(1.0_kind_phys,b)*sqrt(b**2 - 4._kind_phys*a*c))
                 r1 =  q/a
                 r2 =  c/q
            endif
       endif
    endif

    return

  end subroutine roots


  subroutine oneplume( nzm, nzt, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                       wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                       wa, wb, tke, do_condensation, do_precip, plumeheight )
  !**********************************************************************
  ! Calculate a single plume with fixed entrainment
  ! to be used for a dynamic mixing length calculation
  ! By Rachel Storer
  !**********************************************************************
    use physconst,          only: cpair, gravit, zvir

    integer,  intent(in)                 :: nzm, nzt
    real(kind_phys), dimension(nzt), intent(in) :: dzt, iexner_zt,  &
                                            qt, thv, thl
    real(kind_phys), dimension(nzm), intent(in) :: zm, p_zm, iexner_zm, tke

    real(kind_phys), intent(in)                :: wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, &
                                           cwthv, zcb_unset, wa, wb
    logical, intent(in)               :: do_condensation, do_precip

    real(kind_phys), intent(inout)             :: plumeheight

    !local variables

    ! =============================================================================== !
    ! GRID ORIENTATION GENERALIZATION VARIABLES
    ! ------------------------------------------------------------------------------- !
    ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
    ! duplicating code, these variables abstract the vertical loop bounds and slices.
    !
    ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
    ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
    ! kdir          : Directional step (+1 for moving up, -1 for moving down).
    !
    ! STAGGERED GRID INDEXING
    ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
    ! index of the cell center (zt) to the interface (zm) flips depending on whether
    ! memory is loaded top-down or bottom-up. These variables dynamically map them:
    !
    ! kt    : The active thermodynamic cell center associated with the current step.
    !         Upward Sweep:   kt = k - (1-kdir)/2
    !         Downward Sweep: kt = k - (1+kdir)/2
    !
    ! kn    : The NEXT momentum interface in the direction of the current sweep.
    !         Upward Sweep:   kn = k + kdir
    !         Downward Sweep: kn = k - kdir
    !
    ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
    !         kt_up = k - (1-kdir)/2
    !
    ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
    !         kt_dn = k - (1+kdir)/2
    ! =============================================================================== !
    integer :: ksfcm, ktopm, kdir, kt, kn
    integer                     :: k
    real(kind_phys)                    :: thvn, qtn, thln, qcn, thn, qln, qin, qsn, lmixn, zcb, B, wn2, pentexp, pturb, pentw, wp
    real(kind_phys), dimension(nzm)     :: upw, upa, upqt, upthv, upthl, upth, upqs, &
                                   upqc, upql, upqi
    real(kind_phys), dimension(nzt)     :: supqt, supthl
    ! ! fractional entrainment rate
    real(kind_phys), parameter         :: pent = 1.E-3_kind_phys
    ! ! use tke enhanced entrainment
    logical                     :: do_tptke = .false.

    if (zm(1) < zm(nzm)) then
       ksfcm = 1
       ktopm = nzm
       kdir = 1
    else
       ksfcm = nzm
       ktopm = 1
       kdir = -1
    end if

    zcb = zcb_unset

    upw(ksfcm) = 0.5_kind_phys * wmax
    upa(ksfcm) = 0.5_kind_phys * erf( wmax/(sqrt(2._kind_phys)*sigmaw) )

    upqt(ksfcm)  = cwqt * upw(ksfcm) * sigmaqt/sigmaw
    upthv(ksfcm) = cwthv * upw(ksfcm) * sigmathv/sigmaw

    upqt(ksfcm) = qt(ksfcm)+upqt(ksfcm)
    upthv(ksfcm) = thv(ksfcm)+upthv(ksfcm)
    upthl(ksfcm) = upthv(ksfcm) / (1._kind_phys+zvir*upqt(ksfcm))
    upth(ksfcm)  = upthl(ksfcm)

    ! get cloud, lowest momentum level
    if (do_condensation) then
      call condensation_mf(upqt(ksfcm), upthl(ksfcm), p_zm(ksfcm), iexner_zm(ksfcm), &
                           thvn, qcn, thn, qln, qin, qsn, lmixn)
      upthv(ksfcm) = thvn
      upqc(ksfcm)  = qcn
      upql(ksfcm)  = qln
      upqi(ksfcm)  = qin
      upqs(ksfcm)  = qsn
      upth(ksfcm)  = thn
      if (qcn > 0._kind_phys) zcb = zm(ksfcm)
    else
      ! assume no cldliq
      upqc(ksfcm)  = 0._kind_phys
    end if

    do k = ksfcm, ktopm-kdir, kdir
      kt = k - (1-kdir)/2
      kn = k + kdir

      ! get microphysics, autoconversion
      if (do_precip .and. upqc(k) > 0._kind_phys) then
        call precip_mf(upqs(k),upqt(k),upw(k),dzt(kt),abs(zm(kn)-zcb),supqt(kt))
        supthl(kt) = -1._kind_phys*lmixn*supqt(kt)*iexner_zt(kt)/cpair
      else
        supqt(kt)  = 0._kind_phys
        supthl(kt) = 0._kind_phys
      end if
      ! integrate updraft
      if (do_tptke) then
        pturb = (1._kind_phys + clubb_mf_alphturb*sqrt(tke(k))/upw(k))
      else
        pturb = 1._kind_phys
      end if
      pentexp  = exp(-pent*pturb*dzt(kt))
      qtn  = qt(kt) *(1._kind_phys-pentexp ) + upqt (k)*pentexp + supqt(kt)
      thln = thl(kt)*(1._kind_phys-pentexp ) + upthl(k)*pentexp + supthl(kt)

      ! convert source terms to a tendency
      supqt(kt) = supqt(kt)*upw(k)/dzt(kt)
      supthl(kt) = supthl(kt)*upw(k)/dzt(kt)

      ! get cloud, momentum levels
      if (do_condensation) then
        call condensation_mf(qtn, thln, p_zm(kn), iexner_zm(kn), &
                             thvn, qcn, thn, qln, qin, qsn, lmixn)
        if (zcb.eq.zcb_unset .and. qcn > 0._kind_phys) zcb = zm(kn)
      else
        thvn = thln*(1._kind_phys+zvir*qtn)
      end if
      ! get buoyancy
      B=gravit*(0.5_kind_phys*(thvn + upthv(k))/thv(kt)-1._kind_phys)

      ! get wn^2
      wp = wb*pent*pturb
      if (wp==0._kind_phys) then
         wn2 = upw(k)**2._kind_phys+2._kind_phys*wa*B*dzt(kt)
      else
         pentw = exp(-2._kind_phys*wp*dzt(kt))
         wn2 = pentw*upw(k)**2._kind_phys+(1._kind_phys-pentw)*wa*B/wp
      end if

      if (wn2>0._kind_phys) then
        upw(kn)   = sqrt(wn2)
        upthv(kn) = thvn
        upthl(kn) = thln
        upqt(kn)  = qtn
        upqc(kn)  = qcn
        upqs(kn)  = qsn
        upa(kn)   = upa(k)
        upql(kn)  = qln
        upqi(kn)  = qin
        upth(kn)  = thn
      else
        plumeheight = zm(k)
        exit
      end if
    enddo

  end subroutine oneplume


subroutine buoyan_dilute( nzt, nzm, nup, dmpdz, q, t, p, z, pf, &
                          tp, qstp, tl, cape, cin, pblt, lcl, lel, lon, mx, &
                          msg, tpert, landfrac )
!-----------------------------------------------------------------------
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for
!            testing (dmpdp).
! 08/04/05 - Swap to convert dmpdz to dmpdp
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
!
! References:
! Raymond and Blythe (1992) JAS
!
! Author:
! Richard Neale - September 2004
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
! input arguments
!
   integer, intent(in) :: nzt, nzm            ! vertical grid sizes
   integer, intent(in) :: nup           ! number of plumes

   real(kind_phys), intent(in) :: dmpdz(nzt,nup)! Parcel fractional mass entrainment rate (/m) 3D

   real(kind_phys), intent(in) :: q(nzt)        ! spec. humidity
   real(kind_phys), intent(in) :: t(nzt)        ! temperature
   real(kind_phys), intent(in) :: p(nzt)        ! pressure
   real(kind_phys), intent(in) :: z(nzt)        ! height
   real(kind_phys), intent(in) :: pf(nzm)       ! pressure at interfaces
   integer,  intent(in) :: pblt         ! index of pbl depth
   integer,  intent(in) :: msg
   real(kind_phys), intent(in) :: tpert        ! perturbation temperature by pbl processes
   real(kind_phys), intent(in) :: landfrac

! output arguments
   real(kind_phys), intent(out) :: tp(nzt,nup)       ! parcel temperature
   real(kind_phys), intent(out) :: qstp(nzt,nup)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(kind_phys), intent(out) :: tl(nup)          ! parcel temperature at lcl
   real(kind_phys), intent(out) :: cape(nup)        ! convective aval. pot. energy.
   real(kind_phys), intent(out) :: cin (nup)        ! CIN
   integer,  intent(out) :: lcl(nup)         !
   integer,  intent(out) :: lel(nup)         !
   integer,  intent(out) :: lon              ! level of onset of deep convection
   integer,  intent(out) :: mx               ! level of max moist static energy

   !--------------------------Local Variables------------------------------

   ! =============================================================================== !
   ! GRID ORIENTATION GENERALIZATION VARIABLES
   ! ------------------------------------------------------------------------------- !
   ! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
   ! duplicating code, these variables abstract the vertical loop bounds and slices.
   !
   ! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
   ! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
   ! kdir          : Directional step (+1 for moving up, -1 for moving down).
   !
   ! STAGGERED GRID INDEXING
   ! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
   ! index of the cell center (zt) to the interface (zm) flips depending on whether
   ! memory is loaded top-down or bottom-up. These variables dynamically map them:
   !
   ! kt    : The active thermodynamic cell center associated with the current step.
   !         Upward Sweep:   kt = k - (1-kdir)/2
   !         Downward Sweep: kt = k - (1+kdir)/2
   !
   ! kn    : The NEXT momentum interface in the direction of the current sweep.
   !         Upward Sweep:   kn = k + kdir
   !         Downward Sweep: kn = k - kdir
   !
   ! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
   !         kt_up = k - (1-kdir)/2
   !
   ! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
   !         kt_dn = k - (1+kdir)/2
   ! =============================================================================== !
   integer :: ksfct, ktopt, kdir, kn
   integer lelten(nup,mf_num_cin)
   real(kind_phys) capeten(nup,mf_num_cin)     ! provisional value of cape
   real(kind_phys) cinten(nup,mf_num_cin)      ! provisional value of CIN
   real(kind_phys) tv(nzt)
   real(kind_phys) tpv(nzt,nup)
   real(kind_phys) buoy(nzt,nup)
   real(kind_phys) pl(nup)

   real(kind_phys) a1, a2, estp, plexp, hmax, hmn, y
   logical plge600(nup)
   integer knt(nup)
   real(kind_phys) e
   integer i, k, n

   real(kind_phys), parameter :: tiedke_add = 0.5_kind_phys
!-----------------------------------------------------------------------
   if (z(1) < z(nzt)) then
      ksfct = 1
      ktopt = nzt
      kdir = 1
   else
      ksfct = nzt
      ktopt = 1
      kdir = -1
   end if

   do n = 1,mf_num_cin
      do i = 1,nup
         lelten(i,n)  = ksfct
         capeten(i,n) = 0._kind_phys
         cinten (i,n) = 0._kind_phys
      end do
   end do

   lon = ksfct
   mx   = lon
   hmax = 0._kind_phys

   do i = 1,nup
      knt(i) = 0
      lel(i) = ksfct
      cape(i) = 0._kind_phys
      tp(:,i) = t(:)
      qstp(:,i) = q(:)
   end do

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   if (tht_tweaks) then
    tv  (:) = t(:) *(1._kind_phys+q(:)/epsilo)/ (1._kind_phys+q(:))
   else
    tv  (:) = t(:) *(1._kind_phys+1.608_kind_phys*q(:))/ (1._kind_phys+q(:))
   endif
!-tht
   do i = 1,nup
     tpv (:,i) = tv(:)
   end do
   buoy(:,:) = 0._kind_phys

! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
   do k = ksfct, msg-kdir, kdir
       hmn =(cpair+q(k)*cpliq)*t(k)/(1._kind_phys+q(k)) + (1._kind_phys+q(k)/epsilo)/(1._kind_phys+q(k))*gravit*z(k) &
              +(latvap-(cpliq-cpwv)*(t(k)-tmelt))*q(k)
       if ((k - pblt)*kdir <= 0 .and. (k - lon)*kdir >= 0 .and. hmn > hmax) then
          hmax = hmn
          mx = k
       end if
   end do

! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

   do i = 1,nup ! Initialise LCL variables.
      lcl(i) = mx
      tl(i) = t(mx)
      pl(i) = p(mx)
   end do

!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
   call parcel_dilute(nzt, nzm, nup, msg, mx, p, z, t, q, &
                      tpert, tp, tpv, qstp, pl, tl, lcl, &
                      landfrac, dmpdz)

! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,nup
      plge600(i) = pl(i).ge.600._kind_phys ! Just change to always allow buoy calculation.
   end do

! Main buoyancy calculation.
   do k = ksfct, msg-kdir, kdir
      do i=1,nup
         if ((k - mx)*kdir >= 0 .and. plge600(i)) then
          if (tht_tweaks) then
            tv(k) = t(k)* (1._kind_phys+q(k)/epsilo)/ (1._kind_phys+q(k))     !+tht
          else
            tv(k) = t(k)* (1._kind_phys+1.608_kind_phys*q(k))/ (1._kind_phys+q(k)) !orig
          endif
! +0.5K or not? (arbitrary at this point - introduce in parcel_dilute instead? tht)
            buoy(k,i) = tpv(k,i) - tv(k) + tiedke_add  ! +0.5K or not?
         else
            qstp(k,i) = q(k)
            tp(k,i)   = t(k)
            tpv(k,i)  = tv(k)
         endif
      end do
   end do

!-------------------------------------------------------------------------------
! beginning from one below top (first level p>40hPa, msg) check for at most
! num_cin levels of neutral buoyancy (LELten) and compute CAPEten between LCL
   do k = msg-2*kdir, ksfct, -kdir
      do i = 1,nup
         if ((k - lcl(i))*kdir > 0 .and. plge600(i)) then
            if (buoy(k-kdir,i) > 0._kind_phys .and. buoy(k,i) <= 0._kind_phys) then
               knt(i) = min(mf_num_cin,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do

! calculate convective available potential energy (cape).
   do n = 1,mf_num_cin
      do k = msg-kdir, ksfct, -kdir
         do i = 1,nup
            if (plge600(i) .and. (k - mx)*kdir >= 0 .and. (k - lelten(i,n))*kdir < 0) then
               ! Using pf(k) and pf(k+kdir). Since pf is on momentum grid (nzm=nzt+1),
               ! if k ranges from 1 to nzt, k+1 ranges from 2 to nzm (valid).
               ! For top-down (kdir=-1), k-1 ranges from nzt-1 to 0.
               ! The interface pressure above thermo layer k is k - (1-kdir)/2
               kn = k - (1-kdir)/2
               capeten(i,n) = capeten(i,n) + rair*buoy(k,i)*abs(log(pf(kn)/pf(kn+kdir)))
               cinten (i,n) = cinten (i,n) - rair*min(buoy(k,i),0._kind_phys)*abs(log(pf(kn)/pf(kn+kdir)))
            end if
         end do
      end do
   end do

! find maximum cape from all possible tentative capes from one sounding
   do n = 1,mf_num_cin
      do i = 1,nup
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            cin (i) = cinten (i,n) !+tht CIN
            lel(i) = lelten(i,n)
         end if
      end do
   end do

! put lower bound on cape for diagnostic purposes.
   do i = 1,nup
      cape(i) = max(cape(i), 0._kind_phys)
   end do

return
end subroutine buoyan_dilute


 subroutine parcel_dilute (nzt, nzm, nup, msg, klaunch, p, z, t, q, &
                           tpert, tp, tpv, qstp, pl, tl, lcl, &
                           landfrac, dmpdz)
!-tht

! Routine  to determine
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: nzt, nzm
integer, intent(in) :: nup
integer, intent(in) :: msg
integer, intent(in) :: klaunch

real(kind_phys), intent(in)                 :: tpert ! PBL temperature perturbation.
real(kind_phys), intent(in)                 :: landfrac
real(kind_phys), intent(in), dimension(nzt) :: p
real(kind_phys), intent(in), dimension(nzt) :: z
real(kind_phys), intent(in), dimension(nzt) :: t
real(kind_phys), intent(in), dimension(nzt) :: q

real(kind_phys), intent(inout), dimension(nzt,nup) :: tp    ! Parcel temp.
real(kind_phys), intent(inout), dimension(nzt,nup) :: qstp  ! Parcel water vapour (sat value above lcl).
real(kind_phys), intent(inout), dimension(nup)     :: tl    ! Actual temp of LCL.
real(kind_phys), intent(inout), dimension(nup)     :: pl    ! Actual pressure of LCL.
integer,  intent(inout), dimension(nup)     :: lcl   ! Lifting condesation level (first model level with saturation).

real(kind_phys), intent(out), dimension(nzt,nup)   :: tpv   ! Define tpv within this routine.

real(kind_phys), dimension(nzt,nup) :: dmpdz ! Parcel fractional mass entrainment rate (/m) 3D

! =============================================================================== !
! GRID ORIENTATION GENERALIZATION VARIABLES
! ------------------------------------------------------------------------------- !
! To support both top-down (CAM) and bottom-up (CLUBB) grid orientations without
! duplicating code, these variables abstract the vertical loop bounds and slices.
!
! ksfcm / ksfct : Index of the surface for momentum (m) and thermodynamic (t) grids.
! ktopm / ktopt : Index of the model top for momentum (m) and thermodynamic (t) grids.
! kdir          : Directional step (+1 for moving up, -1 for moving down).
!
! STAGGERED GRID INDEXING
! Because momentum (zm) and thermodynamic (zt) grids are staggered, the relative
! index of the cell center (zt) to the interface (zm) flips depending on whether
! memory is loaded top-down or bottom-up. These variables dynamically map them:
!
! kt    : The active thermodynamic cell center associated with the current step.
!         Upward Sweep:   kt = k - (1-kdir)/2
!         Downward Sweep: kt = k - (1+kdir)/2
!
! kn    : The NEXT momentum interface in the direction of the current sweep.
!         Upward Sweep:   kn = k + kdir
!         Downward Sweep: kn = k - kdir
!
! kt_up : The thermodynamic cell center physically ABOVE momentum interface k.
!         kt_up = k - (1-kdir)/2
!
! kt_dn : The thermodynamic cell center physically BELOW momentum interface k.
!         kt_dn = k - (1+kdir)/2
! =============================================================================== !
integer :: ksfct, ktopt, kdir

real(kind_phys) tmix(nzt,nup)        ! Tempertaure of the entraining parcel.
real(kind_phys) qtmix(nzt,nup)       ! Total water of the entraining parcel.
real(kind_phys) qsmix(nzt,nup)       ! Saturated mixing ratio at the tmix.
real(kind_phys) smix(nzt,nup)        ! Entropy of the entraining parcel.
real(kind_phys) xsh2o(nzt,nup)       ! Precipitate lost from parcel.
real(kind_phys) ds_xsh2o(nzt,nup)    ! Entropy change due to loss of condensate.
real(kind_phys) ds_freeze(nzt,nup)   ! Entropy change sue to freezing of precip.
real(kind_phys) dmpdz2d(nzt,nup)     ! variable detrainment rate

real(kind_phys) zl(nup) ! lcl

real(kind_phys) mp(nup)    ! Parcel mass flux.
real(kind_phys) qtp(nup)   ! Parcel total water.
real(kind_phys) sp(nup)    ! Parcel entropy.
real(kind_phys) sp0(nup)   ! Parcel launch entropy.
real(kind_phys) qtp0(nup)  ! Parcel launch total water.
real(kind_phys) mp0(nup)   ! Parcel launch relative mass flux.

real(kind_phys) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(kind_phys) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
real(kind_phys) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(kind_phys) senv       ! Environmental entropy at each grid point.
real(kind_phys) qtenv      ! Environmental total water "   "   ".
real(kind_phys) penv       ! Environmental total pressure "   "   ".
real(kind_phys) zenv
real(kind_phys) tenv       ! Environmental total temperature "   "   ".
real(kind_phys) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(kind_phys) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(kind_phys) dp         ! Layer thickness (center to center)
real(kind_phys) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(kind_phys) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(kind_phys) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(kind_phys) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(kind_phys) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.
real(kind_phys) dmpdz_lnd, dmpdz_mask

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat   ! Number of iterations for condensation/freezing loop.
integer i,k,ii      ! Loop counters.

real(kind_phys) est

if (z(1) < z(nzt)) then
   ksfct = 1
   ktopt = nzt
   kdir = 1
else
   ksfct = nzt
   ktopt = 1
   kdir = -1
end if

nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.

lwmax    = 1.e10_kind_phys   ! tht: don't precipitate
tscool   =-10._kind_phys     ! tht: allow even just mild supercooling?!

qtmix=0._kind_phys
smix=0._kind_phys

qtenv = 0._kind_phys
senv = 0._kind_phys
tenv = 0._kind_phys
penv = 0._kind_phys
zenv = 0._kind_phys

qtp0 = 0._kind_phys
sp0  = 0._kind_phys
mp0 = 0._kind_phys

qtp = 0._kind_phys
sp = 0._kind_phys
mp = 0._kind_phys

new_q = 0._kind_phys
new_s = 0._kind_phys

zl(:)=0._kind_phys

! **** Begin loops ****

do k = ksfct, msg-kdir, kdir
   do i=1,nup

! Initialize parcel values at launch level.

      if (k == klaunch) then
         qtp0(i) = q(k)   ! Parcel launch total water (assuming subsaturated) - OK????.

!+tht: formulate dilution on enthalpy not on entropy
         if (tht_tweaks) then
          sp0(i)  = enthalpy(t(k),p(k),qtp0(i),z(k))  ! Parcel launch enthalpy.
         else
          sp0(i)  = entropy (t(k),p(k),qtp0(i))         ! Parcel launch entropy.
         endif
!-tht
         mp0(i)  = 1._kind_phys       ! Parcel launch relative mass (=1 for dmpdp=0 i.e. undilute).
         smix(k,i)  = sp0(i)
         qtmix(k,i) = qtp0(i)
!+tht: since the function to invert for T is *identical* with sp0(i)=entropy(t), unless there is
! a coding error (likely, given the mess) the result must be t(i,k) (verified 21/2/2014)
         if (tht_tweaks) then
          tmix(k,i) = t(k)
          call qsat_hPa(tmix(k,i),p(k), est, qsmix(k,i))
         else
          tfguess = t(k)
          rcall = 1
          call ientropy (rcall,smix(k,i),p(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
         endif
!-tht
      end if

      if ((k - klaunch)*kdir > 0) then

         dp = -abs(p(k)-p(k-kdir))
         qtenv = 0.5_kind_phys*(q(k)+q(k-kdir))
         tenv  = 0.5_kind_phys*(t(k)+t(k-kdir))
         penv  = 0.5_kind_phys*(p(k)+p(k-kdir))
         zenv  = 0.5_kind_phys*(z(k)+z(k-kdir))

         if (tht_tweaks) then
          senv  = enthalpy(tenv,penv,qtenv,zenv) ! Enthalpy of environment.
         else
          senv  = entropy (tenv,penv,qtenv)      ! Entropy  of environment.
         endif

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*gravit)/(rair*tenv) ! in mb/m since  p in mb.
         dzdp = 1._kind_phys/dpdz                  ! in m/mb
!+tht
! NB: land fudge makes no sense to me - make dmpdz_lnd=dmpdz (as per default code, hard-wired to 1e-3)
        !dmpdp = dmpdz*dzdp
        !dmpdp = dmpdz(i)*dzdp              ! /mb Fractional entrainment 2D
         dmpdp = dmpdz(k,i)*dzdp            ! /mb Fractional entrainment 3D
!-tht

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv
         qtp(i) = qtp(i) - dmpdp*dp*qtenv
         mp(i)  = mp(i)  - dmpdp*dp

! Entrain s and qt to next level.

         smix(k,i)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(k,i) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

         tfguess = tmix(k-kdir,i)
         rcall = 2

         if (tht_tweaks) then
          call ienthalpy(rcall,smix(k,i),p(k),z(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
         else
          call ientropy (rcall,smix(k,i),p(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
         endif

         if (qsmix(k,i) <= qtmix(k,i) .and. qsmix(k-kdir,i) > qtmix(k-kdir,i)) then
            lcl(i) = k
            qxsk   = qtmix(k,i) - qsmix(k,i)
            qxskp1 = qtmix(k-kdir,i) - qsmix(k-kdir,i)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(k-kdir) - qxskp1/dqxsdp
            zl(i)  = z(k-kdir) - qxskp1/dqxsdp *dzdp
            dsdp   = (smix(k,i)  - smix(k-kdir,i))/dp
            dqtdp  = (qtmix(k,i) - qtmix(k-kdir,i))/dp
            slcl   = smix(k-kdir,i)  +  dsdp* (pl(i)-p(k-kdir))
            qtlcl  = qtmix(k-kdir,i) +  dqtdp*(pl(i)-p(k-kdir))

            tfguess = tmix(k,i)
            rcall = 3

            if (tht_tweaks) then
               call ienthalpy(rcall,slcl,pl(i),zl(i),qtlcl,tl(i),qslcl,tfguess)
            else
               call ientropy (rcall,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)
            endif
         endif
!
      end if !  k < klaunch


   end do ! Levels loop
end do ! Columns loop


!   if ( masterproc ) then
!     do k = 1,msg-1
!         do i = 1,nup
!            write(iulog,*) "after, k, nup, dmpdz ", k, i, dmpdz(k,i)
!         end do
!     end do
!   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously
!! provides latent heating to the mixed parcel and so this has to be added back
!! to it. But does this also increase qsmix as well? Also freezing processes


xsh2o = 0._kind_phys
ds_xsh2o = 0._kind_phys
ds_freeze = 0._kind_phys

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy


do k = ksfct, msg-kdir, kdir
   do i=1,nup

! Initialize variables at k=klaunch

      if (k == klaunch) then

! Set parcel values at launch level assume no liquid water.

         tp(k,i)    = tmix(k,i)
         qstp(k,i)  = q(k)
         if (tht_tweaks) then
           tpv(k,i)   =  (tp(k,i) + tpert) * (1._kind_phys+qstp(k,i)/epsilo) / (1._kind_phys+qstp(k,i)) !+tht OK with mx ratio
         else
           tpv(k,i)   =  (tp(k,i) + tpert) * (1._kind_phys+1.608_kind_phys*qstp(k,i)) / (1._kind_phys+qstp(k,i))
         endif
      end if

      if ((k - klaunch)*kdir > 0) then

         if (tht_tweaks) then
           smix(k,i)=entropy(tmix(k,i),p(k),qtmix(k,i)) !+tht make sure to use entropy here
         endif

!----
! Initiate loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.
! Iterate nit_lheat times for s,qt changes.
         do ii=0,nit_lheat-1

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
            xsh2o(k,i) = max (0._kind_phys, qtmix(k,i) - qsmix(k,i) - lwmax)
            ds_xsh2o(k,i) = ds_xsh2o(k-kdir,i) - cpliq * log (tmix(k,i)/tmelt) * max(0._kind_phys,(xsh2o(k,i)-xsh2o(k-kdir,i)))

            if (tmix(k,i) <= tmelt+tscool .and. ds_freeze(k-kdir,i) == 0._kind_phys) then
               ds_freeze(k,i) = (latice/tmix(k,i)) * max(0._kind_phys,qtmix(k,i)-qsmix(k,i)-xsh2o(k,i))
            end if

            if (tmix(k,i) <= tmelt+tscool .and. ds_freeze(k-kdir,i) /= 0._kind_phys) then
               ds_freeze(k,i) = ds_freeze(k-kdir,i)+(latice/tmix(k,i)) * max(0._kind_phys,(qsmix(k-kdir,i)-qsmix(k,i)))
            end if

! Adjust entropy and accordingly to sum of ds (be careful of signs).
            new_s = smix(k,i) + ds_xsh2o(k,i) + ds_freeze(k,i)

! Adjust liquid water and accordingly to xsh2o.
            new_q = qtmix(k,i) - xsh2o(k,i)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(k,i)
            rcall =4
            call ientropy (rcall,new_s, p(k), new_q, tmix(k,i), qsmix(k,i), tfguess)

         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water.

         tp(k,i)    = tmix(k,i)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)
         if (new_q > qsmix(k,i)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(k,i) = qsmix(k,i)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(k,i) = new_q
         end if

         if (tht_tweaks) then
           tpv(k,i) = (tp(k,i)+tpert)* (1._kind_phys+qstp(k,i)/epsilo) / (1._kind_phys+ new_q) !+tht
         else
           tpv(k,i) = (tp(k,i)+tpert)* (1._kind_phys+1.608_kind_phys*qstp(k,i)) / (1._kind_phys+ new_q)
         endif

      end if ! k > klaunch

   end do ! Loop for columns

end do  ! Loop for vertical levels.


return
end subroutine parcel_dilute

!-----------------------------------------------------------------------------------------
real(kind_phys) function entropy(TK,p,qtot)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(kind_phys), intent(in) :: p,qtot,TK
     real(kind_phys) :: qv,qst,e,est,L
     real(kind_phys), parameter :: pref = 1000._kind_phys

L = latvap - (cpliq - cpwv)*(TK-tmelt)         ! T IN CENTIGRADE

call qsat_hPa(TK, p, est, qst)

qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
e = qv*p / (epsilo +qv)

entropy = (cpair + qtot*cpliq)*log( TK/tmelt) - rair*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qst)

end FUNCTION entropy

!-----------------------------------------------------------------------------------------
SUBROUTINE ientropy (rcall,s,p,qt,T,qst,Tfg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg).
! Inverts entropy, pressure and total water qt
! for T and saturated vapor mixing ratio
!

  integer, intent(in) :: rcall
  real(kind_phys), intent(in)  :: s, p, Tfg, qt
  real(kind_phys), intent(out) :: qst, T
  real(kind_phys) :: est
  real(kind_phys) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
  integer :: i

  logical :: converged

  ! Max number of iteration loops.
  integer, parameter :: LOOPMAX = 100
  real(kind_phys), parameter :: EPS = 3.e-8_kind_phys

  converged = .false.

  ! Invert the entropy equation -- use Brent's method
  ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

  T = Tfg                  ! Better first guess based on Tprofile from conv.

  a = Tfg-10    !low bracket
  b = Tfg+10    !high bracket

  fa = entropy(a, p, qt) - s
  fb = entropy(b, p, qt) - s

  c=b
  fc=fb
  tol=0.001_kind_phys

  converge: do i=0, LOOPMAX
     if ((fb > 0.0_kind_phys .and. fc > 0.0_kind_phys) .or. &
          (fb < 0.0_kind_phys .and. fc < 0.0_kind_phys)) then
        c=a
        fc=fa
        d=b-a
        ebr=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if

     tol1=2.0_kind_phys*EPS*abs(b)+0.5_kind_phys*tol
     xm=0.5_kind_phys*(c-b)
     converged = (abs(xm) <= tol1 .or. fb == 0.0_kind_phys)
     if (converged) exit converge

     if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr=fb/fa
        if (a == c) then
           pbr=2.0_kind_phys*xm*sbr
           qbr=1.0_kind_phys-sbr
        else
           qbr=fa/fc
           rbr=fb/fc
           pbr=sbr*(2.0_kind_phys*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_kind_phys))
           qbr=(qbr-1.0_kind_phys)*(rbr-1.0_kind_phys)*(sbr-1.0_kind_phys)
        end if
        if (pbr > 0.0_kind_phys) qbr=-qbr
        pbr=abs(pbr)
        if (2.0_kind_phys*pbr  <  min(3.0_kind_phys*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
           ebr=d
           d=pbr/qbr
        else
           d=xm
           ebr=d
        end if
     else
        d=xm
        ebr=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

     fb = entropy(b, p, qt) - s

  end do converge

  T = b
  call qsat_hPa(T, p, est, qst)

  if (.not. converged) then
     call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

end SUBROUTINE ientropy

! Wrapper for qsat_water that does translation between Pa and hPa
! qsat_water uses Pa internally, so get it right, need to pass in Pa.
! Afterward, set es back to hPa.
subroutine qsat_hPa(t, p, es, qm)
  use wv_saturation, only: qsat_water

  ! Inputs
  real(kind_phys), intent(in) :: t    ! Temperature (K)
  real(kind_phys), intent(in) :: p    ! Pressure (hPa)
  ! Outputs
  real(kind_phys), intent(out) :: es  ! Saturation vapor pressure (hPa)
  real(kind_phys), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

  call qsat_water(t, p*100._kind_phys, es, qm)

  es = es*0.01_kind_phys

end subroutine qsat_hPa

!-----------------------------------------------------------------------------------------
real(kind_phys) function enthalpy(TK,p,qtot,z)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
!
     real(kind_phys), intent(in) :: p,qtot,TK,z
     real(kind_phys) :: qv,qst,e,est,L

L = latvap - (cpliq - cpwv)*(TK-tmelt)

call qsat_hPa(TK, p, est, qst)
qv = min(qtot,qst)                         ! Partition qtot into vapor part only.

 enthalpy = (cpair + qtot*cpliq)* TK         + L*qv + (1._kind_phys+qtot)*gravit*z

return
end FUNCTION enthalpy

!-----------------------------------------------------------------------------------------
 SUBROUTINE ienthalpy (rcall,s,p,z,qt,T,qst,Tfg) !identical with iENTROPY, only function calls swapped
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg).
! Inverts entropy, pressure and total water qt
! for T and saturated vapor mixing ratio
!

  integer, intent(in) :: rcall
  real(kind_phys), intent(in)  :: s, p, z, Tfg, qt
  real(kind_phys), intent(out) :: qst, T
  real(kind_phys) :: est
  real(kind_phys) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
  integer :: i

  logical :: converged

  ! Max number of iteration loops.
  integer, parameter :: LOOPMAX = 100
  real(kind_phys), parameter :: EPS = 3.e-8_kind_phys

  converged = .false.

  ! Invert the entropy equation -- use Brent's method
  ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

  T = Tfg                  ! Better first guess based on Tprofile from conv.

  a = Tfg-10    !low bracket
  b = Tfg+10    !high bracket

  fa = enthalpy(a, p, qt,z) - s
  fb = enthalpy(b, p, qt,z) - s

  c=b
  fc=fb
  tol=0.001_kind_phys

  converge: do i=0, LOOPMAX
     if ((fb > 0.0_kind_phys .and. fc > 0.0_kind_phys) .or. &
          (fb < 0.0_kind_phys .and. fc < 0.0_kind_phys)) then
        c=a
        fc=fa
        d=b-a
        ebr=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if

     tol1=2.0_kind_phys*EPS*abs(b)+0.5_kind_phys*tol
     xm=0.5_kind_phys*(c-b)
     converged = (abs(xm) <= tol1 .or. fb == 0.0_kind_phys)
     if (converged) exit converge

     if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr=fb/fa
        if (a == c) then
           pbr=2.0_kind_phys*xm*sbr
           qbr=1.0_kind_phys-sbr
        else
           qbr=fa/fc
           rbr=fb/fc
           pbr=sbr*(2.0_kind_phys*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_kind_phys))
           qbr=(qbr-1.0_kind_phys)*(rbr-1.0_kind_phys)*(sbr-1.0_kind_phys)
        end if
        if (pbr > 0.0_kind_phys) qbr=-qbr
        pbr=abs(pbr)
        if (2.0_kind_phys*pbr  <  min(3.0_kind_phys*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
           ebr=d
           d=pbr/qbr
        else
           d=xm
           ebr=d
        end if
     else
        d=xm
        ebr=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

     fb = enthalpy(b, p, qt,z) - s

  end do converge

  T = b
  call qsat_hPa(T, p, est, qst)

  if (.not. converged) then
     call endrun('**** ZM_CONV IENTHALPY: Tmix did not converge ****')
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

 end SUBROUTINE ienthalpy

end module clubb_mf
