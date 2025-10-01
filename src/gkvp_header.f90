MODULE GKV_header
!-------------------------------------------------------------------------------
!
!    Header for general use in the fluxtube code
!
!    Notes
!    -----
!      There are some restrictions on numerical parameters:
!         mod( global_nz, nprocz ) = 0, due to z parallelization.
!         mod( global_nv, nprocv ) = 0, due to v parallelization.
!         mod( global_nm+1, nprocm ) = 0, due to m parallelization.
!         nm>=3, due to fft and colli.
!         nzb>=2, due to 4th-oreder z derivative.
!      There are some recommendations on numerical parameters:
!         mod( nxw, nprocw ) = 0, due to w parallelization.
!         mod( global_ny+1, nprocw ) = 0, due to w parallelization.
!         nzb<=nz due to data copy in zfilter.
!
!    Update history
!    --------------
!      gkvp_f0.61 (S. Maeyama, Mar 2021)
!        - equib_type is extended from len=8 to len=15.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)

!--------------------------------------
!  Dimension size (grid numbers)
!--------------------------------------
!  Global simulation domain 
!  in  x, y,z,v,m (0:2*nxw-1,  0:2*nyw-1,-global_nz:global_nz-1,1:2*global_nv,0:global_nm)
!  in kx,ky,z,v,m (   -nx:nx,0:global_ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm)
!
!  integer, parameter :: nxw = 12, nyw = 12
!  integer, parameter :: nx = 7, global_ny = 7 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 2, nyw = 72
!  integer, parameter :: nx = 0, global_ny = 47 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 12, nyw = 12
!  integer, parameter :: nx = 7, global_ny = 7 ! 2/3 de-aliasing rule
  integer, parameter :: nxw = 36, nyw = 36
  integer, parameter :: nx = 23, global_ny = 23 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 72, nyw = 72
!  integer, parameter :: nx = 47, global_ny = 47 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 96, nyw = 96
!  integer, parameter :: nx = 63, global_ny = 64 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 96, nyw = 96
!  integer, parameter :: nx = 63, global_ny = 64 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 108, nyw = 108
!  integer, parameter :: nx = 71, global_ny = 71
!  integer, parameter :: nxw = 120, nyw = 120
!  integer, parameter :: nx = 79, global_ny = 79 ! 2/3 de-aliasing rule
!  integer, parameter :: nxw = 144, nyw = 144
!  integer, parameter :: nx = 95, global_ny = 95

!  integer, parameter :: global_nz = 24, global_nv = 24, global_nm = 15
!  integer, parameter :: global_nz = 48, global_nv = 24, global_nm = 15
!  integer, parameter :: global_nz = 48, global_nv = 48, global_nm = 15
!  integer, parameter :: global_nz = 96, global_nv = 24, global_nm = 15
!  integer, parameter :: global_nz = 96, global_nv = 48, global_nm = 31
!  integer, parameter :: global_nz = 192, global_nv = 24, global_nm = 15
!  integer, parameter :: global_nz = 48, global_nv = 24, global_nm = 31
!  integer, parameter :: global_nz = 48, global_nv = 48, global_nm = 31
!  integer, parameter :: global_nz = 96, global_nv = 48, global_nm = 31
!  integer, parameter :: global_nz = 96, global_nv = 48, global_nm = 31
!  integer, parameter :: global_nz = 96, global_nv = 64, global_nm = 31
  integer, parameter :: global_nz = 64, global_nv = 96, global_nm = 31
!  integer, parameter :: global_nz = 96, global_nv = 96, global_nm = 31
!  integer, parameter :: global_nz = 120, global_nv = 96, global_nm = 31
!   integer, parameter :: global_nz = 144, global_nv = 48, global_nm = 31
!   integer, parameter :: global_nz = 192, global_nv = 48, global_nm = 31
!   integer, parameter :: global_nz = 192, global_nv = 64, global_nm = 31
!   integer, parameter :: global_nz = 252, global_nv = 64, global_nm = 31
!  integer, parameter :: global_nz = 100, global_nv = 48, global_nm = 31
!  integer, parameter :: global_nz = 48, global_nv = 32, global_nm = 31
!  integer, parameter :: global_nz = 96, global_nv = 32, global_nm = 31
!  integer, parameter :: global_nz = 48, global_nv = 96, global_nm = 31

  integer, parameter :: nzb = 3, &  ! the number of ghost grids in z
                        nvb = 2     ! the number of ghost grids in v and m

!--------------------------------------
!  Data distribution for MPI
!--------------------------------------

!  integer, parameter :: nprocw = 1, nprocz = 8, nprocv = 1, nprocm = 2, nprocs = 2
!  integer, parameter :: nprocw = 2, nprocz = 8, nprocv = 1, nprocm = 1, nprocs = 2 ! 32
! integer, parameter :: nprocw = 4, nprocz = 2, nprocv = 2, nprocm = 1, nprocs = 2 ! 32
!  integer, parameter :: nprocw = 2, nprocz = 2, nprocv = 2, nprocm = 2, nprocs = 2 ! 32
!   integer, parameter :: nprocw = 1, nprocz = 4, nprocv = 1, nprocm = 16, nprocs = 2 ! 32
!  integer, parameter :: nprocw = 2, nprocz = 1, nprocv = 1, nprocm = 8, nprocs = 2
!!  integer, parameter :: nprocw = 8, nprocz = 1, nprocv = 1, nprocm = 8, nprocs = 2
!  integer, parameter :: nprocw = 4, nprocz = 2, nprocv = 4, nprocm = 2, nprocs = 2 ! 128
!  integer, parameter :: nprocw = 2, nprocz = 4, nprocv = 4, nprocm = 2, nprocs = 2 ! 128
  integer, parameter :: nprocw = 4, nprocz = 4, nprocv = 4, nprocm = 2, nprocs = 2 ! 256
!  integer, parameter :: nprocw = 1, nprocz = 16, nprocv = 4, nprocm = 2, nprocs = 2 ! 256
!  integer, parameter :: nprocw = 4, nprocz = 8, nprocv = 4, nprocm = 4, nprocs = 2 ! 1024
!  integer, parameter :: nprocw = 4, nprocz = 16, nprocv = 4, nprocm = 4, nprocs = 2 ! 2048
!  integer, parameter :: nprocw = 8, nprocz = 12, nprocv = 4, nprocm = 4, nprocs = 2 ! 3072
!  integer, parameter :: nprocw = 4, nprocz = 32, nprocv = 4, nprocm = 4, nprocs = 2 ! 4096
!  integer, parameter :: nprocw = 16, nprocz = 12, nprocv = 4, nprocm = 4, nprocs = 2 ! 6144
!  integer, parameter :: nprocw = 18, nprocz = 12, nprocv = 8, nprocm = 4, nprocs = 2 ! 13824
!  integer, parameter :: nprocw = 12, nprocz = 24, nprocv = 8, nprocm = 4, nprocs = 2 !
!   integer, parameter :: nprocw = 12, nprocz = 20, nprocv = 8, nprocm = 4, nprocs = 2 ! 15360
!  integer, parameter :: nprocw = 20, nprocz = 12, nprocv = 8, nprocm = 4, nprocs = 2 ! 15360
!  integer, parameter :: nprocw = 16, nprocz = 16, nprocv = 8, nprocm = 4, nprocs = 2 ! 16384

!--------------------------------------
!  Parameters for variable sizes
!--------------------------------------
!  Local simulation domain 
!  in kx,ky,z,v,m (divided in ky,z,v,m) (   -nx:nx,      0:ny,-nz:nz-1,1:2*nv,0:nm)
!  in  x,ky,z,v,m (divided in ky,z,v,m) (0:2*nxw-1,      0:ny,-nz:nz-1,1:2*nv,0:nm)
!  in  y, x,z,v,m (divided in  x,z,v,m) (    0:nyw,0:nxw_size,-nz:nz-1,1:2*nv,0:nm)

  integer, parameter :: nxw_size = (2*nxw-1)/nprocw     ! local allocation size (0:nxw_size)
  integer, parameter :: ny       = global_ny / nprocw   ! local allocation size (0:ny)

  integer, parameter :: nz = global_nz / nprocz,          &
                        nv = global_nv / nprocv,          &
                        nm = (global_nm + 1) / nprocm - 1,&
                        ns = nprocs

  integer, parameter :: nxyz = (2*nx+1)*(ny+1)*(2*nz), &
                        nxy  = (2*nx+1)*(ny+1)

  integer, parameter :: nnx = nxw*2, nny = nyw*2

!--------------------------------------
!  Constants
!--------------------------------------

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, &
                                 ep_0  = 8.8541878128d-12,     & ! vacuum permittivity
                                 mu_0  = 4._DP*pi*1.d-7,       & ! vacuum permeability [N/A^2]
                                 ue    = 1.60217663d-19,       & ! elementary charge [C]
                                 m_e   = 9.1093837d-31,        & ! electron mass [kg]
                                 m_p   = 1.67262192d-27,       & ! proton mass [kg]
                                 twopi = pi * 2._DP,         &
                                 eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )


!--------------------------------------
!  Index Range
!--------------------------------------

! ---- y dimension -------
  integer :: ist_y                           ! local start index of y
  integer :: iend_y                          ! local end   index of y
  integer :: nsize_y                         ! local size of y
  integer :: ist1_y                          ! local start index of y for global start index 1 

  integer :: ist_y_g                         ! global start index of y
  integer :: iend_y_g                        ! global end   index of y

! ---- xw dimension  -------
  integer :: ist_xw                           ! local start index of xw
  integer :: iend_xw                          ! local end   index of xw
  integer :: nsize_xw                         ! local size of xw

  integer :: ist_xw_g                         ! global start index of xw
  integer :: iend_xw_g                        ! global end   index of xw

!--------------------------------------
!  Parameters for time
!--------------------------------------

  real(kind=DP) :: e_limit                           ! elapsed time limit of a job
  real(kind=DP) :: tend                              ! end time
  real(kind=DP) :: dtout_fxv, dtout_ptn, dtout_eng   ! time-spacing for output
  real(kind=DP) :: dtout_dtc                         ! time-spacing for dt control


!--------------------------------------
!  Configuration parameters to be 
!    initialized in init subroutine
!--------------------------------------

!  real(kind=DP), &
!    dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)  :: kvd, kvs
  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm)  :: vdx, vdy, vsy
  real(kind=DP), &
    dimension(-nx:nx,0:ny,-nz:nz-1,0:nm)         :: j0, j1, j2
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: g0, ksq
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: fct_poisson, fct_e_energy
  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: fct_ampere, fct_m_energy
  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm) :: fmx
  real(kind=DP), dimension(-nx:nx)               :: fctgt

!!!%%% Parameters for colli_full %%%
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm) :: xxa
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1) :: nu_h, nu_g, nu_d, nu_p
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm)               :: nu_hs, nu_gs, nu_ds, nu_ps
!!!  real(kind=DP), dimension(1:6,-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1) :: x_tst, y_fld 
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,0:ns-1,1:2) :: c_t0
!!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm,0:ns-1,1:6) :: vfunc
!!!  real(kind=DP), dimension(-nz:nz-1,1:2*nv,0:nm,0:ns-1,1:6) :: vfunc
!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:6) :: jfunc
!!!  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: adbtc
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real(kind=DP), dimension(0:ns-1,0:ns-1) :: ctauiv, calpha, ctheta, cgamma, ceta, cxi

  real(kind=DP), dimension(-nx:nx)          :: kx
  real(kind=DP), dimension(0:ny)            :: ky
  real(kind=DP), dimension(-nz:nz-1)        :: zz, omg
  real(kind=DP), dimension(1:2*nv)          :: vl
  real(kind=DP), dimension(0:nm)            :: mu
  real(kind=DP), dimension(-nz:nz-1,0:nm)   :: vp, mir
  real(kind=DP), dimension(-nz:nz-1)        :: dvp
  real(kind=DP), dimension(-nz:nz-1)        :: dpara, rootg

  complex(kind=DP), dimension(0:ny)         :: ck
  integer, dimension(0:ny)                  :: dj

  real(kind=DP) :: dt_max, dt
  logical :: adapt_dt
                                                     !!! Parameters for the R0 units
  real(kind=DP), dimension(0:ns-1) ::   R0_Ln,  &    ! R0/Lns
                                        R0_Lt,  &    ! R0/Lts
                                           nu,  &    ! collision freq.   
                                         Anum,  &    ! mass number
                                         Znum,  &    ! charge number     
                                          fcs,  &    ! charge-density fraction 
                                          sgn,  &    ! signs of charge   
                                          tau,  &    ! T-ratio
                                         dns1        ! initial perturbation amp.
  real(kind=DP) :: dv, cfsrf, lambda_i, q_0, q_bar, beta, tau_ad, vmax
  real(kind=DP) :: mach, uprime, gamma_e, kxmin_g, kymin_g, tlim_exb
  real(kind=DP) :: Nref, Lref, Tref, Zeff
  integer       :: iFLR, icheck, ibprime, nx0
  real(kind=DP) :: baxfactor

  real(kind=DP) :: courant_num 

  real(kind=DP) :: dref ! fujita added

!--------------------------------------
!  Type of calculation
!--------------------------------------

  character(9)  :: calc_type, &  ! "linear", "lin_freq", "nonlinear"
                   z_bound,   &  ! "zerofixed", "outflow", "mixed"
                   z_filt,    &  ! "on", "off"
                   z_calc,    &  ! "cf4", "up5"
                   col_type,  &  ! "LB", "full", "lorentz"
                   time_advnc    ! "rkg4", "imp_colli", "auto_init"
  real(kind=DP) :: art_diff

  integer :: num_triad_diag 

! fujita added
  integer :: calc_disp           ! if >0, calculate disp. rel. of feedback instability

!--------------------------------------
!  Parameters for numerical settings
!--------------------------------------

  integer :: inum
  logical :: ch_res, init_random
  character(512) :: f_log, f_hst, f_phi, f_fxv, f_cnt
  character(15)  :: equib_type  ! "analytic", "s-alpha", "s-alpha-shift",
                                ! "circ-MHD", "vmec", "eqdsk", "slab"

! --- unit numbers for I/O
  integer, parameter :: inml = 5,  & 
                        olog = 10, &
                        icnt = 20, &
                        ophi = 30, &
                        oAl  = 31, &
                        omom = 32, &
                        otrn = 33, &
                        otri = 34, &
                        ofxv = 40, &
                        ocnt = 50, &
                        odtc = 59, &
                        oeng = 60, &
                        omen = 61, &
                        owes = 62, &
                        owem = 63, &
                        oges = 64, &
                        ogem = 65, &
                        oqes = 66, &
                        oqem = 67, &
                        obln = 68, &
                        ofrq = 69, &
                        odsp = 70, &
                        ocst = 71, &
                        inbz = 14, &
                        ivmc = 15, &
                        omtr = 16, &
                        odbg = 73, &
                        ofzv = 74, &
                        ofbk = 75, &
                        odAt = 76, &
                        ovmc = olog


END MODULE GKV_header
