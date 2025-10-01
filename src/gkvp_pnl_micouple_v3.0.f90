MODULE GKV_pnl

! calculare parallel nonlinear term (PNL)
  
  use GKV_header
  use GKV_mpienv
  use GKV_intgrl, only: intgrl_v0_moment_ms, intgrl_thet, intgrl_v0_moment
  use GKV_fft, only: &
           plan_x_forward, plan_x_backward, & ! for fugaku
           plan_y_forward, plan_y_backward, & ! for fugaku
           plan_xf_y2x, plan_xb_y2x, &
           plan_yf_y2x, plan_yb_y2x, &
           planr_xf_y2x, planr_xb_y2x, &
           planr_yf_y2x, planr_yb_y2x
  use GKV_clock, only: clock_sta, clock_end
  use GKV_micouple, only: micouple_bndry_bound_e
  use GKV_fld, only: fld_emfield_hh
  use GKV_bndry, only: bndry_bound_e

  implicit none

  private
  complex(kind=DP), save, dimension(0:ny,0:2*nxw-1) :: uikx_y2x, uiky_y2x
  real(kind=DP), save :: pnl_max_eachrank

! for fugaku
  integer, parameter :: nbuff = ((2*nz)*(nm+1)-1)/nprocw + 1
  integer, save :: nchunk_zm  = 1, nchunk_yb  = 1, nchunk_xb = 1
  integer, save :: nchunk_yzm = 1, nchunk_xzm = 1
  
  public pnl_term, pnl_max_eachrank, pnl_sum

CONTAINS
!--------------------------------------
  SUBROUTINE pnl_term( ff, psi, chi, dh, cf, ef )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef

    if (trim(calc_type) == "nonlinear") then
       call CALC_PNL_term_y2x( ff, psi, chi, dh, cf, ef )
       !!call CALC_PNL_term_y2zm( ff, psi, dh, cf, ef ) ! fugaku
    else
       ef(:,:,:,:,:) = ( 0._DP, 0._DP )
    end if

    return

  END SUBROUTINE pnl_term

!--------------------------------------
  SUBROUTINE CALC_PNL_term_y2x( ff, psi, chi, dh, cf, ef )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: dfdv
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: pnl_e, pnl_m

    real(kind=DP), allocatable, dimension(:,:,:,:) :: dpdz, dAdt
    real(kind=DP) :: cefv
    integer :: mx, my, iz, iv, im
    
    allocate(dfdv(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dpdz(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
    allocate(dAdt(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))    
    allocate(pnl_e(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(pnl_m(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
!
    call bndry_bound_e( psi )

! M-I coupling part
    if ( trim(z_bound) == "mi_couple"  ) then
      do im = 0, nm
        call micouple_bndry_bound_e ( psi(:,:,:,im) )
      end do
    end if
!M-I coupling part

!...zero reset
    dpdz(:,:,:,:)    = 0._DP
    dAdt(:,:,:,:)    = 0._DP
    dfdv(:,:,:,:,:)  = ( 0._DP, 0._DP ) 
    pnl_e(:,:,:,:,:) = ( 0._DP, 0._DP )   
    pnl_m(:,:,:,:,:) = ( 0._DP, 0._DP )

!...d/dv element
    cefv  = sqrt( Anum(ranks) / tau(ranks) ) / ( 12._DP * dv )

!...make df/dv
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do iv = 1, 2*nv
              dfdv(mx,my,iz,iv,im) = cefv * ( - ff(mx,my,iz,iv+2,im)      &
                                      + 8._DP * ff(mx,my,iz,iv+1,im)      &
                                      - 8._DP * ff(mx,my,iz,iv-1,im)      &
                                              + ff(mx,my,iz,iv-2,im) )
            end do
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
       dfdv(0,0,:,:,:) = ( 0._DP, 0._DP )
    end if
 
    call PNL_es_term_y2x( psi, chi, dfdv, dpdz, pnl_e ) ! electrostatic part

    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              ef(mx,my,iz,iv,im) = ef(mx,my,iz,iv,im) + dref * pnl_e(mx,my,iz,iv,im)               
            end do
          end do
        end do
      end do
    end do
    
    if( beta>0._DP ) then
      call PNL_em_term_y2x( ff, dfdv, dh, cf, ef, dAdt, pnl_m ) ! magnetic part
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ef(mx,my,iz,iv,im) = ef(mx,my,iz,iv,im) + dref * pnl_m(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do 
      end do
    end if

!...debug
    !iz=nz
    !iv=2*nv
    !im=0   
    !do my=0,ny
    !   do mx=0,nx
    !     write(olog, '(2i6, 99ES16.6)' ) mx, my, pnl_e(mx,my,iz,iv,im), pnl_m(mx,my,iz,iv,im)
    !  end do
    !end do
    !write(olog, *)
    !flush(olog)

    !...estimate time step size 
    call pnl_estimate_maxvel_y2x(dpdz,dAdt)

    deallocate( dfdv )
    deallocate( dpdz )
    deallocate( dAdt )
    deallocate( pnl_e )
    deallocate( pnl_m )
    
  END SUBROUTINE CALC_PNL_term_y2x

!--------------------------------------
  SUBROUTINE PNL_es_term_y2x( psi, chi, dfdv, dpdz, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl
    real(kind=DP), dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1), intent(out) :: dpdz
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdz
    real(kind=DP), dimension(-nz:nz-1) :: cefz
    integer :: mx, my, iz, iv, im
    
    allocate(wdz(-nx:nx,0:ny,-nz:nz-1,0:nm))   

    do iz = -nz, nz-1
      cefz(iz)  = 1._DP / ( 12._DP * dpara(iz) )
    end do

    wdz(:,:,:,:) = ( 0._DP, 0._DP )

!...make nonlinear part
    if( beta>0._DP ) then
      call calc_Ez_pssn_y2x( psi, chi, wdz )
    end if

!...add linear part
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx

            wdz(mx,my,iz,im) = cefz(iz) * ( - psi(mx,my,iz+2,im)    &
                                    + 8._DP * psi(mx,my,iz+1,im)    &
                                    - 8._DP * psi(mx,my,iz-1,im)    &
                                            + psi(mx,my,iz-2,im) )  &
                                    + wdz(mx,my,iz,im)
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
       wdz(0,0,:,:) = ( 0._DP, 0._DP )
    end if
    
    call pnl_dft_y2x( dfdv, wdz, dpdz, pnl )

    deallocate( wdz )


  END SUBROUTINE PNL_es_term_y2x


!--------------------------------------
  SUBROUTINE PNL_em_term_y2x( ff, dfdv, dh, cf, ef, dAdt, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv     
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf, ef
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl   
    real(kind=DP), dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1), intent(out) :: dAdt
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: df_st!!, cefX
    !!complex(kind=DP), allocatable, dimension(:,:,:) :: X1, X2
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdt
    complex(kind=DP), allocatable, dimension(:,:,:) :: dpsi
    integer :: mx, my, iz, iv, im

    allocate(df_st(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dpsi(-nx:nx,0:ny,-nz:nz-1))
    allocate(wdt(-nx:nx,0:ny,-nz:nz-1,0:nm))

    df_st(:,:,:,:,:) = ( 0._DP, 0._DP )
    dpsi(:,:,:)      = ( 0._DP, 0._DP )

    ! df w/o dA/dt term
    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              df_st(mx,my,iz,iv,im) = dh(mx,my,iz,iv,im) + cf(mx,my,iz,iv,im) - ef(mx,my,iz,iv,im)
            end do
          end do
        end do
      end do
    end do

    call fld_emfield_hh ( df_st, dpsi )

    ! J_0 * dA/dt
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            wdt(mx,my,iz,im) = j0(mx,my,iz,im) * dpsi(mx,my,iz)
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
       wdt(0,0,:,:) = ( 0._DP, 0._DP )
    end if

    call pnl_dft_y2x( dfdv, wdt, dAdt, pnl )

    deallocate( wdt )
    deallocate( df_st )
    deallocate( dpsi )

!!    deallocate( df_st )
!!    deallocate( X1 )
!!    deallocate( X2 )
!!    deallocate( cefX )

!...old method

    ! df w/o dA/dt term
!!    do im = 0, nm
!!      do iv = 1, 2*nv
!!        do iz = -nz, nz-1
!!          do my = ist_y, iend_y
!!            do mx = -nx, nx
!!              df_st(mx,my,iz,iv,im) = j0(mx,my,iz,im) * sgn(ranks) * Znum(ranks) *                 &
!!                                      sqrt( tau(ranks) / Anum(ranks) ) * vl(iv) *                  &
!!                                      (                                                            &
!!                                      dh(mx,my,iz,iv,im) + cf(mx,my,iz,iv,im) - ef(mx,my,iz,iv,im) &
!!                                      )
!!--!              cefX(mx,my,iz,iv,im) =  ( Znum(ranks) **2 / Anum(ranks) ) *         &
!!--!                                      j0(mx,my,iz,im) ** 2 * ( fmx(iz,iv,im) + dref * ff(mx,my,iz,iv,im) )
!... zeroth approx
!              cefX(mx,my,iz,iv,im) =  ( Znum(ranks) **2 / Anum(ranks) ) *         &
!                                      j0(mx,my,iz,im) ** 2 * fmx(iz,iv,im)
                                      !!(                                                         &
                                      !!vl(iv) ** 2 * fmx(iz,iv,im)  + dref * ff(mx,my,iz,iv,im)  &
                                      !!)
!!            end do
!!          end do
!!        end do
!!      end do
!!    end do
!!    call intgrl_v0_moment_ms( df_st, X1 )
!!    call intgrl_v0_moment_ms( cefX,  X2 )


    
  END SUBROUTINE PNL_em_term_y2x

!--------------------------------------
  SUBROUTINE calc_Ez_pssn_y2x( psi, chi, pssn )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: pssn
    real(kind=DP), dimension(:,:,:,:), allocatable :: dpdx, dpdy, dadx, dady
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: &
                                   wdx1o, wdy1o, wdx2o, wdy2o, &
                                   wdx1e, wdy1e, wdx2e, wdy2e, &
                                   pb1, pb2
    integer :: iv, iprocw
    integer       :: mx, my
    integer, save :: iflg = 0
    
      allocate(dpdx(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dpdy(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dadx(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dady(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(wdx1o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy1o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx1e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy1e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))

      allocate(pb1(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(pb2(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))

      if( iflg == 0 ) then
        iflg = 1
        do mx = 0, nx
          do my = ist_y, iend_y
            uikx_y2x(my,mx) = kx(mx) * ui
            uiky_y2x(my,mx) = ky(my) * ui
          end do
        end do
        do mx = -nx, -1
          do my = ist_y, iend_y
            uikx_y2x(my,2*nxw+mx) = kx(mx) * ui
            uiky_y2x(my,2*nxw+mx) = ky(my) * ui
          end do
        end do
      end if

!$OMP parallel default(none)                      &
!$OMP shared(psi,chi,pssn,dpdx,dpdy,dadx,dady)   &
!$OMP shared(wdx1o,wdy1o,wdx2o,wdy2o) &
!$OMP shared(wdx1e,wdy1e,wdx2e,wdy2e) &
!$OMP shared(pb1,pb2) &
!$OMP private(iprocw)

!$OMP workshare
      pssn(:,:,:,:) = (0._DP, 0._DP)
      wdx1o(:,:,:,:,:) = (0._DP, 0._DP)
      wdy1o(:,:,:,:,:) = (0._DP, 0._DP)
      wdx1e(:,:,:,:,:) = (0._DP, 0._DP)
      wdy1e(:,:,:,:,:) = (0._DP, 0._DP)
      dpdx(:,:,:,:) = 0._DP
      dpdy(:,:,:,:) = 0._DP
      dadx(:,:,:,:) = 0._DP
      dady(:,:,:,:) = 0._DP
!$OMP end workshare


      call pnl_pack_psi_y2x(psi,wdx1o,wdy1o)
!$OMP barrier
      call pnl_transpose_y2x(wdx1o,wdx2o)
      call pnl_transpose_y2x(wdy1o,wdy2o)
      call pnl_pack_psi_y2x(chi,wdx1e,wdy1e)
!$OMP barrier
      call pnl_transpose_y2x(wdx1e,wdx2e)
      call pnl_transpose_y2x(wdy1e,wdy2e)
      call pnl_backwardfft_y2x(wdx2o,dpdx)
      call pnl_backwardfft_y2x(wdy2o,dpdy)
!$OMP barrier
      call pnl_backwardfft_y2x(wdx2e,dadx)
      call pnl_backwardfft_y2x(wdy2e,dady)

      call pnl_pssn_realspcal_y2x( dpdx, dpdy, dadx, dady, pb1 ) 
      call pnl_transpose_x2y( pb1, pb2 )
      call pnl_pssn_unpack_x2y( pb2, pssn ) 
!$OMP end parallel
      
      deallocate(dpdx)
      deallocate(dpdy)
      deallocate(dadx)
      deallocate(dady)
      deallocate(wdx1o)
      deallocate(wdy1o)
      deallocate(wdx2o)
      deallocate(wdy2o)
      deallocate(wdx1e)
      deallocate(wdy1e)
      deallocate(wdx2e)
      deallocate(wdy2e)
      deallocate(pb1)
      deallocate(pb2)

  END SUBROUTINE calc_Ez_pssn_y2x


!--------------------------------------
  SUBROUTINE pnl_dft_y2x( dfdv, wkd1, dpot, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv    
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: wkd1 ! dphi/dz or dA/dt in k-space
    real(kind=DP), dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1), intent(out) :: dpot
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: &
         wc1o, wc2o,  wc1e, wc2e
    complex(kind=DP), dimension(:,:,:,:,:), allocatable ::   &
         wwefo, wwefe, wef2e, wef2o
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: wkd2, wkd3
    integer :: mx, my, iz, iv, im

    allocate(wkd2(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wkd3(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    
    allocate(wc1o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wc2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wc1e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wc2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wwefo(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wwefe(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wef2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
    allocate(wef2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
!    
!$OMP parallel default(none)           &
!$OMP shared(dfdv,pnl)       &
!$OMP shared(wc1o,wc2o,wc1e,wc2e)      &
!$OMP shared(wkd1,wkd2,wkd3,dpot)           &
!$OMP shared(wwefo,wwefe,wef2o,wef2e)  &
!$OMP private(iv)

!$OMP workshare
    pnl(:,:,:,:,:)  = (0._DP, 0._DP)
    wc1o(:,:,:,:,:) = (0._DP, 0._DP)
    wc2o(:,:,:,:,:) = (0._DP, 0._DP)
    wc1e(:,:,:,:,:) = (0._DP, 0._DP)
    wc2e(:,:,:,:,:) = (0._DP, 0._DP)
    wwefe(:,:,:,:,:) = (0._DP, 0._DP)
    wef2e(:,:,:,:,:) = (0._DP, 0._DP)
    wwefo(:,:,:,:,:) = (0._DP, 0._DP)
    wef2o(:,:,:,:,:) = (0._DP, 0._DP)

    wkd2(:,:,:,:,:) = (0._DP, 0._DP)
    wkd3(:,:,:,:,:) = (0._DP, 0._DP)
    dpot(:,:,:,:) = 0._DP
!$OMP end workshare

    call pnl_pack_dpot_y2x( wkd1, wkd2 )
!$OMP barrier
    call pnl_transpose_y2x( wkd2, wkd3 )
!$OMP barrier    
    call pnl_backwardfft_y2x( wkd3, dpot ) ! fftw 

    do iv = 1, 2*nv+4
      if (mod(iv,2) == 1) then ! odd
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2x(wc1e,wc2e)       ! 3
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_transpose_x2y(wwefe,wef2e)     ! 5
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2x(iv,dfdv,wc1o)    ! 1
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_realspcal_y2x(dpot,wc2o,wwefo) ! 3
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_unpack_x2y(iv-4,wef2o,pnl)     ! 5
      else                     ! even
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2x(wc1o,wc2o)       ! 2
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_transpose_x2y(wwefo,wef2o)     ! 4       
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2x(iv,dfdv,wc1e)    ! 2
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_realspcal_y2x(dpot,wc2e,wwefe) ! 4
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_unpack_x2y(iv-4,wef2e,pnl)     ! 6
      end if
!$OMP barrier
    end do
!$OMP end parallel

    deallocate( wkd2 )
    deallocate( wkd3 )
    
    deallocate( wc1o )
    deallocate( wc2o )
    deallocate( wc1e )
    deallocate( wc2e )
 
    deallocate( wwefo )
    deallocate( wwefe )
    deallocate( wef2o )
    deallocate( wef2e )
    
  END SUBROUTINE pnl_dft_y2x


!--------------------------------------
  SUBROUTINE pnl_pack_psi_y2x ( psi, wwdx, wwdy ) 
!--------------------------------------
!     For calculating the nonlinear part of E_para (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdx, wwdy

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1a, w2a
    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1b, w2b
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank
    complex(kind=DP) :: ww

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do my = ist_y, iend_y
        do ix = 0, nx
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = psi(ix,my,iz,im)
          enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = nx+1, 2*nxw-nx-1
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = (0._DP, 0._DP) ! FFTW may destroy input array!
          enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = -nx, -1
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,2*nxw+ix) = psi(ix,my,iz,im)
          enddo
          enddo
        enddo
      enddo

      do my = ist_y, iend_y
        do ix = 0, 2*nxw-1
          do im = 0, nm
          do iz = -nz, nz-1
            ww = w1a(iz,im,my,ix)
            w1a(iz,im,my,ix) = uikx_y2x(my,ix) * ww
            w1b(iz,im,my,ix) = uiky_y2x(my,ix) * ww
          enddo
          enddo
        enddo
      enddo
      call dfftw_execute_dft(planr_xb_y2x(ithd), w1a, w2a)
      call dfftw_execute_dft(planr_xb_y2x(ithd), w1b, w2b)

      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                wwdx(iz,im,ix-ist_xw_g_rank,my,irank) = w2a(iz,im,my,ix)
                wwdy(iz,im,ix-ist_xw_g_rank,my,irank) = w2b(iz,im,my,ix)
              enddo
            enddo
          enddo
        enddo
      enddo
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE pnl_pack_psi_y2x


!--------------------------------------
  SUBROUTINE pnl_pssn_realspcal_y2x ( dpdx, dpdy, dadx, dady, wwef ) !done
!--------------------------------------
!     Calculate Poisson brackets for E x B term calculation (y2x)

    real(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpdx, dpdy, dadx, dady
    complex(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:nyw) :: w3
    real(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:2*nyw-1) :: pbxy
    real(kind=DP) :: cef
    integer :: my, ix, iy, iz, im, irank, ist_y_g_rank, iend_y_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

!NEC$ NOINTERCHANGE
      do iy = 0, 2*nyw-1
!NEC$ NOINTERCHANGE
        do ix = ist_xw, iend_xw
!NEC$ NOINTERCHANGE
          do im = 0, nm
!NEC$ NOINTERCHANGE
            do iz = -nz, nz-1
              pbxy(iz,im,ix,iy) = cef * ( & ! Normalization for 2D Forward FFT
                                          dpdx(iz,im,ix,iy) * dady(iz,im,ix,iy) &
                                        - dadx(iz,im,ix,iy) * dpdy(iz,im,ix,iy) )
            end do 
          end do
        end do
      end do

      call dfftw_execute_dft_r2c(planr_yf_y2x(ithd), pbxy, w3)

      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                wwef(iz,im,ix,my-ist_y_g_rank,irank) = w3(iz,im,ix,my)
              end do
            end do
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_pssn_realspcal_y2x
  
!--------------------------------------
  SUBROUTINE pnl_pssn_unpack_x2y ( wwef, pb )
!--------------------------------------
!     Data unpack for pnl pssn term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: pb

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1, w2
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

    ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1450)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                w2(iz,im,my,ix) = wwef(iz,im,ix-ist_xw_g_rank,my,irank)
              enddo
            end do
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xf_y2x(ithd), w2, w1)

      do my = ist_y, iend_y
        do ix = 0, nx
          do im = 0, nm
            do iz = -nz, nz-1
              pb(ix,my,iz,im) = w1(iz,im,my,ix)
            end do
          end do
        end do
      end do
      do my = ist_y, iend_y
        do ix = -nx, -1
          do im = 0, nm
            do iz = -nz, nz-1
              pb(ix,my,iz,im) = w1(iz,im,my,2*nxw+ix)
            end do
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE pnl_pssn_unpack_x2y
  
!--------------------------------------
  SUBROUTINE pnl_pack_dpot_y2x( win, wout )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: win
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wout

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1, w2
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do my = ist_y, iend_y
        do ix = 0, nx
          do im = 0, nm
            do iz = -nz, nz-1
              w1(iz,im,my,ix) = win(ix,my,iz,im)
            enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = nx+1, 2*nxw-nx-1
          do im = 0, nm
            do iz = -nz, nz-1
              w1(iz,im,my,ix) = (0._DP, 0._DP) ! FFTW may destroy input array!
            enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = -nx, -1
          do im = 0, nm
            do iz = -nz, nz-1
              w1(iz,im,my,2*nxw+ix) = win(ix,my,iz,im)
            enddo
          enddo
        enddo
      enddo

      call dfftw_execute_dft(planr_xb_y2x(ithd), w1, w2)

      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                wout(iz,im,ix-ist_xw_g_rank,my,irank) = w2(iz,im,my,ix)
              enddo
            enddo
          enddo
        enddo
      enddo
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  
  END SUBROUTINE pnl_pack_dpot_y2x

!--------------------------------------
  SUBROUTINE pnl_transpose_y2x ( wwin, wwout ) !done
!--------------------------------------
!     Data transpose for E x B term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwout

!$OMP master
                                           call clock_sta(1420)
                                         ! call fapp_start("nlterm_alltoall1",1420,1)
      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall1",1420,1)
                                           call clock_end(1420)
!$OMP end master

  END SUBROUTINE pnl_transpose_y2x

!--------------------------------------
  SUBROUTINE pnl_backwardfft_y2x ( win, wout ) !done
!--------------------------------------
!     Backward FFT of field(psi,chi) for pnl term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: win
    real(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: wout

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:nyw) :: w3
    real(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:2*nyw-1) :: wtmp
    integer :: ix, my, iz, im, irank, ist_y_g_rank, iend_y_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = win(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
        do ix = ist_xw, iend_xw
          do im = 0, nm
            do iz = -nz, nz-1
              w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
            end do
          end do
        end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, wtmp)
      do my = 0, 2*nyw-1
        do ix = ist_xw, iend_xw
          do im = 0, nm
            do iz = -nz, nz-1
              wout(iz,im,ix,my) = wtmp(iz,im,ix,my)
            end do
          end do
        end do
      end do
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_backwardfft_y2x
  
!--------------------------------------
  SUBROUTINE pnl_pack_dfdv_y2x ( iv, dfdv, wwdv ) !done
!--------------------------------------
!     Data pack for E x B term calculation (y2x)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdv

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1, w2
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank
    complex(kind=DP) :: ww

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master

      do my = ist_y, iend_y
       do ix = 0, nx
        do im = 0, nm
          do iz = -nz, nz-1
            w1(iz,im,my,ix) = dfdv(ix,my,iz,iv,im)
          enddo
        enddo
       enddo
      enddo
      do my = ist_y, iend_y
       do ix = nx+1, 2*nxw-nx-1
        do im = 0, nm
          do iz = -nz, nz-1
            w1(iz,im,my,ix) = (0._DP, 0._DP) ! FFTW may destroy input array!
          enddo
        enddo
       enddo
      enddo
      do my = ist_y, iend_y
       do ix = -nx, -1
        do im = 0, nm
          do iz = -nz, nz-1
            w1(iz,im,my,2*nxw+ix) = dfdv(ix,my,iz,iv,im)
          enddo
        enddo
       enddo
      enddo

      call dfftw_execute_dft(planr_xb_y2x(ithd), w1, w2)

      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                wwdv(iz,im,ix-ist_xw_g_rank,my,irank) = w2(iz,im,my,ix)
              enddo
            enddo
          enddo
        enddo
      enddo
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE pnl_pack_dfdv_y2x

!--------------------------------------
  SUBROUTINE pnl_realspcal_y2x ( dpot, wwdv, wwef ) 
!--------------------------------------

    real(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpot
    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdv
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:nyw) :: w3
    real(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:2*nyw-1) :: dfdv, pnl
    real(kind=DP) :: cef, cs1
    integer :: my, ix, iy, iz, im, irank, ist_y_g_rank, iend_y_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP) ! Normalization for 2D Forward FFT
      cs1 = sgn(ranks) * Znum(ranks) / Anum(ranks)
    
!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = wwdv(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
       do ix = ist_xw, iend_xw
        do im = 0, nm
          do iz = -nz, nz-1
            w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
          end do
        end do
       end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, dfdv)

!NEC$ NOINTERCHANGE
      do iy = 0, 2*nyw-1
!NEC$ NOINTERCHANGE
        do ix = ist_xw, iend_xw
!NEC$ NOINTERCHANGE
          do im = 0, nm
!NEC$ NOINTERCHANGE
            do iz = -nz, nz-1
              pnl(iz,im,ix,iy) = - cef * &
                                   cs1 * dpot(iz,im,ix,iy) * dfdv(iz,im,ix,iy) 
            end do
          end do
        end do
      end do

      call dfftw_execute_dft_r2c(planr_yf_y2x(ithd), pnl, w3)

      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                wwef(iz,im,ix,my-ist_y_g_rank,irank) = w3(iz,im,ix,my)
              end do
            end do
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_realspcal_y2x

                                        
!--------------------------------------
  SUBROUTINE pnl_transpose_x2y ( wwin, wwout )
!--------------------------------------
!     Data transpose for pnl term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwout

!$OMP master
                                           call clock_sta(1440)
                                         ! call fapp_start("nlterm_alltoall2",1440,1)
      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall2",1440,1)
                                           call clock_end(1440)
!$OMP end master

  END SUBROUTINE pnl_transpose_x2y
  
!--------------------------------------
  SUBROUTINE pnl_unpack_x2y ( iv, wwef, pnl ) !done
!--------------------------------------
!     Data unpack for pnl term calculation (y2x)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl
    !NOTE: A noncontiguous subarray as an argument of a subroutine induces a memory copy.
    !      When the subroutine is called in a OpenMP parallel region,
    !      the copied subarray may be treated as a thread-private variable.

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1, w2
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

    ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1450)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                w2(iz,im,my,ix) = wwef(iz,im,ix-ist_xw_g_rank,my,irank)
              enddo
            end do
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xf_y2x(ithd), w2, w1)

      do my = ist_y, iend_y
        do ix = 0, nx
          do im = 0, nm
            do iz = -nz, nz-1
              pnl(ix,my,iz,iv,im) = w1(iz,im,my,ix)
            end do
          end do
        end do
      end do
      do my = ist_y, iend_y
        do ix = -nx, -1
          do im = 0, nm
            do iz = -nz, nz-1
              pnl(ix,my,iz,iv,im) = w1(iz,im,my,2*nxw+ix)
            end do
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE pnl_unpack_x2y


!--------------------------------------
  SUBROUTINE pnl_estimate_maxvel_y2x ( dpdz, dAdt ) 
!--------------------------------------
!     Estimate time step restriction in each MPI processes

    real(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpdz, dAdt

    real(kind=DP) :: cs1, wv_nl
    integer :: ix, iy, iz, im

      pnl_max_eachrank = eps

      cs1 = dref * Znum(ranks) / sqrt( Anum(ranks) * tau(ranks) )
      im=0
!$OMP parallel default(none)                                  &
!$OMP shared(pnl_max_eachrank)           &
!$OMP shared(im,ist_xw,iend_xw,dpdz,dAdt,cs1) &
!$OMP private(iy,ix,iz,wv_nl)

!$OMP do collapse(2) reduction(max:pnl_max_eachrank)
      do iy = 0, 2*nyw-1
         do ix = ist_xw, iend_xw
            do iz = -nz, nz-1
               wv_nl = cs1 * abs( dpdz(iz,im,ix,iy) + dAdt(iz,im,ix,iy) )
               if (pnl_max_eachrank < wv_nl) pnl_max_eachrank = wv_nl              
            end do
         end do
      end do
!$OMP end do

!$OMP end parallel


  END SUBROUTINE pnl_estimate_maxvel_y2x
!!$ end of y2x part

!!$ beginning of y2zm part
!--------------------------------------
  SUBROUTINE CALC_PNL_term_y2zm( ff, psi, dh, cf, ef )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: dfdv
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: pnl_e, pnl_m, df_st

    real(kind=DP), allocatable, dimension(:,:,:) :: dpdz, dAdt
    real(kind=DP) :: cefv
    integer :: mx, my, iz, iv, im
    
    allocate(dfdv(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dpdz(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
    allocate(dAdt(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
    allocate(pnl_e(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(pnl_m(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))

    call bndry_bound_e( psi )

! M-I coupling part
    if ( trim(z_bound) == "mi_couple"  ) then
      do im = 0, nm
        call micouple_bndry_bound_e ( psi(:,:,:,im) )
      end do
    end if
!M-I coupling part

!...zero reset
    dpdz(:,:,:)    = 0._DP
    dAdt(:,:,:)    = 0._DP    
    dfdv(:,:,:,:,:)  = ( 0._DP, 0._DP ) 
    pnl_e(:,:,:,:,:) = ( 0._DP, 0._DP )   
    pnl_m(:,:,:,:,:) = ( 0._DP, 0._DP )

!...d/dv element
    cefv  = sqrt( Anum(ranks) / tau(ranks) ) / ( 12._DP * dv )

!...make df/dv
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do iv = 1, 2*nv
              dfdv(mx,my,iz,iv,im) = cefv * ( - ff(mx,my,iz,iv+2,im)      &
                                      + 8._DP * ff(mx,my,iz,iv+1,im)      &
                                      - 8._DP * ff(mx,my,iz,iv-1,im)      &
                                              + ff(mx,my,iz,iv-2,im) )
            end do
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
       dfdv(0,0,:,:,:) = ( 0._DP, 0._DP )
    end if

    call PNL_es_term_y2zm( psi, dfdv, dpdz, pnl_e ) ! electrostatic part

    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              ef(mx,my,iz,iv,im) = ef(mx,my,iz,iv,im) + dref * pnl_e(mx,my,iz,iv,im)    
            end do
          end do
        end do
      end do
    end do
 
    if( beta>0._DP ) then
      call PNL_em_term_y2zm( ff, dfdv, dh, cf, ef, dAdt, pnl_m ) ! magnetic part
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ef(mx,my,iz,iv,im) = ef(mx,my,iz,iv,im) + dref * pnl_m(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do 
      end do
    end if
   
!...debug
    !iz=nz
    !iv=2*nv
    !im=0   
    !do my=0,ny
    !   do mx=0,nx
    !     write(olog, '(2i6, 99ES16.6)' ) mx, my, pnl_e(mx,my,iz,iv,im), pnl_m(mx,my,iz,iv,im)
    !  end do
    !end do
    !write(olog, *)
    !flush(olog)

!...estimate time step size 
    call pnl_estimate_maxvel_y2zm(dpdz,dAdt)

    deallocate( dfdv )
    deallocate( dpdz )
    deallocate( dAdt )
    deallocate( pnl_e )
    deallocate( pnl_m )
    
  END SUBROUTINE CALC_PNL_term_y2zm

!--------------------------------------
  SUBROUTINE PNL_es_term_y2zm( psi, dfdv, dpdz, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl
    real(kind=DP), dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1), intent(out) :: dpdz

    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
         wc1o, wc2o, wc3o, wc4o,   wc1e, wc2e, wc3e, wc4e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdvo, wwdve,    wwefo, wwefe
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdz1, wdz2, wdz3
    complex(kind=DP), allocatable, dimension(:,:,:) :: wwdz
    real(kind=DP), dimension(-nz:nz-1) :: cefz
    integer :: mx, my, iz, iv, im, iprocw

    allocate(wc1o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc3o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc4o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc1e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc3e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc4e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwdvo(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwdve(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwefo(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwefe(0:global_ny,0:2*nxw-1,0:nbuff-1))    
!      
    allocate(wdz1(-nx:nx,0:ny,-nz:nz-1,0:nm))
    allocate(wdz2(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wdz3(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwdz(0:global_ny,0:2*nxw-1,0:nbuff-1))
    
    do iz = -nz, nz-1
      cefz(iz)  = 1._DP / ( 12._DP * dpara(iz) )
    end do

    wdz1(:,:,:,:) = ( 0._DP, 0._DP )

!...make dphi/dz
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx

            wdz1(mx,my,iz,im) = cefz(iz) * ( - psi(mx,my,iz+2,im)    &
                                     + 8._DP * psi(mx,my,iz+1,im)    &
                                     - 8._DP * psi(mx,my,iz-1,im)    &
                                             + psi(mx,my,iz-2,im) )
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
      wdz1(0,0,:,:) = ( 0._DP, 0._DP )
    end if    
    
!$OMP parallel default(none)           &
!$OMP shared(dfdv,pnl)                 &
!$OMP shared(wc1o,wc2o,wc3o,wc4o)      &
!$OMP shared(wc1e,wc2e,wc3e,wc4e)      &
!$OMP shared(wwdvo,wwdve,wwefo,wwefe)  &
!$OMP shared(wdz1,wdz2,wdz3,wwdz,dpdz) &
!$OMP private(iv,iprocw)
    
!$OMP workshare
      pnl(:,:,:,:,:) = (0._DP, 0._DP)
      wc1o(:,:,:,:)  = (0._DP, 0._DP)
      wc1e(:,:,:,:)  = (0._DP, 0._DP)
      wc2o(:,:,:,:)  = (0._DP, 0._DP)
      wc2e(:,:,:,:)  = (0._DP, 0._DP)      
      wc3o(:,:,:,:)  = (0._DP, 0._DP)
      wc3e(:,:,:,:)  = (0._DP, 0._DP)
      wc4o(:,:,:,:)  = (0._DP, 0._DP)
      wc4e(:,:,:,:)  = (0._DP, 0._DP)
      
      wdz2(:,:,:,:) = (0._DP, 0._DP)
      wdz3(:,:,:,:) = (0._DP, 0._DP)
      wwdz(:,:,:)   = (0._DP, 0._DP)
      dpdz(:,:,:)   = 0._DP

      wwdve(:,:,:) = (0._DP, 0._DP)
      wwdvo(:,:,:) = (0._DP, 0._DP)      
      wwefe(:,:,:) = (0._DP, 0._DP)
      wwefo(:,:,:) = (0._DP, 0._DP)      
!$OMP end workshare
      
    call pnl_pack_dpot_y2zm( wdz1, wdz2 )
!$OMP barrier
    call pnl_transpose_y2zm( wdz2, wdz3 )
!$OMP barrier
    call pnl_unpack_dpot_y2zm( wdz3, wwdz ) ! fftw in x-direction
!$OMP barrier
    call pnl_backwardfft_y2zm( wwdz, dpdz ) ! fftw in y-direction

    do iv = 1, 2*nv+6
      if (mod(iv,2) == 1) then ! odd
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1e,wc2e)        ! 3
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3e,wc4e)        ! 7
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1o)     ! 1
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2o,wwdvo)          ! 3
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpdz,wwdve,wwefe) ! 5
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwefo,wc3o)            ! 5
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4o,pnl)       ! 7
      else                     ! even
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1o,wc2o)        ! 2
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3o,wc4o)        ! 6
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1e)     ! 2
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2e,wwdve)          ! 4
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpdz,wwdvo,wwefo) ! 4
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwefe,wc3e)            ! 6
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4e,pnl)       ! 8
      end if
!$OMP barrier
    end do
!$OMP end parallel

    deallocate( wc1o )
    deallocate( wc2o )
    deallocate( wc3o )
    deallocate( wc4o )
    deallocate( wc1e )
    deallocate( wc2e )
    deallocate( wc3e )
    deallocate( wc4e )

    deallocate( wwdvo )
    deallocate( wwdve )
    deallocate( wwefo )
    deallocate( wwefe )

    deallocate( wdz1 )
    deallocate( wdz2 )
    deallocate( wdz3 )
    deallocate( wwdz )

  END SUBROUTINE PNL_es_term_y2zm

!--------------------------------------
  SUBROUTINE PNL_em_term_y2zm( ff, dfdv, dh, cf, ef, dAdt, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv     
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf, ef
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl   
    real(kind=DP), dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1), intent(out) :: dAdt
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: df_st!!, cefX
    complex(kind=DP), allocatable, dimension(:,:,:) :: dpsi !!X1, X2

    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
         wc1o, wc2o, wc3o, wc4o,   wc1e, wc2e, wc3e, wc4e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdvo, wwdve,    wwefo, wwefe
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdt1, wdt2, wdt3
    complex(kind=DP), allocatable, dimension(:,:,:) :: wwdt   
    
    integer :: mx, my, iz, iv, im, iprocw
    
    allocate(df_st(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dpsi(-nx:nx,0:ny,-nz:nz-1))    
!    allocate(X1(-nx:nx,0:ny,-nz:nz-1))
!    allocate(X2(-nx:nx,0:ny,-nz:nz-1))
!    allocate(cefX(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    
    allocate(wc1o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc3o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc4o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc1e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc3e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc4e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwdvo(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwdve(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwefo(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwefe(0:global_ny,0:2*nxw-1,0:nbuff-1))    
!      
    allocate(wdt1(-nx:nx,0:ny,-nz:nz-1,0:nm))
    allocate(wdt2(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wdt3(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwdt(0:global_ny,0:2*nxw-1,0:nbuff-1))

    df_st(:,:,:,:,:) = ( 0._DP, 0._DP )
    !!cefX(:,:,:,:,:)  = ( 0._DP, 0._DP ) 
   
    ! df w/o dA/dt term
    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              df_st(mx,my,iz,iv,im) = dh(mx,my,iz,iv,im) + cf(mx,my,iz,iv,im) - ef(mx,my,iz,iv,im)
              !df_st(mx,my,iz,iv,im) = j0(mx,my,iz,im) * sgn(ranks) * Znum(ranks) *                 &
              !                        sqrt( tau(ranks) / Anum(ranks) ) * vl(iv) *                  &
              !                        (                                                            &
              !                        dh(mx,my,iz,iv,im) + cf(mx,my,iz,iv,im) - ef(mx,my,iz,iv,im) &
              !                        )
              !cefX(mx,my,iz,iv,im) =  ( Znum(ranks)**2 / Anum(ranks) ) *         &
              !                        j0(mx,my,iz,im) ** 2 *                                    &  
              !                        (                                                         &
              !                        vl(iv) ** 2 * fmx(iz,iv,im)  + dref * ff(mx,my,iz,iv,im) &
              !                        ) 
            end do
          end do
        end do
      end do
    end do    

    call fld_emfield_hh ( df_st, dpsi )
    
    !call intgrl_v0_moment_ms( df_st, X1 )
    !call intgrl_v0_moment_ms( cefX,  X2 )

    ! J_0 * dA/dt
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
             wdt1(mx,my,iz,im) = j0(mx,my,iz,im) * dpsi(mx,my,iz)
            !!wdt1(mx,my,iz,im) = j0(mx,my,iz,im) * beta * X1(mx,my,iz) / ( beta * X2(mx,my,iz) + ksq(mx,my,iz) )
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
      wdt1(0,0,:,:) = ( 0._DP, 0._DP )
    end if
      
!$OMP parallel default(none)           &
!$OMP shared(dfdv,pnl)                 &
!$OMP shared(wc1o,wc2o,wc3o,wc4o)      &
!$OMP shared(wc1e,wc2e,wc3e,wc4e)      &
!$OMP shared(wwdvo,wwdve,wwefo,wwefe)  &
!$OMP shared(wdt1,wdt2,wdt3,wwdt,dAdt) &
!$OMP private(iv,iprocw)

!$OMP workshare
      pnl(:,:,:,:,:) = (0._DP, 0._DP)
      wc1o(:,:,:,:)  = (0._DP, 0._DP)
      wc1e(:,:,:,:)  = (0._DP, 0._DP)
      wc2o(:,:,:,:)  = (0._DP, 0._DP)
      wc2e(:,:,:,:)  = (0._DP, 0._DP)      
      wc3o(:,:,:,:)  = (0._DP, 0._DP)
      wc3e(:,:,:,:)  = (0._DP, 0._DP)
      wc4o(:,:,:,:)  = (0._DP, 0._DP)
      wc4e(:,:,:,:)  = (0._DP, 0._DP)
      
      wdt2(:,:,:,:)  = (0._DP, 0._DP)
      wdt3(:,:,:,:)  = (0._DP, 0._DP)
      wwdt(:,:,:) = (0._DP, 0._DP)
      dAdt(:,:,:) = 0._DP

      wwdve(:,:,:) = (0._DP, 0._DP)
      wwdvo(:,:,:) = (0._DP, 0._DP)      
      wwefe(:,:,:) = (0._DP, 0._DP)
      wwefo(:,:,:) = (0._DP, 0._DP)      
!$OMP end workshare
      
    call pnl_pack_dpot_y2zm( wdt1, wdt2 )
!$OMP barrier
    call pnl_transpose_y2zm( wdt2, wdt3 )
!$OMP barrier
    call pnl_unpack_dpot_y2zm( wdt3, wwdt ) ! fftw in x-direction
!$OMP barrier
    call pnl_backwardfft_y2zm( wwdt, dAdt ) ! fftw in y-direction
      
    do iv = 1, 2*nv+6
      if (mod(iv,2) == 1) then ! odd
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1e,wc2e)        ! 3
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3e,wc4e)        ! 7
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1o)     ! 1
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2o,wwdvo)          ! 3
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dAdt,wwdve,wwefe) ! 5
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwefo,wc3o)            ! 5
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4o,pnl)       ! 7
      else                     ! even
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1o,wc2o)        ! 2
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3o,wc4o)        ! 6
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1e)     ! 2
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2e,wwdve)          ! 4
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dAdt,wwdvo,wwefo) ! 4
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwefe,wc3e)            ! 6
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4e,pnl)       ! 8
      end if
!$OMP barrier
    end do
!$OMP end parallel
    
    deallocate( wc1o )
    deallocate( wc2o )
    deallocate( wc3o )
    deallocate( wc4o )
    deallocate( wc1e )
    deallocate( wc2e )
    deallocate( wc3e )
    deallocate( wc4e )

    deallocate( wwdvo )
    deallocate( wwdve )
    deallocate( wwefo )
    deallocate( wwefe )

    deallocate( wdt1 )
    deallocate( wdt2 )
    deallocate( wdt3 )
    deallocate( wwdt )
    
    deallocate( df_st )
    deallocate( dpsi )
    !deallocate( X1 )
    !deallocate( X2 )
    !deallocate( cefX )
    
  END SUBROUTINE PNL_em_term_y2zm


!--------------------------------------
  SUBROUTINE pnl_pack_dpot_y2zm( win, wout )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: win
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wout

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1410)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              wout(mx,my,ibuff,iprocw) = win(mx,my,iz,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE pnl_pack_dpot_y2zm


!--------------------------------------
  SUBROUTINE pnl_transpose_y2zm ( wc4in, wc4out )
!--------------------------------------
!     Data transpose for pnl term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4in
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4out

!$OMP master
                                           call clock_sta(1420)
                                         ! call fapp_start("nlterm_alltoall1",1420,1)
      call MPI_Alltoall( wc4in,                 &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         wc4out,                &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall1",1420,1)
                                           call clock_end(1420)
!$OMP end master

  END SUBROUTINE pnl_transpose_y2zm

!--------------------------------------
  SUBROUTINE pnl_unpack_dpot_y2zm( wc4, wout )
!--------------------------------------
!     Data unpack for pnl term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wout

    complex(kind=DP), dimension(-nx:nx) :: dpot
    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    integer :: my, ibuff, iprocw, global_my

!$OMP master
                                           call clock_sta(1430)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_yb)
      do ibuff = 0, nbuff-1
        do global_my = 0, global_ny

         !%%% UNPACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          iprocw = global_my / (ny+1)
          my = mod(global_my, ny+1)
          dpot(-nx:nx) = wc4(-nx:nx,my,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Backward x-FFT (kx,ky)->(ky,x) %%%
          w1(0:nx) = dpot(0:nx)
          w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
          w1(2*nxw-nx:2*nxw-1) = dpot(-nx:-1)
          call dfftw_execute_dft(plan_x_backward, w1, w2)
          wout(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_unpack_dpot_y2zm


!--------------------------------------
  SUBROUTINE pnl_backwardfft_y2zm ( win, dpot )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: win
    real(kind=DP), intent(out), &
      dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1) :: dpot

    complex(kind=DP), dimension(0:nyw) :: w3
    integer :: ix, ibuff

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_xb)
      do ibuff = 0, nbuff-1
        do ix = 0, 2*nxw-1

         !%%% Backward y-FFT (ky,x)->(y,x) %%%
          w3(0:global_ny) = win(0:global_ny,ix,ibuff)
          w3(global_ny+1:nyw) = (0._DP, 0._DP) ! FFTW may destroy input array!
          call dfftw_execute_dft_c2r(plan_y_backward, w3, dpot(:,ix,ibuff))
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master
    
  END SUBROUTINE pnl_backwardfft_y2zm

!--------------------------------------
  SUBROUTINE pnl_pack_dfdv_y2zm ( iv, wdv, wc4 )
!--------------------------------------
!
    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wdv
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1410)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              wc4(mx,my,ibuff,iprocw) = wdv(mx,my,iz,iv,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE pnl_pack_dfdv_y2zm

!--------------------------------------
  SUBROUTINE pnl_unpack_y2zm ( wc4, wwdv )
!--------------------------------------
!
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwdv

    complex(kind=DP), dimension(-nx:nx) :: psi
    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    integer :: my, ibuff, iprocw, global_my

!$OMP master
                                           call clock_sta(1430)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_yb)
      do ibuff = 0, nbuff-1
        do global_my = 0, global_ny

         !%%% UNPACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          iprocw = global_my / (ny+1)
          my = mod(global_my, ny+1)
          psi(-nx:nx) = wc4(-nx:nx,my,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Backward x-FFT (kx,ky)->(ky,x) %%%
          w1(0:nx) = psi(0:nx)
          w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
          w1(2*nxw-nx:2*nxw-1) = psi(-nx:-1)
          call dfftw_execute_dft(plan_x_backward, w1, w2)
          wwdv(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_unpack_y2zm

!--------------------------------------
  SUBROUTINE pnl_realspcal_y2zm ( dpot, wwdv, wwef )
!--------------------------------------

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1) :: dpot
    complex(kind=DP), intent(in), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwdv
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwef

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) :: dfdv, pnl
    real(kind=DP) :: cef, cs1
    integer :: ix, iy, ibuff

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP) ! Normalization for 2D Forward FFT
      cs1 = sgn(ranks) * Znum(ranks) / Anum(ranks)

!$OMP master
                                           call clock_sta(1430)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_xb)
      do ibuff = 0, nbuff-1
        do ix = 0, 2*nxw-1

         !%%% Backward y-FFT (ky,x)->(y,x) %%%
          w3(0:global_ny) = wwdv(0:global_ny,ix,ibuff)
          w3(global_ny+1:nyw) = (0._DP, 0._DP) ! FFTW may destroy input array!
          call dfftw_execute_dft_c2r(plan_y_backward, w3, dfdv)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% dpdz*dfdv in (y,x) %%%
          do iy = 0, 2*nyw-1
            pnl(iy) = - cef * & 
                        cs1 * dpot(iy,ix,ibuff) * dfdv(iy)
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Forward y-FFT (y,x)->(ky,x) %%%
          call dfftw_execute_dft_r2c(plan_y_forward, pnl, w3)
          wwef(0:global_ny,ix,ibuff) = w3(0:global_ny)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_realspcal_y2zm

!--------------------------------------
  SUBROUTINE pnl_pack_zm2y ( wwef, wc4 )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwef
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4

    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    complex(kind=DP), dimension(-nx:nx) :: ef
    integer :: mx, my, ibuff, iprocw, global_my

!$OMP master
                                           call clock_sta(1430)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_yb)
      do ibuff = 0, nbuff-1
        do global_my = 0, global_ny

         !%%% Forward x-FFT (ky,x)->(kx,ky) %%%
          w2(0:2*nxw-1) = wwef(global_my,0:2*nxw-1,ibuff) ! FFTW may destroy input array!
          call dfftw_execute_dft(plan_x_forward, w2, w1)
          ef(0:nx) = w1(0:nx)
          ef(-nx:-1) = w1(2*nxw-nx:2*nxw-1)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% PACK: (kx,ky,(z*,m*)*)->(kx,ky*,z*,m*) %%%
          iprocw = global_my / (ny+1)
          my = mod(global_my, ny+1)
          do mx = -nx, nx
            wc4(mx,my,ibuff,iprocw) = ef(mx)
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_pack_zm2y

!--------------------------------------
  SUBROUTINE pnl_transpose_zm2y ( wc4in, wc4out )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4in
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4out

!$OMP master
                                           call clock_sta(1440)
      call MPI_Alltoall( wc4in,                 &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         wc4out,                &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
                                           call clock_end(1440)
!$OMP end master

  END SUBROUTINE pnl_transpose_zm2y

!--------------------------------------
  SUBROUTINE pnl_unpack_zm2y ( iv, wc4, pnl )
!--------------------------------------

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl

    integer :: iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1450)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% UNPACK: (kx,ky,(z*,m*)*)->(kx,ky*,z*,m*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          pnl(:,:,iz,iv,im) = wc4(:,:,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE pnl_unpack_zm2y


!--------------------------------------
  SUBROUTINE pnl_estimate_maxvel_y2zm ( dpdz, dAdt )
!--------------------------------------
!     Estimate time step restriction in each MPI processes

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1) :: dpdz, dAdt

    real(kind=DP) :: cs1, wv_nl
    integer :: ix, iy, ibuff

      pnl_max_eachrank = eps
      
      cs1 = dref * Znum(ranks) / sqrt( Anum(ranks) * tau(ranks) )

!$OMP parallel default(none)                               &
!$OMP shared(pnl_max_eachrank)        &
!$OMP shared(ist_xw,iend_xw,dpdz,dAdt,cs1,vl) &
!$OMP private(iy,ix,ibuff,wv_nl)

!$OMP do collapse(2) reduction(max:pnl_max_eachrank)
        do ibuff = 0, nbuff-1
          do ix = 0, 2*nxw-1
            do iy = 0, 2*nyw-1
              wv_nl = cs1 * abs(dpdz(iy,ix,ibuff) + dAdt(iy,ix,ibuff))
              if (pnl_max_eachrank < wv_nl) pnl_max_eachrank = wv_nl
            end do
          end do
        end do
!$OMP end do

!$OMP end parallel

  END SUBROUTINE pnl_estimate_maxvel_y2zm

!--------------------------------------
  SUBROUTINE pnl_sum ( ff, phi, Al, dpdz, dh, pnint )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: dpdz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al  
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny) :: pnint
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: pnl
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: dfdv
    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi
    
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf, wef, wbf, wpp
    complex(kind=DP), dimension(:,:,:,:), allocatable :: Ezk
    complex(kind=DP), dimension(:,:,:), allocatable :: dAdt    
    complex(kind=DP), dimension(:,:,:), allocatable :: wc3
    complex(kind=DP), dimension(-nx:nx,0:ny) :: wc2
    real(kind=DP), allocatable, dimension(:,:,:,:) :: Ezr
    real(kind=DP) :: cefv
    integer  ::  mx, my, iz, iv, im
    
    allocate(dfdv(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dAdt(-nx:nx,0:ny,-nz:nz-1))
    allocate(Ezk(-nx:nx,0:ny,-nz:nz-1,0:nm))
    allocate(Ezr(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
    allocate(pnl(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(psi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm))
    allocate(chi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm))

    allocate(wc3(-nx:nx,0:ny,-nz:nz-1))

!$OMP parallel workshare
    Ezk(:,:,:,:) = (0._DP, 0._DP)
    dAdt(:,:,:)  = (0._DP, 0._DP)
    pnint(:,:)   = 0._DP
!$OMP end parallel workshare

    cefv  = sqrt( Anum(ranks) / tau(ranks) ) / ( 12._DP * dv )
    if( beta>0._DP ) then
      call fld_emfield_hh ( dh, dAdt )

!$OMP parallel default(none) &
!$OMP shared(psi,chi,phi,Al,j0,ist_y,iend_y) &
!$OMP private(mx,my,iz,im)
!$OMP do collapse(2)
      do im = 0, nm
        do iz = -nz, nz-1
          do my = iend_y, ny
            psi(:,my,iz,im) = (0._DP, 0._DP)
            chi(:,my,iz,im) = (0._DP, 0._DP)
          end do
          do my = ist_y, iend_y
            do mx = -nx, nx
              psi(mx,my,iz,im) = j0(mx,my,iz,im) * phi(mx,my,iz)
              chi(mx,my,iz,im) = j0(mx,my,iz,im) * Al(mx,my,iz)
            end do
          end do
        end do
      end do
!$OMP end do
!$OMP end parallel
     
!...make nonlinear part of Ez
      call calc_Ez_pssn_y2x( psi, chi, Ezk )
    end if

!...make df/dv
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx            
            do iv = 1, 2*nv
              dfdv(mx,my,iz,iv,im) = cefv * ( - ff(mx,my,iz,iv+2,im)      &
                                      + 8._DP * ff(mx,my,iz,iv+1,im)      &
                                      - 8._DP * ff(mx,my,iz,iv-1,im)      &
                                              + ff(mx,my,iz,iv-2,im) )
            end do

            ! note that the sign will be inverted inside the subroutine below
            Ezk(mx,my,iz,im) = dpdz(mx,my,iz,im) + j0(mx,my,iz,im) * dAdt(mx,my,iz) &
                             + Ezk(mx,my,iz,im) 

          end do
        end do
      end do
    end do

    if( rankw==0 ) then
      dfdv(0,0,:,:,:) = (0._DP, 0._DP)
      Ezk(0,0,:,:)    = (0._DP, 0._DP)
    end if
   
    call pnl_dft_y2x( dfdv, Ezk, Ezr, pnl )
 
!$OMP parallel
    do im = 0, nm
!$OMP do
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              pnl(mx,my,iz,iv,im) = dref * pnl(mx,my,iz,iv,im) * tau(ranks) * &
                                    conjg( ff(mx,my,iz,iv,im) ) / fmx(iz,iv,im)
            end do
          end do
        end do
      end do
!$OMP end do nowait
    end do
!$OMP end parallel

    call intgrl_v0_moment ( pnl, wc3 )
    call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
    do my = ist_y, iend_y
      do mx = -nx, nx
        pnint(mx,my) = real( wc2(mx,my), kind=DP )
      end do
    end do
      
    deallocate(dfdv)
    deallocate(dAdt)
    deallocate(Ezk)
    deallocate(Ezr)
    deallocate(pnl)
    deallocate(psi)
    deallocate(chi)
    deallocate(wc3)

  END SUBROUTINE pnl_sum
   
END MODULE GKV_pnl

