MODULE GKV_pnl

! calculare parallel nonlinear term (PNL)
  
  use GKV_header
  use GKV_mpienv
  use GKV_intgrl, only: intgrl_v0_moment_ms, intgrl_thet, intgrl_v0_moment
  use GKV_fft, only: &
           plan_x_forward, plan_x_backward, & ! for fugaku
           plan_y_forward, plan_y_backward    ! for fugaku
  use GKV_clock, only: clock_sta, clock_end
  use GKV_micouple, only: micouple_bndry_bound_e
  use GKV_fld, only: fld_emfield_hh
  use GKV_bndry, only: bndry_bound_e

  implicit none

  private
  real(kind=DP), save :: pnl_max_eachrank
  real(kind=DP), dimension(0:global_ny), save :: gky

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

    real(kind=DP) :: dky
    integer, save :: iflg
    integer :: my
!$  integer :: nthreads, omp_get_num_threads
    data iflg / 0 /
    
    if( iflg == 0 ) then
      iflg = 1
      dky = ky(1) - ky(0)
      do my = 0, global_ny
        gky(my) = dky * real(my, kind=DP)
      end do
!$OMP parallel default(shared)
!$OMP master
!$    nthreads = omp_get_num_threads()
!$    if (nthreads > 1) then
!$      nchunk_zm = ((2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$      nchunk_yb = ((global_ny+1)*nbuff-1) / (nthreads-1) + 1
!$      nchunk_xb = ((2*nxw)*nbuff-1) / (nthreads-1) + 1
!$      nchunk_yzm = ((iend_y-ist_y+1)*(2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$      nchunk_xzm = ((iend_xw-ist_xw+1)*(2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$    end if
!$OMP end master
!$OMP end parallel
    end if

    if (trim(calc_type) == "nonlinear") then
       call CALC_PNL_term_y2zm( ff, psi, chi, dh, cf, ef ) ! fugaku
    else
       ef(:,:,:,:,:) = ( 0._DP, 0._DP )
    end if

    return

  END SUBROUTINE pnl_term

!!$ beginning of y2zm part
!--------------------------------------
  SUBROUTINE CALC_PNL_term_y2zm( ff, psi, chi, dh, cf, ef )
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

    call PNL_es_term_y2zm( psi, chi, dfdv, dpdz, pnl_e ) ! electrostatic part

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
  SUBROUTINE PNL_es_term_y2zm( psi, chi, dfdv, dpdz, pnl )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl
    real(kind=DP), dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1), intent(out) :: dpdz
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdz
    real(kind=DP), dimension(-nz:nz-1) :: cefz
    integer :: mx, my, iz, iv, im
!      
    allocate(wdz(-nx:nx,0:ny,-nz:nz-1,0:nm))
     
    do iz = -nz, nz-1
      cefz(iz)  = 1._DP / ( 12._DP * dpara(iz) )
    end do

    wdz(:,:,:,:) = ( 0._DP, 0._DP )

!...make nonlinear part
    if( beta>0._DP ) then
      call calc_Ez_pssn_y2zm( psi, chi, wdz )
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

    call pnl_dft_y2zm( dfdv, wdz, dpdz, pnl )
    
    deallocate( wdz )

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
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: df_st, cefX
    complex(kind=DP), allocatable, dimension(:,:,:) :: X1, X2
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdt 
    
    integer :: mx, my, iz, iv, im
    
    allocate(df_st(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(X1(-nx:nx,0:ny,-nz:nz-1))
    allocate(X2(-nx:nx,0:ny,-nz:nz-1))
    allocate(cefX(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    
    allocate(wdt(-nx:nx,0:ny,-nz:nz-1,0:nm))

    df_st(:,:,:,:,:) = ( 0._DP, 0._DP )
    cefX(:,:,:,:,:)  = ( 0._DP, 0._DP ) 
   
    ! df w/o dA/dt term
    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              df_st(mx,my,iz,iv,im) = j0(mx,my,iz,im) * sgn(ranks) * Znum(ranks) *                 &
                                      sqrt( tau(ranks) / Anum(ranks) ) * vl(iv) *                  &
                                      (                                                            &
                                      dh(mx,my,iz,iv,im) + cf(mx,my,iz,iv,im) - ef(mx,my,iz,iv,im) &
                                      )
              cefX(mx,my,iz,iv,im) =  ( Znum(ranks)**2 / Anum(ranks) ) *         &
                                      j0(mx,my,iz,im) ** 2 *                                    &  
                                      (                                                         &
                                      vl(iv) ** 2 * fmx(iz,iv,im)  + dref * ff(mx,my,iz,iv,im)  &
                                      ) 
            end do
          end do
        end do
      end do
    end do    

    call intgrl_v0_moment_ms( df_st, X1 )
    call intgrl_v0_moment_ms( cefX,  X2 )

    ! J_0 * dA/dt
    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            wdt(mx,my,iz,im) = j0(mx,my,iz,im) * beta * X1(mx,my,iz) / ( beta * X2(mx,my,iz) + ksq(mx,my,iz) )
          end do
        end do
      end do
    end do

    if( rankw==0 ) then
      wdt(0,0,:,:) = ( 0._DP, 0._DP )
    end if

    call pnl_dft_y2zm( dfdv, wdt, dAdt, pnl )

    deallocate( wdt )
    
    deallocate( df_st )
    deallocate( X1 )
    deallocate( X2 )
    deallocate( cefX )
    
  END SUBROUTINE PNL_em_term_y2zm

!--------------------------------------
  SUBROUTINE calc_Ez_pssn_y2zm( psi, chi, pssn )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: pssn

    real(kind=DP), dimension(:,:,:), allocatable :: dpdx, dpdy, dadx, dady
    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
                         wc1o, wc2o,   wc1e, wc2e,          &
                         pb2, pb3
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdxo, wwdyo, wwefo,      wwdxe, wwdye, wwefe,   pb1
 
      allocate(dpdx(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
      allocate(dpdy(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
      allocate(dadx(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
      allocate(dady(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
      allocate(wc1o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
      allocate(wc2o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
      allocate(wc1e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
      allocate(wc2e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
      allocate(wwdxo(0:global_ny,0:2*nxw-1,0:nbuff-1))
      allocate(wwdyo(0:global_ny,0:2*nxw-1,0:nbuff-1))
      allocate(wwdxe(0:global_ny,0:2*nxw-1,0:nbuff-1))
      allocate(wwdye(0:global_ny,0:2*nxw-1,0:nbuff-1))
      allocate(pb1(0:global_ny,0:2*nxw-1,0:nbuff-1))
      allocate(pb2(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
      allocate(pb3(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    
      wc1o(:,:,:,:) = (0._DP, 0._DP)
      wc2o(:,:,:,:) = (0._DP, 0._DP)
      wc1e(:,:,:,:) = (0._DP, 0._DP)
      wc2e(:,:,:,:) = (0._DP, 0._DP)
      
      pssn(:,:,:,:) = (0._DP, 0._DP)
      dpdx(:,:,:) = 0._DP
      dpdy(:,:,:) = 0._DP
      dadx(:,:,:) = 0._DP
      dady(:,:,:) = 0._DP
    
      call pnl_pack_psi_y2zm(psi,wc1o)
      call pnl_transpose_y2zm(wc1o,wc2o)
      call pnl_pssn_unpack_y2zm(wc2o,wwdxo,wwdyo)
      
      call pnl_pack_psi_y2zm(chi,wc1e)
      call pnl_transpose_y2zm(wc1e,wc2e)
      call pnl_pssn_unpack_y2zm(wc2e,wwdxe,wwdye)
      
      call pnl_backwardfft_y2zm(wwdxo,dpdx)
      call pnl_backwardfft_y2zm(wwdyo,dpdy) 
      call pnl_backwardfft_y2zm(wwdxe,dadx)
      call pnl_backwardfft_y2zm(wwdye,dady)
    
      call pnl_pssn_realspcal_y2zm( dpdx, dpdy, dadx, dady, pb1 )
      call pnl_pack_zm2y(pb1,pb2)
      call pnl_transpose_zm2y(pb2,pb3) 
      call pnl_pssn_unpack_zm2y(pb3,pssn) 
      
      deallocate(dpdx)
      deallocate(dpdy)
      deallocate(dadx)
      deallocate(dady)
      deallocate(wc1o)
      deallocate(wc2o)
      deallocate(wc1e)
      deallocate(wc2e)
      deallocate(wwdxo)
      deallocate(wwdyo)
      deallocate(wwdxe)
      deallocate(wwdye)
      deallocate(pb1)
      deallocate(pb2)
      deallocate(pb3)      
      
  END SUBROUTINE calc_Ez_pssn_y2zm
  
!--------------------------------------
  SUBROUTINE pnl_dft_y2zm( dfdv, wkd1, dpot, pnl )
!--------------------------------------
    
    implicit none
    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pnl    
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dfdv
    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: wkd1 ! dphi/dz or dA/dt in k-space
    real(kind=DP), dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1), intent(out) :: dpot
    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
         wc1o, wc2o, wc3o, wc4o,   wc1e, wc2e, wc3e, wc4e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdvo, wwdve,    wwefo, wwefe
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wkd2, wkd3
    complex(kind=DP), allocatable, dimension(:,:,:) :: wwkd
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
    allocate(wkd2(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wkd3(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwkd(0:global_ny,0:2*nxw-1,0:nbuff-1))
    
!$OMP parallel default(none)           &
!$OMP shared(dfdv,pnl)                 &
!$OMP shared(wc1o,wc2o,wc3o,wc4o)      &
!$OMP shared(wc1e,wc2e,wc3e,wc4e)      &
!$OMP shared(wwdvo,wwdve,wwefo,wwefe)  &
!$OMP shared(wkd1,wkd2,wkd3,wwkd,dpot) &
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

      wkd2(:,:,:,:) = (0._DP, 0._DP)
      wkd3(:,:,:,:) = (0._DP, 0._DP)
      wwkd(:,:,:)   = (0._DP, 0._DP)
      dpot(:,:,:)   = 0._DP

      wwdve(:,:,:) = (0._DP, 0._DP)
      wwdvo(:,:,:) = (0._DP, 0._DP)
      wwefe(:,:,:) = (0._DP, 0._DP)
      wwefo(:,:,:) = (0._DP, 0._DP)
!$OMP end workshare

    call pnl_pack_dpot_y2zm( wkd1, wkd2 )
!$OMP barrier
    call pnl_transpose_y2zm( wkd2, wkd3 )
!$OMP barrier
    call pnl_unpack_dpot_y2zm( wkd3, wwkd ) ! fftw in x-direction
!$OMP barrier
    call pnl_backwardfft_y2zm( wwkd, dpot ) ! fftw in y-direction

    do iv = 1, 2*nv+6
      if (mod(iv,2) == 1) then ! odd
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1e,wc2e)        ! 3
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3e,wc4e)        ! 7
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1o)     ! 1
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2o,wwdvo)          ! 3
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpot,wwdve,wwefe) ! 5
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwefo,wc3o)            ! 5
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4o,pnl)       ! 7
      else                     ! even
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1o,wc2o)        ! 2
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3o,wc4o)        ! 6
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,dfdv,wc1e)     ! 2
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2e,wwdve)          ! 4
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpot,wwdvo,wwefo) ! 4
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

    deallocate( wkd2 )
    deallocate( wkd3 )
    deallocate( wwkd )

      
  END SUBROUTINE pnl_dft_y2zm
  
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
  SUBROUTINE pnl_pack_psi_y2zm( win, wout )
!--------------------------------------

!   array size is different from pnl_pack_dpot_y2zm
    
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: win
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

  END SUBROUTINE pnl_pack_psi_y2zm

  
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
  SUBROUTINE pnl_pssn_unpack_y2zm ( wc4, wwdx, wwdy )
!--------------------------------------
!     Data unpack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwdx, wwdy

    complex(kind=DP), dimension(-nx:nx) :: psi
    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    integer :: my, ibuff, iprocw, global_my

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
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
          w1(0:nx) = ui * kx(0:nx) * psi(0:nx)
          w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
          w1(2*nxw-nx:2*nxw-1) = ui * kx(-nx:-1) * psi(-nx:-1)
          call dfftw_execute_dft(plan_x_backward, w1, w2)
          wwdx(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Backward x-FFT (kx,ky)->(ky,x) %%%
          w1(0:nx) = ui * gky(global_my) * psi(0:nx)
          w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
          w1(2*nxw-nx:2*nxw-1) = ui * gky(global_my) * psi(-nx:-1)
          call dfftw_execute_dft(plan_x_backward, w1, w2)
          wwdy(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_pssn_unpack_y2zm

!--------------------------------------
  SUBROUTINE pnl_pssn_realspcal_y2zm ( dpdx, dpdy, dadx, dady, pb )
!--------------------------------------
!     Calculate Poisson brackets for pnl term calculation (y2zm)

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1) :: dpdx, dpdy, dadx, dady
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: pb

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) :: pbxy
    real(kind=DP) :: cef
    integer :: ix, iy, ibuff

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!$OMP do collapse(2) schedule(dynamic,nchunk_xb)
      do ibuff = 0, nbuff-1
        do ix = 0, 2*nxw-1

         !%%% Poisson brackets in (y,x) %%%
          do iy = 0, 2*nyw-1
            pbxy(iy) = cef * ( & ! Normalization for 2D Forward FFT
                       dpdx(iy,ix,ibuff) * dady(iy,ix,ibuff) &
                     - dpdy(iy,ix,ibuff) * dadx(iy,ix,ibuff) )
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Forward y-FFT (y,x)->(ky,x) %%%
          call dfftw_execute_dft_r2c(plan_y_forward, pbxy, w3)
          pb(0:global_ny,ix,ibuff) = w3(0:global_ny)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_pssn_realspcal_y2zm
  
!--------------------------------------
  SUBROUTINE pnl_pssn_unpack_zm2y ( wc4, pssn )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: pssn

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
          pssn(:,:,iz,im) = wc4(:,:,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE pnl_pssn_unpack_zm2y
  
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
    real(kind=DP), allocatable, dimension(:,:,:) :: Ezr
    !!real(kind=DP), allocatable, dimension(:,:,:,:) :: Ezr    
    real(kind=DP) :: cefv
    integer  ::  mx, my, iz, iv, im
    
    allocate(dfdv(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
    allocate(dAdt(-nx:nx,0:ny,-nz:nz-1))
    allocate(Ezk(-nx:nx,0:ny,-nz:nz-1,0:nm))
    !!!allocate(Ezr(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
    allocate(Ezr(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
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
      call calc_Ez_pssn_y2zm( psi, chi, Ezk )       
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
   
    !!call pnl_dft_y2x( dfdv, Ezk, Ezr, pnl )
    call pnl_dft_y2zm( dfdv, Ezk, Ezr, pnl )
    
!$OMP parallel
    do im = 0, nm
!$OMP do
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              pnl(mx,my,iz,iv,im) = dref * pnl(mx,my,iz,iv,im) * fcs(ranks) * tau(ranks) / Znum(ranks) * &
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
    deallocate(Ezk)
    deallocate(Ezr)
    deallocate(dAdt)
    deallocate(pnl)

  END SUBROUTINE pnl_sum
   
END MODULE GKV_pnl

