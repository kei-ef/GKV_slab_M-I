MODULE GKV_pnl

  use GKV_header
  use GKV_mpienv
  use GKV_fft, only: plan_x_forward, plan_x_backward, &
                     plan_y_forward, plan_y_backward
  use GKV_clock, only: clock_sta, clock_end
  use GKV_micouple, only: rho2R0

  implicit none

  private
  integer, parameter :: nbuff = ((2*nz)*(nm+1)-1)/nprocw + 1
  integer, save :: nchunk_zm  = 1, nchunk_yb  = 1, nchunk_xb = 1
  integer, save :: nchunk_yzm = 1, nchunk_xzm = 1
  
  public PNL_es_part


CONTAINS
!--------------------------------------
  SUBROUTINE PNL_es_part( ff, psi, ef )
!--------------------------------------

    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: wdv1
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: pnl
    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
         wc1o, wc2o, wc3o, wc4o,   wc1e, wc2e, wc3e, wc4e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdvo, wwdve,    wwefo, wwefe
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: wdz1, wdz2, wdz3
    complex(kind=DP), allocatable, dimension(:,:,:) :: wwdz
    real(kind=DP), allocatable, dimension(:,:,:) :: dpdz
    real(kind=DP) :: delta
    integer :: mx, my, iz, iv, im

    allocate(dpdz(0:2*nyw-1,0:2*nxw-1,0:nbuff-1))
    allocate(pnl(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
!
    allocate(wdv1(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm))
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

    if ( trim(z_bound) == "mi_couple"  ) then
      delta = rho2R0
    else   
      delta = 2.d-4
    end if    
    
!... make df/dvz and dpsi/dz

    do im = 0, nm
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do iv = 1, 2*nv
              wdv1(mx,my,iz,iv,im) =        ( - ff(mx,my,iz,iv+2,im)      &
                                      + 8._DP * ff(mx,my,iz,iv+1,im)      &
                                      - 8._DP * ff(mx,my,iz,iv-1,im)      &
                                              + ff(mx,my,iz,iv-2,im) )
            end do

              wdz1(mx,my,iz,im)    =       ( - psi(mx,my,iz+2,im)    &
                                     + 8._DP * psi(mx,my,iz+1,im)    &
                                     - 8._DP * psi(mx,my,iz-1,im)    &
                                             + psi(mx,my,iz-2,im) ) 
          end do
        end do
      end do
    end do
   
!$OMP parallel default(none)                          &
!$OMP shared(wdv1,wdz1,dpdz,pnl)       &
!$OMP shared(wc1o,wc2o,wc3o,wc4o,wc1e,wc2e,wc3e,wc4e) &
!$OMP shared(wdz2,wdz3,wwdz)     &
!$OMP shared(wwdvo,wwdve,wwefo,wwefe)     &
!$OMP private(iv)

!$OMP workshare
      pnl(:,:,:,:,:) = (0._DP, 0._DP)
      wc1o(:,:,:,:)  = (0._DP, 0._DP)
      wc3o(:,:,:,:)  = (0._DP, 0._DP)
      wc1e(:,:,:,:)  = (0._DP, 0._DP)
      wc3e(:,:,:,:)  = (0._DP, 0._DP)

      wwdz(:,:,:) = (0._DP, 0._DP)
      dpdz(:,:,:) = 0._DP
!$OMP end workshare

    call pnl_pack_dpdz( wdz1, wdz2 )
!$OMP barrier
    call pnl_transpose_y2zm( wdz2, wdz3 )
!$OMP barrier    
    call pnl_unpack_dpdz( wdz3, wwdz ) ! fftw in x-direction
!$OMP barrier    
    call pnl_backwardfft( wwdz, dpdz ) ! fftw in y-direction
    
    do iv = 1, 2*nv+6
      if (mod(iv,2) == 1) then ! odd
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1e,wc2e)        ! 3
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3e,wc4e)        ! 7
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,wdv1,wc1o)     ! 1
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2o,wwdvo)          ! 3
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpdz,wwdve,wwefe) ! 5
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwdvo,wc3o)            ! 5
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4o,pnl)       ! 7
      else                     ! even
        if (1+1<=iv .and. iv<=2*nv+1) call pnl_transpose_y2zm(wc1o,wc2o)        ! 2
        if (1+5<=iv .and. iv<=2*nv+5) call pnl_transpose_zm2y(wc3o,wc4o)        ! 6
        if (1  <=iv .and. iv<=2*nv  ) call pnl_pack_dfdv_y2zm(iv,wdv1,wc1e)     ! 2
        if (1+2<=iv .and. iv<=2*nv+2) call pnl_unpack_y2zm(wc2e,wwdve)          ! 4
        if (1+3<=iv .and. iv<=2*nv+3) call pnl_realspcal_y2zm(dpdz,wwdvo,wwefo) ! 4
        if (1+4<=iv .and. iv<=2*nv+4) call pnl_pack_zm2y(wwdve,wc3e)            ! 6
        if (1+6<=iv .and. iv<=2*nv+6) call pnl_unpack_zm2y(iv-6,wc4e,pnl)       ! 8
      end if
!$OMP barrier
    end do
!$OMP end parallel

    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              ef(mx,my,iz,iv,im) = ef(mx,my,iz,iv,im) - delta * pnl(mx,my,iz,iv,im)    
            end do
          end do
        end do
      end do
    end do
    
    deallocate( dpdz )
    deallocate( pnl )

    deallocate( wdv1 )    
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

  END SUBROUTINE PNL_es_part

!--------------------------------------
  SUBROUTINE pnl_pack_dpdz( wdz1, wdz2 )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: wdz1
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wdz2

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
              wdz2(mx,my,ibuff,iprocw) = wdz1(mx,my,iz,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1410)
!$OMP end master
  
  END SUBROUTINE pnl_pack_dpdz


!--------------------------------------
  SUBROUTINE pnl_transpose_y2zm ( wc4in, wc4out )
!--------------------------------------
!     Data transpose for E x B term calculation (y2zm)

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

!!--------------------------------------
!  SUBROUTINE pnl_transpose( wwin, wwout )
!!--------------------------------------
!!     basically the same as exb_transpose_y2x 
!
!    complex(kind=DP), intent(in), &
!      dimension(0:ny,0:nxw_size,-nz:nz-1,0:nm,0:nprocw-1) :: wwin
!    complex(kind=DP), intent(out), &
!      dimension(0:ny,0:nxw_size,-nz:nz-1,0:nm,0:nprocw-1) :: wwout

!!$OMP master
!                                           call clock_sta(1420)
!      call MPI_Alltoall( wwin,                              &
!                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
!                         MPI_DOUBLE_COMPLEX,                &
!                         wwout,                             &
!                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
!                         MPI_DOUBLE_COMPLEX,                &
!                         fft_comm_world,                    &
!                         ierr_mpi )
!                                           call clock_end(1420)
!!$OMP end master
!
!  END SUBROUTINE pnl_transpose

!--------------------------------------
  SUBROUTINE pnl_unpack_dpdz( wc4, wwdz )
!--------------------------------------
!     Data unpack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwdz

    complex(kind=DP), dimension(-nx:nx) :: dpdz
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
          dpdz(-nx:nx) = wc4(-nx:nx,my,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Backward x-FFT (kx,ky)->(ky,x) %%%
          w1(0:nx) = dpdz(0:nx)
          w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
          w1(2*nxw-nx:2*nxw-1) = dpdz(-nx:-1)
          call dfftw_execute_dft(plan_x_backward, w1, w2)
          wwdz(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_unpack_dpdz

  

!--------------------------------------
  SUBROUTINE pnl_backwardfft ( wwdz, dpdz )
!--------------------------------------
!  
    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,-nz:nz-1,0:nm,0:nprocw-1) :: wwdz
    real(kind=DP), intent(inout), &
      dimension(0:2*nyw-1,0:nxw_size,-nz:nz-1,0:nm) :: dpdz

    complex(kind=DP), dimension(0:nyw) :: w3
    integer :: ix, my, iz, im, irank, ist_y_g_rank, iend_y_g_rank

!$OMP master
                                           call clock_sta(1430)
!$OMP end master
!$OMP do collapse(3) schedule(dynamic,nchunk_xzm)
      do im = 0, nm
        do iz = -nz, nz-1
          do ix = ist_xw, iend_xw

           !%%% UNPACK: (x,ky*,z*,m*)->(ky,x*,z*,m*) %%%
            do irank = 0, nprocw-1
              ist_y_g_rank  = (ny+1)*irank
              iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
              do my = ist_y_g_rank, iend_y_g_rank
                w3(my) = wwdz(my-ist_y_g_rank,ix,iz,im,irank)
              end do
            end do
            w3(global_ny+1:nyw) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !%%% Backward y-FFT (ky,x)->(y,x) %%%
            call dfftw_execute_dft_c2r(plan_y_backward, w3, dpdz(:,ix,iz,im))
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          end do
                 
        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE pnl_backwardfft

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
  SUBROUTINE pnl_realspcal_y2zm ( dpdz, wwdv, wwef )
!--------------------------------------

    !!integer, intent(in) :: iv
    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:2*nxw-1,0:nbuff-1) :: dpdz
    complex(kind=DP), intent(in), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwdv
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: wwef

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) :: dfdv, pnle
    real(kind=DP) :: cef, cs1
    integer :: ix, iy, ibuff

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)
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
            pnle(iy) = cef * & ! Normalization for 2D Forward FFT
                       cs1 * dpdz(iy,ix,ibuff) * dfdv(iy) 
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% Forward y-FFT (y,x)->(ky,x) %%%
          call dfftw_execute_dft_r2c(plan_y_forward, pnle, w3)
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
  SUBROUTINE pnl_unpack_zm2y ( iv, wc4, ef )
!--------------------------------------

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    !complex(kind=DP), intent(inout), &
    !  dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: ef
    !NOTE: A noncontiguous subarray as an argument of a subroutine induces a memory copy.
    !      When the subroutine is called in a OpenMP parallel region,
    !      the copied subarray may be treated as a thread-private variable.

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
          ef(:,:,iz,iv,im) = wc4(:,:,ibuff,iprocw)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait
!$OMP master
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE pnl_unpack_zm2y

  
END MODULE GKV_pnl
