MODULE RMHDS_pssn
!-------------------------------------------------------------------------------
!
!    Calculation of Poisson bracket
!
!      created by Maeyama
!      modified by THW for aurora spectral code
!
!... 2023-03-24
!      modified by K. Fujita based on gkvp_exb.f90, 
!                 created by Profs. THW and Maeyama
!
!      FFTW:   1D FFT
!      OpenMP: kx,ky-loop parallelization
!
!-------------------------------------------------------------------------------

!  use RMHDS_header
  use GKV_header
  use GKV_mpienv
  use GKV_clock, only: clock_sta, clock_end
  implicit none
  include "fftw3.f"
  private

  integer(kind=DP), save :: plan_backward_x, plan_forward_x
  integer(kind=DP), save :: plan_backward_y, plan_forward_y
  complex(kind=DP), dimension(0:nny/2,0:nnx-1) :: wwkk
  real(kind=DP), save :: ipb_maxvx, ipb_maxvy
  integer, save :: iflg
  data iflg / 0 /

  public  pssn_brackets, pssn_divgrad, ipb_maxvx, ipb_maxvy

 CONTAINS

!--------------------------------------
SUBROUTINE pssn_brackets( fk, gk, pb, nmw )
!--------------------------------------
!
!    Calculate poisson bracket { fk, gk }
!
  use GKV_header
  use GKV_mpienv

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(in) :: gk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(out) :: pb
  integer, intent(in) :: nmw

  real(kind=DP), dimension(0:2*nyw-1,0:nxw_size) :: dfdx, dfdy, dgdx, dgdy
!!  real(kind=DP), dimension(:,:), allocatable :: dfdx, dfdy, dgdx, dgdy
  complex(kind=DP), dimension(:,:,:), allocatable :: &
               fkyx_dx, fkyx_dy, fkyx_dx_T, fkyx_dy_T, &
               gkyx_dx, gkyx_dy, gkyx_dx_T, gkyx_dy_T, &
               pbyx, pbyx_T

  integer :: im
  integer :: ix, iy, mx, my

!*****
!$OMP master
    call clock_sta(1610)
!$OMP end master
!*****

!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if

    allocate(fkyx_dx(0:ny,0:nxw_size,0:nprocw-1))
    allocate(fkyx_dy(0:ny,0:nxw_size,0:nprocw-1))
    allocate(fkyx_dx_T(0:ny,0:nxw_size,0:nprocw-1))
    allocate(fkyx_dy_T(0:ny,0:nxw_size,0:nprocw-1))

    allocate(gkyx_dx(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dy(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dx_T(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dy_T(0:ny,0:nxw_size,0:nprocw-1))

    allocate(pbyx(0:ny,0:nxw_size,0:nprocw-1))
    allocate(pbyx_T(0:ny,0:nxw_size,0:nprocw-1))

!!!!$OMP PARALLEL DO default(shared) private(im)
    do im = 0, nmw-1

    pb(:,:,im)  = (0._DP, 0._DP)
    pbyx(:,:,:) = (0._DP, 0._DP)

    fkyx_dx(:,:,:) = (0._DP, 0._DP)
    fkyx_dy(:,:,:) = (0._DP, 0._DP)
    gkyx_dx(:,:,:) = (0._DP, 0._DP)
    gkyx_dy(:,:,:) = (0._DP, 0._DP)

    dfdx(:,:) = 0._DP
    dfdy(:,:) = 0._DP
    dgdx(:,:) = 0._DP
    dgdy(:,:) = 0._DP

!... f
    call pssn_pack_y2x( fk(:,:,im), fkyx_dx, fkyx_dy ) ! derivaties in (ky,x)-space
    call pssn_transpose_y2x( fkyx_dx, fkyx_dx_T ) ! transpose
    call pssn_transpose_y2x( fkyx_dy, fkyx_dy_T )
    call pssn_backwardfft_y2x( fkyx_dx_T, fkyx_dy_T, dfdx, dfdy ) ! derivaties in (y,x)-space

!... g
    call pssn_pack_y2x( gk(:,:,im), gkyx_dx, gkyx_dy ) ! derivaties in (ky,x)-space
    call pssn_transpose_y2x( gkyx_dx, gkyx_dx_T ) ! transpose
    call pssn_transpose_y2x( gkyx_dy, gkyx_dy_T )
    call pssn_backwardfft_y2x( gkyx_dx_T, gkyx_dy_T, dgdx, dgdy ) ! derivaties in (y,x)-space

!... { f, g }
    call pssn_realspcal_y2x( dfdx, dfdy, dgdx, dgdy, pbyx ) ! pb in (y,x)-space
    call pssn_transpose_x2y( pbyx, pbyx_T )
    call pssn_unpack_x2y( pbyx_T, pb(:,:,im) ) ! pb in (kx,ky)-space

    call ipb_estimate_max(dfdx,dfdy,dgdx,dgdy)

    end do
!!!!!$OMP END PARALLEL DO

    deallocate(fkyx_dx)
    deallocate(fkyx_dy)
    deallocate(fkyx_dx_T)
    deallocate(fkyx_dy_T)

    deallocate(gkyx_dx)
    deallocate(gkyx_dy)
    deallocate(gkyx_dx_T)
    deallocate(gkyx_dy_T)

    deallocate(pbyx)
    deallocate(pbyx_T)

!*****
!$OMP master
    call clock_end(1610)
!$OMP end master
!*****

    return

END SUBROUTINE pssn_brackets

!-----------------------------------------------
   SUBROUTINE pssn_divgrad( fk, gk, dg, nmw )
!-----------------------------------------------
! 
!  calculate   
!     div( fk * grad(gk) ) 
!
  use GKV_header
  use GKV_mpienv

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(in) :: gk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(out) :: dg
  integer, intent(in) :: nmw

  real(kind=DP), dimension(:,:), allocatable :: fyx
  complex(kind=DP), dimension(:,:,:), allocatable :: &
               fkyx, fkyx_T

  real(kind=DP), dimension(:,:), allocatable :: fdgdx, fdgdy, dgdx, dgdy
  complex(kind=DP), dimension(:,:), allocatable :: ikxfg, ikyfg
  complex(kind=DP), dimension(:,:,:), allocatable :: &
               gkyx_dx, gkyx_dy, gkyx_dx_T, gkyx_dy_T, &
               ikxh, ikxh_T, ikyh, ikyh_T

  integer :: im, ix, iy, mx, my

!*****
!$OMP master
    call clock_sta(1620)
!$OMP end master
!*****

!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if

    allocate(fkyx  (0:ny,0:nxw_size,0:nprocw-1))
    allocate(fkyx_T(0:ny,0:nxw_size,0:nprocw-1))

    allocate(fyx  (0:2*nyw-1,0:nxw_size))
    allocate(fdgdx(0:2*nyw-1,0:nxw_size))
    allocate(fdgdy(0:2*nyw-1,0:nxw_size))
    allocate(dgdx (0:2*nyw-1,0:nxw_size))
    allocate(dgdy (0:2*nyw-1,0:nxw_size))

    allocate(gkyx_dx(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dy(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dx_T(0:ny,0:nxw_size,0:nprocw-1))
    allocate(gkyx_dy_T(0:ny,0:nxw_size,0:nprocw-1))

    allocate(ikxh(0:ny,0:nxw_size,0:nprocw-1))
    allocate(ikxh_T(0:ny,0:nxw_size,0:nprocw-1))
    allocate(ikyh(0:ny,0:nxw_size,0:nprocw-1))
    allocate(ikyh_T(0:ny,0:nxw_size,0:nprocw-1))

    allocate(ikxfg(-nx:nx,0:ny), ikyfg(-nx:nx,0:ny))

!!!!!$OMP DO default(shared) private(im)
    do im = 0, nmw-1

      dg(:,:,im)  = (0._DP, 0._DP)

      fkyx(:,:,:) = (0._DP, 0._DP)

      gkyx_dx(:,:,:) = (0._DP, 0._DP)
      gkyx_dy(:,:,:) = (0._DP, 0._DP)
      dgdx(:,:) = 0._DP
      dgdy(:,:) = 0._DP

      fdgdx(:,:) = 0._DP
      fdgdy(:,:) = 0._DP

      ikxh(:,:,:) = (0._DP, 0._DP)
      ikyh(:,:,:) = (0._DP, 0._DP)

! f
      call fft_pack_y2x( fk(:,:,im), fkyx )
      call pssn_transpose_y2x( fkyx, fkyx_T ) ! transpose
      call fft_backwardfft_y2x( fkyx_T, fyx ) ! f in (y,x)-space

! div g
      call pssn_pack_y2x( gk(:,:,im), gkyx_dx, gkyx_dy ) ! derivaties in (ky,x)-space
      call pssn_transpose_y2x( gkyx_dx, gkyx_dx_T ) ! transpose
      call pssn_transpose_y2x( gkyx_dy, gkyx_dy_T )
      call pssn_backwardfft_y2x( gkyx_dx_T, gkyx_dy_T, dgdx, dgdy ) ! derivaties in (y,x)-space

      do ix = 0, nxw_size
        do iy = 0,  2*nyw-1
          fdgdx(iy,ix) = fyx(iy,ix) * dgdx(iy,ix)
          fdgdy(iy,ix) = fyx(iy,ix) * dgdy(iy,ix)
        end do
      end do

      call fft_realspcal_y2x( fdgdx, ikxh ) ! (y,x) -> (kx,ky)
      call fft_realspcal_y2x( fdgdy, ikyh ) ! (y,x) -> (kx,ky)
 
      call pssn_transpose_x2y( ikxh, ikxh_T )
      call pssn_transpose_x2y( ikyh, ikyh_T )

      call pssn_unpack_x2y( ikxh_T, ikxfg ) ! (kx,ky)
      call pssn_unpack_x2y( ikyh_T, ikyfg ) ! (kx,ky)

      do my = 0, ny
        do mx = -nx, nx
          dg(mx,my,im) = ui * ( kx(mx) * ikxfg(mx,my) + ky(my) * ikyfg(mx,my) )
        end do
      end do

    end do ! im loop
!!!!$OMP END DO


    deallocate(fkyx  )
    deallocate(fkyx_T)
    deallocate(fyx  )
    deallocate(fdgdx)
    deallocate(fdgdy)
    deallocate(dgdx )
    deallocate(dgdy )

    deallocate(gkyx_dx)
    deallocate(gkyx_dy)
    deallocate(gkyx_dx_T)
    deallocate(gkyx_dy_T)

    deallocate(ikxh)
    deallocate(ikxh_T)
    deallocate(ikyh)
    deallocate(ikyh_T)
    deallocate(ikxfg)
    deallocate(ikyfg)

!*****
!$OMP master
    call clock_end(1620)
!$OMP end master
!*****

    return
  END SUBROUTINE pssn_divgrad

!-----------------------------------------------
  SUBROUTINE pssn_pack_y2x( wk, dwdx, dwdy )
!-----------------------------------------------
!  Data pack for 2D nonlinear term calculation (y2x)

  use GKV_header
  implicit none
  complex(kind=DP), intent(in), dimension(-nx:nx,0:ny) :: wk
  complex(kind=DP), intent(inout),  dimension(0:ny,0:nxw_size,0:nprocw-1) :: dwdx, dwdy

!... local variables
  complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
  integer :: ix, my, irank, ist_xw_g_rank, iend_xw_g_rank

    do my = ist_y, iend_y

      !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
      w1(0:nx)             = ui * kx(0:nx) * wk(0:nx,my)
      w1(nx+1:2*nxw-nx-1)  = (0._DP, 0._DP) ! FFTW may destroy input array!
      w1(2*nxw-nx:2*nxw-1) = ui * kx(-nx:-1) * wk(-nx:-1,my)
      call dfftw_execute_dft(plan_backward_x, w1, w2) 

      !%%% PACK: (x,ky*)->(ky,x*) %%%
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do ix = ist_xw_g_rank, iend_xw_g_rank
           dwdx(my,ix-ist_xw_g_rank,irank) = w2(ix)
        enddo
      enddo

      !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
      w1(0:nx) = ui * ky(my) * wk(0:nx,my)
      w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
      w1(2*nxw-nx:2*nxw-1) = ui * ky(my) * wk(-nx:-1,my)
      call dfftw_execute_dft(plan_backward_x, w1, w2)
 
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% PACK: (x,ky*)->(ky,x*) %%%
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do ix = ist_xw_g_rank, iend_xw_g_rank
          dwdy(my,ix-ist_xw_g_rank,irank) = w2(ix)
        enddo
      enddo
             
    enddo

    return
  END SUBROUTINE pssn_pack_y2x

!--------------------------------------
  SUBROUTINE pssn_transpose_y2x ( wwin, wwout )
!--------------------------------------
!   Data transpose for 2d nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwout

    wwout(:,:,:) = ( 0._DP, 0._DP )

      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )

  END SUBROUTINE pssn_transpose_y2x

!--------------------------------------
  SUBROUTINE pssn_backwardfft_y2x ( wwdx, wwdy, dpdx, dpdy )
!--------------------------------------
!   Backward FFT of 2d functions for nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwdx, wwdy
    real(kind=DP), intent(inout), &
      dimension(0:2*nyw-1,0:nxw_size) :: dpdx, dpdy

    complex(kind=DP), dimension(0:nyw) :: w3
    integer :: ix, my, irank, ist_y_g_rank, iend_y_g_rank

    do ix = ist_xw, iend_xw

       !%%% UNPACK: (x,ky*)->(ky,x*) %%%
       do irank = 0, nprocw-1
          ist_y_g_rank  = (ny+1)*irank
          iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
          do my = ist_y_g_rank, iend_y_g_rank
             w3(my) = wwdx(my-ist_y_g_rank,ix,irank)
          end do
       end do
       w3(global_ny+1:nyw) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% Backward y-FFT (ky,x)->(y,x) %%%
       call dfftw_execute_dft_c2r(plan_backward_y, w3, dpdx(:,ix))
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       !%%% UNPACK: (x,ky*)->(ky,x*) %%%
       do irank = 0, nprocw-1
          ist_y_g_rank  = (ny+1)*irank
          iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
          do my = ist_y_g_rank, iend_y_g_rank
             w3(my) = wwdy(my-ist_y_g_rank,ix,irank)
          end do
       end do
       w3(global_ny+1:nyw) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% Backward y-FFT (ky,x)->(y,x) %%%
       call dfftw_execute_dft_c2r(plan_backward_y, w3, dpdy(:,ix))
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end do

  END SUBROUTINE pssn_backwardfft_y2x


!--------------------------------------
  SUBROUTINE fft_backwardfft_y2x ( fin, fout )
!--------------------------------------
!   Backward FFT of 2d functions for nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: fin
    real(kind=DP), intent(out), &
      dimension(0:2*nyw-1,0:nxw_size) :: fout

    complex(kind=DP), dimension(0:nyw) :: ftmp3
    integer :: ix, my, irank, ist_y_g_rank, iend_y_g_rank

    do ix = ist_xw, iend_xw

       !%%% UNPACK: (x,ky*)->(ky,x*) %%%
       do irank = 0, nprocw-1
          ist_y_g_rank  = (ny+1)*irank
          iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
          do my = ist_y_g_rank, iend_y_g_rank
             ftmp3(my) = fin(my-ist_y_g_rank,ix,irank)
          end do
       end do
       ftmp3(global_ny+1:nyw) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% Backward y-FFT (ky,x)->(y,x) %%%
       call dfftw_execute_dft_c2r(plan_backward_y, ftmp3, fout(:,ix))
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end do
!!!!$OMP end do

  END SUBROUTINE fft_backwardfft_y2x

!--------------------------------------
  SUBROUTINE pssn_realspcal_y2x ( dfdx, dfdy, dgdx, dgdy, wwef )
!--------------------------------------
!   Calculate Poisson brackets for 2d nonlinear term calculation (y2x)

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:nxw_size) :: dfdx, dfdy, dgdx, dgdy
    complex(kind=DP), intent(inout), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwef

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) ::  pbxy
    real(kind=DP) :: cef
    integer :: my, ix, iy, irank, ist_y_g_rank, iend_y_g_rank

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

      do ix = ist_xw, iend_xw
 
         !%%% Poisson brackets in (y,x) %%%
         do iy = 0, 2*nyw-1
            pbxy(iy) = cef * ( dfdx(iy,ix) * dgdy(iy,ix) - dfdy(iy,ix) * dgdx(iy,ix) )
         end do
 
         !%%% Forward y-FFT (y,x)->(ky,x) %%%
         call dfftw_execute_dft_r2c(plan_forward_y, pbxy, w3)

         !%%% PACK: (ky,x*)->(x,ky*) %%%
         do irank = 0, nprocw-1
            ist_y_g_rank  = (ny+1)*irank
            iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
            do my = ist_y_g_rank, iend_y_g_rank
               wwef(my-ist_y_g_rank,ix,irank) = w3(my)
            end do
         end do

      end do

    return
  END SUBROUTINE pssn_realspcal_y2x

!--------------------------------------
  SUBROUTINE fft_realspcal_y2x ( fin, fout )
!--------------------------------------
!   Calculate Poisson brackets for 2d nonlinear term calculation (y2x)
!   (y,x) -> (kx,ky)

    real(kind=DP), intent(in), dimension(0:2*nyw-1,0:nxw_size) :: fin
    complex(kind=DP), intent(out), dimension(0:ny,0:nxw_size,0:nprocw-1) :: fout

    complex(kind=DP), dimension(0:nyw) :: ftmp2
    real(kind=DP), dimension(0:2*nyw-1) ::  ftmp1
    real(kind=DP) :: cef
    integer :: my, ix, iy, irank, ist_y_g_rank, iend_y_g_rank

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

      do ix = ist_xw, iend_xw
 
        !%%% Forward y-FFT (y,x)->(ky,x) %%%
        ftmp1(:) = cef * fin(:,ix) 
        call dfftw_execute_dft_r2c(plan_forward_y, ftmp1, ftmp2)

         !%%% PACK: (ky,x*)->(x,ky*) %%%
         do irank = 0, nprocw-1
            ist_y_g_rank  = (ny+1)*irank
            iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
            do my = ist_y_g_rank, iend_y_g_rank
               fout(my-ist_y_g_rank,ix,irank) = ftmp2(my)
            end do
         end do
      end do

    return
  END SUBROUTINE fft_realspcal_y2x

!--------------------------------------
  SUBROUTINE pssn_transpose_x2y ( wwin, wwout )
!--------------------------------------
!  Data transpose for 2d nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwout

      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )

  END SUBROUTINE pssn_transpose_x2y


!--------------------------------------
  SUBROUTINE pssn_unpack_x2y ( wwin, wwout )
!--------------------------------------
!   Data unpack for 2d nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny) :: wwout

    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    integer :: ix, my, irank, ist_xw_g_rank, iend_xw_g_rank

    do my = ist_y, iend_y

       !%%% UNPACK: (ky,x*)->(x,ky*) %%%
       do irank = 0, nprocw-1
          ist_xw_g_rank  = (nxw_size+1)*irank
          iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
          do ix = ist_xw_g_rank, iend_xw_g_rank
             w2(ix) = wwin(my,ix-ist_xw_g_rank,irank)
          enddo
       enddo
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% Forward x-FFT (x,ky)->(kx,ky) %%%
       call dfftw_execute_dft(plan_forward_x, w2, w1)
       wwout(0:nx,my)   = w1(0:nx)
       wwout(-nx:-1,my) = w1(2*nxw-nx:2*nxw-1)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end do

    return

  END SUBROUTINE pssn_unpack_x2y

!--------------------------
  SUBROUTINE fft_pack_y2x( fin, fout )
!-------------------------

  use GKV_header
  implicit none
  complex(kind=DP), intent(in), dimension(-nx:nx,0:ny) :: fin
  complex(kind=DP), dimension(0:ny,0:nxw_size,0:nprocw-1), intent(out) :: fout

  complex(kind=DP), dimension(0:2*nxw-1) :: ftmp1, ftmp2
  integer :: ix, my, irank, ist_xw_g_rank, iend_xw_g_rank

    do my = ist_y, iend_y

      !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
      ftmp1(0:nx)             = fin(0:nx,my)
      ftmp1(nx+1:2*nxw-nx-1)  = (0._DP, 0._DP) ! FFTW may destroy input array!
      ftmp1(2*nxw-nx:2*nxw-1) = fin(-nx:-1,my)
      call dfftw_execute_dft(plan_backward_x, ftmp1, ftmp2) 

      !%%% PACK: (x,ky*)->(ky,x*) %%%
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do ix = ist_xw_g_rank, iend_xw_g_rank
           fout(my,ix-ist_xw_g_rank,irank) = ftmp2(ix)
        enddo
      enddo
    end do

    return
  END SUBROUTINE fft_pack_y2x


!--------------------------------------
  SUBROUTINE exb_fft_pre
!--------------------------------------
!  Initialization of FFT
  complex(kind=DP), dimension(0:2*nxw-1) :: wkx1, wkx2
  complex(kind=DP), dimension(0:nyw)  :: wky1
  real(kind=DP), dimension(0:2*nyw-1) :: wky2

    call dfftw_plan_dft_1d(plan_backward_x,  &
                           (2*nxw),               &
                           wkx1, wkx2,       &
                           FFTW_BACKWARD,    &
                           FFTW_MEASURE)
    call dfftw_plan_dft_c2r_1d(plan_backward_y,  &
                               (2*nyw),               &
                               wky1, wky2,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_r2c_1d(plan_forward_y,   &
                               (2*nyw),               &
                               wky2, wky1,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_forward_x,   &
                           (2*nxw),               &
                           wkx2, wkx1,       &
                           FFTW_FORWARD,     &
                           FFTW_MEASURE)


  END SUBROUTINE exb_fft_pre
  
!--------------------------------------
  SUBROUTINE ipb_estimate_max(dfdx,dfdy,dgdx,dgdy)
!--------------------------------------

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:nxw_size) :: dfdx, dfdy, dgdx, dgdy

    real(kind=DP) :: ipbv_nl
    integer :: ix, iy

    ipb_maxvx = eps
    ipb_maxvy = eps

    do ix = ist_xw, iend_xw
      do iy = 0, 2*nyw-1
        ipbv_nl = abs ( dfdy(iy,ix) )
        if (ipb_maxvx < ipbv_nl) ipb_maxvx = ipbv_nl
      end do
    end do

    do ix = ist_xw, iend_xw
      do iy = 0, 2*nyw-1
        ipbv_nl = abs ( dfdx(iy,ix) )
        if (ipb_maxvy < ipbv_nl) ipb_maxvy = ipbv_nl
      end do
    end do

    do ix = ist_xw, iend_xw
      do iy = 0, 2*nyw-1
        ipbv_nl = abs ( dgdy(iy,ix) )
        if (ipb_maxvx < ipbv_nl) ipb_maxvx = ipbv_nl
      end do
    end do

    do ix = ist_xw, iend_xw
      do iy = 0, 2*nyw-1
        ipbv_nl = abs ( dgdx(iy,ix) )
        if (ipb_maxvy < ipbv_nl) ipb_maxvy = ipbv_nl
      end do
    end do

    return
  END SUBROUTINE ipb_estimate_max

!--------------------------------------
SUBROUTINE nonlinear_sqr( fkxky, fsq, nmw )
!--------------------------------------
!
!    Calculate square of fin 
!
  use GKV_header
  use GKV_mpienv

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(in)  :: fkxky
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nmw-1), intent(out) :: fsq
  integer, intent(in) :: nmw

  real(kind=DP), dimension(0:2*nyw-1,0:nxw_size) :: fyx
  complex(kind=DP), dimension(:,:,:), allocatable :: & 
               fkyx, fkyx_T, fkyxsq, fkyxsq_T
  integer :: im
  integer :: ix, iy, mx, my
  
!%%% Initialize FFTW %%%
  if (iflg == 0) then
    iflg = 1
    call exb_fft_pre
  end if
 
  allocate(fkyx(0:ny,0:nxw_size,0:nprocw-1))
  allocate(fkyx_T(0:ny,0:nxw_size,0:nprocw-1))
  allocate(fkyxsq(0:ny,0:nxw_size,0:nprocw-1))
  allocate(fkyxsq_T(0:ny,0:nxw_size,0:nprocw-1))

  do im = 0, nmw-1

    fsq(:,:,im)    = (0._DP, 0._DP)
    fkyx(:,:,:)    = (0._DP, 0._DP)
    fkyxsq(:,:,:)  = (0._DP, 0._DP)
    fyx (:,:)      = 0._DP
 
!... transform fkxky to fxy
    call sqr_pack_y2x( fkxky(:,:,im), fkyx ) ! f in (ky,x)-space
    call pssn_transpose_y2x( fkyx, fkyx_T )  ! transpose
    call sqr_backwardfft_y2x( fkyx_T, fyx )  ! f in (y,x)-space

    call sqr_realspcal_y2x ( fyx, fkyxsq )     ! get f**2 in k-space

    call pssn_transpose_x2y( fkyxsq, fkyxsq_T )  ! f in 
    call pssn_unpack_x2y( fkyxsq_T, fsq(:,:,im) ) ! b in (kx,ky)-space

  end do

END SUBROUTINE nonlinear_sqr
    
!-----------------------------------------------
  SUBROUTINE sqr_pack_y2x( wk, wkyx )
!-----------------------------------------------
!  Data pack for 2D nonlinear term calculation (y2x)

  use GKV_header
  implicit none
  complex(kind=DP), intent(in), dimension(-nx:nx,0:ny) :: wk
  complex(kind=DP), intent(inout),  dimension(0:ny,0:nxw_size,0:nprocw-1) :: wkyx

!... local variables
  complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
  integer :: ix, my, irank, ist_xw_g_rank, iend_xw_g_rank

    do my = ist_y, iend_y

      !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
      w1(0:nx)             = wk(0:nx,my)
      w1(nx+1:2*nxw-nx-1)  = (0._DP, 0._DP) ! FFTW may destroy input array!
      w1(2*nxw-nx:2*nxw-1) = wk(-nx:-1,my)
      call dfftw_execute_dft(plan_backward_x, w1, w2) 

      !%%% PACK: (x,ky*)->(ky,x*) %%%
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do ix = ist_xw_g_rank, iend_xw_g_rank
           wkyx(my,ix-ist_xw_g_rank,irank) = w2(ix)
        enddo
      enddo
             
    enddo

    return
  END SUBROUTINE sqr_pack_y2x

!--------------------------------------
  SUBROUTINE sqr_backwardfft_y2x ( wkyx, wyx )
!--------------------------------------
!   Backward FFT of 2d functions for nonlinear term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: wkyx
    real(kind=DP), intent(inout), &
      dimension(0:2*nyw-1,0:nxw_size) :: wyx

    complex(kind=DP), dimension(0:nyw) :: w3
    integer :: ix, my, irank, ist_y_g_rank, iend_y_g_rank

    do ix = ist_xw, iend_xw

       !%%% UNPACK: (x,ky*)->(ky,x*) %%%
       do irank = 0, nprocw-1
          ist_y_g_rank  = (ny+1)*irank
          iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
          do my = ist_y_g_rank, iend_y_g_rank
             w3(my) = wkyx(my-ist_y_g_rank,ix,irank)
          end do
       end do
       w3(global_ny+1:nyw) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% Backward y-FFT (ky,x)->(y,x) %%%
       call dfftw_execute_dft_c2r(plan_backward_y, w3, wyx(:,ix))
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end do

  END SUBROUTINE sqr_backwardfft_y2x

!--------------------------------------
  SUBROUTINE sqr_realspcal_y2x ( fyx, fkyxsq )
!--------------------------------------
!   Calculate Poisson brackets for 2d nonlinear term calculation (y2x)

    real(kind=DP), intent(in), &
      dimension(0:2*nyw-1,0:nxw_size) :: fyx
    complex(kind=DP), intent(inout), &
      dimension(0:ny,0:nxw_size,0:nprocw-1) :: fkyxsq

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) ::  fxsq
    real(kind=DP) :: cef
    integer :: my, ix, iy, irank, ist_y_g_rank, iend_y_g_rank

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

!!!!$OMP do default(shared) 
!!!!$OMP private( my, pbxy, w3 )
!!!!$OMP private( irank, ist_y_g_rank, iend_y_g_rank, ix, iy)
      do ix = ist_xw, iend_xw
 
         !%%% Poisson brackets in (y,x) %%%
         do iy = 0, 2*nyw-1
            fxsq(iy) = cef * fyx(iy,ix) ** 2
         end do
 
         !%%% Forward y-FFT (y,x)->(ky,x) %%%
         call dfftw_execute_dft_r2c(plan_forward_y, fxsq, w3)

         !%%% PACK: (ky,x*)->(x,ky*) %%%
         do irank = 0, nprocw-1
            ist_y_g_rank  = (ny+1)*irank
            iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
            do my = ist_y_g_rank, iend_y_g_rank
               fkyxsq(my-ist_y_g_rank,ix,irank) = w3(my)
            end do
         end do

      end do
!!!!$OMP end do

    return
  END SUBROUTINE sqr_realspcal_y2x

  

    
END MODULE RMHDS_pssn
