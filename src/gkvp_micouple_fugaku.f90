MODULE GKV_micouple


!-------------------------------------------------------------------------------
!
!    Modules for the magnetosphere-ionosphere feedback coupling
!
!    --------------
!      for gkvp_f0.57 (T.-H. Watanabe, Sep 2021)
!        - No overlap is assumed for zv data transfer
!        - The current version of pssn_ routines in RMHDS does not 
!          support the domain decomposition for ny. Thus, nprocw should be 1.
!
!    -------------
!      the limitation noted above is resolved. ( K. Fujita )
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
!  use GKV_fld,   only: fld_esfield, fld_emfield_hh, fld_hh2ff
  use GKV_clock, only: clock_sta, clock_end

  use RMHDS_pssn, only: pssn_brackets, pssn_divgrad
  use GKV_fft, only:    plan_x_backward,  plan_y_backward

  implicit none

  private

    complex(kind=DP), save, dimension(-nx:nx,0:ny) :: iphi, icpr, idns, ipsi
    complex(kind=DP), save, dimension(-nx:nx,0:ny) :: didns, qidns
    real(kind=DP),    save, dimension(-nx:nx,0:ny) :: ksqi
    real(kind=DP),    save, dimension(-nx:nx,0:ny) :: fct_jpara, fct_poissoni
 
    real(kind=DP),    save :: valf  !
                            ! rho2R0 is changed to a public variable

! --- parameters for M-I coupling
    real(kind=DP), save :: e0,   &    ! convection electric field
                           idns0,&    ! ionospheric number density
                           mp,   &    ! Pedersen mobility normalized by the ExB mobility
                           mh,   &    ! Hall mobility normalized by the ExB mobility
                           dperp,&    ! diffusion coefficient
                           alpha,&    ! normalized recombination rate : alpha * n0
                           zl,   &    ! field line length
                           rho2L      ! rho_ti / L_MHD
! Fujita added
    character(6)  :: irtrn                     ! function the ionospheric equations are solved for, phi or psi
    real(kind=DP), save :: nu_H                ! hyperviscosity coefficient
    integer :: hvpwr                           ! power of hyperviscosity: k^(2*hvpwr)
    integer :: np                              ! number of spatial points to output 5D delta-f
    integer(kind=DP), allocatable, dimension(:), save :: xp, yp
    integer, parameter :: nbuff = ((2*nz)*(nm+1)-1)/nprocw + 1
!
    integer, save :: nchunk_zm = 1, nchunk_yb  = 1, nchunk_xb = 1
    integer, save :: nchunk_zv = 1, nchunk_yzv = 1, nchunk_yz = 1

! --- unit numbers for I/O
    integer, parameter :: icnt_iono = 120, &
                          ophi_iono = 130, &
                          ocnt_iono = 150, &
                          oeng_idns = 160, &
                          oeng_iphi = 161, &
                          oeng_icpr = 162, &
                          obln_iono = 163, &
                          oroot     = 164

! Table for the plasma dispersion function
    complex(kind=DP), parameter :: coefb(8) = (/ ( -0.017340124574718_DP,  -0.04630639291680_DP  ), &
                                                 ( -0.739916992322501_DP,   0.83951799780998_DP  ), &
                                                 (  5.840628642184073_DP,   0.95360090576437_DP  ), &
                                                 ( -5.583371525286853_DP, -11.20854319126599_DP  ), &
                                                 ( -0.017340124574718_DP,   0.04630639291680_DP  ), &
                                                 ( -0.739916992322501_DP,  -0.83951799780998_DP  ), &
                                                 (  5.840628642184073_DP,  -0.95360090576437_DP  ), &
                                                 ( -5.583371525286853_DP,  11.20854319126599_DP  ) /)
    complex(kind=DP), parameter :: coefc(8) = (/ (  2.237687789201900_DP, -1.625940856173727_DP ), &
                                                 (  1.465234126106004_DP, -1.789620129162444_DP ), &
                                                 (  0.839253981723264_DP, -1.891995045765206_DP ), &
                                                 (  0.273936222628556_DP, -1.941786875844713_DP ), &
                                                 ( -2.237687789201900_DP, -1.625940856173727_DP ), &
                                                 ( -1.465234126106004_DP, -1.789620129162444_DP ), &
                                                 ( -0.839253981723264_DP, -1.891995045765206_DP ), &
                                                 ( -0.273936222628556_DP, -1.941786875844713_DP ) /)

  integer, parameter :: gny = (ny+1)*nprocw-1
  integer, parameter :: nk = min( nx, gny )

!... dispersion relation
  real(kind=DP) :: zet0, dzet
  integer, parameter :: root_max = 50
  real(kind=DP), parameter :: lz = pi

  real(kind=DP), public :: rho2R0
  public   micouple_prep, micouple_init, micouple_close
  public   micouple_literm_zv, micouple_advnc, micouple_bndry_bound_f, micouple_bndry_bound_e
  public   micouple_out_cntrl 
  public   wrt_fxyz

CONTAINS


!--------------------------------------
  SUBROUTINE micouple_prep( mxi, myi, omi, kzi )
!--------------------------------------

    complex(kind=DP), intent(out) :: omi, kzi
    integer, intent(in) :: mxi, myi

    complex(kind=DP) :: idns_egn
    real(kind=DP), allocatable, dimension(:) :: xpr, ypr
    real(kind=DP) :: dnp
    integer ::  mx, my, gmy, is, ip

! additional variables for M-I coupling
    namelist /mip/    irtrn, &   ! phi or psi, passed from I to M
                      np,    &   ! number of sptial points to output f(x,y,z,vz,mu)
                      hvpwr, &   ! power of hyperviscosity: k^(2*hvpwr)
                      nu_H       ! hyperviscosity coefficient


! M-I coupling parameters in the MHD units
    namelist /ionop/  e0,   &   ! convection electric field
                      idns0,&   ! ionospheric number density
                      mp,   &   ! Pedersen mobility normalized by the ExB mobility
                      mh,   &   ! Hall mobility normalized by the ExB mobility
                      dperp,&   ! diffusion coefficient
                      alpha,&   ! normalized recombination rate : alpha * n0
                      zl,   &   ! field line length
                      rho2L     ! rho_ti / L_MHD

      iphi(:,:) = ( 0._DP, 0._DP )
      icpr(:,:) = ( 0._DP, 0._DP )
      idns(:,:) = ( 0._DP, 0._DP )
      ipsi(:,:) = ( 0._DP, 0._DP )

      didns(:,:) = ( 0._DP, 0._DP )
      qidns(:,:) = ( 0._DP, 0._DP )

      do my = ist_y, iend_y
        do mx = -nx, nx
          if (  rankw == 0  .and. mx == 0  .and.  my == 0 ) then
            ksqi(0,0) = 0._DP
          else
            ksqi(mx,my) = 1._DP / ksq(mx,my,-nz)
          end if
        end do
      end do

! Keep electron FLR
!      if( ranks == 0 ) then   ! drift kinetic electron
!        j0(:,:,:,:) = 1._DP
!        j1(:,:,:,:) = 0._DP
!        j2(:,:,:,:) = 0._DP
!        g0(:,:,:)   = 1._DP
!      end if
! Keep electron FLR

      read(inml,nml=mip)

        write( olog, * ) ""
        write( olog, * ) " # Variables for M-I coupling"
        write( olog, * ) " # iono eqs are solved for    = ", trim(irtrn)
        write( olog, * ) " # power of hyperviscosity    = ", hvpwr
        write( olog, * ) " # hyperviscosity coefficient = ", nu_H
        write( olog, * ) ""

      read(inml,nml=ionop)

        write( olog, * ) ""
        write( olog, * ) " # Ionospheric parameters in MHD units"
        write( olog, * ) " # e0           = ", e0
        write( olog, * ) " # idns0        = ", idns0
        write( olog, * ) " # mp           = ", mp
        write( olog, * ) " # dperp        = ", dperp
        write( olog, * ) " # alpha        = ", alpha
        write( olog, * ) " # zl           = ", zl
        write( olog, * ) ""
        write( olog, * ) " # Scaling factors"
        write( olog, * ) " # rho/L_MHD    = ", rho2L
        write( olog, * ) ""

        valf   = 1._DP / sqrt( beta )      ! The bulk ion is assumed to be the reference
        rho2R0 = rho2L * 2._DP * pi / zl
!
        write( olog, * ) " # rho/R0       = ", rho2R0
        write( olog, * ) ""
        write( olog, * ) " # Alfven speed / v_ref = ", valf
        write( olog, * ) " # Alfven transit time tau_A "
        write( olog, * ) " #     = 2*lz*R0 / valf = ", &
                         (zz(0)-zz(-nz))*real(2*nprocz,kind=DP) / valf
        write( olog, * ) ""
        write( olog, * ) " # Particle transit time "
          do is = 0, ns-1
            write( olog, * ) " # is, tau_s = ", is, &
                         (zz(0)-zz(-nz))*real(2*nprocz,kind=DP) / sqrt( tau(is) / Anum(is) )
          end do
        write( olog, * ) ""

! Fujita added
!... Find solutions of the dispersion relation

       if ( inum == 1 ) then
         if( calc_disp > 0 ) then ! calc_disp=0: OFF; =1: GK; =2: MHD 
!---$OMP PARALLEL DO default(shared) private(mx)
           do mx= -nx, nx
             call find_roots( 0, mx, myi, omi, kzi ) 
           end do
!---$OMP PARALLEL END DO
         end if
         call find_roots( 1, mxi, myi, omi, kzi ) 
       end if

! re-scaling
        e0           = e0    * valf / rho2R0
        idns0        = idns0 / valf
        dperp        = dperp * valf / ( rho2R0 * rho2L )
        alpha        = alpha * valf * rho2L / rho2R0

        write( olog, * ) ""
        write( olog, * ) " # Ionospheric parameters in gyrokinetic units"
        write( olog, * ) " # e0           = ", e0
        write( olog, * ) " # idns0        = ", idns0
        write( olog, * ) " # dperp        = ", dperp
        write( olog, * ) " # alpha        = ", alpha
        write( olog, * ) ""

! ... Poisions to write out f(x,y,z,vz,mu) 

      if( np > 0 ) then

        np  = min( np, nnx, nny )
        dnp = 1._DP / dble( np )
 
        allocate( xp(np), yp(np) )
        allocate( xpr(np), ypr(np) )

        write( olog, * ) " # Position of f in fxv/*.fzv.* "
        write( olog, '(a8,2a16)' ) " # index,  x/nnx,  y/nny "
        do ip=1, np

          ! cross section along y = -x/2 + 0.75 line 
          !   to be modifed for arbitary mp/mh ratio
          xpr(ip) = dnp * ( ip - 1 )
          ypr(ip) = - 0.5_DP * xpr(ip) + 0.75_DP 

          ! convert to integer indices
          xp(ip) = nint( xpr(ip) * nnx )
          yp(ip) = nint( ypr(ip) * nny )

          write( olog, '(i8,2ES16.6)' ) ip, xpr(ip), ypr(ip)  

        end do
        write( olog, * ) ""
 
        deallocate( xpr, ypr )
     
      else
        write( olog, * ) " # *.fzv.* files won't be output "
        write( olog, * ) ""
      end if

      flush( olog )


  END SUBROUTINE micouple_prep


!--------------------------------------
  SUBROUTINE micouple_init( phi, Al )
!--------------------------------------

    complex(kind=DP), intent(in), &
    dimension(-nx:nx,0:ny,-nz:nz-1)       :: phi, Al
    complex(kind=DP), dimension(:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:), allocatable :: nw
    real(kind=DP) :: del, dx, dy
    real(kind=DP) :: dt_max_dperp, dt_max_idns0
    real(kind=DP) :: time
    integer :: input_status
    integer  ::  mx, my, iz, iv, im, is
    character(6)   :: crank
    character(3)   :: cnew
    character(1)   :: srank

      write( cnew,  fmt="(i3.3)" ) inum

! dt_max check for dperp and idns0
      dx   = pi / ( kxmin_g * real( nxw, kind=DP ) )
      dy   = pi / ( kymin_g * real( nyw, kind=DP ) )
      del  = min( dx, dy )

      dt_max_dperp = del*del / dperp
      dt_max_idns0 = del*del / ( mp*idns0 )

      dt_max = min( dt_max, dt_max_dperp*0.5_DP, dt_max_idns0*0.5_DP )

      dt     = min( dt, dt_max )

      write( olog, * ) ""
      write( olog, * ) "# Now, dt_max might be changed ... "
      write( olog, * ) "# dt_max_dperp, dt_max_idns0 = ", dt_max_dperp, dt_max_idns0
      write( olog, * ) ""
      write( olog, * ) "# dt,           dt_max       = ", dt, dt_max
      write( olog, * ) ""

! factor for flux correction 

      allocate( wf(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny) )

      iz = -nz
        do im = 0, nm
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iv,im) = fcs(ranks) * vl(iv)**2 * j0(mx,my,iz,im)**2 * fmx(iz,iv,im)
              end do
            end do
          end do
        end do

      call iono_intgrl_v0_moment ( wf, nw, 0 )  ! No sum over species

        do my = ist_y, iend_y
          do mx = -nx, nx

            if ( nw(mx,my) /= 0._DP ) then
              fct_jpara(mx,my) = 1._DP / nw(mx,my)
!             fct_jpara(mx,my) = 1._DP + ksqi(mx,my) * beta * nw(mx,my)
            else
              fct_jpara(mx,my) =  0._DP
            end if

            if( fct_poisson(mx,my,iz) /= 0._DP ) then
              fct_poissoni(mx,my) = 1._DP / fct_poisson(mx,my,iz) 
            else
              fct_poissoni(mx,my) = 0._DP
            end if

          end do
        end do

        if ( rankw == 0 ) then
          fct_jpara(0,0)    = 0._DP
          fct_poissoni(0,0) = 0._DP
        end if

!! debug
!        write( olog, * ) "# fct_japara = "
!        write( olog, fmt="(1p,256e15.7)" ) ( fct_jpara(0,my), my = ist_y, iend_y )
!! debug


!      deallocate( wf )
      deallocate( nw )

! file open
      if( rankg == 0 ) then
    !    open( oeng_iono, file=trim(f_hst)//"eng_iono."//cnew )
        open( oeng_idns, file=trim(f_hst)//"eng_idns."//cnew )
        open( oeng_iphi, file=trim(f_hst)//"eng_iphi."//cnew )
        open( oeng_icpr, file=trim(f_hst)//"eng_icpr."//cnew )
        open( obln_iono, file=trim(f_hst)//"bln_iono."//cnew )
      end if

!      if( rankz == 0 .and. rankv == 0 .and. rankm == 0 .and. ranks == 0 ) then
      if( rankz == 0 ) then

        if ( inum > 1 ) then
          call iono_fileio_open_icnt( trim(f_cnt) )
        end if

          call iono_fileio_open_cnt( trim(f_cnt) )

        if( rankz == 0 .and. rankv == 0 .and. rankm == 0 .and. ranks == 0 ) then
          call iono_fileio_open_phi( trim(f_phi) )
        end if

        if( inum > 1 ) then

          time   = - 1._DP

          do
            call iono_fileio_read_cnt( idns, time, input_status )

            if ( input_status < 0 ) then
              write( olog, * ) &
                 " # end of file of unit=130 is detected --> stop"
              call flush(olog)
              call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
              stop
            end if

            if ( input_status > 0 ) then
              write( olog, * ) &
                 " # input error of unit=130 is detected --> stop"
              call flush(olog)
              call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
              stop
            end if

            if ( time > eps ) exit
          end do

          write( olog, * ) ""
          write( olog, * ) "# ionospheric density read at t = ", time
          write( olog, * ) ""

          iphi(:,:) = phi(:,:,-nz)

          call iono_cpr( iphi, idns, icpr, ipsi )

        else ! inum==1

          iphi(:,:) = phi(:,:,-nz)
          icpr(:,:) = ksq(:,:,-nz) * Al(:,:,-nz) / beta
          call idns_init( iphi, idns, icpr )

        end if

      end if

  END SUBROUTINE micouple_init

!--------------------------------------
  SUBROUTINE idns_init( iphi, idns, icpr )
!--------------------------------------
!    compute the initial ionospheric density

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: iphi, icpr
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: idns

    integer :: mx, my, iz

      iz = -nz

      do my = ist1_y, iend_y
        do mx = -nx, nx

          idns(mx,my) = ( beta * icpr(mx,my) +  mp * idns0 * ksq(mx,my,iz) * iphi(mx,my) ) &
                        / ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )  - ksq(mx,my,iz) * dperp )
        end do
      end do

      if ( rankw == 0 ) then
        my = 0
        do mx = 1, nx
          idns(mx,my) = ( beta * icpr(mx,my) +  mp * idns0 * ksq(mx,my,iz) * iphi(mx,my) ) &
                      / ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )  - ksq(mx,my,iz) * dperp )
          idns(-mx,-my) = conjg( idns(mx,my) )
        end do
        idns(0,0) = ( 0._DP, 0._DP )
      end if


  END SUBROUTINE idns_init





!--------------------------------------
  SUBROUTINE micouple_literm_zv ( ff, psi, im, lf )
!--------------------------------------
!     (z,v)-derivative of ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb)            :: psi
    integer, intent(in) :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)             :: lf

    real(kind=DP), dimension(-nz:nz-1) :: cefz, cefz2
    real(kind=DP) :: cefv, cs1
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1330)
                                         ! call fapp_start("literm_para",1330,1)
!$OMP end master

      cs1    = sgn(ranks) * Znum(ranks) / tau(ranks)
      do iz = -nz, nz-1
        cefz(iz)   = 1._DP / ( 12._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
        cefz2(iz)  = 1._DP / ( 60._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
      end do
      cefv   = 1._DP / ( 12._DP * dv ) * sqrt( tau(ranks) / Anum(ranks) )


      if (trim(z_calc) == "cf4") then

!$OMP do collapse(2) schedule(dynamic,nchunk_zv)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                 - vl(iv) * cefz(iz) * (              &
                     -         ff(mx,my,iz+2,iv)      &
                     + 8._DP * ff(mx,my,iz+1,iv)      &
                     - 8._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )    &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz(iz) * (          &
                         -         psi(mx,my,iz+2)    &
                         + 8._DP * psi(mx,my,iz+1)    &
                         - 8._DP * psi(mx,my,iz-1)    &
                         +         psi(mx,my,iz-2) ) )&
                 - art_diff * (                       &
                     +         ff(mx,my,iz+2,iv)      &
                     - 4._DP * ff(mx,my,iz+1,iv)      &
                     + 6._DP * ff(mx,my,iz  ,iv)      &
                     - 4._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )
            end do
          end do
        end do
      end do
!$OMP end do nowait

! fujita changed the name to keep the original cf4 
!      if (trim(z_calc) == "cf4") then 
      else if (trim(z_calc) == "ondsided") then

        if( rankz /= 0 ) then

!$OMP do collapse(2) schedule(dynamic,nchunk_zv)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz(iz) * (              &
                         -         ff(mx,my,iz+2,iv)      &
                         + 8._DP * ff(mx,my,iz+1,iv)      &
                         - 8._DP * ff(mx,my,iz-1,iv)      &
                         +         ff(mx,my,iz-2,iv) )    &
                     + mir(iz,im) * cefv * (              &
                         -         ff(mx,my,iz,iv+2)      &
                         + 8._DP * ff(mx,my,iz,iv+1)      &
                         - 8._DP * ff(mx,my,iz,iv-1)      &
                         +         ff(mx,my,iz,iv-2) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz(iz) * (          &
                             -         psi(mx,my,iz+2)    &
                             + 8._DP * psi(mx,my,iz+1)    &
                             - 8._DP * psi(mx,my,iz-1)    &
                             +         psi(mx,my,iz-2) ) )&
                     - art_diff * (                       &
                        +         ff(mx,my,iz+2,iv)      &
                         - 4._DP * ff(mx,my,iz+1,iv)      &
                         + 6._DP * ff(mx,my,iz  ,iv)      &
                         - 4._DP * ff(mx,my,iz-1,iv)      &
                         +         ff(mx,my,iz-2,iv) )
                  end do
                end do
              end do
            end do
!$OMP end do nowait
         else ! raknz==0
!*** ORG
          do iv = 1, 2*nv
            do iz = -nz+2, nz-1 

              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz(iz) * (              &
                     -         ff(mx,my,iz+2,iv)      &
                     + 8._DP * ff(mx,my,iz+1,iv)      &
                     - 8._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )    &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz(iz) * (          &
                         -         psi(mx,my,iz+2)    &
                         + 8._DP * psi(mx,my,iz+1)    &
                         - 8._DP * psi(mx,my,iz-1)    &
                         +         psi(mx,my,iz-2) ) )&
                 - art_diff * (                       &
                     +         ff(mx,my,iz+2,iv)      &
                     - 4._DP * ff(mx,my,iz+1,iv)      &
                     + 6._DP * ff(mx,my,iz  ,iv)      &
                     - 4._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )
               end do
             end do
           end do

! one-sided FD with the fourth-order accuracy

        iz = -nz+1 
          do my = ist_y, iend_y
            do mx = -nx, nx
              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                 - vl(iv) * cefz2(iz) * (             &
                     - 1._DP * ff(mx,my,iz+3,iv)      &
                     + 6._DP * ff(mx,my,iz+2,iv)      &
                     +18._DP * ff(mx,my,iz+1,iv)      &
                     +10._DP * ff(mx,my,iz  ,iv)      &
                     -33._DP * ff(mx,my,iz-1,iv) )    &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz2(iz) * (         &
                         - 1._DP * psi(mx,my,iz+3)    &
                         + 6._DP * psi(mx,my,iz+2)    &
                         +18._DP * psi(mx,my,iz+1)    &
                         +10._DP * psi(mx,my,iz  )    &
                         -33._DP * psi(mx,my,iz-1) ) )
            end do
          end do


! one-sided FD with the fourth-order accuracy
        iz = -nz
          do my = ist_y, iend_y
            do mx = -nx, nx
              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                 - vl(iv) * cefz(iz) * (              &
                     -  3._DP * ff(mx,my,iz+4,iv)     &
                     + 16._DP * ff(mx,my,iz+3,iv)     &
                     - 36._DP * ff(mx,my,iz+2,iv)     &
                     + 48._DP * ff(mx,my,iz+1,iv)     &
                     - 25._DP * ff(mx,my,iz  ,iv) )   &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz(iz) * (          &
                         -  3._DP * psi(mx,my,iz+4)   &
                         + 16._DP * psi(mx,my,iz+3)   &
                         - 36._DP * psi(mx,my,iz+2)   &
                         + 48._DP * psi(mx,my,iz+1)   &
                         - 25._DP * psi(mx,my,iz  ) ) )
            end do
          end do

        end do

!*** ORG one block END

! one-sided FD with the third-order accuracy
!        iz = -nz
!          do my = ist_y, iend_y
!            do mx = -nx, nx
!! third 2, -9, 18, -11 / 6
!              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
!                 - vl(iv) * cefz(iz) * 2._DP * (      &
!                     + 2._DP * ff(mx,my,iz+3,iv)      &
!                     - 9._DP * ff(mx,my,iz+2,iv)      &
!                     +18._DP * ff(mx,my,iz+1,iv)      &
!                     -11._DP * ff(mx,my,iz  ,iv) )    &
!                 + mir(iz,im) * cefv * (              &
!                     -         ff(mx,my,iz,iv+2)      &
!                     + 8._DP * ff(mx,my,iz,iv+1)      &
!                     - 8._DP * ff(mx,my,iz,iv-1)      &
!                     +         ff(mx,my,iz,iv-2) )    &
!                 - cs1 * fmx(iz,iv,im) * (            &
!                       vl(iv) * cefz(iz) * 2._DP * (  &
!                         + 2._DP * psi(mx,my,iz+3)    &
!                         - 9._DP * psi(mx,my,iz+2)    &
!                         +18._DP * psi(mx,my,iz+1)    &
!                         -11._DP * psi(mx,my,iz  ) ) )
!            end do
!          end do


! one-sided FD with the second-order accuracy
!        iz = -nz
!          do my = ist_y, iend_y
!            do mx = -nx, nx
!! third 2, -9, 18, -11 / 6
!              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
!                 - vl(iv) * cefz(iz) *  (             &
!                     - 6._DP * ff(mx,my,iz+2,iv)      &
!                     +24._DP * ff(mx,my,iz+1,iv)      &
!                     -18._DP * ff(mx,my,iz  ,iv) )    &
!                 + mir(iz,im) * cefv * (              &
!                     -         ff(mx,my,iz,iv+2)      &
!                     + 8._DP * ff(mx,my,iz,iv+1)      &
!                     - 8._DP * ff(mx,my,iz,iv-1)      &
!                     +         ff(mx,my,iz,iv-2) )    &
!                 - cs1 * fmx(iz,iv,im) * (            &
!                       vl(iv) * cefz(iz) * (          &
!                         - 6._DP * psi(mx,my,iz+2)    &
!                         +24._DP * psi(mx,my,iz+1)    &
!                         -18._DP * psi(mx,my,iz  ) ) )
!            end do
!          end do


!*...org      end do

        end if
!...$OMP end do nowait

      else if (trim(z_calc) == "up5") then
        do iv = 1, 2*nv
          if ( vl(iv) > 0._DP ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_yz)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         - 3._DP * ff(mx,my,iz+2,iv)      &
                         +30._DP * ff(mx,my,iz+1,iv)      &
                         +20._DP * ff(mx,my,iz  ,iv)      &
                         -60._DP * ff(mx,my,iz-1,iv)      &
                         +15._DP * ff(mx,my,iz-2,iv)      &
                         - 2._DP * ff(mx,my,iz-3,iv) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz2(iz) * (         &
                              1._DP * psi(mx,my,iz+3)    &
                           -  9._DP * psi(mx,my,iz+2)    &
                           + 45._DP * psi(mx,my,iz+1)    &
                           - 45._DP * psi(mx,my,iz-1)    &
                           +  9._DP * psi(mx,my,iz-2)    &
                           -  1._DP * psi(mx,my,iz-3) ) )
!                           vl(iv) * cefz(iz) * (          &
!                             -         psi(mx,my,iz+2)    &
!                             + 8._DP * psi(mx,my,iz+1)    &
!                             - 8._DP * psi(mx,my,iz-1)    &
!                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          else
!$OMP do collapse(2) schedule(dynamic,nchunk_yz)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         + 2._DP * ff(mx,my,iz+3,iv)      &
                         -15._DP * ff(mx,my,iz+2,iv)      &
                         +60._DP * ff(mx,my,iz+1,iv)      &
                         -20._DP * ff(mx,my,iz  ,iv)      &
                         -30._DP * ff(mx,my,iz-1,iv)      &
                         + 3._DP * ff(mx,my,iz-2,iv) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz2(iz) * (         &
                              1._DP * psi(mx,my,iz+3)    &
                           -  9._DP * psi(mx,my,iz+2)    &
                           + 45._DP * psi(mx,my,iz+1)    &
                           - 45._DP * psi(mx,my,iz-1)    &
                           +  9._DP * psi(mx,my,iz-2)    &
                           -  1._DP * psi(mx,my,iz-3) ) )
!                           vl(iv) * cefz(iz) * (          &
!                             -         psi(mx,my,iz+2)    &
!                             + 8._DP * psi(mx,my,iz+1)    &
!                             - 8._DP * psi(mx,my,iz-1)    &
!                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          end if
        end do

      else if (trim(z_calc) == "cf6") then

!$OMP do collapse(2) schedule(dynamic,nchunk_zv)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                       - vl(iv) * cefz2(iz) * (           &
                                   ff(mx,my,iz+3,iv)      &
                        -  9._DP * ff(mx,my,iz+2,iv)      &
                        + 45._DP * ff(mx,my,iz+1,iv)      &
                        - 45._DP * ff(mx,my,iz-1,iv)      &
                        +  9._DP * ff(mx,my,iz-2,iv)      &
                        -          ff(mx,my,iz-3,iv) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz2(iz) * (         &
                              1._DP * psi(mx,my,iz+3)     &
                           -  9._DP * psi(mx,my,iz+2)     &
                           + 45._DP * psi(mx,my,iz+1)     &
                           - 45._DP * psi(mx,my,iz-1)     &
                           +  9._DP * psi(mx,my,iz-2)     &
                           -  1._DP * psi(mx,my,iz-3) ) ) &
                     - art_diff * (                       &
                        +         ff(mx,my,iz+2,iv)       &
                        - 4._DP * ff(mx,my,iz+1,iv)       &
                        + 6._DP * ff(mx,my,iz  ,iv)       &
                        - 4._DP * ff(mx,my,iz-1,iv)       &
                        +         ff(mx,my,iz-2,iv) )

!                 - art_diff * (                           &
!                        -          ff(mx,my,iz+3,iv)      &
!                        +  6._DP * ff(mx,my,iz+2,iv)      &
!                        - 15._DP * ff(mx,my,iz+1,iv)      &
!                        + 20._DP * ff(mx,my,iz ,iv)       &
!                        - 15._DP * ff(mx,my,iz-1,iv)      &
!                        +  6._DP * ff(mx,my,iz-2,iv)      &
!                        -          ff(mx,my,iz-3,iv) )   
                  end do
                end do
              end do
            end do
!$OMP end do nowait
      
!...        org
!        write(olog,*) "## up5 is not available for z_bound = mi_couple"
!        call flush(olog)
!        call MPI_Finalize(ierr_mpi)
!        stop

      else

        write(olog,*) "## Illegal choice for z_calc!! ---> stop"
        write(olog,*) " # z_calc: ", trim(z_calc)
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if

!>> hyper viscosity for electrons (fujita added)
      if( ranks==0 ) then
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                 lf(mx,my,iz,iv) = lf(mx,my,iz,iv) &
                                 - ksq(mx,my,iz) ** hvpwr * nu_H * ff(mx,my,iz,iv) 
              end do
            end do
          end do
        end do
      end if
!<<
      

      
!$OMP master
                                         ! call fapp_stop("literm_para",1330,1)
                                           call clock_end(1330)
!$OMP end master


  END SUBROUTINE micouple_literm_zv


!--------------------------------------
  SUBROUTINE micouple_close
!--------------------------------------

! file open
      if( rankg == 0 ) then
        close( oeng_idns )
        close( oeng_iphi )
        close( oeng_icpr )
        close( obln_iono )
      end if

      if( rankz == 0 ) then

        if ( inum > 1 ) then
          call iono_fileio_close_icnt
        end if

        call iono_fileio_close_cnt
        call iono_fileio_close_phi

      end if


  END SUBROUTINE micouple_close


!--------------------------------------
  SUBROUTINE micouple_advnc( hh, ff, phi, Al, dA, istep )
!--------------------------------------
!     This routine should be called only from rankz = 0

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
! fujita added
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: dA
    complex(kind=DP),             &
      dimension(-nx:nx,0:ny) :: ipsi0
    real(kind=DP) :: dpsi
!

    integer, intent(in) :: istep

    real(kind=DP) :: c1, c2, cq, c0
    integer :: mx, my

!$OMP master
    call clock_sta(1630)   
!$OMP end master

      if( rankz /= 0 )  return

      if      ( istep == 1 ) then
        c1   =  0.5_DP
        c2   = -1._DP
        cq   = -2._DP
        c0   =  1._DP
      else if ( istep == 2 ) then
        c1   =  1._DP - sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 3 ) then
        c1   =  1._DP + sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 4 ) then
        c1   =  1._DP / 6._DP
        c2   = -1._DP / 3._DP
        cq   =  0._DP
        c0   =  0._DP
      end if


! time integration using old idns, icpr, and iphi
      call iono_dns( idns, icpr, iphi, didns )

        do my = ist_y, iend_y
          do mx = -nx, nx
            idns(mx,my)  = idns(mx,my) + c1 * didns(mx,my) + c2 * qidns(mx,my)
            qidns(mx,my) =               cq * qidns(mx,my) + c0 * didns(mx,my)
          end do
        end do

!--------------
! phi => jpara
!--------------


      if ( trim(irtrn) == "psi" ) then 

! map the new phi from the GK solution to iono
        iphi(:,:) = phi(:,:,-nz)

! save the old psi
        ipsi0(:,:) = Al(:,:,-nz)

! compute the new icpr and ipsi
        call iono_cpr( iphi, idns, icpr, ipsi )

! reverse mapping to GK
        Al(:,:,-nz) = ipsi(:,:)

! adjust the 1st order moment of ff and hh to icpr
        call iono_correct_flux( ff, Al, icpr, hh )

      else if ( trim(irtrn) == "phi" ) then

!--------------
! jpara => phi
!--------------
!
! map the new Al from the GK solution to iono
        ipsi(:,:) = Al(:,:,-nz)
!
! compute the new icpr and ipsi
        call iono_phi( ipsi, idns, icpr, iphi )
!
! reverse mapping to GK
        phi(:,:,-nz) = iphi(:,:)
!
! adjust the 0th order moment of ff and hh to icpr
        call iono_correct_dns( ff, Al, iphi, hh )

    end if 

!$OMP master
    call clock_end(1630)
!$OMP end master

  END SUBROUTINE micouple_advnc


!--------------------------------------
  SUBROUTINE iono_dns( idns, icpr, iphi, didns )
!--------------------------------------
!     increment within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny)  :: idns, icpr, iphi, didns

    complex(kind=DP), &
      dimension(-nx:nx,0:ny)  :: ipsnb, idnsq

    integer :: mx, my

      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( iphi, idns, ipsnb,   1 )
!!        call nonlinear_sqr( idns, idnsq, 1 ) 
        
      else

        ipsnb(:,:) = ( 0._DP, 0._DP )
!!        idnsq(:,:) = ( 0._DP, 0._DP )
        
      end if

      didns(:,:) = ( icpr(:,:) * beta - 2.d+0 * alpha * idns(:,:) &
                   - ipsnb(:,:) ) * dt
      
! --- reality condition
      if ( rankw == 0 ) then
        my = 0
          do mx = 1, nx
            didns(mx,my) = conjg( didns(-mx,-my) )
          end do
          didns(0,0) = ( 0._DP, 0._DP )
       end if

  END SUBROUTINE iono_dns


!--------------------------------------
  SUBROUTINE iono_cpr( iphi, idns, icpr, ipsi )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: iphi, idns
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: icpr, ipsi

    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: ipsnb, idvgd

    real(kind=DP) :: iomg_abs
    integer :: mx, gmy, my, iz

      iz = -nz

      if( trim(calc_type) == "nonlinear" ) then
        call pssn_brackets( iphi, idns, ipsnb, 1 )
        call pssn_divgrad ( idns, iphi, idvgd, 1 )
!!        ipsnb(:,:) = ( 0._DP, 0._DP )
!!        idvgd(:,:) = ( 0._DP, 0._DP )
      else
        ipsnb(:,:) = ( 0._DP, 0._DP )
        idvgd(:,:) = ( 0._DP, 0._DP )
      end if

      do my = ist_y, iend_y
        do mx = -nx, nx
          icpr(mx,my) = ( - mp * idns0 * ksq(mx,my,iz) * iphi(mx,my)  &
                          + ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )    &
                              - ksq(mx,my,iz) * dperp ) * idns(mx,my) &
                          + mp * idvgd(mx,my) + ipsnb(mx,my)          &
                        ) / beta

!          ipsi(mx,my) = - ksqi(mx,my) * icpr(mx,my)
! The sign of the Ampare's law follows the convenciton of the GK equation
          ipsi(mx,my) = ksqi(mx,my) * icpr(mx,my) * beta

! evaluate div.grad-phi, which appears in idvgd
          !!iomg_abs = abs( ksq(mx,my,iz)*iphi(mx,my) )
          !!if( iomg_abs > icpr_max ) icpr_max = iomg_abs

        end do
      end do



  END SUBROUTINE iono_cpr


!--------------------------------------
  SUBROUTINE iono_phi ( ipsi, idns, icpr, iphi )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: ipsi, idns
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: icpr!!, iphi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny)  :: iphi
    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: ipsnb, idvgd

!! fujita added for iterative calc. of  nonlinear terms
    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: iphi0, iphi_ref, iphi_min
    real(kind=DP) :: diff, diff_max0, diff_max, diff_min
    integer, parameter :: it_max = 500
    integer :: it
!!
    integer :: mx, my, iz

      iz = -nz

! The sign of the Ampare's law follows the convenciton of the GK equation
      icpr(:,:) = ksq(:,:,iz) * ipsi(:,:) / beta

      do my = ist_y, iend_y
        do mx = -nx, nx
          iphi(mx,my)  = ( - beta * icpr(mx,my)                        &
                           + ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )    &
                               - ksq(mx,my,iz) * dperp ) * idns(mx,my) &
                         ) / ( mp * idns0 ) * ksqi(mx,my)
        end do
      end do

      if( trim(calc_type) == "nonlinear" ) then

        iphi0(:,:)    = iphi(:,:)
        
        do it=1, it_max

          iphi_ref(:,:) = iphi(:,:)
          diff_max0     = 0._DP

          call pssn_brackets( iphi, idns, ipsnb, 1 )
          call pssn_divgrad ( idns, iphi, idvgd, 1 )

          do my = ist_y, iend_y
            do mx = -nx, nx
              iphi(mx,my) = iphi0(mx,my) &
                          + ( mp * idvgd(mx,my) + ipsnb(mx,my) ) / ( mp * idns0 )  * ksqi(mx,my)

              diff = abs( iphi(mx,my) / iphi_ref(mx,my) - 1._DP )
              if( abs(iphi_ref(mx,my)) > 0._DP .and. diff>diff_max0 ) diff_max0 = diff

            end do
          end do

          call MPI_Allreduce( diff_max0, diff_max, 1, MPI_DOUBLE_PRECISION, &
                              MPI_MAX, fft_comm_world, ierr_mpi )

          if( diff_max < 1.d-8 ) exit

          if( it==it_max ) then
            write(olog,*) "!! iteration count max is reached", diff_max
            flush(olog)
          end if

        end do

      end if


!... org
!
!      iz = -nz
!
!
!      if( trim(calc_type) == "nonlinear" ) then
!
!        call pssn_brackets( iphi, idns, ipsnb, 1 )
!        call pssn_divgrad ( idns, iphi, idvgd, 1 )
!!!!!        ipsnb(:,:) = ( 0._DP, 0._DP )
!!!!        idvgd(:,:) = ( 0._DP, 0._DP )
!
!      else
!        ipsnb(:,:) = ( 0._DP, 0._DP )
!        idvgd(:,:) = ( 0._DP, 0._DP )
!      end if

! The sign of the Ampare's law follows the convenciton of the GK equation
!      icpr(:,:) = ksq(:,:,iz) * ipsi(:,:) / beta
!
!      do my = ist_y, iend_y
!        do mx = -nx, nx
!          iphi(mx,my)  = ( - beta * icpr(mx,my)                        &
!                           + ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) )    &
!                               - ksq(mx,my,iz) * dperp ) * idns(mx,my) &
!                           + mp * idvgd(mx,my) + ipsnb(mx,my)          &
!                         ) / ( mp * idns0 ) * ksqi(mx,my)
!        end do
!      end do


  END SUBROUTINE iono_phi


!--------------------------------------
  SUBROUTINE iono_correct_flux ( ff, Al, icpr, hh )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny)  :: icpr
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

    complex(kind=DP), dimension(:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:), allocatable :: iflx_old, iflx_new, mcpr
    complex(kind=DP), dimension(:,:,:,:), allocatable :: ff_iono, ff_odd, ff_evn

    integer  ::  mx, my, iz, iv, im

      iz = -nz

! flux given by GK equaitons

      allocate( wf(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( iflx_old(-nx:nx,0:ny) )
      allocate( iflx_new(-nx:nx,0:ny) )
      allocate( mcpr(-nx:nx,0:ny) )

      allocate( ff_iono(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( ff_odd(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( ff_evn(-nx:nx,0:ny,1:2*nv,0:nm) )

!... fujita added
!      allocate( iflx_dif(-nx:nx,0:ny) )



      iz = -nz
!$OMP parallel 
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks) &
                                * sqrt( tau(ranks) / Anum(ranks) ) * vl(iv)

                ff_iono(mx,my,iv,im) = ff(mx,my,iz,iv,im)
              end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


      call iono_intgrl_v0_moment ( wf, iflx_old, 0 )  ! No sum over species

      call iono_intgrl_v0_moment ( wf, mcpr,     1 )  ! Sum over species


        iflx_old(:,:) = iflx_old(:,:) / ( sgn(ranks) * fcs(ranks) )

!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                 write(olog,'(2i6,4e16.6)') mx, my, iflx_old(mx,my), mcpr(mx,my)
!              end do
!           end do
!           write(olog,*)


!! icpr - mcpr
!      if( ranks == 0 ) then     ! supposed to be electrons      
!        iflx_new(:,:) = ( icpr(:,:) - mcpr(:,:) ) / ( sgn(ranks) * fcs(ranks) )
!      else if( ranks == 1 ) then     ! supposed to be ions
!        iflx_new(:,:) = ( 0._DP, 0._DP )
!      end if

! no ion current
      if( ranks == 0 ) then     ! supposed to be electrons      
        iflx_new(:,:) = ( icpr(:,:) - mcpr(:,:) ) / ( sgn(ranks) * fcs(ranks) ) &
                      + iflx_old(:,:)
      else if( ranks == 1 ) then     ! supposed to be ions
        iflx_new(:,:) = iflx_old(:,:)
      end if


      if ( rankw == 0 ) then
        iflx_old(0,0) = ( 0._DP, 0._DP )       !  zero-zero
        iflx_new(0,0) = ( 0._DP, 0._DP )       !  zero-zero
      end if


!------------------------------
! with the odd components of hh
!------------------------------

!---$OMP parallel
      do im = 0, nm
        call iono_symmetry_ff ( ff_iono(:,:,:,im), ff_odd(:,:,:,im), ff_evn(:,:,:,im) )
      end do
!---$OMP end parallel

      iz = -nz

!... correction 
!--$OMP parallel
!      do my = ist_y, iend_y
!        do mx = -nx, nx
!          iflx_dif(mx,my) = iflx_new(mx,my) - iflx_old(mx,my)
!        end do
!      end do
!---$OMP end parallel

!! estimate the correction size
!!      do my = ist_y, iend_y
!!        do mx = -nx, nx
!!          if( abs( iflx_old(mx,my) ) > 1.d-10 ) then
!!            icor_est = abs( iflx_new(mx,my) / iflx_old(mx,my) - 1._DP )
!!            if( icor_max<icor_est ) icor_max = icor_est
!!          end if
!!        end do
!!     end do

!$OMP parallel 
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
!... org
                ff(mx,my,iz,iv,im) = ff_evn(mx,my,iv,im)                      &
                                   + iflx_new(mx,my) * vl(iv) * fmx(iz,iv,im) &
                                   * sqrt( Anum(ranks) / tau(ranks) )
  !  
!           ff(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) &
!                                    + iflx_dif(mx,my) * vl(iv) * fmx(iz,iv,im) &
!                                    * sqrt( Anum(ranks) / tau(ranks) )


             end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

!! debug
!     mx = 0
!     my = 5
!     iz = -nz
!     im = 0
!     write( olog, * ) "# ff, ff_evn, ff_odd "
!     do iv = 1, 2*nv
!       write( olog, fmt="(1p,7e15.7)" ) vl(iv), ff_iono(mx,my,iv,im), ff_evn(mx,my,iv,im), ff_odd(mx,my,iv,im)
!     end do
!     write( olog, * )
!! debug


!------------------------------
!  Reflect the correction to hh
!------------------------------
      call iono_fld_ff2hh( ff, Al, hh )


      deallocate( wf )
      deallocate( iflx_old )
      deallocate( iflx_new )
! fujita added
!      deallocate( iflx_dif )


  END SUBROUTINE iono_correct_flux


!--------------------------------------
  SUBROUTINE iono_correct_dns ( ff, Al, iphi, hh )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny)  :: iphi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh

    complex(kind=DP), dimension(:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:), allocatable :: idns_old, idns_new, mrho
    complex(kind=DP), dimension(:,:,:,:), allocatable :: ff_iono, ff_odd, ff_evn
! fujita added

    integer  ::  mx, my, iz, iv, im

      iz = -nz

! flux given by GK equaitons

      allocate( wf(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( idns_old(-nx:nx,0:ny) )
      allocate( idns_new(-nx:nx,0:ny) )
      allocate( mrho(-nx:nx,0:ny) )

      allocate( ff_iono(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( ff_odd(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( ff_evn(-nx:nx,0:ny,1:2*nv,0:nm) )

      iz = -nz
!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks)

                ff_iono(mx,my,iv,im) = ff(mx,my,iz,iv,im)
              end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


      call iono_intgrl_v0_moment ( wf, idns_old, 0 )  ! No sum over species

      call iono_intgrl_v0_moment ( wf, mrho,     1 )  ! Sum over species


!!! check !!!
        idns_old(:,:) = idns_old(:,:) / ( sgn(ranks) * fcs(ranks) )


      iz = -nz
! no ion current
      if( ranks == 0 ) then     ! supposed to be electrons      
        idns_new(:,:) = ( iphi(:,:) * fct_poissoni(:,:) - mrho(:,:) ) &
                      / ( sgn(ranks) * fcs(ranks) )                   &
                      + idns_old(:,:)
!        idns_new(:,:) = ( iphi(:,:) * fct_poissoni(:,:) ) &
!                      / ( sgn(ranks) * fcs(ranks) )
      else if( ranks == 1 ) then     ! supposed to be ions
        idns_new(:,:) = idns_old(:,:)
!        idns_new(:,:) = ( 0._DP, 0._DP )
      end if

!!! check !!!

      if ( rankw == 0 ) then
        idns_old(0,0) = ( 0._DP, 0._DP )       !  zero-zero
        idns_new(0,0) = ( 0._DP, 0._DP )       !  zero-zero
      end if

!------------------------------
! with the odd components of hh
!------------------------------

!---$OMP parallel
      do im = 0, nm
        call iono_symmetry_ff ( ff_iono(:,:,:,im), ff_odd(:,:,:,im), ff_evn(:,:,:,im) )
      end do
!---$OMP end parallel


      iz = -nz
!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
!!                 ff(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) 
!                ff(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)               &
                ff(mx,my,iz,iv,im) = ff_odd(mx,my,iv,im)             &
                                   + idns_new(mx,my) * fmx(iz,iv,im) 
              end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


!------------------------------
!  Reflect the correction to hh
!------------------------------
      call iono_fld_ff2hh( ff, Al, hh )


      deallocate( wf )
      deallocate( idns_old )
      deallocate( idns_new )



  END SUBROUTINE iono_correct_dns


!--------------------------------------
  SUBROUTINE iono_fld_esfield ( ff, phi )
!--------------------------------------
!     electrostatic field calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv,0:nm) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny) :: phi

    complex(kind=DP), dimension(:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:), allocatable :: nw, ww
    integer  ::  mx, my, iz, iv, im


      allocate( wf(-nx:nx,0:ny,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny) )
      allocate( ww(-nx:nx,0:ny) )

      iz = -nz

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iv,im) = ff(mx,my,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks)
              end do
            end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


      call iono_intgrl_v0_moment ( wf, nw, 1 )    ! Sum over species


      if ( rankw == 0 ) then
          nw(0,0) = ( 0._DP, 0._DP )       !  zero-zero
      end if

       iz = -nz
!$OMP parallel do
          do my = ist_y, iend_y
            do mx = -nx, nx
              phi(mx,my) = nw(mx,my) * fct_poisson(mx,my,iz)
            end do
          end do

      deallocate( wf )
      deallocate( nw )
      deallocate( ww )
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

      if ( rankw == 0 ) then
          phi(0,0) = ( 0._DP, 0._DP )
      end if


  END SUBROUTINE iono_fld_esfield


!--------------------------------------
  SUBROUTINE iono_fld_ff2hh ( ff, Al, hh )
!--------------------------------------
!     ff -> hh

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh

    integer :: mx, my, iz, iv, im


    iz = -nz

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                hh(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)  &
                    + sgn(ranks) * Znum(ranks)  / sqrt( Anum(ranks) * tau(ranks))  &
                    * fmx(iz,iv,im) * vl(iv) * j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


  END SUBROUTINE iono_fld_ff2hh


!--------------------------------------
  SUBROUTINE iono_intgrl_v0_moment ( wf, wn, isw )
!--------------------------------------
!     Calculate the zeroth order velocity moment of wf

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)               :: wn
    integer, intent(in) :: isw

    complex(kind=DP), dimension(:,:), allocatable :: ww
    real(kind=DP) :: cef
    integer :: mx, my, iz, iv, im

      iz = -nz

      cef   = dv * twopi

      allocate( ww(-nx:nx,0:ny) )

      if( rankm == 0 ) then

!---$OMP parallel default(none) &
!---$OMP shared(ww,wf,vp,dvp,cef,ist_y,iend_y) &
!---$OMP private(mx,my,iz,iv,im)
!---$OMP workshare
      ww(:,:) = ( 0._DP, 0._DP )
!---$OMP end workshare

!---$OMP do collapse(2)
          do my = ist_y, iend_y

            do im = 1, nm
              do iv = 1, 2*nv
                do mx = -nx, nx
                  ww(mx,my) = ww(mx,my)                         &
                        + wf(mx,my,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do

          ! for edge compensation
           im = 1
              do iv = 1, 2*nv
                do mx = -nx, nx
                  ww(mx,my) = ww(mx,my)                             &
                        - ( - wf(mx,my,iv,im  ) * vp(iz,im  ) / 12._DP &
                          + ( wf(mx,my,iv,im+1) * vp(iz,im+1)          &
                            - wf(mx,my,iv,im  ) * vp(iz,im  ) * 2._DP  &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)


                end do
              end do
  
          end do
!---$OMP end do
!---$OMP end parallel

      else

!---$OMP parallel default(none) &
!---$OMP shared(ww,wf,vp,dvp,cef,ist_y,iend_y) &
!---$OMP private(mx,my,iz,iv,im)
!---$OMP workshare
      ww(:,:) = ( 0._DP, 0._DP )
!---$OMP end workshare

!---$OMP do collapse(2)
!        do iz = -nz, nz-1
          do my = ist_y, iend_y
  
            do im = 0, nm
              do iv = 1, 2*nv
                do mx = -nx, nx
                  ww(mx,my) = ww(mx,my)                         &
                        + wf(mx,my,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do
  
          end do
!        end do
!---$OMP end do
!---$OMP end parallel

      end if


      if( isw == 0 ) then 
        call MPI_Allreduce( ww, wn, nxy, MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, vel_comm_world, ierr_mpi )    ! No sum over species
      else
        call MPI_Allreduce( ww, wn, nxy, MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, spc_comm_world, ierr_mpi )    ! Sum over species
      end if


      deallocate( ww )


  END SUBROUTINE iono_intgrl_v0_moment


!--------------------------------------
  SUBROUTINE micouple_bndry_bound_f( ff )
!--------------------------------------

   complex(kind=DP), intent(inout), &
     dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff

   complex(kind=DP), dimension(:,:,:,:), allocatable :: zb1, zb2


     allocate( zb1(-nx:nx,0:ny,0:nzb-1,1:2*nv) )
     allocate( zb2(-nx:nx,0:ny,0:nzb-1,1:2*nv) )

!$OMP CRITICAL
       call iono_bndry_zvswap_buffin( ff, zb1 )

       call iono_bndry_zvswap_sendrecv( zb1, zb2 )

       call iono_bndry_zvswap_buffout ( zb2, ff )
!$OMP END CRITICAL

     deallocate( zb1 )
     deallocate( zb2 )


  END SUBROUTINE micouple_bndry_bound_f


!--------------------------------------
  SUBROUTINE iono_bndry_zvswap_buffin( ff, zb1_top )
!--------------------------------------
! This routine should be called after calling the 
! bndry routiens for ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_top

    integer :: iz, iv

      zb1_top(:,:,:,:) = ( 0._DP, 0._DP )

      if( rankz /= nprocz-1 ) return

!$OMP master
                                           call clock_sta(1351)
                                         ! call fapp_start("literm_boundf_bufferin",1351,1)
!$OMP end master

!---$OMP do collapse(2) schedule(dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            zb1_top   (:,:,iz,iv) = ff(:,:, nz-nzb+iz,iv)
          end do
        end do
!---$OMP end do nowait


!---$OMP end master


  END SUBROUTINE iono_bndry_zvswap_buffin


!--------------------------------------
  SUBROUTINE iono_bndry_zvswap_sendrecv ( zb1_top, zb2_top )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_top
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_top

    integer :: slngz, iswp
    integer, dimension(2) :: ireq
    integer, dimension(MPI_STATUS_SIZE,2) :: istatus

      zb2_top(:,:,:,:) = ( 0._DP, 0._DP )

      if( rankz /= nprocz-1 ) return

      slngz  = (2*nx+1)*(ny+1)*(2*nv) * nzb

      if( rankz == nprocz-1 ) then
        iswp   = rank + nprocw * nprocz * ( nprocv - 1 - 2*rankv ) !!* ( nprocs + 1 )
      else
        iswp   = MPI_PROC_NULL
      end if

!!debug
!        write( olog, * ) "# rank, rankv, iswp = ", rank, rankv, iswp
!        flush(olog)
!!debug

      call MPI_irecv( zb2_top,    slngz, MPI_DOUBLE_COMPLEX, iswp, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_isend( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, iswp, 1, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_waitall( 2, ireq, istatus, ierr_mpi )

!debug
!      zb2_top(:,:,:,:) = zb1_top(:,:,:,:)
!debug


  END SUBROUTINE iono_bndry_zvswap_sendrecv


!--------------------------------------
  SUBROUTINE iono_bndry_zvswap_buffout ( zb2_top, ff )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_top
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff

    complex(kind=DP) :: dfdz

    integer :: mx, my, iz, iv


      if( rankz == 0 ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx

! ...org
! 1st order extrapolation -- this method provide the best results in growth rate

!                  ff(mx,my,-nz-1,iv) = ff(mx,my,-nz  ,iv)
!                  ff(mx,my,-nz-2,iv) = ff(mx,my,-nz  ,iv) - ( ff(mx,my,-nz+1,iv) - ff(mx,my,-nz  ,iv) )
!                  ff(mx,my,-nz-3,iv) = ff(mx,my,-nz  ,iv) - ( ff(mx,my,-nz+1,iv) - ff(mx,my,-nz  ,iv) )*2._DP


!... modified: 2nd order
!                   ff(mx,my,-nz-1,iv) = 3._DP * ( ff(mx,my,-nz,iv) - ff(mx,my,-nz+1,iv) ) + ff(mx,my,-nz+2,iv)
!                   ff(mx,my,-nz-2,iv) = 3._DP * ( ff(mx,my,-nz-1,iv) - ff(mx,my,-nz,iv) ) + ff(mx,my,-nz+1,iv)
!                   ff(mx,my,-nz-3,iv) = 3._DP * ( ff(mx,my,-nz-2,iv) - ff(mx,my,-nz-1,iv) ) + ff(mx,my,-nz,iv)

!                   ff(mx,my,-nz-1,iv) = ff(mx,my,-nz,iv)
!                   ff(mx,my,-nz-2,iv) = 3._DP * ( ff(mx,my,-nz-1,iv) - ff(mx,my,-nz,iv) ) + ff(mx,my,-nz+1,iv)
!                   ff(mx,my,-nz-3,iv) = 3._DP * ( ff(mx,my,-nz-2,iv) - ff(mx,my,-nz-1,iv) ) + ff(mx,my,-nz,iv)


!... 4-point Lagrange
                   ff(mx,my,-nz-1,iv) = 4._DP * ff(mx,my,-nz  ,iv) - 6._DP * ff(mx,my,-nz+1,iv) &
                                      + 4._DP * ff(mx,my,-nz+2,iv) - 1._DP * ff(mx,my,-nz+3,iv)
                   ff(mx,my,-nz-2,iv) = 4._DP * ff(mx,my,-nz-1,iv) - 6._DP * ff(mx,my,-nz  ,iv) &
                                      + 4._DP * ff(mx,my,-nz+1,iv) - 1._DP * ff(mx,my,-nz+2,iv)
                   ff(mx,my,-nz-3,iv) = 4._DP * ff(mx,my,-nz-2,iv) - 6._DP * ff(mx,my,-nz-1,iv) &
                                      + 4._DP * ff(mx,my,-nz  ,iv) - 1._DP * ff(mx,my,-nz+1,iv)

!... 5-point Lagrange
!                   ff(mx,my,-nz-1,iv) =  5._DP * ff(mx,my,-nz  ,iv) - 10._DP * ff(mx,my,-nz+1,iv) &
!!                                      + 10._DP * ff(mx,my,-nz+2,iv) -  5._DP * ff(mx,my,-nz+3,iv) + ff(mx,my,-nz+4,iv)
!                   ff(mx,my,-nz-2,iv) =  5._DP * ff(mx,my,-nz-1,iv) - 10._DP * ff(mx,my,-nz  ,iv) &
!                                      + 10._DP * ff(mx,my,-nz+1,iv) -  5._DP * ff(mx,my,-nz+2,iv) + ff(mx,my,-nz+3,iv)
!                   ff(mx,my,-nz-3,iv) =  5._DP * ff(mx,my,-nz-2,iv) - 10._DP * ff(mx,my,-nz-1,iv) &
!                                      + 10._DP * ff(mx,my,-nz  ,iv) -  5._DP * ff(mx,my,-nz+1,iv) + ff(mx,my,-nz+2,iv)



                end do
              end do
          end do
!$OMP end do nowait

      else if ( rankz == nprocz-1 ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do iz = 1, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,nz+iz,iv) = zb2_top(mx,my,nzb-iz,2*nv+1-iv)
                end do
              end do
            end do
          end do
!$OMP end do nowait

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
! This part demands nzb=3
                  ff(mx,my,nz,iv) = ( 4._DP * ( ff(mx,my,nz-1,iv) + ff(mx,my,nz+1,iv) ) &
                                            - ( ff(mx,my,nz-2,iv) + ff(mx,my,nz+2,iv) ) ) / 6._DP
                end do
              end do
          end do
!$OMP end do nowait

      end if


  END SUBROUTINE iono_bndry_zvswap_buffout


!--------------------------------------
  SUBROUTINE micouple_bndry_bound_e ( psi )
!--------------------------------------
! This routine should be called after calling bndry_bound_e

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb) :: psi

    complex(kind=DP) :: dpdz

    integer :: mx, my, iz


      if( rankz == 0 ) then

!$OMP do collapse(2) schedule(dynamic)
              do my = ist_y, iend_y
                do mx = -nx, nx
! 1st order extrapolation -- this method provide the best results in growth rate

!... org

!                  psi(mx,my,-nz-1) =   psi(mx,my,-nz  )
!                  psi(mx,my,-nz-2) =   psi(mx,my,-nz  ) - ( psi(mx,my,-nz+1) - psi(mx,my,-nz  ) )
!                  psi(mx,my,-nz-3) =   psi(mx,my,-nz  ) - ( psi(mx,my,-nz+1) - psi(mx,my,-nz  ) )*2._DP

!... modified: 2nd order
!                  psi(mx,my,-nz-1) = 3._DP * ( psi(mx,my,-nz) - psi(mx,my,-nz+1) ) + psi(mx,my,-nz+2)
!                  psi(mx,my,-nz-2) = 3._DP * ( psi(mx,my,-nz-1) - psi(mx,my,-nz) ) + psi(mx,my,-nz+1)
!                  psi(mx,my,-nz-3) = 3._DP * ( psi(mx,my,-nz-2) - psi(mx,my,-nz-1) ) + psi(mx,my,-nz)

!                  psi(mx,my,-nz-1) = psi(mx,my,-nz)
!                  psi(mx,my,-nz-2) = 3._DP * ( psi(mx,my,-nz-1) - psi(mx,my,-nz) ) + psi(mx,my,-nz+1) 
!                  psi(mx,my,-nz-3) = 3._DP * ( psi(mx,my,-nz-2) - psi(mx,my,-nz-1) ) + psi(mx,my,-nz) 


!                   psi(mx,my,-nz-1) = 2.5_DP * psi(mx,my,-nz) - 2._DP * psi(mx,my,-nz+1) + 0.5_DP * psi(mx,my,-nz+2)
!                   psi(mx,my,-nz-2) = 2.5_DP * psi(mx,my,-nz-1) - 2._DP * psi(mx,my,-nz) + 0.5_DP * psi(mx,my,-nz+1)
!                   psi(mx,my,-nz-3) = 2.5_DP * psi(mx,my,-nz-2) - 2._DP * psi(mx,my,-nz-1) + 0.5_DP * psi(mx,my,-nz)
!... 4-point Lagrange
                   psi(mx,my,-nz-1) = 4._DP * psi(mx,my,-nz  ) - 6._DP * psi(mx,my,-nz+1) &
                                    + 4._DP * psi(mx,my,-nz+2) - 1._DP * psi(mx,my,-nz+3)
                   psi(mx,my,-nz-2) = 4._DP * psi(mx,my,-nz-1) - 6._DP * psi(mx,my,-nz  ) &
                                    + 4._DP * psi(mx,my,-nz+1) - 1._DP * psi(mx,my,-nz+2)
                   psi(mx,my,-nz-3) = 4._DP * psi(mx,my,-nz-2) - 6._DP * psi(mx,my,-nz-1) &
                                    + 4._DP * psi(mx,my,-nz  ) - 1._DP * psi(mx,my,-nz+1)

!... 5-point Lagrange
!                   psi(mx,my,-nz-1) =  5._DP * psi(mx,my,-nz  ) - 10._DP * psi(mx,my,-nz+1) &
!                                    + 10._DP * psi(mx,my,-nz+2) -  5._DP * psi(mx,my,-nz+3) + psi(mx,my,-nz+4)
!                   psi(mx,my,-nz-2) =  5._DP * psi(mx,my,-nz-1) - 10._DP * psi(mx,my,-nz  ) &
!                                    + 10._DP * psi(mx,my,-nz+1) -  5._DP * psi(mx,my,-nz+2) + psi(mx,my,-nz+3)
!                   psi(mx,my,-nz-3) =  5._DP * psi(mx,my,-nz-2) - 10._DP * psi(mx,my,-nz-1) &
!                                    + 10._DP * psi(mx,my,-nz  ) -  5._DP * psi(mx,my,-nz+1) + psi(mx,my,-nz+2)





!! 1st order extrapolation
!                  psi(mx,my,-nz-1) =   psi(mx,my,-nz  ) - ( psi(mx,my,-nz+1) - psi(mx,my,-nz  ) )
!                  psi(mx,my,-nz-2) =   psi(mx,my,-nz  ) - ( psi(mx,my,-nz+1) - psi(mx,my,-nz  ) )*2._DP
!                  psi(mx,my,-nz-3) =   psi(mx,my,-nz  ) - ( psi(mx,my,-nz+1) - psi(mx,my,-nz  ) )*3._DP
!! 2nd order extrapolation
!                  dpdz = ( -3._DP*psi(mx,my,-nz  ) + 4._DP*psi(mx,my,-nz+1) - psi(mx,my,-nz+2) ) * 0.5_DP
!                  psi(mx,my,-nz-1) =   psi(mx,my,-nz  ) - dpdz
!                  psi(mx,my,-nz-2) =   psi(mx,my,-nz  ) - dpdz * 2._DP
!                  psi(mx,my,-nz-3) =   psi(mx,my,-nz  ) - dpdz * 3._DP
!! 2nd order extrapolation
!                  psi(mx,my,-nz-1) = 3._DP * psi(mx,my,-nz  ) &
!                                   - 3._DP * psi(mx,my,-nz+1) &
!                                   + 1._DP * psi(mx,my,-nz+2)
!                  psi(mx,my,-nz-2) = 6._DP * psi(mx,my,-nz  ) &
!                                   - 8._DP * psi(mx,my,-nz+1) &
!                                   + 3._DP * psi(mx,my,-nz+2)
!! 4th order extrapolation  <=  Does NOT work
!                  psi(mx,my,-nz-1) =  5._DP * psi(mx,my,-nz  ) &
!                                   - 10._DP * psi(mx,my,-nz+1) &
!                                   + 10._DP * psi(mx,my,-nz+2) &
!                                   -  5._DP * psi(mx,my,-nz+3) &
!                                   +  1._DP * psi(mx,my,-nz+4)
!                  psi(mx,my,-nz-2) = 15._DP * psi(mx,my,-nz  ) &
!                                   - 40._DP * psi(mx,my,-nz+1) &
!                                   + 45._DP * psi(mx,my,-nz+2) &
!                                   - 24._DP * psi(mx,my,-nz+3) &
!                                   +  5._DP * psi(mx,my,-nz+4)
                end do
              end do
!$OMP end do nowait

      else if ( rankz == nprocz-1 ) then

!$OMP do collapse(2) schedule(dynamic)
            do iz = 1, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  psi(mx,my,nz+iz) = psi(mx,my,nz-iz)
                end do
              end do
            end do
!$OMP end do nowait

!$OMP do collapse(2) schedule(dynamic)
              do my = ist_y, iend_y
                do mx = -nx, nx
!! 1st order
!                  psi(mx,my,nz) = ( 2._DP * ( psi(mx,my,nz-1) + psi(mx,my,nz+1) ) &
!                                          - ( psi(mx,my,nz-2) + psi(mx,my,nz+2) ) ) / 2._DP
! 2nd order => still noisy
                  psi(mx,my,nz) = ( 4._DP * ( psi(mx,my,nz-1) + psi(mx,my,nz+1) ) &
                                          - ( psi(mx,my,nz-2) + psi(mx,my,nz+2) ) ) / 6._DP
                end do
              end do
!$OMP end do nowait

      end if


  END SUBROUTINE micouple_bndry_bound_e


!--------------------------------------
  SUBROUTINE iono_symmetry_ff ( ff_iono, ff_odd, ff_evn )
!--------------------------------------

   complex(kind=DP), intent(inout), &
     dimension(-nx:nx,0:ny,1:2*nv) :: ff_iono, ff_odd, ff_evn

   complex(kind=DP), dimension(:,:,:), allocatable :: zb1, zb2


     allocate( zb1(-nx:nx,0:ny,1:2*nv) )
     allocate( zb2(-nx:nx,0:ny,1:2*nv) )

       call iono_symmetry_ff_buffin( ff_iono, zb1 )

       call iono_symmetry_ff_sendrecv( zb1, zb2 )

       call iono_symmetry_ff_buffout ( ff_iono, zb2, ff_odd, ff_evn )

     deallocate( zb1 )
     deallocate( zb2 )


  END SUBROUTINE iono_symmetry_ff


!--------------------------------------
  SUBROUTINE iono_symmetry_ff_buffin( ff_iono, zb1 )
!--------------------------------------
! This routine should be called for imposing the boundary
! value at hh

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv) :: ff_iono
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,1:2*nv) :: zb1

    integer :: iz, iv

      zb1(:,:,:) = ( 0._DP, 0._DP )

      if( rankz /= 0 ) return

!---$OMP do collapse(2) schedule(dynamic)
        do iv = 1, 2*nv
          zb1(:,:,iv) = ff_iono(:,:,iv)
        end do
!---$OMP end do nowait


!---$OMP end master


  END SUBROUTINE iono_symmetry_ff_buffin


!--------------------------------------
  SUBROUTINE iono_symmetry_ff_sendrecv( zb1, zb2 )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv) :: zb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,1:2*nv) :: zb2

    integer :: slngz, iswp
    integer, dimension(2) :: ireq
    integer, dimension(MPI_STATUS_SIZE,2) :: istatus

      zb2(:,:,:) = ( 0._DP, 0._DP )

      if( rankz /= 0 ) return

      slngz  = (2*nx+1)*(ny+1)*(2*nv)

      if( rankz == 0 ) then
        iswp   = rank + nprocw * nprocz * ( nprocv - 1 - 2*rankv )
      else
        iswp   = MPI_PROC_NULL
      end if

!!debug
!        write( olog, * ) "# rank, rankv, iswp = ", rank, rankv, iswp
!        flush(olog)
!!debug

      call MPI_irecv( zb2, slngz, MPI_DOUBLE_COMPLEX, iswp, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_isend( zb1, slngz, MPI_DOUBLE_COMPLEX, iswp, 1, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_waitall( 2, ireq, istatus, ierr_mpi )

!debug
!      zb2_top(:,:,:) = zb1_top(:,:,:)
!debug


  END SUBROUTINE iono_symmetry_ff_sendrecv


!--------------------------------------
  SUBROUTINE iono_symmetry_ff_buffout ( ff_iono, zb2, ff_odd, ff_evn )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv) :: ff_iono, zb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,1:2*nv) :: ff_odd, ff_evn

    integer :: mx, my, iv


      if( rankz /= 0 ) return

! impose the symmetry in v
!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff_evn(mx,my,iv) = ( ff_iono(mx,my,iv) + zb2(mx,my,2*nv+1-iv) ) * 0.5_DP
                ff_odd(mx,my,iv) = ( ff_iono(mx,my,iv) - zb2(mx,my,2*nv+1-iv) ) * 0.5_DP
              end do
            end do
          end do
!$OMP end do nowait


  END SUBROUTINE iono_symmetry_ff_buffout


!--------------------------------------
  SUBROUTINE iono_fileio_open_icnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(3)   :: cold

    write( crank, fmt="(i6.6)" ) rankg
    write( cold,  fmt="(i3.3)" ) inum-1

    open( icnt_iono, file=path//crank//".cnt_iono."//cold, &
          form="unformatted", status="old", action="read" )

  END SUBROUTINE iono_fileio_open_icnt


!--------------------------------------
  SUBROUTINE iono_fileio_close_icnt
!--------------------------------------

     close( icnt_iono )

  END SUBROUTINE iono_fileio_close_icnt


!--------------------------------------
  SUBROUTINE iono_fileio_open_cnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(3)   :: cnew

    write( crank, fmt="(i6.6)" ) rankg
    write( cnew,  fmt="(i3.3)" ) inum

    open( ocnt_iono, file=path//crank//".cnt_iono."//cnew, &
          form="unformatted" )

  END SUBROUTINE iono_fileio_open_cnt

!--------------------------------------
  SUBROUTINE iono_fileio_close_cnt
!--------------------------------------

     close( ocnt_iono )

  END SUBROUTINE iono_fileio_close_cnt


!--------------------------------------
  SUBROUTINE iono_fileio_open_phi ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( ophi_iono, file=path//crank//"."//srank//".phi_iono."//cnew, &
          form="unformatted" )

  END SUBROUTINE iono_fileio_open_phi

!--------------------------------------
  SUBROUTINE iono_fileio_close_phi
!--------------------------------------

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    close( ophi_iono )

  END SUBROUTINE iono_fileio_close_phi


!--------------------------------------
  SUBROUTINE iono_fileio_read_cnt ( wf, time, istatus )
!--------------------------------------

    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny) :: wf
    real(kind=DP), intent(out) :: time
    integer, optional, intent(out) :: istatus

    integer :: input_status

    read( unit=icnt_iono, iostat=input_status ) time, wf
    if ( present(istatus) ) then
       istatus = input_status
    endif

  END SUBROUTINE iono_fileio_read_cnt


!--------------------------------------
  SUBROUTINE iono_fileio_write_cnt ( wf, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: wf
    real(kind=DP), intent(in) :: time

!    if ( rankz == 0 .and. ranks == 0 .and. vel_rank == 0 ) then
    if ( rankz == 0 ) then
      rewind ocnt_iono
      write( unit=ocnt_iono ) time, wf

      call flush(ocnt_iono)
    end if

  END SUBROUTINE iono_fileio_write_cnt


!--------------------------------------
  SUBROUTINE iono_fileio_write_phi ( idns, iphi, icpr, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: iphi, icpr, idns
    real(kind=DP), intent(in) :: time

    !- OUTPUT binary data phi/"*.phi.*"
    if ( rankz == 0  .and.  ranks == 0 .AND. vel_rank == 0 ) then
       write( unit=ophi_iono ) time, idns, iphi, icpr

       call flush(ophi_iono)

!! debug
!       write(olog,*) "# time = ", time
!       write(olog,fmt="(1p,6e15.7)") iphi(0,ny/2), icpr(0,ny/2), idns(0,ny/2)
!! debug

    end if


  END SUBROUTINE iono_fileio_write_phi


!--------------------------------------
  SUBROUTINE micouple_out_cntrl ( time, id )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    real(kind=DP), save :: tout_iono


      if ( rankz /= 0 ) return

      if( id == 0 ) then

        if ( time == 0._DP ) then
          call iono_wrt ( idns, iphi, icpr, time, id )
        end if

        tout_iono  = ( int( ( time + eps )/dtout_eng ) + 1 ) * dtout_eng
 
      else if( id == 1 ) then

        if ( time >= tout_iono - eps ) then
          call iono_wrt ( idns, iphi, icpr, time, id )
          tout_iono   = tout_iono + dtout_eng
        end if

      else if( id == 2 ) then

          call iono_wrt ( idns, iphi, idns, time, id )

      end if


  END SUBROUTINE micouple_out_cntrl


!--------------------------------------
  SUBROUTINE iono_wrt ( idns, iphi, icpr, time, id )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: idns, iphi, icpr
    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    real(kind=DP), dimension(-nx:nx) :: idns_mode_x, iphi_mode_x, icpr_mode_x
    real(kind=DP), dimension(0:gny) :: idns_mode_y, iphi_mode_y, icpr_mode_y
    real(kind=DP), dimension(0:nk) :: idns_mode_k, iphi_mode_k, icpr_mode_k

    real(kind=DP) :: idns_total_e, icpr_total_e, iphi_total_e
    

      if( id == 2 ) then
        call iono_fileio_write_cnt ( idns, time )

        if ( rankg == 0 ) then
          call flush( oeng_idns )
          call flush( oeng_iphi )
          call flush( oeng_icpr )
          call flush( obln_iono )
        end if

        return
      end if

        call imode_energy ( idns, idns_mode_x, idns_mode_y, idns_mode_k, idns_total_e )
        call imode_energy ( iphi, iphi_mode_x, iphi_mode_y, iphi_mode_k, iphi_total_e )
        call imode_energy ( icpr, icpr_mode_x, icpr_mode_y, icpr_mode_k, icpr_total_e )

      if ( rankg == 0 ) then
        write( unit=oeng_idns, fmt="(f15.8, 2050ES24.15e3)" ) &
               time, idns_total_e, idns_mode_y(0:global_ny)
        write( unit=oeng_iphi, fmt="(f15.8, 2050ES24.15e3)" ) &
               time, iphi_total_e, iphi_mode_y(0:global_ny)
        write( unit=oeng_icpr, fmt="(f15.8, 2050ES24.15e3)" ) &
               time, icpr_total_e, icpr_mode_y(0:global_ny)

      end if
      call iono_fileio_write_phi ( idns, iphi, icpr, time )
      call iono_balance ( idns, iphi, icpr, time )

  END SUBROUTINE iono_wrt

!--------------------------------------
  SUBROUTINE iono_balance ( idns, iphi, icpr, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: idns, iphi, icpr
    real(kind=DP), intent(in) :: time

    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: ipsnb, idvgd

    real(kind=DP) :: ibln_tmp(4), ibln(4)
    integer :: mx, my, iz

!... ionosphere energy balance
      ibln_tmp(:) = 0._DP 

      if( trim(calc_type) == "nonlinear" ) then
        call pssn_brackets( iphi, idns, ipsnb, 1 )
        call pssn_divgrad ( idns, iphi, idvgd, 1 )
      else
        ipsnb(:,:) = ( 0._DP, 0._DP )
        idvgd(:,:) = ( 0._DP, 0._DP )
      end if

      iz = -nz

      do my = ist1_y, iend_y
         do mx = -nx, nx

            ibln_tmp(1) = ibln_tmp(1) + real( beta * icpr(mx,my) * conjg(iphi(mx,my)) )
            ibln_tmp(2) = ibln_tmp(2) + real( - mp * idns0 * ksq(mx,my,iz) * iphi(mx,my) * conjg(iphi(mx,my)) ) &
                                      + real( ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) ) &
                                              ) * idns(mx,my) * conjg(iphi(mx,my))   &
                                            ) 
            ibln_tmp(3) = ibln_tmp(3) + real( ( - ksq(mx,my,iz) * dperp  &
                                              ) * idns(mx,my) * conjg(iphi(mx,my))   &
                                            )
            ibln_tmp(4) = ibln_tmp(4) + real( ( mp * idvgd(mx,my) + ipsnb(mx,my) ) * conjg(iphi(mx,my)) )

         end do
      end do
  
      if( rankw == 0 ) then
        my = 0
         do mx = 1, nx

            ibln_tmp(1) = ibln_tmp(1) + real( beta * icpr(mx,my) * conjg(iphi(mx,my)) )
            ibln_tmp(2) = ibln_tmp(2) + real( - mp * idns0 * ksq(mx,my,iz) * iphi(mx,my) * conjg(iphi(mx,my)) ) &
                                      + real( ( -ui * e0 * ( mp*ky(my) - mh*kx(mx) ) &
                                              ) * idns(mx,my) * conjg(iphi(mx,my))   &
                                            ) 
            ibln_tmp(3) = ibln_tmp(3) + real( ( - ksq(mx,my,iz) * dperp  &
                                              ) * idns(mx,my) * conjg(iphi(mx,my))   &
                                            )
            ibln_tmp(4) = ibln_tmp(4) + real( ( mp * idvgd(mx,my) + ipsnb(mx,my) ) * conjg(iphi(mx,my)) )

         end do
      end if

      call MPI_Reduce( ibln_tmp, ibln, 4, MPI_REAL8, &
                       MPI_SUM, 0, fft_comm_world, ierr_mpi )

      if ( rankg == 0 ) then

        ibln(:) = 2._DP * ibln(:) / beta 
        ibln(1) = -ibln(1) 
        ibln(2) = -ibln(2)
        ibln(4) = -ibln(4)

        write( unit=obln_iono, fmt="(f15.8, 2050ES24.15e3)" ) &
               time,  ibln(1),  & ! Poynting flux
                      ibln(2),  & ! Joule dissipation from linear terms
                      ibln(3),  & ! Collisional diffusion
                      ibln(4)     ! Nonlinear terms ( only divgrad term contributes )

                                  ! JE (1) = jE (2) - Dperp-term (3) + NL (4)


      end if

    END SUBROUTINE iono_balance


!--------------------------------------
  SUBROUTINE imode_energy ( wrk, mode_x, mode_y, mode_k, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny) :: wrk

    real(kind=DP), dimension(08:nx) :: mode_x
    real(kind=DP), dimension(0:nk)  :: mode_k
    real(kind=DP), dimension(0:gny) :: mode_y

    real(kind=DP) :: total_e

! --- local variables

! THW 2023.03.18
    real(kind=DP), dimension(0:nx)        :: mode_xw
    real(kind=DP), dimension(0:ny)        :: mode_yw
! THW 2023.03.18

    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  mx, my

    integer, dimension(0:nk) :: ic

    real(kind=DP) :: dk
    integer       :: ikp


      allocate( wr2(-nx:nx,0:ny) )

      mode_x(:)  = 0._DP
      mode_y(:)  = 0._DP
! THW 2023.03.18
      mode_xw(:)  = 0._DP
      mode_yw(:)  = 0._DP
! THW 2023.03.18
      total_e    = 0._DP

      mode_k(:)  = 0._DP
      ic(:)      = 0

      do my = ist_y, iend_y
        do mx = -nx, nx
          wr2(mx,my) = real( wrk(mx,my) * conjg( wrk(mx,my) )  &
                                                  , kind=DP )
        end do
      end do

! --- ky-modes
      do my = ist1_y, iend_y
         do mx = -nx, nx
            mode_yw(my) = mode_yw(my) + wr2(mx,my) 
         end do
      end do

     if( rankw == 0 ) then
       my = 0
       do mx = 1, nx
         mode_yw(my) = mode_yw(my) + wr2(mx,my)
       end do
     end if

     call MPI_Gather( mode_yw,   ny+1, MPI_DOUBLE_PRECISION, &
                      mode_y,  ny+1, MPI_DOUBLE_PRECISION, &
                      0, fft_comm_world, ierr_mpi )

! --- total in whole k-space
      do my = 0, global_ny
        total_e = total_e + mode_y(my)
      end do

! --- kx-modes
      do my = ist1_y, iend_y
        do mx = 0, nx
          mode_xw(mx) = mode_xw(mx) + wr2(mx,my)
        end do
      end do
  
      if( rankw == 0 ) then
        my = 0
        do mx = 1, nx
          mode_xw(mx) = mode_xw(mx) + wr2(mx,my)
        end do
      end if

      call MPI_Allreduce( mode_xw, mode_x, nx+1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
  


!! --- k_perp-modes
!!      dk = kx(1)
!      dk = ky(1)
!      do my = 1, ny
!        do mx = 1, nx
!          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
!          if( ikp <= nk ) then
!            mode_k(ikp) = mode_k(ikp) + wr2(mx,my) + wr2(-mx,my)
!            ic(ikp)  = ic (ikp) + 2
!          end if
!        end do
!      end do

!      mx = 0
!        do my = 1, ny
!          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
!          if( ikp <= nx ) then
!            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
!            ic(ikp)  = ic (ikp) + 1
!          end if
!        end do

!      my = 0
!        do mx = 1, nx
!          ikp = nint( sqrt( kx(mx)**2 + ky(my)**2 ) / dk )
!          if( ikp <= nx ) then
!            mode_k(ikp) = mode_k(ikp) + wr2(mx,my)
!            ic(ikp)  = ic (ikp) + 1
!          end if
!        end do

!        do ikp = 1, nk
!            mode_k(ikp) = mode_k(ikp) / real( ic(ikp), kind=DP )
!        end do



      deallocate( wr2 )


  END SUBROUTINE imode_energy

!--------------------------------------
   SUBROUTINE find_roots( id, mxin, myi, omg_out, kz_out )
!--------------------------------------
!
!  Note:
!  - Calculation is performed in the MHD unit. 
!  - L_MHD, vA, etc are unity in the MHD normalization
!
   complex(kind=DP), intent(out) :: omg_out, kz_out
   integer, intent(in) :: id 
   integer, intent(in) :: mxin, myi

   complex(kind=DP) :: sol, solmax
   complex(kind=DP), allocatable :: solk(:,:), kz(:,:)
   complex(kind=DP) :: dum0
   complex(kind=DP) :: z_guess(root_max)
   complex(kind=DP) :: sol_tmp, kz_tmp
   real(kind=DP) :: k_unit, kx_pass, ky_pass
   integer, allocatable :: rcnt(:)
   integer :: k, kst, kmax
   integer :: cnt0, root, root_unst, sflg
   integer :: ierr
   character(3) :: cxmd 

   k_unit  = kymin_g / rho2L ! Normalize k by L, not by rho
   z_guess = 0._DP

   if( id == 0 ) then ! find roots for dipersion relation plot

     if( rankg==0 ) then
       if ( mxin >= 0 ) then
         write( cxmd,  fmt="(i3.3)" ) mxin
         open( oroot, file=trim(f_hst)//"dsp_root-p"//cxmd )
       else
         write( cxmd,  fmt="(i3.3)" ) abs(mxin)
         open( oroot, file=trim(f_hst)//"dsp_root-n"//cxmd )
       end if
     end if

     kst     = 1
     kmax    = min( 100, global_ny )

   else if( id == 1 ) then ! find roots of given mxi and myi for initial condition

     if( rankg==0 ) then      
       open( oroot, file=trim(f_hst)//"eigen_mode" )
     end if
     kst     = myi
     kmax    = myi

   end if
   kx_pass = k_unit * mxin
   ky_pass = k_unit

   allocate( solk(kmax,root_max), kz(kmax,root_max), rcnt(kmax) )
   solk(:,:) = ( 0._DP, 0._DP ) ! solutions of omega for given kperp
   kz(:,:)   = ( 0._DP, 0._DP ) ! solutions of kpara for given kperp

!... find roots for ky(1) 
   rcnt(:) = 0
   cnt0    = 0

   call guess_roots( z_guess, kx_pass, ky_pass, cnt0 ) ! get z_guess

   cnt0 = min( cnt0, root_max ) 

! put z_guess and find acutual roots

   do root=1, cnt0
     call solv( z_guess(root), kx_pass, ky_pass, sol, kz_tmp, dum0, sflg )
     if ( sflg == 1 ) then ! root found
        rcnt(kst) = rcnt(kst) + 1
        solk(kst,rcnt(kst)) = sol 
        kz(kst,rcnt(kst))   = kz_tmp
     end if
   end do

   call reorder_roots( solk(kst,:), kz(kst,:), rcnt(kst), cnt0 )
   rcnt(kst)=cnt0

   do k=kst+1, kmax

     ky_pass = k_unit * dble(k)
     rcnt(k) = 0

     do root=1, rcnt(kst)

       z_guess = solk(k-1,root) ! guess z from roots for ky(1)
          
       call solv( z_guess(root), kx_pass, ky_pass, sol, kz_tmp, dum0, sflg )         
       if ( sflg == 1 ) then
         rcnt(k) = rcnt(k) + 1
         solk(k,rcnt(k)) = sol
         kz(k,rcnt(k))   = kz_tmp 
       end if

     end do

     if ( rcnt(k)>1 ) then
       call reorder_roots( solk(k,:), kz(k,:), rcnt(k), cnt0 )
      rcnt(k) = cnt0
     end if

   end do

   ! convert the unit from v_A/L_I to v_A/l 

   solk(:,:) = solk(:,:) * zl
   kz(:,:)   = kz(:,:)   * zl

   ! select physical roots and write out results
   if( rankg==0 ) then

     solmax = ( 0._DP, 0._DP )
     root_unst = 0
   
     write( oroot,* ) " # mx, my, omg_r [vA/l], gamma [vA/l], Re(omg_r*l/pi) [vA], Im(gamma*l/pi) [vA], Re(kz) [1/l], Im(kz) [1/l], -gamma/omg_r"
     do root=1, 10
  
       write( oroot,* ) " # root = ", root 
       do k=kst, kmax

         if ( rcnt(k)<root ) cycle

         if ( real(solk(k,root)) < 0._DP ) then
           write( oroot,'(a6, 2i6,8ES16.6)' ) &
           " # ", nint( kx_pass/twopi ), k, solk(k,root), solk(k,root)/pi, kz(k,root), -aimag(solk(k,root))/real(solk(k,root))
         else
           write( oroot,'(2i6,8ES16.6)' )     &
                  nint( kx_pass/twopi ), k, solk(k,root), solk(k,root)/pi, kz(k,root), -aimag(solk(k,root))/real(solk(k,root))

           if( aimag( solk(k,root) - solmax ) > 0._DP ) then
             root_unst = root
             solmax    = solk(k,root)
           end if

         end if

       end do

       write( oroot,* )
       write( oroot,* )

     end do

     write( oroot,* ) " # most unstable root number is ", root_unst

     flush( oroot )
     close( oroot )
  
   end if

   if( id==1 ) then ! return most unstable root after converting to values in gkv unit
    
     call MPI_Bcast( root_unst,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 

     omg_out = solk(myi,root_unst) * valf / twopi 
     kz_out  = kz(myi,root_unst) / twopi

   end if

   return

   END SUBROUTINE find_roots

!----------------------------------------------
!     solver in Newton method
!     
!     Calculate frequency, omega (named z here)
!     and parallel wavenumber, kz 
!
!----------------------------------------------
  SUBROUTINE solv( zom0, kx, ky, zom, kzout, egnf, sflg )

      complex(kind=DP), intent(in) ::  zom0
      complex(kind=DP), intent(out) :: kzout
      real(kind=DP)    :: kx, ky
      complex(kind=DP) :: zom, kz, egnf

      complex(kind=DP) :: f, df, z, znew

      real(kind=DP)    :: eps, eps0

      integer      isw, i, sflg
 
!!      data     eps0     /  1.d-7 /
      sflg =1
      eps0 = 1.d-6

        z          = zom0
        znew       = zom0

        i          = 0

   10   continue
          i          = i + 1

          call  fun ( z,  f, kx, ky, kz, egnf )
          call dfunr( z, df, kx, ky, egnf )

          znew       = z - f / df
          eps        = abs( znew - z ) / abs( z )
          z          = znew
          if( i .ge. 30 ) then
             sflg = 0
             goto 999
          end if

         if( eps .gt. eps0 )   goto 10

  999   continue

        zom = z
        kzout = kz

      return

  END SUBROUTINE solv

!----------------------------------------------
  SUBROUTINE solv_kz( z, kx, ky, kzout )
!----------------------------------------------
!
! Give the frequency omega ( named z here )
! and get consistent kz 
!
      real(kind=DP)    :: kx, ky
      complex(kind=DP) :: egnf

      complex(kind=DP) :: f, df, z, kz, kznew, kzout

      real(kind=DP)    :: eps, eps0

      integer      isw, i
!!      data     eps0     /  1.d-7 /
      eps0 = 1.d-6

        kz     = z
        kznew  = z

        i      = 0

   10   continue
          i          = i + 1

          call  fun_kz( z,  f, kx, ky, kz )
          call dfun_kz( z, df, kx, ky, kz )

          kznew  = kz - f / df
          eps    = abs( kznew - kz ) / abs( kz )
          kz     = kznew
          if( i .ge. 30 ) then
             goto 999
          end if

         if( eps .gt. eps0 )   goto 10

  999   continue

        kzout = kz

      return

  END SUBROUTINE solv_kz

!----------------------------------------------
!     functions to be solved
!
!----------------------------------------------
  SUBROUTINE fun( z, f, kx, ky, kzout, egnf )

    use GKV_math, only: math_g0
      complex(kind=DP), intent(in)  :: z
      complex(kind=DP), intent(out) :: f, kzout
      complex(kind=DP) :: kz, chiz, zres, cothZ
      real(kind=DP)    :: facG, bs, Gam0i, KME0

      real(kind=DP), intent(in) :: kx, ky
      complex(kind=DP) :: egnf

      real(kind=DP)    :: kperp

        kperp =  sqrt( kx**2 + ky**2 ) 
        KME0  =  e0 * ( mp * ky - mh * kx )

        if ( calc_disp > 1 ) then ! MHD 

           kz = z
           chiz = 1._DP

        else ! FLR effect is included

           bs    = ( kperp * rho2L )**2
           call math_g0( bs, Gam0i )
           facG = 1._DP - Gam0i 
           call solv_kz ( z, kx, ky, kz )
           chiz = ( bs / facG ) * ( kz / z )  

        end if

        cothZ =  cosh( ui * kz * zl ) / sinh( ui * kz * zl )
        zres = -mp * idns0 * chiz * cothZ                          ! magnetosphere impedance
        f   = ( z + ui * 2.d+0 * alpha ) * ( 1._DP + zres ) &
               + ui * kperp**2 * dperp &
               - KME0 

        kzout = kz
!
!   MHD case: 
!   cotZ =  cosh( ui*z*zl ) / sinh( ui*z*zl )
!   res  = -5.0d0 * cotZ!! 
!! 
     return

  END SUBROUTINE fun

!===============================
  SUBROUTINE fun_kz( z, f, kx, ky, kz )
!==============================
!
!   Eight-pole Pade approximation is used
!   for the plasme dispersion function, funcZ
!
    use GKV_math, only: math_g0

      complex(kind=DP), intent(in)  :: z
      complex(kind=DP), intent(out) :: f
      complex(kind=DP) :: kz, zeta
      complex(kind=DP) :: funcZ, facZ
      real(kind=DP), intent(in) :: kx, ky
      real(kind=DP) :: kperp
      real(kind=DP) :: Gam0i, facG, bs, vte
      integer :: i

        kperp =  sqrt( kx**2 + ky**2 ) ! ( k * L_I )**2
        bs    = ( kperp * rho2L )**2   ! ( k * rho )**2
        call math_g0( bs, Gam0i )

        facG = 1._DP - Gam0i

        vte   = sqrt( beta / Anum(0) ) ! electron thermal velocity relative to vA 
        zeta  = z / kz / vte / sqrt( 2._DP )

        funcZ = ( 0._DP, 0._DP ) 
        do i=1, 8
          funcZ = funcZ + coefb(i) / ( zeta - coefc(i) )
        end do
        
        facZ  = 1._DP + zeta * funcZ  

! Gyrokinetic dispersion relation 
        f  = z**2 - kz**2 * bs * ( 1._DP / facZ + 1._DP / facG )

     return

  END SUBROUTINE fun_kz


!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfunr( z, df, kx, ky, egnf )

      complex(kind=DP) :: z, df, kz
      real(kind=DP)    :: kx, ky

      complex(kind=DP) :: egnf

      complex(kind=DP) :: ff1, ff2, z1, z2
      real(kind=DP)    :: dx, d

!!      data     dx  /  1.d-10 /
      dx = 1.d-8

        d          = dble(z) * dx
        z1         = z + dcmplx( d*0.5d0, 0.d0 ) ! Omega+dOmega/2
        z2         = z - dcmplx( d*0.5d0, 0.d0 ) ! Omega-dOmega/2

        call fun( z1, ff1, kx, ky, kz, egnf )
        call fun( z2, ff2, kx, ky, kz, egnf )

        df         = ( ff1 - ff2 ) / d

      return

  END SUBROUTINE dfunr

!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfun_kz( z, df, kx, ky, kz )

      complex(kind=DP), intent(in)  :: z, kz
      complex(kind=DP), intent(out) :: df
      real(kind=DP), intent(in)    :: kx, ky
  
      complex(kind=DP) :: ff1, ff2, kz1, kz2
      real(kind=DP)    :: dx, d
 
!      data     dx  /  1.d-10 /
      dx = 1.d-8

        d    = dble(kz) * dx
        kz1  = kz + dcmplx( d*0.5d0, 0.d0 ) ! kz+dkz/2
        kz2  = kz - dcmplx( d*0.5d0, 0.d0 ) ! kz-dkz/2

        call fun_kz( z, ff1, kx, ky, kz1 )
        call fun_kz( z, ff2, kx, ky, kz2 )

        df   = ( ff1 - ff2 ) / d

      return

  END SUBROUTINE dfun_kz

!--------------------------------------
   SUBROUTINE guess_roots( z_guess, kx_pass, ky_pass, rcnt )
!--------------------------------------

!  Search each vicinity of multiple roots of f=0

     implicit none
     complex(kind=DP), intent(out) :: z_guess(root_max)
     real(kind=DP),intent(in) :: kx_pass, ky_pass
     integer,intent(out) :: rcnt

     complex(kind=DP)  :: f, z
     complex(kind=DP) :: dum0, kz
     real(kind=DP) :: fr, fi
     real(kind=DP) :: z1, z2, dz1, dz2, z1max, z2min 
     integer :: sgnr0, sgni0, sgnr, sgni
     integer, parameter :: imax = 200
     integer :: i, j

     z1max = 8._DP * pi / zl    ! upper boundary of the real freq.
     z2min = - 0.01_DP * z1max  ! lower boundary of the growth rate
     dz1   = z1max / dble(imax) 
     dz2   = 0.01_DP * dz1 

     rcnt = 0
     z1 = 0.0d0

     do i=1,imax
        z1 = z1 + dz1
        if (z1 > z1max) exit
        z2 = z2min
        do j=1, 2*imax
           z2 = z2 + dz2
           z  = z1 + ui*z2 
           call fun( z, f, kx_pass, ky_pass, kz, dum0 )
           if ( f /= f ) cycle ! exclude f = NaN
           fr = real(f)
           fi = aimag(f)
           if( i==1 .and. j==1 ) then
              sgnr0 = int(fr/abs(fr))
              sgni0 = int(fi/abs(fi))
           else
              sgnr = int(fr/abs(fr))
              sgni = int(fi/abs(fi))
              if(sgnr==sgnr0 .and. sgni/=sgni0) then
                 cycle
              else if (sgnr/=sgnr0 .and. sgni/=sgni0) then
                 sgnr0 = sgnr
                 sgni0 = sgni
                 rcnt  = rcnt + 1
                 if(rcnt > root_max) exit
                 z_guess(rcnt)  = z
                 exit
              end if
           end if
        end do
     end do

     return

   END SUBROUTINE guess_roots

!--------------------------------------
   SUBROUTINE reorder_roots( solk, kz, rcnt, ucnt)
!--------------------------------------

!  Arrange the roots (sol) in order of Re(Omega) from smaller to larger
!  then remove duplicated roots and count the number of unique ones (ucnt).

   implicit none
   complex(kind=DP), intent(inout) :: solk(root_max), kz(root_max)
   complex(kind=DP) :: sol_unique(rcnt), kz_unique(rcnt)
   complex(kind=DP) :: sol_tmp, kz_tmp

   real(kind=DP) :: diff_lim
   integer, intent(in)  :: rcnt
   integer, intent(out) :: ucnt ! number of unique solution
   integer :: i, j, flg

   diff_lim = 0.01_DP * pi / zl ! if |A-B| < diff_lim, A and B are considered the same

      do i=1,rcnt-1
         do j=i+1,rcnt
  
            if( real(solk(i) - solk(j)) > 0.0d0 ) then  
                 
               sol_tmp   = solk(i)
               kz_tmp    = kz(i)  
    
               solk(i) = solk(j)
               solk(j) = sol_tmp

               kz(i)   = kz(j)
               kz(j)   = kz_tmp

            end if

         end do
      end do

! screen out non-unique solutions 
      ucnt=0
      do i=1,rcnt-1
         flg=0
         do j=i+1,rcnt
            if( abs(solk(i) - solk(j)) < diff_lim ) then 
               flg=1
            end if
         end do
         if(flg==0) then
            ucnt=ucnt+1
            sol_unique(ucnt) = solk(i)
            kz_unique(ucnt)  = kz(i)
         end if
      end do

      solk(:) = ( 0._DP, 0._DP ) 
      kz(:)   = ( 0._DP, 0._DP )
      do i=1,ucnt
         solk(i) = sol_unique(i)
         kz  (i) = kz_unique(i)
      end do

      return
   END SUBROUTINE reorder_roots



!--------------------------------------
  SUBROUTINE wrt_fxyz ( ff, time )
!--------------------------------------

!   fout: f on sampled (np) field lines 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP), intent(in) :: time
 
    real(kind=DP), allocatable, dimension(:,:,:,:) :: fout
    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
                      wc1o, wc2o, wc1e, wc2e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                      wwo, wwe
    real(kind=DP), allocatable, dimension(:,:,:)   :: wr1o, wr1e
    real(kind=DP), allocatable, dimension(:,:,:,:) :: wr2o, wr2e

    integer :: iv, im
! for debug
    integer :: iz


!$OMP master
                                           call clock_sta(1800)
!$OMP end master

    if( np <= 0 ) return

    allocate( fout (np,-nz:nz-1,1:2*nv,0:nm) )

    allocate(wc1o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2o(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc1e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wc2e(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1))
    allocate(wwo(0:global_ny,0:2*nxw-1,0:nbuff-1))
    allocate(wwe(0:global_ny,0:2*nxw-1,0:nbuff-1))
!    allocate(wr1o(0:2*nxw-1,0:2*nyw-1,0:nbuff-1))
!    allocate(wr1e(0:2*nxw-1,0:2*nyw-1,0:nbuff-1))
!    allocate(wr2o(0:2*nxw-1,0:2*nyw-1,0:nbuff-1,0:nprocw-1))
!    allocate(wr2e(0:2*nxw-1,0:2*nyw-1,0:nbuff-1,0:nprocw-1))
    allocate(wr1o(np,0:2*nyw-1,0:nbuff-1))
    allocate(wr1e(np,0:2*nyw-1,0:nbuff-1))
    allocate(wr2o(np,0:2*nyw-1,0:nbuff-1,0:nprocw-1))
    allocate(wr2e(np,0:2*nyw-1,0:nbuff-1,0:nprocw-1))


!$OMP parallel default(shared) &
!$OMP private(iv)
!$OMP workshare
      wc1o(:,:,:,:) = (0._DP, 0._DP)
      wc1e(:,:,:,:) = (0._DP, 0._DP)
!$OMP end workshare

!$OMP barrier
      do iv = 1, 2*nv+5
        if (mod(iv,2) == 1) then ! odd
          if (1+4<=iv .and. iv<=2*nv+4) call fout_transpose_back(wr1o,wr2o)       ! 5,7,9,...
          if (1+5<=iv .and. iv<=2*nv+5) call fout_reduction(iv-5,wr2e,fout)       ! 7,9,11,...
          if (1+1<=iv .and. iv<=2*nv+1) call fout_transpose_y2zm(wc1e,wc2e)       ! 3,5,7,...
          if (1  <=iv .and. iv<=2*nv  ) call fout_pack_ff_y2zm(iv,ff,wc1o)        ! 1,3,5,... 
          if (1+2<=iv .and. iv<=2*nv+2) call fout_unpack_y2zm(wc2o,wwo)           ! 3,5,7,... 
          if (1+3<=iv .and. iv<=2*nv+3) call fout_realspcal_y2zm(wwe,wr1e)        ! 5,7,9,...
        else                     ! even
          if (1+4<=iv .and. iv<=2*nv+4) call fout_transpose_back(wr1e,wr2e)       ! 6,8,10,...
          if (1+5<=iv .and. iv<=2*nv+5) call fout_reduction(iv-5,wr2o,fout)       ! 6,8,10,...
          if (1+1<=iv .and. iv<=2*nv+1) call fout_transpose_y2zm(wc1o,wc2o)       ! 2,4,6,...
          if (1  <=iv .and. iv<=2*nv  ) call fout_pack_ff_y2zm(iv,ff,wc1e)        ! 2,4,6,...
          if (1+2<=iv .and. iv<=2*nv+2) call fout_unpack_y2zm(wc2e,wwe)           ! 4,6,8,...
          if (1+3<=iv .and. iv<=2*nv+3) call fout_realspcal_y2zm(wwo,wr1o)        ! 4,6,8,...
        end if 
!$OMP barrier
      end do
!$OMP end parallel

    if(rankw==0) then
      write( ofzv ) time, fout 
   
      call flush(ofzv) 
      ! debug
      !do iz=-nz, nz-1
      !  do iv = 1, 2*nv
      !     write(odbg,'(2i6,99ES16.6)') (2*nz)*rankz+iz+nz, (2*nv)*rankv+iv, fout(6,iz,iv,0:nm) 
      !  end do
      !  write(odbg,*)
      !end do
      !write(odbg,*)
      !write(odbg,*)

    end if

    deallocate( fout )
    
    deallocate(wc1o)
    deallocate(wc2o)
    deallocate(wc1e)
    deallocate(wc2e)
    deallocate(wwo)
    deallocate(wwe)
    deallocate(wr1o)
    deallocate(wr2o)
    deallocate(wr1e)
    deallocate(wr2e)

!$OMP master
                                           call clock_end(1800)
!$OMP end master



  END SUBROUTINE wrt_fxyz

!--------------------------------------
  SUBROUTINE fout_pack_ff_y2zm ( iv, fin, wc4 )
!--------------------------------------
!     Data pack for E x B term calculation (y2zm)

    integer, intent(in) :: iv
    !complex(kind=DP), intent(in), &
    !  dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: fin
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: fin
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              wc4(mx,my,ibuff,iprocw) = fin(mx,my,iz,iv,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait

  END SUBROUTINE fout_pack_ff_y2zm

!--------------------------------------
  SUBROUTINE fout_transpose_y2zm ( wc4in, wc4out )
!--------------------------------------
!     Data transpose for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4in
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4out

!$OMP master
                                           call clock_sta(1820)
!$OMP end master

!$OMP master
  
      call MPI_Alltoall( wc4in,                 &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         wc4out,                &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )

!$OMP end master

!$OMP master
                                           call clock_end(1820)
!$OMP end master

  END SUBROUTINE fout_transpose_y2zm

!--------------------------------------
  SUBROUTINE fout_unpack_y2zm ( wc4, ww )
!--------------------------------------
!     Data unpack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
!    complex(kind=DP), intent(out), &
!      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: ww
    complex(kind=DP), intent(out), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: ww

    complex(kind=DP), dimension(-nx:nx) :: psi
    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
    integer :: my, ibuff, iprocw, global_my

!$OMP master
                                           call clock_sta(1810)
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
!          ww(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
          ww(global_my,0:2*nxw-1,ibuff) = w2(0:2*nxw-1)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait

!$OMP master
                                           call clock_end(1810)
!$OMP end master

  END SUBROUTINE fout_unpack_y2zm

!--------------------------------------
  SUBROUTINE fout_realspcal_y2zm ( ww, wout )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(0:global_ny,0:2*nxw-1,0:nbuff-1) :: ww
    real(kind=DP), intent(out), &
      dimension(np,0:2*nyw-1,0:nbuff-1) :: wout
!      dimension(0:2*nxw-1,0:2*nyw-1,0:nbuff-1) :: wout

    complex(kind=DP), dimension(0:nyw) :: w3
    real(kind=DP), dimension(0:2*nyw-1) :: ftmp
    real(kind=DP) :: cef
    integer :: ix, iy, ip, ibuff

!$OMP master
                                           call clock_sta(1811)
!$OMP end master

!      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)

!$OMP do collapse(2) schedule(dynamic,nchunk_xb)
      do ibuff = 0, nbuff-1
!        do ix = 0, 2*nxw-1
         do ip =1, np 
         !%%% Backward y-FFT (ky,x)->(y,x) %%%
!          w3(0:global_ny) = ww(0:global_ny,ix,ibuff)
          w3(0:global_ny) = ww(0:global_ny,xp(ip),ibuff)
          w3(global_ny+1:nyw) = (0._DP, 0._DP) ! FFTW may destroy input array!
          call dfftw_execute_dft_c2r(plan_y_backward, w3, ftmp)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !%%% f in (x,y) %%%
          do iy = 0, 2*nyw-1
            !!fxy(ix,iy,iv,ibuff) = cef * ftmp(iy) 
            wout(ip,iy,ibuff) = ftmp(iy) 
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!$OMP end do nowait

!$OMP master
                                           call clock_end(1811)
!$OMP end master


  END SUBROUTINE fout_realspcal_y2zm


!--------------------------------------
  SUBROUTINE fout_transpose_back( wrin, wrout )
!--------------------------------------

!    real(kind=DP), intent(in), &
!      dimension(0:2*nxw-1,0:2*nyw-1,0:nbuff-1) :: wrin
!    real(kind=DP), intent(out), &
!      dimension(0:2*nxw-1,0:2*nyw-1,0:nbuff-1,0:nprocw-1) :: wrout
    real(kind=DP), intent(in), &
      dimension(np,0:2*nyw-1,0:nbuff-1) :: wrin
    real(kind=DP), intent(out), &
      dimension(np,0:2*nyw-1,0:nbuff-1,0:nprocw-1) :: wrout



!$OMP master

       call MPI_Allgather( wrin,                 &
                           np*2*nyw*nbuff, &
                           MPI_DOUBLE_PRECISION,    &
                           wrout,                &
                           np*2*nyw*nbuff, &
                           MPI_DOUBLE_PRECISION,    &
                           fft_comm_world,        &
                           ierr_mpi )

!$OMP end master



  END SUBROUTINE fout_transpose_back


!--------------------------------------
  SUBROUTINE fout_reduction( iv, fin, fout )
!--------------------------------------

!    real(kind=DP), intent(in), &
!      dimension(0:2*nxw-1,0:2*nyw-1,0:nbuff-1,0:nprocw-1) :: fin
!    real(kind=DP), intent(out), &
!      dimension(np,-nz:nz-1,1:2*nv,0:nm) :: fout
    real(kind=DP), intent(in), &
      dimension(np,0:2*nyw-1,0:nbuff-1,0:nprocw-1) :: fin
    real(kind=DP), intent(out), &
      dimension(np,-nz:nz-1,1:2*nv,0:nm) :: fout
    integer, intent(in) :: iv

    integer :: im, iz, izm, ibuff, iprocw, ip

!$OMP master
                                           call clock_sta(1830)
!$OMP end master

    do im = 0, nm
      do iz=-nz,nz-1
        izm = (2*nz)*im + (iz + nz)
        ibuff = mod(izm, nbuff)
        iprocw = izm / nbuff
        do ip=1, np
           !!fout( ip,iz,iv,im ) = fin(xp(ip),yp(ip),ibuff,iprocw) 
           fout( ip,iz,iv,im ) = fin(ip,yp(ip),ibuff,iprocw) 
        end do
        !write(odbg,'(5i8,999es16.6)') iz, iv, im, ibuff, iprocw, fout(6,iz,iv,im), fin(6,yp(6),ibuff,iprocw) 
      end do
    end do


!$OMP master
                                           call clock_end(1830)
!$OMP end master



  END SUBROUTINE fout_reduction


END MODULE GKV_micouple
