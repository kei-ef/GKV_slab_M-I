Module MOD_CLOCK
!
! COPIED FROM GKV
!
   use GKV_header
   use MOD_MPI

   implicit none
 
   private
   public   clock_timer, clock_etime, clock_sta, clock_end, clock_reset

   real(kind=DP), dimension(1:2000), save :: sss, eee, elt
   integer, dimension(1:2000), save :: ccount

CONTAINS
!--------------------------------------
  SUBROUTINE clock_timer( isw )
!--------------------------------------

    integer, intent(in)  :: isw
!!    integer, intent(out) :: iflg

    real(kind=DP), save  :: sss0, eee0
    real(kind=DP)        :: ttotl

      if( isw == 0 ) then

        if( myrank == 0 ) then
          call clock_etime ( sss0 )
        end if

        call MPI_Bcast( sss0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr )

        !!iflg = 0

        sss(1:2000) = 0._DP
        eee(1:2000) = 0._DP
        elt(1:2000) = 0._DP
        ccount(1:2000) = 0

      else if( isw == 2 ) then

        if( myrank == 0 ) then
          call clock_etime ( eee0 )
        end if

        call MPI_Bcast( eee0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr )

        !!iflg = 0
        
        if( myrank==0) then

           ttotl   = eee0 - sss0
           write( olog, * ) ""
           write( olog, * ) " ######### elapsed time (sec) and call count #########"
           write( olog, '(a22,f15.5,i15)' ) " total               = ", ttotl
           write( olog, '(a22,f15.5,i15)' ) " fxv                 = ", elt(10), ccount(10)
           write( olog, '(a22,f15.5,i15)' ) "  - read             = ", elt(12), ccount(12)
           write( olog, '(a22,f15.5,i15)' ) "  - write            = ", elt(15), ccount(15)
           write( olog, '(a22,f15.5,i15)' ) " pot (magnetosphere) = ", elt(20), ccount(20)
           write( olog, '(a22,f15.5,i15)' ) "  - write            = ", elt(25), ccount(25)
           write( olog, '(a22,f15.5,i15)' ) " iono                = ", elt(30), ccount(30)
           write( olog, '(a22,f15.5,i15)' ) "  - write            = ", elt(35), ccount(35)
           write( olog, '(a22,f15.5,i15)' ) " fft                 = ", elt(100), ccount(100)
           write( olog, * ) ""
        
         end if

      end if


  END SUBROUTINE clock_timer

!--------------------------------------
  SUBROUTINE clock_etime( ttt )
!--------------------------------------

    real(kind=DP) :: ttt


      ttt = MPI_Wtime()


  END SUBROUTINE clock_etime


!--------------------------------------
  SUBROUTINE clock_sta( id )
!--------------------------------------

    integer, intent(in) :: id


      call clock_etime( sss(id) )


  END SUBROUTINE clock_sta


!--------------------------------------
  SUBROUTINE clock_end( id )
!--------------------------------------

    integer, intent(in) :: id


      call clock_etime( eee(id) )

      elt(id) = elt(id) + eee(id) - sss(id)

      ccount(id)=ccount(id)+1


  END SUBROUTINE clock_end

!--------------------------------------
  SUBROUTINE clock_reset
!--------------------------------------


      elt(10:2000) = 0._DP
      ccount(10:2000) = 0


  END SUBROUTINE clock_reset

END Module MOD_CLOCK


