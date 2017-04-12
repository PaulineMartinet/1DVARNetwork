SUBROUTINE RTIASI_rttvi( &
     KNPF,   &
     IERR)

! Description: Initialisation for RTIASI.
!
! Owner:
!     Andrew Collard
!
! Version Date     Comment
! ------- -------- -------
! 1.1     29/07/97 Original code adapted from RTIASI driver (IASRTM)  
!                  (Marco matricardi, ECMWF).  ADC
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code"
!
!     Module used:


IMPLICIT NONE

!  Subroutine arguments:
INTEGER, INTENT(IN):: KNPF    ! Number of processed profiles.
INTEGER, INTENT(OUT):: IERR   ! Error Flag

!       Local parameters:

INTEGER, PARAMETER :: NSATID = 1  ! Number of satellites to be processed.

!       Local scalars:
INTEGER            :: I

!       Local direct arrays:
!     End of program arguments


INTEGER :: KIDSAT(NSATID)        ! Array of satellite  
                                 ! identification indices.

!-----End of header-----------------------------------------------------

DO I=1,NSATID
  KIDSAT(I)=1
END DO

!-----------------------------------------------------------------------
!          1.    INITIALISING ROUTINE FOR RTIASI
!-----------------------------------------------------------------------

CALL RTIASI_RTIAI(IERR,KIDSAT,NSATID,KNPF)


END SUBROUTINE RTIASI_RTTVI
