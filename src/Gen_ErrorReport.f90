!+ Routine to handle error reports

SUBROUTINE Gen_ErrorReport( &
  NameOfRoutine,          & ! in
  Message,                & ! in
  ErrorStatus)              ! in

! Description:
!   Writes out fatal and warning error messages. Execution is terminated
!   in the event of a fatal error.
!
! Method:
!   If ErrorStatus = 0 Var_MessageReport called. 
!   If ErrorStatus < 0 a Warning message is output.
!   If ErrorStatus > 0 a Fatal Error message is output and program is 
!                                                       aborted.
!
!   The output message is constructed from NameOfRoutine and Message(:). 
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0       03/06/98 Stripped down version of OPS routine with the same
!                    name for NWPSAF 1DVar.  A. Collard  Met Office.
!
! Code Description:
!   Language:       Fortran 90.
!   Software Standards: GTDP 8
!
! Declarations:
!
! Modules used:


IMPLICIT NONE

INCLUDE 'Gen_MessageReport.interface'

!* Subroutine arguments

CHARACTER(LEN=*), INTENT(IN) :: NameOfRoutine ! Calling this one
CHARACTER(LEN=*), INTENT(IN) :: Message(:)    ! Message to output
INTEGER, INTENT(IN)         :: ErrorStatus   ! Input error code

!* End of Subroutine arguments

! Local variables

INTEGER                         :: lines         ! number of message lines
CHARACTER (LEN=120), ALLOCATABLE :: LocalMessage(:)

!- End of header

lines = SIZE(Message, 1)

ALLOCATE (LocalMessage(lines+1))

LocalMessage(2:)=Message(1:)

IF (ErrorStatus < 0) THEN
      LocalMessage(1)="WARNING"   

  CALL Gen_MessageReport(       &
    NameOfRoutine,              &
    LocalMessage)

ELSE IF (ErrorStatus > 0) THEN
      LocalMessage(1)="FATAL ERROR"

  CALL Gen_MessageReport(       &
    NameOfRoutine,              &
    LocalMessage)

  STOP

ELSE ! 0 is OK, just print message
   CALL Gen_MessageReport ( &
      NameOfRoutine,        & ! in
      Message)          ! in, optional
END IF

DEALLOCATE (LocalMessage)

END SUBROUTINE Gen_ErrorReport
