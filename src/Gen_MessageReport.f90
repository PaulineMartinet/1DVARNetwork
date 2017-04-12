!+ Routine to handle general messages

SUBROUTINE Gen_MessageReport( &
  NameOfRoutine,          & ! in
  Message)                  ! in

! Description:
!   Writes out messages. 
!
! Method:
!
!   The output message is construced from NameOfRoutine and Message(:). 
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0       03/06/98 Stripped down version of OPS subroutine with the
!                    same name.  A. Collard.   Met Office.
!
! Code Description:
!   Language:       Fortran 90.
!   Software Standards: GTDP 8
!

IMPLICIT NONE

!* Subroutine arguments

CHARACTER(LEN=*), INTENT(IN) :: NameOfRoutine ! Calling this one
CHARACTER(LEN=*), INTENT(IN) :: Message(:)    ! Message to output

!* End of Subroutine arguments

! Local scalars:

 Integer                         :: lines         ! number of message lines
 Integer                         :: line          ! current  line of message
!INTEGER :: UnitOut

!- End of header
lines = SIZE(Message, 1)

WRITE(*,'(A)') TRIM(NameOfRoutine)
DO line = 1, lines
  IF (Message(line) /= "") THEN    
    WRITE(*,'(A)') TRIM(Message(line))
  END IF
END DO

END SUBROUTINE Gen_MessageReport
