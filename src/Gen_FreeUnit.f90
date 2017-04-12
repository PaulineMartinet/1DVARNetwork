!+ Frees a fortran I/O unit number assigned by Gen_GetUnit.

SUBROUTINE Gen_FreeUnit( &
  unit,               & ! InOut
  ErrorStatus )         ! Out  

! Description:
!   Frees a fortran I/O unit number previously assigned by Gen_GetUnit.
!
! Method:
!   The module GenMod_FortranIO keeps track of which unit numbers have been 
!   assigned.
!
! Owner: John Bray
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       21/07/97 Creation from VarMod_FortranIO 1.8 JB
! 1.2       04/09/97 ErrorStatus INTENT(OUT). JB
! 1.4       12/01/98 ErrorStatus OK if no problem. JB
!
! Code Description:
!   Language:       Fortran 90.
!   Software Standards: GTDP 8
!
! Parent module: GenMod_FortranIO
!
! Declarations:

! Modules used:

USE NWPSAFMod_Params, ONLY :  &
    LB, &
    UB, &
    UnitsInUse, &
    StatusOK,                           & ! Good ErrorStatus value.
    StatusWarning                         ! Warning error ErrorStatus value.

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'

! Subroutine arguments

INTEGER, INTENT(INOUT) :: unit        ! Unit number supplied
INTEGER, INTENT(OUT)   :: ErrorStatus ! Error status flag

! Local constants

CHARACTER (LEN=*), PARAMETER :: RoutineName = "Gen_FreeUnit"

!- End of header

IF (unit >= LB .AND. unit < UB) THEN
  UnitsInUse(unit) = .FALSE.
  unit = 0
  ErrorStatus = StatusOK
ELSE
  ErrorStatus = StatusWarning
  CALL Gen_ErrorReport(                                        &
     RoutineName,                                          &
     (/ 'Attempt to free an illegal unit number' /),          &
     ErrorStatus )
END IF

END SUBROUTINE Gen_FreeUnit
