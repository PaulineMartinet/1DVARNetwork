!+ Opens a file with next free unit number.

SUBROUTINE NWPSAF_OpenFile( &
  file_pathname, & ! in
  WhatAccess,    & ! in
  WhatAction,    & ! in
  WhatStatus,    & ! in
  WhatForm,      & ! inout
  file_unit,     & ! out
  RecordLength   ) ! optional in

!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2004, EUMETSAT, All Rights Reserved.
!
! Description: Opens the file with specified pathname and next free unit number.
!              Various checks are done to check if the file exists, is of the
!              correct type for the inputs etc.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     09/11/99 Original version modified from ATOVS_OpenFile (A.J. Smith)
!                                                         A. Collard
!                                                         Met Office
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,   &
    DebugMode,     &
    StatusWarning, &
    StatusFatal,   &
    StatusOk

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_GetUnit.interface'
INCLUDE 'Gen_MessageReport.interface'

! Subroutine arguments:
CHARACTER (LEN=*), INTENT(IN) :: file_pathname
CHARACTER (LEN=*), INTENT(IN) :: WhatAccess
CHARACTER (LEN=*), INTENT(IN) :: WhatAction
CHARACTER (LEN=*), INTENT(IN) :: WhatForm
CHARACTER (LEN=*), INTENT(IN) :: WhatStatus
INTEGER, INTENT(OUT) :: file_unit
INTEGER, OPTIONAL, INTENT(IN) :: RecordLength

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_OpenFile"

! Local scalars:
INTEGER :: unit_status
INTEGER :: open_status
INTEGER :: pathsize
INTEGER :: ErrStatRep
LOGICAL :: okay_to_open
LOGICAL :: file_exists
CHARACTER(LEN=7)   :: file_format
CHARACTER(LEN=11)  :: ActualForm
CHARACTER(LEN=160) :: ErrorMessage(3)
CHARACTER(LEN=160) :: Message(2)

!-------------------------------------------------------------------------------

!Initialise variables
!--------------------
  unit_status = 0
  open_status = 0
  ErrStatRep  = StatusOk
  ErrorMessage(:) = ' '
  Message(:) = ' '
  okay_to_open = .TRUE.
  ActualForm = Trim(WhatForm)

!Check inputs
!------------
  IF (Trim(WhatAccess)=='DIRECT' .AND. .NOT.(PRESENT(RecordLength))) THEN
    ErrStatRep = StatusFatal
    pathsize = MIN(LEN(file_pathname),120)
    WRITE( UNIT=ErrorMessage(1),FMT='(3A)' ) &
      'Error opening ',file_pathname(1:pathsize)
    ErrorMessage(2) = 'Direct Access requires RECL to be specified'
    CALL Gen_ErrorReport( RoutineName, &
                          ErrorMessage, &
                          ErrStatRep )
  ENDIF

  pathsize = MIN(LEN(file_pathname),160)
  ErrorMessage(1) = file_pathname(1:pathsize)

  IF (TRIM(WhatForm) == 'FORMATTED') THEN
     File_format = 'YES'
  ELSE
     File_format='NO'
  END IF

  INQUIRE( FILE=file_pathname, &
           EXIST=file_exists, &
           FORMATTED=file_format, &
           NUMBER=file_unit )

  IF ( Trim(WhatStatus) == 'OLD' .AND. .NOT. file_exists ) THEN
    File_Unit = 0
    ErrStatRep = StatusWarning
    okay_to_open=.FALSE.
    pathsize = MIN(LEN(file_pathname),120)
    WRITE( UNIT=ErrorMessage(1),FMT='(3A)' ) &
      'Error: required file ',file_pathname(1:pathsize), ' does not exist'
  ELSEIF ( file_unit >=0 ) THEN
    ErrStatRep = StatusWarning
    okay_to_open=.FALSE.
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I4)') &
      'File already open and connected to unit',file_unit
!  ELSEIF ( file_format == 'NO' .AND. Trim(WhatForm) == 'FORMATTED' ) THEN
!    ActualForm = 'UNFORMATTED'
!    ErrStatRep = StatusWarning
!    ErrorMessage(2) = &
!      'Requested formatted input/output which is not permitted for this file'
!    ErrorMessage(3) = &
!      'Ignoring and opening as unformatted...'
!  ELSEIF ( file_format == 'YES' .AND. Trim(WhatForm) == 'UNFORMATTED' ) THEN
!    ActualForm = 'FORMATTED'
!    ErrStatRep = StatusWarning
!    ErrorMessage(2) = &
!      'Requested unformatted input/output which is not permitted for this file'
!    ErrorMessage(3) = &
!      'Ignoring and opening as formatted...'
  ENDIF

  OpenIt: IF ( okay_to_open ) THEN

  !Get unit number and write messages
  !----------------------------------
    CALL Gen_GetUnit(file_unit,unit_status)
    IF ( GeneralMode >= DebugMode ) THEN
      pathsize = MIN(LEN(file_pathname),80)
      WRITE( UNIT=Message(1),FMT='(A)' ) file_pathname(1:pathsize)
        
      IF ( unit_status /= StatusOk ) THEN
        ErrorMessage(1) = Message(1)
        ErrorMessage(2) = '...no unit available'
        ErrStatRep = StatusWarning
      ELSE
        WRITE(Message(2), FMT='(5A,I4)') 'Action=', TRIM(WhatAction), &
             ' Form=', TRIM(ActualForm), ' Assigned unit', file_unit
      ENDIF
    ENDIF

  !Open the file
  !-------------
    IF ( ErrStatRep == StatusOk ) THEN

       IF (PRESENT(RecordLength)) THEN
          OPEN( UNIT   = file_unit, &
               FILE   = file_pathname,&
               IOSTAT = open_status, &
               ACCESS = WhatAccess, &
               ACTION = WhatAction, &
               STATUS = WhatStatus, &
               RECL   = RecordLength, &
               FORM   = ActualForm )
       ELSE
         OPEN( UNIT   = file_unit, &
               FILE   = file_pathname,&
               IOSTAT = open_status, &
               ACCESS = WhatAccess, &
               ACTION = WhatAction, &
               STATUS = WhatStatus, &
               FORM   = ActualForm )
       ENDIF

       IF ( open_status /= 0 ) THEN
          ErrorMessage(1) = Message(1)
          ErrorMessage(2) = '...failed to open'        
          ErrStatRep = StatusWarning
       ENDIF

    ENDIF

  ENDIF OpenIt

!Report errors
!-------------
  IF ( ErrStatRep /= StatusOk ) THEN
    CALL Gen_ErrorReport( RoutineName, &
                          ErrorMessage, &
                          ErrStatRep )
  ELSE
    IF ( GeneralMode >= DebugMode ) THEN
      WRITE(Message(2), FMT=*) TRIM(Message(2)), '...Opened'
      CALL Gen_MessageReport( RoutineName,Message(1:2) )
    ENDIF
  ENDIF

END SUBROUTINE NWPSAF_OpenFile
