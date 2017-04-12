    SUBROUTINE NWPSAF_Read_ControlData ( &
         RT_Params) ! out 

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
!------------------------------------------------------------------------------
! Description: Reads file containing main control data for NWPSAF 1DVar
!
! Method:
!   Reads control data and makes sure the correct sections of the structures
!   defined in NWPSAFMod_Params and NWPSAFMod_RTModel are assigned.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.0      22/10/99  Original Version 
! 1.1      09/11/00  Added FirstOb and LastOb.  ADC 
! 1.2      10/01/01  Revised for standalone NWPSAF 1DVar.
! 1.3      01/03/01  Added ATOVS Roger Saunders
! 1.4      16/04/02  Fixed small bug where MaxCloudChans was mistyped 
!                    MaxRetChans.   ADC
! 2.3      28/05/02  Changed to Namelist version.   A. Collard.
! 3.0.1    15/07/03  Add cloud retrievals (from M. Szyndel). A. Collard.
! 3.0.4    02/03/04  Changes for more generic instrument treatment. A. Collard.
! 3.0.5    29/03/04  Remove SoundingType_Text.  Add in check for forcing
!                    use of Eqn_101.                         A. Collard.
! 3.1.0    05/04/05  Added support for RTTOV8.               E. Pavelin.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date     Comment
! 13      30/01/09 Added RTTOV9. E. Pavelin.
! 25      20/02/12 Added RTTOV10. P. Weston.
! 28      22/02/12 Added whether Lqtotal is 0 or 1 for cloud liquid water
!                  retrievals.
! 31      18/01/13 Added RTTOV11. P. Weston.
!
! Code Description:
!   Language:           Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY:         &
     StatusOK,                    &
     StatusFatal,                 &
     StatusWarning,               &
     FirstOb, LastOb,             &
     Minimisation_Method,         &
     Newtonian,                   &
     Additional_Cost_Function,    &
     No_Additional_Cost_Function, &
     Perform1DVar,                &
     Force_Eqn_101,               &
     Control,                     &
     Lqtotal,                     &
     Read_CLW_Background,         &
     Output_Dir

USE NWPSAFMod_RTModel, ONLY : &
     RTParams_Type,         &
     RTModelToUse,          &
     RTTOV,                 &
     RTTOV8,                &
     RTTOV9,                &
     RTTOV10,               &
     RTTOV11,               &
     RTIASI,                &
     Gastropod

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
!xxx INCLUDE 'NWPSAF_SetUpInstrumentData.interface'
INCLUDE 'NWPSAF_OpenFile.interface'

! Subroutine arguments:

TYPE(RTParams_type), INTENT(OUT) :: RT_Params        !Info for the RT Model

! Declarations

CHARACTER(LEN=70) :: ErrorMessage(2)  ! Message for Gen_ErrorReport
CHARACTER(LEN=*), PARAMETER :: Filename = 'ControlData.NL'
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NWPSAF_Read_ControlData'
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER (LEN=5) :: Status = "OLD"

INTEGER :: ErrStatRep      = StatusOK
INTEGER :: fileunit              ! I/O unit number
INTEGER :: ReadStatus
!INTEGER :: OpenStatus

!-----------------------------------------------------------------------------

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

ErrorMessage(:) = ''

!--------------------------------------------------------
!1. Open the Control File and Read in Header Data
!--------------------------------------------------------
 
CALL NWPSAF_OpenFile( Trim(FileName),      & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     fileunit )             ! out

! Read Namelist:

IF (FileUnit <= 0) THEN
   Errormessage(1) = 'Error opening Control Namelist'
   ErrStatRep = StatusFatal
   CALL Gen_ErrorReport( RoutineName,  &
        ErrorMessage, &
        ErrStatRep    )
END IF

Read( fileunit, nml=Control, IOSTAT=ReadStatus) 

IF (ReadStatus /= 0) THEN
   Errormessage(1) = ' Error Reading Control Namelist'
   Errormessage(2) = ' You may need to remove comments if compiling with F90'
   ErrStatRep = StatusFatal
   CALL Gen_ErrorReport( RoutineName,  &
        ErrorMessage, &
        ErrStatRep    )
END IF

CLOSE(fileunit)

!------------------------------------------------------
! 2. Check some of the variables from the namelist
!------------------------------------------------------

! First observation to be processed:
!------
IF (FirstOb < 0 ) THEN
  Errormessage(1) = ' Error Reading FirstOb'
  ErrStatRep = StatusFatal
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
END IF

! Last observation to be processed:
!------
IF (LastOb < FirstOb ) THEN
  Errormessage(1) = ' Error Reading LastOb'
  ErrStatRep = StatusFatal
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
END IF

!------------------------------------------------------
! 3. Check some more of the variables from the namelist
!------------------------------------------------------

! Minimisation method:
!------
IF (Minimisation_Method < 0 .OR. Minimisation_Method > 2) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' Minimisation Method',Minimisation_Method,' is not allowed'
  ErrStatRep = StatusFatal
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
ELSE IF (Minimisation_Method == 0) THEN
  Errormessage(1) = &
       ' Minimisation Method 0 chosen - turning off Minimisation'
  Perform1DVar = .FALSE.
  ErrStatRep = StatusWarning
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
END IF

IF ( Force_Eqn_101 .AND. Minimisation_Method /= Newtonian ) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' If Force_Eqn_101 is .TRUE., Minimisation_Method must be 1'
  WRITE(Errormessage(1),FMT=*)  &
       ' Resetting Minimisation_Method'
  ErrStatRep = StatusWarning
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
  Minimisation_Method = Newtonian
END IF

IF ( Force_Eqn_101 .AND. &
     Additional_Cost_Function /= No_Additional_Cost_Function ) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' If Force_Eqn_101 is .TRUE., Additional_Cost_Function must be 1'
  WRITE(Errormessage(1),FMT=*)  &
       ' Resetting Additional_Cost_Function'
  ErrStatRep = StatusWarning
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )
  Additional_Cost_Function = No_Additional_Cost_Function
END IF

! Check the type of Mw cloud liquid water retrieval

IF (Lqtotal < 0 .OR. Lqtotal > 1) THEN
WRITE(Errormessage(1),FMT=*)  &
       ' Lqtotal has to be 1 or 0'
  ErrStatRep = StatusFatal
  CALL Gen_ErrorReport( RoutineName,  &
                        ErrorMessage, &
                        ErrStatRep    )

ENDIF

! RTModel
!---------

SELECT CASE (Trim(RTModelToUse))
  CASE ('RTIASI')
     RT_Params % RTModel = RTIASI
  CASE ('RTTOV')
     RT_Params % RTModel = RTTOV
  CASE ('RTTOV8')
     RT_Params % RTModel = RTTOV8
  CASE ('RTTOV9')
     RT_Params % RTModel = RTTOV9
  CASE ('RTTOV10')
     RT_Params % RTModel = RTTOV10
  CASE ('RTTOV11')
     RT_Params % RTModel = RTTOV11
  CASE ('Gastropod')
     RT_Params % RTModel = Gastropod
  CASE DEFAULT
     WRITE(Errormessage(1),FMT=*)  &
          ' RT Model ',Trim(RTModelToUse),' is not supported'
     ErrorMessage(2) = 'Only RTIASI, RTTOV, RTTOV8, RTTOV9, RTTOV10, RTTOV11 and Gastropod are allowed'
     ErrStatRep = StatusFatal
     CALL Gen_ErrorReport( RoutineName,  &
          ErrorMessage, &
          ErrStatRep    )
END SELECT

END SUBROUTINE NWPSAF_Read_ControlData
