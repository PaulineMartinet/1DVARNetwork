Subroutine NWPSAF_RTTOV7_Interface ( &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! inout
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrorCode)                 ! out
       
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
! Description: Interface between NWPSAF 1DVar code and RTTOV7
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     12/11/01 Original version - split out from Fastmodel_Interface
! 2.2     22/01/02 Set CLW in RTTOV to -1 (disables all CLW calculations)
!                  Add option to initialise a subset of channels.
! 2.3     22/05/02 Change absolute channel number from rank 2 to rank 1 
!                  array.                                     A. Collard.
! 3.0.1   15/07/03 Add cloudy radiances (M. Szyndel)          A. Collard.
! 3.0.4   03/03/04 Use unmodified RTTOV subroutines.          A. Collard.
! 3.0.5   29/03/04 Change of warning message.                 A. Collard.
! 3.0.6   18/06/04 Include tolerances in soft limits.         A. Collard.
! 3.0.7   08/01/07 Changed minimum soft limits for u & v to 100m/s. E. Pavelin.
!
! Code Description:
!   Language:		Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     StatusWarning,      & 
     GeneralMode,        &
     VerboseMode,        &
     DebugMode,          &
     Retrieved_Elements, &
     Ret_FirstQ,         &
     Ret_LastQ,          &
     Ret_q2,             &
     CloudyRetrieval,    &
     Cloud_Min_Pressure

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type,            &
     FastmodelMode_CleanUp,    &
     FastmodelMode_Initialise, &
     FastmodelMode_Forward,    &
     FastmodelMode_Gradient,   &
     GuessProf, &
     BackGrProf, &
     ProfSize,  &
     Prof_FirstT, &
     Prof_FirstQ, &
     Prof_FirstO3, &
     Prof_LastT, &
     Prof_LastQ, &
     Prof_LastO3, &
     Prof_T2, &
     Prof_q2, &
     Prof_Tstar, &
     Prof_pstar, &
     Prof_uwind, &
     Prof_vwind, &
     Prof_CTP, &
     Prof_CloudCover, &
     Num_RTlevels, &
     Num_Profs, &
     Soft_Limits, &
     ! The following are RTTOV specific
     Num_Press_Levels,   &
     Num_Atm_Prof_Vars,  &
     Num_Surf_Vars,      &
     Num_Skin_Vars,      &
     Num_Cloud_Vars,     &
     Max_Num_Sats

USE MOD_CPARAM, ONLY : &
     jpch, &
     jpnsat

USE MOD_TOVCHN, ONLY: &
     ChanConstants        ! Calibration Constants

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_MessageReport.interface'

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrorCode


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV7_Interface"

REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv
REAL, PARAMETER ::  O3_MassMixToPPMV  = 6.035045e5  ! kg/kg -> ppmv

! Local variables:

INTEGER :: I, IChan1, IChan2, Variable_Number
INTEGER :: ErrStatRep                 ! error status for Gen_ErrorReport
INTEGER :: First_Instr, Last_Instr
INTEGER :: IPrint = 0
INTEGER :: Instrument
INTEGER :: Radiance_Or_BT = 2
INTEGER :: NumInstChans
INTEGER :: SatIndex

CHARACTER(LEN=80) :: ErrorMessage(4)  ! Message for Gen_ErrorReport
CHARACTER(LEN=80) :: Message(2)       ! Message for Gen_MessageReport

! Arrays used in interface with both RT Models
!-----------

INTEGER :: First_Chan_Pos !) Position of first and last channels for  
INTEGER :: Last_Chan_Pos  !) current instrument in RT_Params % TotalBTs
INTEGER :: KPROF(UsedChans % NumChans)                        
INTEGER :: Inst_Chans (UsedChans % NumChans)
INTEGER :: RTModel_RTSurfaceType(Num_Profs)
REAL :: RTModel_Jacobian (ProfSize, UsedChans % NumChans)
REAL :: RTModel_SatZenithAngle(Num_Profs)          
REAL :: RTModel_SolarZenithAngle(Num_Profs)        
! REAL :: Surf_Emiss(UsedChans % NumChans * Num_Profs)           
REAL :: Surf_Emiss_IO(UsedChans % NumChans * Num_Profs)           
REAL :: Surf_Emiss_K(UsedChans % NumChans * Num_Profs)           

REAL :: RTM_TotalRads(UsedChans % NumChans * Num_Profs)          
REAL :: RTM_TotalBTs(UsedChans % NumChans * Num_Profs)                     
 
REAL, POINTER :: RTProf(:)

! Arrays used in RTTOV interface (should be set to zero size if RTTOV is
! not to be used) 
!--------

REAL :: Atmos_Prof(Num_Press_Levels, Num_Atm_Prof_Vars, Num_Profs)
REAL :: Surf_Vars(Num_Surf_Vars, Num_Profs)            
REAL :: Skin_Vars(Num_Skin_Vars, Num_Profs)         
REAL :: Cloud_Vars(Num_Cloud_Vars, Num_Profs)                
INTEGER :: IFail(Max_Num_Sats)                        
REAL :: OCast_Rads(UsedChans % NumChans, 2*Num_Press_Levels+2) 
REAL :: OCast_Rads_at_CldTop(UsedChans % NumChans ) 
REAL :: Trans_to_Space(UsedChans%NumChans, Num_Press_Levels)
REAL :: Surf_to_Space_Trans(UsedChans % NumChans)
REAL :: Atmos_Prof_K(Num_Press_Levels, Num_Atm_Prof_Vars, &
     UsedChans % NumChans * Num_Profs )
REAL :: Surf_Vars_K(Num_Surf_Vars, UsedChans % NumChans * Num_Profs)   
REAL :: Skin_Vars_K(Num_Skin_Vars, UsedChans % NumChans * Num_Profs)
REAL :: Cloud_Vars_K(Num_Cloud_Vars, UsedChans % NumChans * Num_Profs)
INTEGER :: Valid_Channels(jpch,jpnsat)

! Variables used in initialising RTTOV interface 
!--------

INTEGER :: Actual_Num_Instr
INTEGER :: Max_Parallel_Profiles
INTEGER :: Max_Channels
INTEGER :: Max_Channels_Used
INTEGER, POINTER :: First_Instr_Chan
INTEGER, POINTER :: Num_Instr_Chan 
INTEGER :: Last_Instr_Chan 
REAL, ALLOCATABLE :: RTTOV_Pressures(:)
REAL, ALLOCATABLE :: Minimum_Temperature(:)
REAL, ALLOCATABLE :: Maximum_Temperature(:)
REAL, ALLOCATABLE :: Minimum_Humidity(:)
REAL, ALLOCATABLE :: Maximum_Humidity(:)
REAL, ALLOCATABLE :: Minimum_Ozone(:)
REAL, ALLOCATABLE :: Maximum_Ozone(:)
LOGICAL :: Cloud_Switch

!---------------------------------------------------------------

IF ( GeneralMode >= VerboseMode ) THEN
   WRITE(UNIT=Message(1),FMT='("RTTOV7 Called in MODE = ",i3)') &
        Fastmodel_Mode
   CALL Gen_MessageReport( RoutineName,Message(1:1) )
END IF

! Error code set to zero here

ErrorCode = 0
ErrorMessage(:)=' '

! Check that Num_Profs is set correctly (it's just too complicated to 
! have it set to anything but one).

IF (Num_Profs /= 1) THEN
   ErrorMessage(1)='Num_Profs should be set to 1'
   WRITE( UNIT=ErrorMessage(2),FMT='(A,I5)' ) &
        'Currently it is ',Num_Profs
   ErrStatRep = StatusFatal
   CALL Gen_ErrorReport( RoutineName,  &
        ErrorMessage, &
        ErrStatRep    )
END IF

! Start calls to RTTOV7

RTTOV_FastmodelMode : IF &
     ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
   IF ( GeneralMode >= VerboseMode ) THEN
      Message(1) = 'Initialising RTTOV'
      CALL Gen_MessageReport( RoutineName,Message(1:1) )
   ENDIF
   ALLOCATE(RTTOV_Pressures(Num_RTLevels))
   ALLOCATE(Minimum_Temperature(Num_RTLevels))
   ALLOCATE(Maximum_Temperature(Num_RTLevels))
   ALLOCATE(Minimum_Humidity(Num_RTLevels))
   ALLOCATE(Maximum_Humidity(Num_RTLevels))
   ALLOCATE(Minimum_Ozone(Num_RTLevels))
   ALLOCATE(Maximum_Ozone(Num_RTLevels))
   Valid_Channels(:,:)=0

   IF (RT_Params % Num_Instruments > jpnsat) THEN
      ErrorMessage(1)='RT_Params % Num_Instruments > jpnsat'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I5,A,I5)' ) &
           'Currently they are ',RT_Params % Num_Instruments,' and ',jpnsat
      ErrorMessage(3)='Either reduce the number of instruments required in'// &
           ' the observations file ' 
      ErrorMessage(4)='or increase jpnsat in MOD_CPARAMS.f90'
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
    

   DO I = 1, RT_Params % Num_Instruments
      First_Instr_Chan => RT_Params % First_Channel_for_Instrument(I)
      Num_Instr_Chan   => RT_Params % NumChans(I)
      IF (Num_Instr_Chan > 0) Actual_Num_Instr = I
      Last_Instr_Chan  =  First_Instr_Chan + Num_Instr_Chan - 1
      Valid_Channels(1:Num_Instr_Chan,I) = RT_Params % &
           Absolute_Channel_Number(First_Instr_Chan:Last_Instr_Chan)
      IF (jpch > Num_Instr_Chan) &
           Valid_Channels(Num_Instr_Chan+1:jpch,I) = 0
   END DO

   IF (RT_Params % Num_Instruments < jpnsat) &
        Valid_Channels(:,RT_Params % Num_Instruments+1:jpnsat) = 0

   ! Reset number of instruments to prevent some unnecessary initialisations
   ! (this isn't the perfect solution yet)
   RT_Params % Num_Instruments = Actual_Num_Instr

   CALL RTTVI ( & 
        ErrorCode,                      & !out
        Max_Parallel_Profiles,          & !out
        Max_Num_Sats,                   & !out
        Num_Press_Levels,               & !out
        Max_Channels,                   & !out
        Max_Channels_Used,              & !out
        Num_Atm_Prof_Vars,              & !out
        Num_Surf_Vars,                  & !out
        Num_Skin_Vars,                  & !out
        Num_Cloud_Vars,                 & !out
        RT_Params % Num_Instruments,    & !in
        RT_Params % SeriesChoice(1:Actual_Num_Instr),    & !in
        RT_Params % PlatformChoice(1:Actual_Num_Instr),  & !in
        RT_Params % SubTypeChoice(1:Actual_Num_Instr),   & !in
        RT_Params % NumChans(1:Actual_Num_Instr),        & !inout
        RTTOV_Pressures(:),             & !out
        Minimum_Temperature(:),         & !out
        Maximum_Temperature(:),         & !out
        Minimum_Humidity(:),            & !out
        Maximum_Humidity(:),            & !out
        Minimum_Ozone(:),               & !out
        Maximum_Ozone(:),               & !out
        Valid_Channels(:,:),            & !inout
        RT_Params % RTCoeffs(:))          !inout

   IF (ErrorCode /= 0) THEN 
      ErrorMessage(1)='Error in RTTOV_RTTVI'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
           'Error is Code ',ErrorCode
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
   
   IF (Max_Parallel_Profiles /= Num_Profs) THEN
      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
           'Mismatch in number of parallel profiles on '// &
           'initialising RTTOV'
      WRITE( UNIT=ErrorMessage(2),FMT=* ) &
           'Required = ',Num_Profs,', from RTTOV=',Max_Parallel_Profiles
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
   
   IF (Num_Press_Levels /= Num_RTLevels) THEN
      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
           'Mismatch in number of pressure levels on '// &
           'initialising RTTOV'
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
    
   ! Set up soft limits

   ALLOCATE( Soft_Limits % Minimum(ProfSize))
   ALLOCATE( Soft_Limits % Maximum(ProfSize))

   ! Default values:
   Soft_Limits % Minimum(:) = 0.
   Soft_Limits % Maximum(:) = 1.e10

   ! Profile minima (include allowable tolerances)
   Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = Minimum_Temperature(:)-0.5
   Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = LOG(Minimum_Humidity(:))-0.5
   Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = &
        Minimum_Ozone(:)*O3_MassMixToPPMV*0.8
   Soft_Limits % Minimum(Prof_T2) = Minimum_Temperature(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_Q2) = LOG(Minimum_Humidity(Num_RTLevels))-0.5
   Soft_Limits % Minimum(Prof_PStar) = 300.
   Soft_Limits % Minimum(Prof_UWind) = -100.
   Soft_Limits % Minimum(Prof_VWind) = -100.
   Soft_Limits % Minimum(Prof_TStar) = Minimum_Temperature(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_CTP) = Cloud_Min_Pressure
   Soft_Limits % Minimum(Prof_CloudCover) = 0.

   ! Profile maxima
   Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = Maximum_Temperature(:)+0.5
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = LOG(Maximum_Humidity(:))+0.5
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
        Maximum_Ozone(:)*O3_MassMixToPPMV*1.2
   Soft_Limits % Maximum(Prof_T2) = Maximum_Temperature(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_Q2) = LOG(Maximum_Humidity(Num_RTLevels))+0.5
   ! (This is the RTTOV hard limit:) 
   Soft_Limits % Maximum(Prof_PStar) = 1200.
   Soft_Limits % Maximum(Prof_UWind) = 100.
   Soft_Limits % Maximum(Prof_VWind) = 100.
   Soft_Limits % Maximum(Prof_TStar) = Maximum_Temperature(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_CTP) = 1200.
   Soft_Limits % Maximum(Prof_CloudCover) = 1.0

   ! These are the RTTOV hard limits:
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
        MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),-2.99578)
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
        MIN(Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3),20.)
   Soft_Limits % Maximum(Prof_Q2) = MIN(Soft_Limits % Maximum(Prof_Q2),-2.99578)

   DO I=1,ProfSize
      IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Mismatch in RTTOV soft limits'
         WRITE( UNIT=ErrorMessage(2),FMT=* )  I, &
              Soft_Limits % Minimum(I), &
              Soft_Limits % Maximum(I)
         ErrStatRep = StatusWarning
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END IF
   END DO

   ! Deallocate some arrays

   DEALLOCATE(RTTOV_Pressures)
   DEALLOCATE(Minimum_Temperature)
   DEALLOCATE(Maximum_Temperature)
   DEALLOCATE(Minimum_Humidity)
   DEALLOCATE(Maximum_Humidity)
   DEALLOCATE(Minimum_Ozone)
   DEALLOCATE(Maximum_Ozone)

ELSE IF (FastModel_Mode == FastModelMode_CleanUp) THEN

   CALL CLEANUP

ELSE

   SatIndex = RT_Params % SatIndex
   First_Instr = RT_Params % SatID(SatIndex) % First_Instr
   Last_Instr = RT_Params % SatID(SatIndex) % Last_Instr

   !---------------------------------------------------------------------
   !  Put profile into the RTModel vector remembering to change the 
   !  humidity variable accordingly.  (This is currently only used
   !  when the Jacobian units are being converted, but has more utility
   !  if the entire profile vector is passed as one.)
   !---------------------------------------------------------------------

   IF (WhichProf == BackGrProf) THEN
      RTProf => RT_Params % RTBack
   ELSE IF (WhichProf == GuessProf) THEN
      RTProf => RT_Params % RTGuess
   ELSE
      WRITE( UNIT=ErrorMessage(1),FMT='(A,I2,A)' )  &
           'Incorrect profile (',WhichProf,') specified'
      ErrStatRep = StatusWarning
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF

   !  Unpack the RTProf Profile Vector into the input arrays used for 
   ! RTTOV:
   !------
   
   ! First do the atmospheric profile variables:
   DO Variable_Number = 1, Num_Atm_Prof_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is temperature
         Atmos_Prof(:,Variable_Number,1) = RTProf(Prof_FirstT : Prof_LastT)
      CASE (2)  ! This is humidity (in kg/kg)(convert from LOG)
         Atmos_Prof(:,Variable_Number,1) = EXP(RTProf(Prof_FirstQ:Prof_LastQ)) 
      CASE (3)  ! This is Ozone
         Atmos_Prof(:,Variable_Number,1) = &
              RTProf(Prof_FirstO3 : Prof_LastO3) / O3_MassMixToPPMV
      CASE (4)  ! This is Cloud Liquid Water (-1.0 == disabled for now).
         Atmos_Prof(:, Variable_Number, 1) = -1.0
      CASE DEFAULT
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Not expecting more than four atmospheric'// &
              'profile elements'
         ErrStatRep = StatusFatal
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END SELECT
   END DO
   ! Now do the Surface Variables
   DO Variable_Number = 1, Num_Surf_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is surface temperature
         Surf_Vars(Variable_Number, 1) = RTProf(Prof_T2)
      CASE (2)  ! This is surface humidity (kg/kg) (convert from LOG)
         Surf_Vars(Variable_Number, 1) = EXP(RTProf(Prof_Q2))
      CASE (3)  ! This is surface pressure
         Surf_Vars(Variable_Number, 1) = RTProf(Prof_PStar)
      CASE (4)  ! This is U wind
         Surf_Vars(Variable_Number, 1) = RTProf(Prof_UWind)
      CASE (5)  ! This is V wind
         Surf_Vars(Variable_Number, 1) = RTProf(Prof_VWind)
      CASE DEFAULT
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Not expecting more than five surface variables'
         ErrStatRep = StatusFatal
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END SELECT
   END DO
   ! Now do the Skin Variables
   DO Variable_Number = 1, Num_Skin_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is skin temperature
         Skin_Vars(Variable_Number, 1) = RTProf(Prof_TStar)
      CASE (2,3,4,5,6)  ! Fastem parameters.  Set to zero for now.
         Skin_Vars(Variable_Number, 1) = 0.0
      CASE DEFAULT 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Not expecting more than six skin variables'
         ErrStatRep = StatusFatal
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END SELECT
   END DO
   ! Now do the Cloud Variables
   DO Variable_Number = 1, Num_Cloud_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is cloud top pressure (make sure is strictly within the 
                ! atmosphere as otherwise RTTOV may crash!)
         Cloud_Vars(Variable_Number, 1) = RTProf(Prof_CTP)
         Cloud_Vars(Variable_Number, 1) = &
              MIN(Cloud_Vars(Variable_Number, 1),RTProf(Prof_PStar))
         Cloud_Vars(Variable_Number, 1) = &
              MAX(Cloud_Vars(Variable_Number, 1), &
              MINVAL(RT_Params % Pressure_Pa(:) / 100.))
      CASE (2)  ! This is cloud fraction (make sure is strictly within the 
                ! range 0-1 as otherwise RTTOV may crash!)
         Cloud_Vars(Variable_Number, 1) = RTProf(Prof_CloudCover)
         Cloud_Vars(Variable_Number, 1) = MIN(Cloud_Vars(Variable_Number, 1),1.0)
         Cloud_Vars(Variable_Number, 1) = MAX(Cloud_Vars(Variable_Number, 1),0.0)
       CASE DEFAULT 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Not expecting more than two cloud variables'
         ErrStatRep = StatusFatal
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END SELECT
   END DO
   ! End of profile assignments
   
   ! Use default surface emissivities for now
   ! Surf_Emiss(:)=0.

   ! No cloud processing at present
   Cloud_Switch = .FALSE.
   
   RTModel_SatZenithAngle(1) = RT_Params % SatZenithAngle
   RTModel_SolarZenithAngle(1) = RT_Params % SolarZenAngle
   RTModel_RTSurfaceType(1) = RT_Params % RTSurfaceType
   KPROF(1:UsedChans % NumChans) = 1
   
   IF (GeneralMode >= DebugMode) IPrint = 1
   SELECT CASE ( Fastmodel_Mode )
     CASE( FastmodelMode_Forward )
        Instrument_Loop : DO Instrument = First_Instr, Last_Instr
           NumInstChans = &
                COUNT(RT_Params % Instrument_Number( &
                SatIndex,UsedChans % Channels) == Instrument)
           IChan2 = 1
           First_Chan_Pos = 0
           Last_Chan_Pos = 0
           DO IChan1 = 1, UsedChans % NumChans
              IF (RT_Params % Instrument_Number(SatIndex, &
                   UsedChans % Channels(IChan1)) == Instrument) THEN
                 Inst_Chans(IChan2) = UsedChans % Channels(IChan1) - &
                      RT_Params % First_Channel_for_Instrument(Instrument) + 1
                 IF (First_Chan_Pos == 0) First_Chan_Pos = IChan1
                 Last_Chan_Pos = IChan1
                 IChan2 = IChan2 + 1
              END IF
           END DO

           IF (First_Chan_Pos > 0) THEN

              Surf_Emiss_IO(1:NumInstChans) = 0.

              CALL RTTOV( &
                   Num_Profs,                            & ! in (Should be 1!)
                   Num_Press_Levels,                     & ! in
                   RT_Params % Pressure_Pa(:) / 100.,    & ! in
                   RTModel_SatZenithAngle(:),            & ! in
                   RTModel_SolarZenithAngle(:),          & ! in
                   RTModel_RTSurfaceType(:),             & ! in
                   Instrument,                           & ! in
                   NumInstChans,                         & ! in
                   Inst_Chans(1:NumInstChans),           & ! in
                   KPROF(:),                             & ! in
                   Atmos_Prof(:,:,:),                    & ! in
                   Surf_Vars(:,:),                       & ! in
                   Skin_Vars(:,:),                       & ! in
                   Cloud_Vars(:,:),                      & ! in
                   Surf_Emiss_IO(1:NumInstChans),        & ! inout
                   IFail(:),                             & ! out
                   RTM_TotalRads(1:NumInstChans),        & ! out
                   RTM_TotalBTs(1:NumInstChans),         & ! out
                   OCast_Rads(1:NumInstChans,:),         & ! out
                   OCast_Rads_at_CldTop(1:NumInstChans), & ! out
                   Trans_to_Space(1:NumInstChans,:),     & ! out
                   Surf_to_Space_Trans(1:NumInstChans),  & ! out
                   Cloud_Switch)                           ! in
              
              ErrorCode=MAXVAL(IFail)
              ! Assume soft limts have been set in this routine, so only
              ! report a catestrophic error!
              IF (ErrorCode < 20) ErrorCode=0

              !  Unpack the required output arrays from RTTOV into the 
              !  RT_Param structure
              !------

              RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
                   RTM_TotalRads(1:NumInstChans) 
              RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                   RTM_TotalBTs(1:NumInstChans)
              IF (CloudyRetrieval) THEN
                 RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      1:2*Num_Press_Levels+2) = &
                      OCast_Rads(1:UsedChans % NumChans,1:2*Num_Press_Levels+2)
                 RT_Params % tc1(First_Chan_Pos:Last_Chan_Pos) = &
                      ChanConstants(Instrument) % &
                          tc1(Inst_Chans(1:NumInstChans))
                 RT_Params % tc2(First_Chan_Pos:Last_Chan_Pos) = &
                      ChanConstants(Instrument) % &
                          tc2(Inst_Chans(1:NumInstChans))
                 RT_Params % bcon1(First_Chan_Pos:Last_Chan_Pos) = &
                      ChanConstants(Instrument) % &
                          bcon1(Inst_Chans(1:NumInstChans))
                 RT_Params % bcon2(First_Chan_Pos:Last_Chan_Pos) = &
                      ChanConstants(Instrument) % &
                           bcon2(Inst_Chans(1:NumInstChans))
              END IF

           END IF

        END DO Instrument_Loop  

     CASE ( FastmodelMode_Gradient )
        Instrument_LoopK : DO Instrument = First_Instr, Last_Instr
           NumInstChans = &
                COUNT(RT_Params % Instrument_Number( &
                SatIndex,UsedChans % Channels) == Instrument)
           IChan2 = 1
           First_Chan_Pos = 0
           Last_Chan_Pos = 0
           DO IChan1 = 1, UsedChans % NumChans
              IF (RT_Params % Instrument_Number(SatIndex, &
                   UsedChans % Channels(IChan1)) == Instrument) THEN
                 Inst_Chans(IChan2) = UsedChans % Channels(IChan1) - &
                      RT_Params % First_Channel_for_Instrument(Instrument) + 1
                 IF (First_Chan_Pos == 0) First_Chan_Pos = IChan1
                 Last_Chan_Pos = IChan1
                 IChan2 = IChan2 + 1
              END IF
           END DO

           IF (First_Chan_Pos > 0) THEN

              Surf_Emiss_IO(1:NumInstChans) = 0.
              OCast_Rads(:,:) = 0.

              CALL RTTOVK( &
                   Num_Profs,                & ! in (This should be 1!)
                   Num_Press_Levels,                             & ! in
                   RT_Params % Pressure_Pa(:) / 100.,            & ! in
                   RTModel_SatZenithAngle(:),                    & ! in
                   RTModel_SolarZenithAngle(:),                  & ! in
                   RTModel_RTSurfaceType(:),                     & ! in
                   Instrument,                                   & ! in
                   NumInstChans,                                 & ! in
                   Inst_Chans(1:NumInstChans),                   & ! in
                   KPROF(:),                                     & ! in
                   Atmos_Prof_K(:,:,1:NumInstChans),             & ! out
                   Surf_Vars_K(:,1:NumInstChans),                & ! out
                   Skin_Vars_K(:,1:NumInstChans),                & ! out
                   Cloud_Vars_K(:,1:NumInstChans),               & ! out
                   Surf_Emiss_K(1:NumInstChans),                 & ! out
                   Atmos_Prof(:,:,:),                            & ! in
                   Surf_Vars(:,:),                               & ! in
                   Skin_Vars(:,:),                               & ! in 
                   Cloud_Vars(:,:),                              & ! in
                   Surf_Emiss_IO(1:NumInstChans),                & ! inout
                   RTM_TotalRads(1:NumInstChans),                & ! out
                   RTM_TotalBTs(1:NumInstChans),                 & ! out
                   Radiance_Or_BT,                               & ! in
                   Cloud_Switch,                                 & ! in
                   IFail(:),                                     & ! out
                   OCast_Rads(1:NumInstChans,:))                   ! inout

              ErrorCode=MAXVAL(IFail)
              ! Assume soft limts have been set in this routine, so only
              ! report a catestrophic error!
              IF (ErrorCode < 20) ErrorCode=0

              !  Unpack the required output arrays from RTTOV into the 
              !  RT_Param structure
              !------

              RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
                   RTM_TotalRads(1:NumInstChans) 
              RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                   RTM_TotalBTs(1:NumInstChans)
              IF (CloudyRetrieval) &
                   RT_Params % CloudyRadiances(1:UsedChans % NumChans, &
                   1:2*Num_Press_Levels+2) = OCast_Rads &
                   (1:UsedChans % NumChans,1:2*Num_Press_Levels+2)

 
              !  Unpack the output K arrays from RTTOV into RTModel_Jacobian:
              !------
              
              ! First do the atmospheric profile variables:
              DO Variable_Number = 1, Num_Atm_Prof_Vars
                 SELECT CASE (Variable_Number)
                   CASE (1)  ! This is temperature
                      RTModel_Jacobian(Prof_FirstT:Prof_LastT, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Atmos_Prof_K(:,Variable_Number, 1:NumInstChans)
                   CASE (2)  ! This is humidity
                      RTModel_Jacobian(Prof_FirstQ:Prof_LastQ, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Atmos_Prof_K(:,Variable_Number, 1:NumInstChans)
                   CASE (3)  ! This is Ozone
                      RTModel_Jacobian(Prof_FirstO3:Prof_LastO3, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Atmos_Prof_K(:,Variable_Number, 1:NumInstChans) &
                           * O3_MassMixToPPMV
                   CASE (4)  ! (Don't worry about cloud liquid water for now)
                   CASE DEFAULT
                      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                           'Not expecting more than four atmospheric'// &
                           'profile elements'
                      ErrStatRep = StatusFatal
                      CALL Gen_ErrorReport( RoutineName,  &
                           ErrorMessage, &
                           ErrStatRep    )
                 END SELECT
              END DO
              ! Now do the Surface Variables
              DO Variable_Number = 1, Num_Surf_Vars
                 SELECT CASE (Variable_Number)
                   CASE (1)  ! This is surface temperature
                      RTModel_Jacobian(Prof_T2, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Surf_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (2)  ! This is surface humidity (kg/kg)
                      RTModel_Jacobian(Prof_Q2, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Surf_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (3)  ! This is surface pressure
                      RTModel_Jacobian(Prof_PStar, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Surf_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (4)  ! This is U wind
                      RTModel_Jacobian(Prof_UWind, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Surf_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (5)  ! This is V wind
                      RTModel_Jacobian(Prof_VWind, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Surf_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE DEFAULT
                      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                           'Not expecting more than five surface variables'
                      ErrStatRep = StatusFatal
                      CALL Gen_ErrorReport( RoutineName,  &
                           ErrorMessage, &
                           ErrStatRep    )
                 END SELECT
              END DO
              ! Now do the Skin Variables
              DO Variable_Number = 1, Num_Skin_Vars
                 SELECT CASE (Variable_Number)
                   CASE (1)  ! This is skin temperature
                      RTModel_Jacobian(Prof_TStar, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Skin_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (2,3,4,5,6)  ! Fastem parameters.  (ignore)
                   CASE DEFAULT 
                      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                           'Not expecting more than six skin variables'
                      ErrStatRep = StatusFatal
                      CALL Gen_ErrorReport( RoutineName,  &
                           ErrorMessage, &
                           ErrStatRep    )
                 END SELECT
              END DO
              ! Now do the Cloud Variables
              DO Variable_Number = 1, Num_Cloud_Vars
                 SELECT CASE (Variable_Number)
                   CASE (1)  ! This is cloud top pressure
                      RTModel_Jacobian(Prof_CTP,First_Chan_Pos:Last_Chan_Pos) = &
                           Cloud_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE (2)  ! This is cloud fraction
                      RTModel_Jacobian(Prof_CloudCover, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           Cloud_Vars_K(Variable_Number, 1:NumInstChans)
                   CASE DEFAULT 
                      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                           'Not expecting more than two cloud variables'
                      ErrStatRep = StatusFatal
                      CALL Gen_ErrorReport( RoutineName,  &
                           ErrorMessage, &
                           ErrStatRep    )
                 END SELECT
              END DO
           END IF
           ! End of Jacobian assignments

        END DO Instrument_LoopK

        
        ! Pick out the Jacobian elements to be used by the minimisation
        !------------
        RT_Params % H_matrix_T = RTModel_Jacobian(Retrieved_Elements,:)

        !-------------------------------------------------------------
        ! Water Vapour Jacobians must be converted from 
        ! kg/kg to log(kg/kg)
        !-------------------------------------------------------------
        IF (Ret_FirstQ > 0) THEN
           DO I = Ret_FirstQ, Ret_LastQ
              IF (Retrieved_Elements(I) > 0) &
                   RT_Params % H_matrix_T(I,:) = &
                   RT_Params % H_matrix_T(I,:) * &
                   Atmos_Prof(Retrieved_Elements(I)-Prof_FirstQ+1,2,1)
           END DO
        END IF

        !-------------------------------------------------------------
        ! This is the surface humidity Jacobian
        !-------------------------------------------------------------
        IF (Ret_q2 > 0) THEN
           IF (Retrieved_Elements(Ret_q2) > 0) &
                RT_Params % H_matrix_T(Ret_q2,:) = &
                RT_Params % H_matrix_T(Ret_q2,:) * &
                Surf_Vars(2,1)
        END IF
        RT_Params % H_matrix = TRANSPOSE(RT_Params % H_matrix_T)

     CASE DEFAULT
        WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
             'Incorrect Fastmodel Mode in Fastmodel Interface'
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &
             ErrStatRep    )
   END SELECT
   
ENDIF RTTOV_FastmodelMode

End Subroutine NWPSAF_RTTOV7_Interface
