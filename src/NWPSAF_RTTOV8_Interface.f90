Subroutine NWPSAF_RTTOV8_Interface ( &
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
! Description: Interface between NWPSAF 1DVar code and RTTOV8
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     05/04/05 Original version - based on RTTOV7_Interface.  E. Pavelin.
! 1.1     08/12/06 Removed duplicate DebugMode reference. E. Pavelin.
! 1.2     04/01/07 Changed size of profiles_k to work with RTTOV8.7.
! 1.2              Removed tabs. E. Pavelin
! 1.3     04/01/07 Made more arrays allocatable. E. Pavelin.
! 1.4     08/01/07 Changed minimum soft limit for u & v to -100m/s. E. Pavelin.
!
! Code Description:
!   Language:		Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type, &
     MaxChanUsed

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     StatusWarning,      & 
     GeneralMode,        &
     OperationalMode,    &
     ProductionMode,     &
     DiagnosticMode,     &
     DebugMode,          &
     VerboseMode,        &
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
     Num_ProfElementsUsed, &
     Soft_Limits, &
     ! The following are RTTOV specific
     Num_Press_Levels,   &
     Num_Atm_Prof_Vars,  &
     Num_Surf_Vars,      &
     Num_Skin_Vars,      &
     Num_Cloud_Vars

USE NWPSAFMod_RTTOV8, ONLY: &
     RTTOV8_Coef

USE rttov_const, ONLY : &
     gas_id_watervapour, &
     gas_id_ozone,      &
     errorstatus_warning
     
USE rttov_types, ONLY : &
     profile_Type,      &
     radiance_Type,     &
     transmission_Type, &
     jprb

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_MessageReport.interface'

!----- The following are RTTOV8 specific
INCLUDE 'rttov_setup.interface'
INCLUDE 'rttov_setupchan.interface'
INCLUDE 'rttov_direct.interface'
INCLUDE 'rttov_k.interface'


!----- Subroutine arguments:
INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans   ! Instrument channel selection
INTEGER, INTENT(OUT)                    :: ErrorCode   ! Returned error code


!----- Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV8_Interface"

REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv
REAL, PARAMETER ::  O3_MassMixToPPMV  = 6.035045e5  ! kg/kg -> ppmv

INTEGER, PARAMETER :: Err_Unit = 0        ! RTTOV8 error unit (STDOUT)
INTEGER, PARAMETER :: jpnsat = 20         ! Formerly used by RTTOV7 as max sensors

!----- Local variables:

INTEGER :: I, IChan1, IChan2, Variable_Number
INTEGER :: ErrStatRep                 ! error status for Gen_ErrorReport
INTEGER :: First_Instr, Last_Instr
INTEGER :: IPrint = 0
INTEGER :: Instrument
INTEGER :: NumInstChans
INTEGER :: SatIndex

CHARACTER(LEN=80) :: ErrorMessage(4)  ! Message for Gen_ErrorReport
CHARACTER(LEN=80) :: Message(2)       ! Message for Gen_MessageReport

! Arrays used in interface with both RT Models
!-----------

INTEGER :: First_Chan_Pos !) Position of first and last channels for  
INTEGER :: Last_Chan_Pos  !) current instrument in RT_Params % TotalBTs
INTEGER :: Inst_Chans (UsedChans % NumChans)
REAL :: RTModel_Jacobian (ProfSize, UsedChans % NumChans)
 
REAL, POINTER :: RTProf(:)

! Arrays used in RTTOV interface (should be set to zero size if RTTOV is
! not to be used) 
!--------

Type(profile_Type) :: profiles(Num_Profs)
Type(profile_Type), ALLOCATABLE :: profiles_k(:)
REAL :: OCast_Rads(UsedChans % NumChans, 2*Num_Press_Levels+2) 
INTEGER :: Valid_Channels(MaxChanUsed,jpnsat)
INTEGER :: RTInstrument(3,RT_Params%Num_Instruments)

! Variables used in initialising RTTOV interface 
!--------

INTEGER :: Actual_Num_Instr
INTEGER, POINTER :: First_Instr_Chan
INTEGER, POINTER :: Num_Instr_Chan 
INTEGER :: Last_Instr_Chan 
INTEGER :: Verbosity_Level
LOGICAL :: Cloud_Switch

! Variables needed for RTTOV8 interface
!--------
INTEGER :: nfrequencies
INTEGER :: nchannels
INTEGER :: nbtout
INTEGER :: errorstatus(RT_Params%Num_Instruments)
INTEGER :: errorstatus2(Num_Profs)
INTEGER :: nchan(1)
INTEGER :: channels(UsedChans % NumChans)
REAL(kind=jprb), ALLOCATABLE :: surfem(:)           
LOGICAL, ALLOCATABLE :: calcemiss(:)
LOGICAL :: switchrad
Type(transmission_Type) :: transmission
Type(transmission_Type) :: transmission_k
Type(radiance_Type)     :: radiancedata
INTEGER :: lprofiles(UsedChans % NumChans)
INTEGER, ALLOCATABLE :: polarisations(:,:)
REAL(kind=jprb),ALLOCATABLE :: emissivity(:)
REAL(kind=jprb), ALLOCATABLE :: emissivity_k(:)
INTEGER :: wv_pos, o3_pos
!---------------------------------------------------------------

IF ( GeneralMode >= VerboseMode ) THEN
   WRITE(UNIT=Message(1),FMT='("RTTOV8 Called in MODE = ",i3)') &
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

! Start calls to RTTOV8

RTTOV_FastmodelMode : IF &
     ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
!----------------- RTTOV8 INITIALISATION SECTION ----------------------------
   IF ( GeneralMode >= VerboseMode ) THEN
      Message(1) = 'Initialising RTTOV8'
      CALL Gen_MessageReport( RoutineName,Message(1:1) )
   ENDIF

   Valid_Channels(:,:)=0

!----- The following block of code is almost obsolete now
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
    

!----- Extract the channel numbers for each instrument
!----- These are the channel numbers in the the RT_Params state vector etc.
!----- RT_Params % NumChans(:) contains the number of channels for each instrument.
   DO I = 1, RT_Params % Num_Instruments
      First_Instr_Chan => RT_Params % First_Channel_for_Instrument(I)
      Num_Instr_Chan   => RT_Params % NumChans(I)
      IF (Num_Instr_Chan > 0) Actual_Num_Instr = I
      Last_Instr_Chan  =  First_Instr_Chan + Num_Instr_Chan - 1
      Valid_Channels(1:Num_Instr_Chan,I) = RT_Params % &
           Absolute_Channel_Number(First_Instr_Chan:Last_Instr_Chan)
      IF (MaxChanUsed > Num_Instr_Chan) &
           Valid_Channels(Num_Instr_Chan+1:MaxChanUsed,I) = 0
   END DO

!----- Set remaining elements of Valid_Channels to zero
   IF (RT_Params % Num_Instruments < jpnsat) &
        Valid_Channels(:,RT_Params % Num_Instruments+1:jpnsat) = 0

   ! Reset number of instruments to prevent some unnecessary initialisations
   ! (this isn't the perfect solution yet)
   RT_Params % Num_Instruments = Actual_Num_Instr
   
!----- Set up RTTOV8 Instrument definition
   RTInstrument(1,:) = RT_Params % SeriesChoice(1:Actual_Num_Instr)
   RTInstrument(2,:) = RT_Params % PlatformChoice(1:Actual_Num_Instr)
   RTInstrument(3,:) = RT_Params % SubTypeChoice(1:Actual_Num_instr)

   !----- Translate verbosity mode into RTTOV8 Verbosity_Level
   SELECT CASE (GeneralMode)
     CASE (OperationalMode)
       Verbosity_Level = 1
     CASE (ProductionMode)
       Verbosity_Level = 1
     CASE (DiagnosticMode)
       Verbosity_Level = 2
     CASE (DebugMode)
       Verbosity_Level = 2
     CASE (VerboseMode)
       Verbosity_Level = 3
     CASE DEFAULT
       Verbosity_Level = 3
   END SELECT

   !----- Allocate space for RTTOV8 coefficients
   ALLOCATE( RTTOV8_Coef(RT_Params % Num_Instruments) )
  
   !----- General setup for RTTOV8.
   !----- This loads the coefficient files and sets up some basic parameters.
   CALL RTTOV_SETUP ( &
       errorstatus,                    & !out
       Err_Unit,                       & !in
       Verbosity_Level,                & !in
       RT_Params%Num_Instruments,      & !in
       RTTOV8_Coef(:),                 & !out
       RTInstrument,                   & !in
       Valid_Channels )                  !in

   ErrorCode = MAXVAL(errorstatus)

   IF (ErrorCode /= 0) THEN 
      ErrorMessage(1)='Error in RTTOV_SETUP'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
           'Error is Code ',ErrorCode
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF

   Num_Press_Levels = RTTOV8_Coef(1)%nlevels
   
   IF ( Num_Press_Levels /= Num_RTLevels ) THEN
      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
           'Mismatch in number of pressure levels on '// &
           'initialising RTTOV8'
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

   wv_pos = RTTOV8_Coef(1) % fmv_gas_pos( gas_id_watervapour )
   o3_pos = RTTOV8_Coef(1) % fmv_gas_pos( gas_id_ozone )

   ! Profile minima (include allowable tolerances)
   Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = RTTOV8_Coef(1)%lim_prfl_tmin(:)-0.5
   Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = &
        LOG(RTTOV8_Coef(1)%lim_prfl_gmin(:,wv_pos)/H2O_MassMixToPPMV)-0.5
   Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = &
        RTTOV8_Coef(1)%lim_prfl_gmin(:,o3_pos)*0.8
   Soft_Limits % Minimum(Prof_T2) = RTTOV8_Coef(1)%lim_prfl_tmin(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_Q2) = &
        LOG(RTTOV8_Coef(1)%lim_prfl_gmin(Num_RTLevels,wv_pos) / &
        H2O_MassMixToPPMV)-0.5
   Soft_Limits % Minimum(Prof_PStar) = 300.
   Soft_Limits % Minimum(Prof_UWind) = -100.
   Soft_Limits % Minimum(Prof_VWind) = -100.
   Soft_Limits % Minimum(Prof_TStar) = RTTOV8_Coef(1)%lim_prfl_tmin(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_CTP) = Cloud_Min_Pressure
   Soft_Limits % Minimum(Prof_CloudCover) = 0.

   ! Profile maxima
   Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = RTTOV8_Coef(1)%lim_prfl_tmax(:)+0.5
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
        LOG(RTTOV8_Coef(1)%lim_prfl_gmax(:,wv_pos)/H2O_MassMixToPPMV)+0.5
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
        RTTOV8_Coef(1)%lim_prfl_gmax(:,o3_pos)*1.2
   Soft_Limits % Maximum(Prof_T2) = RTTOV8_Coef(1)%lim_prfl_tmax(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_Q2) = &
        LOG(RTTOV8_Coef(1)%lim_prfl_gmax(Num_RTLevels,wv_pos) / &
        H2O_MassMixToPPMV)+0.5
   ! (This is the RTTOV hard limit:) 
   Soft_Limits % Maximum(Prof_PStar) = 1200.
   Soft_Limits % Maximum(Prof_UWind) = 100.
   Soft_Limits % Maximum(Prof_VWind) = 100.
   Soft_Limits % Maximum(Prof_TStar) = RTTOV8_Coef(1)%lim_prfl_tmax(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_CTP) = 1200.
   Soft_Limits % Maximum(Prof_CloudCover) = 1.0

   ! These are the RTTOV hard limits:
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
        MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),-2.99578)
   !Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
   !     MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),-0.9862)
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
        MIN(Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3),20.)
   Soft_Limits % Maximum(Prof_Q2) = MIN(Soft_Limits % Maximum(Prof_Q2),-2.99578)
   !Soft_Limits % Maximum(Prof_Q2) = MIN(Soft_Limits % Maximum(Prof_Q2),-0.9862)

   ! Check for profile mismatches
   DO I=1,ProfSize
      IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Mismatch in RTTOV8 soft limits'
         WRITE( UNIT=ErrorMessage(2),FMT=* )  I, &
              Soft_Limits % Minimum(I), &
              Soft_Limits % Maximum(I)
         ErrStatRep = StatusWarning
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END IF
   END DO

!------------- END OF INITIALISATION SECTION ------------------------

ELSE IF (FastModel_Mode == FastModelMode_CleanUp) THEN

   DEALLOCATE( RTTOV8_Coef )

ELSE   ! If not doing initialise or cleanup...

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

   !----- Allocate profile arrays
   ALLOCATE( profiles(1)%p(RTTOV8_Coef(1)%nlevels) )
   ALLOCATE( profiles(1)%t(RTTOV8_Coef(1)%nlevels) )
   ALLOCATE( profiles(1)%q(RTTOV8_Coef(1)%nlevels) )
   ALLOCATE( profiles(1)%o3(RTTOV8_Coef(1)%nlevels) )
   ALLOCATE( profiles(1)%clw(RTTOV8_Coef(1)%nlevels) )

   !  Unpack the RTProf Profile Vector into the input arrays used for 
   ! RTTOV:
   !------   

   profiles(1) % nlevels = Num_Press_Levels
   profiles(1) % ozone_data = .TRUE.
   profiles(1) % co2_data = .FALSE.
   profiles(1) % clw_data = .FALSE.
   profiles(1) % p(:) = RT_Params%Pressure_Pa(:)/100.0

   Num_Cloud_Vars    = 2  ! Set these values (used to be set by RTTVI)
   Num_Atm_Prof_Vars = 4
   Num_Surf_Vars     = 5 
   Num_Skin_Vars     = 1
      
   ! First do the atmospheric profile variables:
   DO Variable_Number = 1, Num_Atm_Prof_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is temperature
         profiles(1)%t(:) = RTProf(Prof_FirstT : Prof_LastT)
      CASE (2)  ! This is humidity (in ppmv) (convert from LOG kg/kg)
         profiles(1)%q(:) = EXP(RTProf(Prof_FirstQ:Prof_LastQ)) * H2O_MassMixToPPMV
      CASE (3)  ! This is Ozone (ppmv)
         profiles(1)%o3(:) = RTProf(Prof_FirstO3 : Prof_LastO3)
      CASE (4)  ! This is Cloud Liquid Water (-1.0 == disabled for now).
         profiles(1)%clw(:) = -1.0
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
         profiles(1)%s2m%t = RTProf(Prof_T2)
      CASE (2)  ! This is surface humidity (ppmv) (convert from LOG kg/kg)
         profiles(1)%s2m%q = EXP(RTProf(Prof_Q2)) * H2O_MassMixToPPMV
      CASE (3)  ! This is surface pressure
         profiles(1)%s2m%p = RTProf(Prof_PStar)
      CASE (4)  ! This is U wind
         profiles(1)%s2m%u = RTProf(Prof_UWind)
      CASE (5)  ! This is V wind
         profiles(1)%s2m%v = RTProf(Prof_VWind)
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
	 profiles(1)%skin%t = RTProf(Prof_TStar)
      CASE (2,3,4,5,6)  ! Fastem parameters.  Set to zero for now.
         profiles(1)%skin%fastem(:) = 0.0
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
         profiles(1)%ctp = RTProf(Prof_CTP)
         profiles(1)%ctp = MIN(profiles(1)%ctp,RTProf(Prof_PStar)*1._JPRB)
         profiles(1)%ctp = MAX(profiles(1)%ctp, &
              MINVAL(RT_Params % Pressure_Pa(:) / 100._JPRB))
      CASE (2)  ! This is cloud fraction (make sure is strictly within the 
                ! range 0-1 as otherwise RTTOV may crash!)
         profiles(1)%cfraction = RTProf(Prof_CloudCover)
         profiles(1)%cfraction = MIN(profiles(1)%cfraction,1.0_JPRB)
         profiles(1)%cfraction = MAX(profiles(1)%cfraction,0.0_JPRB)
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
      
   
   !! No cloud processing at present (calculates extra radiances)
   !Cloud_Switch = .FALSE.
   ! Switch on cloud processing (to get downcld)
   Cloud_Switch = .TRUE.
   
   profiles(1)%zenangle = RT_Params % SatZenithAngle
   !RTModel_SolarZenithAngle(1) = RT_Params % SolarZenAngle
   profiles(1)%skin%surftype = RT_Params % RTSurfaceType
        
   IF (GeneralMode >= DebugMode) IPrint = 1
   SELECT CASE ( Fastmodel_Mode )
     CASE( FastmodelMode_Forward )
        Instrument_Loop : DO Instrument = First_Instr, Last_Instr
	   ! Count no. of channels for current instrument
           NumInstChans = &
                COUNT(RT_Params % Instrument_Number( &
                SatIndex,UsedChans % Channels) == Instrument)

           nchan(:) = NumInstChans

           Call RTTOV_SETUPCHAN( &
               Num_Profs,                &  ! in 
	       nchan,                    &  ! in
	       RTTOV8_Coef(Instrument),  &  ! in
	       nfrequencies,             &  ! out
	       nchannels,                &  ! out
	       nbtout )                     ! out

           !----- Allocate radiance arrays
           !radiancedata%lcloud = .TRUE.
           ALLOCATE( radiancedata%clear(nchannels))
           ALLOCATE( radiancedata%clear_out(nchannels))
           ALLOCATE( radiancedata%cloudy(nchannels))
           ALLOCATE( radiancedata%total(nchannels))
           ALLOCATE( radiancedata%total_out(nchannels))
           ALLOCATE( radiancedata%out(nchannels))
           ALLOCATE( radiancedata%out_clear(nchannels))
           ALLOCATE( radiancedata%bt(nchannels))
           ALLOCATE( radiancedata%bt_clear(nchannels))
           ALLOCATE( radiancedata%upclear(nchannels))
           ALLOCATE( radiancedata%dnclear(nchannels))
           ALLOCATE( radiancedata%reflclear(nchannels))
           ALLOCATE( radiancedata%overcast(RTTOV8_Coef(1)%nlevels,nchannels))
           ALLOCATE( radiancedata%downcld(RTTOV8_Coef(1)%nlevels,nchannels))
           !----- Allocate transmission arrays
           ALLOCATE( transmission%tau_surf(nchannels) )
           ALLOCATE( transmission%tau_layer(RTTOV8_Coef(1)%nlevels,nchannels) )
           ALLOCATE( transmission%od_singlelayer(RTTOV8_Coef(1)%nlevels,nchannels) )

           ALLOCATE( polarisations(nchannels,3) )
           ALLOCATE( emissivity(nchannels) )

           ALLOCATE( surfem(nchannels) )
           ALLOCATE( calcemiss(nchannels) )
           ! Use default surface emissivities for now
           surfem(:) = 0.0
           calcemiss(:) = .TRUE.

            Call RTTOV_SETUPINDEX( &
               nchan,                     &  ! in
               Num_Profs,                 &  ! in
               nfrequencies,              &  ! in
               nchannels,                 &  ! in
               nbtout,                    &  ! in
               RTTOV8_Coef(Instrument),   &  ! in
               surfem(1:nchannels),       &  ! in
               lprofiles(1:NumInstChans),        &  ! out
               channels(1:NumInstChans),         &  ! out
               polarisations(1:nchannels, :),    &  ! out
               emissivity(1:nchannels) )            ! out

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

           !Surf_Emiss_IO(:) = 0.0
	   
           IF (First_Chan_Pos > 0) THEN

	      CALL rttov_direct( &
	           errorstatus2,                          & ! out
               nfrequencies,                          & ! in
               nchannels,                             & ! in
               nbtout,                                & ! in
               Num_Profs,                             & ! in
               Inst_Chans(1:NumInstChans),            & ! in
               polarisations(1:nchannels, :),         & ! in
               lprofiles(1:NumInstChans),             & ! in
               profiles(1),                           & ! in
               RTTOV8_Coef(Instrument),               & ! in
               Cloud_Switch,                          & ! in
               calcemiss(1:nchannels),                & ! in
               emissivity(1:nchannels),               & ! inout
               transmission,                          & ! out
               radiancedata )                           ! out

              ErrorCode=MAXVAL(errorstatus2)
              ! Assume soft limts have been set in this routine, so only
              ! report a catestrophic error!
              IF (ErrorCode <= errorstatus_warning) ErrorCode=0

              IF (ErrorCode > errorstatus_warning) THEN 
                 ErrorMessage(1)='Error in RTTOV_DIRECT'
                 WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                   'Error is Code ',ErrorCode
                 ErrStatRep = StatusFatal
                 CALL Gen_ErrorReport( RoutineName,  &
                    ErrorMessage, &
                    ErrStatRep    )
              END IF
 
              !  Unpack the required output arrays from RTTOV into the 
              !  RT_Param structure
              !------

              RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
	           radiancedata%total_out(1:NumInstChans) 
              RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                   radiancedata%out(1:NumInstChans)

              IF (CloudyRetrieval) THEN
                 RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      1:RTTOV8_Coef(1)%nlevels) = &
		      TRANSPOSE(radiancedata%overcast(:,1:NumInstChans))
                 RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      RTTOV8_Coef(1)%nlevels+1:2*RTTOV8_Coef(1)%nlevels) = &
		      TRANSPOSE(radiancedata%downcld(:,1:NumInstChans))
                 RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      2*RTTOV8_Coef(1)%nlevels+1) = radiancedata%upclear(1:NumInstChans)
                 RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      2*RTTOV8_Coef(1)%nlevels+2) = radiancedata%reflclear(1:NumInstChans)
		      
                 RT_Params % tc1(First_Chan_Pos:Last_Chan_Pos) = &
                      RTTOV8_Coef(Instrument) % &
                          ff_bco(Inst_Chans(1:NumInstChans))
                 RT_Params % tc2(First_Chan_Pos:Last_Chan_Pos) = &
                      RTTOV8_Coef(Instrument) % &
                          ff_bcs(Inst_Chans(1:NumInstChans))
                 RT_Params % bcon1(First_Chan_Pos:Last_Chan_Pos) = &
                      RTTOV8_Coef(Instrument) % fc_planck_c1 * &
		      RTTOV8_Coef(Instrument) % ff_cwn(UsedChans%Channels)**3
                 RT_Params % bcon2(First_Chan_Pos:Last_Chan_Pos) = &
                      RTTOV8_Coef(Instrument) % fc_planck_c2 * &
		      RTTOV8_Coef(Instrument) % ff_cwn(UsedChans%Channels)
              END IF
	      	      
           END IF
           !----- Deallocate radiance arrays
           DEALLOCATE( radiancedata%clear)
           DEALLOCATE( radiancedata%clear_out)
           DEALLOCATE( radiancedata%cloudy)
           DEALLOCATE( radiancedata%total)
           DEALLOCATE( radiancedata%total_out)
           DEALLOCATE( radiancedata%out)
           DEALLOCATE( radiancedata%out_clear)
           DEALLOCATE( radiancedata%bt)
           DEALLOCATE( radiancedata%bt_clear)
           DEALLOCATE( radiancedata%upclear)
           DEALLOCATE( radiancedata%dnclear)
           DEALLOCATE( radiancedata%reflclear)
           DEALLOCATE( radiancedata%overcast)
           DEALLOCATE( radiancedata%downcld)
           !----- Deallocate transmission arrays
           DEALLOCATE( transmission%tau_surf)
           DEALLOCATE( transmission%tau_layer)
           DEALLOCATE( transmission%od_singlelayer)
           DEALLOCATE( emissivity)
           DEALLOCATE( polarisations )
           DEALLOCATE( surfem )
           DEALLOCATE( calcemiss )	  	  
        END DO Instrument_Loop  

     CASE ( FastmodelMode_Gradient )
        Instrument_LoopK : DO Instrument = First_Instr, Last_Instr
           NumInstChans = &
                COUNT(RT_Params % Instrument_Number( &
                SatIndex,UsedChans % Channels) == Instrument)

           nchan(:) = NumInstChans
           Call RTTOV_SETUPCHAN( &
               Num_Profs,                &  ! in 
	       nchan,                    &  ! in
	       RTTOV8_Coef(Instrument),  &  ! in
	       nfrequencies,             &  ! out
	       nchannels,                &  ! out
	       nbtout )                     ! out
	  	  
           ALLOCATE( profiles_k(nchannels) )	  	  
           DO i=1,nchannels
             ALLOCATE( profiles_k(i)%p(RTTOV8_Coef(1)%nlevels) )
             ALLOCATE( profiles_k(i)%t(RTTOV8_Coef(1)%nlevels) )
             ALLOCATE( profiles_k(i)%q(RTTOV8_Coef(1)%nlevels) )
             ALLOCATE( profiles_k(i)%o3(RTTOV8_Coef(1)%nlevels) )
             ALLOCATE( profiles_k(i)%clw(RTTOV8_Coef(1)%nlevels) )
           ENDDO

           !----- Allocate radiance arrays
           !radiancedata%lcloud = .TRUE.
           ALLOCATE( radiancedata%clear(nchannels))
           ALLOCATE( radiancedata%clear_out(nchannels))
           ALLOCATE( radiancedata%cloudy(nchannels))
           ALLOCATE( radiancedata%total(nchannels))
           ALLOCATE( radiancedata%total_out(nchannels))
           ALLOCATE( radiancedata%out(nchannels))
           ALLOCATE( radiancedata%out_clear(nchannels))
           ALLOCATE( radiancedata%bt(nchannels))
           ALLOCATE( radiancedata%bt_clear(nchannels))
           ALLOCATE( radiancedata%upclear(nchannels))
           ALLOCATE( radiancedata%dnclear(nchannels))
           ALLOCATE( radiancedata%reflclear(nchannels))
           ALLOCATE( radiancedata%overcast(RTTOV8_Coef(1)%nlevels,nchannels))
           ALLOCATE( radiancedata%downcld(RTTOV8_Coef(1)%nlevels,nchannels))
           !----- Allocate transmission arrays
           ALLOCATE( transmission%tau_surf(nchannels) )
           ALLOCATE( transmission%tau_layer(RTTOV8_Coef(1)%nlevels,nchannels) )
           ALLOCATE( transmission%od_singlelayer(RTTOV8_Coef(1)%nlevels,nchannels) )
           ALLOCATE( transmission_k%tau_surf(nchannels) )
           ALLOCATE( transmission_k%tau_layer(RTTOV8_Coef(1)%nlevels,nchannels) )
           ALLOCATE( transmission_k%od_singlelayer(RTTOV8_Coef(1)%nlevels,nchannels) )
           ALLOCATE( polarisations(nchannels,3) )
           ALLOCATE( emissivity(nchannels) )
           ALLOCATE( emissivity_k(nchannels) )

           ALLOCATE( surfem(nchannels) )
           ALLOCATE( calcemiss(nchannels) )
           ! Use default surface emissivities for now
           surfem(:) = 0.0
           calcemiss(:) = .TRUE.

           Call RTTOV_SETUPINDEX( &
               nchan,                     &  ! in
               Num_Profs,                 &  ! in
               nfrequencies,              &  ! in
               nchannels,                 &  ! in
               nbtout,                    &  ! in
               RTTOV8_Coef(Instrument),   &  ! in
               surfem(1:nchannels),                    &  ! in
               lprofiles(1:NumInstChans),              &  ! out
               channels(1:NumInstChans),               &  ! out
               polarisations(1:nchannels, :),          &  ! out
               emissivity(1:nchannels) )                  ! out

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
              transmission_k % tau_surf(:)         = 0.0
              transmission_k % tau_layer(:,:)      = 0.0
              transmission_k % od_singlelayer(:,:) = 0.0
              emissivity_k(:)                      = 0.0   

              !Surf_Emiss_IO(1:NumInstChans) = 0.
              OCast_Rads(:,:) = 0.

              switchrad = .TRUE.    ! Select BT calculations
              CALL RTTOV_K( &
	            errorstatus2,                          & ! out
                nfrequencies,                          & ! in
                nchannels,                             & ! in
                nbtout,                                & ! in
                Num_Profs,                             & ! in
                Inst_Chans(1:NumInstChans),            & ! in
                polarisations(1:nchannels, :),         & ! in
                lprofiles(1:NumInstChans),             & ! in
                profiles(1),                           & ! in
                RTTOV8_Coef(Instrument),               & ! in
                Cloud_Switch,                          & ! in
                switchrad,                             & ! in
                calcemiss(1:nchannels),                & ! in
                emissivity(1:nchannels),               & ! inout
                profiles_k(1:nchannels),               & ! inout
                emissivity_k(1:nchannels),             & ! inout
                transmission,                          & ! inout
                transmission_k,                        & ! inout
                radiancedata )                           ! inout
		   
	      
              ErrorCode=MAXVAL(errorstatus2)
              ! Assume soft limits have been set in this routine, so only
              ! report a catestrophic error!
              IF (ErrorCode <= errorstatus_warning) ErrorCode=0

              IF (ErrorCode > errorstatus_warning) THEN 
                 ErrorMessage(1)='Error in RTTOV_K'
                 WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                   'Error is Code ',ErrorCode
                 ErrStatRep = StatusFatal
                 CALL Gen_ErrorReport( RoutineName,  &
                    ErrorMessage, &
                    ErrStatRep    )
              END IF

              !  Unpack the required output arrays from RTTOV into the 
              !  RT_Param structure
              !------
              RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
                   radiancedata%total_out(1:NumInstChans)
              RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                   radiancedata%out(1:NumInstChans)

              IF (CloudyRetrieval) THEN
                   RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      1:RTTOV8_Coef(1)%nlevels) = &
		      TRANSPOSE(radiancedata%overcast(:,1:NumInstChans))
                   RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      RTTOV8_Coef(1)%nlevels+1:2*RTTOV8_Coef(1)%nlevels) = &
		      TRANSPOSE(radiancedata%downcld(:,1:NumInstChans))
                   RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      2*RTTOV8_Coef(1)%nlevels+1) = radiancedata%upclear(1:NumInstChans)
                   RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                      2*RTTOV8_Coef(1)%nlevels+2) = &
		      radiancedata%reflclear(1:NumInstChans)
              ENDIF
 
              !  Unpack the output K arrays from RTTOV into RTModel_Jacobian:
              !------
              
              ! First do the atmospheric profile variables:
              DO Variable_Number = 1, Num_Atm_Prof_Vars
                 SELECT CASE (Variable_Number)
                   CASE (1)  ! This is temperature
		      DO i=First_Chan_Pos,Last_Chan_Pos
                        RTModel_Jacobian(Prof_FirstT:Prof_LastT, i) = &
			   profiles_k(i-First_Chan_Pos+1)%t(:)
		      ENDDO
                   CASE (2)  ! This is humidity
		      DO i=First_Chan_Pos,Last_Chan_Pos
                        RTModel_Jacobian(Prof_FirstQ:Prof_LastQ, i) = &
			   profiles_k(i-First_Chan_Pos+1)%q(:) !/ H2O_MassMixToPPMV
		      ENDDO
                   CASE (3)  ! This is Ozone
		      DO i=First_Chan_Pos,Last_Chan_Pos
                        RTModel_Jacobian(Prof_FirstO3:Prof_LastO3, i) = &
			   profiles_k(i-First_Chan_Pos+1)%o3(:)
		      ENDDO
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
                           profiles_k(1:NumInstChans)%s2m%t
                   CASE (2)  ! This is surface humidity (kg/kg)
                      RTModel_Jacobian(Prof_Q2, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           profiles_k(1:NumInstChans)%s2m%q !/ H2O_MassMixToPPMV
                   CASE (3)  ! This is surface pressure
                      RTModel_Jacobian(Prof_PStar, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           profiles_k(1:NumInstChans)%s2m%p
                   CASE (4)  ! This is U wind
                      RTModel_Jacobian(Prof_UWind, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           profiles_k(1:NumInstChans)%s2m%u
                   CASE (5)  ! This is V wind
                      RTModel_Jacobian(Prof_VWind, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           profiles_k(1:NumInstChans)%s2m%v
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
                           profiles_k(1:NumInstChans)%skin%t
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
                           profiles_k(1:NumInstChans)%ctp
                   CASE (2)  ! This is cloud fraction
                      RTModel_Jacobian(Prof_CloudCover, &
                           First_Chan_Pos:Last_Chan_Pos) = &
                           profiles_k(1:NumInstChans)%cfraction
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
           !----- Deallocate radiance arrays
           DEALLOCATE( radiancedata%clear)
           DEALLOCATE( radiancedata%clear_out)
           DEALLOCATE( radiancedata%cloudy)
           DEALLOCATE( radiancedata%total)
           DEALLOCATE( radiancedata%total_out)
           DEALLOCATE( radiancedata%out)
           DEALLOCATE( radiancedata%out_clear)
           DEALLOCATE( radiancedata%bt)
           DEALLOCATE( radiancedata%bt_clear)
           DEALLOCATE( radiancedata%upclear)
           DEALLOCATE( radiancedata%dnclear)
           DEALLOCATE( radiancedata%reflclear)
           DEALLOCATE( radiancedata%overcast)
           DEALLOCATE( radiancedata%downcld)
           !----- Deallocate transmission arrays
           DEALLOCATE( transmission%tau_surf)
           DEALLOCATE( transmission%tau_layer)
           DEALLOCATE( transmission%od_singlelayer)
           DEALLOCATE( transmission_k%tau_surf)
           DEALLOCATE( transmission_k%tau_layer )
           DEALLOCATE( transmission_k%od_singlelayer)
           DEALLOCATE( polarisations )
           DEALLOCATE( emissivity )
           DEALLOCATE( emissivity_k )
           DEALLOCATE( surfem )
           DEALLOCATE( calcemiss )
           DO i=1,nchannels
             DEALLOCATE( profiles_k(i)%p )
             DEALLOCATE( profiles_k(i)%t )
             DEALLOCATE( profiles_k(i)%q )
             DEALLOCATE( profiles_k(i)%o3 )
             DEALLOCATE( profiles_k(i)%clw )
           ENDDO
           DEALLOCATE( profiles_k )
        END DO Instrument_LoopK

        
        ! Pick out the Jacobian elements to be used by the minimisation
        !------------
        RT_Params % H_matrix_T = RTModel_Jacobian(Retrieved_Elements,:)

        !-------------------------------------------------------------
        ! Water Vapour Jacobians must be converted from 
        ! kg/kg to log(kg/kg)
        ! dy/d(ln q) = dy/dq * q/dy   (dy = 1K)
        !-------------------------------------------------------------

        IF (Ret_FirstQ > 0) THEN
           DO I = Ret_FirstQ, Ret_LastQ
              IF (Retrieved_Elements(I) > 0) &
                   RT_Params % H_matrix_T(I,:) = &
                   RT_Params % H_matrix_T(I,:) * &
                   profiles(1)%q(Retrieved_Elements(I)-Prof_FirstQ+1) !/ H2O_MassMixToPPMV

           END DO
        END IF

        !-------------------------------------------------------------
        ! This is the surface humidity Jacobian
        !-------------------------------------------------------------
        IF (Ret_q2 > 0) THEN
           IF (Retrieved_Elements(Ret_q2) > 0) &
                RT_Params % H_matrix_T(Ret_q2,:) = &
                RT_Params % H_matrix_T(Ret_q2,:) * &
                profiles(1)%s2m%q !/ H2O_MassMixToPPMV
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


   IF( (FastModel_Mode .eq. FastModelMode_Forward) .OR. &
       (FastModel_Mode .eq. FastModelMode_Gradient) ) THEN
      !----- Deallocate profile arrays
      DEALLOCATE( profiles(1)%p )
      DEALLOCATE( profiles(1)%t )
      DEALLOCATE( profiles(1)%q )
      DEALLOCATE( profiles(1)%o3 )
      DEALLOCATE( profiles(1)%clw )
   END IF
   
ENDIF RTTOV_FastmodelMode

End Subroutine NWPSAF_RTTOV8_Interface
