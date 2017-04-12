Subroutine NWPSAF_RTTOV11_Interface ( &
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
! Description: Interface between NWPSAF 1DVar code and RTTOV10
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 31      18/01/13 Original version, based on RTTOV10_interface P. Weston.
! 32      28/01/14 Changes to allow for HDF5 coefficients and channel 
!                  selections to be used. P. Weston
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
     OperationalMode,    &
     ProductionMode,     &
     DiagnosticMode,     &
     DebugMode,          &
     Retrieved_Elements, &
     Ret_FirstQ,         &
     Ret_LastQ,          &
     Ret_q2,             &
     CloudyRetrieval,    &
     Cloud_Min_Pressure, &
     RTsea,              &
     Use_EmisAtlas,      &
     Atlas_Dir,          &
     Atlas_Ver,          &
     MwClwRetrieval,     &
     Lqtotal,            &
     Ozone_Present

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
     Prof_FirstCLW, &
     Prof_LastT, &
     Prof_LastQ, &
     Prof_LastO3, &
     Prof_LastCLW, &
     Prof_T2, &
     Prof_q2, &
     Prof_Tstar, &
     Prof_pstar, &
     Prof_uwind, &
     Prof_vwind, &
     Prof_CTP, &
     Prof_CloudCover, &
     Prof_CLW,    &
     Num_RTlevels, &
     Num_Profs, &
     Soft_Limits, &
     UseModelLevels, &
     ! The following are RTTOV specific
     Num_Press_Levels,   &
     Num_Atm_Prof_Vars,  &
     Num_Surf_Vars,      &
     Num_Skin_Vars,      &
     Num_Cloud_Vars

USE NWPSAFMod_RTTOV11, ONLY: &
     RTTOV11_Coefs,            &
     RTTOV11_Opts,             &
     RTTOV11_Chanprof,         &
     jpch,                     &
     jpnsat

USE NWPSAFMod_Constants, ONLY :  & 
     epsilon_1000

USE rttov_const, ONLY :   &
     gas_id_watervapour,  &
     gas_id_ozone,        &
     errorstatus_fatal,   &
     sensor_id_mw,        &
     sensor_id_po

USE rttov_types, ONLY: &
     profile_type,      &
     radiance_type,     &
     transmission_type, &
     rttov_emissivity

#ifdef _RTTOV_HDF
use rttov_hdf_mod
#endif

Use parkind1, Only : jpim,jprb

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_MessageReport.interface'
INCLUDE 'rttov_alloc_rad.interface'
INCLUDE 'rttov_alloc_transmission.interface'
#include <rttov_setup_emis_atlas.interface>
#include <rttov_dealloc_coefs.interface>
INCLUDE 'rttov_deallocate_emis_atlas.interface'
INCLUDE 'rttov_direct.interface'
INCLUDE 'rttov_get_emis.interface'
INCLUDE 'rttov_k.interface'
#include <rttov_read_coefs.interface>
#include <rttov_init_coefs.interface>
INCLUDE 'NWPSAF_Layers_to_LWP.interface'

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrorCode


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV11_Interface"

REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv
REAL, PARAMETER ::  O3_MassMixToPPMV  = 6.035045e5  ! kg/kg -> ppmv

! Local variables:

INTEGER :: I, IChan1, IChan2, Variable_Number, J
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
INTEGER :: KPROF(UsedChans % NumChans)                        
INTEGER :: Inst_Chans (UsedChans % NumChans)
REAL :: RTModel_Jacobian (ProfSize, UsedChans % NumChans)
TYPE(rttov_emissivity) :: Surf_Emiss(UsedChans % NumChans * Num_Profs) 
TYPE(rttov_emissivity) :: Surf_Emiss_K(UsedChans % NumChans * Num_Profs) 
 
REAL, POINTER :: RTProf(:)

! Arrays used in RTTOV interface (should be set to zero size if RTTOV is
! not to be used) 
!--------

TYPE(profile_type)  :: Profiles(Num_Profs)
TYPE(profile_type)  :: Profiles_K(UsedChans % NumChans * Num_Profs)
TYPE(radiance_type) :: Radiance
TYPE(radiance_type) :: Radiance_K
REAL :: OCast_Rads(UsedChans % NumChans, 2*Num_Press_Levels+2) 
INTEGER(kind=jpim) :: Valid_Channels(jpch,jpnsat)
INTEGER(kind=jpim) :: RTInstrument(3,RT_Params%Num_Instruments)
INTEGER :: wv_pos, o3_pos
LOGICAL :: calcemiss(UsedChans % NumChans * Num_Profs)
Type(transmission_type) :: transmission
Type(transmission_type) :: transmission_K
INTEGER(kind=jpim) :: ErrorStatus(RT_Params%Num_Instruments)

! Variables used in initialising RTTOV interface 
!--------

INTEGER :: Actual_Num_Instr
INTEGER, POINTER :: First_Instr_Chan
INTEGER, POINTER :: Num_Instr_Chan 
INTEGER :: Last_Instr_Chan 
INTEGER :: Nchans(RT_Params%Num_Instruments)
LOGICAL :: Cloud_Switch
INTEGER, PARAMETER :: ASW_ALLOCATE = 1
INTEGER, PARAMETER :: ASW_DEALLOCATE = 0
LOGICAL            :: verbosity

! Variables used for cloud liquid water retrieval
REAL :: LWP_calc
REAL :: cloud_structure(Num_RTlevels)
INTEGER         :: QtotalOption
INTEGER         :: Num_WetLevels ! Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
REAL            :: wt(Ret_LastQ - Ret_FirstQ + 1) 
REAL            :: wqtotal(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: wq(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: wql(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: wpress(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: esat(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: SESAT
REAL            :: SDlnes_DT
REAL            :: Dlnes_DT(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: qsat(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: epsilon = epsilon_1000/1000.0
REAL            :: Dqsat_dT(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: delta_eps
REAL            :: RTModel_Jacobian_Add(ProfSize, UsedChans % NumChans)
REAL            :: clw_profiles(Num_RTlevels)
!---------------------------------------------------------------

IF ( GeneralMode >= VerboseMode ) THEN
   WRITE(UNIT=Message(1),FMT='("RTTOV11 Called in MODE = ",i3)') &
        Fastmodel_Mode
   CALL Gen_MessageReport( RoutineName,Message(1:1) )
END IF

! Error code set to zero here

ErrorStatus(:) = 0
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

! Start calls to RTTOV11

! ****** INITIALISE ******
RTTOV_FastmodelMode : IF &
     ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
   IF ( GeneralMode >= VerboseMode ) THEN
      Message(1) = 'Initialising RTTOV11'
      CALL Gen_MessageReport( RoutineName,Message(1:1) )
   ENDIF

   Valid_Channels(:,:)=0

! Check that there aren't too many instruments
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
      ! Note: Channel numbers will be zero if not defined in obsfile.
      ! In this case need to not pass Valid_Channels to RTTOV_READ_COEFS
      Nchans(I) = 0
      DO J = 1, Num_Instr_Chan
         IF( Valid_Channels(J,I) > 0 ) THEN
            Nchans(I) = Nchans(I) + 1
         END IF
      END DO
   END DO

   IF (RT_Params % Num_Instruments < jpnsat) &
        Valid_Channels(:,RT_Params % Num_Instruments+1:jpnsat) = 0

   ! Reset number of instruments to prevent some unnecessary initialisations
   RT_Params % Num_Instruments = Actual_Num_Instr

   !----- Set up RTTOV11 Instrument definition
   RTInstrument(1,:) = RT_Params % SeriesChoice(1:Actual_Num_Instr)
   RTInstrument(2,:) = RT_Params % PlatformChoice(1:Actual_Num_Instr)
   RTInstrument(3,:) = RT_Params % SubTypeChoice(1:Actual_Num_instr)

   !----- Translate verbosity mode into RTTOV11 Verbosity flag
   SELECT CASE (GeneralMode)
     CASE (OperationalMode)
       Verbosity = .False.
     CASE (ProductionMode)
       Verbosity = .False.
     CASE (DiagnosticMode)
       Verbosity = .True.
     CASE (DebugMode)
       Verbosity = .True.
     CASE (VerboseMode)
       Verbosity = .True.
     CASE DEFAULT
       Verbosity = .True.
   END SELECT

   Allocate(RTTOV11_Coefs(RT_Params % Num_Instruments))
   Allocate(RTTOV11_Opts(RT_Params % Num_Instruments))

   ! Initialise options structure

   RTTOV11_Opts(:) % config % do_checkinput     = .true.
   RTTOV11_Opts(:) % config % apply_reg_limits  = .false.
   RTTOV11_Opts(:) % config % verbose           = Verbosity
   RTTOV11_Opts(:) % rt_all % switchrad         = .true.
   RTTOV11_Opts(:) % rt_all % addrefrac         = .false.
   RTTOV11_Opts(:) % rt_all % use_q2m           = .true.
   RTTOV11_Opts(:) % rt_mw % do_lambertian      = .false.
   RTTOV11_Opts(:) % rt_mw % fastem_version     = 5
   IF (MwClwRetrieval) THEN
     RTTOV11_Opts(:) % rt_mw % clw_data         = .true.
   ELSE
     RTTOV11_Opts(:) % rt_mw % clw_data         = .false.
   END IF
   RTTOV11_Opts(:) % rt_ir % addsolar           = .false.
   RTTOV11_Opts(:) % rt_ir % do_nlte_correction = .false.
   RTTOV11_Opts(:) % rt_ir % addaerosl          = .false.
   RTTOV11_Opts(:) % rt_ir % addclouds          = .false.
   RTTOV11_Opts(:) % rt_ir % user_aer_opt_param = .false.
   RTTOV11_Opts(:) % rt_ir % user_cld_opt_param = .false.
   RTTOV11_Opts(:) % rt_ir % cldstr_threshold   = 0.001_jprb
   RTTOV11_Opts(:) % rt_ir % ozone_data         = .false.
   RTTOV11_Opts(:) % rt_ir % co2_data           = .false.
   RTTOV11_Opts(:) % rt_ir % n2o_data           = .false.
   RTTOV11_Opts(:) % rt_ir % co_data            = .false.
   RTTOV11_Opts(:) % rt_ir % ch4_data           = .false.
   RTTOV11_Opts(:) % rt_ir % pc % addpc         = .false.
   RTTOV11_Opts(:) % rt_ir % pc % addradrec     = .false.
   RTTOV11_Opts(:) % rt_ir % pc % ipcbnd        = -1_jpim
   RTTOV11_Opts(:) % rt_ir % pc % ipcreg        = -1_jpim
   RTTOV11_Opts(:) % interpolation % addinterp  = .false.
   RTTOV11_Opts(:) % interpolation % lgradp     = .false.
   RTTOV11_Opts(:) % interpolation % spacetop   = .true.

#ifdef _RTTOV_HDF
   Call open_hdf(.false., ErrorCode)
   IF (ErrorCode /= 0) THEN 
      ErrorMessage(1)='Error in OPEN_HDF'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorCode
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
#endif

   DO I=1,RT_Params%Num_Instruments
     ! Read in RTTOV 11 coefficients
     IF (Nchans(I) > 0) THEN
        CALL RTTOV_READ_COEFS( &
             ErrorStatus(I),                           & ! out
             RTTOV11_Coefs(I),                         & ! out
             RTTOV11_Opts(I),                          & ! in
             channels = Valid_Channels(1:Nchans(I),I), & ! in
             instrument = RTInstrument(1:Nchans(I),I))   ! in
     ELSE
        CALL RTTOV_READ_COEFS( &
             ErrorStatus(I),                           & ! out
             RTTOV11_Coefs(I),                         & ! out
             RTTOV11_Opts(I),                          & ! in
             instrument = RTInstrument(1:Nchans(I),I))   ! in
     END IF

     ErrorCode = ErrorStatus(I)
     IF (ErrorCode /= 0) THEN 
        ErrorMessage(1)='Error in RTTOV_READ_COEFS'
        WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
             'Error is Code ',ErrorCode
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &
             ErrStatRep    )
     END IF

     Call RTTOV_INIT_COEFS( &
          ErrorStatus(I),   & ! out
          RTTOV11_Opts(I),  & ! in
          RTTOV11_Coefs(I)  ) ! inout

     ErrorCode = ErrorStatus(I)
     IF (ErrorCode /= 0) THEN 
        ErrorMessage(1)='Error in RTTOV_INIT_COEFS'
        WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
             'Error is Code ',ErrorCode
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &
             ErrStatRep    )
     END IF

   END DO   

#ifdef _RTTOV_HDF
   call close_hdf(ErrorCode)
   IF (ErrorCode /= 0) THEN 
      ErrorMessage(1)='Error in OPEN_HDF'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
           'Error is Code ',ErrorCode
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   END IF
#endif

   ! Determine whether ozone is present in the background file and coefficients
   IF (Ozone_Present .AND. &
      ( ANY(RTTOV11_Coefs(1) % coef % fmv_gas_id(:) == gas_id_ozone))) THEN
     RTTOV11_Opts(:) % rt_ir % ozone_data = .true.
   ELSE
     RTTOV11_Opts(:) % rt_ir % ozone_data = .false.
   END IF

   ! Determine whether pressure level interpolation is needed for
   ! each instrument.
   IF ( ANY(RTTOV11_Coefs(:)%coef%nlevels /= Num_RTLevels) ) &
     UseModelLevels=.TRUE.
   DO i=1,RT_Params%Num_Instruments
     IF(RTTOV11_Coefs(i)%coef%nlevels /= Num_RTLevels) THEN
       RTTOV11_Opts(i) % interpolation % addinterp = .TRUE.
       RTTOV11_Opts(i) % config % apply_reg_limits = .TRUE.
       RTTOV11_Opts(i) % interpolation % lgradp = .TRUE.
       IF(GeneralMode >= DebugMode) Then
         Write(*,*)'N.B.: Using profile interpolation for instrument ',i
       END IF
     END IF
   END DO
    
   ! Set up soft limits

   ALLOCATE( Soft_Limits % Minimum(ProfSize))
   ALLOCATE( Soft_Limits % Maximum(ProfSize))

   ! Default values:
   Soft_Limits % Minimum(:) = 0.
   Soft_Limits % Maximum(:) = 1.e10

   ! Default values:
   Soft_Limits % Minimum(:) = 0.
   Soft_Limits % Maximum(:) = 1.e10

   wv_pos = RTTOV11_Coefs(1)%coef % fmv_gas_pos( gas_id_watervapour )
   IF(RTTOV11_Opts(1) % rt_ir % ozone_data) o3_pos = &
        RTTOV11_Coefs(1) % coef % fmv_gas_pos( gas_id_ozone )

   ! Profile minima (include allowable tolerances)

   IF ( .NOT. UseModelLevels ) THEN

     ! Mininum temperature
     Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmin(:)-0.5
     ! Minimum water vapour
     Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = &
          LOG(RTTOV11_Coefs(1) % coef % lim_prfl_gmin(:,wv_pos)/ &
          H2O_MassMixToPPMV)-0.5
     ! Minimum ozone (if present)
     IF(RTTOV11_Opts(1) % rt_ir % ozone_data) THEN
        Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = &
             RTTOV11_Coefs(1) % coef % lim_prfl_gmin(:,o3_pos)*0.8
     Else
        Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = -9999.99
     END IF   
     ! Minimum surface variables
     Soft_Limits % Minimum(Prof_T2) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmin(Num_RTLevels)-0.5
     Soft_Limits % Minimum(Prof_Q2) = &
          LOG(RTTOV11_Coefs(1) % coef % lim_prfl_gmin(Num_RTLevels,wv_pos)/ &
          H2O_MassMixToPPMV)-0.5
     Soft_Limits % Minimum(Prof_PStar) = 300.
     Soft_Limits % Minimum(Prof_UWind) = -100.
     Soft_Limits % Minimum(Prof_VWind) = -100.
     Soft_Limits % Minimum(Prof_TStar) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmin(Num_RTLevels)-0.5
     Soft_Limits % Minimum(Prof_CTP) = Cloud_Min_Pressure
     Soft_Limits % Minimum(Prof_CloudCover) = 0.

     ! Profile maxima

     ! Maximum temperature
     Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmax(:)+0.5
     ! Maximum water vapour
     Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
          LOG(RTTOV11_Coefs(1) % coef % lim_prfl_gmax(:,wv_pos)/ &
          H2O_MassMixToPPMV)+0.5
     ! Maximum ozone (if present)
     IF(RTTOV11_Opts(1) % rt_ir % ozone_data) Then
        Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
             RTTOV11_Coefs(1) % coef % lim_prfl_gmax(:,o3_pos)*1.2
     Else
        Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = 99999.99
     END IF   
     ! Maximum surface variables
     Soft_Limits % Maximum(Prof_T2) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmax(Num_RTLevels)+0.5
     Soft_Limits % Maximum(Prof_Q2) = &
          LOG(RTTOV11_Coefs(1) % coef % lim_prfl_gmax(Num_RTLevels,wv_pos)/ &
          H2O_MassMixToPPMV)+0.5
     ! (This is the RTTOV hard limit:) 
     Soft_Limits % Maximum(Prof_PStar) = 1200.
     Soft_Limits % Maximum(Prof_UWind) = 100.
     Soft_Limits % Maximum(Prof_VWind) = 100.
     Soft_Limits % Maximum(Prof_TStar) = &
          RTTOV11_Coefs(1) % coef % lim_prfl_tmax(Num_RTLevels)+0.5
     Soft_Limits % Maximum(Prof_CTP) = 1200.
     Soft_Limits % Maximum(Prof_CloudCover) = 1.0

     ! These are the RTTOV hard limits:
     Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
          MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),-2.99578)

     Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
          MIN(Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3),20.)

     Soft_Limits % Maximum(Prof_Q2) = &
          MIN(Soft_Limits % Maximum(Prof_Q2),-2.99578)


     ! Check for profile mismatches
     DO I=1,ProfSize
        IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN 
           WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                'Mismatch in RTTOV11 soft limits'
           WRITE( UNIT=ErrorMessage(2),FMT=* )  I, &
                Soft_Limits % Minimum(I), &
                Soft_Limits % Maximum(I)
           ErrStatRep = StatusWarning
           CALL Gen_ErrorReport( RoutineName,  &
                ErrorMessage, &
                ErrStatRep    )
        END IF
     END DO

   END IF

ELSE IF (FastModel_Mode == FastModelMode_CleanUp) THEN

   DO I=1,RT_Params % Num_Instruments
     Call RTTOV_DEALLOC_COEFS( ErrorCode,    & ! out
                               RTTOV11_Coefs(I))  ! in

     IF (ErrorCode /= 0) THEN 
        ErrorMessage(1)='Error in RTTOV_DEALLOC_COEFS'
        WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
             'Error is Code ',ErrorCode
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &
             ErrStatRep    )
     END IF
   END DO

   IF (ASSOCIATED(RTTOV11_Coefs)) Deallocate(RTTOV11_Coefs)
   IF (ASSOCIATED(RTTOV11_Opts)) Deallocate(RTTOV11_Opts)
   IF (ASSOCIATED(RTTOV11_Chanprof)) Deallocate(RTTOV11_Chanprof)

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

   ! Allocate profile arrays (Don't use RTTOV_ALLOCPROF)
   Allocate( Profiles(1) % p(Num_RTLevels) )
   Allocate( Profiles(1) % t(Num_RTLevels) )
   Allocate( Profiles(1) % q(Num_RTLevels) )
   Allocate( Profiles(1) % o3(Num_RTLevels) )
   Allocate( Profiles(1) % clw(Num_RTLevels) )
   Profiles(1) % Be = 0    ! These need to be initialised for RTTOV to work
   Profiles(1) % cosbk = 0

   !  Unpack the RTProf Profile Vector into the input arrays used for 
   ! RTTOV:
   !------   
   profiles(1) % nlevels = Num_RTLevels
   profiles(1) % nlayers = Num_RTLevels - 1
   profiles(1) % idg = 1 ! Ice scheme - not used
   profiles(1) % ish = 1 ! Crystal shape - not used

   profiles(1) % p(:) = RT_Params%Pressure_Pa(:)/100.0

   Num_Cloud_Vars    = 2  ! Set these values (used to be set by RTTOV7)
   Num_Atm_Prof_Vars = 4
   Num_Surf_Vars     = 5 
   Num_Skin_Vars     = 1

   ! First do the atmospheric profile variables:
   DO Variable_Number = 1, Num_Atm_Prof_Vars
      SELECT CASE (Variable_Number)
      CASE (1)  ! This is temperature
         profiles(1)%t(:) = RTProf(Prof_FirstT : Prof_LastT)
      CASE (2)  ! This is humidity (in ppmv) (convert from LOG kg/kg)
         profiles(1)%q(:) = EXP(RTProf(Prof_FirstQ:Prof_LastQ)) * &
              H2O_MassMixToPPMV
      CASE (3)  ! This is Ozone (ppmv)
         IF(RTTOV11_Opts(1) % rt_ir % ozone_data) Then
            profiles(1)%o3(:) = RTProf(Prof_FirstO3 : Prof_LastO3)
         End If
      CASE (4)  ! This is Cloud Liquid Water 
         IF (MwClwRetrieval) THEN
            profiles(1)%clw(:) = RTProf(Prof_FirstCLW : Prof_LastCLW)
         ELSE
            profiles(1)%clw(:)=-1 !(-1.0 == disabled for now).
         END IF
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

   profiles(1)%s2m%wfetc = 100000.0 ! Wind fetch
   profiles(1)%skin%fastem(:) = 0.0 ! Initialise FASTEM coefs to zero

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
              MINVAL(RT_Params % Pressure_Pa(:) / 100.*1._JPRB))
      CASE (2)  ! This is cloud fraction (make sure is strictly within the 
                ! range 0-1 as otherwise RTTOV may crash!)
         profiles(1)%cfraction = RTProf(Prof_CloudCover)
         profiles(1)%cfraction = MIN(profiles(1)%cfraction,1._JPRB)
         profiles(1)%cfraction = MAX(profiles(1)%cfraction,0._JPRB)
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

   ! Don't do full cloudy calculation
   Cloud_Switch = .FALSE.
   
   Profiles(1)%zenangle = RT_Params % SatZenithAngle
   Profiles(1)%azangle = RT_Params % SatAzimAngle
   Profiles(1)%sunzenangle = RT_Params % SolarZenAngle
   Profiles(1)%sunazangle = RT_Params % SatSolarAzimAngle
   Profiles(1)%date = RT_Params % Date
   Profiles(1)%elevation = RT_Params % Elevation / 1000._JPRB
   Profiles(1)%latitude = RT_Params % Latitude
   Profiles(1)%longitude = RT_Params % Longitude
   Profiles(1)%skin%surftype = RT_Params % RTSurfaceType
   Profiles(1)%skin%watertype = 1 ! Default to ocean water for now
   Profiles(1)%snow_frac = 0.0 ! Default to no snow for now
   KPROF(1:UsedChans % NumChans) = 1
   
   IF (GeneralMode >= DebugMode) IPrint = 1
   SELECT CASE ( Fastmodel_Mode )
! ****** FORWARD MODEL CALL ******
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

           IF (First_Chan_Pos == 0) CYCLE

           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance,                   & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_ALLOCATE,               & ! in
                                 init=.TRUE.                 ) ! in

           IF (ErrorCode /= 0) THEN 
              ErrorMessage(1)='Error in RTTOV_ALLOC_RAD'
              WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                   'Error is Code ',ErrorCode
              ErrStatRep = StatusFatal
              CALL Gen_ErrorReport( RoutineName,  &
                   ErrorMessage, &
                   ErrStatRep    )
           END IF

           !----- Allocate transmission arrays
           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission, &
                profiles(1) % nlayers, &
                NumInstChans, &
                ASW_ALLOCATE, &
                init = .TRUE. )

           Allocate(RTTOV11_Chanprof(NumInstChans))
           RTTOV11_Chanprof(1:NumInstChans)%chan = &
                 Inst_Chans(1:NumInstChans)
           RTTOV11_Chanprof(1:NumInstChans)%prof = KPROF*1_JPIM

           ! Initialise emissivity arrays
           IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN

             IF (RTTOV11_Coefs(Instrument)%coef%id_sensor == sensor_id_mw &
                 .OR. RTTOV11_Coefs(Instrument)%coef%id_sensor == &
                 sensor_id_po)  THEN ! MW
               CALL rttov_setup_emis_atlas( &
                    ErrorCode,                      & ! out
                    RTTOV11_Opts(Instrument),       & ! in
                    Profiles(1) % Date(2),          & ! in
                    RTTOV11_Coefs(Instrument),      & ! in
                    path=Atlas_Dir,                 & ! in
                    mw_atlas_ver = Atlas_Ver        ) ! in
             ELSE ! IR
               CALL rttov_setup_emis_atlas( &
                    ErrorCode,                      & ! out
                    RTTOV11_Opts(Instrument),       & ! in
                    Profiles(1) % Date(2),          & ! in
                    RTTOV11_Coefs(Instrument),      & ! in
                    path=Atlas_Dir,                 & ! in
                    ir_atlas_read_std=.TRUE.        ) ! in
             END IF

             IF (ErrorCode >= errorstatus_fatal) THEN 
               ErrorMessage(1)='Error in rttov_atlas_setup'
               WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                 'Error is Code ',ErrorCode
               ErrStatRep = StatusFatal
               CALL Gen_ErrorReport( RoutineName,  &
                    ErrorMessage, &
                    ErrStatRep    )
             END IF                

             Surf_Emiss(1:NumInstChans)%emis_in = 0.
             Surf_Emiss(1:NumInstChans)%emis_out = 0.

             CALL rttov_get_emis( &
                  ErrorCode,                      & ! out
                  RTTOV11_Opts(Instrument),       & ! in
                  RTTOV11_Chanprof,               & ! in
                  Profiles,                       & ! in
                  RTTOV11_Coefs(Instrument),      & ! in
                  Emissivity = Surf_Emiss%emis_in ) ! out

             IF (ErrorCode >= errorstatus_fatal) THEN 
               ErrorMessage(1)='Error in rttov_get_emis'
               WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                 'Error is Code ',ErrorCode
               ErrStatRep = StatusFatal
               CALL Gen_ErrorReport( RoutineName,  &
                 ErrorMessage, &
                 ErrStatRep    )
             END IF

             Calcemiss(1:NumInstChans) = .FALSE.

           ELSE

             Surf_Emiss(1:NumInstChans)%emis_in = 0.
             Surf_Emiss(1:NumInstChans)%emis_out = 0.
             Calcemiss(1:NumInstChans) = .TRUE.

           END IF

           CALL RTTOV_DIRECT( &
                ErrorStatus(Instrument),   & ! out
                RTTOV11_Chanprof,          & ! in
                RTTOV11_Opts(Instrument),  & ! in
                Profiles,                  & ! in
                RTTOV11_Coefs(Instrument), & ! in
                transmission,              & ! inout
                Radiance,                  & ! inout
                calcemis=calcemiss,        & ! in
                emissivity=Surf_Emiss )      ! inout

           ErrorCode = ErrorStatus(Instrument)
           IF (ErrorCode >= errorstatus_fatal) THEN 
              ErrorMessage(1)='Error in RTTOV_DIRECT'
              WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                'Error is Code ',ErrorCode
              ErrStatRep = StatusFatal
              CALL Gen_ErrorReport( RoutineName,  &
                 ErrorMessage, &
                 ErrStatRep    )
           END IF

           Deallocate(RTTOV11_Chanprof)

           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission, &
                profiles(1) % nlayers, &
                NumInstChans,&
                ASW_DEALLOCATE)


           !  Unpack the required output arrays from RTTOV into the 
           !  RT_Param structure
           !------

           RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
                Radiance%total(1:NumInstChans) 
           RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                Radiance%bt(1:NumInstChans)

           IF (CloudyRetrieval) THEN
              RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                   1:profiles(1) % nlayers) = &
	           TRANSPOSE(Radiance%overcast(:,1:NumInstChans))
              RT_Params % tc1(First_Chan_Pos:Last_Chan_Pos) = &
                   RTTOV11_Coefs(Instrument)%coef % &
                   ff_bco(Inst_Chans(1:NumInstChans))
              RT_Params % tc2(First_Chan_Pos:Last_Chan_Pos) = &
                   RTTOV11_Coefs(Instrument)%coef % &
                   ff_bcs(Inst_Chans(1:NumInstChans))
              RT_Params % bcon1(First_Chan_Pos:Last_Chan_Pos) = &
                   RTTOV11_Coefs(Instrument)%coef % fc_planck_c1 * &
                   RTTOV11_Coefs(Instrument)%coef % &
                   ff_cwn(UsedChans%Channels)**3
              RT_Params % bcon2(First_Chan_Pos:Last_Chan_Pos) = &
                   RTTOV11_Coefs(Instrument)%coef % fc_planck_c2 * &
                   RTTOV11_Coefs(Instrument)%coef % &
                   ff_cwn(UsedChans%Channels)
           END IF

           RT_Params % RTEmissivity = Surf_Emiss%emis_out

           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance,                   & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_DEALLOCATE   )            ! in

           IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN
             CALL rttov_deallocate_emis_atlas( RTTOV11_Coefs(Instrument) )
           END IF

        END DO Instrument_Loop  

        ! Deallocate profile arrays
        Deallocate( Profiles(1) % p )
        Deallocate( Profiles(1) % t )
        Deallocate( Profiles(1) % q )
        Deallocate( Profiles(1) % o3 )
        Deallocate( Profiles(1) % clw )

     CASE ( FastmodelMode_Gradient )
!****** GRADIENT CALL ******
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

           IF (First_Chan_Pos == 0) CYCLE

           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance,                   & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_ALLOCATE,               & ! in
                                 init=.TRUE.                 ) ! in
           ! Allocate and initialise (to 0) radiance_k
           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance_K,                 & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_ALLOCATE,               & ! in
                                 init=.TRUE.                 ) ! in

           IF (ErrorCode /= 0) THEN 
              ErrorMessage(1)='Error in RTTOV_ALLOC_RAD'
              WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                   'Error is Code ',ErrorCode
              ErrStatRep = StatusFatal
              CALL Gen_ErrorReport( RoutineName,  &
                   ErrorMessage, &
                   ErrStatRep    )
           END IF

           ! Initialise radiance_k % bt to 1K as RTTOV11_Opts%rt_all%switchrad
           ! is set to true
           Radiance_k % bt = 1.0

           OCast_Rads(:,:) = 0.

           ! Allocate and initialise Jacobian arrays
           Do i=1,NumInstChans
              Allocate( Profiles_K(i) % p(Num_RTLevels) )
              Profiles_K(i) % p(:) = 0.
              Allocate( Profiles_K(i) % t(Num_RTLevels) )
              Profiles_K(i) % t(:) = 0.
              Allocate( Profiles_K(i) % q(Num_RTLevels) )
              Profiles_K(i) % q(:) = 0.
              Allocate( Profiles_K(i) % o3(Num_RTLevels) )
              Profiles_K(i) % o3(:) = 0.
              Allocate( Profiles_K(i) % clw(Num_RTLevels) )
              Profiles_K(i) % clw(:) = 0.
              Profiles_K(i) % s2m % t = 0.
              Profiles_K(i) % s2m % q = 0.
              Profiles_K(i) % s2m % p = 0.
              Profiles_K(i) % s2m % u = 0.
              Profiles_K(i) % s2m % v = 0.
              Profiles_K(i) % skin % t = 0.
              Profiles_K(i) % ctp = 0.
              Profiles_K(i) % cfraction = 0.
           End Do

           !----- Allocate transmission arrays
           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission, &
                profiles(1) % nlayers, &
                NumInstChans,&
                ASW_ALLOCATE)
           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission_K, &
                profiles(1) % nlayers, &
                NumInstChans,&
                ASW_ALLOCATE)

           transmission_k%tau_levels(:,:)=0.0
           transmission_k%tau_total(:)=0.0
           Allocate(RTTOV11_Chanprof(NumInstChans))
           RTTOV11_Chanprof(1:NumInstChans)%chan = &
                 Inst_Chans(1:NumInstChans)
           RTTOV11_Chanprof(1:NumInstChans)%prof = KPROF*1_JPIM

           ! Initialise emissivity arrays
           IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN

             IF (RTTOV11_Coefs(Instrument)%coef%id_sensor == sensor_id_mw &
                 .OR. RTTOV11_Coefs(Instrument)%coef%id_sensor == &
                 sensor_id_po)  THEN ! MW
               CALL rttov_setup_emis_atlas( &
                    ErrorCode,                      & ! out
                    RTTOV11_Opts(Instrument),       & ! in
                    Profiles(1) % Date(2),          & ! in
                    RTTOV11_Coefs(Instrument),      & ! in
                    path=Atlas_Dir,                 & ! in
                    mw_atlas_ver = Atlas_Ver        ) ! in
             ELSE ! IR
               CALL rttov_setup_emis_atlas( &
                    ErrorCode,                      & ! out
                    RTTOV11_Opts(Instrument),       & ! in
                    Profiles(1) % Date(2),          & ! in
                    RTTOV11_Coefs(Instrument),      & ! in
                    path=Atlas_Dir,                 & ! in
                    ir_atlas_read_std=.TRUE.        ) ! in
             END IF

             IF (ErrorCode >= errorstatus_fatal) THEN 
               ErrorMessage(1)='Error in rttov_atlas_setup'
               WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                 'Error is Code ',ErrorCode
               ErrStatRep = StatusFatal
               CALL Gen_ErrorReport( RoutineName,  &
                    ErrorMessage, &
                    ErrStatRep    )
             END IF                

             Surf_Emiss(1:NumInstChans)%emis_in = 0.
             Surf_Emiss(1:NumInstChans)%emis_out = 0.

             CALL rttov_get_emis( &
                  ErrorCode,                      & ! out
                  RTTOV11_Opts(Instrument),       & ! in
                  RTTOV11_Chanprof,               & ! in
                  Profiles,                       & ! in
                  RTTOV11_Coefs(Instrument),      & ! in
                  Emissivity = Surf_Emiss%emis_in ) ! out

             IF (ErrorCode >= errorstatus_fatal) THEN 
               ErrorMessage(1)='Error in rttov_get_emis'
               WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                 'Error is Code ',ErrorCode
               ErrStatRep = StatusFatal
               CALL Gen_ErrorReport( RoutineName,  &
                 ErrorMessage, &
                 ErrStatRep    )
             END IF

             Calcemiss(1:NumInstChans) = .FALSE.

           ELSE

             Surf_Emiss(1:NumInstChans)%emis_in = 0.
             Surf_Emiss(1:NumInstChans)%emis_out = 0.
             Calcemiss(1:NumInstChans) = .TRUE.

           END IF

           Surf_Emiss_K(1:NumInstChans)%emis_in = 0.
           Surf_Emiss_K(1:NumInstChans)%emis_out = 0.

           CALL RTTOV_K( &
                ErrorStatus(Instrument),              & ! out
                RTTOV11_Chanprof,                     & ! in
                RTTOV11_Opts(Instrument),             & ! in
                Profiles,                             & ! in
                Profiles_K,                           & ! inout
                RTTOV11_Coefs(Instrument),            & ! in
                transmission,                         & ! inout
                transmission_K,                       & ! inout
                Radiance,                             & ! inout
                Radiance_K,                           & ! inout
                calcemiss,                            & ! in
                Surf_Emiss,                           & ! inout
                Surf_Emiss_K )                          ! inout

           ErrorCode = ErrorStatus(Instrument)
           IF (ErrorCode >= errorstatus_fatal) THEN 
              ErrorMessage(1)='Error in RTTOV_K'
              WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
                'Error is Code ',ErrorCode
              ErrStatRep = StatusFatal
              CALL Gen_ErrorReport( RoutineName,  &
                 ErrorMessage, &
                 ErrStatRep    )
           END IF

           Deallocate(RTTOV11_Chanprof)

           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission, &
                profiles(1) % nlayers, &
                NumInstChans,&
                ASW_DEALLOCATE)
           CALL rttov_alloc_transmission( &
                ErrorCode, &
                transmission_k, &
                profiles(1) % nlayers, &
                NumInstChans,&
                ASW_DEALLOCATE)

           IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN
             CALL rttov_deallocate_emis_atlas( RTTOV11_Coefs(Instrument) )
           END IF


           !  Unpack the required output arrays from RTTOV into the 
           !  RT_Param structure
           !------
           RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
                radiance%total(1:NumInstChans)
           RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
                radiance%bt(1:NumInstChans)

           IF (CloudyRetrieval) THEN
              RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
                   1:profiles(1) % nlayers) = &
                        TRANSPOSE(radiance%overcast(:,1:NumInstChans))
! The following commented out as they are not used
!             RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
!                  RTTOV11_Coefs(1)%coef%nlevels+1:&
!                  2*RTTOV11_Coefs(1)%coef%nlevels) = &
!                  TRANSPOSE(radiance%downcld(:,1:NumInstChans))
!             RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
!                  2*RTTOV11_Coefs(1)%coef%nlevels+1) = &
!                  radiance%upclear(1:NumInstChans)
!             RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
!                  2*RTTOV11_Coefs(1)%coef%nlevels+2) = &
!                  radiance%reflclear(1:NumInstChans)
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
                CASE (2)  ! This is humidity (kg/kg)
                   DO i=First_Chan_Pos,Last_Chan_Pos
                     RTModel_Jacobian(Prof_FirstQ:Prof_LastQ, i) = &
                          profiles_k(i-First_Chan_Pos+1)%q(:) !* &
                          !H2O_MassMixToPPMV
		   ENDDO
                CASE (3)  ! This is Ozone
		   IF(RTTOV11_Opts(Instrument) % rt_ir % ozone_data) Then
                     DO i=First_Chan_Pos,Last_Chan_Pos
                       RTModel_Jacobian(Prof_FirstO3:Prof_LastO3, i) = &
		            profiles_k(i-First_Chan_Pos+1)%o3(:)
		     ENDDO
                   END IF
                CASE (4)  
                   ! (Cloud liquid water retrieval is handled as 
                   ! LWP retrieval or qtotal retrieval)
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
                        profiles_k(1:NumInstChans)%s2m%q !* H2O_MassMixToPPMV
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
                   RTModel_Jacobian(Prof_CTP,First_Chan_Pos:Last_Chan_Pos) &
                        = profiles_k(1:NumInstChans)%ctp
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
  
           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance,                   & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_DEALLOCATE   )            ! in
           Call RTTOV_ALLOC_RAD( ErrorCode,                  & ! out
                                 NumInstChans,               & ! in
                                 Radiance_K,                 & ! in/out
                                 profiles(1) % nlayers,      & ! in
                                 ASW_DEALLOCATE   )            ! in

           ! For cloud liquid water profiles, we retrieve liquid water path 
           ! (LWP) or qtotal

           !--------------------------------------------------------------------
           ! 2.5  Calculate total liquid water path from the specific liquid
           !      water profile to create the minimalisation variable LWP_calc
           !      in the state vector. NWPSAF_Layers_to_LWP integrates up the
           !      water on the fixed pressure levels.
           !--------------------------------------------------------------------
        
           IF(MwClwRetrieval .AND. Lqtotal==0) THEN ! LWP retrieval
           
              clw_profiles =  profiles(1)%clw(:)
              CALL NWPSAF_Layers_to_LWP( clw_profiles,      & ! In
                                     RT_params,         &
                                     LWP_calc           ) ! Out

              ! Get the cloud structure
              CALL NWPSAF_CloudStructure(RT_params, cloud_structure(:))
              !-----------------------------------------------------------------
              ! 2.6 Compute dTb/dLWP = dLWP                                    
              !-----------------------------------------------------------------
              DO i=First_Chan_Pos,Last_Chan_Pos
                 CALL NWPSAF_LayerK_to_LWPK (profiles_k(i)%clw(:),         &
                                           cloud_structure,              &
                                           RTModel_Jacobian(Prof_CLW, i) )
              END DO
           END IF
           Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1

           IF(MwClwRetrieval .AND. Lqtotal==1) THEN ! qtotal retrieval
           
              QtotalOption=0
              Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
              wt(1:Num_WetLevels )= RT_Params % RTguess(Prof_LastT - &
                   Num_WetLevels + 1:Prof_LastT)
              wqtotal(1:Num_WetLevels)= EXP(RT_Params % RTguess(Prof_LastQ - &
                   Num_WetLevels + 1:Prof_LastQ)) + RT_Params % RTguess(&
                   Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)
              CALL NWPSAF_Qtot_to_q_ql(wqtotal,     &
                                     wq,          &
                                     wql,         &
                                     wt,          &
                                     RT_Params,   &
                                     QtotalOption )
           
              DO j=1,Num_WetLevels
                 ! Kmat= dTB/dln(qtotal)=qtotal dTB/dqtotal
                 DO i=First_Chan_Pos,Last_Chan_Pos
                    RTModel_Jacobian(Prof_LastQ - Num_WetLevels + 1:Prof_LastQ,&
                         i) = profiles_k(i-First_Chan_Pos+1)%q(Prof_LastT - &
                         Num_WetLevels+1:Prof_LastT) *wq(j) + &
                         profiles_k(i-First_Chan_Pos+1)%clw(Prof_LastT - &
                         Num_WetLevels + 1:Prof_LastT)  * wql(j)           
                 END DO
              END DO
           
              QtotalOption=2
              CALL NWPSAF_Qtot_to_q_ql(wqtotal,     &
                                     wq,          &
                                     wql,         &
                                     wt,          &
                                     RT_Params,   &
                                     QtotalOption )
           
              wpress(:)=RT_Params % Pressure_Pa(Prof_LastT - Num_WetLevels + &
                   1:Prof_LastT)
              delta_eps=1./epsilon -1.
              DO  i=1,Num_WetLevels
                 CALL NWPSAF_svp(wt(i),SESAT) ! SESAT is in Pa.
                 esat(i)     = SESAT
                 CALL NWPSAF_svp_deriv(wt(i),SDlnes_DT)
                 Dlnes_DT(i) = SDlnes_DT
              ENDDO
           
              WHERE (esat(:) > wPress(:))
                 qsat(:)     = 1.0
                 Dqsat_DT(:) = 0.0 
              ELSEWHERE
                 qsat(:)     = epsilon / ( wpress(:)/esat(:) - (1.-epsilon))
                 Dqsat_dT(:) = qsat(:) * Dlnes_DT(:) *(1.+delta_eps*qsat(:))
              ENDWHERE
           
           
              RTModel_Jacobian_Add(:, :) = 0.0
           
              DO j=1,Num_WetLevels
                 DO i=First_Chan_Pos, Last_Chan_Pos
                    RTModel_Jacobian_Add(Prof_LastQ- Num_WetLevels + &
                         1:Prof_LastQ, i)=     &
                         profiles_k(i-First_Chan_Pos+1)%q(Prof_LastT - &
                         Num_WetLevels+1:Prof_LastT) *wq(j) * dqsat_dT(j) + &
                         profiles_k(i-First_Chan_Pos+1)%clw(Prof_LastT - &
                         Num_WetLevels+1:Prof_LastT) *wql(j) * dqsat_dT(j)
                 END DO
              END DO
       
              ! Add this to the temperature jacobian
              DO i=First_Chan_Pos,Last_Chan_Pos
                 RTModel_Jacobian(Prof_LastT - Num_WetLevels + 1:Prof_LastT,&
                      i) = RTModel_Jacobian(Prof_LastT - Num_WetLevels + &
                      1:Prof_LastT, i) + RTModel_Jacobian_Add(Prof_LastQ - &
                      Num_WetLevels + 1:Prof_LastQ, i)
              END DO
           END IF

           ! Deallocate Jacobian arrays
           DO i=1,NumInstChans
              DEALLOCATE( Profiles_K(i) % p )
              DEALLOCATE( Profiles_K(i) % t )
              DEALLOCATE( Profiles_K(i) % q )
              DEALLOCATE( Profiles_K(i) % o3 )
              DEALLOCATE( Profiles_K(i) % clw )
           END DO

           ! End of Jacobian assignments

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
                   profiles(1)%q(Retrieved_Elements(I)-Prof_FirstQ+1) !/ &
                   !H2O_MassMixToPPMV

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

        ! Deallocate profile arrays
        DEALLOCATE( Profiles(1) % p )
        DEALLOCATE( Profiles(1) % t )
        DEALLOCATE( Profiles(1) % q )
        DEALLOCATE( Profiles(1) % o3 )
        DEALLOCATE( Profiles(1) % clw )

   CASE DEFAULT
        WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
             'Incorrect Fastmodel Mode in Fastmodel Interface'
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &
             ErrStatRep    )
   END SELECT



ENDIF RTTOV_FastmodelMode

End Subroutine NWPSAF_RTTOV11_Interface
