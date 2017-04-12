Subroutine NWPSAF_RTIASI_Interface( &
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
! Description: Interface between NWPSAF 1DVar code and RTIASI.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     14/11/01 Original version - split out from Fastmodel_Interface
! 2.2     02/05/02 Added in error message in case the user attempts to 
!                  initialise a subset of channels.
! 3.0.5   15/04/04 Set up profile limits.                      A. Collard.
! 3.0.6   18/06/04 Include tolerances in soft limits.          A. Collard.
!
! Code Description:
!   Language:		Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_Constants, ONLY : &
     RTIASI_ProfMin_Temperature, &
     RTIASI_ProfMax_Temperature, &
     RTIASI_ProfMin_Humidity, &
     RTIASI_ProfMax_Humidity, &
     RTIASI_ProfMin_Ozone, &
     RTIASI_ProfMax_Ozone

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     StatusWarning,      &
     GeneralMode,        &
     VerboseMode,        &
     DebugMode,          &
     Retrieved_Elements, &
     Ret_FirstQ,         &
     Ret_LastQ,          &
     Ret_q2

USE NWPSAFMod_RTmodel, ONLY : &
     Num_RTLevels,             &
     RTParams_Type,            &
     FastmodelMode_CleanUp,    &
     FastmodelMode_Initialise, &
     FastmodelMode_Forward,    &
     FastmodelMode_Gradient,   &
     ProfSize,                 &
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
     Num_Profs, &
     GuessProf, &
     Soft_Limits, &
     BackGrProf

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
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTIASI_Interface"

REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv

! Local variables:

INTEGER :: I, Channel
INTEGER :: ErrStatRep                 ! error status for Gen_ErrorReport
INTEGER :: IPrint = 0
INTEGER :: NumSatIDs = 1
CHARACTER(LEN=80) :: ErrorMessage(2)  ! Message for Gen_ErrorReport
CHARACTER(LEN=80) :: Message(2)       ! Message for Gen_MessageReport

! Arrays used in interface with both RT Models
!-----------

INTEGER, POINTER :: RTModel_UsedChans (:)
INTEGER :: RTModel_RTSurfaceType(Num_Profs)
REAL :: RTModel_Profile (ProfSize, Num_Profs)
REAL, ALLOCATABLE :: RTModel_Jacobian (:,:)
REAL :: RTModel_SatZenithAngle(Num_Profs)         
REAL :: RTModel_SatAzimuthAngle(Num_Profs)        

REAL, POINTER :: RTModel_TotalRadiances(:)
REAL, POINTER :: RTModel_TotalBTs(:)                     

!---------------------------------------------------------------

! Error code is not used at present - set to zero here

ErrorCode = 0
ErrorMessage(:)=' '

NULLIFY(RTModel_UsedChans)
NULLIFY(RTModel_TotalRadiances) 
NULLIFY(RTModel_TotalBTs) 

RTModel_FastmodelMode : IF &
     ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
   IF ( GeneralMode >= VerboseMode ) THEN
      Message(1) = 'Initiailising RTIASI'
      CALL Gen_MessageReport( RoutineName,Message(1:1) )
   ENDIF

   IF ( ANY(RT_Params % Absolute_Channel_Number(:) /= 0)) THEN
      ErrorMessage(1) = 'Cannot initialise a subset of channels with'
      ErrorMessage(2) = 'this version of RTIASI'
      ErrStatRep = StatusWarning
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
   ENDIF

   CALL RTIASI_RTTVI( &
        Num_Profs, & ! in
        ErrorCode)   ! out

   ! Set up soft limits
   
   ALLOCATE( Soft_Limits % Minimum(ProfSize))
   ALLOCATE( Soft_Limits % Maximum(ProfSize))

   ! Default values:
   Soft_Limits % Minimum(:) = 0.
   Soft_Limits % Maximum(:) = 1.e10

   ! Profile minima
   Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = RTIASI_ProfMin_Temperature(:)-0.5
   Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = RTIASI_ProfMin_Humidity(:)-0.5
   Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = RTIASI_ProfMin_Ozone(:)*0.8
   Soft_Limits % Minimum(Prof_T2) = RTIASI_ProfMin_Temperature(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_Q2) = RTIASI_ProfMin_Humidity(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_PStar) = 300.
   Soft_Limits % Minimum(Prof_UWind) = 0.
   Soft_Limits % Minimum(Prof_VWind) = 0.
   Soft_Limits % Minimum(Prof_TStar) = RTIASI_ProfMin_Temperature(Num_RTLevels)-0.5

   ! Profile maxima
   Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = RTIASI_ProfMax_Temperature(:)+0.5
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = RTIASI_ProfMax_Humidity(:)+0.5
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = RTIASI_ProfMax_Ozone(:)*1.2
   Soft_Limits % Maximum(Prof_T2) = RTIASI_ProfMax_Temperature(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_Q2) = RTIASI_ProfMax_Humidity(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_PStar) = 1200.
   Soft_Limits % Maximum(Prof_UWind) = 100.
   Soft_Limits % Maximum(Prof_VWind) = 100.
   Soft_Limits % Maximum(Prof_TStar) = RTIASI_ProfMax_Temperature(Num_RTLevels)+0.5

   DO I=1,ProfSize
      IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Mismatch in RTIASI soft limits'
         WRITE( UNIT=ErrorMessage(2),FMT=* )  I, &
              Soft_Limits % Minimum(I), &
              Soft_Limits % Maximum(I)
         ErrStatRep = StatusWarning
         CALL Gen_ErrorReport( RoutineName,  &
              ErrorMessage, &
              ErrStatRep    )
      END IF
   END DO

ELSE IF (FastModel_Mode /= FastModelMode_CleanUp) THEN

   NumSatIDs = 1
   !---------------------------------------------------------------------
   !  Put Profile into the RTIASI vector remembering to change the 
   !  humidity variable accordingly.
   !---------------------------------------------------------------------

   IF (WhichProf == BackGrProf) THEN
      RTModel_Profile(:, RT_Params % SatIndex) = RT_Params % RTBack
   ELSE IF (WhichProf == GuessProf) THEN
      RTModel_Profile(:, RT_Params % SatIndex) = RT_Params % RTGuess
   ELSE
      WRITE( UNIT=ErrorMessage(1),FMT='(A,I2,A)' )  &
           'Incorrect profile (',WhichProf,') specified'
      ErrStatRep = StatusWarning
      CALL Gen_ErrorReport( RoutineName,  &
           ErrorMessage, &
           ErrStatRep    )
    END IF
   
   RTModel_Profile(Prof_Firstq:Prof_Lastq,RT_Params % SatIndex) =  &
        EXP(RTModel_Profile(Prof_Firstq:Prof_LastQ, RT_Params % SatIndex)) * &
        H2O_MassMixToPPMV
   RTModel_Profile(Prof_q2, RT_Params % SatIndex) = &
        EXP(RTModel_Profile(Prof_q2, RT_Params % SatIndex)) * &
        H2O_MassMixToPPMV
   
   RTModel_TotalRadiances => &
        RT_Params % TotalRadiances(1:UsedChans % NumChans)
   RTModel_TotalBTs => &
        RT_Params % TotalBTs(1:UsedChans % NumChans) 

   RTModel_UsedChans => UsedChans % Channels(1:UsedChans % NumChans)
   RTModel_SatZenithAngle(:) = RT_Params % SatZenithAngle
   RTModel_SatAzimuthAngle(:) = RT_Params % SatAzimAngle
   RTModel_RTSurfaceType = RT_Params % RTSurfaceType
   IF (GeneralMode >= DebugMode) IPrint = 1
   SELECT CASE ( Fastmodel_Mode )
     CASE( FastmodelMode_Forward )

        CALL RTIASI_Direct( &
             RTModel_Profile(:,:),        & ! in
             RTModel_UsedChans(:),        & ! in
             RTModel_SatZenithAngle(:),   & ! in
             RTModel_SatAzimuthAngle(:),  & ! in
             RTModel_RTSurfaceType(:),    & ! in
             Num_Profs,                   & ! in (This had better be 1!)
             UsedChans % NumChans,        & ! in
             ProfSize,                    & ! in 
             IPrint,                      & ! in
             RTModel_TotalRadiances(:),   & ! out
             RTModel_TotalBTs(:) )          ! out

     CASE ( FastmodelMode_Gradient )

        ALLOCATE(RTModel_Jacobian(ProfSize, UsedChans % NumChans * Num_Profs))
        CALL RTIASI_K( & 
             RTModel_Profile(:,:),       & ! in
             RTModel_UsedChans(:),       & ! in
             RTModel_SatZenithAngle(:),  & ! in
             RTModel_SatAzimuthAngle(:), & ! in
             RTModel_RTSurfaceType(:),   & ! in
             Num_Profs,                  & ! in (This had better be 1!)
             UsedChans % NumChans,       & ! in
             ProfSize,                   & ! in
             RT_Params % SatIndex,       & ! in
             IPrint,                     & ! in
             RTModel_TotalRadiances(:),  & ! out
             RTModel_TotalBTs(:),        & ! out
             RTModel_Jacobian(:,:) )       ! out

        RT_Params % H_matrix_T = RTModel_Jacobian(Retrieved_Elements,:)
        DEALLOCATE(RTModel_Jacobian)

        !-----------------------------------------------------------------
        ! Water Vapour Jacobians must be converted from ppmv to log(kg/kg)
        !-----------------------------------------------------------------
        IF (Ret_FirstQ > 0) THEN
           DO I = Ret_FirstQ, Ret_LastQ
              IF (Retrieved_Elements(I) > 0) THEN
                 RT_Params % H_matrix_T(I,:) = &
                      RT_Params % H_matrix_T(I,:) * &
                      RTModel_Profile(Retrieved_Elements(I),1)
                 
                 !-----------------------------------------------------------
                 ! The following picks up those pesky RTIASI Jacobian spikes 
                 ! (allowing them to be ingored for the iteration). 
                 !-----------------------------------------------------------
                 DO Channel=1, UsedChans % NumChans
                    IF (ABS(RT_Params % H_matrix_T(I,Channel)) > 5.) THEN
                       WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
                            'Bad Jacobian for Channel/Prof Element ='
                       WRITE( UNIT=ErrorMessage(2),FMT='(3I4)' )  &
                            UsedChans % Channels (Channel),Channel,I
                       !xxx                         ErrorCode = 9000
                       ErrStatRep = StatusWarning
                       CALL Gen_ErrorReport( RoutineName,  &
                            ErrorMessage, &
                            ErrStatRep    )
                       ! Zero Jacobian for entire channel
                       RT_Params % H_matrix_T(:,Channel) = 0.
                    END IF
                 END DO
              END IF
           END DO
        END IF
        
        !-----------------------------------------------------------------
        ! This is the surface humidity Jacobian
        !-----------------------------------------------------------------
        IF (Ret_q2 > 0) THEN
           IF (Retrieved_Elements(Ret_q2) > 0) &
                RT_Params % H_matrix_T(Ret_q2,:) = &
                RT_Params % H_matrix_T(Ret_q2,:) * &
                RTModel_Profile(Retrieved_Elements(Ret_q2),1) 
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
   
   NULLIFY(RTModel_UsedChans)
   NULLIFY(RTModel_TotalRadiances) 
   NULLIFY(RTModel_TotalBTs) 
END IF RTModel_FastmodelMode

End Subroutine NWPSAF_RTIASI_Interface
