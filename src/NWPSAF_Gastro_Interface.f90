Subroutine NWPSAF_Gastro_Interface( &
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
! Description: Interface between NWPSAF 1DVar code and Fastmodels
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     16/10/03 Original version for Gastropod                 V. Sherlock.
! 3.0.3   26/02/04 Incorporated into main NWPSAF_1DVar development.  A. Collard.
! 3.0.4   05/03/04 Replace Plevels_RTmodel_hPa with RT_Params % Pressure_Pa/100
!                  and set up profile limits.                      A. Collard.
! 3.0.6   18/06/04 Include tolerances in soft limits.              A. Collard.
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
     Gastro_ProfMin_Temperature, &
     Gastro_ProfMax_Temperature, &
     Gastro_ProfMin_Humidity, &
     Gastro_ProfMax_Humidity, &
     Gastro_ProfMin_Ozone, &
     Gastro_ProfMax_Ozone

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     StatusWarning,      &
     GeneralMode,        &
     VerboseMode,        &
     Retrieved_Elements, &
     CloudyRetrieval

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
     Soft_Limits, &
     GuessProf, &
     BackGrProf

USE GastroMod_Hardware, ONLY : &
     imepkind

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_MessageReport.interface'
INCLUDE 'Gastro_Initialise.interface'
INCLUDE 'Gastro_1DVar.interface'

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrorCode


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Gastro_Interface"

INTEGER, PARAMETER :: DNAV=4  ! 4 input profile variables (P,T,Q,O3) for Gastro

! Local variables:

INTEGER :: ErrStatRep                 ! error status for Gen_ErrorReport
INTEGER :: NumSatIDs = 1
CHARACTER(LEN=80) :: ErrorMessage(2)  ! Message for Gen_ErrorReport
CHARACTER(LEN=80) :: Message(2)       ! Message for Gen_MessageReport

INTEGER :: I, ILEV
INTEGER :: NLEV
INTEGER :: NLEV_GIP
REAL :: Ps

! Arrays used in interface with both RT Models
!-----------

INTEGER, POINTER :: RTModel_UsedChans (:)
REAL :: RTModel_Profile (ProfSize, Num_Profs)
REAL, ALLOCATABLE :: RTModel_Jacobian (:,:)

REAL, POINTER :: RTModel_TotalBTs(:)                     

INTEGER :: RTModel_RTSurfaceType(Num_Profs)
REAL :: RTModel_SatAzimuthAngle(Num_Profs)        

REAL(kind=imepkind) :: RTModel_SatZenithAngle(Num_Profs)         
REAL(kind=imepkind) :: Gastro_Ts

REAL, ALLOCATABLE :: Gastro_Input_Profile(:,:)
REAL, ALLOCATABLE :: Gastro_Jacobian(:,:,:)
REAL, ALLOCATABLE :: Gastro_TsJacobian(:)

!---------------------------------------------------------------

! Error code is not used at present - set to zero here

IF (CloudyRetrieval) THEN
   WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
        'Gastropod does not support cloudy retrievals'
   ErrStatRep = StatusFatal
   CALL Gen_ErrorReport( RoutineName,  &
        ErrorMessage, &                  
        ErrStatRep    )
END IF

ErrorCode = 0
ErrorMessage(:)=' '

NULLIFY(RTModel_UsedChans)
NULLIFY(RTModel_TotalBTs) 

RTModel_FastmodelMode : IF &
     ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
   IF ( GeneralMode >= VerboseMode ) THEN
      Message(1) = 'Initialising Gastropod'
      CALL Gen_MessageReport( RoutineName,Message(1:1) )
   ENDIF

   CALL Gastro_Initialise( ) 

   ! Set up soft limits
   
   ALLOCATE( Soft_Limits % Minimum(ProfSize))
   ALLOCATE( Soft_Limits % Maximum(ProfSize))

   ! Default values:
   Soft_Limits % Minimum(:) = 0.
   Soft_Limits % Maximum(:) = 1.e10

   ! Profile minima
   Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = Gastro_ProfMin_Temperature(:)-0.5
   Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = Gastro_ProfMin_Humidity(:)-0.5
   Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = Gastro_ProfMin_Ozone(:)*0.8
   Soft_Limits % Minimum(Prof_T2) = Gastro_ProfMin_Temperature(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_Q2) = Gastro_ProfMin_Humidity(Num_RTLevels)-0.5
   Soft_Limits % Minimum(Prof_PStar) = 300.
   Soft_Limits % Minimum(Prof_UWind) = 0.
   Soft_Limits % Minimum(Prof_VWind) = 0.
   Soft_Limits % Minimum(Prof_TStar) = Gastro_ProfMin_Temperature(Num_RTLevels)-0.5

   ! Profile maxima
   Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = Gastro_ProfMax_Temperature(:)+0.5
   Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = Gastro_ProfMax_Humidity(:)+0.5
   Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = Gastro_ProfMax_Ozone(:)*1.2
   Soft_Limits % Maximum(Prof_T2) = Gastro_ProfMax_Temperature(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_Q2) = Gastro_ProfMax_Humidity(Num_RTLevels)+0.5
   Soft_Limits % Maximum(Prof_PStar) = 1200.
   Soft_Limits % Maximum(Prof_UWind) = 100.
   Soft_Limits % Maximum(Prof_VWind) = 100.
   Soft_Limits % Maximum(Prof_TStar) = Gastro_ProfMax_Temperature(Num_RTLevels)+0.5

   DO I=1,ProfSize
      IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN 
         WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
              'Mismatch in Gastropod soft limits'
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

   !---------------------------------------------------------------------
   ! Pack the profile into the form required for Gastropod.

   ! 1. Determine the number of levels.

   Ps=RTModel_Profile(Prof_pstar,RT_Params % SatIndex) 

   NLEV=Num_RTlevels
   DO ilev=1, Num_RTLevels-1
     IF ( ( Ps .gt. RT_Params % Pressure_Pa(ilev) / 100. ) .and. &
        ( Ps .le. RT_Params % Pressure_Pa(ilev+1) / 100. ) ) THEN
	NLEV=ilev
	EXIT
     ENDIF
   ENDDO

   NLEV_GIP=NLEV+1
   ALLOCATE (Gastro_Input_Profile(DNAV, NLEV_GIP)) 
   
   ilev=NLEV_GIP
   
   Gastro_Input_Profile(1,1)=RTModel_Profile(Prof_pstar,RT_Params % SatIndex) 
   Gastro_Input_Profile(2,1)=RTModel_Profile(Prof_T2,RT_Params % SatIndex) 
   Gastro_Input_Profile(3,1)=EXP(RTModel_Profile(Prof_q2,RT_Params % SatIndex)) 
   Gastro_Input_Profile(4,1)=RTModel_Profile(Prof_LastO3,RT_Params % SatIndex) 
   
   DO ilev=NLEV, 1, -1
      Gastro_Input_Profile(1, NLEV_GIP-ilev+1)=&
           RT_Params % Pressure_Pa(ilev) / 100.
      Gastro_Input_Profile(2, NLEV_GIP-ilev+1)= & 
			RTModel_Profile(ilev,RT_Params % SatIndex) 
      Gastro_Input_Profile(3, NLEV_GIP-ilev+1)= & 
		EXP(RTModel_Profile(Num_RTlevels+ilev,RT_Params % SatIndex)) 
      Gastro_Input_Profile(4, NLEV_GIP-ilev+1)= & 
		RTModel_Profile(2*Num_RTlevels+ilev,RT_Params % SatIndex) 
   ENDDO

   Gastro_Ts = RTModel_Profile(Prof_Tstar,RT_Params % SatIndex)

   RTModel_TotalBTs => &
        RT_Params % TotalBTs(1:UsedChans % NumChans) 

   RTModel_UsedChans => UsedChans % Channels(1:UsedChans % NumChans)
   RTModel_SatZenithAngle(:) = RT_Params % SatZenithAngle
   RTModel_SatAzimuthAngle(:) = RT_Params % SatAzimAngle
   RTModel_RTSurfaceType = RT_Params % RTSurfaceType

   SELECT CASE ( Fastmodel_Mode )
     CASE( FastmodelMode_Forward )

        CALL Gastro_1DVar( &
             RTModel_SatZenithAngle(:),   & ! in
	     NLEV_GIP,                 & 
	     Gastro_Input_Profile,        & 
	     Gastro_Ts,                   & 
             UsedChans % NumChans,        & ! in
             RTModel_UsedChans(:),        & ! in
             RTModel_TotalBTs(:) )          ! out

        DEALLOCATE(Gastro_Input_Profile) 

     CASE ( FastmodelMode_Gradient )

        ALLOCATE(Gastro_Jacobian(DNAV,NLEV_GIP, UsedChans % NumChans))
        ALLOCATE(Gastro_TsJacobian(UsedChans % NumChans))

        CALL Gastro_1DVar( &
             RTModel_SatZenithAngle(:),   & ! in 
	     NLEV_GIP,                    & 
	     Gastro_Input_Profile,        & 
	     Gastro_Ts,                   & 
             UsedChans % NumChans,        & ! in
             RTModel_UsedChans(:),        & ! in
             RTModel_TotalBTs(:),         & ! out
             Gastro_Jacobian(:,:,:),        & ! out
	     Gastro_TsJacobian(:))          ! out

        ALLOCATE(RTModel_Jacobian(ProfSize, UsedChans % NumChans))
	RTModel_Jacobian(:,:)=0.0

	! Pack Gastro Jacobians into the RTModel_Jacobian array
        ! Convert Water Vapour Jacobians from kg/kg to log(kg/kg)
   
        RTModel_Jacobian(Prof_T2,:)=Gastro_Jacobian(2,1,:)
        RTModel_Jacobian(Prof_q2,:)=Gastro_Jacobian(3,1,:)*&
				    Gastro_Input_Profile(3,1)
        RTModel_Jacobian(Prof_Tstar,:)=Gastro_TsJacobian(:)
   
        DO ilev=NLEV, 1, -1
           RTModel_Jacobian(ilev,:)=Gastro_Jacobian(2, NLEV_GIP-ilev+1,:) 
           RTModel_Jacobian(Num_RTlevels+ilev,:)=&
	               Gastro_Jacobian(3,NLEV_GIP-ilev+1,:)*&
	               Gastro_Input_Profile(3,NLEV_GIP-ilev+1) 
           RTModel_Jacobian(2*Num_RTlevels+ilev,:)=&
		       Gastro_Jacobian(4, NLEV_GIP-ilev+1,:) 
        ENDDO

        RT_Params % H_matrix_T = RTModel_Jacobian(Retrieved_Elements,:)
        DEALLOCATE(RTModel_Jacobian)

        RT_Params % H_matrix = TRANSPOSE(RT_Params % H_matrix_T)

        DEALLOCATE(Gastro_Input_Profile)
        DEALLOCATE(Gastro_Jacobian)
        DEALLOCATE(Gastro_TsJacobian)

     CASE DEFAULT
        WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
             'Incorrect Fastmodel Mode in Fastmodel Interface'
        ErrStatRep = StatusFatal
        CALL Gen_ErrorReport( RoutineName,  &
             ErrorMessage, &                  
             ErrStatRep    )
   END SELECT
   
   NULLIFY(RTModel_UsedChans)
   NULLIFY(RTModel_TotalBTs) 
END IF RTModel_FastmodelMode

End Subroutine NWPSAF_Gastro_Interface
