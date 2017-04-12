Subroutine NWPSAF_RTTOV11_Brute_Force_Jacobians ( &
     NumInstChans ,  & ! number of channels for the current instrument
     First_chan_pos  , & ! premier canal pour l'instrument 
     Last_chan_pos,  & ! dernier canal pour l'instrument 
     UsedChans,        & ! in %channels for the computation of brute force
     Profiles_K  ,         &
     RT_Params,       & ! inout
     ErrorCode)                 ! out

! This routine is intended to compute brute force Jacobians instead of
!calling RTTOV_K while this routine is not ready.
! 27/10/2014  P. MARTINET

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
     BruteForceProf, &
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

USE rttov_types, ONLY: &
     profile_type,      &
     radiance_type,     &
     transmission_type, &
     rttov_emissivity

Use parkind1, Only : jpim,jprb

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'Gen_MessageReport.interface'
INCLUDE 'rttov_alloc_rad.interface'
INCLUDE 'rttov_alloc_transmission.interface'
#include <rttov_dealloc_coefs.interface>
INCLUDE 'rttov_direct.interface'
INCLUDE 'rttov_k.interface'
#include <rttov_read_coefs.interface>
#include <rttov_init_coefs.interface>
INCLUDE 'NWPSAF_Layers_to_LWP.interface'

! Subroutine arguments:

TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans ! channels to compute the brute force jacobians
TYPE(profile_type),INTENT(INOUT)  :: Profiles_K(UsedChans % NumChans * Num_Profs)
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(OUT)                    :: ErrorCode
INTEGER, INTENT(IN) :: First_Chan_Pos !) Position of first and last channels for  
INTEGER, INTENT(IN) :: Last_Chan_Pos  !) current instrument in RT_Params % TotalBTs
INTEGER, INTENT(IN) :: NumInstChans

!Local variables
INTEGER :: Fastmodel_Mode       !   Mode in which fastmodel is called
REAL :: BT_original_profiles(NumInstChans)       ! BT differences
REAL :: BT_perturbed_profiles(Num_Atm_Prof_Vars*Num_RTlevels,NumInstChans)       ! BT differences
INTEGER :: Ichan
INTEGER :: Ilevel,Variable_Number ,j
REAL :: delta_profile
REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv
INTEGER::WhichProf   ! Which profile to use 
REAL, POINTER :: TEST(:)
REAL :: profile_inc(Num_RTlevels*Num_Atm_Prof_Vars)

TEST => RT_Params % RTGuess
RT_Params % RTGuess_BF(Prof_FirstT:Prof_CloudCover) = TEST(Prof_FirstT:Prof_CloudCover) !INITIALISATION OF THE BRUTE FORCE CALCULATION WITH
WhichProf = BruteForceProf !guess profile is the one modified in the 1D-Var
Fastmodel_Mode = FastmodelMode_Forward

!=====LOOP OVER CHANNELS AND PROFILES LEVELS TO PERTURB THE PROFILE
DO Variable_Number = 1, Num_Atm_Prof_Vars
	SELECT CASE(Variable_Number)
      		CASE (1)  ! This is temperature
			DO Ilevel =1,Num_RTlevels
  !       		 	RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1) = RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1) + 1 ! 1K perturbation
        		 	RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1) = RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1) &
& - 0.005*RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1) ! 1K perturbation
				profile_inc(Prof_FirstT+Ilevel-1)= - 0.005*RT_Params % RTGuess_BF(Prof_FirstT+Ilevel-1)
				
				!write(*,*) 'computation jacobian Temperature', Ilevel
				CALL NWPSAF_Fastmodel_interface(     &
     				Fastmodel_Mode,          & ! in
     				RT_Params,               & ! in
     				WhichProf,               & ! in
     				UsedChans,               & ! in
     				ErrorCode)                 ! out			

				RT_Params % RTGuess_BF(Prof_FirstT:Prof_CloudCover) = TEST(Prof_FirstT:Prof_CloudCover) ! reinitialization of the profile 
				BT_perturbed_profiles(Prof_FirstT+Ilevel-1,:)=RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos)

			ENDDO
		 CASE (2)  ! This is humidity (in ppmv) (convert from LOG kg/kg)
			DO Ilevel =1,Num_RTlevels
        		 	RT_Params % RTGuess_BF(Prof_FirstQ+Ilevel-1) = RT_Params % RTGuess_BF(Prof_FirstQ+Ilevel-1) &
& - 0.01*RT_Params % RTGuess_BF(Prof_FirstQ+Ilevel-1) ! 1K perturbation
				profile_inc(Prof_FirstQ+Ilevel-1)= - 0.01*RT_Params % RTGuess_BF(Prof_FirstQ+Ilevel-1)

				CALL NWPSAF_Fastmodel_interface(     &
     				Fastmodel_Mode,          & ! in
     				RT_Params,               & ! in
     				WhichProf,               & ! in
     				UsedChans,               & ! in
     				ErrorCode)                 ! out

				RT_Params % RTGuess_BF(Prof_FirstT:Prof_CloudCover) = TEST(Prof_FirstT:Prof_CloudCover) ! reinitialization of the profile 

				BT_perturbed_profiles(Prof_FirstQ+Ilevel-1,:)=RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos)
				!DO Ichan=1,UsedChans % NumChans
				!	Jacobians(Prof_FirsTQ+Ilevel-1,Ichan)=(BT_perturbed_profiles(Ichan)-BT_original_profiles(Ichan))/EXP(0.001)
				!ENDDO
			ENDDO
	END SELECT
END DO
write(*,*) 'END OF JACOBIANS CALCULATION=============='

!=====CALL TO THE DIRECT MODEL TO COMPUTE DIRECT BTs========
WhichProf = GuessProf !guess profile is the one modified in the 1D-Var
Fastmodel_Mode = FastmodelMode_Forward
CALL NWPSAF_Fastmodel_interface(     &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! in
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrorCode)                 ! out

BT_original_profiles=RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos)

!===============================================

!===COMPUTATION OF THE JACOBIANS========================
DO Ilevel =1,Num_RTlevels
	DO Ichan=1, NumInstChans
		Profiles_K(Ichan)%t(Ilevel)=(BT_perturbed_profiles(Prof_FirstT+Ilevel-1,Ichan)-BT_original_profiles(Ichan))/&
&profile_inc(Prof_FirstT+Ilevel-1)
		Profiles_K(Ichan)%q(Ilevel)=(BT_perturbed_profiles(Prof_FirsTQ+Ilevel-1,Ichan)-BT_original_profiles(Ichan))/&
&profile_inc(Prof_FirstQ+Ilevel-1)
	ENDDO
ENDDO

IF(NumInstChans > 5) THEN

	OPEN(UNIT=3000,FILE=&
	&'/home/francesco/Documenti/1DVAR/1DVar_v1.0_bf/&
	&Jacobians.dat',FORM='FORMATTED',ACCESS='SEQUENTIAL')

	DO Ichan=1,NumInstChans
		DO Ilevel =1,Num_RTlevels
			write(3000,*) Profiles_K(Ichan)%t(Ilevel)
			!write(3000,*) Profiles_K(Ichan)%q(Ilevel)
		ENDDO
	ENDDO

 	CLOSE(3000)

	OPEN(UNIT=3000,FILE=&
	&'/home/francesco/Documenti/1DVAR/1DVar_v1.0_bf/&
	&Pressure.dat',FORM='FORMATTED',ACCESS='SEQUENTIAL')

	DO Ilevel =1,Num_RTlevels
		write(3000,*) RT_Params%Pressure_Pa(Ilevel)/100.0
	ENDDO

	CLOSE(3000)


ENDIF

End Subroutine NWPSAF_RTTOV11_Brute_Force_Jacobians





