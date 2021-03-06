!+ Find a solution to the satellite sounding inverse problem

Subroutine NWPSAF_Minimize( &
     Obs,             & ! inout
     BTObserved,      & ! in
     Bmatrix,         & ! in
     BmatrixInverse,  & ! in
     R_matrix,        & ! in
     obnumber,        & ! in
     UsedChans,       & ! in
     RT_Params,       & ! inout
     ErrorCode, &
     converged       ) ! out
  
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
! Description:
!
!   find the most probable atmospheric state vector by minimising a cost
!   function through a series of iterations. If a solution exists, the
!   iterations will converge when the iterative increments are acceptably
!   small. A limit on the total number of iterations allowed, is imposed.
!
!   Two formulations are used.
!
!   The first uses Eqn. 101 in Rodgers (1976), this is ONLY used where:
!     1) The length of the observation vector is less than the length 
!        of the state vector, 
!     2) Where no additional cost function terms are provided and 
!     3) where Newtonian minimisation is desired.
!   The form of the solution is then:
!
!   x_(n+1) = xb + Wn.(ym-y(xn) - H'.(xb-xn))
!   where: x is an atmospheric state vector, subscripted b=background,n=nth
!           iteration
!          Wn = B.Hn'.(Hn.B.Hn'+R)^-1
!          B is the background error covariance matrix
!          R is the combined forward model and ob error covariance matrix
!
!   Delta_x = (x_(n+1) - x_n) is checked for convergence after each iteration
!
!   The second formulation is used in all other cases.  The solution is 
!
!   x_(n+1) = xb + U^-1 V
!   where U=(B^-1 + H^T R^-1 H + J2 + gamma I)
!         V=H^T R^-1 [(y-y(x_n))+H(x_n-xb)] + gamma (x_n-xb) - J1
!   and   J_extra=J0+J1.(x-xb)+(x-xb)^T.J2.(x-xb) is the additional cost 
!         function
!         gamma = is the multiplier used in the Marquardt-Levenberg routine
!                                    (set to zero for simple inverse Hessian).
!
!   When gamma and J_extra are zero this is simply Rogers (1976), Eqn. 100.
!
! Reference: 
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     23/11/99 Original code based on ATOVS_Minimize.  A. Collard.
! 1.2     07/01/02 Fix bug so that now the cost function is evaluated
!                  with the DeltaBT for the current rather than last 
!                  iteration.                  A. Collard
! 2.2     19/02/02 Make sure Guess is correctly updated if the profile
!                  variables are reset due to reaching limits.  
!                  Also clean up the convergence tests.    A. Collard.
! 2.3     23/02/02 Modified calls to NWPSAF_Calculate_Cost_Function.
!                  Minor modifications to ensure that the first guess
!                  does not *have* to equal the background.
!                  Change to Obs % JCost_Gradient calculation.
!                                                          A. Collard.
! 3.0.1   15/07/03 Add cloud cover code (from M. Szyndel). A. Collard.
! 3.0.2   04/02/04 Pass RT_Params%Guess to NWPSAF_Calculate_Cost_Function
!                  to be used in additional cost function calc. ADC.
! 3.0.4   04/03/04 Rename NWPSAF_MaxIterations, MaxIterations.  Pass presures
!                  to NWPSAF_CheckIteration.                 A. Collard.
! 3.0.5   02/04/04 Remove superfluous code rescaling JCost_Gradient and 
!                  improve gradient test.
!                                                          A. Collard.
! 3.0.6   21/06/04 Remove pressure argument to NWPSAF_CheckIteration.
!                                                          A. Collard.
!
! Ticket Data      Comment
! ------ ----       -------
! 28     22/02/12  Added option for cloud liquid water retrieval TR Sreerekha
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type, &
    R_Eigenvectors, &
    Analysis_Error_Covariance, &
    Prop_Measurement_Noise

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type

USE NWPSAFMod_Params, ONLY :    &
    GeneralMode,              &
    DebugMode,                &
    VerboseMode,              &
!!!fradea
    Read_CLW_Background,      &
!!!fradea 
    MaxIterations,            &
    Additional_Cost_Function, &
    No_Additional_Cost_Function, &
    Minimisation_Method,      &
    Newtonian,                &    
    Marquardt_Levenberg,      &
    DeltaJ,                   &
    SmallJCost_Gradient,      &
    Allow_Eqn_101,            &
    Force_Eqn_101,            &
    Retrieved_Elements,       &
    CloudyRetrieval,          &
    MwClwRetrieval,           &
    IWVRetrieval,           & !PM change
    Lqtotal,                  &
    FileUnit_MinimisationLog,   &
    FileUnit_BTMinimisationLog, &
    FileUnit_AMatrix, FileUnit_AmMatrix, &
    Ret_FirstQ, &
    Ret_LastQ

USE NWPSAFMod_RTmodel, ONLY : &
     Num_RTlevels, &  
     RTParams_Type, &
     Num_ProfElementsUsed, &
     Prof_CloudCover, &
     Prof_CTP, &
     Prof_CLW, &
     Prof_IWV, & !PM change
     FastmodelMode_Gradient, &
     GuessProf, &
     Prof_FirstCLW, &
     Prof_LastCLW, &
     Prof_FirstT, &
     Prof_LastT, &
     Prof_FirstQ, &
     Prof_LastQ

 Use parkind1, Only : jpim,jprb

IMPLICIT NONE

INCLUDE 'NWPSAF_Calculate_Cost_Function.interface'
INCLUDE 'NWPSAF_CheckIteration.interface'
INCLUDE 'NWPSAF_Fastmodel_Interface.interface'
INCLUDE 'NWPSAF_Minimize_100.interface'
INCLUDE 'NWPSAF_Minimize_100ML.interface'
INCLUDE 'NWPSAF_Minimize_101.interface'
INCLUDE 'NWPSAF_RMatrix_ChanSelect.interface'
INCLUDE 'NWPSAF_SatMatInv.interface'
INCLUDE 'NWPSAF_Layers_to_LWP.interface'

! Subroutine arguments:
REAL,    INTENT(IN)  :: BTObserved(:)        ! Observation BT's
REAL,    INTENT(IN)  :: Bmatrix(:,:)         ! B-matrix
REAL,    INTENT(IN)  :: BmatrixInverse(:,:)  ! inverse of B-matrix
INTEGER, INTENT(IN)  :: obnumber             ! Profile number
INTEGER, INTENT(OUT) :: ErrorCode
TYPE(Ob_type)           , INTENT(INOUT) :: Obs       ! Observed/Retrieval data 
TYPE(R_Matrix_Type)        , INTENT(IN) :: R_Matrix  ! R-Matrix
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans ! Chans to be used
TYPE(RTParams_Type)     , INTENT(INOUT) :: RT_Params ! RT Model Data
LOGICAL, INTENT(OUT) :: Converged      
! Local variables:
INTEGER :: Fastmodel_Mode       !   Mode in which fastmodel is called
INTEGER :: Work_Dimension   ! Work array dimension
INTEGER :: Iterate  ! loop counter
INTEGER :: Status   ! from NWPSAF_SatMatInv
REAL :: Back(Num_ProfElementsUsed)          ! background profile
REAL :: Guess(Num_ProfElementsUsed)         ! guess profile 
REAL :: Guess_CTP                           ! guess CTP
REAL :: Guess_CldFrac                       ! guess Cloud Fraction
REAL :: LWP_back
REAL :: IWV_back !PM change
REAL :: DeltaBT(UsedChans % NumChans)       ! BT differences
REAL :: Delta_Profile(Num_ProfElementsUsed) ! profile increments 
REAL :: Gamma                      !Used in Marquart-Levenberg minimisation
REAL :: Old_Gamma                  !Used in Marquart-Levenberg minimisation
REAL :: JCost        ! Cost function
REAL :: JCost_Gradient(Num_ProfElementsUsed) ! Cost function gradient
REAL :: JOld         ! Stored Cost function
REAL :: LWP_FirstGuess
REAL :: IWV_FirstGuess !PM change
REAL :: cloud_structure(Num_RTlevels)
REAL :: q_structure(Num_RTlevels) !PM change
!REAL (kind=jprb)           :: clw(Num_RTlevels)

! splitting of total water into water vapor and cloud water
REAL            :: LWP_split
REAL            :: LWP_splitmod
REAL            :: LWP_ratio
INTEGER         :: QtotalOption
INTEGER         :: Num_WetLevels
REAL            :: wt(Ret_LastQ - Ret_FirstQ + 1) ! Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
REAL            :: wqtotal(Ret_LastQ - Ret_FirstQ + 1)
REAL, PARAMETER ::  H2O_MassMixToPPMV = 1.6078e6    ! kg/kg -> ppmv
REAL            :: wq(Ret_LastQ - Ret_FirstQ + 1)
REAL            :: wql(Ret_LastQ - Ret_FirstQ + 1)


TYPE(R_Matrix_Type) :: R_SubMatrix  ! Channel-selected R-matrix

LOGICAL :: Eqn_101 = .FALSE.   
LOGICAL :: Out_of_range      
LOGICAL :: Profile_Variables_Reset

INTEGER :: WhichProf
INTEGER :: i

!-----------------------------------------------------------------------------

!-----------------------
!1. Initialize variables
!-----------------------

ErrorCode = 0

Out_of_range = .FALSE.
Converged    = .FALSE.

! N.B. Some of the atmospheric profile elements used in the RT calculations are
! either not used in minimisation or have different units. 

Back = RT_Params % RTBack(Retrieved_Elements)


IF (CloudyRetrieval) THEN

   !Set first guess profile = background profile 
   !leaving CTP and cloud fraction guess unchanged

   Guess_CTP = RT_Params % RTguess(Prof_CTP)
   Guess_CldFrac = RT_Params % RTguess(Prof_CloudCover)

   RT_Params % RTguess(Prof_CTP) = Guess_CTP
   RT_Params % RTguess(Prof_CloudCover) = Guess_CldFrac

END IF

IF (MwClwRetrieval.AND. Lqtotal==1 ) THEN
   LWP_split = 0.0
   LWP_splitmod =0.0
   QtotalOption=1
   Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
   wt(1:Num_WetLevels )= RT_Params % RTguess(Prof_LastT - Num_WetLevels + 1:Prof_LastT)
   !RT_Params % RTguess(humidity) is in Log(kg/kg), hence take exponential to get units in kg/kg
   ! Multiply humidity units in kg/kg by H2O_MassMixToPPMV to get humidity units in ppmv 
   !(not necessary here but when input to rttov this has to be done)
   ! clw is in units of kg/kg anyway
   ! add both in the same units to get wqtotal

!   wqtotal(1:Num_WetLevels)=EXP(RT_Params % RTguess(Prof_LastQ- Num_WetLevels + 1:Prof_LastQ)) + &
 !       RT_Params % RTguess(Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)
   wqtotal(1:Num_WetLevels)=RT_Params % RTguess(Prof_LastQ- Num_WetLevels + 1:Prof_LastQ) + &
        RT_Params % RTguess(Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)
    CALL NWPSAF_Qtot_to_q_ql(wqtotal, wq, wql, wt,RT_Params, QtotalOption)
   RT_Params % RTguess(Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)= wql
   !clw(:)=RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW)
   CALL NWPSAF_Layers_to_LWP(RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW),RT_Params,LWP_split)
   IF (LWP_split > 0.5) THEN
      LWP_ratio = 0.5/LWP_split
      RT_Params % RTguess(Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)= &
           RT_Params % RTguess(Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)*LWP_ratio             
      !clw(:) = RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW)
      CALL NWPSAF_Layers_to_LWP(RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW),RT_Params,LWP_splitmod)
   ENDIF
 END IF

IF (MwClwRetrieval.AND. Lqtotal==0 ) THEN
      CALL NWPSAF_Layers_to_LWP(RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW),RT_Params, LWP_back)
      RT_Params % RTBack(Prof_ClW) = LWP_back
END IF

!===========MODIFICATION PAULINE MARTINET================================
IF (IWVRetrieval ) THEN
      CALL NWPSAF_Layers_to_LWP(RT_Params % RTguess(Prof_FirstQ:Prof_LastQ),RT_Params, IWV_back)
      RT_Params % RTBack(Prof_IWV) = IWV_back
      IWV_FirstGuess = IWV_back! RT_Params % RTguess(Prof_ClW) !=0.1 ! kg/m**2
      RT_Params % RTguess(Prof_IWV) = IWV_FirstGuess
     RT_Params % RT1stguess(Prof_IWV) = IWV_FirstGuess 
      !q_structure(:)=0
      !q_structure(:) = &
        !RT_Params % RTBack(Prof_FirstQ:Prof_LastQ)/RT_Params % RTBack(Prof_IWV)
        !RT_Params % RTguess(Prof_FirstQ:Prof_LastQ) = IWV_FirstGuess * q_structure(:)
END IF
!================================================================


IF(MwClwRetrieval .AND. Lqtotal==0) THEN
   !!!fradea
   IF (Read_CLW_Background == 1) THEN 
       LWP_FirstGuess = LWP_back! RT_Params % RTguess(Prof_ClW) !=0.1 ! kg/m**2i
   ELSE
       LWP_FirstGuess = 0.1! RT_Params % RTguess(Prof_ClW) !=0.1 ! kg/m**2
   ENDIF
   !!!fradea
   RT_Params % RTguess(Prof_ClW) = LWP_FirstGuess
   RT_Params % RT1stguess(Prof_ClW) = LWP_FirstGuess
   cloud_structure(:) = 0.0
   CALL NWPSAF_LWP_to_Layers(LWP_firstguess, RT_Params,cloud_structure )
   !The cloud liquid water profile is modified when using this LWP_firstGuess
END IF

Guess = RT_Params % RTGuess(Retrieved_Elements)

Delta_Profile(:) = Guess - Back 

! Check that the first guess is within bounds acceptable to the RT model.
! If not, RTGuess is updated but NOT RTBack (i.e., only the linearisation  
! point and not the a priori information is changed).


!!!fradea
CALL NWPSAF_CheckIteration( &
     RT_Params % RTGuess,      & ! inout
     Delta_Profile,            & ! inout
     Profile_Variables_Reset,  & ! out
     Out_of_Range,             & ! out
     Bmatrix              )      ! in
!!!fradea


! If profile variables are reset, we need to reset Guess again
IF (Profile_Variables_Reset) Guess(:) = RT_Params % RTGuess(Retrieved_Elements)

ALLOCATE(Analysis_Error_Covariance(Num_ProfElementsUsed,Num_ProfElementsUsed))
ALLOCATE(Prop_Measurement_Noise(Num_ProfElementsUsed,Num_ProfElementsUsed))

!------------------------------------------------------------------------
! 1.1. R_SubMatrix is allocated and set up.  It is the observational
!      + forward model error covariance matrix after channel selection.
!------------------------------------------------------------------------

CALL NWPSAF_RMatrix_ChanSelect( &
     UsedChans, &   ! in
     R_Matrix,  &   ! in
     R_SubMatrix)   ! out

!-----------------
!2. Iterative Loop
!-----------------
! The loop is exited on one of three conditions:
! 1) The improvement in the cost function is sufficiently
!    small to imply convergence at an acceptable solution.
! 2) The maximum number of allowed iterations has been reached. In most
!    cases, one of the above criteria will have occurred.
! 3) An error occurs (or the profile goes out of bounds)

! First get H matrix for background situation.

WhichProf = GuessProf
Fastmodel_Mode = FastmodelMode_Gradient
ALLOCATE(RT_Params % H_matrix(UsedChans % NumChans, Num_ProfElementsUsed))
ALLOCATE(RT_Params % H_matrix_T(Num_ProfElementsUsed,UsedChans % NumChans))

CALL NWPSAF_Fastmodel_interface(     &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! in
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrorCode)                 ! out

DeltaBT = BTObserved(UsedChans % Channels(1 : UsedChans % NumChans)) - &
     RT_Params % TotalBTs(1 : UsedChans % NumChans)

IF (GeneralMode >= DebugMode) THEN
   WRITE(FileUnit_BTMinimisationLog,*) '------------------------------'
   WRITE(FileUnit_BTMinimisationLog,*) 'ObNumber =',obnumber
   WRITE(FileUnit_BTMinimisationLog,*) 'Iteration =',0  
   WRITE(FileUnit_BTMinimisationLog,*) 'Gamma =',0
   WRITE(FileUnit_BTMinimisationLog,*) 'Brightness Temperature Difference:'
   WRITE(FileUnit_BTMinimisationLog,*) DeltaBT
   WRITE(FileUnit_BTMinimisationLog,*) 
END IF

! Set a flag if the Rogers Eqn. 101 form is to be used

IF (UsedChans % NumChans < Num_ProfElementsUsed .AND. &
     Additional_Cost_Function == No_Additional_Cost_Function .AND. &
     Minimisation_Method == Newtonian .AND. Allow_Eqn_101 ) THEN
   Eqn_101 = .TRUE.
ELSEIF (Force_Eqn_101) THEN
   Eqn_101 = .TRUE.
ELSE
   Eqn_101 = .FALSE.
END IF


! Set up the dimension of the extra work array used in the minimisation
! routines.

IF (R_SubMatrix % RType == R_Eigenvectors) THEN
   Work_Dimension = R_SubMatrix % Num_Elements
ELSE
   Work_Dimension = UsedChans % NumChans
ENDIF

! Calculate cost function so that convergence can be tested after first 
! iteration (also used in Marquardt-Levenberg algorithm)

Delta_Profile = Guess - Back

CALL NWPSAF_Calculate_Cost_Function( &
     BMatrixInverse,         & ! in
     R_SubMatrix,            & ! inout
     Delta_Profile,          & ! in
     DeltaBT,                & ! in
     Work_Dimension,         & ! in
     Eqn_101,                & ! in
     RT_Params % RTGuess,    & ! in
     JCost,                  & ! out
     Status)                   ! out
JOld = JCost

! If Marquardt-Levenberg is being used, calculate the cost-function
! for the first guess and set up Gamma and some more work arrays.

IF (Minimisation_Method == Marquardt_Levenberg) THEN
   Gamma = 1.e-4
ELSE
   Gamma = 0.
   Old_Gamma = 1.
END IF

Iteration_Loop: DO Iterate = 1, MaxIterations

   IF ( GeneralMode >= VerboseMode ) WRITE(*,*) 'ITER=',Iterate

   !2.1) Calculate the new profile
   !----
   
   Delta_Profile = Guess - Back

   ! The observation+forward model error is specified one of three ways
   ! (full matrix, band diagonal and eigenvalues/eigenvectors).  The 
   ! following subroutines determines in which form this error covariance is 
   ! given, calculates H^T.R^-1.H and H^T.O^-1.(y-y(x)) and from there the U 
   ! matrix and V vectors.  Finally a new profile increment is determined by 
   ! solving the appropriate retrieval equation using Cholesky decomposition.

   IF (Eqn_101) THEN
      CALL NWPSAF_Minimize_101(      &
           Bmatrix,                & ! in
           DeltaBT,                & ! in  
           UsedChans % NumChans,   & ! in
           Work_Dimension,         & ! in
           RT_Params % H_Matrix,   & ! in
           RT_Params % H_Matrix_T, & ! in
           R_SubMatrix,            & ! in
           Delta_Profile,          & ! inout
           Status           )        ! out
   ELSE IF (Minimisation_Method == Marquardt_Levenberg) THEN
      Old_Gamma = Gamma
      CALL NWPSAF_Minimize_100ML( &
           BTObserved,             & ! in 
           RT_Params,              & ! inout
           BmatrixInverse,         & ! in
           Back,                   & ! in
           UsedChans,              & ! in
           Work_Dimension,         & ! in
           R_SubMatrix,            & ! inout
           DeltaBT,                & ! inout  
           Delta_Profile,          & ! inout
           Gamma,                  & ! inout
           JCost,                  & ! inout
           Status,                 & ! out
           Out_of_Range,           & ! out   
           Bmatrix          )        ! in
   ELSE 
      CALL NWPSAF_Minimize_100( &
           BmatrixInverse,         & ! in
           DeltaBT,                & ! in  
           UsedChans % NumChans,   & ! in
           Work_Dimension,         & ! in
           RT_Params % H_Matrix,   & ! in
           RT_Params % H_Matrix_T, & ! in
           R_SubMatrix,            & ! inout
           Delta_Profile,          & ! inout
           Status           )        ! out
   END IF
   
   IF ( Status /= 0 ) THEN
      Converged = .FALSE.
      EXIT Iteration_Loop
   END IF

   IF (Iterate == 1 .AND. ObNumber == 1 .AND. &
        GeneralMode >= DebugMode) THEN
      IF (Eqn_101) THEN
         write (*,*) 'Analysis_Error_Covariance not calculated for Eqn_101'
      ELSE
         CALL NWPSAF_SatMatInv(            &
              Num_ProfElementsUsed,      & ! in
              Num_ProfElementsUsed,      & ! in
              Analysis_Error_Covariance, & ! inout
              Status)                      ! out
         WRITE(FileUnit_AMatrix,*) '------------------------------'
         WRITE(FileUnit_AMatrix,*) 'Analysis Error Covariance:'
         WRITE(FileUnit_AMatrix,*) &
              'Calculated from the Background for Observation 1'
         WRITE(FileUnit_AMatrix,*) 'Number of Elements = ',Num_ProfElementsUsed
         WRITE(FileUnit_AMatrix,*) Analysis_Error_Covariance
         WRITE(FileUnit_AMatrix,*) 
         WRITE(FileUnit_AMatrix,*) '------------------------------'
         WRITE(FileUnit_AMatrix,*) '------------------------------'
      END IF
   END IF
    
   Guess(:) = Back(:) + Delta_Profile(:)
   RT_Params % RTGuess(Retrieved_Elements) = Guess(:)

   ! The first call of RTTOV uses the guess profile rather than the background profile.
   ! If the two are different, the first returned BTs will NOT correspond to the background!

   !2.2) QC new profile and check for convergence
   !----

!!!fradea
   CALL NWPSAF_CheckIteration( &
        RT_Params % RTGuess,      & ! inout
        Delta_Profile,            & ! inout
        Profile_Variables_Reset,  & ! out
        Out_of_Range,             & ! out
        Bmatrix              )      ! in
!!!fradea



   IF (Profile_Variables_Reset) THEN
      ! Need to reset Guess again
      Guess(:) = RT_Params % RTGuess(Retrieved_Elements)
      IF ( GeneralMode >= VerboseMode ) WRITE(*,*) &
           'Profile Variables Reset'
   END IF

   IF (GeneralMode >= DebugMode) THEN
      IF (Iterate == 1) THEN
         WRITE(FileUnit_MinimisationLog,*) '------------------------------'
         WRITE(FileUnit_MinimisationLog,*) 'ObNumber =',obnumber
         WRITE(FileUnit_MinimisationLog,*) 'Iteration =',0  
         WRITE(FileUnit_MinimisationLog,*) 'Gamma =',0.
         WRITE(FileUnit_MinimisationLog,*) 'BackGround Profile:'
         WRITE(FileUnit_MinimisationLog,*) RT_Params % RTBack(:)
         WRITE(FileUnit_MinimisationLog,*) 
         WRITE(FileUnit_MinimisationLog,*) 'Cost Function =',JCost
      END IF
      WRITE(FileUnit_MinimisationLog,*) '------------------------------'
      WRITE(FileUnit_MinimisationLog,*) 'ObNumber =',obnumber
      WRITE(FileUnit_MinimisationLog,*) 'Iteration =',Iterate  
      WRITE(FileUnit_MinimisationLog,*) 'Gamma =',Gamma
      WRITE(FileUnit_MinimisationLog,*) 'Profile:'
      WRITE(FileUnit_MinimisationLog,*) RT_Params % RTGuess
      WRITE(FileUnit_MinimisationLog,*) 
   END IF

!   IF (Out_of_Range) THEN
!        IF ( GeneralMode >= DebugMode ) WRITE(*,*) &
!             'Iteration out of range - exiting minimisation'
!        EXIT Iteration_Loop
!   END IF

   IF(MwClwRetrieval .AND. Lqtotal==0) THEN

      cloud_structure(:) = 0.0
      CALL NWPSAF_LWP_to_Layers(RT_Params % RTguess(Prof_ClW), RT_Params,cloud_structure )
      !The cloud liquid water profile is modified when using this LWP_firstGuess
   END IF

!===MODIFICATION PM========================================
   IF(IWVRetrieval) THEN
     q_structure(:) = 0.0
     q_structure(:) = &
        RT_Params % RTBack(Prof_FirstQ:Prof_LastQ)/RT_Params % RTBack(Prof_IWV)
        RT_Params % RTguess(Prof_FirstQ:Prof_LastQ) = RT_Params % RTguess(Prof_IWV) * q_structure(:)
     ! The cloud liquid water profile is modified when using this LWP_firstGuess
   END IF

   Fastmodel_Mode = FastmodelMode_Gradient
   WhichProf = GuessProf
   CALL NWPSAF_Fastmodel_interface(     &
        Fastmodel_Mode,          & ! in
        RT_Params,               & ! in
        WhichProf,               & ! in
        UsedChans,               & ! in
        ErrorCode)                 ! out


   DeltaBT = BTObserved(UsedChans % Channels(1 : UsedChans % NumChans)) - &
        RT_Params % TotalBTs(1 : UsedChans % NumChans)
   
   IF (GeneralMode >= DebugMode) THEN
      WRITE(FileUnit_BTMinimisationLog,*) '------------------------------'
      WRITE(FileUnit_BTMinimisationLog,*) 'ObNumber =',obnumber
      WRITE(FileUnit_BTMinimisationLog,*) 'Iteration =',Iterate  
      WRITE(FileUnit_BTMinimisationLog,*) 'Gamma =',Gamma
      WRITE(FileUnit_BTMinimisationLog,*) 'Brightness Temperature Difference:'
      WRITE(FileUnit_BTMinimisationLog,*) DeltaBT
      WRITE(FileUnit_BTMinimisationLog,*) 
   END IF

   
   !Check for convergence. All profile increments must satisfy the given
   !criteria for convergence to occur. If the increments are out of
   !sensible range then end minimization.
   

   CALL NWPSAF_Calculate_Cost_Function( &
        BMatrixInverse,         & ! in          
        R_SubMatrix,            & ! inout
        Delta_Profile,          & ! in          
        DeltaBT,                & ! in      
        Work_Dimension,         & ! in          
        Eqn_101,                & ! in          
        RT_Params % RTGuess,    & ! in
        JCost,                  & ! out         
        Status,                 & ! out         
        H_Matrix_T = RT_Params % H_Matrix_T, & ! optional in     
        JCost_Gradient = JCost_Gradient)       ! optional out

 
   IF (Eqn_101) THEN 
      ! JCost_Gradient is not currently produced in this case.
      Obs% Jcost_Gradient = 0.
   ELSE
      ! Calculate a scalar cost function gradient using the B-Matrix as a
      ! metric.  Ignore those components where the B-matrix is assumed to
      ! be large, i.e., CTP and CloudCover (see NWPSAF_InitBMatrix).
      IF (CloudyRetrieval) THEN
         WHERE(Retrieved_Elements == Prof_CTP .OR. &
              Retrieved_Elements == Prof_CloudCover) JCost_Gradient = 0.0
      END IF

      Obs% Jcost_Gradient = &
           DOT_PRODUCT(JCost_Gradient,MATMUL(BMatrix,JCost_Gradient))
   END IF

   IF (GeneralMode >= DebugMode) &
        WRITE(FileUnit_MinimisationLog,*) 'Cost Function =',JCost, &
        'Cost Function Gradient =',Obs% Jcost_Gradient

   ! Here we have the three convergence tests:
   ! 1) Cost function doesn't change much
   ! 2) Gamma not increasing (Marquardt-Levenberg only)
   ! 3) The cost function gradient isn't large
   ! 


   IF (ABS(JCost-Jold)/JCost < DeltaJ   &
        .AND. Gamma/Old_Gamma < 1.01      &
        .AND. ( Minimisation_Method /= Marquardt_Levenberg .OR. & 
        ABS(Obs % Jcost_Gradient) < SmallJCost_Gradient * JCost**2) ) THEN
      IF (Profile_Variables_Reset) THEN
         JCost=-JCost
         !!!fradea
         Converged = .FALSE.
         write(*,*)'profile variables reset'
         !!!fradea
      ELSE
         Converged = .TRUE.
         IF (GeneralMode >= VerboseMode) THEN
            WRITE(*,*) 'Convergence Criteria:'
            WRITE(*,*) 'JCost, JOld:',JCost, JOld
            WRITE(*,*) 'Gamma, Old_Gamma:',Gamma, Old_Gamma
            WRITE(*,*) 'Obs % Jcost_Gradient, SmallJCost_Gradient * JCost**2:',&
                 Obs % Jcost_Gradient, SmallJCost_Gradient * JCost**2
         END IF
      END IF
     
      EXIT Iteration_Loop
   END IF

   JOld = JCost

END DO Iteration_Loop

        
!--------------------------------
!3. Final checks and calculations
!--------------------------------


!3.1) Output Error Covariance Matrix to the end of the Minimisation log.
!----

IF (GeneralMode >= DebugMode) THEN
   IF (Eqn_101) THEN
      write (*,*) 'Analysis_Error_Covariance not calculated for Eqn_101'
   ELSE
      CALL NWPSAF_SatMatInv(            &
           Num_ProfElementsUsed,      &
           Num_ProfElementsUsed,      &
           Analysis_Error_Covariance, &
           Status)
      Prop_Measurement_Noise = &
           MATMUL(Prop_Measurement_Noise, Analysis_Error_Covariance)
      Prop_Measurement_Noise = &
           MATMUL(Analysis_Error_Covariance, Prop_Measurement_Noise)
      WRITE(FileUnit_AMatrix,*) '------------------------------'
      WRITE(FileUnit_AMatrix,*) 'Analysis Error Covariance:'
      WRITE(FileUnit_AMatrix,*) 'Observation Number = ',obnumber
      WRITE(FileUnit_AMatrix,*) 'Number of Elements = ',Num_ProfElementsUsed
      WRITE(FileUnit_AMatrix,*) Analysis_Error_Covariance
      WRITE(FileUnit_AMatrix,*) 
      WRITE(FileUnit_AMatrix,*) '------------------------------'
      WRITE(FileUnit_AMatrix,*) '------------------------------'
      WRITE(FileUnit_AmMatrix,*) '------------------------------'
      WRITE(FileUnit_AmMatrix,*) 'Propagated Measurement Noise Error Covariance:'
      WRITE(FileUnit_AmMatrix,*) 'Observation Number = ',obnumber
      WRITE(FileUnit_AmMatrix,*) 'Number of Elements = ',Num_ProfElementsUsed
      WRITE(FileUnit_AmMatrix,*) Prop_Measurement_Noise
      WRITE(FileUnit_AmMatrix,*) 
      WRITE(FileUnit_AmMatrix,*) '------------------------------'
      WRITE(FileUnit_AmMatrix,*) '------------------------------'
   END IF
END IF

!3.2) For Convergence
!----
!Calculate retrieval BTs and finish

IF ( Converged ) THEN

    Obs% Jcost = JCost / REAL(UsedChans % NumChans)

    IF ( GeneralMode >= VerboseMode ) THEN
       WRITE(*,'(A,F8.3,A,F8.3,A,F8.3)') &
            'Cost function: J/N = ', Obs% Jcost
       WRITE(*,'(A,F8.3)') &
            'Gradient: ', Obs% Jcost_Gradient
    END IF
    
ELSE 

   WRITE(*,*) '*********** Minimisation Failed ***********'
   IF (GeneralMode >= DebugMode) WRITE(FileUnit_MinimisationLog,*) &
       '*********** Minimisation Failed ***********'

END IF


Obs% NIter = Iterate

!-----------------------------------------
!4. Interpolate increments to model levels (currently not used)
!-----------------------------------------
!Note: the following routine might require some alteration before use.

DEALLOCATE(R_SubMatrix % Matrix)
DEALLOCATE(RT_Params % H_matrix)
DEALLOCATE(RT_Params % H_matrix_T)
DEALLOCATE(Analysis_Error_Covariance)
DEALLOCATE(Prop_Measurement_Noise)

End Subroutine NWPSAF_Minimize
