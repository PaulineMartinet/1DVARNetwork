!+Construct a liquid water profile from the up-dated LWP
Subroutine NWPSAF_LWP_to_Layers(&
     LWP_new, & !in
     RT_Params, & !inout
     cloud_structure)!inout

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
! Construct a liquid water profile from the up-dated LWP
! calculated from the NWPSAF_Minimize routine (which uses total column). 
!
! History:
!
! Ticket Date     Comment
! ------- -------- -------
!   28   21/02/12 Original code based on SSMI_LWP_to_Layers.  TR Sreerekha.
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY :  &
    DebugMode

USE NWPSAFMod_ObsInfo, ONLY :        &
     Ob_Type

USE NWPSAFMod_RTmodel, ONLY : &
    Num_RTlevels, &
    RTParams_Type, &
    Prof_FirstCLW, &
    Prof_LastCLW

IMPLICIT NONE

! Subroutine arguments:
TYPE(RTParams_Type)     , INTENT(INOUT) :: RT_Params ! RT Model Data
!TYPE(Ob_type)           , INTENT(INOUT) :: Obs       ! Observed/Retrieval data
REAL                :: cloud_structure(Num_RTlevels)
REAL                :: LWP_new

!local
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_LWP_to_Layers"


!----------------End of header -----------------------------------------

CALL NWPSAF_CloudStructure(RT_params, cloud_structure(:))

RT_Params % RTguess(Prof_FirstCLW:Prof_LastCLW) = LWP_new * cloud_structure(:)

End Subroutine NWPSAF_LWP_to_Layers
