!+ Split total water content (qtotal) into q and ql
Subroutine NWPSAF_Qtot_to_q_ql(&
     qtotal, & !In
     q,      & !Out
     ql,     & !Out
     t,      & !In
     RT_Params, & !In
     QtotOption) !In          
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
! IF QtotOption=1 : Split total water content (qtotal) into 
!   water vapor content (q) and (ql) cloud liquid water content.
! If QtotOption=0 : Compute derivatives: (q) =dq/dqtotal and (ql) = dql/dqtotal
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
!   28    24/02/12 Original code based on SSMI_Qtot_to_q_ql.f90  TR Sreerekha.
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

!!$USE SSMIMod_Params, ONLY :  &
!!$    Wet_levels,             &
!!$    Pres_level,             &
!!$    RH_1,                   &
!!$    RH_2,                   &
!!$    SplitQtotal,            &
!!$    DebugMode,              & 
!!$    DebugMode2,             & 
!!$    Reference_Pressures,    &
!!$    TopQ,                   &
!!$    epsilon_1000
    
USE NWPSAFMod_Params, ONLY :    &
     DebugMode, &
     Ret_FirstQ, &
     Ret_LastQ

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type, &
     Prof_FirstT, &
     Prof_LastT

USE NWPSAFMod_Constants, ONLY :  & 
    epsilon_1000,       &
    RH_1,                   &
    RH_2,                   &
    SplitQtotal

IMPLICIT NONE

!Subroutine Arguments
REAL                :: qtotal(Ret_LastQ - Ret_FirstQ + 1) ! Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
REAL                :: q(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: ql(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: t(Ret_LastQ - Ret_FirstQ + 1)
INTEGER             :: QtotOption
TYPE(RTParams_Type)     , INTENT(INOUT) :: RT_Params ! RT Model Data

!local
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_qtot_to_q_ql" 
INTEGER                      ::  i ! loop index
REAL                         :: epsilon = epsilon_1000/1000.0

REAL                :: qsat(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: RH_qtotal(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: qvfix(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: qexcess(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: wpress(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: esat(Ret_LastQ - Ret_FirstQ + 1)
REAL                :: SESAT
INTEGER             :: Num_WetLevels

!----------------End of header -----------------------------------------

!Compute saturated specific humidity profile (qsat) for Wet_levels only
 Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
 wpress(:)=RT_Params % Pressure_Pa(Prof_LastT - Num_WetLevels + 1:Prof_LastT)

 DO i=1,Num_WetLevels
    CALL NWPSAF_svp(t(i),SESAT) ! SESAT is in Pa.
    esat(i)=SESAT
 ENDDO

 WHERE (esat(:) > wPress(:))
    qsat(:) = 1.0
 ELSEWHERE
    qsat(:)= epsilon / ( wpress(:)/esat(:) - (1.-epsilon))
 ENDWHERE

 RH_qtotal = qtotal / qsat

 IF (QtotOption == 1) THEN

    ! Split qtotal into q and ql 

    qvfix = qsat * (RH_1 + SplitQtotal * (RH_2 - RH_1))      
    qexcess = qtotal - RH_1 * qsat

    WHERE ( RH_qtotal < RH_1 ) 
       q  = qtotal
    ENDWHERE

    WHERE ( RH_qtotal >= RH_1 .AND. RH_qtotal < RH_2 )
        q= RH_1 *qsat + SplitQtotal* qexcess
    ENDWHERE

    WHERE ( RH_qtotal >= RH_2 ) 
       q=qvfix
    ENDWHERE

    ql = qtotal - q


 ELSE IF (QtotOption == 0) THEN
  
    ! Compute derivates
    ! q = dq/dqtotal and ql = dql/dqtotal

    WHERE ( RH_qtotal < RH_1 )
        q  = 1.0       
        ql = 0.0
    ENDWHERE
    WHERE ( RH_qtotal >= RH_1 .AND. RH_qtotal < RH_2 ) 
        q = SplitQtotal
        ql = (1. - SplitQtotal)
    ENDWHERE
    WHERE ( RH_qtotal >= RH_2 )
        q  = 0.0       
        ql = 1.0
    ENDWHERE

 ELSE IF (QtotOption == 2) THEN

    ! Compute derivatives of q and ql with respect to qsat
    !   q = dq/dqsat  and ql = dql/dqsat

    WHERE ( RH_qtotal(:) < RH_1 )
        q(:)  = 0.0
    ENDWHERE
    WHERE ( RH_qtotal(:) >= RH_1 .AND. RH_qtotal(:) < RH_2 )
        q(:) = (1. - SplitQtotal) *  RH_1
    ENDWHERE
    WHERE ( RH_qtotal(:) >= RH_2 )
        q(:)  = (RH_1 + SplitQtotal * (RH_2 - RH_1))
    ENDWHERE

    ql(:) = -q(:)

 ELSE

    WRITE(6,*) ' QtotOption not well defined'
    WRITE(6,*) ' QtotOption =', QtotOption
    STOP

 ENDIF

End Subroutine NWPSAF_Qtot_to_q_ql 
