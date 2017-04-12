!+ Convert cloud liquid water to total liquid water path
SUBROUTINE NWPSAF_Layers_to_LWP( &
     clw,  & ! in
     RT_Params, & !in
     lwp)    ! out
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
! To sum all cloud liquid water from RTTOV levels to Total Liquid Water
! Path as required by descent routine. 
!
! History:
!
! Ticket    Date     Comment
! ------    ----     -------
! 28     22/02/12 Original code based on SSMI_Layers_to_LWP.  TR Sreerekha.
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Constants, ONLY :  & 
    gravity

USE NWPSAFMod_RTmodel, ONLY : &
    Num_RTlevels, &
    RTParams_Type

!Use parkind1, Only : jpim,jprb

IMPLICIT NONE

! Imported Type Definitions:

REAL       :: clw(Num_RTlevels) ! cloud liquid water on RTTOV levels
REAL            :: lwp             ! Liquid water path
TYPE(RTParams_Type), INTENT(INOUT):: RT_Params     ! Info for RT Model

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Layers_to_LWP"
INTEGER                      :: i               ! loop counter

!---- End of header --------------------------------------------------------

! Initialise output

  lwp = 0.0

! Loop over pressure levels

! Note that here LWP does not include 2m level clw.
  DO i=2,Num_RTlevels
     lwp = lwp +(RT_Params % Pressure_Pa(i)-RT_Params % Pressure_Pa(i-1)) &
             *(clw(i)+clw(i-1))
  
  ENDDO
  
  lwp=(0.5*lwp)/gravity
  
  
END SUBROUTINE NWPSAF_Layers_to_LWP
