!+Compute Liquid Water Path sensitivity 
SUBROUTINE NWPSAF_LayerK_to_LWPK( &
     clw_k,  & ! In
     cloud_structure, & !in
     lwp_k)    ! Out
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
! Compute Liquid Water Path sensitivity as required by minimization routine.
!       dTB/dLWP
!
! History:
!
! Ticket Date     Comment
! ------- -------- -------
! 28     22/02/12 Original code based on SSMI_LayerK_to_LWPK.  TR Sreerekha.
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

USE NWPSAFMod_Constants, ONLY :  & 
    gravity

USE NWPSAFMod_RTmodel, ONLY : &
    Num_RTlevels

Use parkind1, Only : jpim,jprb
! Imported Type Definitions:

REAL (kind=jprb)           :: clw_k(Num_RTlevels) ! cloud liquid water on RTTOV levels
REAL            :: cloud_structure(Num_RTlevels)
REAL(kind=jprb)            :: lwp_k            ! Liquid water path

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_LayerK_to_LWPK"
INTEGER                      :: i               ! loop counter

!---- End of header --------------------------------------------------------

! Initialise output
write(*,*) 'initialise output'
  lwp_k= 0.0

! Loop over pressure levels

  DO i=2, Num_RTlevels    
     lwp_k = lwp_k +  0.5* (clw_k(i)+clw_k(i-1)) * Cloud_Structure(i)
  ENDDO



END SUBROUTINE NWPSAF_LayerK_to_LWPK
