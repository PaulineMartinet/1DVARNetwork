!+ Calculate saturation vapour pressure
SUBROUTINE NWPSAF_svp (T,             & !In
                     Sat_Vap_Pres)    !Out
   
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
! Calculate saturation vapour pressure in Pa, given the temperature (K)
!
! History:
!
! Ticket Date     Comment
! ------- -------- -------
! 28     21/02/12 Original code based on SSMI_svp.  TR Sreerekha.
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

 USE NWPSAFMod_satvapor, ONLY :  &
     awater,                 &
     bwater,                 &
     aice,                   &
     bice,                   &
     Tmelt,                  &
     Cteten,                 &
     bothphases                       


 IMPLICIT NONE

! Arguments
 REAL            :: T             ! Temperature in Kelvin
 REAL            :: Sat_Vap_Pres  ! Pa       

! Continuous
 REAL            :: aconst
 REAL            :: bconst

!  Use Tetens' formula:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.
! Define constants to be used depending on phase (water or ice)

  aconst = awater
  bconst = bwater

  IF (bothphases.AND. T < Tmelt) THEN
        aconst = aice  
        bconst = bice  
  ENDIF

! Compute Sat_Vap_Pres in hPa         

  Sat_Vap_Pres = Cteten*exp ( aconst*(T-Tmelt)/(T-bconst) ) 


END SUBROUTINE NWPSAF_svp
 
