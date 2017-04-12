!+ Calculate derivative of ln(saturation vapour pressure) w.r.t T 
SUBROUTINE NWPSAF_svp_deriv (T,             & !In
                           dlnes_dT)        !Out
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
! Calculates derivative of natural logarithm of saturation vapour pressure 
! with respect to temperature.
! T is in Kelvin, es in Pa.
!
! Calculate saturation vapour pressure in Pa, given the temperature (K)
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     24/02/12 Original code based on SSMI_svp_deriv.  TR Sreerekha.
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

!Arguments
REAL            :: T  ! Temperature in Kelvin
REAL            :: dlnes_dT

REAL            :: aconst
REAL            :: bconst
REAL            :: aw           ! T-bconst

!  Use Tetens' formula:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.
!
!  Define constants to be used depending on phase (water or ice)

aconst = awater
bconst = bwater

IF (T < Tmelt.AND.bothphases) THEN
   aconst = aice
   bconst = bice
ENDIF

! Compute Dlnes_DT

aw = T-bconst
dlnes_dT = aconst*( Tmelt-bconst )/( aw*aw )


END SUBROUTINE NWPSAF_svp_deriv
 
