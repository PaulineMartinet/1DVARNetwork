!+ Routine that calculates saturation water vapor pressure
Module NWPSAFMod_satvapor
!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2004, EUMETSAT, All Rights Reserved
! 
! Description:
!  Store's parameters for Tetens's formula that computes
!  saturation water vapor pressure.
!
!
! History:
! Ticket  Date     Comment
! 28      24/02/12 Based on SSMIMod_satvapor.  TR Sreerekha
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Declarations:
!

! Teten's formula parameters:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.

REAL,    PARAMETER  ::  awater = 17.269
REAL,    PARAMETER  ::  bwater = 35.86
REAL,    PARAMETER  ::  aice   = 21.875
REAL,    PARAMETER  ::  bice   = 7.66
REAL,    PARAMETER  ::  Tmelt  = 273.16
REAL,    PARAMETER  ::  Cteten = 610.78
LOGICAL, PARAMETER  ::  bothphases=.TRUE.


END MODULE NWPSAFMod_satvapor
