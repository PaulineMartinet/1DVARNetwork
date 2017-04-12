!+ routines and constants for RTModel

MODULE NWPSAFMod_RTTOV8

!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2009, EUMETSAT, All Rights Reserved.
!
! Description: contains information for RTTOV-8
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 13      30/01/09 Original module.      E. Pavelin.
!                  !
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

USE rttov_types, ONLY : &
    rttov_coef

! Coefficient structures for RTTOV8
TYPE(rttov_coef), POINTER :: RTTOV8_coef(:)

END MODULE NWPSAFMod_RTTOV8



