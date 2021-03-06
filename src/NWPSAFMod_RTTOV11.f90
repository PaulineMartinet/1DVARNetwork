!+ routines and constants for RTModel

MODULE NWPSAFMod_RTTOV11

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
! Description: contains information for RTTOV-11
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 31      18/01/13 Original module, based on RTTOV10 module. P. Weston.
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

USE rttov_types, ONLY : &
    rttov_coefs,        &
    rttov_options,      &
    rttov_chanprof

! Coefficient structure for RTTOV11
TYPE(rttov_coefs), POINTER :: RTTOV11_Coefs(:)

! Options list for RTTOV11
TYPE(rttov_options), POINTER :: RTTOV11_Opts(:)

! Channels and profiles structure for RTTOV11
TYPE(rttov_chanprof), POINTER :: RTTOV11_Chanprof(:)

Integer, PARAMETER :: jpch   = 10000    ! Max. no. of channels
Integer, PARAMETER :: jpnsat = 30       ! Max no sensors to be used


END MODULE NWPSAFMod_RTTOV11



