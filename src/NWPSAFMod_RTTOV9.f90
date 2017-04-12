!+ routines and constants for RTModel

MODULE NWPSAFMod_RTTOV9

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
! Description: contains information for RTTOV-9
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 13      30/01/09 Original module.      E. Pavelin.
! 13      11/02/09 Added ozone_present.  E. Pavelin.
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

USE rttov_types, ONLY : &
    rttov_coef,         &
    rttov_coef_scatt_ir, &
    rttov_optpar_ir

! Coefficient structures for RTTOV9
TYPE(rttov_coef), POINTER :: RTTOV9_coef(:)
TYPE(rttov_coef_scatt_ir), POINTER :: RTTOV9_Coef_Scatt_IR(:)
TYPE(rttov_optpar_ir), POINTER :: optp(:)

! Whether or not to interpolate levels for each instrument
LOGICAL, POINTER :: addinterp(:)

! Whether or not ozone coefficients are present
LOGICAL          :: ozone_present

Integer, PARAMETER :: jpch   = 10000    ! Max. no. of channels
Integer, PARAMETER :: jpnsat = 30       ! Max no sensors to be used


END MODULE NWPSAFMod_RTTOV9



