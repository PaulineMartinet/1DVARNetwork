Subroutine NWPSAF_Fastmodel_Interface( &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! inout
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrorCode)                 ! out
       
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
! Description: Interface between NWPSAF 1DVar code and Fastmodels
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     03/06/99 Original code started. A. Collard.
! 1.1     03/02/00 Reorder so that the fastmodel to be used is chosen
!                  first and then the access mode.
! 1.2     24/03/00 Add test to RTIASI Jacobian for spikes.
! 1.3     30/01/01 Changed to interface block format.
! 1.4     21/03/01 Addition of RTTOV Code.  A. Collard.
! 1.5     04/10/01 Changes to add in RTTOV7.  A. Collard.
! 1.6     12/11/01 Tidy up and add code to allow processing of ATOVS as
!                  a single instrument (RTTOV7 requires HIRS,AMSU-A 
!                  & AMSU-B separately)
! 3.0.3   26/02/04 Add Gastropod call (V. Sherlock)     A. Collard.
! 3.1.0   05/04/05 Add RTTOV8 calls                     E. Pavelin.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date     Comment
! 13      29/01/09 Added RTTOV9 Calls. E. Pavelin.
! 25      20/02/12 Added RTTOV10 calls. P. Weston.
! 31      18/01/12 Added RTTOV11 calls. P. Weston.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type,            &
     RTIASI,                   &
     RTTOV,                    &
     RTTOV8,                   &
     RTTOV9,                   &
     RTTOV10,                  &
     RTTOV11,                  &
     Gastropod 
 
USE NWPSAFMod_Params, ONLY : &
     StatusFatal
       
IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'
INCLUDE 'NWPSAF_Gastro_Interface.interface'
INCLUDE 'NWPSAF_RTTOV7_Interface.interface'
INCLUDE 'NWPSAF_RTTOV8_Interface.interface'
INCLUDE 'NWPSAF_RTTOV9_Interface.interface'
INCLUDE 'NWPSAF_RTTOV10_Interface.interface'
INCLUDE 'NWPSAF_RTTOV11_Interface.interface'
INCLUDE 'NWPSAF_RTIASI_Interface.interface'

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params    ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf    ! Which profiel to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrorCode


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Fastmodel_Interface"

! Local variables:

INTEGER :: ErrStatRep                 ! error status for Gen_ErrorReport
CHARACTER(LEN=80) :: Message(2)       ! Message for Gen_MessageReport

!---------------------------------------------------------------

! Error code is not used at present - set to zero here

ErrorCode = 0

SELECT CASE (RT_Params % RTModel)

   !-----------------------------------------------------------
   ! This is the interface with RTTOV8
   !----------------------------------------------------------
   CASE(RTTOV8)
      CALL NWPSAF_RTTOV8_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with RTTOV7
   !----------------------------------------------------------
   CASE(RTTOV)
      CALL NWPSAF_RTTOV7_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with RTTOV9
   !----------------------------------------------------------
   CASE(RTTOV9)
      CALL NWPSAF_RTTOV9_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with RTTOV10
   !----------------------------------------------------------
   CASE(RTTOV10)
      CALL NWPSAF_RTTOV10_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with RTTOV11
   !----------------------------------------------------------
   CASE(RTTOV11)
      CALL NWPSAF_RTTOV11_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with RTIASI
   !----------------------------------------------------------
   CASE(RTIASI)
      CALL NWPSAF_RTIASI_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! This is the interface with Gastropod 
   !----------------------------------------------------------
   CASE(Gastropod)
      CALL NWPSAF_Gastro_Interface( &
           Fastmodel_Mode,          & ! in
           RT_Params,               & ! inout
           WhichProf,               & ! in
           UsedChans,               & ! in
           ErrorCode)                 ! out
   !-----------------------------------------------------------
   ! Else we have a problem
   !----------------------------------------------------------
   CASE DEFAULT
      Message(1) = 'Required RT Model is not supported'
      ErrStatRep = StatusFatal
      CALL Gen_ErrorReport( RoutineName, &
           Message(1:1), &
           ErrStatRep )
END SELECT

  
End Subroutine NWPSAF_Fastmodel_Interface
