Subroutine NWPSAF_Gastro_Interface( &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! inout
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrStatRep)                ! out
       
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
! Description: Dummy interface for when Gastropod is not required.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 3.3     26/02/04 Original version.   A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_Params, ONLY : &
     StatusFatal

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type

IMPLICIT NONE

INCLUDE 'Gen_ErrorReport.interface'

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  ! Mode in which 
                                                           ! the fastmodel is 
                                                           ! to be run
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrStatRep

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Gastro_Interface"

! Local variables
CHARACTER(LEN=80) :: ErrorMessage(1)  ! Message for Gen_ErrorReport
INTEGER :: Tmp                        ! Used to stop compiler warnings!
REAL :: TmpR                        ! Used to stop compiler warnings!

! Stop compiler warnings about variables not being used.

Tmp=FastModel_Mode
RT_Params % RTModel = Tmp
Tmp = WhichProf
TmpR = UsedChans % NumChans

! Fatal error message:

WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
     'Gastropod not supported'
ErrStatRep = StatusFatal
CALL Gen_ErrorReport( RoutineName,  &
     ErrorMessage, &                  
     ErrStatRep    )

End Subroutine NWPSAF_Gastro_Interface
