SUBROUTINE RTIASI_RTIAI(KERR,KIDSAT,KNSAT,KNPF)

!     Description:
!     Initialisation for IASI radiative transfer routine, RTIASI.
!     To be called before first call to RTIASI.
!
!     Method:
!     Simple
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            10/01/1991  Original code (RTTOV). J.R.Eyre. ECMWF.
!     2            12/11/1999  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_PRFCON,ONLY : &
!     Imported parameters:
 JPLEV   ,   & ! Number of pressure levels
!     Imported scalars:
  NLEVW  ,   & ! Upper level for water vapour transmittance calculation
!                  ! (assumed transmittance=1. above this)
  UPPRES ,   & ! Upper level for transmittance calculation
!     Imported arrays:
  XPRES  ,   & ! Standard pressure levels for transmittance (and, currently,
!                  ! radiative transfer) calculation; from top down; in mb=hpa
  XPRES2 ,   & ! XPRES**2; in mb**2
  DPRES  ,   & ! Intervals between standard pressure levels; in mb
  DPP        ! DPRES*XPRES


USE RTIASI_GEOCON,ONLY : &
!     Imported scalars:
 HTSAT   ,   & ! Satellite altitude in km
 RATOE       ! Ratio satellite-orbit/earth radii


IMPLICIT NONE

!     SUBROUTINE RTIASI_arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) ::  KNSAT          ! Number of satellites to be
                                      ! processed.
INTEGER,INTENT(IN) :: KNPF            ! Number of processed profiles.
INTEGER,INTENT(IN) ::  KIDSAT (KNSAT) ! Array of satellite
                                      ! identification indices.
 
!       Scalar arguments with intent out:
INTEGER,INTENT(OUT)::  KERR           ! Error flag

!     End of subroutine arguments


!       Local scalars:
INTEGER ::  J,JSAT,ISAT
INTEGER ::  ILEVW          ! Top level for water vap. transmittance calc.
INTEGER ::  IPLEV          ! Number of standard pressure levels.
INTEGER ::  INDMAX         ! Max. value of satellite index.
REAL    ::  TOPPRES        ! Top pressure used for transmittance calc.
REAL    ::  REARTH         ! Radius of earth

!       Local arrays:
REAL    ::  ZPRES  (JPLEV) ! Standard 43 pressure levels
REAL    ::  ZHTSAT (KNSAT) ! Eight of satellite orbit


!-----End of header-----------------------------------------------------


DATA ZPRES       /.1,.29,.69,1.42,2.611,4.407,6.95,10.37,14.81, &
20.4,27.26,35.51,45.29,56.73,69.97,85.18,102.05,122.04, &
143.84,167.95,194.36,222.94,253.71,286.6,321.5,358.28, &
396.81,436.95,478.54,521.46,565.54,610.6,656.43,702.73, &
749.12,795.09,839.95,882.8,922.46,957.44,985.88,1005.43,1013.25/
DATA TOPPRES     /0.004985/
DATA ILEVW       /1/
DATA REARTH      /6370.949/
DATA IPLEV       /43/
DATA INDMAX      /1/


!-----------------------------------------------------------------------
!         1.   SET UP PROFILE CONSTANTS.
!-----------------------------------------------------------------------

KERR=0

  DO ISAT=1,KNSAT
    ZHTSAT(ISAT)=826.99
  ENDDO


!-----1.1  Set up pressure level constants------------------------------

UPPRES=TOPPRES
  IF(IPLEV /= JPLEV)THEN
    KERR=1
    RETURN
  ENDIF
    DO J=1,IPLEV
      XPRES(J)=ZPRES(J)
      XPRES2(J)=XPRES(J)*XPRES(J)
    ENDDO
DPRES(1)=XPRES(1)-UPPRES
DPP(1)=DPRES(1)*UPPRES
  DO J=2,IPLEV
    IF(J == 2)THEN
      DPRES(J)=XPRES(J-1)-UPPRES
    ELSE IF(J > 2)THEN
      DPRES(J)=XPRES(J-1)-XPRES(J-2)
    END IF
    DPP(J)=DPRES(J)*XPRES(J-1)
  ENDDO
!-----------------------------------------------------------------------


!-----1.2  Set up upper level humidity constants------------------------

NLEVW=ILEVW

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!         2.   SET UP DATA FOR NEW SATELLITE.
!-----------------------------------------------------------------------


!-----2.1  Set up satellite geometry constants--------------------------

!     (sat index: IASI=1)
DO JSAT=1,KNSAT
  ISAT=KIDSAT(JSAT)
  IF (ISAT > INDMAX) THEN
    KERR=11
    RETURN
  ENDIF
  HTSAT(JSAT)=ZHTSAT(ISAT)
  RATOE(JSAT)=(REARTH+HTSAT(JSAT))/REARTH
END DO
!-----------------------------------------------------------------------


!-----2.2  Set up satellite-specific data for IASI----------------------

CALL RTIASI_IASCF(KERR,KIDSAT,KNSAT,KNPF)

!-----------------------------------------------------------------------


RETURN
END SUBROUTINE RTIASI_RTIAI
!     Initialise satellite-dependent data for IASI.
SUBROUTINE RTIASI_IASCF(KERR,KIDSAT,KNSAT,KNPF)

!     Description:
!     Initialise satellite-dependent data and coefficients
!     for IASI radiative transfert routines.
!
!     Method:
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            10/12/1990  Original code(RTTOV). J.R.Eyre. ECMWF.
!     2            12/11/1999  Store data for IASI only. Marco Matricardi. ECMWF
!     2.1          06/09/2001  Make code flexible so that unformatted data can
!                              be read in independent of default precision for
!                              the code.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE NWPSAFMod_Params, ONLY : &
     Coeffs_Dir

USE RTIASI_TAUCFN, ONLY : &
!     Imported parameters:
  JPCOFM  ,  & ! Max. number of mixed gas coefficients
  JPCOFW  ,  & ! Max. number of water vapour coefficients
  JPCOFO  ,  & ! Max. number of ozone coefficients
!     Imported arrays
  CFMPOS  ,  & ! Mixed gases coefficients(positive transmittances)
  CFWPOS  ,  & ! Water vapour coefficients(single regression scheme)
  CFWPOS1 ,  & ! Water vapour coefficients(optically thin regime)
  CFWPOS2 ,  & ! Water vapour coefficients(optically thick regime)
  CFO        ! Ozone coefficients


USE RTIASI_IASCHN,ONLY : &
!     Imported parameters:
  JPNSAT  ,  & ! Max. number of satellite coefficients stored
  JPCH    ,  & ! Max. number of IASI channels
!     Imported arrays:
  WVNUM   ,  & ! Wavenumber of IASI channel in cm**-1
  GAMMA   ,  & ! "Gamma factor" transmittance corrections
  EMSCOEF    ! Regression coefficients to compute emissivity over water

USE RTIASI_EMISIR,ONLY : &
!     Imported arrays:
  EMSIR      ! Infrared emissivities

USE RTIASI_PRFREF,ONLY : &
!     Imported parameters:
  JPLEV   ,  & ! Number of pressure levels
!     Imported arrays:
  TREF    ,  & ! Reference temperature profile in K
  TREFO   ,  & ! Reference temperature for ozone in K
  WREF    ,  & ! Reference water vapour volume mixing ratio in ppmv
  WREFO   ,  & ! Reference water vapour for ozone in ppmv
  OREF       ! Reference ozone volume mixing ratio in ppmv


IMPLICIT NONE

INCLUDE 'NWPSAF_OpenFile.interface'

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN)   :: KNSAT          ! Number of processed satellites

!       Array arguments with intent in:
INTEGER,INTENT(IN)   :: KIDSAT (KNSAT) ! Array of satellite indices

INTEGER,INTENT(IN)   :: KNPF           ! Number of processed profiles.

!       Scalar arguments with intent out:
INTEGER,INTENT(INOUT)  :: KERR           ! Error label

!     End of subroutine arguments


!       Local parameters:
INTEGER,PARAMETER    ::  NBYTES_SP = SELECTED_REAL_KIND ( 6, 37)
INTEGER,PARAMETER    ::  NBYTES_DP = SELECTED_REAL_KIND (14,307)

!       Local scalars:
INTEGER :: JSAT,JK,KK,JCH,I,JI,JJ,J,ICHAN,IUEMS,KEMS
INTEGER :: ICH          ! Number of channel coefficients to be read
INTEGER :: ICOFM        ! Number of mixed gases coefficients to be read
INTEGER :: ICOFW        ! Number of water vap. coefficients to be read
INTEGER :: ICOFO        ! Number of ozone coefficients to be read
INTEGER :: ILEV         ! Number of pressure levels
INTEGER :: IUMIXPOS     ! File unit number
!INTEGER :: IUMIXNEG     ! File unit number
INTEGER :: IUH2O        ! File unit number
INTEGER :: IUH2O1       ! File unit number
INTEGER :: IUH2O2       ! File unit number
INTEGER :: IUOZ         ! File unit number
INTEGER :: IUWAOPC      ! File unit number
INTEGER :: IUGAM        ! File unit number
INTEGER :: IRC          ! Record number
INTEGER :: ISAT
INTEGER :: Precision    ! Number of bytes in single precision real
REAL    :: W

!       Local arrays

REAL(KIND=NBYTES_SP), ALLOCATABLE :: CFMPOS_IN  (:,:)
REAL(KIND=NBYTES_SP), ALLOCATABLE :: CFWPOS_IN  (:,:) 
REAL(KIND=NBYTES_SP), ALLOCATABLE :: CFWPOS1_IN (:,:) 
REAL(KIND=NBYTES_SP), ALLOCATABLE :: CFWPOS2_IN (:,:) 
REAL(KIND=NBYTES_SP), ALLOCATABLE :: CFO_IN     (:,:) 
REAL(KIND=NBYTES_SP), ALLOCATABLE :: TREF_IN  (:)
REAL(KIND=NBYTES_SP), ALLOCATABLE :: TREFO_IN (:) 
REAL(KIND=NBYTES_SP), ALLOCATABLE :: WREF_IN  (:)
REAL(KIND=NBYTES_SP), ALLOCATABLE :: WREFO_IN (:)
REAL(KIND=NBYTES_SP), ALLOCATABLE :: OREF_IN  (:)

REAL(KIND=NBYTES_DP) :: EMSCOEF_IN(9)  ! Input emissivity coefficients
                         ! *Could be allocatable but only 9 elements*

! Local variables for NWPSAF_OpenFile

CHARACTER(LEN=80) :: Filename
CHARACTER(LEN=10) :: Access = "DIRECT"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "UNFORMATTED"
CHARACTER(LEN=3)  :: Status = "OLD"
INTEGER :: RecordLength

!-----End of header-----------------------------------------------------


DATA ICH      /8461/
DATA ICOFM    /14/
DATA ICOFW    /14/
DATA ICOFO    /14/
DATA ILEV     /43/
!DATA IUMIXPOS /13/
!DATA IUMIXNEG /14/
!DATA IUH2O    /15/
!DATA IUH2O1   /16/
!DATA IUH2O2   /17/
!DATA IUOZ     /18/
!DATA IUWAOPC  /19/
!DATA IUGAM    /20/


!-----------------------------------------------------------------------
! 
! 0.1 Check current precision and set up variables for I/O
!
!-----------------------------------------------------------------------

IF (KIND(1.0) == NBYTES_SP) THEN 
   Precision = 4
ELSE IF (KIND(1.0) == NBYTES_DP) THEN
   Precision = 8
   ALLOCATE( CFMPOS_IN  (JPCOFM,JPLEV) ) 
   ALLOCATE( CFWPOS_IN  (JPCOFW,JPLEV) ) 
   ALLOCATE( CFWPOS1_IN (JPCOFW,JPLEV) ) 
   ALLOCATE( CFWPOS2_IN (JPCOFW,JPLEV) ) 
   ALLOCATE( CFO_IN     (JPCOFO,JPLEV) ) 
   ALLOCATE( TREF_IN  (JPLEV) ) 
   ALLOCATE( TREFO_IN (JPLEV) ) 
   ALLOCATE( WREF_IN  (JPLEV) ) 
   ALLOCATE( WREFO_IN (JPLEV) ) 
   ALLOCATE( OREF_IN  (JPLEV) ) 
ELSE               ! This would be a little unusual
   WRITE(*,*) 'WARNING: Unexpected precision'
   Precision = 0
   ALLOCATE( CFMPOS_IN  (JPCOFM,JPLEV) ) 
   ALLOCATE( CFWPOS_IN  (JPCOFW,JPLEV) ) 
   ALLOCATE( CFWPOS1_IN (JPCOFW,JPLEV) ) 
   ALLOCATE( CFWPOS2_IN (JPCOFW,JPLEV) ) 
   ALLOCATE( CFO_IN     (JPCOFO,JPLEV) ) 
   ALLOCATE( TREF_IN  (JPLEV) ) 
   ALLOCATE( TREFO_IN (JPLEV) ) 
   ALLOCATE( WREF_IN  (JPLEV) ) 
   ALLOCATE( WREFO_IN (JPLEV) ) 
   ALLOCATE( OREF_IN  (JPLEV) ) 
END IF

ALLOCATE(CFMPOS(JPCOFM,JPLEV,JPCH,JPNSAT))
ALLOCATE(CFWPOS(JPCOFW,JPLEV,JPCH,JPNSAT)) 
ALLOCATE(CFWPOS1(JPCOFW,JPLEV,JPCH,JPNSAT))
ALLOCATE(CFWPOS2(JPCOFW,JPLEV,JPCH,JPNSAT))
ALLOCATE(CFO(JPCOFO,JPLEV,JPCH,JPNSAT))

!-----------------------------------------------------------------------
!        1.   OPEN AND READ IASI MIXED GAS TRANS COEFFS AND REF PROFILE.
!             Positive transmittances.
!-----------------------------------------------------------------------
IF (ICH /= JPCH) THEN
  KERR=3
  RETURN
ENDIF

IF (ILEV /= JPLEV) THEN
  KERR=4
  RETURN
ENDIF

IF (ICOFM /= JPCOFM) THEN
  KERR=5
  RETURN
ENDIF

!-----! f77 w/s---------------------------------------------------------
!      OPEN (IUMIXPOS,STATUS='OLD',ACCESS='DIRECT',RECL=JPLEV*JPCOFM,
!-----! f90 w/s---------------------------------------------------------
!OPEN (IUMIXPOS,STATUS='OLD',ACCESS='DIRECT',RECL=4*JPLEV*JPCOFM, &
!FILE='IASIMIXCOEFv1.dat')

Filename = Trim(Coeffs_Dir) // 'IASIMIXCOEFv1.dat'
RecordLength = 4*JPLEV*JPCOFM

CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUMIXPOS,            & ! out
     RecordLength )         ! optional in


IF (Precision == 4) THEN
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      IF (ISAT < 1.OR.ISAT > JPNSAT) THEN
         KERR=9
         RETURN
      ENDIF
            
      DO JK=1,ICH
         IRC=JK
         IF (JK == 1) THEN
            WRITE(*,'(/)')
            WRITE(*,'(A32)')'READING MIXED GASES COEFFICIENTS'
         END IF
!---------Note that different reference profiles are read for water-----
!         vapour and ozone respectively.
         READ (UNIT=IUMIXPOS,REC=IRC) &
              ((CFMPOS(JI,JJ,JK,JSAT),JI=1,ICOFM),JJ=1,ILEV)
      ENDDO

      WRITE(*,'(A38)')'MIXED GASES COEFFICIENTS READ:POSITIVE'
      
      IRC=IRC+1
      READ (UNIT=IUMIXPOS,REC=IRC) (TREF(J),J=1,ILEV), &
           (WREF(J), J=1,ILEV),(TREFO(J),J=1,ILEV),(WREFO(J),J=1,ILEV), &
           (OREF(J),J=1,ILEV)
      
      WRITE(*,'(/)')
      WRITE(*,'(A22)')'REFERENCE PROFILE READ'
      WRITE(*,'(/)')
      
      WRITE(*,'(6A13)')'TREFWV','WREFWV','TREFOZ','WREFOZ','OREF','LEVEL'
      WRITE(*,'(/)')
      
      DO KK=1,43
         WRITE(*,'(5(F13.5),I13)')TREF(KK),WREF(KK),TREFO(KK), &
              WREFO(KK),OREF(KK),KK
      END DO
      
      WRITE(*,'(/)')     
   ENDDO
ELSE
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      IF (ISAT < 1.OR.ISAT > JPNSAT) THEN
         KERR=9
         RETURN
      ENDIF


      DO JK=1,ICH
         IRC=JK
         IF (JK == 1) THEN
            WRITE(*,'(/)')
            WRITE(*,'(A32)')'READING MIXED GASES COEFFICIENTS'
         END IF
!---------Note that different reference profiles are read for water-----
!         vapour and ozone respectively.
         READ (UNIT=IUMIXPOS,REC=IRC) &
              ((CFMPOS_IN(JI,JJ),JI=1,ICOFM),JJ=1,ILEV)
         CFMPOS(:,:,JK,JSAT) = CFMPOS_IN(:,:)
      ENDDO

      WRITE(*,'(A38)')'MIXED GASES COEFFICIENTS READ:POSITIVE'
      
      IRC=IRC+1
      READ (UNIT=IUMIXPOS,REC=IRC) (TREF_IN(J),J=1,ILEV), &
           (WREF_IN(J), J=1,ILEV),(TREFO_IN(J),J=1,ILEV), &
           (WREFO_IN(J),J=1,ILEV), &
           (OREF_IN(J),J=1,ILEV)
      
      WRITE(*,'(/)')
      WRITE(*,'(A22)')'REFERENCE PROFILE READ'
      WRITE(*,'(/)')
      
      WRITE(*,'(6A13)')'TREFWV','WREFWV','TREFOZ','WREFOZ','OREF','LEVEL'
      WRITE(*,'(/)')
      
      DO KK=1,43
         WRITE(*,'(5(F13.5),I13)')TREF_IN(KK),WREF_IN(KK),TREFO_IN(KK), &
              WREFO_IN(KK),OREF_IN(KK),KK
      END DO
      
      WRITE(*,'(/)')
      
   ENDDO
   
   TREF  (:) = TREF_IN  (:)  
   TREFO (:) = TREFO_IN (:)  
   WREF  (:) = WREF_IN  (:)  
   WREFO (:) = WREFO_IN (:)  
   OREF  (:) = OREF_IN  (:)  

ENDIF
!-----------------------------------------------------------------------


IF (ICOFW /= JPCOFW) THEN
  KERR=6
  RETURN
ENDIF



!-----------------------------------------------------------------------
!        2.   OPEN AND READ IASI WATER VAPOUR TRANSMITTANCE COEFFS.
!-------------Positive transmittances:optically thin regime ------------
!-----------------------------------------------------------------------

!-----! f77 w/s---------------------------------------------------------
!     OPEN (IUH2O1,STATUS='OLD',ACCESS='DIRECT',RECL=JPLEV*JPCOFW,
!-----! f90 w/s---------------------------------------------------------
!OPEN (IUH2O1,STATUS='OLD',ACCESS='DIRECT',RECL=4*JPLEV*JPCOFW, &
!FILE='IASIWVCOEFRG1v1.dat')

Filename = Trim(Coeffs_Dir) // 'IASIWVCOEFRG1v1.dat'
RecordLength = 4*JPLEV*JPCOFW
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUH2O1,              & ! out
     RecordLength )         ! optional in

IF (Precision == 4) THEN
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O1,REC=IRC) &
              ((CFWPOS1(JI,JJ,JK,JSAT),JI=1,ICOFW),JJ=1,ILEV)
      ENDDO
      WRITE(*,'(A52)')'WATER VAPOUR COEFFICIENTS READ:OPTICALLY THIN REGIME'
   ENDDO
ELSE
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O1,REC=IRC) &
              ((CFWPOS1_IN(JI,JJ),JI=1,ICOFW),JJ=1,ILEV)
         CFWPOS1(:,:,JK,JSAT) = CFWPOS1_IN(:,:)
      ENDDO
      WRITE(*,'(A52)')'WATER VAPOUR COEFFICIENTS READ:OPTICALLY THIN REGIME'
   ENDDO
END IF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!        3.   OPEN AND READ IASI WATER VAPOUR TRANSMITTANCE COEFFS.
!-------------Positive transmittances:optically thick regime -----------
!-----------------------------------------------------------------------

!-----! f77 w/s---------------------------------------------------------
!     OPEN (IUH2O2,STATUS='OLD',ACCESS='DIRECT',RECL=JPLEV*JPCOFW,
!-----! f90 w/s---------------------------------------------------------
!OPEN (IUH2O2,STATUS='OLD',ACCESS='DIRECT',RECL=4*JPLEV*JPCOFW, &
!FILE='IASIWVCOEFRG2v1.dat')

Filename = Trim(Coeffs_Dir) // 'IASIWVCOEFRG2v1.dat'
RecordLength = 4*JPLEV*JPCOFW
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUH2O2,              & ! out
     RecordLength )         ! optional in

IF (Precision == 4) THEN
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O2,REC=IRC) &
              ((CFWPOS2(JI,JJ,JK,JSAT),JI=1,ICOFW),JJ=1,ILEV)
      ENDDO
      WRITE(*,'(A53)')'WATER VAPOUR COEFFICIENTS READ:OPTICALLY THICK REGIME'
   ENDDO
ELSE
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O2,REC=IRC) &
              ((CFWPOS2_IN(JI,JJ),JI=1,ICOFW),JJ=1,ILEV)
         CFWPOS2(:,:,JK,JSAT) = CFWPOS2_IN(:,:)
      ENDDO
      WRITE(*,'(A53)')'WATER VAPOUR COEFFICIENTS READ:OPTICALLY THICK REGIME'
   ENDDO
END IF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!        4.   OPEN AND READ IASI WATER VAPOUR TRANSMITTANCE COEFFS.
!-------------Positive transmittances:single regression scheme----------
!-----------------------------------------------------------------------

!-----! f77 w/s---------------------------------------------------------
!     OPEN (IUH2O,STATUS='OLD',ACCESS='DIRECT',RECL=JPLEV*JPCOFW,
!-----! f90 w/s---------------------------------------------------------
!OPEN (IUH2O ,STATUS='OLD',ACCESS='DIRECT',RECL=4*JPLEV*JPCOFW, &
!FILE='IASIWVCOEFv1.dat')

Filename = Trim(Coeffs_Dir) // 'IASIWVCOEFv1.dat'
RecordLength = 4*JPLEV*JPCOFW
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUH2O,               & ! out
     RecordLength )         ! optional in

IF (Precision == 4) THEN
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O,REC=IRC) &
              ((CFWPOS(JI,JJ,JK,JSAT),JI=1,ICOFW),JJ=1,ILEV)
      ENDDO
      WRITE(*,'(A55)')'WATER VAPOUR COEFFICIENTS READ:SINGLE REGRESSION SCHEME'
   ENDDO
ELSE
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A33)')'READING WATER VAPOUR COEFFICIENTS'
         END IF
         
         READ (UNIT=IUH2O,REC=IRC) &
              ((CFWPOS_IN(JI,JJ),JI=1,ICOFW),JJ=1,ILEV)
         CFWPOS(:,:,JK,JSAT) = CFWPOS_IN(:,:)
      ENDDO
      WRITE(*,'(A55)')'WATER VAPOUR COEFFICIENTS READ:SINGLE REGRESSION SCHEME'
   ENDDO
END IF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!        5.   OPEN AND READ IASI OZONE TRANSMITTANCE COEFFICIENTS.
!-----------------------------------------------------------------------
IF (ICOFO /= JPCOFO) THEN
  KERR=7
  RETURN
ENDIF

!-----! f77 w/s---------------------------------------------------------
!     OPEN (IUOZ,STATUS='OLD',ACCESS='DIRECT',RECL=JPLEV*JPCOFO,
!-----! f90 w/s---------------------------------------------------------
!OPEN (IUOZ,STATUS='OLD',ACCESS='DIRECT',RECL=4*JPLEV*JPCOFO, &
!     FILE='IASIOZCOEFv1.dat')

Filename = Trim(Coeffs_Dir) // 'IASIOZCOEFv1.dat'
RecordLength = 4*JPLEV*JPCOFO

CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUOZ,                & ! out
     RecordLength )         ! optional in

IF (Precision == 4) THEN
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A26)')'READING OZONE COEFFICIENTS'
         END IF
         
         READ (UNIT=IUOZ,REC=IRC) &
              ((CFO(JI,JJ,JK,JSAT),JI=1,ICOFO),JJ=1,ILEV)
      ENDDO
      WRITE(*,'(A23)')'OZONE COEFFICIENTS READ'
   ENDDO
ELSE 
   DO JSAT=1,KNSAT
      ISAT=KIDSAT(JSAT)
      DO JK=1,ICH
         IRC=JK
         
         IF(JK == 1)THEN
            WRITE(*,'(/)')
            WRITE(*,'(A26)')'READING OZONE COEFFICIENTS'
         END IF
         
         READ (UNIT=IUOZ,REC=IRC) &
              ((CFO_IN(JI,JJ),JI=1,ICOFO),JJ=1,ILEV)
         CFO(:,:,JK,JSAT) = CFO_IN(:,:)
      ENDDO
      WRITE(*,'(A23)')'OZONE COEFFICIENTS READ'
   ENDDO
END IF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!        6.   STORE IASI GAMMAS AND CHANNEL WAVENUMBERS
!-----------------------------------------------------------------------

WRITE(*,'(/)')
WRITE(*,'(A14)')'READING GAMMAS'
!OPEN(IUGAM,FILE='IASIGAMMAS.dat')
Filename = Trim(Coeffs_Dir) // 'IASIGAMMAS.dat'
Access='SEQUENTIAL'
Form='FORMATTED'
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUGAM)                 ! out

  DO JSAT=1,KNSAT
    DO JCH=1,JPCH
      READ(IUGAM,*)W,GAMMA(JCH,JSAT)     ! Gammas are currently set to 1
      WVNUM(JCH)=645.+0.25*(JCH-1)
    ENDDO
  ENDDO
  WRITE(*,'(A11)')'GAMMAS READ'

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!        7.   STORE REGRESSION COEFFICIENTS TO COMPUTE IASI CHANNELS
!             EMISSIVITIES
!-----------------------------------------------------------------------
!        Note:NWPSAF_EMSOW_COEF.dat is used to compute ocean water emissivities
!             NWPSAF_EMSFW_COEF.dat is used to compute fresh water emissivities

WRITE(*,'(/)')
WRITE(*,'(A55)')'READING REGRESSION COEFFICIENTS TO COMPUTE EMISSIVITIES'
!OPEN(IUWAOPC,FILE='NWPSAF_EMSOW_COEF.dat', &
!FORM='UNFORMATTED')
Filename = Trim(Coeffs_Dir) // 'NWPSAF_EMSOW_COEF.dat'
Access='SEQUENTIAL'
Form='UNFORMATTED'
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUWAOPC)               ! out

IF (Precision == 8) THEN
   DO I=1,JPCH
      READ(IUWAOPC)(EMSCOEF(J,I),J=1,9)
   ENDDO
ELSE
   DO I=1,JPCH
      READ(IUWAOPC)(EMSCOEF_IN(J),J=1,9)
      EMSCOEF(:,I)=EMSCOEF_IN(:)
   ENDDO
END IF

WRITE(*,'(A28)')'REGRESSION COEFFICIENTS READ'

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        7.1   This is imported in from SETPRO.
!              ADC 22/12/99
!-----------------------------------------------------------------------
WRITE(*,'(/)')
WRITE(*,'(A20)')'READING EMISSIVITIES'
!OPEN(IUEMS,FILE='IASIEMS0.dat')

Filename = Trim(Coeffs_Dir) // 'IASIEMS0.dat'
Form   = "FORMATTED"
CALL NWPSAF_OpenFile( Trim(filename), & ! in
     Access,              & ! in
     Action,              & ! in
     Status,              & ! in
     Form,                & ! inout
     IUEMS)                 ! out

KEMS=KNPF*JPCH

    DO I=1,KEMS
      READ(IUEMS,'(I4,1X,F10.8)')ICHAN,EMSIR(I)
    END DO

WRITE(*,'(A17)')'EMISSIVITIES READ'


!-----------------------------------------------------------------------
!        8.   DEALLOCATE ARRAYS AND CLOSE ALL INPUT FILES.
!-----------------------------------------------------------------------

IF (Precision /= 4) THEN 
   DEALLOCATE( CFMPOS_IN  ) 
   DEALLOCATE( CFWPOS_IN  ) 
   DEALLOCATE( CFWPOS1_IN ) 
   DEALLOCATE( CFWPOS2_IN ) 
   DEALLOCATE( CFO_IN     ) 
   DEALLOCATE( TREF_IN  ) 
   DEALLOCATE( TREFO_IN ) 
   DEALLOCATE( WREF_IN  ) 
   DEALLOCATE( WREFO_IN ) 
   DEALLOCATE( OREF_IN  ) 
END IF

CLOSE(IUMIXPOS)
!CLOSE(IUMIXNEG)
CLOSE(IUH2O)
CLOSE(IUH2O1)
CLOSE(IUH2O2)
CLOSE(IUOZ)
CLOSE(IUGAM)
CLOSE(IUWAOPC)
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_IASCF
!     Set up input profile to RTIASI
SUBROUTINE RTIASI_SETPRO(PROF,KLENPF,KNPF)

!     Description:
!     Read the input profile to RTIASI. Also check that input values
!     are resonable. If unreasonable value are found a warning message
!     is written to an output file. Read user supplied emissivities.
!
!     Method:
!     Simple
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            12/11/1999  Original code. Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_PRFCON,ONLY : &
!     Imported parameters:
 JPLEV     ! Number of pressure levels

IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF               ! Number of profiles processed.
INTEGER,INTENT(IN) :: KLENPF             ! Dimension of input profile to
REAL,INTENT(IN)    :: PROF (KLENPF,KNPF) ! Input profile to RTIASI
                                         ! RTIASI.

!     End of subroutine arguments


!       Local parameters:

CHARACTER (LEN=50), PARAMETER :: WARNING = &
     '!!!!!!!!!!!!!!!!!!!!!WARNING !!!!!!!!!!!!!!!!!!!!'

CHARACTER (LEN=300), PARAMETER :: ERRSURFT    =            &
     '('' ' // WARNING   // &
     ' '' ,/,'' SURFACE TEMPERATURE VALUE IS OUTSIDE '  // &
     'REASONABLE RANGE'',/,'' ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: ERRSKINT    =            &
     '('' ' // WARNING   // &
     ' '' ,/,'' SKIN    TEMPERATURE VALUE IS OUTSIDE '  // &
     'REASONABLE RANGE'',/,'' ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: ERRSURFP    =            &
     '('' ' // WARNING   // &
     ' '' ,/,'' SURFACE PRESSURE    VALUE IS OUTSIDE '  // &
     'REASONABLE RANGE'',/,'' ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: ERRSURFH    =            &
     '('' ' // WARNING   // &
     ' '' ,/,'' SURFACE HUMIDITY    VALUE IS OUTSIDE '  // &
     'REASONABLE RANGE'',/,'' ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: ERRSURFW    =            &
     '('' ' // WARNING   // &
     ' '' ,/,'' SURFACE WIND SPEED  VALUE IS OUTSIDE '  // &
     'REASONABLE RANGE'',/,'' ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MINTEMP     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' TEMPERATURE IS OUTSIDE THE RANGE OF VALUES '      // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MINIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MAXTEMP     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' TEMPERATURE IS OUTSIDE THE RANGE OF VALUES '      // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MAXIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MINWVAP     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' WATER VAPOUR IS OUTSIDE THE RANGE OF VALUES '     // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MINIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MAXWVAP     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' WATER VAPOUR IS OUTSIDE THE RANGE OF VALUES '     // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MAXIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MINOZON     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' OZONE IS OUTSIDE THE RANGE OF VALUES '            // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MINIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=300), PARAMETER :: MAXOZON     =                &
     '('' ' // WARNING // ' '' ,/,''LEVEL '',I3,'           // &
     ' '' OZONE IS OUTSIDE THE RANGE OF VALUES '            // &
     ' USED TO COMPUTE THE REGRESSION COEFFICIENTS '',/,'   // &
     ' ''MAXIMUM  VALUE= '',F15.5,/,''ACTUAL VALUE = '',F15.5)'

CHARACTER (LEN=3)   :: CPROF

REAL, PARAMETER :: Tol=1.001

!       Local scalars:
INTEGER :: I,KP,JJ,KC

!       Local arrays:
REAL    :: MINW (JPLEV)        ! Minimum value of humidity used to
                               ! compute the regression coefficients.
REAL    :: MAXW (JPLEV)        ! Maximum value of humidity used to
                               ! compute the regression coefficients.
REAL    :: MINT (JPLEV)        ! Minimum value of temperature used to
                               ! compute the regression coefficients.
REAL    :: MAXT (JPLEV)        ! Maximum value of temperature used to
                               ! compute the regression coefficients.
REAL    :: MINO (JPLEV)        ! Minimum value of ozone used to
                               ! compute the regression coefficients.
REAL    :: MAXO (JPLEV)        ! Maximum value of ozone used to
                               ! compute the regression coefficients.
!-----End of header-----------------------------------------------------


DATA MINW/4.361205,5.409970,5.635613,5.429358,4.849816,4.623666, &
4.274895,3.963845,3.741639,3.690460,3.549987,2.994473,2.739244, &
1.916146,2.411677,1.643729,1.651366,2.598531,2.077630,1.727889, &
2.220125,3.201714,3.211240,3.177834,1.446612,0.3097142,2.192869, &
2.783983,5.390884,8.050692,8.050692,21.16817,41.51120,100.7937, &
127.084,113.4186,112.6335,142.723,134.5775,135.7767,148.3392, &
76.34663,42.33958/

DATA MAXW/5.867325,6.233867,6.175500,6.046309,5.741012,5.703165, &
5.844031,5.830869,5.369009,5.397994,5.592806,4.852201,4.593871, &
4.456130,4.333647,4.025418,3.889585,47.88945,115.3653,220.1605, &
373.2881,594.2026,1024.177,1496.372,2253.174,3406.157,4853.446, &
6701.576,8711.748,10387.42,11990.57,14121.00,16615.18,18891.13, &
21266.94,23771.49,26735.84,28639.89,32190.38,36126.02,37362.10, &
37787.14,38028.82/

DATA MINT/190.9450,215.5050,212.8999,198.4644,198.4496,204.5798, &
189.3978,186.9405,185.8896,186.1667,183.1527,183.6486,184.5035, &
183.7715,186.2423,185.1949,186.2353,189.8716,189.5201,189.3538, &
191.6970,196.3361,201.8183,204.3071,207.7184,208.7572,211.7302, &
216.1444,220.3337,222.7711,226.1762,230.0194,233.5735,235.5786, &
236.9657,237.2966,235.5742,234.810,234.7339,235.1377,236.1328, &
233.8223,232.2004/

DATA MAXT/288.1028,296.6161,310.956,309.3185,300.2095,289.6117, &
290.7231,285.7654,264.9886,261.3101,257.8293,253.1519,251.5889, &
254.8312,250.4318,243.4454,242.3198,241.5646,240.8297,242.2247, &
243.8277,245.8346,253.0286,253.6444,259.7492,266.1804,271.1689, &
273.4873,274.3341,275.3654,279.0348,282.4702,284.9223,289.6340, &
294.6081,299.0551,303.4185,306.0771,306.4893,308.8958,310.5137, &
311.5381,311.9639/

DATA MINO/0.6800000,1.249224,1.409621,1.413208,1.414658,1.415650, &
1.415996,1.416218,1.416394,1.416506,1.360321,0.5584311,0.3678274, &
0.1907501,0.1429634,0.1146927,4.6684265E-02,7.8972578E-03, &
6.9316402E-03,5.8482438E-03,4.9787536E-03,4.7907233E-03, &
8.0226958E-03,7.8008175E-03,1.0673031E-02,5.3272247E-03, &
7.5416565E-03,4.4134855E-03,4.1854382E-03,3.3203363E-03, &
1.8514395E-03,1.6137362E-03,1.5238980E-03,1.3898611E-03, &
4.8571825E-04,4.8571825E-04,4.8571825E-04,0.0000000E+00, &
0.0000000E+0,0.0000000E+00,0.0000000E+00,0.0000000E+00, &
0.0000000E+00/

DATA MAXO/8.950000,9.244218,9.340413,9.371043,9.385607,9.392838, &
9.396525,9.403474,9.433378,8.851357,8.767143,6.276594,6.064847, &
5.385192,3.544674,2.912089,2.250107,1.991318,1.680061,1.226486, &
0.8982038,0.8089733,0.5956802,0.4170179,0.3237114,0.2100811, &
0.1716347,0.1342554,0.1286316,0.1264586,0.1167932,0.1103659, &
0.1134226,0.1211433,0.1026807,0.1048491,9.9005103E-02, &
9.4605148E-02,9.2247009E-02,9.0234280E-02,8.9395046E-02, &
8.8834248E-02,8.8595532E-02/

!DATA IUEMS/21/

MAXT=MAXT*Tol
MAXW=MAXW*Tol
MAXO=MAXO*Tol
MINT=MINT/Tol
MINW=MINW/Tol
MINO=MINO/Tol

DO KP=1,KNPF

  WRITE(CPROF,'(I3)')KP

    DO JJ=1,3
      IF(CPROF(JJ:JJ).EQ.' ')THEN
        KC=JJ
      END IF
    END DO

  WRITE(*,'(/)')
  WRITE(*,'(A22,I3)')'READING INPUT PROFILE ',KP



!-----------------------------------------------------------------------
!          1.    STORE INPUT PROFILE VARIABLES
!-----------------------------------------------------------------------


!      PROF(1:43,KP)      ! Temperature
!      PROF(44:86,KP)     ! Humidity
!      PROF(89:129,KP)    ! Ozone
!      PROF(130,KP)       ! Total column cloud liquid water (Not used)
!      PROF(131,KP)       ! Surface air temperature
!      PROF(132,KP)       ! Surface air humidity
!      PROF(133,KP)       ! Surface pressure
!      PROF(134,KP)       ! U Surface wind
!      PROF(135,KP)       ! V Surface wind
!      PROF(136,KP)       ! Skin temperature
!      PROF(137,KP)       ! Cloud top pressure
!      PROF(138,KP)       ! Effective cloud coverage


!-------Input file------------------------------------------------------
!  OPEN(10,FILE='PROFIN'//CPROF(KC+1:3)//'.dat')
!-----------------------------------------------------------------------


!-------Output files----------------------------------------------------
!  OPEN(11,FILE='WARNMESS'//CPROF(KC+1:3)//'.out')
!-----------------------------------------------------------------------



    DO I=1,43
!      READ(10,*)PRES(I),PROF(I,KP),PROF(I+43,KP),PROF(I+86,KP)


!-----Check if actual pressure levels match standard levels-------------

!        IF(ABS(PRES(I)-XPRES(I)) > 0.1 ) THEN
!          WRITE(*,'(/)')
!          WRITE(*,FMT=PRESERR)I,XPRES(I),PRES(I)
!          WRITE(*,'(/)') &
!             'ERROR IN PRESSURE LEVEL:SEE FILE WARNMESS'// &
!             CPROF(KC+1:3)//'.out'
!          STOP
!        END IF

!-----Check if temperature profile is out of bounds---------------------

        IF(PROF(I,KP) < MINT(I))THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=MINTEMP)I,MINT(I),PROF(I,KP)
          WRITE(*,'(/)')
        ELSE IF(PROF(I,KP) > MAXT(I))THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=MAXTEMP)I,MAXT(I),PROF(I,KP)
          WRITE(*,'(/)')
        END IF

!-----Check if humidity profile is out of bounds------------------------

        IF(PROF(I+43,KP) < MINW(I))THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=MINWVAP)I,MINW(I),PROF(I+43,KP)
          WRITE(*,'(/)')
        ELSE IF(PROF(I+43,KP) > MAXW(I))THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=MAXWVAP)I,MAXW(I),PROF(I+43,KP)
          WRITE(*,'(/)')
        END IF

!-----Check if ozone profile is out of bounds---------------------------

        IF(PROF(I+86,KP) < MINO(I))THEN
          WRITE(*,'(/)')
          WRITE(*,MINOZON)I,MINO(I),PROF(I+86,KP)
          WRITE(*,'(/)')
        ELSE IF (PROF(I+86,KP) > MAXO(I))THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=MAXOZON)I,MAXO(I),PROF(I+86,KP)
          WRITE(*,'(/)')
        END IF
     ENDDO

     !          DO I=1,9
     !            READ(10,*)PROF(129+I,KP)
     !          END DO
     
     !-----Check if surface temperature is reasonable------------------------

     IF (PROF(131,KP) > 400. )THEN
        WRITE(*,'(/)')
        WRITE(*,FMT=ERRSURFT)PROF(131,KP)
        !WRITE(*,'(/)') &
        !     'SURFACE TEMPERATURE TOO LOW:SEE FILE WARNMESS'// &
        !     CPROF(KC+1:3)//'.out'
        STOP
     ELSE IF (PROF(131,KP) < 150.) THEN
        WRITE(*,'(/)')
        WRITE(*,FMT=ERRSURFT)PROF(131,KP)
        !          WRITE(*,'(/)') &
        !               'SURFACE TEMPERATURE TOO HIGH:SEE FILE WARNMESS'// & 
        !               CPROF(KC+1:3)//'.out'
        STOP
     END IF
     
!-----Check if skin temperature is reasonable---------------------------

        IF  (PROF(136,KP) > 400. ) THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=ERRSKINT)PROF(136,KP)
!          WRITE(*,'(/)') &
!               'SKIN TEMPERATURE TOO LOW:SEE FILE WARNMESS'// &
!               CPROF(KC+1:3)//'.out'
          STOP
        ELSE IF (PROF(136,KP) < 150.) THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=ERRSKINT)PROF(136,KP)
!          WRITE(*,'(/)') &
!               'SKIN TEMPERATURE TOO HIGH:SEE FILE WARNMESS'// &
!               CPROF(KC+1:3)//'.out'
          STOP
        END IF

!-----Check if surface pressure is reasonable---------------------------

        IF  (PROF(133,KP) < 400.) THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=ERRSURFP)PROF(133,KP)
!          WRITE(*,'(/)') &
!               'SURFACE PRESSURE TOO LOW:SEE FILE WARNMESS'// &
!               CPROF(KC+1:3)//'.out'
          STOP
          ELSE IF ( PROF(133,KP) > 1100.) THEN
          WRITE(*,'(/)')
          WRITE(*,FMT=ERRSURFP)PROF(133,KP)
!          WRITE(*,'(/)') &
!               'SURFACE PRESSURE TOO HIGH:SEE FILE WARNMESS'// &
!               CPROF(KC+1:3)//'.out'
          STOP
          END IF

!-----Check if surface humidity is reasonable---------------------------

         IF  (PROF(132,KP) < 0. ) THEN
           WRITE(*,'(/)')
           WRITE(*,FMT=ERRSURFH)PROF(132,KP)
!           WRITE(*,'(/)') &
!                'SURFACE HUMIDITY TOO LOW:SEE FILE WARNMESS'// &
!                CPROF(KC+1:3)//'.out'
           STOP
         ELSE IF ( PROF(132,KP) > 85000.) THEN
           WRITE(*,'(/)')
           WRITE(*,FMT=ERRSURFH)PROF(132,KP)
!           WRITE(*,'(/)') &
!                'SURFACE HUMIDITY TOO HIGH:SEE FILE WARNMESS'// &
!                CPROF(KC+1:3)//'.out'
           STOP
         END IF

!-----Check if surface wind speed is reasonable-------------------------

         IF (SQRT(PROF(134,KP)**2+PROF(135,KP)**2) > 100. ) THEN
           WRITE(*,'(/)')
           WRITE(*,FMT=ERRSURFW)SQRT(PROF(134,KP)**2 + &
                PROF(135,KP)**2)
!           WRITE(*,'(/)') &
!                'SURFACE WIND SPEED TOO HIGH:SEE FILE WARNMESS'// &
!                CPROF(KC+1:3)//'.out'
           STOP
         END IF

  WRITE(*,'(A18)')'INPUT PROFILE READ'

END DO              ! End profile loop


!-----------------------------------------------------------------------
!        2.   STORE EMISSIVITIES FOR IASI CHANNELS
!-----------------------------------------------------------------------

! This has been moced to IASCF.  ADC 22/12/99

!-----------------------------------------------------------------------


RETURN
END SUBROUTINE RTIASI_SETPRO
!     Fast radiative transfer model for IASI.
SUBROUTINE RTIASI(PROF,KNPF,KLENPF,PANGL,PANGS,KSURF, &
KCHAN,KPROF,KNCHPF,PRAD,PTB,KSAT, &
TAUM,TAUW,TAUO,TAU,TAUSFC,EMS,IPRINT)

!     Description:
!     Compute multi-channel IASI radiances and brightness
!     temperatures for many profiles.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     3) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            23/08/1990  Original Code (RTTOV). J.R.Eyre. ECMWF.
!     2            12/11/1999  Marco Matricardi. ECMWF.
!                              i)  New fast transmittance scheme to predict IASI
!                                  radiances.
!                              ii) Fast emissivity model to compute sea surface
!                                  infrared emissivities.
!                              iii)New definitionm of average layer temperature
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPLEV  ,  & ! Number of pressure levels
  JPPF   ,  & ! Max. number of profiles to be processed
  JPCOFM ,  & ! Max. number of mixed gas coefficients
  JPCOFW ,  & ! Max. number of water vapour coefficients
  JPCOFO    ! Max. number of ozone coefficients

IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN):: IPRINT
INTEGER, INTENT(IN):: KNPF           ! Number of processed profiles
INTEGER, INTENT(IN):: KNCHPF         ! Number of processed radiances
INTEGER, INTENT(IN):: KLENPF         ! Lenght of input profile to RTIASI
INTEGER, INTENT(IN):: KSAT           ! Input satellite index

!       Array arguments with intent in:
INTEGER, INTENT(IN):: KSURF (KNPF)   ! Array of surface types indices
INTEGER, INTENT(IN):: KCHAN (KNCHPF) ! Array of channel indices
INTEGER, INTENT(IN):: KPROF (KNCHPF) ! Array of profile indices
REAL,INTENT(IN)    :: PROF           & ! Input profile to RTIASI
(KLENPF,KNPF)

!       Array arguments with intent out:
REAL,INTENT(OUT)   :: PRAD   (KNCHPF)       ! Computed radiance
REAL,INTENT(OUT)   :: PTB    (KNCHPF)       ! Computed br. temperature
REAL,INTENT(OUT)   :: TAU    (KNCHPF,JPLEV) ! Tau from each level to space
REAL,INTENT(OUT)   :: TAUSFC (KNCHPF)       ! Tau from surface to space
REAL,INTENT(OUT)   :: TAUM   (KNCHPF,JPLEV) ! Fixed gases transmittances
REAL,INTENT(OUT)   :: TAUW   (KNCHPF,JPLEV) ! Water vapour  transmittances
REAL,INTENT(OUT)   :: TAUO   (KNCHPF,JPLEV) ! Ozone transmittances
REAL,INTENT(OUT)   :: EMS    (KNCHPF)       ! Emissivity for IASI channels


!     End of subroutine arguments


!       Local arrays:
REAL :: XXM  (JPCOFM,JPLEV,JPPF)  ! Functions of prof. for mix. gas calc.
REAL :: XXW  (JPCOFW,JPLEV,JPPF)  ! Functions of prof. for wat. va. calc.
REAL :: XXO  (JPCOFO,JPLEV,JPPF)  ! Functions of prof. for ozone    calc.
REAL :: XXW1 (JPCOFW,JPLEV,JPPF)  ! Functions of prof.: level 1, wa. vap.
REAL :: XXO1 (JPCOFO,JPLEV,JPPF)  ! Functions of prof.: level 1, ozone
REAL :: OPDPMP  (KNCHPF,JPLEV)    ! Mixed gases optical depth.
REAL :: OPDPO   (KNCHPF,JPLEV)    ! Ozone optical depth
REAL :: OPDPWP1 (KNCHPF,JPLEV)    ! Wa. vap. op. depth :weak absorption
REAL :: OPDPWP2 (KNCHPF,JPLEV)    ! Wa. vap. op. depth :strong absorption
REAL :: OPDPWP  (KNCHPF,JPLEV)    ! Wa. vap. op. depth :single  scheme
REAL :: B       (KNCHPF,JPLEV)    ! Planck functions for temp. profiles
REAL :: BCG     (KNCHPF,JPLEV)    ! Curtis-Godson layer temperature
REAL :: BA      (KNCHPF)          ! Planck funtions for surface air temp.
REAL :: BS      (KNCHPF)          ! Planck funtions for surface skin temp.
REAL :: RADOV   (KNCHPF,JPLEV)    ! Overcast radiances for cloud at each
                                  ! standard pressure level.
REAL :: RADO    (KNCHPF)          ! Overcast radiances for given cloud-top
                                  ! pressures.
REAL :: BDT     (KNCHPF,JPLEV)    ! Stores upwelling radiation from
                                  ! atmosphere above each level.
REAL :: BDTR    (KNCHPF,JPLEV)    ! Stores downwelling radiation at each
                                  ! level from atmosphere above.
REAL :: FCLD    (KNCHPF)          ! Effective fractional cloud cover for
                                  ! each chan * prof.
REAL :: RADCL   (KNCHPF)          ! Clear column radiances
REAL :: TBCL    (KNCHPF)          ! Clear column brightness temp. in K
REAL :: PANGL   (KNPF)            ! Satellite sounding Nadir angle
REAL :: PANGS   (KNPF)            ! Solar zenith angle.
REAL :: PWWR    (JPLEV,JPPF)      ! Pressure weighted water vapour ratio


!-----End of header-----------------------------------------------------


!-----Note: all radiances are in mW/m**2/sr/cm**-1----------------------


IF( IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(59A)') &
'SETTING UP PROFILE VARIABLES FOR TRANSMITTANCE CALCULATIONS'
ENDIF

!-----------------------------------------------------------------------
!       1.  SET UP PROFILE VARIABLES.
!-----------------------------------------------------------------------
CALL RTIASI_PRFIN(PROF,KNPF,KLENPF,PANGL,PANGS,KSURF,KSAT, &
                 XXM,XXW,XXW1,XXO,XXO1,PWWR,IPRINT)

IF( IPRINT == 1)THEN
WRITE(*,'(24A)')'PROFILE VARIABLES SET UP'
WRITE(*,'(/)')
WRITE(*,'(26A)')'CALCULATING OPTICAL DEPTHS'
ENDIF

!-----------------------------------------------------------------------
!       2.  CALCULATE LAYER OPTICAL DEPTHS.
!-----------------------------------------------------------------------
CALL RTIASI_OPDEP(KCHAN,KPROF,KNCHPF,KSAT,XXM,XXW,XXW1,XXO, &
XXO1,OPDPMP,OPDPO,OPDPWP1,OPDPWP2,OPDPWP)

IF( IPRINT == 1)THEN
WRITE(*,'(25A)')'OPTICAL DEPTHS CALCULATED'

WRITE(*,'(/)')
WRITE(*,'(26A)')'CALCULATING TRANSMITTANCES'
ENDIF

!-----------------------------------------------------------------------
!       3.  CALCULATE TRANSMITTANCES ON RADIATIVE TRANSFER MODEL LEVELS.
!-----------------------------------------------------------------------
CALL RTIASI_RTTAU(KCHAN,KPROF,KNCHPF,KSAT,TAU,TAUSFC,PWWR,OPDPMP, &
OPDPO,OPDPWP1,OPDPWP2,OPDPWP,TAUM,TAUW,TAUO &
)

IF( IPRINT == 1)THEN
WRITE(*,'(25A)')'TRANSMITTANCES CALCULATED'

WRITE(*,'(/)')
WRITE(*,'(23A)')'INTEGRATING RT EQUATION'
WRITE(*,'(/)')
ENDIF

!-----------------------------------------------------------------------
!       4.  INTEGRATE RADIATIVE TRANSFER EQUATION.
!-----------------------------------------------------------------------
CALL RTIASI_RTINT(KCHAN,KPROF,KNCHPF,PRAD,PTB, &
EMS,TAU,TAUSFC,B,BCG,BA,BS,RADOV,RADO,BDT,BDTR,FCLD,RADCL,TBCL, &
PANGL,KNPF,KSAT,IPRINT)

IF( IPRINT == 1)THEN
WRITE(*,'(38A)')'RADIATIVE TRANSFER EQUATION INTEGRATED'
ENDIF


RETURN
END SUBROUTINE RTIASI
!     Fast radiative transfer model for IASI.
SUBROUTINE RTIASIK(PROF,PROF_D,KNPF,KLENPF,PANGL,PANGS,KSURF, &
KCHAN,KPROF,KNCHPF,KSAT,PTB_D,PRAD_D,IPRINT)

!     Description:
!     K of subroutine RTIASI
!     Compute multi-channel IASI radiances and brightness
!     temperatures for many profiles.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2)See: User's manual for RTIASI (Available from EUMETSAT)
!     3)See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            23/08/1990  Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999  Marco Matricardi. ECMWF.
!                              i)  New fast transmittance scheme to predict IASI
!                                  radiances.
!                              ii) Fast emissivity model to compute sea surface
!                                  infrared emissivities.
!                              iii)New definitionm of average layer temperature
!     2.1          04/01/2000  TB_D and PRAD_D defined via module.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct module used:

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPLEV             ,  & ! Number of pressure levels
  JPPF              ,  & ! Max. number of profiles to be processed
  JPCOFM            ,  & ! Max. number of mixed gas coefficients
  JPCOFW            ,  & ! Max. number of water vapour coefficients
  JPCOFO            ,  & ! Max. number of ozone coefficients
  JPCHPF               ! Max. number of processed radiances

!     K module used:

USE RTIASI_IRCLDK, ONLY : &
!     Imported arrays:
  FRACPC               ! Fraction of standardd pressure level interval
                       ! by which cloud is above level NLEVCD.

USE RTIASI_EMISIRK,ONLY : &
!     Imported tangent linear arrays:
  EMSIR                ! Matrix of partial derivatives with respect
                       ! to the emissivities

USE RTIASI_SURFK, ONLY : &
!     Imported arrays:
  FRACPS        ! Fraction of standard pressure level interval by
                ! which surf is above level nlevsf.

USE RTIASI_PRFVARK,ONLY : &
!     Imported arrays:
  TEMPW  ,  & ! Curtis-Godson Temperature profile in K
  TEMP   ,  & ! Temperature profile in K
  WMIX   ,  & ! Water vapour volume mixing ratio in ppmv
  OMIX   ,  & ! Ozone volume mixing ratio in ppmv
  TA     ,  & ! Surface air temperature in K
  WMIXS  ,  & ! Water vapour surface volume mixing ratio in ppmv
  TS     ,  & ! Skin temperature IN K
  SURFP  ,  & ! Surface pressure in mb=hPa
  CLDP   ,  & ! Cloud-top pressure in mb=hPa
  CLDF   ,  & ! Fractional (ir) cloud cover
  CLW    ,  & ! Total column cloud liquid water in mm
  SURFWU ,  & ! U component of surface wind in m/sec
  SURFWV    ! V component of surface wind in m/sec

IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN):: IPRINT
INTEGER, INTENT(IN):: KNPF                ! Number of processed profiles.
INTEGER, INTENT(IN):: KNCHPF              ! Number of processed radiances.
INTEGER, INTENT(IN):: KLENPF              ! Lenght of input profile
                                          ! to RTIASI.
INTEGER, INTENT(IN):: KSAT                ! Input satellite index.

!       Array arguments with intent in:
INTEGER, INTENT(IN):: KSURF (KNPF)        ! Array of surface types indices
INTEGER, INTENT(IN):: KCHAN (KNCHPF)      ! Array of channel indices.
INTEGER, INTENT(IN):: KPROF (KNCHPF)      ! Array of profile indices.

!       Direct array arguments with intent out:
REAL,INTENT(IN)    :: PROF_D (KLENPF,KNPF)! Input profile to RTIASI.

!       K array arguments with intent out:
REAL,INTENT(OUT)   :: PROF (KLENPF,KNCHPF)! Matrix of partial derivatives
                                          ! with respect to the input
                                          ! profile


!     End of subroutine arguments




!       Direct Local arrays:
REAL :: XXM_D     (JPCOFM,JPLEV,JPPF) ! Func. of prof for mix. gas calc.
REAL :: XXW_D     (JPCOFW,JPLEV,JPPF) ! Func. of prof. for wa. vap. calc.
REAL :: XXO_D     (JPCOFO,JPLEV,JPPF) ! Func. of profile for ozone calc.
REAL :: XXW1_D    (JPCOFW,JPLEV,JPPF) ! Func. of profile: level 1,wa. vap.
REAL :: XXO1_D    (JPCOFO,JPLEV,JPPF) ! Func. of profile: level 1, ozone
REAL :: OPDPMP_D  (KNCHPF,JPLEV)      ! Mixed gases optical depth.
REAL :: OPDPO_D   (KNCHPF,JPLEV)      ! Ozone optical depth
REAL :: OPDPWP1_D (KNCHPF,JPLEV)      ! Wa. vap. op. depth :weak absorpt.
REAL :: OPDPWP2_D (KNCHPF,JPLEV)      ! Wa. vap. op. depth :strong absorpt
REAL :: OPDPWP_D  (KNCHPF,JPLEV)      ! Wa. vap. op. depth:single scheme
REAL :: TAU_D     (KNCHPF,JPLEV)      ! Trans. from each level to space
REAL :: TAUSFC_D  (KNCHPF)            ! Trans. from surface to space
REAL :: TAUM_D    (KNCHPF,JPLEV)      ! Fixed gases transmittances
REAL :: TAUW_D    (KNCHPF,JPLEV)      ! Water vapour transmittances
REAL :: TAUO_D    (KNCHPF,JPLEV)      ! Ozone transmittances
REAL :: B_D       (KNCHPF,JPLEV)      ! Planck functions for temp. prof.
REAL :: BCG_D     (KNCHPF,JPLEV)
REAL :: BA_D      (KNCHPF)            ! Planck func. for surface air temp.
REAL :: BS_D      (KNCHPF)            ! Planck func. for surf. skin temp.
REAL :: RADOV_D   (KNCHPF,JPLEV)      ! Overcast radiances for cloud at
                                      ! each standard pressure level.
REAL :: RADO_D    (KNCHPF)            ! Overcast radiances for given
                                      ! cloud-top pressures.
REAL :: BDT_D     (KNCHPF,JPLEV)      ! Stores upwelling radiation from
                                      ! atmosphere above each level.
REAL :: BDTR_D    (KNCHPF,JPLEV)      ! Stores downwelling radiation at
                                      ! each level from atmosphere above.
REAL :: FCLD_D    (KNCHPF)            ! Effective fractional cloud cover
                                      ! for each chan * prof.
REAL :: RADCL_D   (KNCHPF)            ! Clear column radiances
REAL :: TBCL_D    (KNCHPF)            ! Clear column brightness temp.
REAL :: EMS_D     (KNCHPF)            ! Emissivity for IASI channels
REAL :: PANGL     (KNPF)              ! Satellite sounding Nadir angle
REAL :: PANGS     (KNPF)              ! Solar zenith angle.
REAL :: PWWR_D    (JPLEV,JPPF)
REAL :: PRAD_D    (KNCHPF)            ! Computed radiance
REAL :: PTB_D     (KNCHPF)            ! Computed brigh. temperature

!       K local arrays:
REAL :: PRAD      (KNCHPF)            ! Computed radiance
REAL :: PTB       (KNCHPF)            ! Computed brigh. temperature
REAL :: EMS       (KNCHPF)            ! Emissivity for IASI channel
REAL :: XXM       (JPCOFM,JPLEV,JPCHPF) ! Func. of prof for mix. gas calc
REAL :: XXW       (JPCOFW,JPLEV,JPCHPF) ! Func. of prof for wa. vap. calc
REAL :: XXO       (JPCOFO,JPLEV,JPCHPF) ! Func. of prof for ozone calc.
REAL :: XXW1      (JPCOFW,JPLEV,JPCHPF) ! Func. of prof: level 1, wa. vap
REAL :: XXO1      (JPCOFO,JPLEV,JPCHPF) ! Func. of prof: level 1, ozone
REAL :: OPDPMP    (KNCHPF,JPLEV)      ! Mixed gases optical depth.
REAL :: OPDPO     (KNCHPF,JPLEV)      ! Ozone optical depth
REAL :: OPDPWP1   (KNCHPF,JPLEV)      ! Wa. vap. op. depth :strong absorpt
REAL :: OPDPWP2   (KNCHPF,JPLEV)      ! Wa. vap. op. depth :strong absorpt
REAL :: OPDPWP    (KNCHPF,JPLEV)      ! Wa. vap. op. depth.:single scheme
REAL :: TAU       (KNCHPF,JPLEV)      ! Trans. from each level to space
REAL :: TAUSFC    (KNCHPF)            ! Trans. from surface to space
REAL :: TAUM      (KNCHPF,JPLEV)
REAL :: TAUW      (KNCHPF,JPLEV)
REAL :: TAUO      (KNCHPF,JPLEV)
REAL :: B         (KNCHPF,JPLEV)      ! Planck functions for temp. prof.
REAL :: BCG       (KNCHPF,JPLEV)
REAL :: BA        (KNCHPF)            ! Planck func. for surf. air temp.
REAL :: BS        (KNCHPF)            ! Planck func. for surf. skin temp.
REAL :: RADOV     (KNCHPF,JPLEV)      ! Overcast radiances for cloud at
                                      ! each standard pressure level.
REAL :: RADO      (KNCHPF)            ! Overcast radiances for given
                                      ! cloud-top pressures.
REAL :: BDT       (KNCHPF,JPLEV)      ! Stores upwelling radiation from
                                      ! atmosphere above each level.
REAL :: BDTR      (KNCHPF,JPLEV)      ! Stores downwelling radiation at
                                      ! each level from atmosphere above.
REAL :: FCLD      (KNCHPF)            ! Effective fractional cloud cover
                                      ! for each chan * prof.
REAL :: RADCL     (KNCHPF)            ! Clear column radiances
REAL :: PWWR      (JPLEV,JPCHPF)


!-----End of header-----------------------------------------------------


!-----All radiances in mW/m**2/sr/cm**-1--------------------------------


!-----------------------------------------------------------------------
!        0.   REPEAT DIRECT CALCULATION
!-----------------------------------------------------------------------

IF(IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(43A)')'K  OF RTIASI:REPEATING DIRECT CALCULATIONS'

WRITE(*,'(/)')
WRITE(*,'(59A)') &
'SETTING UP PROFILE VARIABLES FOR TRANSMITTANCE CALCULATIONS'
ENDIF


CALL RTIASI_PRFIN(PROF_D,KNPF,KLENPF,PANGL,PANGS,KSURF,KSAT, &
                 XXM_D,XXW_D,XXW1_D,XXO_D,XXO1_D,PWWR_D,IPRINT)

IF(IPRINT == 1)THEN
WRITE(*,'(24A)')'PROFILE VARIABLES SET UP'

WRITE(*,'(/)')
WRITE(*,'(26A)')'CALCULATING OPTICAL DEPTHS'
ENDIF


CALL RTIASI_OPDEP(KCHAN,KPROF,KNCHPF,KSAT,XXM_D,XXW_D,XXW1_D,XXO_D, &
XXO1_D,OPDPMP_D,OPDPO_D,OPDPWP1_D,OPDPWP2_D,OPDPWP_D)

IF(IPRINT == 1)THEN
WRITE(*,'(25A)')'OPTICAL DEPTHS CALCULATED'

WRITE(*,'(/)')
WRITE(*,'(26A)')'CALCULATING TRANSMITTANCES'
ENDIF

CALL RTIASI_RTTAU(KCHAN,KPROF,KNCHPF,KSAT,TAU_D,TAUSFC_D,PWWR_D, &
OPDPMP_D,OPDPO_D,OPDPWP1_D,OPDPWP2_D,OPDPWP_D, &
TAUM_D,TAUW_D,TAUO_D)

IF(IPRINT == 1)THEN
WRITE(*,'(25A)')'TRANSMITTANCES CALCULATED'

WRITE(*,'(/)')
WRITE(*,'(23A)')'INTEGRATING RT EQUATION'
WRITE(*,'(/)')
ENDIF

CALL RTIASI_RTINT(KCHAN,KPROF,KNCHPF,PRAD_D,PTB_D, &
EMS_D,TAU_D,TAUSFC_D,B_D,BCG_D,BA_D,BS_D,RADOV_D,RADO_D,BDT_D, &
BDTR_D,FCLD_D,RADCL_D,TBCL_D,PANGL,KNPF,KSAT,IPRINT)

IF(IPRINT == 1)THEN
WRITE(*,'(38A)')'RADIATIVE TRANSFER EQUATION INTEGRATED'

WRITE(*,'(/)')
WRITE(*,'(28A)')'DIRECT CALCULATIONS REPEATED'
WRITE(*,'(/)')
ENDIF



!-----------------------------------------------------------------------
!       4.  INTEGRATE RADIATIVE TRANSFER EQUATION.
!-----------------------------------------------------------------------


IF(IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(47A)')'----------K              CALCULATIONS----------'
WRITE(*,'(/)')
WRITE(*,'(29A)')'INTEGRATING K  OF RT EQUATION'
WRITE(*,'(/)')
ENDIF


!-----Initialize variables----------------------------------------------

RADCL(:)     =0
RADO(:)      =0
FCLD(:)      =0
RADOV(:,:)   =0
B(:,:)       =0
BDT(:,:)     =0
BDTR(:,:)    =0
BCG(:,:)     =0
BS(:)        =0
BA(:)        =0


SURFWU(:)    =0
SURFWV(:)    =0
TEMP(:,:)    =0
WMIX(:,:)    =0
OMIX(:,:)    =0
TA(:)        =0
WMIXS(:)     =0
TS(:)        =0
SURFP(:)     =0
CLDP(:)      =0
CLDF(:)      =0
CLW(:)       =0
FRACPC(:)    =0
FRACPS(:)    =0

OPDPO  (:,:) =0
OPDPMP (:,:) =0
OPDPWP1(:,:) =0
OPDPWP2(:,:) =0
OPDPWP (:,:) =0
PWWR   (:,:) =0
XXM(:,:,:)   =0
XXW(:,:,:)   =0
XXO(:,:,:)   =0
XXW1(:,:,:)  =0
XXO1(:,:,:)  =0
EMS(:)       =0
TAU(:,:)     =0
TAUM(:,:)    =0
TAUW(:,:)    =0
TAUO(:,:)    =0
TAUSFC(:)    =0
TEMPW(:,:)   =0
EMSIR(:)     =0


PRAD(:)      =1.
PTB (:)      =1.
!-----------------------------------------------------------------------



!CALL RTINTK(KCHAN,KPROF,KNCHPF,PRAD_D,PRAD,PTB_D,PTB,EMS_D, &
!EMS,TAU_D,TAU,TAUSFC_D,TAUSFC,B_D,B,BCG_D,BCG,BA_D,BA,BS_D,BS, &
!RADOV_D,RADOV,RADO_D,RADO,BDT_D,BDT,BDTR_D,BDTR,FCLD_D,FCLD, &
!RADCL_D,RADCL,TBCL_D,TBCL,PANGL,KNPF,KSAT,IPRINT)
CALL RTIASI_RTINTK(KCHAN,KPROF,KNCHPF,PRAD_D,PRAD,PTB_D,PTB,EMS_D, &
EMS,TAU_D,TAU,TAUSFC_D,TAUSFC,B_D,B,BCG_D,BCG,BA_D,BA,BS_D,BS, &
RADOV_D,RADOV,RADO_D,RADO,BDT,BDTR_D,BDTR,FCLD_D,FCLD, &
RADCL_D,RADCL,PANGL,KNPF,KSAT,IPRINT)

IF(IPRINT == 1)THEN
  WRITE(*,'(44A)')'K OF RADIATIVE TRANSFER EQUATION INTEGRATED'
ENDIF


!-----------------------------------------------------------------------
!       3.  CALCULATE TRANSMITTANCES ON RADIATIVE TRANSFER MODEL LEVELS.
!-----------------------------------------------------------------------

IF(IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(27A)')'K CALCULATIONS'
WRITE(*,'(32A)')'CALCULATING K OF TRANSMITTANCES'
ENDIF


CALL RTIASI_RTTAUK(KCHAN,KPROF,KNCHPF,KSAT,TAU_D,TAU, &
TAUSFC,PWWR_D,OPDPMP_D,OPDPMP,OPDPO_D,OPDPO, &
OPDPWP1_D,OPDPWP1,OPDPWP2_D,OPDPWP2,OPDPWP_D,OPDPWP, &
TAUM_D,TAUM,TAUW_D,TAUW,TAUO_D,TAUO)

IF(IPRINT == 1)THEN
  WRITE(*,'(30A)')'K OF TRANSMITTANCES CALCULATED'
ENDIF

!-----------------------------------------------------------------------
!       2.  CALCULATE LAYER OPTICAL DEPTHS.
!-----------------------------------------------------------------------

IF(IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(27A)')'K CALCULATIONS'
WRITE(*,'(32A)')'CALCULATING K OF OPTICAL DEPTHS'
ENDIF


CALL RTIASI_OPDEPK(KCHAN,KNCHPF,KSAT,XXM,XXW, &
XXW1,XXO,XXO1,OPDPMP,OPDPO, &
OPDPWP1,OPDPWP2,OPDPWP)

IF(IPRINT == 1)THEN
  WRITE(*,'(31A)')'K OF OPTICAL DEPTHS CALCULATED'
ENDIF


!-----------------------------------------------------------------------
!       1.  SET UP PROFILE VARIABLES.
!-----------------------------------------------------------------------

IF(IPRINT == 1)THEN
WRITE(*,'(/)')
WRITE(*,'(27A)')'K CALCULATIONS'
WRITE(*,'(62A)') &
'SETTING UP PROFILE VARIABLES FOR AD TRANSMITTANCE CALCULATIONS'
ENDIF

CALL RTIASI_PRFINK(PROF,KNPF,KNCHPF,KPROF,KLENPF, &
                 XXM_D,XXM,XXW_D,XXW,XXW1,XXO_D,XXO, &
                 XXO1,PWWR_D,PWWR)



IF(IPRINT == 1)THEN
  WRITE(*,'(24A)')'PROFILE VARIABLES SET UP'
ENDIF



RETURN
END SUBROUTINE RTIASIK
!     Set up profile variables for radiative transfer calculations.
SUBROUTINE RTIASI_PRFIN(PROF,KNPF,KLENPF,PANGL,PANGS,KSURF,KSAT, &
                 XXM,XXW,XXW1,XXO,XXO1,PWWR,IPRINT)

!     Description:
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     2) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/1990   Original code(RTTOV). J.R.EYRE. ECMWF
!     2           12/11/1999   Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOCON,ONLY : &
!     Imported arrays:
  RATOE         ! Ratio satellite-orbit/earth radii


USE RTIASI_GEOPTH, ONLY : &
!     Imported arrays:
  XPATH  ,      & ! Secant of viewing path angle at surface
  XPATH1 ,      & ! XPATH - 1
  SQTPTH ,      & ! Sqrt of XPATH
  XPATHS        ! Secant of solar zenith angle at surface


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NLEVSF ,      & ! Index of nearest standard pressure level
                ! at/below surface.
  FRACPS ,      & ! Fraction of standard pressure level interval by
                ! which surf is above level nlevsf.
  NSTYPE        ! Surface type index; 1=sea, 2=land.


USE RTIASI_PRFVAR,ONLY : &
!     Imported parameters:
  JPLEV  ,      & ! Number of pressure levels
  JPPF   ,      & ! Max. number of profiles to be processed
!     Imported arrays:
  TEMP   ,      & ! Temperature profile in K
  WMIX   ,      & ! Water vapour volume mixing ratio in ppmv
  OMIX   ,      & ! Ozone volume mixing ratio in ppmv
  TA     ,      & ! Surface air temperature in K
  WMIXS  ,      & ! Water vapour surface volume mixing ratio in ppmv
  TS     ,      & ! Skin temperature IN K
  SURFP  ,      & ! Surface pressure in mb=hPa
  CLDP   ,      & ! Cloud-top pressure in mb=hPa
  CLDF   ,      & ! Fractional (ir) cloud cover
  CLW    ,      & ! Total column cloud liquid water in mm
  SURFWU ,      & ! U component of surface wind in m/sec
  SURFWV        ! V component of surface wind in m/sec


USE RTIASI_IRCLD, ONLY : &
!     Imported arrays:
  NLEVCD    , & ! Index of nearest standard pressure level
                ! at/below cloud top.
  FRACPC        ! Fraction of standardd pressure level interval
                ! by which cloud is above level NLEVCD.


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCOFM ,      & ! Max. number of mixed gas coefficients
  JPCOFW ,      & ! Max. number of water vapour coefficients
  JPCOFO          ! Max. number of ozone coefficients


IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: IPRINT
INTEGER,INTENT(IN) :: KNPF                       ! Number of processed
                                                 ! profiles.
INTEGER,INTENT(IN) :: KLENPF                     ! Lenght of input profile
                                                 ! to RTIASI.
INTEGER,INTENT(IN) :: KSAT                       ! Input satellite index

!       Array arguments with intent in:
INTEGER,INTENT(IN) :: KSURF (KNPF)               ! Array of surface types
                                                 ! indices.
REAL,INTENT(IN)    :: PANGL (KNPF)               ! Satellite sounding
                                                 ! Nadir angle.
REAL,INTENT(IN)    :: PANGS (KNPF)               ! Solar zenith angle.
REAL,INTENT(IN)    :: PROF  (KLENPF,KNPF)        ! Input profile to RTIASI

!       Array arguments with intent out:
REAL,INTENT(OUT)   :: XXM   (JPCOFM,JPLEV,JPPF)  ! Functions of profile
                                                 ! for fixed gas
                                                 ! calculations.
REAL,INTENT(OUT)   :: XXW   (JPCOFW,JPLEV,JPPF)    ! Functions of profile
                                                 ! for water vapour
                                                 ! calculations.
REAL,INTENT(OUT)   :: XXO   (JPCOFO,JPLEV,JPPF)  ! Functions of profile
                                                 ! for ozone calculations.
REAL,INTENT(OUT)   :: XXW1  (JPCOFW,JPLEV,JPPF)  ! Functions of profile:
                                                 ! level 1, water vapour.
REAL,INTENT(OUT)   :: XXO1  (JPCOFO,JPLEV,JPPF)  ! Functions of profile:
                                                 ! level 1, ozone.
REAL,INTENT(OUT)   :: PWWR  (JPLEV,JPPF)         ! Functions of profile:
                                                 ! pressure weighted water
                                                 ! vapour ratio.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,JL
REAL    :: ZA                     ! Zenith angle
REAL    :: DTR                    ! Degs to rads = PI/180
REAL    :: ZANGMX                 ! Maximum value for solar zenith angle
REAL    :: ZPTHMX


!-----End of header-----------------------------------------------------


DATA DTR/0.0174533/
DATA ZANGMX/89./
DATA ZPTHMX/1.0E6/


!-----------------------------------------------------------------------
!         1.   SETS UP PROFILE-DEPENDENT 'KNOWN' QUANTITIES.
!-----------------------------------------------------------------------



!-----1.1  Satellite and solar geometry.--------------------------------

DO J=1,KNPF

   ZA=SIN(PANGL(J)*DTR)*RATOE(KSAT)
   XPATH(J)=1./SQRT(1.-ZA*ZA)           ! View angle path factor
   SQTPTH(J)=SQRT(XPATH(J))             ! .. and its square root
   XPATH1(J)=XPATH(J)-1.                ! .. and path factor - 1

     IF(PANGS(J) <= ZANGMX) THEN
       XPATHS(J)=1./COS(PANGS(J)*DTR)   ! Solar path factor
     ELSE                               ! if sun near/below horizon,
       XPATHS(J)=ZPTHMX                 ! set to very large value
     ENDIF

ENDDO
!-----------------------------------------------------------------------



!-----1.2  Sets up surface type info.-----------------------------------

DO J=1,KNPF
 NSTYPE(J)=KSURF(J)
ENDDO
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         2.   UNPACK PROFILE VECTOR.
!-----------------------------------------------------------------------
DO  J=1,KNPF

  DO JL=1,43
   TEMP(JL,J)=PROF(JL,J)           ! Temperature
  ENDDO

    DO JL=44,86
     WMIX(JL-43,J)=PROF(JL,J)      ! W. vapour volume mixing ratio [ppmv]
    ENDDO

      DO JL=87,129
       OMIX(JL-86,J)=PROF(JL,J)    ! Ozone volume mixing ratio [ppmv]
      END DO

ENDDO

        DO  J=1,KNPF
          TA(J)       =PROF(131,J) ! Surface temperature in K
          WMIXS(J)    =PROF(132,J) ! Not currently used
          TS(J)       =PROF(136,J) ! Skin temperature in K
          SURFP(J)    =PROF(133,J) ! Surface pressure in mb=hPa
          CLDP(J)     =PROF(137,J) ! Cloud-top pressure in mb=hPa
          CLDF(J)     =PROF(138,J) ! Fractional (ir) cloud cover
          CLW(J)      =PROF(130,J) ! Not currently used
          SURFWU(J)   =PROF(134,J) ! U component of surface wind in m/sec
          SURFWV(J)   =PROF(135,J) ! V component of surface wind in m/sec
        ENDDO


!-----------------------------------------------------------------------
!         3.   SET RELATED PROFILE PARAMETERS.
!-----------------------------------------------------------------------


!-----3.1  Set surface pressure parameters.-----------------------------

CALL RTIASI_PRSLEV(1,SURFP,NLEVSF,FRACPS,KNPF,IPRINT)
!-----------------------------------------------------------------------


!-----3.2  Set cloud pressure parameters.-------------------------------

CALL RTIASI_PRSLEV(2,CLDP,NLEVCD,FRACPC,KNPF,IPRINT)
!-----------------------------------------------------------------------


!-----3.3  Set microwave cloud parameters (not yet coded)---------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!         4.   SETS UP PROFILE VARIABLES FOR TRANSMITTANCE CALC.
!-----------------------------------------------------------------------

CALL RTIASI_PRFTAU(KNPF,XXM,XXW,XXW1,XXO,XXO1,PWWR)
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_PRFIN
!     Set up profile variables for radiative transfer calculations.
!SUBROUTINE RTIASI_PRFINK(PROF,KNPF,KNCHPF,KPROF,KLENPF, &
!                 XXM_D,XXM,XXW_D,XXW,XXW1_D,XXW1,XXO_D,XXO,XXO1_D, &
!                 XXO1,PWWR_D,PWWR)
SUBROUTINE RTIASI_PRFINK(PROF,KNPF,KNCHPF,KPROF,KLENPF, &
                 XXM_D,XXM,XXW_D,XXW,XXW1,XXO_D,XXO, &
                 XXO1,PWWR_D,PWWR)

!     Description:
!     K of subroutine PRFIN
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     2) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/1990   Original code. J.R.EYRE. ECMWF
!     2           12/11/1999   Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct module used:

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCOFM            , & ! Max. number of mixed gas coefficients
  JPCOFW            , & ! Max. number of water vapour coefficients
  JPCOFO            , & ! Max. number of ozone coefficients
  JPLEV             , & ! Number of levels
  JPPF              , & ! Max. number of processed profiles
  JPCHPF              ! Max. number of processed radiances

USE RTIASI_IRCLD, ONLY : &
!     Imported arrays:
  NLEVCD_D=>NLEVCD     ! Index of nearest standard pressure level
                       ! at/below cloud top.

USE RTIASI_SURF, ONLY : &
!     Imported arrays:
       NLEVSF_D =>NLEVSF  ! Index of nearest standard pressure level
                          ! at/below surface.
!     K module used:


USE RTIASI_PRFVARK,ONLY : &
!     Imported arrays:
  TEMP              , & ! Temperature profile in K
  WMIX              , & ! Water vapour volume mixing ratio in ppmv
  OMIX              , & ! Ozone volume mixing ratio in ppmv
  TA                , & ! Surface air temperature in K
  WMIXS             , & ! Water vapour surface volume mixing ratio in ppmv
  TS                , & ! Skin temperature IN K
  SURFP             , & ! Surface pressure in mb=hPa
  CLDP              , & ! Cloud-top pressure in mb=hPa
  CLDF              , & ! Fractional (ir) cloud cover
  CLW               , & ! Total column cloud liquid water in mm
  SURFWU            , & ! U component of surface wind in m/sec
  SURFWV              ! V component of surface wind in m/sec


USE RTIASI_SURFK, ONLY : &
!     Imported arrays:
  FRACPS              ! Fraction of standard pressure level interval by
                      ! which surf is above level nlevsf.

USE RTIASI_IRCLDK, ONLY : &
!     Imported arrays:
  FRACPC              ! Fraction of standardd pressure level interval
                      ! by which cloud is above level NLEVCD.


IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF                      ! Number of processed
                                                ! profiles.
INTEGER,INTENT(IN) :: KLENPF                    ! Lenght of input profile
                                                ! to RTIASI.
INTEGER, INTENT(IN):: KNCHPF                    ! Number of processed
                                                ! radiances.


!       Array arguments with intent in:
INTEGER, INTENT(IN):: KPROF (KNCHPF)            ! Array of profile indices


!       K arguments with intent out:
REAL,INTENT(OUT)   :: PROF  (KLENPF,KNCHPF)     ! Input profile to RTIASI


!       K array arguments with intent in:
REAL,INTENT(INOUT) :: XXM    (JPCOFM,JPLEV,JPCHPF) ! Functions of profile
                                                  ! for fixed gas
                                                  ! calculations.
REAL,INTENT(INOUT) :: XXW    (JPCOFW,JPLEV,JPCHPF) ! Functions of profile
                                                  ! for water vapour
                                                  ! calculations.
REAL,INTENT(INOUT) :: XXO    (JPCOFO,JPLEV,JPCHPF) ! Functions of profile
                                                  ! for ozone calculations.
REAL,INTENT(INOUT) :: XXW1   (JPCOFW,JPLEV,JPCHPF) ! Functions of profile:
                                                  ! level 1, water vapour.
REAL,INTENT(INOUT) :: XXO1   (JPCOFO,JPLEV,JPCHPF) ! Functions of profile:
                                                  ! level 1, ozone.
REAL,INTENT(INOUT) :: PWWR   (JPLEV,JPCHPF)        ! Functions of profile:
                                                  ! pressure weighted water

!       Array arguments with intent in:
REAL,INTENT(IN)   :: XXM_D  (JPCOFM,JPLEV,JPPF) ! Functions of profile
                                                ! for fixed gas
                                                ! calculations.
REAL,INTENT(IN)   :: XXW_D  (JPCOFW,JPLEV,JPPF) ! Functions of profile
                                                ! for water vapour
                                                ! calculations.
REAL,INTENT(IN)   :: XXO_D  (JPCOFO,JPLEV,JPPF) ! Functions of profile
                                                ! for ozone calculations.
REAL,INTENT(IN)   :: PWWR_D (JPLEV,JPPF)        ! Functions of profile:
                                                ! pressure weighted water
                                                ! vapour ratio.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,JL


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   SETS UP PROFILE-DEPENDENT 'KNOWN' QUANTITIES.
!-----------------------------------------------------------------------



!-----1.1  Satellite and solar geometry.--------------------------------
!           No Tangent Linear required
!-----------------------------------------------------------------------



!-----1.2  Sets up surface type info.-----------------------------------
!           No Tangent Linear required
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         4.   SETS UP PROFILE VARIABLES FOR TRANSMITTANCE CALC.
!-----------------------------------------------------------------------

!CALL PRFTAUK(KNPF,KNCHPF,KPROF,XXM_D,XXM,XXW_D,XXW,XXW1_D,XXW1, &
!              XXO_D,XXO,XXO1_D,XXO1,PWWR_D,PWWR)
CALL RTIASI_PRFTAUK(KNPF,KNCHPF,KPROF,XXM_D,XXM,XXW_D,XXW,XXW1, &
              XXO_D,XXO,XXO1,PWWR_D,PWWR)
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         3.   SET RELATED PROFILE PARAMETERS.
!-----------------------------------------------------------------------


!-----3.3  Set microwave cloud parameters (not yet coded)---------------
!-----------------------------------------------------------------------


!-----3.2  Set cloud pressure parameters.-------------------------------

!CALL PRSLEVK(2,CLDP,NLEVCD_D,FRACPC,KNPF,KNCHPF,KPROF)
CALL RTIASI_PRSLEVK(CLDP,NLEVCD_D,FRACPC,KNPF,KNCHPF,KPROF)
!-----------------------------------------------------------------------



!-----3.1  Set surface pressure parameters.-----------------------------

!CALL PRSLEVK(1,SURFP,NLEVSF_D,FRACPS,KNPF,KNCHPF,KPROF)
CALL RTIASI_PRSLEVK(SURFP,NLEVSF_D,FRACPS,KNPF,KNCHPF,KPROF)
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!         2.   UNPACK PROFILE VECTOR.
!-----------------------------------------------------------------------



DO  J=1,KNCHPF
  PROF(135,J)=SURFWV(J) ! V component of surface wind in m/sec
  PROF(134,J)=SURFWU(J) ! U component of surface wind in m/sec
  PROF(130,J)=CLW(J)    ! Not currently used
  PROF(138,J)=CLDF(J)   ! Fractional (ir) cloud cover
  PROF(137,J)=CLDP(J)   ! Cloud-top pressure in mb=hPa
  PROF(133,J)=SURFP(J)  ! Surface pressure in mb=hPa
  PROF(136,J)=TS(J)     ! Skin temperature in K
  PROF(132,J)=WMIXS(J)  ! Not currently used
  PROF(131,J)=TA(J)     ! Surface temperature in K
ENDDO


DO  J=1,KNCHPF

  DO JL=87,129
    PROF(JL,J)=OMIX(JL-86,J)
  ENDDO

    DO JL=44,86
      PROF(JL,J)=WMIX(JL-43,J)
    ENDDO

      DO JL=1,43
        PROF(JL,J)=TEMP(JL,J)
      ENDDO
ENDDO



RETURN
END SUBROUTINE RTIASI_PRFINK
!     Locate given pressures on array of fixed levels.
SUBROUTINE RTIASI_PRSLEV(ISWITCH,PRES,KLEV,PFRAC,KNPF,IPRINT)

!     Description:
!     Finds the indices and fractional levels locating
!     given pressures on an array of fixed pressure levels.
!     Exception: if pressure greater than highest value in fixed
!     levels, set KLEV=last level, and then PFRAC is negative.
!
!     Method.
!     Linear interpolation in pressure.
!     1) See: User's manual for RTIASI (Available from EUMETSAT)
!     2) See: ECMWF TECH MEM 176.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code(RTTOV). J.R.EYRE. ECMWF
!     2            12/11/1999  Move to FORTRAN 90. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_PRFCON,ONLY : &
!     Imported parameters:
 JPLEV   ,   & ! Number of pressure levels
!     Imported arrays:
  XPRES  ,   & ! Standard pressure levels for transmittance (and, currently,
!                  ! radiative treansfer) calculation; from top down; in mb=hpa
  DPRES      ! Intervals between standard pressure levels; in mb


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: IPRINT
INTEGER, INTENT(IN) :: KNPF         ! Number of profiles to be processed.

!       Array arguments with intent in:
REAL,INTENT(IN)     :: PRES  (KNPF) ! Input array of pressures
                                    ! to be located.

!       Array arguments with intent out:
INTEGER,INTENT(OUT) :: KLEV  (KNPF) ! Output array of level indices.
REAL,INTENT(OUT)    :: PFRAC (KNPF) ! Output array of fractions by which
                                    ! given pressures are between level
                                    ! KLEV and adjacent level of lower
                                    ! pressure.

!     End of subroutine arguments


!      Local scalars:
INTEGER :: J,JL,ILEV,ISWITCH


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   FIND NEAREST LEVEL AND CALCULATE FRACTIONAL PRESSURE.
!-----------------------------------------------------------------------

DO J=1,KNPF
  DO JL=JPLEV,2,-1
    IF (PRES(J) > XPRES(JL)) THEN
      ILEV=MIN((JL+1),JPLEV)
      KLEV(J)=ILEV
        IF(IPRINT == 1) THEN
        IF(ISWITCH == 2)THEN
          WRITE(*,'(A31,I2)')'NEAREST LEVEL TO TOP CLOUD IS: ',KLEV(J)
        ELSE IF(ISWITCH == 1)THEN
          WRITE(*,'(A29,I2)')'NEAREST LEVEL TO SURFACE IS: ',KLEV(J)
        END IF
        ENDIF
      PFRAC(J)=(XPRES(ILEV)-PRES(J))/DPRES(ILEV)
      EXIT
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_PRSLEV
!     Locate given pressures on array of fixed levels.
!SUBROUTINE PRSLEVK(ISWITCH,PRES,KLEV,PFRAC,KNPF,KNCHPF,KPROF)
SUBROUTINE RTIASI_PRSLEVK(PRES,KLEV,PFRAC,KNPF,KNCHPF,KPROF)

!     Description:
!     K of subroutine PRSLEV
!     Finds the indices and fractional levels locating
!     given pressures on an array of fixed pressure levels.
!     Exception: if pressure greater than highest value in fixed,
!     levels, set KLEV=last level, and then PFRAC is negative.
!
!     Method.
!     Linear interpolation in pressure.
!     1) See: User's manual for RTIASI (Available from EUMETSAT)
!     2) See: ECMWF TECH MEM 176.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code. J.R.EYRE. ECMWF
!     2            12/11/1999  Move to FORTRAN 90. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_PRFCON,ONLY : &
!     Imported arrays:
  DPRES      ! Intervals between standard pressure levels; in mb

IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KNPF             ! Number of profiles to be processed.
INTEGER, INTENT(IN) :: KNCHPF           ! Number of processed radiances.

!       Array arguments with intent in:
INTEGER,INTENT(IN)  :: KLEV  (KNPF)     ! Input array of level indices.
INTEGER, INTENT(IN) :: KPROF (KNCHPF)   ! Array of profile indices.

!       K array arguments with intent in:
REAL,INTENT(INOUT)  :: PFRAC (KNCHPF)   ! Output array of fractions by which
                                        ! given pressures are between level
                                        ! KLEV and adjacent level of lower
                                        ! pressure.

!       K array arguments with intent out:
REAL,INTENT(INOUT)  :: PRES  (KNCHPF)   ! Input array of pressures
                                        ! to be located.

!     End of subroutine arguments


!      Local scalars:
INTEGER :: J,ILEV,JPK


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   FIND NEAREST LEVEL AND CALCULATE FRACTIONAL PRESSURE.
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  JPK=KPROF(J)
  ILEV=KLEV(JPK)
  PRES(J)=PRES(J)-PFRAC(J)/DPRES(ILEV)
ENDDO

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_PRSLEVK
!     Store profile variables for transmittance calc.
SUBROUTINE RTIASI_PRFTAU(KNPF,XXM,XXW,XXW1,XXO,XXO1,PWWR)

!     Description:
!     Calculate and store the profile variables required
!     in subsequent transmittance calculations.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            09/04/91.   Original code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOPTH, ONLY : &
!     Imported arrays:
  XPATH      ! Secant of viewing path angle at surface


USE RTIASI_PRFVAR,ONLY : &
!     Imported parameters:
  JPLEV   ,  & ! Number of pressure levels
  JPPF    ,  & ! Max. number of profiles to be processed
!     Imported arrays:
  TEMP    ,  & ! Temperature profile in K
  WMIX    ,  & ! Water vapour volume mixing ratio in ppmv
  OMIX         ! Ozone volume mixing ratio in ppmv

USE RTIASI_PRFCON,ONLY : &
!     Imported arrays:
  DPP        ! DPRES*XPRES


USE RTIASI_PRFREF,ONLY : &
!     Imported arrays:
  TREF    ,  & ! Reference temperature profile in K
  TREFO   ,  & ! Reference temperature for ozone in K
  WREF    ,  & ! Reference water vapour volume mixing ratio in ppmv
  WREFO   ,  & ! Reference water vapour for ozone in ppmv
  OREF       ! Reference ozone volume mixing ratio in ppmv


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCOFM  ,  & ! Max. number of mixed gas coefficients
  JPCOFW  ,  & ! Max. number of water vapour coefficients
  JPCOFO     ! Max. number of ozone coefficients


USE RTIASI_PREDVAR,ONLY: &
!     Profile dependent variables used to set up
!     IASI predictors.
!     Imported parameters:
  RATT, &
  RATO, &
  RATW, &
  PWTR, &
  PWOR, &
  PWORTR, &
  DELT, &
  DELTO, &
  RATWO, &
  RATTO, &
  PUWOR


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF                    ! Number of processed
                                              ! profiles.
!       Array arguments with intent out:
REAL,INTENT(OUT)   :: XXM(JPCOFM,JPLEV,JPPF)  ! Functions of profile for
                                              ! fixed gas calculations.
REAL,INTENT(OUT)   :: XXW(JPCOFW,JPLEV,JPPF)  ! Functions of profile for
                                              ! wat. vapour calculations.
REAL,INTENT(OUT)   :: XXO(JPCOFO,JPLEV,JPPF)  ! Functions of profile for
                                              ! ozone calculations.
REAL,INTENT(OUT)   :: XXW1(JPCOFW,JPLEV,JPPF) ! Functions of profile:
                                              ! level 1, water vapour.
REAL,INTENT(OUT)   :: XXO1(JPCOFO,JPLEV,JPPF) ! Functions of profile:
                                              ! level 1, ozone.
REAL,INTENT(OUT)   :: PWWR(JPLEV,JPPF)        ! Functions of profile:
                                              ! pressure weighted
                                              ! water vapour ratio.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: JK,JL
REAL    :: DPPRT
REAL    :: DPPWP
REAL    :: DPPWR
REAL    :: OU
REAL    :: OUR
REAL    :: DPPO
REAL    :: DPPOR
REAL    :: DPPRORT

!       Local arrays:
REAL    :: SUM1(JPPF)
REAL    :: SUM2(JPPF)
REAL    :: SUM3(JPPF)
REAL    :: SUM4(JPPF)
REAL    :: SUM5(JPPF)
REAL    :: SUM6(JPPF)
REAL    :: SUM7(JPPF)
REAL    :: SUM8(JPPF)


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!        1.   INTERPOLATE PROFILE TO TRANSMITTANCE P-LEVELS.
!-----------------------------------------------------------------------
!          Not required yet.  input, trans and rt p-levels all same



!-----------------------------------------------------------------------
!        2.   CALCULATE VARIABLES FOR TRANSMITTANCE CALCULATION.
!-----------------------------------------------------------------------

SUM1(:)=0.
SUM2(:)=0.
SUM3(:)=0.
SUM4(:)=0.
SUM5(:)=0.
SUM6(:)=0.
SUM7(:)=0.
SUM8(:)=0.


DO JL=1,JPLEV
  DO JK=1,KNPF
    IF(JL == 1)THEN

!-----------SET UP PROFILE VARIABLES FOR MIXED GASES--------------------

      RATT(JL,JK)=TEMP(JL,JK)/TREF(JL)
      DPPRT=0.

      SUM1(JK)=SUM1(JK)+DPPRT
      PWTR(JL,JK)=SUM1(JK)

!-----------SET UP PROFILE VARIABLES FOR WATER VAPOUR-------------------

      RATW(JL,JK)=WMIX(JL,JK)/WREF(JL)
      DELT(JL,JK)=TEMP(JL,JK)-TREF(JL)

      DPPWP=(DPP(JL)*WMIX(JL,JK))
      DPPWR=(DPP(JL)*WREF(JL))

      SUM2(JK)=SUM2(JK)+DPPWP
      SUM3(JK)=SUM3(JK)+DPPWR
      PWWR(JL,JK)=SUM2(JK)/SUM3(JK)

!-----------SET UP PROFILE VARIABLES FOR OZONE--------------------------

      RATO(JL,JK)=OMIX(JL,JK)/OREF(JL)
      DELTO(JL,JK)=TEMP(JL,JK)-TREFO(JL)
      RATWO(JL,JK)=WMIX(JL,JK)/WREFO(JL)
      RATTO(JL,JK)=TEMP(JL,JK)/TREFO(JL)

      DPPO=DPP(JL)*OMIX(JL,JK)
      DPPOR=DPP(JL)*OREF(JL)
      OU=OMIX(JL,JK)
      OUR=OREF(JL)
      DPPRORT=0.

      SUM4(JK)=SUM4(JK)+OU
      SUM5(JK)=SUM5(JK)+OUR
      PUWOR(JL,JK)=SUM4(JK)/SUM5(JK)

      SUM6(JK)=SUM6(JK)+DPPO
      SUM7(JK)=SUM7(JK)+DPPOR
      PWOR(JL,JK)=SUM6(JK)/SUM7(JK)

      SUM8(JK)=SUM8(JK)+DPPRORT
      PWORTR(JL,JK)=SUM8(JK)

!-----------------------------------------------------------------------

    ELSE IF(JL > 1)THEN

!-----------SET UP PROFILE VARIABLES FOR MIXED GASES--------------------

      RATT(JL,JK)=(TEMP(JL,JK)+TEMP(JL-1,JK))/(TREF(JL)+TREF(JL-1) &
     )
      DPPRT=DPP(JL)*RATT(JL-1,JK)

      SUM1(JK)=SUM1(JK)+DPPRT
      PWTR(JL,JK)=SUM1(JK)

!-----------SET UP PROFILE VARIABLES FOR WATER VAPOUR-------------------

      RATW(JL,JK)=(WMIX(JL,JK)+WMIX(JL-1,JK))/(WREF(JL)+WREF(JL-1) &
      )
      DELT(JL,JK)=0.5*(TEMP(JL,JK)+TEMP(JL-1,JK)-TREF(JL)-TREF(JL- &
      1))

      DPPWP=DPP(JL)*0.5*(WMIX(JL,JK)+WMIX(JL-1,JK))
      DPPWR=DPP(JL)*0.5*(WREF(JL)+WREF(JL-1))

      SUM2(JK)=SUM2(JK)+DPPWP
      SUM3(JK)=SUM3(JK)+DPPWR
      PWWR(JL,JK)=SUM2(JK)/SUM3(JK)

!-----------SET UP PROFILE VARIABLES FOR OZONE--------------------------

      RATO(JL,JK)=(OMIX(JL,JK)+OMIX(JL-1,JK))/(OREF(JL)+OREF(JL-1) &
      )
      DELTO(JL,JK)=0.5*(TEMP(JL,JK)+TEMP(JL-1,JK)-TREFO(JL)-TREFO( &
      JL-1))
      RATWO(JL,JK)=(WMIX(JL,JK)+WMIX(JL-1,JK))/(WREFO(JL)+WREFO(JL &
      -1))
      RATTO(JL,JK)=(TEMP(JL,JK)+TEMP(JL-1,JK))/(TREFO(JL)+TREFO(JL &
      -1))

      DPPO=DPP(JL)*0.5*(OMIX(JL,JK)+OMIX(JL-1,JK))
      DPPOR=DPP(JL)*0.5*(OREF(JL)+OREF(JL-1))
      OU=0.5*(OMIX(JL,JK)+OMIX(JL-1,JK))
      OUR=0.5*(OREF(JL)+OREF(JL-1))
      DPPRORT=DPP(JL)*DELTO(JL-1,JK)*RATO(JL-1,JK)

      SUM4(JK)=SUM4(JK)+OU
      SUM5(JK)=SUM5(JK)+OUR
      PUWOR(JL,JK)=SUM4(JK)/SUM5(JK)

      SUM6(JK)=SUM6(JK)+DPPO
      SUM7(JK)=SUM7(JK)+DPPOR
      PWOR(JL,JK)=SUM6(JK)/SUM7(JK)

      SUM8(JK)=SUM8(JK)+DPPRORT
      PWORTR(JL,JK)=SUM8(JK)
    END IF
  ENDDO                                    ! End of profile loop
ENDDO                                      ! End of level loop



DO JL=1,JPLEV
  DO JK=1,KNPF

!---------PREDICTORS FOR MIXED GAS TRANSMITTANCES-----------------------

    XXM(1,JL,JK)=XPATH(JK)
    XXM(2,JL,JK)=XPATH(JK)*XPATH(JK)
    XXM(3,JL,JK)=XPATH(JK)*RATT(JL,JK)
    XXM(4,JL,JK)=XPATH(JK)*RATT(JL,JK)*RATT(JL,JK)
    XXM(5,JL,JK)=RATT(JL,JK)
    XXM(6,JL,JK)=RATT(JL,JK)*RATT(JL,JK)
    XXM(7,JL,JK)=XPATH(JK)*PWTR(JL,JK)
    XXM(8,JL,JK)=XPATH(JK)*PWTR(JL,JK)/RATT(JL,JK)
    XXM(9,JL,JK)=SQRT(XPATH(JK))
    XXM(10,JL,JK)=SQRT(XPATH(JK))*PWTR(JL,JK)**0.25

!-----------------------------------------------------------------------



!---------PREDICTORS FOR WATER VAPOUR TRANSMITTANCE---------------------

    IF(JL == 1)THEN      !LAYER SPECIFIC PREDICTORS ARE USED FOR LAYER 1

      XXW1(1,JL,JK)=XPATH(JK)*RATW(JL,JK)
      XXW1(2,JL,JK)=1

      XXW(1,JL,JK)=0.
      XXW(2,JL,JK)=0.
      XXW(3,JL,JK)=0.
      XXW(4,JL,JK)=0.
      XXW(5,JL,JK)=0.
      XXW(6,JL,JK)=0.
      XXW(7,JL,JK)=0.
      XXW(8,JL,JK)=0.
      XXW(9,JL,JK)=0.
      XXW(10,JL,JK)=0.
      XXW(11,JL,JK)=0.
      XXW(12,JL,JK)=0.
      XXW(13,JL,JK)=0.
      XXW(14,JL,JK)=0.
    ELSE
      XXW(1,JL,JK)=XPATH(JK)*RATW(JL,JK)
      XXW(2,JL,JK)=(XPATH(JK)*RATW(JL,JK))**0.5
      XXW(3,JL,JK)=XPATH(JK)*RATW(JL,JK)**2/PWWR(JL,JK)
      XXW(4,JL,JK)=XPATH(JK)*RATW(JL,JK)*DELT(JL,JK)
      XXW(5,JL,JK)=(XPATH(JK)*RATW(JL,JK))**2
      XXW(6,JL,JK)=(XPATH(JK)*RATW(JL,JK))**0.5*DELT(JL,JK)
      XXW(7,JL,JK)=(XPATH(JK)*RATW(JL,JK))**0.25
      XXW(8,JL,JK)=(XPATH(JK)*RATW(JL,JK))**0.5*RATW(JL,JK)/PWWR(JL,JK)
      XXW(9,JL,JK)=(XPATH(JK)*RATW(JL,JK))**3
      XXW(10,JL,JK)=RATW(JL,JK)
      XXW(11,JL,JK)=XPATH(JK)*RATW(JL,JK)*DELT(JL,JK)*ABS(DELT(JL,JK))
      XXW(12,JL,JK)=XPATH(JK)*RATW(JL,JK)**2
      XXW(13,JL,JK)=(XPATH(JK)*PWWR(JL,JK))**4
      XXW(14,JL,JK)=(XPATH(JK)*PWWR(JL,JK))**2
    END IF

!-----------------------------------------------------------------------



!---------PREDICTORS FOR OZONE TRANSMITTANCE----------------------------

    IF(JL == 1)THEN      !LAYER SPECIFIC PREDICTORS ARE USED FOR LAYER 1

      XXO1(1,JL,JK)=XPATH(JK)*RATO(JL,JK)
      XXO1(2,JL,JK)=1

      XXO(1,JL,JK)=0.
      XXO(2,JL,JK)=0.
      XXO(3,JL,JK)=0.
      XXO(4,JL,JK)=0.
      XXO(5,JL,JK)=0.
      XXO(6,JL,JK)=0.
      XXO(7,JL,JK)=0.
      XXO(8,JL,JK)=0.
      XXO(9,JL,JK)=0.
    ELSE
      XXO(1,JL,JK)=XPATH(JK)*RATO(JL,JK)
      XXO(2,JL,JK)=(XPATH(JK)*RATO(JL,JK))**0.5
      XXO(3,JL,JK)=XPATH(JK)*RATO(JL,JK)*DELTO(JL,JK)
      XXO(4,JL,JK)=(XPATH(JK)*RATO(JL,JK))**2
      XXO(5,JL,JK)=(XPATH(JK)*RATO(JL,JK))**0.5*DELTO(JL,JK)
      XXO(6,JL,JK)=XPATH(JK)*RATO(JL,JK)*RATO(JL,JK)*PUWOR(JL,JK)
      XXO(7,JL,JK)=SQRT(XPATH(JK)*RATO(JL,JK))*RATO(JL,JK)/PUWOR(JL,JK)
      XXO(8,JL,JK)=XPATH(JK)*RATO(JL,JK)*PWOR(JL,JK)/PUWOR(JL,JK)
      XXO(9,JL,JK)=XPATH(JK)*RATO(JL,JK)*SQRT(XPATH(JK)*PUWOR(JL,JK))
      XXO(10,JL,JK)=XPATH(JK)*XPATH(JK)*RATO(JL,JK)*PWORTR(JL,JK)
    END IF
  ENDDO                                     ! End of level loop
ENDDO                                       ! End of profile loop

RETURN
END SUBROUTINE RTIASI_PRFTAU
!     Store profile variables for transmittance calc.

!SUBROUTINE PRFTAUK(KNPF,KNCHPF,KPROF,XXM_D,XXM,XXW_D,XXW,XXW1_D, &
!                   XXW1,XXO_D,XXO,XXO1_D,XXO1,PWWR_D,PWWR)
SUBROUTINE RTIASI_PRFTAUK(KNPF,KNCHPF,KPROF,XXM_D,XXM,XXW_D,XXW, &
                   XXW1,XXO_D,XXO,XXO1,PWWR_D,PWWR)

!     Description:
!     K of subroutine  PRFTAU
!     Calculate and store the profile variables required
!     in subsequent transmittance calculations.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            09/04/91.   Original code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct module used:


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCHPF                ! Max. number of processed radiances

USE RTIASI_GEOPTH, ONLY : &
!     Imported arrays:
  XPATH              ! Secant of viewing path angle at surface

USE RTIASI_PRFVAR,ONLY : &
!     Imported parameters:
  JPLEV             , & ! Number of pressure levels
  JPPF              , & ! Max. number of profiles to be processed
!     Imported arrays:
  WMIX_D   =>WMIX   , & ! Water vapour volume mixing ratio in ppmv
  OMIX_D   =>OMIX       ! Ozone volume mixing ratio in ppmv

USE RTIASI_PRFCON,ONLY : &
!     Imported arrays:
  DPP                 ! DPRES*XPRES


USE RTIASI_PRFREF,ONLY : &
!     Imported arrays:
  TREF              , & ! Reference temperature profile in K
  TREFO             , & ! Reference temperature for ozone in K
  WREF              , & ! Reference water vapour volume mixing ratio in ppmv
  WREFO             , & ! Reference water vapour for ozone in ppmv
  OREF                ! Reference ozone volume mixing ratio in ppmv


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCOFM            , & ! Max. number of mixed gas coefficients
  JPCOFW            , & ! Max. number of water vapour coefficients
  JPCOFO              ! Max. number of ozone coefficients


USE RTIASI_PREDVAR,ONLY : &
!     Profile dependent variables used to set up
!     IASI predictors
!     Imported parameters:
  RATT_D   =>RATT   , &
  RATO_D   =>RATO   , &
  PWOR_D   =>PWOR   , &
  PWORTR_D =>PWORTR , &
  DELT_D   =>DELT   , &
  DELTO_D  =>DELTO  , &
  PUWOR_D  =>PUWOR


!    K module used:

USE RTIASI_PRFVARK,ONLY : &
!     Imported arrays:
  TEMP              , & ! Temperature profile in K
  WMIX              , & ! Water vapour volume mixing ratio in ppmv
  OMIX                  ! Ozone volume mixing ratio in ppmv


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF                         ! Number of processed
INTEGER, INTENT(IN):: KNCHPF                       ! Number of processed
                                                   ! radiances.

!       Array arguments with intent in:
INTEGER, INTENT(IN):: KPROF  (KNCHPF)              ! Array of profile
                                                   ! indices.
!       Direct array arguments with intent in:
REAL,INTENT(IN)    :: XXM_D  (JPCOFM,JPLEV,JPPF)   ! Functions of prof for
                                                   ! fixed gas calc.
REAL,INTENT(IN)    :: XXW_D  (JPCOFW,JPLEV,JPPF)   ! Functions of prof for
                                                   ! wat. vapour calc.
REAL,INTENT(IN)    :: XXO_D  (JPCOFO,JPLEV,JPPF)   ! Functions of prof for
                                                   ! ozone calcs.
REAL,INTENT(IN)    :: PWWR_D (JPLEV,JPPF)          ! Functions of profile:
                                                   ! pressure weighted

!       K array with intent in:
REAL,INTENT(INOUT) :: XXM    (JPCOFM,JPLEV,JPCHPF) ! Functions of prof for
                                                   ! fixed gas calc.
REAL,INTENT(INOUT) :: XXW    (JPCOFW,JPLEV,JPCHPF) ! Functions of prof for
                                                   ! wat. vapour calc.
REAL,INTENT(INOUT) :: XXO    (JPCOFO,JPLEV,JPCHPF) ! Functions of prof for
                                                   ! ozone calc.
REAL,INTENT(INOUT) :: XXW1   (JPCOFW,JPLEV,JPCHPF) ! Functions of profile:
                                                   ! level 1, water vapour
REAL,INTENT(INOUT) :: XXO1   (JPCOFO,JPLEV,JPCHPF) ! Functions of profile:
                                                   ! level 1, ozone.
REAL,INTENT(INOUT) :: PWWR   (JPLEV,JPCHPF)        ! Functions of profile:
                                                   ! pressure weighted
                                                   ! water vapour ratio.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: JK,JL,JPK

!       Local tangent linear scalars:
REAL    :: DPPWP
REAL    :: DPPWR
REAL    :: DPPO
REAL    :: DPPOR
REAL    :: OU
REAL    :: OUR
REAL    :: DPPRT
REAL    :: DPPRORT

!       Local arrays:
REAL    :: SUM2_D  (JPLEV,JPPF)
REAL    :: SUM3_D  (JPLEV,JPPF)
REAL    :: SUM4_D  (JPLEV,JPPF)
REAL    :: SUM5_D  (JPLEV,JPPF)
REAL    :: SUM6_D  (JPLEV,JPPF)
REAL    :: SUM7_D  (JPLEV,JPPF)
REAL    :: DPPWP_D (JPLEV,JPPF)
REAL    :: DPPWR_D (JPLEV,JPPF)
REAL    :: DPPO_D  (JPLEV,JPPF)
REAL    :: DPPOR_D (JPLEV,JPPF)
REAL    :: OU_D    (JPLEV,JPPF)
REAL    :: OUR_D   (JPLEV,JPPF)

!       Local tangent linear arrays:
REAL    :: SUM1    (JPCHPF)
REAL    :: SUM2    (JPCHPF)
REAL    :: SUM3    (JPCHPF)
REAL    :: SUM4    (JPCHPF)
REAL    :: SUM5    (JPCHPF)
REAL    :: SUM6    (JPCHPF)
REAL    :: SUM7    (JPCHPF)
REAL    :: SUM8    (JPCHPF)
REAL    :: RATT    (JPLEV,JPCHPF)
REAL    :: RATW    (JPLEV,JPCHPF)
REAL    :: PWTR    (JPLEV,JPCHPF)
REAL    :: RATWO   (JPLEV,JPCHPF)
REAL    :: RATTO   (JPLEV,JPCHPF)
REAL    :: RATO    (JPLEV,JPCHPF)
REAL    :: DELT    (JPLEV,JPCHPF)
REAL    :: PWOR    (JPLEV,JPCHPF)
REAL    :: PUWOR   (JPLEV,JPCHPF)
REAL    :: DELTO   (JPLEV,JPCHPF)
REAL    :: PWORTR  (JPLEV,JPCHPF)




!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!        1.   INTERPOLATE PROFILE TO TRANSMITTANCE P-LEVELS.
!-----------------------------------------------------------------------
!          Not required yet.  input, trans and rt p-levels all same





!-----------------------------------------------------------------------
!        2.   CALCULATE VARIABLES FOR TRANSMITTANCE CALCULATION.
!-----------------------------------------------------------------------


!-----Initialize variables----------------------------------------------

SUM2_D(:,:)=0.
SUM3_D(:,:)=0.
SUM4_D(:,:)=0.
SUM5_D(:,:)=0.
SUM6_D(:,:)=0.
SUM7_D(:,:)=0.

SUM1(:)    =0.
SUM2(:)    =0.
SUM3(:)    =0.
SUM4(:)    =0.
SUM5(:)    =0.
SUM6(:)    =0.
SUM7(:)    =0.
SUM8(:)    =0.

RATO(:,:)  =0.
PWORTR(:,:)=0.
PUWOR(:,:) =0.
PWOR(:,:)  =0.
DELTO(:,:) =0.
RATW(:,:)  =0.
DELT(:,:)  =0.
PWTR(:,:)  =0.
RATT(:,:)  =0.
DPPRORT    =0.
DPPOR      =0.
DPPO       =0.
OUR        =0.
OU         =0.
DPPWR      =0.
DPPWP      =0.
DPPRT      =0.
RATTO(:,:) =0.
RATWO(:,:) =0
!-----------------------------------------------------------------------





!-----Compute direct variables------------------------------------------

DO JL=1,JPLEV
  DO JK=1,KNPF
    IF(JL == 1)THEN

      DPPWP_D(JL,JK)=(DPP(JL)*WMIX_D(JL,JK))
      DPPWR_D(JL,JK)=(DPP(JL)*WREF(JL))
      SUM2_D(JL,JK)=DPPWP_D(JL,JK)
      SUM3_D(JL,JK)=DPPWR_D(JL,JK)

      DPPO_D(JL,JK)=DPP(JL)*OMIX_D(JL,JK)
      DPPOR_D(JL,JK)=DPP(JL)*OREF(JL)
      OU_D(JL,JK)=OMIX_D(JL,JK)
      OUR_D(JL,JK)=OREF(JL)
      SUM4_D(JL,JK)=OU_D(JL,JK)
      SUM5_D(JL,JK)=OUR_D(JL,JK)
      SUM6_D(JL,JK)=DPPO_D(JL,JK)
      SUM7_D(JL,JK)=DPPOR_D(JL,JK)

    ELSE IF(JL > 1)THEN

      DPPWP_D(JL,JK)=DPP(JL)*0.5*(WMIX_D(JL,JK)+WMIX_D(JL-1,JK))
      DPPWR_D(JL,JK)=DPP(JL)*0.5*(WREF(JL)+WREF(JL-1))
      SUM2_D(JL,JK)=SUM2_D(JL-1,JK)+DPPWP_D(JL,JK)
      SUM3_D(JL,JK)=SUM3_D(JL-1,JK)+DPPWR_D(JL,JK)

      DPPO_D(JL,JK)=DPP(JL)*0.5*(OMIX_D(JL,JK)+OMIX_D(JL-1,JK))
      DPPOR_D(JL,JK)=DPP(JL)*0.5*(OREF(JL)+OREF(JL-1))
      OU_D(JL,JK)=0.5*(OMIX_D(JL,JK)+OMIX_D(JL-1,JK))
      OUR_D(JL,JK)=0.5*(OREF(JL)+OREF(JL-1))
      SUM4_D(JL,JK)=SUM4_D(JL-1,JK)+OU_D(JL,JK)
      SUM5_D(JL,JK)=SUM5_D(JL-1,JK)+OUR_D(JL,JK)
      SUM6_D(JL,JK)=SUM6_D(JL-1,JK)+DPPO_D(JL,JK)
      SUM7_D(JL,JK)=SUM7_D(JL-1,JK)+DPPOR_D(JL,JK)
    END IF
  ENDDO
ENDDO

!-----------------------------------------------------------------------




!-----Predictors for ozone transmittance--------------------------------

DO JL=JPLEV,1,-1
  DO JK=1,KNCHPF
     JPK=KPROF(JK)

    IF(JL == 1)THEN      !Layer specific predictors are used for layer 1

      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*XXO1(1,JL,JK)
      XXO1(2,JL,JK)=0.

      XXO(1,JL,JK)=0.
      XXO(2,JL,JK)=0.
      XXO(3,JL,JK)=0.
      XXO(4,JL,JK)=0.
      XXO(5,JL,JK)=0.
      XXO(6,JL,JK)=0.
      XXO(7,JL,JK)=0.
      XXO(8,JL,JK)=0.
      XXO(9,JL,JK)=0.
      XXO(10,JL,JK)=0.

    ELSE


      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*XPATH(JPK)* &
      PWORTR_D(JL,JPK)*XXO(10,JL,JK)
      PWORTR(JL,JK)=PWORTR(JL,JK)+XPATH(JPK)*XXO_D(1,JL,JPK)* &
      XXO(10,JL,JK)

      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*SQRT(XPATH(JPK)* &
      PUWOR_D(JL,JPK))*XXO(9,JL,JK)
      PUWOR(JL,JK)=PUWOR(JL,JK)+XXO_D(1,JL,JPK)*0.5*XPATH(JPK)* &
      XXO(9,JL,JK)/SQRT(XPATH(JPK)*PUWOR_D(JL,JPK))

      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*XXO(8,JL,JK)* &
      PWOR_D(JL,JPK)/PUWOR_D(JL,JPK)
      PWOR(JL,JK)=PWOR(JL,JK)+XXO_D(1,JL,JPK)*XXO(8,JL,JK)/ &
      PUWOR_D(JL,JPK)
      PUWOR(JL,JK)=PUWOR(JL,JK)-XXO_D(1,JL,JPK)*PWOR_D(JL,JPK)* &
      XXO(8,JL,JK)/PUWOR_D(JL,JPK)**2

      RATO(JL,JK)=RATO(JL,JK)+1.5*XXO_D(2,JL,JPK)*XXO(7,JL,JK)/ &
      PUWOR_D(JL,JPK)
      PUWOR(JL,JK)=PUWOR(JL,JK)-XXO_D(7,JL,JPK)*XXO(7,JL,JK)/ &
      PUWOR_D(JL,JPK)

      RATO(JL,JK)=RATO(JL,JK)+2*XXO_D(1,JL,JPK)*XXO(6,JL,JK)* &
      PUWOR_D(JL,JPK)
      PUWOR(JL,JK)=PUWOR(JL,JK)+XXO(6,JL,JK)*XXO_D(4,JL,JPK)/ &
      XPATH(JPK)

      DELTO(JL,JK)=DELTO(JL,JK)+XXO_D(2,JL,JPK)*XXO(5,JL,JK)
      RATO(JL,JK)=RATO(JL,JK)+ DELTO_D(JL,JPK)*0.5*XPATH(JPK)* &
      XXO(5,JL,JK)/XXO_D(2,JL,JPK)

      RATO(JL,JK)=RATO(JL,JK)+2*XPATH(JPK)*XXO(4,JL,JK)* &
      XXO_D(1,JL,JPK)

      DELTO(JL,JK)=DELTO(JL,JK)+XXO_D(1,JL,JPK)*XXO(3,JL,JK)
      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*XXO(3,JL,JK)* &
      DELTO_D(JL,JPK)

      RATO(JL,JK)=RATO(JL,JK)+0.5*XPATH(JPK)*XXO(2,JL,JK)/ &
      XXO_D(2,JL,JPK)

      RATO(JL,JK)=RATO(JL,JK)+XPATH(JPK)*XXO(1,JL,JK)

    ENDIF

!-----Predictors for water vapour transmittance-------------------------

    IF(JL == 1)THEN

      RATW(JL,JK)=RATW(JL,JK)+XPATH(JPK)*XXW1(1,JL,JK)
      XXW1(2,JL,JK)=0.

      XXW(1,JL,JK)=0.
      XXW(2,JL,JK)=0.
      XXW(3,JL,JK)=0.
      XXW(4,JL,JK)=0.
      XXW(5,JL,JK)=0.
      XXW(6,JL,JK)=0.
      XXW(7,JL,JK)=0.
      XXW(8,JL,JK)=0.
      XXW(9,JL,JK)=0.
      XXW(10,JL,JK)=0.
      XXW(11,JL,JK)=0.
      XXW(12,JL,JK)=0.
      XXW(13,JL,JK)=0.
      XXW(14,JL,JK)=0.

    ELSE

      PWWR(JL,JK)=PWWR(JL,JK)+2*(XPATH(JPK)*PWWR_D(JL,JPK))* &
      XPATH(JPK)*XXW(14,JL,JK)

      PWWR(JL,JK)=PWWR(JL,JK)+4*(XPATH(JPK)*PWWR_D(JL,JPK))**3* &
      XPATH(JPK)*XXW(13,JL,JK)

      RATW(JL,JK)=RATW(JL,JK)+2*XXW_D(1,JL,JPK)*XXW(12,JL,JK)

      IF ( DELT_D(JL,JPK) >= 0. ) THEN
        RATW(JL,JK)=RATW(JL,JK)+XPATH(JPK)*XXW(11,JL,JK)* &
        DELT_D(JL,JPK)*DELT_D(JL,JPK)
        DELT(JL,JK)=DELT(JL,JK)+XXW_D(1,JL,JPK)*2* &
        DELT_D(JL,JPK)*XXW(11,JL,JK)
      ELSE IF ( DELT_D(JL,JPK) < 0. ) THEN
        RATW(JL,JK)=RATW(JL,JK)+XPATH(JPK)*XXW(11,JL,JK)* &
        DELT_D(JL,JPK)*ABS(DELT_D(JL,JPK))
        DELT(JL,JK)=DELT(JL,JK)+XXW_D(1,JL,JPK)*XXW(11,JL,JK)* &
        ABS(DELT_D(JL,JPK))
        DELT(JL,JK)=DELT(JL,JK)-XXW_D(1,JL,JPK)*DELT_D(JL,JPK)* &
        XXW(11,JL,JK)
      ENDIF

     RATW(JL,JK)=RATW(JL,JK)+XXW(10,JL,JK)

      RATW(JL,JK)=RATW(JL,JK)+3*XXW_D(1,JL,JPK)**2*XPATH(JPK)* &
      XXW(9,JL,JK)

      RATW(JL,JK)=RATW(JL,JK)+1.5*XXW_D(2,JL,JPK)*XXW(8,JL,JK)/ &
      PWWR_D(JL,JPK)
      PWWR(JL,JK)=PWWR(JL,JK)-XXW_D(2,JL,JPK)**3*XXW(8,JL,JK)/ &
      (XPATH(JPK)*PWWR_D(JL,JPK)**2)

      RATW(JL,JK)=RATW(JL,JK)+XPATH(JPK)*XXW(7,JL,JK)*0.25/ &
      XXW_D(7,JL,JPK)**3

      DELT(JL,JK)=DELT(JL,JK)+XXW_D(2,JL,JPK)*XXW(6,JL,JK)
      RATW(JL,JK)=RATW(JL,JK)+0.5*DELT_D(JL,JPK)*XPATH(JPK)* &
      XXW(6,JL,JK)/XXW_D(2,JL,JPK)

      RATW(JL,JK)=RATW(JL,JK)+2*XXW_D(1,JL,JPK)*XXW(5,JL,JK)* &
      XPATH(JPK)

      DELT(JL,JK)=DELT(JL,JK)+XXW_D(1,JL,JPK)*XXW(4,JL,JK)
      RATW(JL,JK)=RATW(JL,JK)+XXW(4,JL,JK)*XPATH(JPK)* &
      DELT_D(JL,JPK)

      RATW(JL,JK)=RATW(JL,JK)+(2*XXW_D(1,JL,JPK)*XXW(3,JL,JK)/ &
      PWWR_D(JL,JPK))
      PWWR(JL,JK)=PWWR(JL,JK)-XXW_D(3,JL,JPK)*XXW(3,JL,JK)/ &
      PWWR_D(JL,JPK)

      RATW(JL,JK)=RATW(JL,JK)+0.5*XPATH(JPK)*XXW(2,JL,JK)/ &
      XXW_D(2,JL,JPK)

      RATW(JL,JK)=RATW(JL,JK)+XPATH(JPK)*XXW(1,JL,JK)

    ENDIF

!-----Predictors for mixed gas transmittances---------------------------


    IF (JL /= 1)THEN
      PWTR(JL,JK)=PWTR(JL,JK)+(XPATH(JPK)**2)*0.25*XXM(10,JL,JK)/ &
      (XXM_D(10,JL,JPK))**3
    ELSE IF (JL == 1 )THEN
      XXM(10,JL,JK)=0.
    ENDIF

    XXM(9,JL,JK)=0.

    PWTR(JL,JK)=PWTR(JL,JK)+XPATH(JPK)*XXM(8,JL,JK)/RATT_D(JL,JPK)
    RATT(JL,JK)=RATT(JL,JK)-XXM_D(8,JL,JPK)*XXM(8,JL,JK)/ &
    XXM_D(5,JL,JPK)

    PWTR(JL,JK)=PWTR(JL,JK)+XPATH(JPK)*XXM(7,JL,JK)

    RATT(JL,JK)=RATT(JL,JK)+2*XXM(6,JL,JK)*XXM_D(5,JL,JPK)

    RATT(JL,JK)=RATT(JL,JK)+XXM(5,JL,JK)

    RATT(JL,JK)=RATT(JL,JK)+2*XXM(4,JL,JK)*XXM_D(3,JL,JPK)

    RATT(JL,JK)=RATT(JL,JK)+XPATH(JPK)*XXM(3,JL,JK)

    XXM(2,JL,JK)=0.

    XXM(1,JL,JK)=0.

  ENDDO
ENDDO



DO JL=JPLEV,1,-1
  DO JK=1,KNCHPF
     JPK=KPROF(JK)
    IF(JL == 1)THEN

!-----Set up profile variables for ozone--------------------------------


      SUM8(JK)=SUM8(JK)+PWORTR(JL,JK)
      PWORTR(JL,JK)=0.

      DPPRORT=DPPRORT+SUM8(JK)

      SUM6(JK)=SUM6(JK)+PWOR(JL,JK)/SUM7_D(JL,JPK)
      SUM7(JK)=SUM7(JK)+SUM6_D(JL,JPK)*PWOR(JL,JK)/ &
      SUM7_D(JL,JPK)**2
      PWOR(JL,JK)=0

      DPPOR=DPPOR+SUM7(JK)
      DPPO=DPPO+SUM6(JK)

      SUM4(JK)=SUM4(JK)+PUWOR(JL,JK)/SUM5_D(JL,JPK)
      SUM5(JK)=SUM5(JK)+SUM4_D(JL,JPK)*PUWOR(JL,JK)/ &
      SUM5_D(JL,JPK)**2
      PUWOR(JL,JK)=0.

      OUR=OUR+SUM5(JK)
      OU=OU+SUM4(JK)

      DPPRORT=0.
      OUR=0.


      OMIX(JL,JK)=OMIX(JL,JK)+OU
      OU=0.
      DPPOR=0.
      OMIX(JL,JK)=OMIX(JL,JK)+DPP(JL)*DPPO
      DPPO=0.


      TEMP(JL,JK)=TEMP(JL,JK)+RATTO(JL,JK)/TREFO(JL)
      RATTO(JL,JK)=0.

      WMIX(JL,JK)=WMIX(JL,JK)+RATWO(JL,JK)/WREFO(JL)
      RATWO(JL,JK)=0.

      TEMP(JL,JK)=TEMP(JL,JK)+DELTO(JL,JK)
      DELTO(JL,JK)=0.

      OMIX(JL,JK)=OMIX(JL,JK)+RATO(JL,JK)/OREF(JL)
      RATO(JL,JK)=0.



!-----Set up profile variables for water vapour-------------------------


      SUM2(JK)=SUM2(JK)+PWWR(JL,JK)/SUM3_D(JL,JPK)
      SUM3(JK)=SUM3(JK)-SUM2_D(JL,JPK)*PWWR(JL,JK)/ &
      SUM3_D(JL,JPK)**2
      PWWR(JL,JK)=0.

      DPPWR=DPPWR+SUM3(JK)

      DPPWP=DPPWP+SUM2(JK)

      DPPWR=0.

      WMIX(JL,JK)=WMIX(JL,JK)+DPP(JL)*DPPWP
      DPPWP=0.


      TEMP(JL,JK)=TEMP(JL,JK)+DELT(JL,JK)
      DELT(JL,JK)=0.

      WMIX(JL,JK)=WMIX(JL,JK)+RATW(JL,JK)/WREF(JL)
      RATW(JL,JK)=0.

!-----Set up profile variables for mixed gases--------------------------

      SUM1(JK)=SUM1(JK)+PWTR(JL,JK)
      PWTR(JL,JK)=0.

      DPPRT=DPPRT+SUM1(JK)
      DPPRT=0.

      TEMP(JL,JK)=TEMP(JL,JK)+RATT(JL,JK)/TREF(JL)
      RATT(JL,JK)=0.


    ELSE IF(JL > 1)THEN

!-----Set up profile variables for ozone--------------------------------



      SUM8(JK)=SUM8(JK)+PWORTR(JL,JK)
      PWORTR(JL,JK)=0.

      DPPRORT=DPPRORT+SUM8(JK)

      SUM6(JK)=SUM6(JK)+PWOR(JL,JK)/SUM7_D(JL,JPK)
      SUM7(JK)=SUM7(JK)+SUM6_D(JL,JPK)*PWOR(JL,JK)/ &
      SUM7_D(JL,JPK)**2
      PWOR(JL,JK)=0.

      DPPOR=DPPOR+SUM7(JK)

      DPPO=DPPO+SUM6(JK)

      SUM4(JK)=SUM4(JK)+PUWOR(JL,JK)/SUM5_D(JL,JPK)
      SUM5(JK)=SUM5(JK)+SUM4_D(JL,JPK)*PUWOR(JL,JK)/ &
      SUM5_D(JL,JPK)**2
      PUWOR(JL,JK)=0.

      OUR=OUR+SUM5(JK)

      OU=OU+SUM4(JK)

      DELTO(JL-1,JK)=DELTO(JL-1,JK)+DPP(JL)*DPPRORT* &
      RATO_D(JL-1,JPK)
      RATO(JL-1,JK)=RATO(JL-1,JK)+DPP(JL)*DELTO_D(JL-1,JPK)* &
      DPPRORT
      DPPRORT=0.
      OUR=0.

      OMIX(JL,JK)=OMIX(JL,JK)+0.5*OU
      OMIX(JL-1,JK)=OMIX(JL-1,JK)+0.5*OU
      OU=0.

      DPPOR=0.

      OMIX(JL,JK)=OMIX(JL,JK)+DPP(JL)*0.5*DPPO
      OMIX(JL-1,JK)=OMIX(JL-1,JK)+DPP(JL)*0.5*DPPO
      DPPO=0

      TEMP(JL,JK)=TEMP(JL,JK)+RATTO(JL,JK)/(TREFO(JL)+TREFO(JL &
      -1))

     TEMP(JL-1,JK)=TEMP(JL-1,JK)+RATTO(JL,JK)/(TREFO(JL)+TREFO(JL &
      -1))
      RATTO(JL,JK)=0.


       WMIX(JL,JK)=WMIX(JL,JK)+RATWO(JL,JK)/(WREFO(JL)+WREFO(JL &
      -1))
      WMIX(JL-1,JK)=WMIX(JL-1,JK)+RATWO(JL,JK)/(WREFO(JL)+WREFO(JL &
      -1))
      
      RATWO(JL,JK)=0.

      TEMP(JL,JK)=TEMP(JL,JK)+0.5*DELTO(JL,JK)
      TEMP(JL-1,JK)=TEMP(JL-1,JK)+0.5*DELTO(JL,JK)
      DELTO(JL,JK)=0.


       OMIX(JL,JK)=OMIX(JL,JK)+RATO(JL,JK)/(OREF(JL)+OREF(JL-1) &
      )
      OMIX(JL-1,JK)=OMIX(JL-1,JK)+RATO(JL,JK)/(OREF(JL)+OREF(JL-1) &
      )
       RATO(JL,JK)=0




!-----Set up profile variables for water vapour-------------------------


      SUM2(JK)=SUM2(JK)+PWWR(JL,JK)/SUM3_D(JL,JPK)
      SUM3(JK)=SUM3(JK)-SUM2_D(JL,JPK)*PWWR(JL,JK)/ &
      SUM3_D(JL,JPK)**2
      PWWR(JL,JK)=0.

      DPPWR=DPPWR+SUM3(JK)

      DPPWP=DPPWP+SUM2(JK)

      DPPWR=0.


      WMIX(JL,JK)=WMIX(JL,JK)+DPP(JL)*0.5*DPPWP
      WMIX(JL-1,JK)=WMIX(JL-1,JK)+DPP(JL)*0.5*DPPWP
      DPPWP=0.

      TEMP(JL,JK)=TEMP(JL,JK)+0.5*DELT(JL,JK)
      TEMP(JL-1,JK)=TEMP(JL-1,JK)+0.5*DELT(JL,JK)
      DELT(JL,JK)=0.

      WMIX(JL,JK)=WMIX(JL,JK)+RATW(JL,JK)/(WREF(JL)+WREF(JL-1))
      WMIX(JL-1,JK)=WMIX(JL-1,JK)+RATW(JL,JK)/ &
      (WREF(JL)+WREF(JL-1))
      RATW(JL,JK)=0

!-----Set up profile variables for mixed gases--------------------------

      SUM1(JK)=SUM1(JK)+PWTR(JL,JK)
      PWTR(JL,JK)=0.

      DPPRT=DPPRT+SUM1(JK)

      RATT(JL-1,JK)=RATT(JL-1,JK)+DPP(JL)*DPPRT
      DPPRT=0

      TEMP(JL,JK)=TEMP(JL,JK)+RATT(JL,JK)/(TREF(JL)+TREF(JL-1))
      TEMP(JL-1,JK)=TEMP(JL-1,JK)+RATT(JL,JK)/ &
      (TREF(JL)+TREF(JL-1))

      RATT(JL,JK)=0.


    ENDIF
  ENDDO
ENDDO


RETURN
END SUBROUTINE RTIASI_PRFTAUK
!     Calculate optical depths.
SUBROUTINE RTIASI_OPDEP(KCHAN,KPROF,KNCHPF,KSAT,XXM,XXW,XXW1,XXO, &
XXO1,OPDPMP,OPDPO,OPDPWP1,OPDPWP2,OPDPWP)

!     Description:
!     Calculate optical depths for a number of channels
!     and profiles from every pressure level to space.
!
!     Method:
!     Layer optical depths are computed for ozone,mixed gases and water
!     wapour.For mixed gases one set of optical depth are computed to
!     account for positive transmittances.For water vapour
!     three sets are computed to account
!     for three different regimes of positive transmittances.One set
!     is computed for ozone.
!     Note that these are fictious layer optical depths in that they
!     come from ratioing layer to space transmittances that have been
!     shifted by a constant value (currently 0.06) to ensure we are
!     dealing with positive ransmittances.The effective layer to space
!     transmittances are recovered using the algorithm described
!     in subroutine RTTAU.
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2)See: User's manual for RTIASI (Available from EUMETSAT)
!     3)See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Moved to F90. Also see "Method".
!                              Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOPTH, ONLY : &
!     Imported parameters:
  JPPF         ! Max no. profiles

USE RTIASI_PRFCON,ONLY : &
!     Imported parameters:
 JPLEV    ,    & ! Number of pressure levels
!     Imported scalars:
  NLEVW        ! Upper level for water vapour transmittance calculation
!                    ! (assumed transmittance=0. above this)


USE RTIASI_TAUCFN, ONLY : &
!     Imported parameters:
  JPCOFM  ,    & ! Max. number of mixed gas coefficients
  JPCOFW  ,    & ! Max. number of water vapour coefficients
  JPCOFO  ,    & ! Max. number of ozone coefficients
!     Imported arrays
  CFMPOS  ,    & ! Mixed gases coeffs(positive transmittances)
  CFWPOS  ,    & ! Water vapour coeffs(single regression scheme)
  CFWPOS1 ,    & ! Water vapour coeffs(optically thin regime)
  CFWPOS2 ,    & ! Water vapour coeffs(optically thick regime)
  CFO          ! Ozone coeffs


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KSAT                        ! Input satellite index
INTEGER, INTENT(IN) :: KNCHPF                      ! Input number of rad.
                                                   ! (channels * profiles)

!       Array arguments with intent in:
INTEGER, INTENT(IN) :: KCHAN   (KNCHPF)            ! Input array of ch.
                                                   ! indices for each rad.
INTEGER, INTENT(IN) :: KPROF   (KNCHPF)            ! Input array of prof.
                                                   ! indices for each rad.
REAL,INTENT(IN)     :: XXM     (JPCOFM,JPLEV,JPPF) ! Functions of prof.
                                                   ! for fixed gas calc.
REAL,INTENT(IN)     :: XXW     (JPCOFW,JPLEV,JPPF) ! Functions of prof.
                                                   ! for wat. vapour calc.
REAL,INTENT(IN)     :: XXO     (JPCOFO,JPLEV,JPPF) ! Functions of profile
                                                   ! for ozone calc.
REAL,INTENT(IN)     :: XXW1    (JPCOFW,JPLEV,JPPF) ! Functions of profile:
                                                   ! level 1, wat. vapour.
REAL,INTENT(IN)     :: XXO1    (JPCOFO,JPLEV,JPPF) ! Functions of profile:
                                                   ! level 1, ozone.

!       Array arguments with intent out:
REAL,INTENT(OUT)    :: OPDPMP  (KNCHPF,JPLEV)      ! Mixed gases op. depth
REAL,INTENT(OUT)    :: OPDPO   (KNCHPF,JPLEV)      ! Ozone optical depth.
REAL,INTENT(OUT)    :: OPDPWP1 (KNCHPF,JPLEV)      ! Wat. vapour op. depth
                                                   ! :weak absorption.
REAL,INTENT(OUT)    :: OPDPWP2 (KNCHPF,JPLEV)      ! Wat. vapour op. depth
                                                   ! :strong absorption.
REAL,INTENT(OUT)    :: OPDPWP  (KNCHPF,JPLEV)      ! Wat. vapour op. depth
                                                   ! :single regr. scheme

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,JL,ICH,IPF

!       Local arrays:
REAL    :: ZOPDP (KNCHPF,JPLEV)
REAL    :: OPDPA (KNCHPF)


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE GASEOUS LAYER OPTICAL DEPTHS.
!-----------------------------------------------------------------------



!-----1.1  Uniformly mixed gases:positive transmittances----------------
!
OPDPA(1:KNCHPF)=0.

DO J=1,KNCHPF
  ICH=KCHAN(J)
  IPF=KPROF(J)
    DO JL=1,JPLEV
          ZOPDP(J,JL)= &
               CFMPOS(1,JL,ICH,KSAT)*XXM(1,JL,IPF)+ &
               CFMPOS(2,JL,ICH,KSAT)*XXM(2,JL,IPF)+ &
               CFMPOS(3,JL,ICH,KSAT)*XXM(3,JL,IPF)+ &
               CFMPOS(4,JL,ICH,KSAT)*XXM(4,JL,IPF)+ &
               CFMPOS(5,JL,ICH,KSAT)*XXM(5,JL,IPF)+ &
               CFMPOS(6,JL,ICH,KSAT)*XXM(6,JL,IPF)+ &
               CFMPOS(7,JL,ICH,KSAT)*XXM(7,JL,IPF)+ &
               CFMPOS(8,JL,ICH,KSAT)*XXM(8,JL,IPF)+ &
               CFMPOS(9,JL,ICH,KSAT)*XXM(9,JL,IPF)+ &
               CFMPOS(10,JL,ICH,KSAT)*XXM(10,JL,IPF)

    ENDDO
ENDDO

DO JL=1,JPLEV
  DO J=1,KNCHPF
    OPDPA(J)=OPDPA(J)+ZOPDP(J,JL)
    OPDPMP(J,JL)=ZOPDP(J,JL)
  ENDDO
ENDDO
!-----------------------------------------------------------------------


!-----1.2  Water vapour:positive transmittances regime 1----------------

OPDPA(1:KNCHPF)=0.

DO J=1,KNCHPF
  ICH=KCHAN(J)
  IPF=KPROF(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
        ZOPDP(J,JL)=CFWPOS1(1,JL,ICH,KSAT)*XXW1(1,JL,IPF)+ &
        CFWPOS1(2,JL,ICH,KSAT)*XXW1(2,JL,IPF)
      ELSE
          ZOPDP(J,JL)= &
               CFWPOS1(1,JL,ICH,KSAT)*XXW(1,JL,IPF)+ &
               CFWPOS1(2,JL,ICH,KSAT)*XXW(2,JL,IPF)+ &
               CFWPOS1(3,JL,ICH,KSAT)*XXW(3,JL,IPF)+ &
               CFWPOS1(4,JL,ICH,KSAT)*XXW(4,JL,IPF)+ &
               CFWPOS1(5,JL,ICH,KSAT)*XXW(5,JL,IPF)+ &
               CFWPOS1(6,JL,ICH,KSAT)*XXW(6,JL,IPF)+ &
               CFWPOS1(7,JL,ICH,KSAT)*XXW(7,JL,IPF)+ &
               CFWPOS1(8,JL,ICH,KSAT)*XXW(8,JL,IPF)+ &
               CFWPOS1(9,JL,ICH,KSAT)*XXW(9,JL,IPF)+ &
               CFWPOS1(10,JL,ICH,KSAT)*XXW(10,JL,IPF)+ &
               CFWPOS1(11,JL,ICH,KSAT)*XXW(11,JL,IPF)+ &
               CFWPOS1(12,JL,ICH,KSAT)*XXW(12,JL,IPF)+ &
               CFWPOS1(13,JL,ICH,KSAT)*XXW(13,JL,IPF)+ &
               CFWPOS1(14,JL,ICH,KSAT)*XXW(14,JL,IPF)

      END IF
    ENDDO
ENDDO


DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    OPDPA(J)=OPDPA(J)+ZOPDP(J,JL)
    OPDPWP1(J,JL)=ZOPDP(J,JL)
  ENDDO
ENDDO
!-----------------------------------------------------------------------


!-----1.3  Water vapour:positive transmittances regime 2----------------

OPDPA(1:KNCHPF)=0.

DO J=1,KNCHPF
  ICH=KCHAN(J)
  IPF=KPROF(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
        ZOPDP(J,JL)=CFWPOS2(1,JL,ICH,KSAT)*XXW1(1,JL,IPF)+ &
        CFWPOS2(2,JL,ICH,KSAT)*XXW1(2,JL,IPF)
      ELSE IF(JL /= 1)THEN
          ZOPDP(J,JL)= &
               CFWPOS2(1,JL,ICH,KSAT)*XXW(1,JL,IPF)+ &
               CFWPOS2(2,JL,ICH,KSAT)*XXW(2,JL,IPF)+ &
               CFWPOS2(3,JL,ICH,KSAT)*XXW(3,JL,IPF)+ &
               CFWPOS2(4,JL,ICH,KSAT)*XXW(4,JL,IPF)+ &
               CFWPOS2(5,JL,ICH,KSAT)*XXW(5,JL,IPF)+ &
               CFWPOS2(6,JL,ICH,KSAT)*XXW(6,JL,IPF)+ &
               CFWPOS2(7,JL,ICH,KSAT)*XXW(7,JL,IPF)+ &
               CFWPOS2(8,JL,ICH,KSAT)*XXW(8,JL,IPF)+ &
               CFWPOS2(9,JL,ICH,KSAT)*XXW(9,JL,IPF)+ &
               CFWPOS2(10,JL,ICH,KSAT)*XXW(10,JL,IPF)+ &
               CFWPOS2(11,JL,ICH,KSAT)*XXW(11,JL,IPF)+ &
               CFWPOS2(12,JL,ICH,KSAT)*XXW(12,JL,IPF)+ &
               CFWPOS2(13,JL,ICH,KSAT)*XXW(13,JL,IPF)+ &
               CFWPOS2(14,JL,ICH,KSAT)*XXW(14,JL,IPF)
      END IF
    ENDDO
ENDDO

DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    OPDPA(J)=OPDPA(J)+ZOPDP(J,JL)
    OPDPWP2(J,JL)=ZOPDP(J,JL)
 ENDDO
ENDDO
!-----------------------------------------------------------------------


!-----1.3  Water vapour:single regression scheme -----------------------

OPDPA(1:KNCHPF)=0.

DO J=1,KNCHPF
  ICH=KCHAN(J)
  IPF=KPROF(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
        ZOPDP(J,JL)=CFWPOS(1,JL,ICH,KSAT)*XXW1(1,JL,IPF)+ &
        CFWPOS(2,JL,ICH,KSAT)*XXW1(2,JL,IPF)
      ELSE IF(JL /= 1)THEN
          ZOPDP(J,JL)= &
               CFWPOS(1,JL,ICH,KSAT)*XXW(1,JL,IPF)+ &
               CFWPOS(2,JL,ICH,KSAT)*XXW(2,JL,IPF)+ &
               CFWPOS(3,JL,ICH,KSAT)*XXW(3,JL,IPF)+ &
               CFWPOS(4,JL,ICH,KSAT)*XXW(4,JL,IPF)+ &
               CFWPOS(5,JL,ICH,KSAT)*XXW(5,JL,IPF)+ &
               CFWPOS(6,JL,ICH,KSAT)*XXW(6,JL,IPF)+ &
               CFWPOS(7,JL,ICH,KSAT)*XXW(7,JL,IPF)+ &
               CFWPOS(8,JL,ICH,KSAT)*XXW(8,JL,IPF)+ &
               CFWPOS(9,JL,ICH,KSAT)*XXW(9,JL,IPF)+ &
               CFWPOS(10,JL,ICH,KSAT)*XXW(10,JL,IPF)+ &
               CFWPOS(11,JL,ICH,KSAT)*XXW(11,JL,IPF)+ &
               CFWPOS(12,JL,ICH,KSAT)*XXW(12,JL,IPF)+ &
               CFWPOS(13,JL,ICH,KSAT)*XXW(13,JL,IPF)+ &
               CFWPOS(14,JL,ICH,KSAT)*XXW(14,JL,IPF)
      END IF
    ENDDO
ENDDO

DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    OPDPA(J)=OPDPA(J)+ZOPDP(J,JL)
    OPDPWP(J,JL)=ZOPDP(J,JL)
 ENDDO
ENDDO
!-----------------------------------------------------------------------



!-----1.5  Ozone--------------------------------------------------------

OPDPA(1:KNCHPF)=0.

DO J=1,KNCHPF
  ICH=KCHAN(J)
  IPF=KPROF(J)
    DO JL=1,JPLEV
      IF(JL == 1)THEN
        ZOPDP(J,JL)= &
        CFO(1,JL,ICH,KSAT)*XXO1(1,JL,IPF)+ &
        CFO(2,JL,ICH,KSAT)*XXO1(2,JL,IPF)
      ELSE IF(JL /= 1)THEN
        ZOPDP(J,JL)= &
               CFO(1,JL,ICH,KSAT)*XXO(1,JL,IPF)+ &
               CFO(2,JL,ICH,KSAT)*XXO(2,JL,IPF)+ &
               CFO(3,JL,ICH,KSAT)*XXO(3,JL,IPF)+ &
               CFO(4,JL,ICH,KSAT)*XXO(4,JL,IPF)+ &
               CFO(5,JL,ICH,KSAT)*XXO(5,JL,IPF)+ &
               CFO(6,JL,ICH,KSAT)*XXO(6,JL,IPF)+ &
               CFO(7,JL,ICH,KSAT)*XXO(7,JL,IPF)+ &
               CFO(8,JL,ICH,KSAT)*XXO(8,JL,IPF)+ &
               CFO(9,JL,ICH,KSAT)*XXO(9,JL,IPF)+ &
               CFO(10,JL,ICH,KSAT)*XXO(10,JL,IPF)

      END IF
    ENDDO
 ENDDO

DO JL=1,JPLEV
  DO J=1,KNCHPF
    OPDPA(J)=OPDPA(J)+ZOPDP(J,JL)
    OPDPO(J,JL)=ZOPDP(J,JL)
  ENDDO
ENDDO
!-----------------------------------------------------------------------


RETURN
END SUBROUTINE RTIASI_OPDEP
!     Calculate optical depths.
!SUBROUTINE OPDEPK(KCHAN,KPROF,KNCHPF,KSAT,XXM_D,XXM,XXW_D,XXW, &
!XXW1_D,XXW1,XXO_D,XXO,XXO1_D,XXO1,OPDPMP_D,OPDPMP,OPDPO_D,OPDPO, &
!OPDPWP1_D,OPDPWP1,OPDPWP2_D,OPDPWP2,OPDPWP_D,OPDPWP)
SUBROUTINE RTIASI_OPDEPK(KCHAN,KNCHPF,KSAT,XXM,XXW, &
XXW1,XXO,XXO1,OPDPMP,OPDPO, &
OPDPWP1,OPDPWP2,OPDPWP)

!     Description:
!     K of subroutine OPDEP.
!     Calculate optical depths for a number of channels
!     and profiles from every pressure level to space.
!
!     Method:
!     Layer optical depths are computed for ozone,mixed gases and water
!     wapour.For mixed gases one set of optical depth are computed to
!     account for positive transmittances.For water vapour
!     three sets are computed to account
!     for three different regimes of positive transmittances.One set
!     is computed for ozone.
!     Note that these are fictious layer optical depths in that they
!     come from ratioing layer to space transmittances that have been
!     shifted by a constant value (currently 0.06) to ensure we are
!     dealing with positive ransmittances.The effective layer to space
!     transmittances are recovered using the algorithm described
!     in subroutine RTTAU.
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2)See: User's manual for RTIASI (Available from EUMETSAT)
!     3)See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Moved to F90. Also see "Method".
!                              Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOPTH, ONLY : &
!     Imported parameters:
  JPCHPF       ! Max no. profiles


USE RTIASI_PRFCON,ONLY : &
!     Imported parameters:
  JPLEV,     & ! Number of pressure levels
!     Imported scalars:
  NLEVW        ! Upper level for water vapour transmittance calculation
!                  ! (assumed transmittance=0. above this)



USE RTIASI_TAUCFN, ONLY : &
!     Imported parameters:
  JPCOFM,    & ! Max. number of mixed gas coefficients
  JPCOFW,    & ! Max. number of water vapour coefficients
  JPCOFO,    & ! Max. number of ozone coefficients
!     Imported arrays
  CFMPOS,    & ! Mixed gases coeffs(positive transmittances)
  CFWPOS,    & ! Water vapour coeffs(single reg. scheme)
  CFWPOS1,   & ! Water vapour coeffs(positive transmittances)
  CFWPOS2,   & ! Water vapour coeffs(positive transmittances)
  CFO        ! Ozone coeffs


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KSAT                           ! Input sat. index.
INTEGER, INTENT(IN) :: KNCHPF                         ! Input number of
                                                      ! radiances
                                                      ! (channels * profiles).

!       Array arguments with intent in:
INTEGER, INTENT(IN) :: KCHAN(KNCHPF)                  ! Input array of
                                                      ! channel indices
                                                      ! for each radiance

!       Tangent linear array arguments with intent in:
REAL,INTENT(INOUT)  :: XXM      (JPCOFM,JPLEV,JPCHPF) ! Functions of prof
                                                      ! for fixed gas
                                                      ! calculations.
REAL,INTENT(INOUT)  :: XXW      (JPCOFW,JPLEV,JPCHPF) ! Functions of prof
                                                      ! for wat. vapour
                                                      ! calculations.
REAL,INTENT(INOUT)  :: XXO      (JPCOFO,JPLEV,JPCHPF) ! Functions of prof
                                                      ! for ozone
                                                      ! calculations.
REAL,INTENT(INOUT)  :: XXW1     (JPCOFW,JPLEV,JPCHPF) ! Functions of prof:
                                                      ! level 1, wat. vap.
REAL,INTENT(INOUT)  :: XXO1     (JPCOFO,JPLEV,JPCHPF) ! Functions of prof:
                                                      ! level 1, ozone.

!       Tangent linear array arguments with intent out:
REAL,INTENT(INOUT)  :: OPDPMP   (KNCHPF,JPLEV)        ! Mix gases op depth
REAL,INTENT(INOUT)  :: OPDPO    (KNCHPF,JPLEV)        ! Ozone op depth.
REAL,INTENT(INOUT)  :: OPDPWP1  (KNCHPF,JPLEV)        ! Wat. vap op depth
                                                      ! :weak absorption.
REAL,INTENT(INOUT)  :: OPDPWP2  (KNCHPF,JPLEV)        ! Wat. vap op depth
                                                      ! :strong absorption
REAL,INTENT(INOUT)  :: OPDPWP   (KNCHPF,JPLEV)        ! Wat. vap op depth
                                                      ! :sing. reg. scheme

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,JL,ICH


!       Tangent linera local arrays:
REAL    :: ZOPDP (KNCHPF,JPLEV)
REAL    :: OPDPA (KNCHPF)

OPDPA(:)   =0.
ZOPDP(:,:) =0.



!-----End of header-----------------------------------------------------



!-----------------------------------------------------------------------
!         2.   Liquid water (microwave channels).
!              Not yet coded.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE GASEOUS LAYER OPTICAL DEPTHS.
!-----------------------------------------------------------------------

!-----1.6  Ozone--------------------------------------------------------

DO JL=1,JPLEV
  DO J=1,KNCHPF
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPO(J,JL)
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPA(J)
  ENDDO
ENDDO


DO J=1,KNCHPF
  ICH=KCHAN(J)
    DO JL=1,JPLEV
      IF(JL == 1)THEN
        XXO1(1,JL,J)=XXO1(1,JL,J)+CFO(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO1(2,JL,J)=XXO1(2,JL,J)+CFO(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      ELSE IF(JL /= 1)THEN
        XXO(10,JL,J)=XXO(10,JL,J)+CFO(10,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(9,JL,J) =XXO(9,JL,J)+CFO(9,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(8,JL,J) =XXO(8,JL,J)+CFO(8,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(7,JL,J) =XXO(7,JL,J)+CFO(7,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(6,JL,J) =XXO(6,JL,J)+CFO(6,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(5,JL,J) =XXO(5,JL,J)+CFO(5,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(4,JL,J) =XXO(4,JL,J)+CFO(4,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(3,JL,J) =XXO(3,JL,J)+CFO(3,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(2,JL,J) =XXO(2,JL,J)+CFO(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXO(1,JL,J) =XXO(1,JL,J)+CFO(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      END IF
    ENDDO
ENDDO

ZOPDP(:,:)=0.


!-----1.5  Water vapour:single regression scheme -----------------------

DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPWP(J,JL)
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPA(J)
 ENDDO
ENDDO


DO J=1,KNCHPF
  ICH=KCHAN(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
        XXW1(1,JL,J)=XXW1(1,JL,J)+CFWPOS(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW1(2,JL,J)=XXW1(2,JL,J)+CFWPOS(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      ELSE IF(JL /= 1)THEN
        XXW(14,JL,J)=XXW(14,JL,J)+CFWPOS(14,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(13,JL,J)=XXW(13,JL,J)+CFWPOS(13,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(12,JL,J)=XXW(12,JL,J)+CFWPOS(12,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(11,JL,J)=XXW(11,JL,J)+CFWPOS(11,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(10,JL,J)=XXW(10,JL,J)+CFWPOS(10,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(9,JL,J)=XXW(9,JL,J)+CFWPOS(9,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(8,JL,J)=XXW(8,JL,J)+CFWPOS(8,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(7,JL,J)=XXW(7,JL,J)+CFWPOS(7,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(6,JL,J)=XXW(6,JL,J)+CFWPOS(6,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(5,JL,J)=XXW(5,JL,J)+CFWPOS(5,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(4,JL,J)=XXW(4,JL,J)+CFWPOS(4,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(3,JL,J)=XXW(3,JL,J)+CFWPOS(3,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(2,JL,J)=XXW(2,JL,J)+CFWPOS(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(1,JL,J)=XXW(1,JL,J)+CFWPOS(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      END IF
    ENDDO
ENDDO

ZOPDP(:,:)=0.

!-----1.4  Water vapour:positive transmittances optically thick regime--

DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPWP2(J,JL)
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPA(J)
 ENDDO
ENDDO


DO J=1,KNCHPF
  ICH=KCHAN(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
        XXW1(1,JL,J)=XXW1(1,JL,J)+CFWPOS2(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW1(2,JL,J)=XXW1(2,JL,J)+CFWPOS2(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      ELSE IF(JL /= 1)THEN
        XXW(14,JL,J)=XXW(14,JL,J)+CFWPOS2(14,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(13,JL,J)=XXW(13,JL,J)+CFWPOS2(13,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(12,JL,J)=XXW(12,JL,J)+CFWPOS2(12,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(11,JL,J)=XXW(11,JL,J)+CFWPOS2(11,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(10,JL,J)=XXW(10,JL,J)+CFWPOS2(10,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(9,JL,J)=XXW(9,JL,J)+CFWPOS2(9,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(8,JL,J)=XXW(8,JL,J)+CFWPOS2(8,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(7,JL,J)=XXW(7,JL,J)+CFWPOS2(7,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(6,JL,J)=XXW(6,JL,J)+CFWPOS2(6,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(5,JL,J)=XXW(5,JL,J)+CFWPOS2(5,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(4,JL,J)=XXW(4,JL,J)+CFWPOS2(4,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(3,JL,J)=XXW(3,JL,J)+CFWPOS2(3,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(2,JL,J)=XXW(2,JL,J)+CFWPOS2(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW(1,JL,J)=XXW(1,JL,J)+CFWPOS2(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      END IF
    ENDDO
ENDDO

ZOPDP(:,:)=0.

!-----1.3  Water vapour:positive transmittances optically thin regime---


DO JL=NLEVW,JPLEV
  DO J=1,KNCHPF
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPWP1(J,JL)
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPA(J)
  ENDDO
ENDDO


DO J=1,KNCHPF
  ICH=KCHAN(J)
    DO JL=NLEVW,JPLEV
      IF(JL == 1)THEN
      XXW1(1,JL,J)=XXW1(1,JL,J)+CFWPOS1(1,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
        XXW1(2,JL,J)=XXW1(2,JL,J)+CFWPOS1(2,JL,ICH,KSAT)* &
        ZOPDP(J,JL)
      ELSE IF(JL /= 1)THEN
        XXW(1,JL,J)=XXW(1,JL,J)+CFWPOS1(1,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(2,JL,J)=XXW(2,JL,J)+CFWPOS1(2,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(3,JL,J)=XXW(3,JL,J)+CFWPOS1(3,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(4,JL,J)=XXW(4,JL,J)+CFWPOS1(4,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(5,JL,J)=XXW(5,JL,J)+CFWPOS1(5,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(6,JL,J)=XXW(6,JL,J)+CFWPOS1(6,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(7,JL,J)=XXW(7,JL,J)+CFWPOS1(7,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(8,JL,J)=XXW(8,JL,J)+CFWPOS1(8,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(9,JL,J)=XXW(9,JL,J)+CFWPOS1(9,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(10,JL,J)=XXW(10,JL,J)+CFWPOS1(10,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(11,JL,J)=XXW(11,JL,J)+CFWPOS1(11,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(12,JL,J)=XXW(12,JL,J)+CFWPOS1(12,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(13,JL,J)=XXW(13,JL,J)+CFWPOS1(13,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
        XXW(14,JL,J)=XXW(14,JL,J)+CFWPOS1(14,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      END IF
    ENDDO
ENDDO

ZOPDP(:,:)=0.


!-----1.2  Uniformly mixed gases:positive transmittances----------------

DO JL=1,JPLEV
  DO J=1,KNCHPF
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPMP(J,JL)
    ZOPDP(J,JL)=ZOPDP(J,JL)+OPDPA(J)
  ENDDO
ENDDO


DO J=1,KNCHPF
  ICH=KCHAN(J)
    DO JL=1,JPLEV
      XXM(1,JL,J)=XXM(1,JL,J)+CFMPOS(1,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(2,JL,J)=XXM(2,JL,J)+CFMPOS(2,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(3,JL,J)=XXM(3,JL,J)+CFMPOS(3,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(4,JL,J)=XXM(4,JL,J)+CFMPOS(4,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(5,JL,J)=XXM(5,JL,J)+CFMPOS(5,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(6,JL,J)=XXM(6,JL,J)+CFMPOS(6,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(7,JL,J)=XXM(7,JL,J)+CFMPOS(7,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(8,JL,J)=XXM(8,JL,J)+CFMPOS(8,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(9,JL,J)=XXM(9,JL,J)+CFMPOS(9,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
      XXM(10,JL,J)=XXM(10,JL,J)+CFMPOS(10,JL,ICH,KSAT)* &
       ZOPDP(J,JL)
    ENDDO
ENDDO

RETURN
END SUBROUTINE RTIASI_OPDEPK
!     Calculate transmittance on levels of radiative transfer model.

SUBROUTINE RTIASI_RTTAU(KCHAN,KPROF,KNCHPF,KSAT,TAU,TAUSFC,PWWR, &
     OPDPMP,OPDPO,OPDPWP1,OPDPWP2,OPDPWP,&
     TAUM,TAUW,TAUO)

!     Description:
!     Interpolate optical depths on to levels of *rt* model
!     (which, at present, entails only surface transmittance, as
!     other optical depths are on *rt* levels) and to convert
!     optical depths to transmittances.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     3) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Moved to F90. Modified for IASI. 
!                              Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_IASCHN,ONLY : &
!     Imported arrays:
       GAMMA      ! "Gamma factor" transmittance corrections


USE RTIASI_SURF , ONLY : &
!     Imported parameters:
       JPPF   , &   ! Max no. profiles
!     Imported arrays:
       NLEVSF , &   ! Index of nearest standard pressure level
!                   ! at/below surface.
       FRACPS       ! Fraction of standard pressure level interval by
                    ! which surf is above level nlevsf.


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
       JPLEV      ! Number of pressure levels


  IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
  INTEGER, INTENT(IN) :: KNCHPF                  ! Number of processed 
                                                 ! radiances.
  INTEGER, INTENT(IN) :: KSAT                    ! Input satellite index.
  INTEGER, INTENT(IN) :: KPROF   (KNCHPF)        ! Array of prof. indices
  INTEGER, INTENT(IN) :: KCHAN   (KNCHPF)        ! Array of ch. indices

!      Array arguments with intent in:
  REAL,INTENT(IN)     :: OPDPMP  (KNCHPF,JPLEV)  ! Mixed gases optical
                                                 ! depth.
  REAL,INTENT(IN)     :: OPDPO   (KNCHPF,JPLEV)  ! Ozone optical depth.
  REAL,INTENT(IN)     :: OPDPWP1 (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                 ! :weak absorption.
  REAL,INTENT(IN)     :: OPDPWP2 (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                 ! :strong absorption.
  REAL,INTENT(IN)     :: OPDPWP  (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                 ! :single reg. scheme
  REAL,INTENT(IN)     :: PWWR    (JPLEV,JPPF)    ! Pressure weighted 
                                                 ! water vapour ratio.

!       Array argumetnts with intent out:
  REAL,INTENT(OUT)    :: TAU     (KNCHPF,JPLEV)  ! Total transmittance.
  REAL,INTENT(OUT)    :: TAUSFC  (KNCHPF)        ! Surface transmittance.
  REAL,INTENT(OUT)    :: TAUM    (KNCHPF,JPLEV)  ! Fixed gases tran.s
  REAL,INTENT(OUT)    :: TAUW    (KNCHPF,JPLEV)  ! Water vapour trans.
  REAL,INTENT(OUT)    :: TAUO    (KNCHPF,JPLEV)  ! Ozone transmittance.

!     End of subroutine arguments



!       Local scalars:
  INTEGER :: JL,ICH,IPF,ISF,J
  REAL    :: TAULST
  REAL    :: ZZ
  REAL    :: ZA

!       Local arrays:
  REAL    :: OPDP(KNCHPF,JPLEV)


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   INTERPOLATE OPTICAL DEPTHS TO *RT* MODEL LEVELS.
!-----------------------------------------------------------------------
!              Not required (yet), except for ...



!-----------------------------------------------------------------------
!         2.   CONVERT OPTICAL DEPTHS TO TRANSMITTANCES.
!-----------------------------------------------------------------------

  DO J=1,KNCHPF

!-------2.1 Mixed gases-------------------------------------------------

     TAULST=1.
     DO JL=1,JPLEV
        IF(JL == 1)THEN
           
!-------------Level 1 transmittance is recovered------------------------

           OPDP(J,JL)=OPDPMP(J,JL)
           OPDP(J,JL)=MIN(LOG(1.06),OPDP(J,JL))
           TAUM(J,JL)=TAULST*EXP(OPDP(J,JL))-0.06
        ELSE IF(JL /= 1)THEN

!-------------This condition means that we are in a region of positive--
!             transmittances and the appropriate coefficients should
!             be used.

           IF(-ALOG(TAUM(J,JL-1)+0.06) < 2.8)THEN

!---------------The layer transmittance is recovered using the----------
!               fictiuous optical depth.

              OPDP(J,JL)=OPDPMP(J,JL)

!---------------If the layer optical depth is equal to zero, the--------
!               transmittance is set to a small negative value.

              IF (OPDP(J,JL) == 0..AND.TAUM(J,JL-1) < 0.06) THEN
                 OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
              END IF

!-----------------Check for reasonable value----------------------------

              IF(TAULST > 0.OR.TAULST < -0.06)THEN
                 OPDP(J,JL)=MIN(OPDP(J,JL),0.)
              ELSE IF(TAULST > -0.06.AND.TAULST < 0.)THEN
                 OPDP(J,JL)=MAX(0.,OPDP(J,JL))
                 IF (OPDP(J,JL) >= LOG(0.06/(TAULST+0.06))) THEN
                    OPDP(J,JL)=0.
                 END IF
              END IF

!---------------If the layer optical depth is equal to zero, the--------
!               transmittance is set to a small negative value.

              IF (OPDP(J,JL) == 0..AND.TAUM(J,JL-1) < 0.06) THEN
                 OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
              END IF

!-------------This condition means that we are in a region of negative--
!             transmittances:the transmittance is set equal to -0.0003--

           ELSE IF(-ALOG(TAUM(J,JL-1)+0.06) > 2.8)THEN
              OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
           END IF
           ZZ = EXP(OPDP(J,JL)) + (EXP(OPDP(J,JL))*0.06/TAULST) - &
                (0.06/TAULST)
           TAUM(J,JL)=TAULST*ZZ
        END IF
        IF (TAUM(J,JL) <= 1.0E-05) THEN
           TAUM(J,JL)=-0.0003
        END IF
        TAULST=TAUM(J,JL)
     ENDDO                                      ! End of level loop


!-------2.2 Water vapour------------------------------------------------

     TAULST=1.
     DO JL=1,JPLEV
        IF(JL == 1)THEN
           OPDP(J,JL) = OPDPWP1(J,JL)
           OPDP(J,JL) =MIN(LOG(1.06),OPDP(J,JL))
           TAUW(J,JL)=(TAULST*EXP(OPDP(J,JL)))-0.06
        ELSE IF(JL /= 1)THEN

!-------------This condition means that we are in a region of positive -
!             transmittances and the appropriate coefficients should be-
!             used.

           IF(-ALOG(TAUW(J,JL-1)+0.06) < 2.775)THEN

!---------------Appropriate coefficients are used when this condition---
!               is fullfilled. 

              IF(-ALOG(TAUW(J,JL-1)+0.06)*PWWR(JL,KPROF(J)) < 0.8)THEN
                 OPDP(J,JL)=OPDPWP1(J,JL)

!---------------Appropriate coefficients are used when this condition is 
!               fullfilled:strong water vapour absorption.

              ELSE IF(-ALOG(TAUW(J,JL-1)+0.06)*PWWR(JL,KPROF(J))> 2) &
                   THEN
                 OPDP(J,JL)=OPDPWP2(J,JL)

              ELSE IF (-ALOG(TAUW(J,JL-1)+0.06)*PWWR(JL,KPROF(J)) &
                   >= 0.8 .OR. &
                   -ALOG(TAUW(J,JL-1)+0.06)*PWWR(JL,KPROF(J)) <= 2) THEN
                 OPDP(J,JL)=OPDPWP(J,JL)

              END IF

!-------------If the layer optical depth is equal to zero, the ---------
!             transmittance is set to a small negative value (if the
!             optical depth is equal to zero  negative tranmsittances  
!             should be encountered).

              IF (OPDP(J,JL) == 0..AND.TAUW(J,JL-1) < 0.06) THEN
                 OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
              END IF

!---------------Check for reasonable value------------------------------

              IF(TAULST > 0.OR.TAULST < -0.06)THEN
                 OPDP(J,JL)=MIN(OPDP(J,JL),0.)
              ELSE IF(TAULST > -0.06.AND.TAULST < 0.)THEN
                 OPDP(J,JL)=MAX(0.,OPDP(J,JL))
                 IF (OPDP(J,JL) >= LOG(0.06/(TAULST+0.06))) THEN
                    OPDP(J,JL)=0.
                 END IF
              END IF

!-------------If the layer optical depth is equal to zero, the----------
!             transmittance is set to a small negative value 

              IF (OPDP(J,JL) == 0..AND.TAUW(J,JL-1) < 0.06) THEN
                 OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
              END IF

!-------------This condition means that we are in a region of negative--
!             transmittances :the transmittance is set equal to -0.0003-

           ELSE IF(-ALOG(TAUW(J,JL-1)+0.06) > 2.775)THEN
              OPDP(J,JL)=ALOG((-0.0003+0.06)/(TAULST+0.06))
           END IF
           ZZ = EXP(OPDP(J,JL)) + (EXP(OPDP(J,JL))*0.06/TAULST) - &
                (0.06/TAULST)
           TAUW(J,JL)=TAULST*ZZ
        END IF
        IF (TAUW(J,JL) <= 1.0E-05) THEN
           TAUW(J,JL)=-0.0003
        END IF
        TAULST=TAUW(J,JL)
     ENDDO                                       !End of level loop 


!-------2.3 Ozone-------------------------------------------------------

     TAULST=1.
     DO JL=1,JPLEV
        IF(JL == 1)THEN
           OPDP(J,JL)=OPDPO(J,JL)
           OPDP(J,JL)=MIN(LOG(1.06),OPDP(J,JL))
           TAUO(J,JL)=TAULST*EXP(OPDP(J,JL))-0.06
        ELSE IF(JL /= 1)THEN
           OPDP(J,JL)=OPDPO(J,JL)

!---------------Check for reasonable value------------------------------

           IF(TAULST > 0.OR.TAULST < -0.06)THEN
              OPDP(J,JL)=MIN(OPDP(J,JL),0.)
           ELSE IF(TAULST > -0.06.AND.TAULST < 0.)THEN
              OPDP(J,JL)=MAX(0.,OPDP(J,JL))
              IF (OPDP(J,JL) >= LOG(0.06/(TAULST+0.06))) THEN
                 OPDP(J,JL)=0.
              END IF
           END IF
!-----------------------------------------------------------------------

           ZZ = EXP(OPDP(J,JL)) + (EXP(OPDP(J,JL))*0.06/TAULST) - &
                (0.06/TAULST)
           TAUO(J,JL)=TAULST*ZZ
        END IF
        IF (TAUO(J,JL) <= 1.0E-05) THEN
           TAUO(J,JL)=-0.0003
        END IF
        TAULST=TAUO(J,JL)
     ENDDO                                       !End of level loop

!-----------------------------------------------------------------------

     DO JL=1,JPLEV
        TAU(J,JL)=TAUM(J,JL)*TAUW(J,JL)*TAUO(J,JL)
!         TAU(J,JL)=TAUM(J,JL)*TAUW(J,JL)
!         TAU(J,JL)=TAUW(J,JL)
     ENDDO
  END DO                                !End of channel/profile loop

!-----End of section 2--------------------------------------------------



!-----------------------------------------------------------------------
!         3.   PRINT TRANSMITTANCES FOR SELECTED CHANNEL
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
!         4.   APPLY GAMMA CORRECTIONS 
!-----------------------------------------------------------------------

  DO J=1,KNCHPF
     ICH=KCHAN(J)
     DO JL=1,JPLEV
        IF (TAU(J,JL) > 0) THEN
           TAU(J,JL)=TAU(J,JL)**GAMMA(ICH,KSAT)
        END IF
        IF (TAU(J,JL) < 0) THEN
           TAU(J,JL)=-(ABS(TAU(J,JL))**GAMMA(ICH,KSAT))
        END IF
     END DO
  END DO

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         5.   CALCULATE TRANSMITTANCE AT SURFACE.
!-----------------------------------------------------------------------
  DO J=1,KNCHPF
     IPF=KPROF(J)
     ISF=NLEVSF(IPF)
     ZA=-ALOG(TAU(J,ISF)+0.06)+FRACPS(IPF) * &
          (-ALOG(TAU(J,ISF-1)+0.06)+ALOG(TAU(J,ISF)+0.06))
     TAUSFC(J)=EXP(-ZA)-0.06
  ENDDO



  RETURN
END SUBROUTINE RTIASI_RTTAU
!     Calculate transmittance on levels of radiative transfer model.
SUBROUTINE RTIASI_RTTAUK(KCHAN,KPROF,KNCHPF,KSAT,TAU_D,TAU, &
     TAUSFC,PWWR_D,OPDPMP_D,OPDPMP,OPDPO_D,OPDPO,OPDPWP1_D, &
     OPDPWP1,OPDPWP2_D,OPDPWP2,OPDPWP_D,OPDPWP, &
     TAUM_D,TAUM,TAUW_D,TAUW,TAUO_D,TAUO)

!     Description:
!     K of subroutine RTTAU
!     Interpolate optical depths on to levels of *rt* model
!     (which, at present, entails only surface transmittance, as
!     other optical depths are on *rt* levels) and to convert
!     optical depths to transmittances.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     3) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Moved to F90. Modified for IASI.
!                              Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct module used:

  USE RTIASI_IASCHN,ONLY : &
!     Imported arrays:
       GAMMA                ! "Gamma factor" transmittance corrections


  USE RTIASI_SURF, ONLY : &
!     Imported parameters:
       JPPF            , &   ! Max no. profiles
!     Imported arrays:
       NLEVSF_D =>NLEVSF, &  ! Index of nearest standard pressure level
                             ! at/below surface.
       FRACPS_D =>FRACPS     ! Fraction of standard pressure level interval by
                             ! which surf is above level nlevsf.


  USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
       JPLEV            ! Number of pressure levels


!     Tangenl linear module used:

  USE RTIASI_SURFK, ONLY : &
!     Imported arrays:
       FRACPS               ! Fraction of standard pressure level interval by
                            ! which surf is above level nlevsf.


  IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
  INTEGER, INTENT(IN) :: KNCHPF                    ! Number of processed 
                                                   ! radiances*profiles.
  INTEGER, INTENT(IN) :: KSAT                      ! Input satellite index.

!       Direct array arguments with intent in:
  INTEGER, INTENT(IN) :: KPROF      (KNCHPF)       ! Array of prof. indices
  INTEGER, INTENT(IN) :: KCHAN      (KNCHPF)       ! Array of ch. indices
  REAL,INTENT(IN)     :: OPDPMP_D   (KNCHPF,JPLEV) ! Mixed gases optical
                                                   ! depth.
  REAL,INTENT(IN)     :: OPDPO_D   (KNCHPF,JPLEV)  ! Ozone optical depth.
  REAL,INTENT(IN)     :: OPDPWP1_D (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :weak absorption.
  REAL,INTENT(IN)     :: OPDPWP2_D (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :strong absorption.
  REAL,INTENT(IN)     :: OPDPWP_D  (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :single reg. scheme
  REAL,INTENT(IN)     :: PWWR_D    (JPLEV,JPPF)    ! Pressure weighted 
                                                   ! water vapour ratio.
  REAL,INTENT(IN)     :: TAU_D     (KNCHPF,JPLEV)  ! Total transmittance.
  REAL,INTENT(IN)     :: TAUM_D    (KNCHPF,JPLEV)  ! Fixed gases trans.
  REAL,INTENT(IN)     :: TAUW_D    (KNCHPF,JPLEV)  ! Water vapour trans.
  REAL,INTENT(IN)     :: TAUO_D    (KNCHPF,JPLEV)  ! Ozone transmittance.

!       K array arguments with intent out:
  REAL,INTENT(INOUT)  :: OPDPMP    (KNCHPF,JPLEV)  ! Mixed gases optical
                                                   ! depth.
  REAL,INTENT(INOUT)  :: OPDPO     (KNCHPF,JPLEV)  ! Ozone optical depth.
  REAL,INTENT(INOUT)  :: OPDPWP1   (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :weak absorption.
  REAL,INTENT(INOUT)  :: OPDPWP2   (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :strong absorption.
  REAL,INTENT(INOUT)  :: OPDPWP    (KNCHPF,JPLEV)  ! Wa. vap. optical depth
                                                   ! :single reg. scheme

!       K array argumetnts with intent in:
  REAL,INTENT(INOUT)  :: TAU       (KNCHPF,JPLEV)  ! Total transmittance.
  REAL,INTENT(INOUT)  :: TAUSFC    (KNCHPF)        ! Surface transmittance.
  REAL,INTENT(INOUT)  :: TAUM      (KNCHPF,JPLEV)  ! Fixed gases trans.
  REAL,INTENT(INOUT)  :: TAUW      (KNCHPF,JPLEV)  ! Water vapour trans.
  REAL,INTENT(INOUT)  :: TAUO      (KNCHPF,JPLEV)  ! Ozone transmittance

!     End of subroutine arguments



!       Local scalars:
  INTEGER :: JL,ICH,IPF,ISF,J
  REAL    :: ZA,ZA_D
  
!       Local K scalars:
  REAL    :: TAULSTM,TAULSTW,TAULSTO
  REAL    :: ZZ

!       Local arrays:
  REAL    :: OPDPM_D   (KNCHPF,JPLEV)
  REAL    :: OPDPW_D   (KNCHPF,JPLEV)
  REAL    :: OPDPOZ_D  (KNCHPF,JPLEV)
  REAL    :: OPDPM1_D   (JPLEV)
  REAL    :: OPDPW1_D   (JPLEV)
  REAL    :: OPDPOZ1_D  (JPLEV)
  REAL    :: OPDPM2_D   (JPLEV)
  REAL    :: OPDPW2_D   (JPLEV)
  REAL    :: OPDPM3_D   (JPLEV)
  REAL    :: OPDPW3_D   (JPLEV)
  REAL    :: TAULSTM_D (KNCHPF,JPLEV+1)
  REAL    :: TAULSTW_D (KNCHPF,JPLEV+1)
  REAL    :: TAULSTO_D (KNCHPF,JPLEV+1)
  REAL    :: ZZM_D     (KNCHPF,JPLEV)
  REAL    :: ZZW_D     (KNCHPF,JPLEV)
  REAL    :: ZZO_D     (KNCHPF,JPLEV)

!       Local K arrays:
  REAL    :: OPDPM     (KNCHPF,JPLEV)
  REAL    :: OPDPW     (KNCHPF,JPLEV)
  REAL    :: OPDPOZ    (KNCHPF,JPLEV)



!-----End of header-----------------------------------------------------



!-----Initialize variables----------------------------------------------

  ZA          =0.
  ZZ          =0.
  TAULSTO     =0.
  TAULSTW     =0.
  TAULSTM     =0.
  OPDPOZ(:,:) =0
  OPDPW(:,:)  =0
  OPDPM(:,:)  =0
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         1.   INTERPOLATE OPTICAL DEPTHS TO *RT* MODEL LEVELS.
!-----------------------------------------------------------------------
!              Not required (yet), except for ...



!-----------------------------------------------------------------------
!         4.   CALCULATE TRANSMITTANCE AT SURFACE.
!-----------------------------------------------------------------------
  DO J=1,KNCHPF
     IPF=KPROF(J)
     ISF=NLEVSF_D(IPF)
     
     ZA_D=-ALOG(TAU_D(J,ISF)+0.06)+FRACPS_D(IPF)* &
          (-ALOG(TAU_D(J,ISF-1)+0.06)+ALOG(TAU_D(J,ISF)+0.06))

     ZA=ZA-TAUSFC(J)*EXP(-ZA_D)
     TAU(J,ISF)=TAU(J,ISF)-ZA/(TAU_D(J,ISF)+0.06)
     FRACPS(J)=FRACPS(J)+ZA * &
          (-ALOG(TAU_D(J,ISF-1)+0.06)+ALOG(TAU_D(J,ISF)+0.06))
     TAU(J,ISF-1)=TAU(J,ISF-1)- FRACPS_D(IPF)*ZA/ &
          (TAU_D(J,ISF-1)+0.06)
     TAU(J,ISF)=TAU(J,ISF)+FRACPS_D(IPF)*ZA/(TAU_D(J,ISF)+0.06)
     ZA=0.

!-----------------------------------------------------------------------
!         3.   APPLY GAMMA CORRECTIONS
!-----------------------------------------------------------------------

     ICH=KCHAN(J)
     DO JL=1,JPLEV
        IF (TAU_D(J,JL) > 0) THEN
           TAU(J,JL)=GAMMA(ICH,KSAT) * &
                TAU_D(J,JL)**(GAMMA(ICH,KSAT)-1)*TAU(J,JL)
        END IF
        IF (TAU_D(J,JL) < 0) THEN
           TAU(J,JL)=+GAMMA(ICH,KSAT) * &
                ABS(TAU_D(J,JL))**(GAMMA(ICH,KSAT)-1)*TAU(J,JL)
        END IF
     END DO

!-----------------------------------------------------------------------




!-----REPEAT DIRECT CALCULATIONS----------------------------------------


!-----Ozone-------------------------------------------------------------

     DO JL=1,JPLEV
        IF(JL == 1)THEN
           OPDPOZ_D(J,JL)=OPDPO_D(J,JL)
           TAULSTO_D(J,JL)=1
        ELSE IF(JL /= 1)THEN
           OPDPOZ_D(J,JL)=OPDPO_D(J,JL)
           OPDPOZ1_D(JL)=OPDPOZ_D(J,JL)
           IF(TAULSTO_D(J,JL) > 0.OR.TAULSTO_D(J,JL) < -0.06)THEN
              OPDPOZ_D(J,JL)=MIN(OPDPOZ_D(J,JL),0.)
           ELSE IF(TAULSTO_D(J,JL) > -0.06 .AND. &
                TAULSTO_D(J,JL) < 0.) THEN
              OPDPOZ_D(J,JL)=MAX(0.,OPDPOZ_D(J,JL))
              IF (OPDPOZ_D(J,JL) >= &
                   LOG(0.06/(TAULSTO_D(J,JL)+0.06))) THEN
                 OPDPOZ_D(J,JL)=0.
              END IF
           ENDIF
           ZZO_D(J,JL)=EXP(OPDPOZ_D(J,JL)) + &
                (EXP(OPDPOZ_D(J,JL))*0.06/TAULSTO_D(J,JL))- &
                (0.06/TAULSTO_D(J,JL))
        END IF
        TAULSTO_D(J,JL+1)=TAUO_D(J,JL)
     ENDDO
     
!-----------------------------------------------------------------------

!-----Water Vapour------------------------------------------------------

     DO JL=1,JPLEV
        IF(JL == 1)THEN
           OPDPW_D(J,JL) = OPDPWP1_D(J,JL)
           TAULSTW_D(J,JL)=1.
        ELSE IF(JL /= 1)THEN
           IF(-ALOG(TAUW_D(J,JL-1)+0.06) < 2.775)THEN
              IF(-ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J)) < 0.8)THEN
                 OPDPW_D(J,JL) = OPDPWP1_D(J,JL)
              ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J))> 2) THEN
                 OPDPW_D(J,JL) = OPDPWP2_D(J,JL)
              ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J))>= 0.8 .OR. &
                   -ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J)) <=2)THEN
                 OPDPW_D(J,JL) = OPDPWP_D(J,JL)
              END IF
              
              OPDPW1_D(JL)=OPDPW_D(J,JL)
              IF (OPDPW_D(J,JL) == 0..AND.TAUW_D(J,JL-1) < 0.06) THEN
                 OPDPW_D(J,JL)=ALOG((-0.0003+0.06)/ &
                      (TAULSTW_D(J,JL)+0.06))
              END IF
              
              OPDPW2_D(JL)=OPDPW_D(J,JL)
              IF(TAULSTW_D(J,JL) > 0.OR.TAULSTW_D(J,JL) < -0.06)THEN
                 OPDPW_D(J,JL)=MIN(OPDPW_D(J,JL),0.)
              ELSE IF(TAULSTW_D(J,JL) > -0.06 .AND. &
                   TAULSTW_D(J,JL) < 0.)THEN
                 OPDPW_D(J,JL)=MAX(0.,OPDPW_D(J,JL))
                 IF (OPDPW_D(J,JL) >= &
                      LOG(0.06/(TAULSTW_D(J,JL)+0.06))) THEN
                    OPDPW_D(J,JL)=0.
                 END IF
              END IF
              
              OPDPW3_D(JL)=OPDPW_D(J,JL)
              IF (OPDPW_D(J,JL) == 0..AND.TAUW_D(J,JL-1) < 0.06) THEN
                 OPDPW_D(J,JL) = &
                      ALOG((-0.0003+0.06)/(TAULSTW_D(J,JL)+0.06))
              END IF
              
           ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06) > 2.775)THEN
              OPDPW_D(J,JL) = &
                   ALOG((-0.0003+0.06)/(TAULSTW_D(J,JL)+0.06))
           END IF
           ZZW_D(J,JL)=EXP(OPDPW_D(J,JL))+(EXP(OPDPW_D(J,JL))*0.06/ &
                TAULSTW_D(J,JL))-(0.06/TAULSTW_D(J,JL))
        END IF
        TAULSTW_D(J,JL+1)=TAUW_D(J,JL)
     ENDDO
     
!-----------------------------------------------------------------------

!-----Fixed gases-------------------------------------------------------

     DO JL=1,JPLEV
        IF(JL == 1)THEN
           OPDPM_D(J,JL)  = OPDPMP_D(J,JL)
           TAULSTM_D(J,JL)=1.
        ELSE IF(JL /= 1)THEN
           IF(-ALOG(TAUM_D(J,JL-1)+0.06) < 2.8)THEN
              OPDPM_D(J,JL)  = OPDPMP_D(J,JL)
              
              OPDPM1_D(JL)=OPDPW_D(J,JL)
              IF (OPDPM_D(J,JL) == 0. .AND. &
                   TAUM_D(J,JL-1) < 0.06)THEN
                 OPDPM_D(J,JL)=ALOG((-0.0003+0.06)/ &
                      (TAULSTM_D(J,JL)+0.06))
              ENDIF
              
              OPDPM2_D(JL)=OPDPW_D(J,JL)
              IF(TAULSTM_D(J,JL) > 0.OR.TAULSTM_D(J,JL) < -0.06)THEN
                 OPDPM_D(J,JL)=MIN(OPDPM_D(J,JL),0.)
              ELSE IF(TAULSTM_D(J,JL) > -0.06 .AND. &
                   TAULSTM_D(J,JL) < 0.)THEN
                 OPDPM_D(J,JL)=MAX(0.,OPDPM_D(J,JL))
                 IF (OPDPM_D(J,JL) >= &
                      LOG(0.06/(TAULSTM_D(J,JL)+0.06))) THEN
                    OPDPM_D(J,JL)=0.
                 END IF
              END IF
              
              OPDPM3_D(JL)=OPDPW_D(J,JL)
              IF (OPDPM_D(J,JL) == 0..AND.TAUM_D(J,JL-1) < 0.06) THEN
                 OPDPM_D(J,JL) = &
                      ALOG((-0.0003+0.06)/(TAULSTM_D(J,JL)+0.06))
              ENDIF
              
           ELSE IF(-ALOG(TAUM_D(J,JL-1)+0.06) > 2.8)THEN
              OPDPM_D(J,JL) = &
                   ALOG((-0.0003+0.06)/(TAULSTM_D(J,JL)+0.06))
           ENDIF
           ZZM_D(J,JL)=EXP(OPDPM_D(J,JL))+(EXP(OPDPM_D(J,JL))*0.06/ &
                TAULSTM_D(J,JL))-(0.06/TAULSTM_D(J,JL))
        ENDIF
        TAULSTM_D(J,JL+1)=TAUM_D(J,JL)
     ENDDO

     DO JL=1,JPLEV
        TAUM(J,JL)=TAUM(J,JL)+TAU(J,JL)*TAUW_D(J,JL)*TAUO_D(J,JL)
        TAUW(J,JL)=TAUW(J,JL)+TAU(J,JL)*TAUM_D(J,JL)*TAUO_D(J,JL)
        TAUO(J,JL)=TAUO(J,JL)+TAU(J,JL)*TAUM_D(J,JL)*TAUW_D(J,JL)
     ENDDO



!-----K CALCULATIONS----------------------------------------------------

!-----Ozone-------------------------------------------------------------

     DO JL=JPLEV,1,-1
        TAUO(J,JL)=TAUO(J,JL)+TAULSTO
        TAULSTO=0
        IF (TAUO_D(J,JL) <= 1.0E-05) THEN
           TAUO(J,JL)=0.
        END IF
        IF(JL == 1)THEN
           OPDPOZ(J,JL)=OPDPOZ(J,JL)+EXP(OPDPOZ_D(J,JL))*TAUO(J,JL)
           IF(OPDPOZ_D(J,JL) > LOG(1.06) ) THEN
              OPDPOZ(J,JL)=0.
           ENDIF
           OPDPO(J,JL)=OPDPO(J,JL)+OPDPOZ(J,JL)
        ELSE IF(JL /= 1)THEN
           
           TAULSTO=TAULSTO+TAUO(J,JL)*ZZO_D(J,JL)
           ZZ=ZZ+TAUO(J,JL)*TAULSTO_D(J,JL)
           
           OPDPOZ(J,JL)=OPDPOZ(J,JL)+EXP(OPDPOZ_D(J,JL))*ZZ
           TAULSTO=TAULSTO+0.06*ZZ/TAULSTO_D(J,JL)**2
           OPDPOZ(J,JL)=OPDPOZ(J,JL)+ &
                EXP(OPDPOZ_D(J,JL))*ZZ*0.06/TAULSTO_D(J,JL)
           TAULSTO=TAULSTO- &
                EXP(OPDPOZ_D(J,JL))*0.06*ZZ/TAULSTO_D(J,JL)**2
           ZZ=0
           
           IF(TAULSTO_D(J,JL) > 0.OR.TAULSTO_D(J,JL) < -0.06)THEN
              IF ( OPDPOZ1_D(JL) > 0. ) THEN
!xxx                OPDPOZ(J,JL-1)=OPDPOZ(J,JL-1)+ OPDPOZ(J,JL)
                 OPDPOZ(J,JL) = 0.
              END IF
           ELSE IF(TAULSTO_D(J,JL) > -0.06 .AND. &
                TAULSTO_D(J,JL) < 0.) THEN
              IF (OPDPOZ1_D(JL) >= &
                   LOG(0.06/(TAULSTO_D(J,JL)+0.06))) THEN
!xxx                OPDPOZ(J,JL-1)=OPDPOZ(J,JL-1)+OPDPOZ(J,JL)
                 OPDPOZ(J,JL) = 0.
              END IF
              IF ( OPDPOZ1_D(JL) < 0. ) THEN
!xxx                 OPDPOZ(J,JL-1)=OPDPOZ(J,JL-1)+OPDPOZ(J,JL)
                 OPDPOZ(J,JL) = 0.
              END IF
           ENDIF
           OPDPO(J,JL)=OPDPO(J,JL)+OPDPOZ(J,JL)
        ENDIF
     ENDDO

!-----WATER VAPOUR----------------------------------------------------

     DO JL=JPLEV,1,-1
        TAUW(J,JL)=TAUW(J,JL)+TAULSTW
        TAULSTW=0.
        IF (TAUW_D(J,JL) <= 1.0E-05) THEN
           TAUW(J,JL)=0.
        END IF
        IF(JL == 1)THEN
           OPDPW(J,JL)=OPDPW(J,JL)+TAUW(J,JL)*EXP(OPDPW_D(J,JL))
           IF (OPDPW_D(J,JL) > LOG(1.06) ) THEN
              OPDPW(J,JL)=0.
           ENDIF
           OPDPWP1(J,JL)=OPDPWP1(J,JL)+OPDPW(J,JL)
        ELSE IF(JL /= 1)THEN
           
           ZZ=ZZ+TAUW(J,JL)*TAULSTW_D(J,JL)
           TAULSTW=TAULSTW+TAUW(J,JL)*ZZW_D(J,JL)
           
           OPDPW(J,JL)=OPDPW(J,JL)+ZZ*EXP(OPDPW_D(J,JL))
           TAULSTW=TAULSTW+0.06*ZZ/TAULSTW_D(J,JL)**2
           OPDPW(J,JL)=OPDPW(J,JL)+EXP(OPDPW_D(J,JL))*ZZ*0.06/ &
                TAULSTW_D(J,JL)
           TAULSTW=TAULSTW-EXP(OPDPW_D(J,JL))*0.06* &
                ZZ/TAULSTW_D(J,JL)**2
           ZZ=0
           
           IF(-ALOG(TAUW_D(J,JL-1)+0.06) < 2.775)THEN

!-------------If the layer optical depth is equal to zero, the----------
!             transmittance is set to a small negative value
              IF (OPDPW3_D(JL) == 0..AND.TAUW_D(J,JL-1) < 0.06) THEN
                 TAULSTW=TAULSTW-((TAULSTW_D(J,JL)+0.06)/ &
                      (-0.0003+0.06))*((-0.0003+0.06)/ &
                      (TAULSTW_D(J,JL)+0.06)**2)*OPDPW(J,JL)
              ENDIF

!---------------Check for reasonable value------------------------------
              IF(TAULSTW_D(J,JL) > 0.OR.TAULSTW_D(J,JL) < -0.06)THEN
                 IF ( OPDPW2_D(JL) > 0. ) THEN
!xxx              OPDPW(J,JL-1)=OPDPW(J,JL-1)+OPDPW(J,JL)
                    OPDPW(J,JL) = 0.
                 END IF
              ELSE IF(TAULSTW_D(J,JL) > -0.06 .AND. &
                   TAULSTW_D(J,JL) < 0.)THEN
                 IF (OPDPW2_D(JL) >= &
                      LOG(0.06/(TAULSTW_D(J,JL)+0.06))) THEN
!xxx            OPDPW(J,JL-1)=OPDPW(J,JL-1)+OPDPW(J,JL)
                    OPDPW(J,JL) = 0.
                 END IF
                 IF ( OPDPW2_D(JL) < 0. ) THEN
!xxx            OPDPW(J,JL-1)=OPDPW(J,JL-1)+OPDPW(J,JL)
                    OPDPW(J,JL) = 0.
                 ENDIF
              ENDIF
          
!-------------If the layer optical depth is equal to zero, the ---------
!             transmittance is set to a small negative value (if the
!             optical depth is equal to zero  negative tranmsittances
!             should be encountered).

              IF (OPDPW1_D(JL) == 0..AND.TAUW_D(J,JL-1) < 0.06) THEN
                 TAULSTW=TAULSTW-((TAULSTW_D(J,JL)+0.06)/ &
                      (-0.0003+0.06))*((-0.0003+0.06)/ &
                      (TAULSTW_D(J,JL)+0.06)**2)*OPDPW(J,JL)
              ENDIF

              IF(-ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J))>=0.8 .AND. &
                   -ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J))<=2)THEN
                 OPDPWP(J,JL)=OPDPWP(J,JL)+OPDPW  (J,JL)
              ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06)* &
                   PWWR_D(JL,KPROF(J)) < 2)THEN
                 OPDPWP1(J,JL)=OPDPWP1(J,JL)+OPDPW  (J,JL)
              ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06) * &
                   PWWR_D(JL,KPROF(J))> 0.8) THEN
                 OPDPWP2(J,JL)=OPDPWP2(J,JL)+OPDPW  (J,JL)
              ENDIF


           ELSE IF(-ALOG(TAUW_D(J,JL-1)+0.06) > 2.775)THEN
              TAULSTW=TAULSTW-((TAULSTW_D(J,JL)+0.06)/ &
                   (-0.0003+0.06))*((-0.0003+0.06)/ &
                   (TAULSTW_D(J,JL)+0.06)**2)*OPDPW(J,JL)
           ENDIF
              
        ENDIF
     ENDDO

!-----Fixed gases------------------------------------------------------

     DO JL=JPLEV,1,-1
        TAUM(J,JL)=TAUM(J,JL)+TAULSTM
        TAULSTM=0
        IF (TAUM_D(J,JL) <= 1.0E-05) THEN
           TAUM(J,JL)=0.
        END IF
        IF(JL == 1)THEN
           OPDPM(J,JL)=OPDPM(J,JL)+TAUM(J,JL)*EXP(OPDPM_D(J,JL))
           IF( OPDPM_D(J,JL) > LOG(1.06) ) THEN
              OPDPM(J,JL)=0.
           ENDIF
           OPDPMP(J,JL)=OPDPMP(J,JL)+OPDPM(J,JL)
        ELSE IF(JL /= 1)THEN

           ZZ=ZZ+TAUM(J,JL)*TAULSTM_D(J,JL)
           TAULSTM=TAULSTM+TAUM(J,JL)*ZZM_D(J,JL)
           
           OPDPM(J,JL)=OPDPM(J,JL)+ZZ*EXP(OPDPM_D(J,JL))
           TAULSTM=TAULSTM+0.06*ZZ/TAULSTM_D(J,JL)**2
           OPDPM(J,JL)=OPDPM(J,JL)+EXP(OPDPM_D(J,JL))*ZZ*0.06/ &
                TAULSTM_D(J,JL)
           TAULSTM=TAULSTM - EXP(OPDPM_D(J,JL))*0.06*ZZ/ &
                TAULSTM_D(J,JL)**2
           ZZ=0
           IF(-ALOG(TAUM_D(J,JL-1)+0.06) < 2.8)THEN
              IF (OPDPM3_D(JL) == 0..AND.TAUM_D(J,JL-1) < 0.06) THEN
                 TAULSTM=TAULSTM-((TAULSTM_D(J,JL)+0.06)/ &
                      (-0.0003+0.06))*((-0.0003+0.06)/ &
                      (TAULSTM_D(J,JL)+0.06)**2)*OPDPM(J,JL)
              ENDIF
              IF(TAULSTM_D(J,JL) > 0.OR.TAULSTM_D(J,JL) < -0.06)THEN
                 IF ( OPDPM2_D(JL) > 0 ) THEN
!xxx               OPDPM(J,JL-1)=OPDPM(J,JL-1)+OPDPM(J,JL)
                    OPDPM(J,JL) = 0.
                 END IF
              ELSE IF(TAULSTM_D(J,JL) > -0.06 .AND. &
                   TAULSTM_D(J,JL) < 0.)THEN
                 IF (OPDPM2_D(JL) >= &
                      LOG(0.06/(TAULSTM_D(J,JL)+0.06))) THEN
!xxx             OPDPM(J,JL-1)=OPDPM(J,JL-1)+OPDPM(J,JL)
                    OPDPM(J,JL) = 0.
                 END IF
                 IF (OPDPM2_D(JL) < 0.) THEN
!xxx           OPDPM(J,JL-1)=OPDPM(J,JL-1)+OPDPM(J,JL)
                    OPDPM(J,JL) = 0.
                 END IF
              END IF
              IF (OPDPM1_D(JL) == 0. .AND. &
                   TAUM_D(J,JL-1) < 0.06) THEN
                 TAULSTM=TAULSTM-((TAULSTM_D(J,JL)+0.06)/ &
                      (-0.0003+0.06))*((-0.0003+0.06)/ &
                      (TAULSTM_D(J,JL)+0.06)**2)*OPDPM(J,JL)
              ENDIF
              OPDPMP(J,JL)=OPDPMP(J,JL)+OPDPM(J,JL)
           ELSE IF(-ALOG(TAUM_D(J,JL-1)+0.06) > 2.8)THEN
              TAULSTM=TAULSTM-((TAULSTM_D(J,JL)+0.06)/ &
                   (-0.0003+0.06))*((-0.0003+0.06)/ &
                   (TAULSTM_D(J,JL)+0.06)**2)*OPDPM(J,JL)
           ENDIF
        END IF
        
     ENDDO
  ENDDO



  RETURN
END SUBROUTINE RTIASI_RTTAUK

!     Performs integration of radiative transfer equation in RTIASI.
SUBROUTINE RTIASI_RTINT(KCHAN,KPROF,KNCHPF,PRAD,PTB, &
EMS,TAU,TAUSFC,B,BCG,BA,BS,RADOV,RADO,BDT,BDTR,FCLD,RADCL,TBCL, &
PANGL,KNPF,KSAT,IPRINT)

!     Description
!     Perform integration of radiative transfer equation
!     in RTIASI suite, calculating radiances and brightness
!     temperatures.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2)See: User's manual for RTIASI (Available from EUMETSAT)
!     3)See J.R.Eyre: ECMWF TECH. MEMO 176
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            25/06/91.   Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Moved to F90. New definition of average layer
!                              temperature. Emissivities are fully accounted for
!                              using a fast model developed by M.Matricardi.
!                              ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used

USE RTIASI_IRCLD, ONLY : &
!     Imported arrays:
  NLEVCD ,   & ! Index of nearest standard pressure level
             ! at/below cloud top.
  FRACPC     ! Fraction of standardd pressure level interval
             ! by which cloud is above level NLEVCD.


USE RTIASI_PRFVAR,ONLY : &
!     Imported parameters:
  JPLEV  ,   & ! Number of pressure levels
!     Imported arrays:
  TEMPW  ,   & ! Curtis-Godson Temperature profile in K
  TEMP   ,   & ! Temperature profile in K
  TA     ,   & ! Surface air temperature in K
  TS     ,   & ! Skin temperature IN K
  CLDF         ! Fractional (ir) cloud cover


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NLEVSF ,   & ! Index of nearest standard pressure level
             ! at/below surface.
  FRACPS     ! Fraction of standard pressure level interval by
             ! which surf is above level nlevsf.

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCH       ! Max. number of IASI channels


IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: IPRINT
INTEGER, INTENT(IN) :: KNPF                  ! Number of processed profiles.
INTEGER, INTENT(IN) :: KNCHPF                ! Input number of radiances
                                             ! (channels*profiles).
INTEGER, INTENT(IN) :: KSAT

!       Array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN  (KNCHPF)       ! Input array of channel
                                             ! indices for each radiance.
INTEGER,INTENT(IN)  :: KPROF  (KNCHPF)       ! Input array of profile
                                             ! indices for each radiance.
REAL,INTENT(INOUT)  :: BCG    (KNCHPF,JPLEV)
REAL,INTENT(IN)     :: TAU    (KNCHPF,JPLEV) ! Transmittances from each
                                             ! level to space (chans*profs,
                                             ! standard pressure levels).
REAL,INTENT(IN)     :: TAUSFC (KNCHPF)       ! Transmittances from surface
                                             ! to space.
REAL,INTENT(IN)     :: PANGL  (KNPF)         ! Nadir angle.

!       Array arguments with intent out:
REAL,INTENT(OUT)    :: EMS    (KNCHPF)       ! Surface emissivities.
REAL,INTENT(OUT)    :: PRAD   (KNCHPF)       ! Output array of radiances.
REAL,INTENT(OUT)    :: PTB    (KNCHPF)       ! Output array of brightness
                                             ! temperatures.
REAL,INTENT(OUT)    :: B      (KNCHPF,JPLEV) ! Planck functions for
                                             ! temperatures profiles.
REAL,INTENT(OUT)    :: BA     (KNCHPF)       ! Planck funtions for surface
                                             ! air temperatures.
REAL,INTENT(OUT)    :: BS     (KNCHPF)       ! Planck funtions for surface
                                             ! skin temperatures.
REAL,INTENT(OUT)    :: RADOV  (KNCHPF,JPLEV) ! Overcast radiances for cloud
                                             ! at each standard pressure
                                             ! level.
REAL,INTENT(OUT)    :: RADO   (KNCHPF)       ! Overcast radiances for
                                             ! given cloud-top pressures.
REAL,INTENT(OUT)    :: BDT    (KNCHPF,JPLEV) ! Stores upwelling radiation
                                             ! from atmosphere above each
                                             ! level.
REAL,INTENT(OUT)    :: BDTR   (KNCHPF,JPLEV) ! Stores downwelling radiation
                                             ! at each level from atmosphere
                                             ! above.
REAL,INTENT(OUT)    :: FCLD   (KNCHPF)       ! Effective fractional cloud
                                             ! cover.
REAL,INTENT(OUT)    :: RADCL  (KNCHPF)       ! Clear column radiances.
REAL,INTENT(OUT)    :: TBCL   (KNCHPF)       ! Clear column brightness
                                             ! temperatures in K.


!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,JLEV,ILEVWRT,JLEV1,ILEV1,IPROF,ICLD
REAL :: ZBDT
REAL :: ZBDTR
REAL :: ZTAUMN
REAL :: ZFRAC

!       Local arrays:
REAL :: ZT(KNCHPF)
REAL :: ZTCG(KNCHPF)


!-----End of header-----------------------------------------------------


DATA ZTAUMN/1.0E-8/
TEMPW(:,:) =300.0


!-----------------------------------------------------------------------
!         1.   COMPUTE GURTIS-GODSON LAYER WEIGHTED TEMPERATURES
!-----------------------------------------------------------------------

IF(IPRINT==1)THEN
WRITE(*,'(A44)')'COMPUTING GODSON WEIGHTED LAYER TEMPERATURES'
ENDIF

DO J=1,KNPF
  ILEVWRT=NLEVSF(J)
    IF (FRACPS(J) < 0.) THEN
      ILEVWRT=ILEVWRT+1
    END IF
  CALL RTIASI_LAYERS(J,ILEVWRT,PANGL,KNPF)
ENDDO

IF(IPRINT==1)THEN
WRITE(*,'(A43)')'GODSON WEIGHTED LAYER TEMPERATURES COMPUTED'
WRITE(*,'(/)')
ENDIF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         2.   SET UP PLANCK FUNCTIONS, ETC.
!-----------------------------------------------------------------------

!-----2.1  Convert all temperatures to planck functions-----------------



DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    ZTCG(J)=TEMPW(JLEV,KPROF(J))
  ENDDO

!-----Compute the planck function for layer jlev using the Curtis-Godson
!     weighted temperature for that layer
!-----------------------------------------------------------------------
  CALL RTIASI_PLNCX(BCG(1:KNCHPF,JLEV),KCHAN,ZTCG,KNCHPF)
ENDDO


DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    ZT(J)=TEMP(JLEV,KPROF(J))
  ENDDO
  CALL RTIASI_PLNCX(B(1:KNCHPF,JLEV),KCHAN,ZT,KNCHPF)
ENDDO


!-----Surface air and skin temperatures---------------------------------

DO J=1,KNCHPF
  ZT(J)=TA(KPROF(J))
END DO
CALL RTIASI_PLNCX(BA,KCHAN,ZT,KNCHPF)

DO J=1,KNCHPF
  ZT(J)=TS(KPROF(J))
END DO
CALL RTIASI_PLNCX(BS,KCHAN,ZT,KNCHPF)
!-----------------------------------------------------------------------




IF(IPRINT==1)THEN
WRITE(*,'(A22)')'COMPUTING EMISSIVITIES'
ENDIF


!-----2.2  Set up surface emissivities.---------------------------------

CALL RTIASI_EMISS(KCHAN,KPROF,KNCHPF,EMS,PANGL,KNPF,KSAT)
!-----------------------------------------------------------------------

IF(IPRINT==1)THEN
WRITE(*,'(A21)')'EMISSIVITIES COMPUTED'
WRITE(*,'(/)')
ENDIF





!-----2.3  Set up cloud fractional coverage-----------------------------
!          Currently: same for all IASI channels

DO J=1,KNCHPF
  IF(KCHAN(J) <= JPCH) THEN
    FCLD(J)=CLDF(KPROF(J))
   ELSE
    FCLD(J)=0.
  ENDIF
ENDDO
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
!         3.   ATMOSPHERIC CONTRIBUTION.
!-----------------------------------------------------------------------

!-----LAYER ABOVE TOP PRESSURE LEVEL------------------------------------

DO J=1,KNCHPF
  BDT(J,1)=B(J,1)*(1.-TAU(J,1))      ! Upward
  BDTR(J,1)=BDT(J,1)/TAU(J,1)        ! Downward
ENDDO




!-----LAYERS BETWEEN STANDARD PRESSURE LEVELS---------------------------

DO JLEV=2,JPLEV
  JLEV1=JLEV-1
    DO J=1,KNCHPF


      ZBDT=BCG(J,JLEV-1)* &
      (TAU(J,JLEV1)-TAU(J,JLEV))
      BDT(J,JLEV)=BDT(J,JLEV1)+ZBDT
        IF (ABS(TAU(J,JLEV)) > ZTAUMN .AND. &
             ABS(TAU(J,JLEV1)) > ZTAUMN) THEN
          ZBDTR=BDTR(J,JLEV1)+ &
          ZBDT/(TAU(J,JLEV)*TAU(J,JLEV1))
        ELSE
          ZBDTR=BDTR(J,JLEV1)
        ENDIF
      BDTR(J,JLEV)=ZBDTR
    ENDDO
ENDDO


!-----NEAR-SURFACE LAYER + ADD UPWARD AND DOWNWARD PARTS----------------

DO J=1,KNCHPF
  ILEV1=NLEVSF(KPROF(J))-1
    IF (FRACPS(KPROF(J)) < 0.) THEN
      ILEV1=ILEV1+1
    END IF


  ZBDT=BCG(J,ILEV1)*(TAU(J,ILEV1)-TAUSFC(J))

    IF (ABS(TAUSFC(J)) > ZTAUMN .AND. ABS(TAU(J,ILEV1)) > ZTAUMN) THEN
      ZBDTR=BDTR(J,ILEV1)+ZBDT/(TAU(J,ILEV1)*TAUSFC(J))
    ELSE
      ZBDTR=BDTR(J,ILEV1)
    ENDIF
  ZBDT=BDT(J,ILEV1)+ZBDT
  RADCL(J)=ZBDT+ZBDTR*(1.-EMS(J))*TAUSFC(J)*TAUSFC(J)
ENDDO


!-----End of section 3--------------------------------------------------




!-----------------------------------------------------------------------
!         4.   SURFACE CONTRIBUTION.
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  RADCL(J)=RADCL(J)+EMS(J)*TAUSFC(J)*BS(J)
ENDDO
!           *** SOLAR REFLECTION BY SURFACE COULD BE ADDED HERE


!-----------------------------------------------------------------------
!         5.   CALCULATE OVERCAST RADIANCES.
!-----------------------------------------------------------------------

!----- 5.1  Add cloud emission to give overcast rads at all levels.-----
!           (and store them all for possible use outside
!           rtiasi, e.g. in cloud first-guess routine.)

DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    RADOV(J,JLEV)=BDT(J,JLEV)+B(J,JLEV)*TAU(J,JLEV)
  ENDDO
ENDDO
!               *** THIS ASSUMES BLACK (NON-REFLECTIVE) CLOUD.
!                   NON-BLACK CLOUD COULD BE ADDED AT THIS POINT
!                   QUITE EASILY, INCLUDING REFLECTION OF
!                   DOWNWELLING RADIATION STORED IN BDTR.
!-----------------------------------------------------------------------



!-----5.2  Interpolate to given cloud-top pressures---------------------

DO J=1,KNCHPF
  IPROF=KPROF(J)
  ICLD=NLEVCD(IPROF)
  ZFRAC=FRACPC(IPROF)
  RADO(J)=RADOV(J,ICLD)*(1.-ZFRAC)+RADOV(J,ICLD-1)*ZFRAC
ENDDO


!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         6.   CALCULATE (PARTLY) CLOUDY RADIANCES.
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  PRAD(J)=RADCL(J)*(1.-FCLD(J))+RADO(J)*FCLD(J)
ENDDO



!-----------------------------------------------------------------------
!         7.   CONVERT RADIANCES TO BRIGHTNESS TEMPERATURES.
!-----------------------------------------------------------------------


CALL RTIASI_BRIGV(PTB,PRAD,KCHAN,KNCHPF)       !    Cloudy
CALL RTIASI_BRIGV(TBCL,RADCL,KCHAN,KNCHPF)     !    Clear


RETURN
END SUBROUTINE RTIASI_RTINT
!     Performs integration of radiative transfer equation in RTIASI.
!SUBROUTINE RTINTK(KCHAN,KPROF,KNCHPF,PRAD_D,PRAD,PTB_D,PTB, &
!EMS_D,EMS,TAU_D,TAU,TAUSFC_D,TAUSFC,B_D,B,BCG_D,BCG,BA_D,BA,BS_D, &
!BS,RADOV_D,RADOV,RADO_D,RADO,BDT_D,BDT,BDTR_D,BDTR,FCLD_D,FCLD, &
!RADCL_D,RADCL,TBCL_D,TBCL,PANGL,KNPF,KSAT,IPRINT)
SUBROUTINE RTIASI_RTINTK(KCHAN,KPROF,KNCHPF,PRAD_D,PRAD,PTB_D,PTB, &
EMS_D,EMS,TAU_D,TAU,TAUSFC_D,TAUSFC,B_D,B,BCG_D,BCG,BA_D,BA,BS_D, &
BS,RADOV_D,RADOV,RADO_D,RADO,BDT,BDTR_D,BDTR,FCLD_D,FCLD, &
RADCL_D,RADCL,PANGL,KNPF,KSAT,IPRINT)

!     Description
!     K of subroutine RTINT
!     Perform integration of radiative transfer equation
!     in RTIASI suite, calculating radiances and brightness
!     temperatures.
!
!     Method:
!     1)See: Marco Matricardi and Roger Saunders:Fast radiative transfer model
!       for simulation of infrared atmospheric sounding interferometer
!       radiances. Applied Optics, 38, pp. 5679-5691 (1999).
!     2) See: User's manual for RTIASI (Available from EUMETSAT)
!     3) See: ECMWF Technical Memorandum 176 (Available from ECMWF)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            25/06/91.   Original Code. J.R.Eyre. ECMWF.
!     2            12/11/1999. Modified to simulate radiances for IASI.
!                              Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct module used

USE RTIASI_IRCLD, ONLY : &
!     Imported arrays:
  NLEVCD_D=>NLEVCD, &  ! Index of nearest standard pressure level
                       ! at/below cloud top.
  FRACPC_D=>FRACPC     ! Fraction of standardd pressure level interval
                       ! by which cloud is above level NLEVCD.


USE RTIASI_PRFVAR,ONLY : &
!     Imported parameters:
  JPLEV , &               ! Number of pressure levels
!     Imported arrays:
  TEMPW_D =>TEMPW ,    & ! Curtis-Godson temperature profile
  TEMP_D  =>TEMP  ,    & ! Temperature profile in K
  TA_D    =>TA    ,    & ! Surface air temperature in K
  TS_D    =>TS           ! Skin temperature IN K


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NLEVSF_D=>NLEVSF,    & ! Index of nearest standard pressure level
                       ! at/below surface.
  FRACPS_D=>FRACPS     ! Fraction of standard pressure level interval by
                       ! which surf is above level nlevsf.


USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCH                 ! Max. number of IASI channels


!     Tangent linear module used:


USE RTIASI_PRFVARK,ONLY : &
!     Imported arrays:
  TEMPW           ,    & ! Curtis-Godson temperature profile
  TEMP            ,    & ! Temperature profile in K
  CLDF            ,    & ! Fractional (ir) cloud cover
  TA              ,    & ! Surface air temperature in K
  TS                     ! Skin temperature IN K


USE RTIASI_IRCLDK, ONLY : &
!     Imported arrays:
  FRACPC               ! Fraction of standardd pressure level interval
                       ! by which cloud is above level NLEVCD.


IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: IPRINT
INTEGER, INTENT(IN) :: KNPF                    ! Number of processed prof.
INTEGER, INTENT(IN) :: KNCHPF                  ! Input number of radianc.
                                               ! (channels*profiles).
INTEGER, INTENT(IN) :: KSAT

!       Direct array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN    (KNCHPF)       ! Input array of channel
                                               ! indices for each radiance
INTEGER,INTENT(IN)  :: KPROF    (KNCHPF)       ! Input array of profile
                                               ! indices for each radiance
REAL,INTENT(IN)     :: PANGL    (KNPF)         ! Nadir angle.
REAL,INTENT(IN)     :: TAU_D    (KNCHPF,JPLEV) ! Transmittances from each
                                               ! lev. to top (chans*profs,
                                               ! standard pressure levels)
REAL,INTENT(IN)     :: TAUSFC_D (KNCHPF)       ! Transmittances from
                                               ! surface to space.
REAL,INTENT(IN)     :: PRAD_D   (KNCHPF)       ! Output array of radiances
REAL,INTENT(IN)     :: PTB_D    (KNCHPF)       ! Output array of bright.
                                               ! temperatures.
REAL,INTENT(IN)     :: EMS_D    (KNCHPF)       ! Surface emissivities.
REAL,INTENT(IN)     :: B_D      (KNCHPF,JPLEV) ! Planck functions for
                                               ! temperatures profiles.
REAL,INTENT(IN)     :: BCG_D    (KNCHPF,JPLEV)
REAL,INTENT(IN)     :: BA_D     (KNCHPF)       ! Planck funtions for surf
                                               ! air temperatures.
REAL,INTENT(IN)     :: BS_D     (KNCHPF)       ! Planck funtions for surf
                                               ! skin temperatures.
REAL,INTENT(IN)     :: RADOV_D  (KNCHPF,JPLEV) ! Overcast radiances for
                                               ! cloud at each standard
                                               ! pressure level.
REAL,INTENT(IN)     :: RADO_D   (KNCHPF)       ! Overcast radiances for
                                               ! given cloud-top pressures
REAL,INTENT(IN)     :: BDTR_D   (KNCHPF,JPLEV) ! Stores downwelling
                                               ! radiation at each level
                                               ! from atmosphere above.
REAL,INTENT(IN)     :: FCLD_D   (KNCHPF)       ! Effective fractional
                                               ! cloud cover.
REAL,INTENT(IN)     :: RADCL_D  (KNCHPF)       ! Clear column radiances.

!       K array arguments with intent in:
REAL,INTENT(INOUT)  :: EMS      (KNCHPF)       ! Surface emissivities.
REAL,INTENT(INOUT)  :: PRAD     (KNCHPF)       ! Output array of radiances
REAL,INTENT(INOUT)  :: PTB      (KNCHPF)       ! Output array of bright.
                                               ! temperatures.
REAL,INTENT(INOUT)  :: B        (KNCHPF,JPLEV) ! Planck functions for
                                               ! temperatures profiles.
REAL,INTENT(INOUT)  :: BCG      (KNCHPF,JPLEV)
REAL,INTENT(INOUT)  :: BA       (KNCHPF)       ! Planck funtions for surf
                                               ! air temperatures.
REAL,INTENT(INOUT)  :: BS       (KNCHPF)       ! Planck funtions for surf
                                               ! skin temperatures.
REAL,INTENT(INOUT)  :: RADOV    (KNCHPF,JPLEV) ! Overcast radiances for
                                               ! cloud at each standard
                                               ! pressure level.
REAL,INTENT(INOUT)  :: RADO     (KNCHPF)       ! Overcast radiances for
                                               ! given cloud-top pressures
REAL,INTENT(INOUT)  :: BDT      (KNCHPF,JPLEV) ! Stores upwelling
                                               ! radiation from atmosph.
                                               ! above each level.
REAL,INTENT(INOUT)  :: BDTR     (KNCHPF,JPLEV) ! Stores downwelling
                                               ! radiation at each level
                                               ! from atmosphere above.
REAL,INTENT(INOUT)  :: FCLD     (KNCHPF)       ! Effective fractional
                                               ! cloud cover.
REAL,INTENT(INOUT)  :: RADCL    (KNCHPF)       ! Clear column radiances.

!       K array arguments with intent out:
REAL,INTENT(INOUT)  :: TAU      (KNCHPF,JPLEV) ! Transmittances from each
                                               ! lev. to top (chans*profs,
                                               ! standard pressure levels)
REAL,INTENT(INOUT)  :: TAUSFC   (KNCHPF)       ! Transmittances from
                                               ! surface to space.


!     End of subroutine arguments


!       Direct local scalars:
INTEGER :: J,JLEV,JLEV1,ILEV1,IPROF,ICLD
REAL    :: ZBDT_D
REAL    :: ZBDTR_D
REAL    :: ZTAUMN
REAL    :: ZFRAC_D

!       K local scalars:
REAL    :: ZBDT
REAL    :: ZBDTR

!       Direct local arrays:
REAL    :: ZT_D   (KNCHPF)
REAL    :: ZTCG_D (KNCHPF)
REAL    :: ZIC

!       K local arrays:
REAL    :: ZT     (KNCHPF)
REAL    :: ZTCG   (KNCHPF)



!-----End of header-----------------------------------------------------

ZBDT  =0.
ZBDTR =0.
FCLD(:)=0


ZTAUMN = 1.0E-8

!-----------------------------------------------------------------------
!         7.   CONVERT RADIANCES TO BRIGHTNESS TEMPERATURES.
!-----------------------------------------------------------------------

!      PRAD(:)=0
!      RADCL(:)=0.

CALL RTIASI_BRIGVK(PTB_D,PTB,PRAD_D,PRAD,KCHAN,KNCHPF)       !    CLOUDY
!     CALL BRIGVK(TBCL_D,TBCL,RADCL_D,RADCL,KCHAN,KNCHPF)   !    CLEAR



!-----------------------------------------------------------------------
!         6.   CALCULATE (PARTLY) CLOUDY RADIANCES.
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  RADCL(J)=RADCL(J)+PRAD(J)*(1.-FCLD_D(J))
  RADO(J) =RADO(J)+PRAD(J)*FCLD_D(J)
  FCLD(J) =FCLD(J)+PRAD(J)*(RADO_D(J)-RADCL_D(J))
ENDDO


!-----------------------------------------------------------------------
!         5.   CALCULATE OVERCAST RADIANCES.
!-----------------------------------------------------------------------

!-----5.2  Interpolate to given cloud-top pressures---------------------

DO J=1,KNCHPF
  IPROF=KPROF(J)
  ICLD=NLEVCD_D(IPROF)
  ZFRAC_D=FRACPC_D(IPROF)
  RADOV(J,ICLD)=RADOV(J,ICLD)+RADO(J)*(1.-ZFRAC_D)
  RADOV(J,ICLD-1)=RADOV(J,ICLD-1)+RADO(J)*ZFRAC_D
  FRACPC(J)=FRACPC(J)+ &
  (RADOV_D(J,ICLD-1)-RADOV_D(J,ICLD))*RADO(J)
ENDDO



!----- 5.1  Add cloud emission to give overcast rads at all levels.-----
!           (and store them all for possible use outside
!           rtiasi, e.g. in cloud first-guess routine.)

DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    BDT(J,JLEV)=BDT(J,JLEV)+RADOV(J,JLEV)
    B(J,JLEV)=B(J,JLEV)+RADOV(J,JLEV)*TAU_D(J,JLEV)
    TAU(J,JLEV)=TAU(J,JLEV)+RADOV(J,JLEV)*B_D(J,JLEV)
  ENDDO
ENDDO
!               *** THIS ASSUMES BLACK (NON-REFLECTIVE) CLOUD.
!                   NON-BLACK CLOUD COULD BE ADDED AT THIS POINT
!                   QUITE EASILY, INCLUDING REFLECTION OF
!                   DOWNWELLING RADIATION STORED IN BDTR.
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!         4.   SURFACE CONTRIBUTION.
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  EMS(J)=EMS(J)+RADCL(J)*TAUSFC_D(J)*BS_D(J)
  TAUSFC(J)=TAUSFC(J)+RADCL(J)*EMS_D(J)*BS_D(J)
  BS(J)=BS(J)+RADCL(J)*EMS_D(J)*TAUSFC_D(J)
ENDDO
!           *** SOLAR REFLECTION BY SURFACE COULD BE ADDED HERE



!-----------------------------------------------------------------------
!         3.   ATMOSPHERIC CONTRIBUTION.
!-----------------------------------------------------------------------


!-----NEAR-SURFACE LAYER + ADD UPWARD AND DOWNWARD PARTS----------------

DO J=1,KNCHPF

  ILEV1=NLEVSF_D(KPROF(J))-1
    IF (FRACPS_D(KPROF(J)) < 0.) THEN
      ILEV1=ILEV1+1
    END IF

  ZBDT_D=BCG_D(J,ILEV1)*(TAU_D(J,ILEV1)-TAUSFC_D(J))
    IF (ABS(TAUSFC_D(J)) > ZTAUMN .AND. ABS(TAU_D(J,ILEV1)) > ZTAUMN) &
         THEN
      ZIC=1./(TAU_D(J,ILEV1)*TAUSFC_D(J))
      ZBDTR_D=BDTR_D(J,ILEV1)+ZBDT_D*ZIC
    ELSE
      ZBDTR_D=BDTR_D(J,ILEV1)
    ENDIF


  ZBDT=ZBDT+RADCL(J)
  ZBDTR=ZBDTR+TAUSFC_D(J)**2*(1.-EMS_D(J))*RADCL(J)
  TAUSFC(J)=TAUSFC(J)+2.*TAUSFC_D(J)*ZBDTR_D*(1.-EMS_D(J))*RADCL(J)
  EMS(J)=EMS(J)-TAUSFC_D(J)**2*ZBDTR_D*RADCL(J)
  BDT(J,ILEV1)=BDT(J,ILEV1)+ZBDT

  IF (ABS(TAUSFC_D(J)) > ZTAUMN .AND. ABS(TAU_D(J,ILEV1)) > ZTAUMN) THEN
     BDTR(J,ILEV1)=BDTR(J,ILEV1)+ZBDTR
     ZBDT=ZBDT+ZBDTR*ZIC
     TAU(J,ILEV1)=TAU(J,ILEV1)-(ZBDT_D*ZIC)*ZBDTR*TAUSFC_D(J)*ZIC
     TAUSFC(J)=TAUSFC(J)-TAU_D(J,ILEV1)*ZBDTR*ZIC*ZBDT_D*ZIC
     ZBDTR=0.
  ELSE
     BDTR(J,ILEV1)=BDTR(J,ILEV1)+ZBDTR
     ZBDTR=0.
  ENDIF

  BCG(J,ILEV1)=BCG(J,ILEV1)+ZBDT*(TAU_D(J,ILEV1)-TAUSFC_D(J))
  TAU(J,ILEV1)=TAU(J,ILEV1)+ZBDT*BCG_D(J,ILEV1)
  TAUSFC(J)=TAUSFC(J)-ZBDT*BCG_D(J,ILEV1)
  ZBDT=0.

ENDDO

!-----LAYERS BETWEEN STANDARD PRESSURE LEVELS---------------------------

DO JLEV=JPLEV,2,-1
  JLEV1=JLEV-1
    DO J=1,KNCHPF
      ZBDT_D=BCG_D(J,JLEV-1)* &
      (TAU_D(J,JLEV1)-TAU_D(J,JLEV))

    ZBDTR=ZBDTR+BDTR(J,JLEV)

      IF (ABS(TAU_D(J,JLEV)) > ZTAUMN .AND. &
           ABS(TAU_D(J,JLEV1)) > ZTAUMN) THEN
        ZIC=1./(TAU_D(J,JLEV)*TAU_D(J,JLEV1))
        BDTR(J,JLEV1)=BDTR(J,JLEV1)+ZBDTR
        ZBDT=ZBDT+ZBDTR*ZIC
     TAU(J,JLEV)=TAU(J,JLEV)-(ZBDT_D*ZIC)*ZBDTR*TAU_D(J,JLEV1)*ZIC
      TAU(J,JLEV1)=TAU(J,JLEV1)-TAU_D(J,JLEV)*ZBDTR*ZIC*ZBDT_D*ZIC
        ZBDTR=0.
      ELSE
        BDTR(J,JLEV1)=BDTR(J,JLEV1)+ZBDTR
        ZBDTR=0.
      ENDIF
    BDT(J,JLEV1)=BDT(J,JLEV1)+BDT(J,JLEV)
    ZBDT=ZBDT+BDT(J,JLEV)

    BCG(J,JLEV-1)=BCG(J,JLEV-1)+ZBDT*(TAU_D(J,JLEV1)-TAU_D(J,JLEV) &
)
    TAU(J,JLEV1)=TAU(J,JLEV1)+ZBDT*BCG_D(J,JLEV-1)
    TAU(J,JLEV)=TAU(J,JLEV)-ZBDT*BCG_D(J,JLEV-1)
    ZBDT=0.
  ENDDO
ENDDO

!-----LAYER ABOVE TOP PRESSURE LEVEL------------------------------------

DO J=1,KNCHPF
  BDT(J,1)=BDT(J,1)+BDTR(J,1)/TAU_D(J,1)
  TAU(J,1)=TAU(J,1)-BDTR_D(J,1)*BDTR(J,1)/TAU_D(J,1)

  B(J,1)=B(J,1)+(1.-TAU_D(J,1))*BDT(J,1)
  TAU(J,1)=TAU(J,1)-B_D(J,1)*BDT(J,1)
ENDDO

!-----------------------------------------------------------------------
!         2.   SET UP PLANCK FUNCTIONS, ETC.
!-----------------------------------------------------------------------

!-----2.3  Set up cloud fractional coverage-----------------------------
!          Currently: same for all IASI channels

DO J=1,KNCHPF
  IF(KCHAN(J) <= JPCH) THEN
    CLDF(J)=CLDF(J)+FCLD(J)
  ELSE
    FCLD(J)=0.
  ENDIF
ENDDO

IF(IPRINT==1)THEN
WRITE(*,'(A22)')'COMPUTING EMISSIVITIES'
ENDIF

!-----2.2  Set up surface emissivities.---------------------------------

CALL RTIASI_EMISSK(KCHAN,KPROF,KNCHPF,EMS,PANGL,KNPF,KSAT)
!-----------------------------------------------------------------------

IF(IPRINT==1)THEN
WRITE(*,'(A21)')'EMISSIVITIES COMPUTED'
WRITE(*,'(/)')
ENDIF

!-----2.1  Convert all temperatures to planck functions-----------------


DO J=1,KNCHPF
  ZT_D(J)=TS_D(KPROF(J))
ENDDO

ZT(:)=0.

CALL RTIASI_PLNCXK(BS_D,BS,KCHAN,ZT_D,ZT,KNCHPF)

DO J=1,KNCHPF
  TS(J)=TS(J)+ZT(J)
ENDDO


DO J=1,KNCHPF
  ZT_D(J)=TA_D(KPROF(J))
ENDDO

ZT(:)=0.

CALL RTIASI_PLNCXK(BA_D,BA,KCHAN,ZT_D,ZT,KNCHPF)

DO J=1,KNCHPF
  TA(J)=TA(J)+ZT(J)
ENDDO

ZT(:)=0.


DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    ZT_D(J)=TEMP_D(JLEV,KPROF(J))
  ENDDO

  CALL RTIASI_PLNCXK(B_D(1:KNCHPF,JLEV),B(1:KNCHPF,JLEV),KCHAN,ZT_D, &
               ZT,KNCHPF)
  DO J=1,KNCHPF
    TEMP(JLEV,J)=TEMP(JLEV,J)+ZT(J)
    ZT(J)=0.
  ENDDO
ENDDO



ZTCG(:)  =0.



DO JLEV=1,JPLEV
  DO J=1,KNCHPF
    ZTCG_D(J)=TEMPW_D(JLEV,KPROF(J))
  ENDDO

  CALL RTIASI_PLNCXK(BCG_D(1:KNCHPF,JLEV),BCG(1:KNCHPF,JLEV),KCHAN, &
               ZTCG_D,ZTCG,KNCHPF)

  DO J=1,KNCHPF
     TEMPW(JLEV,J)=TEMPW(JLEV,J)+ZTCG(J)
     ZTCG(J)=0.
  ENDDO
ENDDO



!-----------------------------------------------------------------------
!         1.   COMPUTE GURTIS-GODSON LAYER WEIGHTED TEMPERATURES
!-----------------------------------------------------------------------

IF(IPRINT==1)THEN
WRITE(*,'(A66)') &
'COMPUTING CURTIS-GODSON TANGENT LINEAR WEIGHTED LAYER TEMPERATURES'
ENDIF


  CALL RTIASI_LAYERSK(PANGL,KNPF,KNCHPF,KPROF)

IF(IPRINT==1)THEN
WRITE(*,'(A65)') &
'TANGENT LINEAR CURTIS-GODSON WEIGHTED LAYER TEMPERATURES COMPUTED'
WRITE(*,'(/)')
ENDIF

!-----------------------------------------------------------------------



RETURN
END SUBROUTINE RTIASI_RTINTK
!     Compute sea surface emissivity for radiative transfer calculations
SUBROUTINE RTIASI_EMISS(KCHAN,KPROF,KNCHPF,PEMS,PANGL,KNPF,KSAT)

!     Description:
!     Set up surface emissivities for all channels and all
!     profiles.
!
!     Method:
!     See K.MASUDA,T.TAKASHIMA,Y.TAKAYAMA.
!     Emissivity of pure water and sea waters for the
!     sea surface in the infrared window regions.
!     Remote Sensing of Environment 24:313-329 (1988)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            21/03/99.   Original Code. M. Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOCON,ONLY : &
!     Imported scalars:
 RATOE       ! Ratio satellite-orbit/earth radii


USE RTIASI_EMISIR,ONLY : &
!     Imported arrays:
  EMSIR      ! Infrared emissivities


USE RTIASI_PRFVAR,ONLY : &
!     Imported arrays:
  SURFWU,    & ! U component of surface wind in m/sec
  SURFWV     ! V component of surface wind in m/sec


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NSTYPE     ! Surface type index; 1=sea, 2=land.


USE RTIASI_IASCHN,ONLY : &
!     Imported parameters:
  JPCH,      & ! Max. number of IASI channels
!     Imported arrays:
  EMSCOEF    ! Regression coefficients used to compute emissivities


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNCHPF        ! Input number of radiances
                                    ! (chans * profs)
INTEGER,INTENT(IN) :: KNPF          ! Number of processed profiles
INTEGER,INTENT(IN) :: KSAT

!       Array arguments with intent in:
INTEGER,INTENT(IN) :: KCHAN(KNCHPF) ! Input array of channel indices
                                    ! for each radiance.
INTEGER,INTENT(IN) :: KPROF(KNCHPF) ! Input array of profile indices
                                    ! for each radiance.
REAL,INTENT(IN)    :: PANGL(KNPF)   ! Zenith angle

!       Array arguments with intent out:
REAL,INTENT(OUT)   :: PEMS(KNCHPF)  ! Emissivities

!     End of subroutine arguments


!       Local parameters:


!       Local scalars:
INTEGER :: J,ICHAN,IPROF,ISTYPE
REAL :: WINDSP         ! Surface wind speed.
REAL :: ZANGLE         ! Local zenith angle
REAL :: TETA           ! Nadir angle
REAL :: ZSIN
REAL :: PI

!       Local arrays:
REAL :: A(JPCH)
REAL :: B(JPCH)
REAL :: CC(JPCH)


!-----End of header-----------------------------------------------------

DATA PI/3.14159265359/

!-----------------------------------------------------------------------
!         1.   SETS UP SURFACE EMISSIVITIES.
!-----------------------------------------------------------------------



DO J=1,KNCHPF                           !Loop over channels.
  ICHAN=KCHAN(J)
  IPROF=KPROF(J)
  ISTYPE=NSTYPE(IPROF)

  TETA=(PANGL(IPROF)/180.)*PI           !Convert nadir angle to radiants.
  ZSIN=SIN(TETA)*RATOE(KSAT)
  ZANGLE=180.*ASIN(ZSIN)/PI             !Zenith angle

  WINDSP=SQRT(SURFWU(IPROF)**2+         & ! Surface wind speed.
  SURFWV(IPROF)**2)


    IF(ISTYPE == 1)THEN                 !Obtain ir emissivities over sea.

      IF(EMSIR(ICHAN) == 0. ) THEN
        A(ICHAN)=EMSCOEF(1,ICHAN)+EMSCOEF(2,ICHAN)*WINDSP+ &
                 EMSCOEF(3,ICHAN)*WINDSP**2

        B(ICHAN)=EMSCOEF(4,ICHAN)+EMSCOEF(5,ICHAN)*WINDSP+ &
                 EMSCOEF(6,ICHAN)*WINDSP**2


        CC(ICHAN)=EMSCOEF(7,ICHAN)+EMSCOEF(8,ICHAN)*WINDSP

        PEMS(J)= A(ICHAN)+ (B(ICHAN)-A(ICHAN))* &
                 Exp(( (EMSCOEF(9,ICHAN)-60.)**2. &
                    - (ZANGLE-EMSCOEF(9,ICHAN))**2. )/ &
                    CC(ICHAN))


      ELSE IF ( EMSIR(ICHAN) /= 0. ) THEN
        PEMS(J)=EMSIR(ICHAN)
      END IF
    ELSE IF (ISTYPE == 2)THEN           ! Obtain ir emissivities over land
      PEMS(J)=0.9
    ENDIF


ENDDO                                   !End of channel loop


RETURN
END SUBROUTINE RTIASI_EMISS
!     Compute sea surface emissivity for radiative transfer calculations
SUBROUTINE RTIASI_EMISSK(KCHAN,KPROF,KNCHPF,PEMS,PANGL,KNPF,KSAT)

!     Description:
!     Set up surface emissivities for all channels and all
!     profiles.
!
!     Method:
!     See K.MASUDA,T.TAKASHIMA,Y.TAKAYAMA.
!     Emissivity of pure water and sea waters for the
!     sea surface in the infrared window regions.
!     Remote Sensing of Environment 24:313-329 (1988)
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            21/03/99.   Original Code. M. Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_GEOCON,ONLY : &
!     Imported scalars:
 RATOE       ! Ratio satellite-orbit/earth radii


USE RTIASI_EMISIR,ONLY : &
!     Imported arrays:
  EMSIR_D=>EMSIR      ! Infrared emissivities


USE RTIASI_PRFVAR,ONLY : &
!     Imported arrays:
  SURFWU_D=>SURFWU, & ! U component of surface wind in m/sec
  SURFWV_D=>SURFWV  ! V component of surface wind in m/sec


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NSTYPE_D=>NSTYPE  ! Surface type index; 1=sea, 2=land.


USE RTIASI_IASCHN,ONLY : &
!     Imported parameters:
  JPCH,      & ! Max. number of IASI channels
!     Imported arrays:
  EMSCOEF

!     Tangent linear module used:


USE RTIASI_PRFVARK,ONLY : &
!     Imported tangent linear arrays:
  SURFWU,    & ! U component of surface wind in m/sec
  SURFWV     ! V component of surface wind in m/sec


USE RTIASI_EMISIRK,ONLY : &
!     Imported tangent linear arrays:
  EMSIR      ! Infrared emissivities


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNCHPF        ! Input number of radiances
                                    ! (chans * profs)
INTEGER,INTENT(IN) :: KNPF          ! Number of processed profiles
INTEGER,INTENT(IN) :: KSAT

!       Array arguments with intent in:
INTEGER,INTENT(IN) :: KCHAN(KNCHPF) ! Input array of channel indices
                                    ! for each radiance.
INTEGER,INTENT(IN) :: KPROF(KNCHPF) ! Input array of profile indices
                                    ! for each radiance.
REAL,INTENT(IN)    :: PANGL(KNPF)   ! Zenith angle

!       Array arguments with intent inout:
REAL,INTENT(INOUT) :: PEMS(KNCHPF)  ! Emissivities

!     End of subroutine arguments

!       Local scalars:
INTEGER :: J,ICHAN,IPROF,ISTYPE
REAL :: WINDSP_D         ! Surface wind speed.
REAL :: EXPF
REAL :: FAC
REAL :: ZANGLE         ! Local zenith angle
REAL :: TETA           ! Nadir angle
REAL :: ZSIN
REAL :: PI

!       Local arrays:
REAL :: A_D(JPCH)
REAL :: B_D(JPCH)
REAL :: CC_D(JPCH)

!       Local tangent linear arrays:
REAL :: A(JPCH)
REAL :: B(JPCH)
REAL :: CC(JPCH)

!       Local tangent linear scalars:
REAL :: WINDSP         ! Surface wind speed.


!-----End of header-----------------------------------------------------

DATA PI/3.14159265359/

!-----------------------------------------------------------------------
!         1.   SETS UP SURFACE EMISSIVITIES.
!-----------------------------------------------------------------------

!      SURFWU(:)=0.
!      SURFWV(:)=0.
!      EMSIR(:)=0.

A(:)   =0.
B(:)   =0.
CC(:)  =0.
WINDSP =0.


DO J=1,KNCHPF                           !Loop over channels.
  ICHAN=KCHAN(J)
  IPROF=KPROF(J)
  ISTYPE=NSTYPE_D(IPROF)
  ISTYPE=1

  TETA=(PANGL(IPROF)/180.)*PI           !Convert nadir angle to radiants.
  ZSIN=SIN(TETA)*RATOE(KSAT)
  ZANGLE=180.*ASIN(ZSIN)/PI

  WINDSP_D=SQRT(SURFWU_D(IPROF)**2+          & ! Direct surface wind speed.
  SURFWV_D(IPROF)**2)




    IF(ISTYPE == 1)THEN                 !Obtain ir emissivities over sea.

      IF(EMSIR_D(ICHAN) == 0. ) THEN
        A_D(ICHAN)=EMSCOEF(1,ICHAN)+EMSCOEF(2,ICHAN)*WINDSP_D+ &
                 EMSCOEF(3,ICHAN)*WINDSP_D**2

        B_D(ICHAN)=EMSCOEF(4,ICHAN)+EMSCOEF(5,ICHAN)*WINDSP_D+ &
                 EMSCOEF(6,ICHAN)*WINDSP_D**2


        CC_D(ICHAN)=EMSCOEF(7,ICHAN)+EMSCOEF(8,ICHAN)*WINDSP_D

        EXPF=Exp(( (EMSCOEF(9,ICHAN)-60.)**2. &
                    - (ZANGLE-EMSCOEF(9,ICHAN))**2. )/ &
                    CC_D(ICHAN))
        FAC=-( (EMSCOEF(9,ICHAN)-60.)**2. &
                    - (ZANGLE-EMSCOEF(9,ICHAN))**2. )/ &
                    CC_D(ICHAN)**2

        A(ICHAN)=A(ICHAN)+ PEMS(J)
        B(ICHAN)=B(ICHAN)+PEMS(J)*EXPF
        CC(ICHAN)=CC(ICHAN)+B_D(ICHAN)*EXPF*PEMS(J)*FAC
        A(ICHAN)=A(ICHAN)-PEMS(J)*EXPF
        CC(ICHAN)=CC(ICHAN)-A_D(ICHAN)*EXPF*PEMS(J)*FAC

        WINDSP=WINDSP+EMSCOEF(8,ICHAN)*CC(ICHAN)
        WINDSP=WINDSP+(EMSCOEF(5,ICHAN)+2*EMSCOEF(6,ICHAN)* &
               WINDSP_D)*B(ICHAN)
        WINDSP=WINDSP+(EMSCOEF(2,ICHAN)+2*EMSCOEF(3,ICHAN)* &
               WINDSP_D)*A(ICHAN)


      ELSE IF ( EMSIR_D(ICHAN) /= 0. ) THEN
        EMSIR(ICHAN)=EMSIR(ICHAN)+PEMS(J)
      END IF
    ELSE IF (ISTYPE == 2)THEN           ! Obtain ir emissivities over land
      PEMS(J)=0.
    ENDIF

    IF ( SURFWU_D(IPROF) /= 0. .AND. SURFWV_D(IPROF) /= .0 ) THEN
      SURFWU(J)=SURFWU(J)+SURFWU_D(IPROF)*WINDSP/WINDSP_D
      SURFWV(J)=SURFWV(J)+SURFWV_D(IPROF)*WINDSP/WINDSP_D
      WINDSP=0.
    ELSE
      SURFWU(J)=SURFWU(J)+WINDSP/SQRT(2.)
      SURFWV(J)=SURFWV(J)+WINDSP/SQRT(2.)
      WINDSP=0.
    ENDIF


ENDDO                                   !End of channel loop


RETURN
END SUBROUTINE RTIASI_EMISSK
!     Calculate a  weighted integrated quantity using simpsons rule
SUBROUTINE RTIASI_INTEG(MXDIM,NREC,NP,A,B,ARX,ARY,MODY,SWITCH,ARW1,MOW1, &
ARW2,MOW2,XINT)

!     Description:
!     To calculate a  weighted integrated quantity using
!     Simpsons rule
!
!     Method:
!     A NP*2 panel Simpsons rule calculation
!     is performed between the limits A and B for the
!     quantity ARY*ARW1*ARW2 and the result returned as
!     XINT. The function values are interpolated using
!     appropriate schemes specified by mode.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            27/02/91    Original Code. D.P. Edwards
!     2            12/11/1999  Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


 IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER ,INTENT(IN) :: MXDIM 
INTEGER ,INTENT(IN) :: NREC  
INTEGER ,INTENT(IN) :: NP    
INTEGER ,INTENT(IN) :: MODY  
INTEGER ,INTENT(IN) :: MOW1  
INTEGER ,INTENT(IN) :: MOW2  
REAL    ,INTENT(IN) :: A     
REAL    ,INTENT(IN) :: B     
LOGICAL ,INTENT(IN) :: SWITCH

!       Array arguments with intent in:
REAL    ,INTENT(IN) :: ARX  (MXDIM)
REAL    ,INTENT(IN) :: ARY  (MXDIM)
REAL    ,INTENT(IN) :: ARW1 (MXDIM)
REAL    ,INTENT(IN) :: ARW2 (MXDIM)

!       Array scalars with intent out:
REAL    ,INTENT(OUT):: XINT


!     End of subroutine arguments


!       Local scalars:
INTEGER :: I
REAL    :: HINC
REAL    :: SS
REAL    :: S2
REAL    :: S4
REAL    :: SY
REAL    :: SA
REAL    :: SB
REAL    :: ARG
REAL    :: H
REAL    :: SW1
REAL    :: SW2


!-----End of header-----------------------------------------------------

!-----Simpsons rule integration pannel width----------------------------
HINC = (B - A)/FLOAT(2*NP)


!-----Argument for function value at lower integration boundary---------
IF (.NOT. SWITCH) THEN
  CALL RTIASI_INTER(MXDIM,NREC,MODY,A,ARX,ARY,SY,H)
  SA = SY
ELSE
  CALL RTIASI_INTER(MXDIM,NREC,MODY,A,ARX,ARY,SY,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW1,A,ARX,ARW1,SW1,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW2,A,ARX,ARW2,SW2,H)
  SA = SY*SW1*SW2
ENDIF
!-----------------------------------------------------------------------


!-----Argument for function value at upper integration boundary---------
IF (.NOT. SWITCH) THEN
  CALL RTIASI_INTER(MXDIM,NREC,MODY,B,ARX,ARY,SY,H)
  SB = SY
ELSE
  CALL RTIASI_INTER(MXDIM,NREC,MODY,B,ARX,ARY,SY,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW1,B,ARX,ARW1,SW1,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW2,B,ARX,ARW2,SW2,H)
  SB = SY*SW1*SW2
ENDIF
!-----------------------------------------------------------------------


S4 = 0.0
DO I=1,NP
  ARG = A + HINC*FLOAT(2*I - 1) ! Argument for function values at even
                                ! integration points.

    IF (.NOT. SWITCH) THEN
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      SS = SY
    ELSE
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW1,ARG,ARX,ARW1,SW1,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW2,ARG,ARX,ARW2,SW2,H)
      SS = SY*SW1*SW2
    ENDIF
  S4 = S4 + SS              ! Sum function values a even integ. points
ENDDO


S2 = 0.0
DO I=1,NP-1
  ARG = A + HINC*FLOAT(2*I) ! Argument for function values at
                            ! odd integration points.
    IF (.NOT. SWITCH) THEN
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      SS = SY
    ELSE
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW1,ARG,ARX,ARW1,SW1,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW2,ARG,ARX,ARW2,SW2,H)
      SS = SY*SW1*SW2
    ENDIF
  S2 = S2 + SS              ! Sum function values at odd integ. points
ENDDO

XINT = HINC/3.0 * (SA + SB + 4.0*S4 + 2.0*S2)          ! Simpsons rule

RETURN
END SUBROUTINE RTIASI_INTEG
!     Calculate a  weighted integrated quantity using simpsons rule
SUBROUTINE RTIASI_INTEGTL(MXDIM,NREC,NP,A_D,A,B_D,B,ARX,ARY_D,ARY,MODY, &
SWITCH,ARW1,MOW1,ARW2,MOW2,XINT_D,XINT)

!     Description:
!     Tangent linear of subroutine INTEG
!     To calculate a  weighted integrated quantity using
!     Simpsons rule
!
!     Method:
!     A NP*2 panel Simpsons rule calculation
!     is performed between the limits A and B for the
!     quantity ARY*ARW1*ARW2 and the result returned as
!     XINT. The function values are interpolated using
!     appropriate schemes specified by mode.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            27/02/91    Original Code. D.P. Edwards
!     2            12/11/1999  Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


 IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: MXDIM
INTEGER,INTENT(IN) :: NREC
INTEGER,INTENT(IN) :: NP
INTEGER,INTENT(IN) :: MODY
INTEGER,INTENT(IN) :: MOW1
INTEGER,INTENT(IN) :: MOW2
REAL,INTENT(IN)    :: A_D
REAL,INTENT(IN)    :: B_D
LOGICAL,INTENT(IN) :: SWITCH

!       Scalar tangent linear arguments with intent in:
REAL,INTENT(IN)    :: A
REAL,INTENT(IN)    :: B

!       Array arguments with intent in:
REAL,INTENT(IN)    :: ARX   (MXDIM)
REAL,INTENT(IN)    :: ARY_D (MXDIM)
REAL,INTENT(IN)    :: ARW1  (MXDIM)
REAL,INTENT(IN)    :: ARW2  (MXDIM)

!       Array tangent linear arguments with intent in:
REAL,INTENT(IN)    :: ARY(MXDIM)

!       Array scalars with intent out:
REAL,INTENT(OUT  ) :: XINT_D

!       Tangent linear scalars with intent out:
REAL,INTENT(OUT  ) :: XINT


!     End of subroutine arguments


!       Local scalars:
INTEGER :: I
REAL    :: HINC_D
REAL    :: SS_D
REAL    :: S2_D
REAL    :: S4_D
REAL    :: SY_D
REAL    :: SA_D
REAL    :: SB_D
REAL    :: ARG
REAL    :: H_D
REAL    :: H1_D
REAL    :: SW1
REAL    :: SW2
REAL    :: SY2_D
REAL    :: SY4_D
REAL    :: SYA_D
REAL    :: SYB_D

!      Local tangent linear scalars:
REAL    :: HINC
REAL    :: SS
REAL    :: S2
REAL    :: S4
REAL    :: SY
REAL    :: SA
REAL    :: SB
REAL    :: H



!-----End of header-----------------------------------------------------


!-----Simpsons rule integration pannel width----------------------------
HINC_D = (B_D - A_D)/FLOAT(2*NP)
HINC = (B - A)/FLOAT(2*NP)


!-----Argument for function value at lower integration boundary---------
IF (.NOT. SWITCH) THEN
  CALL RTIASI_INTER(MXDIM,NREC,4,A_D,ARX,ARY_D,SY_D,H_D)
  CALL RTIASI_INTER(MXDIM,NREC,1,A_D,ARX,ARY_D,SYA_D,H1_D)
  CALL RTIASI_INTERTL(MXDIM,NREC,MODY,A_D,ARX,ARY, SY ,SYA_D,H )
  SA_D = SY_D
  SA   = SY
ELSE
  CALL RTIASI_INTER(MXDIM,NREC,MODY,A,ARX,ARY,SY,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW1,A,ARX,ARW1,SW1,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW2,A,ARX,ARW2,SW2,H)
  SA = SY*SW1*SW2
ENDIF
!-----------------------------------------------------------------------


!-----Argument for function value at upper integration boundary---------
IF (.NOT. SWITCH) THEN
  CALL RTIASI_INTER(MXDIM,NREC,4,B_D,ARX,ARY_D,SY_D,H_D)
  CALL RTIASI_INTER(MXDIM,NREC,1,B_D,ARX,ARY_D,SYB_D,H1_D)
  CALL RTIASI_INTERTL(MXDIM,NREC,MODY,B_D,ARX,ARY ,SY ,SYB_D,H )
  SB_D = SY_D
  SB  = SY
ELSE
  CALL RTIASI_INTER(MXDIM,NREC,MODY,B,ARX,ARY,SY,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW1,B,ARX,ARW1,SW1,H)
  CALL RTIASI_INTER(MXDIM,NREC,MOW2,B,ARX,ARW2,SW2,H)
  SB = SY*SW1*SW2
ENDIF
!-----------------------------------------------------------------------


S4_D = 0.
S4   = 0.
DO I=1,NP
  ARG = A_D + HINC_D*FLOAT(2*I - 1) ! Argument for function values at even
                                    ! integration points.

    IF (.NOT. SWITCH) THEN
      CALL RTIASI_INTER(MXDIM,NREC,4,ARG,ARX,ARY_D,SY_D,H_D)
      CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SY4_D,H1_D)
      CALL RTIASI_INTERTL(MXDIM,NREC,MODY,ARG,ARX,ARY ,SY ,SY4_D ,H)
      SS = SY
      SS_D = SY_D
    ELSE
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW1,ARG,ARX,ARW1,SW1,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW2,ARG,ARX,ARW2,SW2,H)
      SS = SY*SW1*SW2
    ENDIF
  S4 = S4 + SS                      ! Sum function values at even
                                    ! integration points.
  S4_D=S4_D+SS_D
ENDDO


S2 = 0.0
S2_D=0
DO I=1,NP-1
  ARG = A_D + HINC_D*FLOAT(2*I)     ! Argument for function values at
                                    ! odd integration points.
    IF (.NOT. SWITCH) THEN
      CALL RTIASI_INTER(MXDIM,NREC,4,ARG,ARX,ARY_D,SY_D,H_D)
      CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SY2_D,H1_D)
      CALL RTIASI_INTERTL(MXDIM,NREC,MODY,ARG,ARX,ARY ,SY ,SY2_D ,H)
      SS = SY
      SS_D=SY_D
    ELSE
      CALL RTIASI_INTER(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW1,ARG,ARX,ARW1,SW1,H)
      CALL RTIASI_INTER(MXDIM,NREC,MOW2,ARG,ARX,ARW2,SW2,H)
      SS = SY*SW1*SW2
    ENDIF
  S2 = S2 + SS                      ! Sum function values at odd
                                    ! integration points.
  S2_D=S2_D+SS_D
ENDDO

XINT_D = HINC_D/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D)
XINT   = HINC/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D) &
      + HINC_D/3.0 * (SA + SB + 4.0*S4 + 2.0*S2)

RETURN
END SUBROUTINE RTIASI_INTEGTL
!     Calculate a  weighted integrated quantity using simpsons rule
!SUBROUTINE INTEGK(MXDIM,NREC,NP,A_D,A,B_D,B,ARX,ARY_D,ARY,MODY, &
!SWITCH,ARW1,MOW1,ARW2,MOW2,XINT_D,XINT)
SUBROUTINE RTIASI_INTEGK(MXDIM,NREC,NP,A_D,A,B_D,B,ARX,ARY_D,ARY,MODY, &
SWITCH,XINT_D,XINT)

!     Description:
!     K of subroutine INTEG
!     To calculate a  weighted integrated quantity using
!     Simpsons rule
!
!     Method:
!     A NP*2 panel Simpsons rule calculation
!     is performed between the limits A and B for the
!     quantity ARY*ARW1*ARW2 and the result returned as
!     XINT. The function values are interpolated using
!     appropriate schemes specified by mode.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            27/02/91    Original Code. D.P. Edwards
!     2            96/11/05    D.P. Edwards
!     3            19/1/1999   Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


 IMPLICIT NONE

!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: MXDIM
INTEGER,INTENT(IN) :: NREC
INTEGER,INTENT(IN) :: NP
INTEGER,INTENT(IN) :: MODY
REAL,INTENT(IN)    :: A_D
REAL,INTENT(IN)    :: B_D
LOGICAL,INTENT(IN) :: SWITCH

!       Scalar tangent linear arguments with intent in:
REAL,INTENT(INOUT)    :: A
REAL,INTENT(INOUT)    :: B

!       Array arguments with intent in:
REAL,INTENT(IN)    :: ARX(MXDIM)
REAL,INTENT(IN)    :: ARY_D(MXDIM)

!       Array tangent linear arguments with intent in:
REAL,INTENT(INOUT)    :: ARY(MXDIM,2)

!       Array scalars with intent out:
REAL,INTENT(OUT  ) :: XINT_D

!       Tangent linear scalars with intent out:
REAL,INTENT(IN) :: XINT


!     End of subroutine arguments


!       Local scalars:
INTEGER :: I
REAL    :: HINC_D
REAL    :: SS_D
REAL    :: S2_D
REAL    :: S4_D
REAL    :: SY_D
REAL    :: SA_D
REAL    :: SB_D
REAL    :: ARG
REAL    :: H_D
REAL    :: SY2_D
REAL    :: SY4_D
REAL    :: SYA_D
REAL    :: SYB_D

!      Local tangent linear scalars:
REAL    :: HINC
REAL    :: SS
REAL    :: S2
REAL    :: S4
REAL    :: SY
REAL    :: SA
REAL    :: SB
REAL    :: H



!-----End of header-----------------------------------------------------



HINC_D = (B_D - A_D)/FLOAT(2*NP)

  IF (.NOT. SWITCH) THEN
    CALL RTIASI_INTER(MXDIM,NREC,4,A_D,ARX,ARY_D,SY_D,H_D)
    SA_D = SY_D
  ENDIF

  IF (.NOT. SWITCH) THEN
    CALL RTIASI_INTER(MXDIM,NREC,4,B_D,ARX,ARY_D,SY_D,H_D)
    SB_D = SY_D
  ENDIF

S4_D = 0.
DO I=1,NP
  ARG = A_D + HINC_D*FLOAT(2*I - 1)
    IF (.NOT. SWITCH) THEN
      CALL RTIASI_INTER(MXDIM,NREC,4,ARG,ARX,ARY_D,SY_D,H_D)
      SS_D = SY_D
    ENDIF
  S4_D=S4_D+SS_D
ENDDO

S2_D=0
  DO I=1,NP-1
    ARG = A_D + HINC_D*FLOAT(2*I)
      IF (.NOT. SWITCH) THEN
        CALL RTIASI_INTER(MXDIM,NREC,4,ARG,ARX,ARY_D,SY_D,H_D)
        SS_D=SY_D
      ENDIF
    S2_D=S2_D+SS_D
ENDDO

XINT_D = HINC_D/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D)

HINC=0.
SA=0.
SB=0.
S4=0.
S2=0.
SS=0.
SY=0.
H=0.

HINC=HINC+XINT/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D)
SA=SA+HINC_D/3.0 *XINT
SB=SB+HINC_D/3.0 *XINT
S4=S4+HINC_D/3.0 * (4.0*XINT)
S2=S2+ HINC_D/3.0 * (2.0*XINT)

DO I=NP-1,1,-1
  ARG = A_D + HINC_D*FLOAT(2*I)     ! Argument for function values at
                                    ! odd integration points.
  SS=SS+S2
    IF (.NOT. SWITCH) THEN
      SY=SY+SS
      SS=0.
      CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SY2_D,H_D)
      CALL RTIASI_INTERK(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,SY2_D,H)
    ENDIF
ENDDO

DO I=NP,1,-1
  ARG = A_D + HINC_D*FLOAT(2*I - 1) ! Argument for function values at even
                                    ! integration points.
  SS=SS+S4
    IF (.NOT. SWITCH) THEN
      SY=SY+SS
      SS=0.
      CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SY4_D,H_D)
      CALL RTIASI_INTERK(MXDIM,NREC,MODY,ARG,ARX,ARY,SY,SY4_D,H)
    ENDIF
ENDDO

IF (.NOT. SWITCH) THEN
  SY=SY+SB
  ARG=B_D
  CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SYB_D,H_D)
  CALL RTIASI_INTERK(MXDIM,NREC,MODY,B_D,ARX,ARY ,SY ,SYB_D,H )
ENDIF

IF (.NOT. SWITCH) THEN
  SY=SY+SA
  ARG=A_D
  CALL RTIASI_INTER(MXDIM,NREC,1,ARG,ARX,ARY_D,SYA_D,H_D)
  CALL RTIASI_INTERK(MXDIM,NREC,MODY,A_D,ARX,ARY ,SY ,SYA_D,H )
ENDIF

B=B+HINC/FLOAT(2*NP)
A=A-HINC/FLOAT(2*NP)

RETURN
END SUBROUTINE RTIASI_INTEGK
!     Interpolation routine
SUBROUTINE RTIASI_INTER(MXDIM,NREC,MODE,ARG,ARX,ARY,SS,H)

!     Description:
!     To interpolate at the point ARG between the array ARX
!     and the corresponding function values ARY.
!
!     Method:
!     Exponential interpolation is used for pressure and
!     number density, linear interpolation for temperature.
!     If the mode number is <0, negative numbers are set = 0.
!
!     Owner:
!     Marco Matricardi.ECMWF
!
!     History:
!     Version      Date        Comment
!     1            15/08/94    Original code. 4.0 D.P. Edwards
!     2            12/11/1999  Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


USE RTIASI_Arrays, ONLY : &
     XL

IMPLICIT NONE


!     Subroutine arguments:

!        Scalars with intent in:
INTEGER,INTENT(IN)  :: MXDIM      ! Array dimension
INTEGER,INTENT(IN)  :: NREC       ! No. of elements in arrays ARX and ARY
INTEGER,INTENT(IN)  :: MODE       ! Interpolation mode
REAL,INTENT(IN)     :: ARG        ! Interpolation argument

!       Scalars with intent out:
REAL,INTENT(OUT)    :: SS         ! Interpolated function value at ARG
REAL,INTENT(OUT)    :: H          ! Gradient or scale height value

!       Arrays with intent in:
REAL,INTENT(IN)     :: ARX(MXDIM) ! X value array
REAL,INTENT(IN)     :: ARY(MXDIM) ! Y value function array

!     End of subroutine arguments


!       Local scalars:
INTEGER :: IR,NO,NOO,KS,KF,K,KL,KD
REAL    :: AA
REAL    :: BB
!REAL    :: XL
REAL    :: BIGNUMBER


!-----End of header----------------------------------------------------


DATA BIGNUMBER/1.0E+20/


IF (ARG >= ARX(1) .AND. ARG <= ARX(NREC)) THEN
  DO IR=1,NREC-1
    IF (ARG >= ARX(IR) .AND. ARG < ARX(IR+1)) THEN
      KL = IR
    END IF
  ENDDO
      IF (ARG == ARX(NREC)) THEN
        KL = NREC - 1
      END IF
ELSEIF (ARG < ARX(1)) THEN
  KL = 1
ELSEIF (ARG > ARX(NREC)) THEN
  KL = NREC - 1
ENDIF

!-----INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KL TO KL+1-----


!-----Exponential interpolation-----------------------------------------

IF (ABS(MODE) == 1. .OR. ABS(MODE) == 4.) THEN
      IF (ARY(KL+1) == ARY(KL)) THEN
        H = BIGNUMBER
        SS=ARY(KL)
      ELSE
        H = -(ARX(KL+1) - ARX(KL))/LOG(ARY(KL+1)/ARY(KL))
        SS = ARY(KL)*EXP(-(ARG-ARX(KL))/H)
      ENDIF
    IF (ABS(MODE) == 4.) THEN
      SS = 1.0/SS
    END IF
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF

!-----Linear interpolation----------------------------------------------

ELSEIF (ABS(MODE) == 2 .OR. ABS(MODE) == 5) THEN
  H = (ARY(KL+1) - ARY(KL))/(ARX(KL+1) - ARX(KL))
  SS = ARY(KL) + H*(ARG - ARX(KL))
    IF (ABS(MODE) == 5) THEN
      SS = 1.0/SS
    END IF
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF

!-----Linear-log interpolation------------------------------------------

ELSEIF (ABS(MODE) == 7) THEN
  H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(ARX(KL+1)-ARX(KL))
  SS = EXP(LOG(ARY(KL)) + H*(ARG - ARX(KL)))
    IF (MODE < 0 .AND. SS < 0.0) THEN
      SS = 0.0
    END IF

!-----Log-log interpolation---------------------------------------------

ELSEIF (ABS(MODE) == 8) THEN
  H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(LOG(ARX(KL+1))-LOG(ARX(KL)))
  SS = EXP(LOG(ARY(KL)) + H*(LOG(ARG) - LOG(ARX(KL))))
    IF (MODE < 0 .AND. SS < 0.0) THEN
      SS = 0.0
    END IF

!-----Logarithmic interpolation-----------------------------------------

ELSEIF (ABS(MODE) == 3) THEN
  AA = ARX(KL+1)/ARX(KL)
  BB = ARG/ARX(KL)
    IF (AA == BB) THEN
      SS = ARY(KL+1)
    ELSE
      H = (ARY(KL+1) - ARY(KL))/LOG(AA)
      SS = ARY(KL) + H*LOG(BB)
    ENDIF
      IF (MODE < 0 .AND. SS < 0.0)THEN
        SS = 0.0
      END IF

!-----Lagrangian interpolation------------------------------------------

ELSEIF (ABS(MODE) == 6) THEN

  NO = 4      ! Number of data point points to interpolate between

!-----FIND DATA POINTS BETWEEN WHICH TO INTERPOLATE---------------------

  NOO = NO


  DO KD=1,10000
    IF (ARG < ARX(1)) THEN
      NOO = 2
      KS = 1
      KF = 1 + NOO - 1
    ELSEIF (ARG > ARX(NREC)) THEN
      NOO = 2
      KS = NREC - NOO + 1
      KF = NREC
    ELSE
      IF (MOD(NOO,2) == 0) THEN
        KS = KL - 0.5*NOO + 1
        KF = KL + 0.5*NOO
      ELSE
        KS = KL - 0.5*(NOO - 1) + 1
        KF = KL + 0.5*(NOO + 1)
      ENDIF
        IF (KS < 1) THEN
          KS = 1
        END IF
          IF (KF > NREC) THEN
            KF = NREC
          END IF
    ENDIF

!-----INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KS TO KF-------

  SS = 0.0
    DO K=KS,KF
      SS = SS + XL(MXDIM,KS,KF,K,ARG,ARX)*ARY(K) &
      /XL(MXDIM,KS,KF,K,ARX(K),ARX)
    ENDDO
  H = NOO

!-----IF INTERPOLATION HAS OVERSHOT, REDUCE ORDER AND TRY AGAIN---------

  IF (ARG < ARX(1) .OR. ARG > ARX(NREC)) THEN
    IF (SS < 0.0) THEN
      IF (NOO == 2) THEN
        SS = 0.0
        EXIT
      ELSE
        NOO = NOO - 1
        CYCLE
      ENDIF
    ENDIF
  ELSE
    IF (((ARY(KL) <= ARY(KL+1)) .AND. &
    (SS < ARY(KL) .OR. SS > ARY(KL+1))) .OR. &
    ((ARY(KL) > ARY(KL+1)) .AND. &
    (SS > ARY(KL) .OR. SS < ARY(KL+1)))) THEN
      NOO = NOO-1
      CYCLE
    ENDIF
  ENDIF
END DO
ENDIF

!-----------------------------------------------------------------------
 RETURN
 END SUBROUTINE RTIASI_INTER
!     Interpolation routine
SUBROUTINE RTIASI_INTERTL(MXDIM,NREC,MODE,ARG,ARX,ARY,SS,SS_D,H)

!     Description:
!     TL of subroutine INTER
!     To interpolate at the point ARG between the array ARX
!     and the corresponding function values ARY.
!
!     Method:
!     Exponential interpolation is used for pressure and
!     number density, linear interpolation for temperature.
!     If the mode number is <0, negative numbers are set = 0.
!
!     Owner:
!     Marco Matricardi.ECMWF
!
!     History:
!     Version      Date        Comment
!     1            15/08/94    Original code. 4.0 D.P. Edwards
!     2            12/11/1999  Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

USE RTIASI_Arrays, ONLY : &
     XL

IMPLICIT NONE


!     Subroutine arguments:

!        Scalars with intent in:
INTEGER,INTENT(IN)  :: MXDIM      ! Array dimension
INTEGER,INTENT(IN)  :: NREC       ! Number of elements in arrays ARX and ARY
INTEGER,INTENT(IN)  :: MODE       ! Interpolation mode
REAL   ,INTENT(IN)  :: ARG        ! Interpolation argument
REAL   ,INTENT(IN)  :: SS_D

!       Scalars with intent out:
REAL,INTENT(OUT)    :: SS         ! Interpolated function value at ARG
REAL,INTENT(OUT)    :: H          ! Gradient or scale height value

!       Arrays with intent in:
REAL,INTENT(IN)     :: ARX(MXDIM) !X value array
REAL,INTENT(IN)     :: ARY(MXDIM) !Y value function array

!     End of subroutine arguments


!       Local scalars:
INTEGER :: IR,NO,NOO,KS,KF,K,KL,KD
REAL :: AA
REAL :: BB
!REAL :: XL
REAL :: BIGNUMBER


!-----End of header----------------------------------------------------


DATA BIGNUMBER/1.0E+20/


IF (ARG >= ARX(1) .AND. ARG <= ARX(NREC)) THEN
  DO IR=1,NREC-1
    IF (ARG >= ARX(IR) .AND. ARG < ARX(IR+1)) THEN
      KL = IR
    END IF
  ENDDO
      IF (ARG == ARX(NREC)) THEN
        KL = NREC - 1
      END IF
ELSEIF (ARG < ARX(1)) THEN
  KL = 1
ELSEIF (ARG > ARX(NREC)) THEN
  KL = NREC - 1
ENDIF

!-----INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KL TO KL+1-----


!-----Exponential interpolation-----------------------------------------

IF (ABS(MODE) == 1. .OR. ABS(MODE) == 4.) THEN
      IF (ARY(KL+1) == ARY(KL)) THEN
        H = BIGNUMBER
        SS=ARY(KL)
      ELSE
        H = -(ARX(KL+1) - ARX(KL))/LOG(ARY(KL+1)/ARY(KL))
        SS = ARY(KL)*EXP(-(ARG-ARX(KL))/H)
      ENDIF
    IF (ABS(MODE) == 4.) THEN
      SS = 1.0/SS
    END IF
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF

!-----Linear interpolation----------------------------------------------

ELSEIF (ABS(MODE) == 2 .OR. ABS(MODE) == 5) THEN
  H = (ARY(KL+1) - ARY(KL))/(ARX(KL+1) - ARX(KL))
  SS = ARY(KL) + H*(ARG - ARX(KL))
    IF (ABS(MODE) == 5) THEN
      SS=-SS/(SS_D**2)
    END IF
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF

!-----Linear-log interpolation------------------------------------------

ELSEIF (ABS(MODE) == 7) THEN
  H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(ARX(KL+1)-ARX(KL))
  SS = EXP(LOG(ARY(KL)) + H*(ARG - ARX(KL)))
    IF (MODE < 0 .AND. SS < 0.0) THEN
      SS = 0.0
    END IF

!-----Log-log interpolation---------------------------------------------

ELSEIF (ABS(MODE) == 8) THEN
  H=(LOG(ARY(KL+1))-LOG(ARY(KL)))/(LOG(ARX(KL+1))-LOG(ARX(KL)))
  SS = EXP(LOG(ARY(KL)) + H*(LOG(ARG) - LOG(ARX(KL))))
    IF (MODE < 0 .AND. SS < 0.0) THEN
      SS = 0.0
    END IF

!-----Logarithmic interpolation-----------------------------------------

ELSEIF (ABS(MODE) == 3) THEN
  AA = ARX(KL+1)/ARX(KL)
  BB = ARG/ARX(KL)
    IF (AA == BB) THEN
      SS = ARY(KL+1)
    ELSE
      H = (ARY(KL+1) - ARY(KL))/LOG(AA)
      SS = ARY(KL) + H*LOG(BB)
    ENDIF
      IF (MODE < 0 .AND. SS < 0.0)THEN
        SS = 0.0
      END IF

!-----Lagrangian interpolation------------------------------------------

ELSEIF (ABS(MODE) == 6) THEN

  NO = 4      ! Number of data point points to interpolate between

!-----FIND DATA POINTS BETWEEN WHICH TO INTERPOLATE---------------------

  NOO = NO


  DO KD=1,10000
    IF (ARG < ARX(1)) THEN
      NOO = 2
      KS = 1
      KF = 1 + NOO - 1
    ELSEIF (ARG > ARX(NREC)) THEN
      NOO = 2
      KS = NREC - NOO + 1
      KF = NREC
    ELSE
      IF (MOD(NOO,2) == 0) THEN
        KS = KL - 0.5*NOO + 1
        KF = KL + 0.5*NOO
      ELSE
        KS = KL - 0.5*(NOO - 1) + 1
        KF = KL + 0.5*(NOO + 1)
      ENDIF
        IF (KS < 1) THEN
          KS = 1
        END IF
          IF (KF > NREC) THEN
            KF = NREC
          END IF
    ENDIF

!-----INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KS TO KF-------

  SS = 0.0
    DO K=KS,KF
      SS = SS + XL(MXDIM,KS,KF,K,ARG,ARX)*ARY(K) &
      /XL(MXDIM,KS,KF,K,ARX(K),ARX)
    ENDDO
  H = NOO

!-----IF INTERPOLATION HAS OVERSHOT, REDUCE ORDER AND TRY AGAIN---------

  IF (ARG < ARX(1) .OR. ARG > ARX(NREC)) THEN
    IF (SS < 0.0) THEN
      IF (NOO == 2) THEN
        SS = 0.0
        EXIT
      ELSE
        NOO = NOO - 1
        CYCLE
      ENDIF
    ENDIF
  ELSE
    IF (((ARY(KL) <= ARY(KL+1)) .AND. &
    (SS < ARY(KL) .OR. SS > ARY(KL+1))) .OR. &
    ((ARY(KL) > ARY(KL+1)) .AND. &
    (SS > ARY(KL) .OR. SS < ARY(KL+1)))) THEN
      NOO = NOO-1
      CYCLE
    ENDIF
  ENDIF
END DO
ENDIF

!-----------------------------------------------------------------------
 RETURN
 END SUBROUTINE RTIASI_INTERTL
!     Compute Lagrangian interpolation coefficients
FUNCTION XL(MXDIM,KS,KF,K,ARG,ARR)

!     Description:
!     XL - To compute Lagrangian interpolation coefficients
!
!     Method:
!     Simple
!
!     Owner:
!     Marco Matricardi.ECMWF
!
!     History:
!     Version      Date        Comment
!     1            01/01/89    Original code. D.P. Edwards.
!     2            12/11/1999  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code"


IMPLICIT NONE



!     Function arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: MXDIM       ! Array dimension.
INTEGER,INTENT(IN) :: KS          ! Lower limit of lagrangian sum.
INTEGER,INTENT(IN) :: KF          ! Upper limit of lagrangian sum.
REAL   ,INTENT(IN) :: ARG         ! Interpolation argument.

!       Array arguments with intent in:
REAL   ,INTENT(IN) :: ARR(MXDIM)  ! Array to interpolate between.

!     End of function arguments


!       Local scalars
INTEGER :: J
INTEGER :: K
REAL    :: PROD
REAL    :: XL


!-----End of header----------------------------------------------------

PROD = 1.0
  DO J=KS,KF
    IF (J /= K) THEN
      PROD = PROD*(ARG - ARR(J))
    ENDIF
  ENDDO
XL = PROD

RETURN
END FUNCTION XL





!     Interpolation routine
SUBROUTINE RTIASI_INTERK(MXDIM,NREC,MODE,ARG,ARX,ARY,SS,SS_D,H)

!     Description:
!     To interpolate at the point ARG between the array ARX
!     and the corresponding function values ARY.
!
!     Method:
!     Exponential interpolation is used for pressure and
!     number density, linear interpolation for temperature.
!     If the mode no is negative numbers are set = 0.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            15/08/94    Original code. 4.0 D.P. Edwards
!     2            96/11/05    D.P. Edwards
!     3            19/1/1999   Marco Matricardi. ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


IMPLICIT NONE


!     Subroutine arguments:

!        Scalars with intent in:
INTEGER,INTENT(IN)  :: MXDIM
INTEGER,INTENT(IN)  :: NREC
INTEGER,INTENT(IN)  :: MODE
REAL,INTENT(IN)     :: ARG
REAL,INTENT(IN)     :: SS_D

!       Scalars with intent out:
REAL,INTENT(INOUT)    :: SS
REAL,INTENT(INOUT)    :: H

!       Arrays with intent in:
REAL,INTENT(IN)     :: ARX(MXDIM)
REAL,INTENT(INOUT)  :: ARY(MXDIM)

!     End of subroutine arguments


!       Local scalars:
INTEGER :: IR,KL
REAL    :: BIGNUMBER



!-----End of header----------------------------------------------------


DATA BIGNUMBER/1.0E+20/

!      ARY(:)=0.
!      H=0.

IF (ARG >= ARX(1) .AND. ARG <= ARX(NREC)) THEN
  DO IR=1,NREC-1
    IF (ARG >= ARX(IR) .AND. ARG < ARX(IR+1)) THEN
      KL = IR
    END IF
  ENDDO
      IF (ARG == ARX(NREC)) THEN
        KL = NREC - 1
      END IF
ELSEIF (ARG < ARX(1)) THEN
  KL = 1
ELSEIF (ARG > ARX(NREC)) THEN
  KL = NREC - 1
ENDIF

!-----INTERPOLATE FUNCTION VALUE AT ARG FROM DATA POINTS KL TO KL+1-----


!-----Exponential interpolation-----------------------------------------

IF (ABS(MODE) == 1. .OR. ABS(MODE) == 4.) THEN
      IF (ARY(KL+1) == ARY(KL)) THEN
        H = BigNumber
        SS=ARY(KL)
      ELSE
        H = -(ARX(KL+1) - ARX(KL))/LOG(ARY(KL+1)/ARY(KL))
        SS = ARY(KL)*EXP(-(ARG-ARX(KL))/H)
      ENDIF
    IF (ABS(MODE) == 4.) THEN
      SS = 1.0/SS
    END IF
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF

!-----Linear interpolation----------------------------------------------

ELSEIF (ABS(MODE) == 2 .OR. ABS(MODE) == 5) THEN
      IF (MODE < 0 .AND. SS < 0.0) THEN
        SS = 0.0
      END IF
    IF (ABS(MODE) == 5) THEN
      SS = -SS/(SS_D**2)
    ENDIF

  ARY(KL)=ARY(KL)+SS
  H=H+SS*(ARG - ARX(KL))
  ARY(KL+1)=ARY(KL+1)+H/(ARX(KL+1) - ARX(KL))
  ARY(KL)=ARY(KL)-H/(ARX(KL+1) - ARX(KL))

  H=0.
  SS=0.

ENDIF

!-----------------------------------------------------------------------
 RETURN
 END SUBROUTINE RTIASI_INTERK




!     Produce Curtis-Godson weighted average layer  parameters
SUBROUTINE RTIASI_LAYERS(NPROF,ILEVWRT,PANGL,KNPF)

!     Description:
!     Produce Curtis-Godson weighted average layer  parameters
!
!     Method:
!     The atmosphere is divided in a number of layers
!     (currently 42) and the input  atmospheric profile data
!     is used to produce Curtis-Godson weighted average
!     layer temperatures
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            96/07/15    Version 4.0  D.P. Edwards
!     2            12/11/1999  To compute Curtis Godson weighted average layer
!                              temperatures to be used by RTIASI in radiance calculations.
!                              calculations. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:

USE RTIASI_PRFVAR,ONLY : &
!     Imported arrays:
  TEMPW   ,   & ! Curtis-Godson Temperature profile in K
  TEMP    ,   & ! Temperature profile in K
  WMIX    ,   & ! Water vapour volume mixing ratio in ppmv
  TA      ,   & ! Surface air temperature in K
  WMIXS   ,   & ! Water vapour surface volume mixing ratio in ppmv
  SURFP         ! Surface pressure in mb=hPa


USE RTIASI_PRFCON,ONLY : &
!     Imported arrays:
  XPRES       ! Standard pressure levels for transmittance (and, currently,
!                   ! radiative treansfer) calculation; from top down; in mb=hpa


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF

!       Array arguments with intent in:
REAL   ,INTENT(IN) :: PANGL(KNPF)

!     End of subroutine arguments


!       Local parameters:
INTEGER,PARAMETER :: MAXLEV=50

!       Local scalars:
INTEGER           :: I,IR,NPTS,IPATH
INTEGER           :: NP
INTEGER           :: NR
INTEGER           :: NOLAY
INTEGER           :: ILEVWRT
INTEGER           :: NPROF
REAL              :: EQRAD           ! Earth equatoral radius
REAL              :: ALPHA           ! Earth excentricity
REAL              :: RTOD            ! Radian to degree conversion factor
REAL              :: XLAT            ! Latitude
REAL              :: ZERHT           ! Eight of first pressure level
REAL              :: REARTH          ! Earth radius at XLAT
REAL              :: ARG
REAL              :: DIF
REAL              :: SS
REAL              :: A,B,TP

!       Local arrays:
REAL              :: XHL    (MAXLEV) ! Layer lower height [km]
REAL              :: XHU    (MAXLEV) ! Layer upper height [km]
REAL              :: XPL    (MAXLEV) ! Layer lower pressure [mb]
REAL              :: XPU    (MAXLEV) ! Layer upper pressure [mb]
REAL              :: SINREF (MAXLEV) ! Sine of local zenith angle
REAL              :: ALT    (MAXLEV)
REAL              :: PRES   (MAXLEV) ! Profile level pressures [mb]
REAL              :: AIRNO  (MAXLEV)
REAL              :: WATER  (MAXLEV) ! Profile volume mixing ratio [ppmv]
REAL              :: DUS    (MAXLEV) ! Data point air densities [cm-3]
REAL              :: HUS    (MAXLEV) ! Height of pressure levels.
REAL              :: TEMPP  (MAXLEV) ! Temperature profile [K]


!-----End of header-----------------------------------------------------


DATA EQRAD /6378.388/
DATA ALPHA /3.367E-3/
DATA RTOD  /0.01745326/
DATA XLAT  /45./
DATA ZERHT /0./
DATA NP    /10/


NR=ILEVWRT
NOLAY=NR-1


!-----UNPACK INPUT PROFILE----------------------------------------------
DO I=1,NR
  IF(I == 1)THEN
    PRES(I)=SURFP(NPROF)
    WATER(I)=WMIXS(NPROF)
    TEMPP(I)=TA(NPROF)
  ELSE
    PRES(I)=XPRES(NR-I+1)
    WATER(I)=WMIX(NR-I+1,NPROF)
    TEMPP(I)=TEMP(NR-I+1,NPROF)
  END IF
ENDDO
!-----------------------------------------------------------------------


!-----SET UP LAYER LOWER AND UPPER PRESSURES----------------------------
DO I=1,NOLAY
  IF(I == 1)THEN
    XPL(I)=SURFP(NPROF)
    XPU(I)=XPRES(NR-I)
  ELSE
    XPL(I)=XPRES(NR-I+1)
    XPU(I)=XPRES(NR-I)
  END IF
ENDDO
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     CALCULATE THE EARTH RADIUS AT THE THIS VALUE OF XLAT
!-----------------------------------------------------------------------

!-----A reference value of 45 degrees is assumed for the latitude-------


DIF = (1.0 - ALPHA)**2
ARG = RTOD*ABS(XLAT)
REARTH = SQRT(EQRAD**2*DIF/(SIN(ARG)**2 + DIF*COS(ARG)**2))

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!     CALCULATE AIR DENSITY
!-----------------------------------------------------------------------

DO IR=1,NR
  DUS(IR)=6.022E+26*1.E-06*1.E+02*PRES(IR)/(8.3143E+03*TEMPP(IR))
END DO

!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
!     CALCULATE HEIGHTS OF INPUT PRESSURE LEVELS
!-----------------------------------------------------------------------

CALL RTIASI_PRESHT(NR,NP,XLAT,ZERHT,REARTH,PRES,DUS,WATER,HUS,MAXLEV)

!-----------------------------------------------------------------------

NPTS = NR

DO IR=1,NPTS
  ALT(IR) =  HUS(IR)
  AIRNO(IR) = DUS(IR)
END DO

DO I=1,NOLAY
  XHL(I)=ALT(I)
  XHU(I)=ALT(I+1)
  SINREF(I)=SIN(PANGL(NPROF))
END DO




!-----------------------------------------------------------------------
!     CALCULATE CURTIS-GODSON WEIGHTED GAS PARAMETERS FOR EACH LAYER
!-----------------------------------------------------------------------


DO IPATH=1,NOLAY
  SS = SINREF(IPATH)
  A  = XHL(IPATH)
  B  = XHU(IPATH)
  CALL RTIASI_RAYTCE(NPTS,NP,A,B,ALT,TEMPP,PRES,AIRNO,REARTH,SS, &
              TP,MAXLEV)

  TEMPW(-IPATH+1+NOLAY,NPROF)=TP
END DO

!-----------------------------------------------------------------------

 RETURN
 END SUBROUTINE RTIASI_LAYERS
!     Produce Curtis-Godson weighted average layer  parameters
SUBROUTINE RTIASI_LAYERSK(PANGL,KNPF,KNCHPF,KPROF)

!     Description:
!     K of subroutine LAYERS
!     Produce Curtis-Godson weighted average layer  parameters
!
!     Method:
!     The atmosphere is divided in a number of layers
!     (currently 42) and the input  atmospheric profile data
!     is used to produce Curtis-Godson weighted average
!     layer temperatures
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            12/11/1999  To compute Curtis Godson weighted average layer
!                              temperatures to be used by RTIASI in radiance
!                              calculations. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Direct Module used:

USE RTIASI_PRFVAR,ONLY : &
!     Imported arrays:
  TEMP_D   =>TEMP   ,  & ! Temperature profile in K
  WMIX_D   =>WMIX   ,  & ! Water vapour volume mixing ratio in ppmv
  TA_D     =>TA     ,  & ! Surface air temperature in K
  WMIXS_D  =>WMIXS  ,  & ! Water vapour surface volume mixing ratio in ppmv
  SURFP_D  =>SURFP       ! Surface pressure in mb=hPa


USE RTIASI_PRFCON,ONLY : &
!     Imported arrays:
  XPRES                ! Standard pressure levels for transmittance (and,
                       ! currently,radiative transfer) calculation;
                       ! from top down; in mb=hpa


USE RTIASI_SURF, ONLY : &
!     Imported arrays:
  NLEVSF_D =>NLEVSF ,  & ! Index of nearest standard pressure level
                       ! at/below surface.
  FRACPS_D =>FRACPS    ! Fraction of standard pressure level interval by
                       ! which surf is above level nlevsf.

!     K module used:

USE RTIASI_PRFVARK,ONLY : &
!     Imported arrays:
  TEMPW             ,  & ! Curtis-Godson Temperature profile in K
  TEMP              ,  & ! Temperature profile in K
  WMIX              ,  & ! Water vapour volume mixing ratio in ppmv
  TA                ,  & ! Surface air temperature in K
  WMIXS             ,  & ! Water vapour surface volume mixing ratio in ppmv
  SURFP                  ! Surface pressure in mb=hPa


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF
INTEGER,INTENT(IN) :: KNCHPF

!       Array arguments with intent in:
REAL      ,INTENT(IN) :: PANGL(KNPF)
INTEGER   ,INTENT(IN) :: KPROF(KNCHPF)

!     End of subroutine arguments


!       Local parameters:
INTEGER,PARAMETER :: MAXLEV=50

!       Local scalars:
INTEGER           :: I,J,IR,NPTS,IPATH
INTEGER           :: ILEVWRT
INTEGER           :: NP
INTEGER           :: NR
INTEGER           :: NOLAY
REAL              :: EQRAD          ! Earth equatoral radius
REAL              :: ALPHA          ! Earth excentricity
REAL              :: RTOD           ! Radian to degree conversion factor
REAL              :: XLAT           ! Latitude
REAL              :: ZERHT          ! Eight of first pressure level
REAL              :: REARTH         ! Earth radius at XLAT
REAL              :: ARG
REAL              :: DIF

!       Tangent linear local scalars:
REAL              :: A(MAXLEV,KNCHPF),B(MAXLEV,KNCHPF)

!       Local arrays:
REAL              :: XHL_D    (MAXLEV,KNPF) ! Layer lower height [km]
REAL              :: XHU_D    (MAXLEV,KNPF) ! Layer upper height [km]
REAL              :: XPL_D    (MAXLEV,KNPF) ! Layer lower pressure [mb]
REAL              :: XPU_D    (MAXLEV,KNPF) ! Layer upper pressure [mb]
REAL              :: SINREF   (MAXLEV,KNPF) ! Sine of local zenith angle
REAL              :: ALT_D    (MAXLEV,KNPF)
REAL              :: AIRNO_D  (MAXLEV,KNPF)
REAL              :: WATER_D  (MAXLEV,KNPF) ! Vol. mixing ratio [ppmv]
REAL              :: DUS_D    (MAXLEV,KNPF) ! Air densities [cm-3]
REAL              :: HUS_D    (MAXLEV,KNPF) ! Height of pressure levels.
REAL              :: HUS_DD   (MAXLEV)
REAL              :: TEMPP_D  (MAXLEV,KNPF) ! Temperature profile [K]
REAL              :: PRES_D   (MAXLEV,KNPF) ! Pressure values.
REAL              :: A_D      (MAXLEV,KNPF)
REAL              :: B_D      (MAXLEV,KNPF)
REAL              :: TP       (MAXLEV,KNCHPF)


!      K local arrays:
REAL              :: PRES     (MAXLEV,KNCHPF) ! Pressure values.
REAL              :: WATER    (MAXLEV,KNCHPF) ! volume mixing ratio [ppmv]
REAL              :: TEMPP    (MAXLEV,KNCHPF) ! Temperature profile [K]
REAL              :: XPL      (MAXLEV,KNCHPF) ! Layer lower pressure [mb]
REAL              :: XPU      (MAXLEV,KNCHPF) ! Layer upper pressure [mb]
REAL              :: DUS      (MAXLEV,KNCHPF) ! Air densities [cm-3]
REAL              :: XHL      (MAXLEV,KNCHPF) ! Layer lower height [km]
REAL              :: XHU      (MAXLEV,KNCHPF) ! Layer upper height [km]
REAL              :: ALT      (MAXLEV,KNCHPF)
REAL              :: AIRNO    (MAXLEV,KNCHPF)
REAL              :: HUS      (MAXLEV,KNCHPF) ! Height of pressure levels.


!-----End of header-----------------------------------------------------


DATA EQRAD      /6378.388/
DATA ALPHA      /3.367E-3/
DATA RTOD       /0.01745326/
DATA XLAT       /45./
DATA ZERHT      /0./
DATA NP         /10/

TP(:,:)    =0.
A(:,:)     =0.
B(:,:)     =0.
TEMPP(:,:) =0.
PRES(:,:)  =0.
AIRNO(:,:) =0.


!-----------------------------------------------------------------------
!     CALCULATE THE EARTH RADIUS AT THE THIS VALUE OF XLAT
!-----------------------------------------------------------------------

!-----A reference value of 45 degrees is assumed for the latitude-------

DIF = (1.0 - ALPHA)**2
ARG = RTOD*ABS(XLAT)
REARTH = SQRT(EQRAD**2*DIF/(SIN(ARG)**2 + DIF*COS(ARG)**2))

!-----------------------------------------------------------------------



!-----REPEAT DIRECT CALCULATIONS----------------------------------------

DO J=1,KNPF
  ILEVWRT=NLEVSF_D(J)
  IF (FRACPS_D(J) < 0.) THEN
    ILEVWRT=ILEVWRT+1
  END IF

  NR=ILEVWRT
  NOLAY=NR-1


!-----UNPACK INPUT PROFILE----------------------------------------------

  DO I=1,NR
    IF(I == 1)THEN
      PRES_D  (I,J) =SURFP_D (J)
      WATER_D (I,J) =WMIXS_D (J)
      TEMPP_D (I,J) =TA_D    (J)
    ELSE
      PRES_D  (I,J) =XPRES   (NR-I+1)
      WATER_D (I,J) =WMIX_D  (NR-I+1,J)
      TEMPP_D (I,J) =TEMP_D  (NR-I+1,J)
    END IF
  ENDDO
!-----------------------------------------------------------------------

!-----SET UP LAYER LOWER AND UPPER PRESSURES----------------------------

  DO I=1,NOLAY
    IF(I == 1)THEN
      XPL_D (I,J) =SURFP_D (J)
      XPU_D (I,J) =XPRES   (NR-I)
    ELSE
      XPL_D (I,J) =XPRES   (NR-I+1)
      XPU_D (I,J) =XPRES   (NR-I)
    END IF
  ENDDO
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!     CALCULATE AIR DENSITY
!-----------------------------------------------------------------------

  DO IR=1,NR
    DUS_D(IR,J)=6.022E+26*1.E-06*1.E+02*PRES_D(IR,J)/ &
            (8.3143E+03*TEMPP_D(IR,J))
  ENDDO

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     CALCULATE HEIGHTS OF INPUT PRESSURE LEVELS
!-----------------------------------------------------------------------

  CALL RTIASI_PRESHT(NR,NP,XLAT,ZERHT,REARTH,PRES_D(:,J),DUS_D(:,J), &
            WATER_D(:,J),HUS_DD,MAXLEV)


  HUS_D   (:,J)=HUS_DD   (:)

  NPTS = NR

  DO IR=1,NPTS
    ALT_D   (IR,J)   = HUS_D(IR,J)
    AIRNO_D (IR,J)   = DUS_D(IR,J)
  END DO

  DO I=1,NOLAY
    XHL_D  (I,J)    = ALT_D     (I,J)
    XHU_D  (I,J)    = ALT_D     (I+1,J)
    SINREF (I,J)    = SIN       (PANGL (J))
  END DO


!-----------------------------------------------------------------------
!     CALCULATE CURTIS-GODSON WEIGHTED GAS PARAMETERS FOR EACH LAYER
!-----------------------------------------------------------------------

  DO IPATH=1,NOLAY
    A_D(IPATH,J)  = XHL_D(IPATH,J)
    B_D(IPATH,J)  = XHU_D(IPATH,J)
  ENDDO
!-----------------------------------------------------------------------

ENDDO



!-----K CALCULATIONS----------------------------------------------------

!-----------------------------------------------------------------------
!     CALCULATE CURTIS-GODSON WEIGHTED GAS PARAMETERS FOR EACH LAYER
!-----------------------------------------------------------------------

DO J=1,KNCHPF
  ILEVWRT=NLEVSF_D(KPROF(J))

  IF (FRACPS_D(KPROF(J)) < 0.) THEN
    ILEVWRT=ILEVWRT+1
  END IF

  NR=ILEVWRT
  NOLAY=NR-1


  DO IPATH=NOLAY,1,-1
    TP(IPATH,J)=TEMPW(-IPATH+1+NOLAY,J)
    TEMPW(-IPATH+1+NOLAY,J)=0.
  ENDDO

ENDDO


CALL RTIASI_RAYTCEK(NPTS,NP,A_D,A,B_D, &
       B,ALT_D,TEMPP_D,TEMPP, &
       PRES_D,PRES,AIRNO_D,AIRNO, &
       REARTH,SINREF,TP,MAXLEV,NOLAY,KNCHPF,KNPF,KPROF)


WATER(:,:) =0
HUS(:,:)   =0
DUS(:,:)   =0
ALT(:,:)   =0
XHU(:,:)   =0
XHL(:,:)   =0
XPL(:,:)   =0


DO J=1,KNCHPF

  ILEVWRT=NLEVSF_D(KPROF(J))

  IF (FRACPS_D(KPROF(J)) < 0.) THEN
    ILEVWRT=ILEVWRT+1
  END IF

  NR=ILEVWRT
  NOLAY=NR-1

  DO IPATH=NOLAY,1,-1
    XHU(IPATH,J)=XHU(IPATH,J)+B(IPATH,J)
    XHL(IPATH,J)=XHL(IPATH,J)+A(IPATH,J)
  ENDDO


  DO I=NOLAY,1,-1
    ALT(I+1,J)=ALT(I+1,J)+XHU(I,J)
    XHU(I,J)=0.
    ALT(I,J)  =ALT(I,J)  +XHL(I,J)
    XHL(I,J)=0.
  ENDDO

  NPTS = NR

  DO IR=NPTS,1,-1
    DUS(IR,J)=DUS(IR,J)+AIRNO(IR,J)
    AIRNO(IR,J)=0.
    HUS(IR,J)=HUS(IR,J)+ALT(IR,J)
    ALT(IR,J)=0.
  ENDDO

ENDDO



!-----------------------------------------------------------------------
!     CALCULATE HEIGHTS OF INPUT PRESSURE LEVELS
!-----------------------------------------------------------------------

CALL RTIASI_PRESHTK(NR,NP,XLAT,ZERHT,REARTH,PRES_D,PRES, &
     DUS_D,DUS,WATER_D,WATER, &
     HUS_D,HUS,MAXLEV,KNPF,KNCHPF,KPROF)



HUS(:,:)=0.



DO J=1,KNCHPF
  ILEVWRT=NLEVSF_D(KPROF(J))

  IF (FRACPS_D(KPROF(J)) < 0.) THEN
    ILEVWRT=ILEVWRT+1
  END IF

  NR=ILEVWRT
  NOLAY=NR-1


!-----------------------------------------------------------------------
!     CALCULATE AIR DENSITY
!-----------------------------------------------------------------------

  DO IR=NR,1,-1
    PRES(IR,J)=PRES(IR,J)+6.022E+26*1.E-06*1.E+02*DUS(IR,J)/ &
           (8.3143E+03*TEMPP_D(IR,KPROF(J)))
    TEMPP(IR,J)=TEMPP(IR,J)-(6.022E+26*1.E-06*1.E+02* &
    PRES_D(IR,KPROF(J))/(8.3143E+03*TEMPP_D(IR,KPROF(J))**2))* &
    DUS(IR,J)
    DUS(IR,J)=0.
  ENDDO

!-----SET UP LAYER LOWER AND UPPER PRESSURES----------------------------

  DO I=NOLAY,1,-1
    IF(I == 1)THEN
      XPU(I,J)=0.
      SURFP(J)=SURFP(J)+XPL(I,J)
      XPL(I,J)  =0.
    ELSE
      XPL(I,J)  =0.
      XPU(I,J)  =0.
    ENDIF
  ENDDO


!-----UNPACK INPUT PROFILE----------------------------------------------

  DO I=NR,1,-1
    IF(I == 1)THEN
      TA(J)=TA(J)+TEMPP(I,J)
      TEMPP(I,J)=0.
      WMIXS(J)=WMIXS(J)+WATER(I,J)
      WATER(I,J)=0.
      SURFP(J)=SURFP(J)+PRES(I,J)
      PRES(I,J)=0.
    ELSE
      TEMP(NR-I+1,J)=TEMP(NR-I+1,J)+TEMPP(I,J)
      TEMPP(I,J)=0.
      WMIX(NR-I+1,J)=WMIX(NR-I+1,J)+WATER(I,J)
      WATER(I,J)=0.
    END IF
  ENDDO
!-----------------------------------------------------------------------

ENDDO


RETURN
END SUBROUTINE RTIASI_LAYERSK
!     Calculate the heights of pressure levels
SUBROUTINE RTIASI_PRESHT(NR,NP,XLAT,ZERHT,REARTH,PUS,DUS,GASUS,HUS,MAXLEV &
)

!     Description:
!     To calculate the heights of pressure levels
!
!     Method.
!     A calculation is performed based on the hydrostatic
!     equation to obtain the heights of pressure levels.
!     virtual temperatures are used to account for the
!     presence of water vapour. the variation of gravity
!     with latitude and height is also considered.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            95/08/16    Version 4.0  D.P. Edwards
!     2            12/11/1999  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


 IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER ,INTENT(IN) :: NR
INTEGER ,INTENT(IN) :: NP
INTEGER ,INTENT(IN) :: MAXLEV
REAL    ,INTENT(IN) :: XLAT
REAL    ,INTENT(IN) :: ZERHT
REAL    ,INTENT(IN) :: REARTH

!       Arrays arguments with intent in:
REAL    ,INTENT(IN) :: PUS(MAXLEV)
REAL    ,INTENT(IN) :: DUS(MAXLEV)
REAL    ,INTENT(IN) :: GASUS(MAXLEV)

!       Array arguments with intent out:
REAL    ,INTENT(OUT):: HUS(MAXLEV)

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,IR
REAL    :: RTOD
REAL    :: AVOG
REAL    :: AIRMOL
REAL    :: H2OMOL
REAL    :: GRED
REAL    :: GRAV0
REAL    :: ARG
REAL    :: A1
REAL    :: ZZ
REAL    :: H0
LOGICAL :: SWITCH

!       Local arrays:
REAL    :: ARX(MAXLEV)
REAL    :: ARY(MAXLEV)
REAL    :: DUM(MAXLEV)


!-----End of header-----------------------------------------------------


DATA RTOD   /0.01745326/
DATA AVOG   /6.022E26/
DATA AIRMOL /28.964/
DATA H2OMOL /18.0115/

DUM(:)=1.

!-----VALUE OF GRAVITY AT EARTHS SURFACE AT LATITUDE XLAT [m.s-2]-------

ARG = COS(2.0*RTOD*XLAT)

!-----Planet specific:- calculation of GRAV0----------------------------

GRAV0 = 9.80616 - 0.025928*ARG + 0.000069*ARG*ARG

!-----Assume gravity falls off as g(h) = g(0) - 2*g(0)*h/R(earth)-------
!                                      = g(0) - GRED*h
GRED = 2.0E-3*GRAV0/REARTH
!-----------------------------------------------------------------------

DO IR=1,NR
  J = NR - IR + 1
  ARX(J) = PUS(IR)
  ARY(J)  = DUS(IR)*(AIRMOL* &
  (1.0E6 - GASUS(IR))+ H2OMOL*GASUS(IR))/AVOG
ENDDO

HUS(1) = ZERHT                     ! Height of the lowest pressure level

DO IR = 2,NR
  SWITCH = .FALSE.
  CALL RTIASI_INTEG &
  (MAXLEV,NR,NP,PUS(IR),PUS(IR-1), &
  ARX,ARY,4,SWITCH,DUM,2,DUM,2,ZZ)

  ZZ = ZZ*1.0E2                    ! Integral in units [m2.s-1]
  A1 = GRAV0/GRED
  H0 = 1.0E3*HUS(IR-1)
  HUS(IR) = 1.0E-3*(A1 -           & ! Height of this pressure level
  SQRT(A1*A1 - H0*(2.0*A1 - H0) - 2.0*ZZ/GRED))
ENDDO

RETURN
END SUBROUTINE RTIASI_PRESHT
!     Calculate the heights of pressure levels
SUBROUTINE RTIASI_PRESHTK(NR,NP,XLAT,ZERHT,REARTH,PUS_D,PUS,DUS_D,DUS, &
GASUS_D,GASUS,HUS_D,HUS,MAXLEV,KNPF,KNCHPF,KPROF)

!     Description:
!     K of subroutine PRESHT
!     To calculate the heights of pressure levels
!
!     Method.
!     A calculation is performed based on the hydrostatic
!     equation to obtain the heights of pressure levels.
!     virtual temperatures are used to account for the
!     presence of water vapour. the variation of gravity
!     with latitude and height is also considered.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            12/11/1999  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


 IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN) :: KNPF
INTEGER,INTENT(IN) :: KNCHPF
INTEGER,INTENT(IN) :: NR
INTEGER,INTENT(IN) :: NP
INTEGER,INTENT(IN) :: MAXLEV
REAL,INTENT(IN)    :: XLAT
REAL,INTENT(IN)    :: ZERHT
REAL,INTENT(IN)    :: REARTH

!       Arrays arguments with intent in:
INTEGER,INTENT(IN) :: KPROF   (KNCHPF)
REAL,INTENT(IN)    :: PUS_D   (MAXLEV,KNPF)
REAL,INTENT(IN)    :: DUS_D   (MAXLEV,KNPF)
REAL,INTENT(IN)    :: GASUS_D (MAXLEV,KNPF)

!      K array arguments with intent out:
REAL,INTENT(INOUT) :: PUS     (MAXLEV,KNCHPF)
REAL,INTENT(INOUT) :: DUS     (MAXLEV,KNCHPF)
REAL,INTENT(INOUT) :: GASUS   (MAXLEV,KNCHPF)

!       Array arguments with intent out:
REAL,INTENT(OUT)   :: HUS_D   (MAXLEV,KNPF)

!       K array arguments with intent in:
REAL,INTENT(INOUT) :: HUS     (MAXLEV,KNCHPF)

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,IR,JJ,I
REAL    :: RTOD
REAL    :: AVOG
REAL    :: AIRMOL
REAL    :: H2OMOL
REAL    :: GRED
REAL    :: GRAV0
REAL    :: ARG
REAL    :: A1
REAL    :: ZZI_D(KNPF)
REAL    :: HINC_D
REAL    :: A_D
REAL    :: B_D
REAL    :: SA_D
REAL    :: SY_D
REAL    :: SY2_D
REAL    :: SB_D
REAL    :: S4_D
REAL    :: SY4_D
REAL    :: SYA_D
REAL    :: SYB_D
REAL    :: H_D
REAL    :: SS_D
REAL    :: S2_D
REAL    :: XINT_D
LOGICAL :: SWITCH

!       K local scalars:
REAL    :: ZZ
REAL    :: H0
REAL    :: H
REAL    :: ZZDUM
REAL    :: HINC
REAL    :: SA
REAL    :: SB
REAL    :: S4
REAL    :: S2
REAL    :: SS
REAL    :: SY

!       Local arrays:
REAL    :: ARX   (MAXLEV,KNPF)
REAL    :: ARY_D (MAXLEV,KNPF)
REAL    :: DUM   (MAXLEV)
REAL    :: H0_D  (MAXLEV,KNPF)
REAL    :: ZZ_D  (MAXLEV,KNPF)

!       K local arrays:
REAL    :: ARY   (MAXLEV,KNCHPF)


!-----End of header-----------------------------------------------------



DATA RTOD    /0.01745326/
DATA AVOG    /6.022E26/
DATA AIRMOL  /28.964/
DATA H2OMOL  /18.0115/

SWITCH=.FALSE.

DUM(1:MAXLEV)=1.


!-----VALUE OF GRAVITY AT EARTHS SURFACE AT LATITUDE XLAT [m.s-2]-------

ARG = COS(2.0*RTOD*XLAT)

!-----Planet specific:- calculation of GRAV0----------------------------

GRAV0 = 9.80616 - 0.025928*ARG + 0.000069*ARG*ARG

!-----Assume gravity falls off as g(h) = g(0) - 2*g(0)*h/R(earth)-------
!                                      = g(0) - GRED*h
GRED = 2.0E-3*GRAV0/REARTH
!-----------------------------------------------------------------------

A1      = GRAV0/GRED



!-----REPEAT DIRECT CALCULATIONS----------------------------------------

DO JJ=1,KNPF
  HUS_D(1,JJ)=ZERHT

  DO IR=1,NR
    J = NR - IR + 1
    ARX(J,JJ) = PUS_D(IR,JJ)
    ARY_D(J,JJ)  = DUS_D(IR,JJ)*(AIRMOL* &
    (1.0E6 - GASUS_D(IR,JJ))+ H2OMOL*GASUS_D(IR,JJ))/AVOG
  ENDDO

    DO IR = 2,NR

      CALL RTIASI_INTEGTL &
      (MAXLEV,NR,NP,PUS_D(IR,JJ),PUS_D(IR,JJ),PUS_D(IR-1,JJ), &
      PUS_D(IR-1,JJ),ARX(:,JJ),ARY_D(:,JJ),DUM,4,SWITCH,DUM, &
      2,DUM,2,ZZI_D(JJ),ZZDUM)


      ZZ_D(IR,JJ)= ZZI_D(JJ)
      ZZ_D(IR,JJ)= ZZ_D(IR,JJ)*1.0E2

      H0_D(IR,JJ)    = 1.0E3*HUS_D(IR-1,JJ)

      HUS_D(IR,JJ) = 1.0E-3*(A1 -          & ! Height of this pressure level
      SQRT(A1*A1 - H0_D(IR,JJ)*(2.0*A1 - H0_D(IR,JJ)) - &
      2.0*ZZ_D(IR,JJ)/GRED))
    ENDDO
ENDDO


!-----------------------------------------------------------------------




!-----DO K CALCULATIONS-------------------------------------------------

ARY(:,:)=0

DO J=1,KNCHPF

  ZZ=0.
  H0=0.
  ZZDUM=0.

  DO IR = NR,2,-1


    H0=H0+1.0E-3*(-0.5*(-2*A1*HUS(IR,J)+2*H0_D(IR,KPROF(J))* &
    HUS(IR,J))/SQRT(A1*A1 - H0_D(IR,KPROF(J))* &
    (2.0*A1 - H0_D(IR,KPROF(J))) - 2.0*ZZ_D(IR,KPROF(J))/GRED))

    ZZ=ZZ+1.0E-3*(-0.5*(-2*HUS(IR,J)/GRED)/ &
    SQRT(A1*A1 - H0_D(IR,KPROF(J))*(2.0*A1 - H0_D(IR,KPROF(J))) - &
    2.0*ZZ_D(IR,KPROF(J))/GRED))

    HUS(IR,J)=0.


    HUS(IR-1,J)=HUS(IR-1,J)+1.0E3*H0
    H0=0.

    ZZ=ZZ*1.0E2

!----------------------------------------------------------------------------

    HINC_D = (PUS_D(IR-1,KPROF(J))-PUS_D(IR,KPROF(J)))/FLOAT(2*NP)
    A_D=PUS_D(IR,KPROF(J))
    B_D=PUS_D(IR-1,KPROF(J))


    IF (.NOT. SWITCH) THEN
      SY_D=ARY_D(NR-IR+1,KPROF(J))
      SY_D=1/SY_D
      SA_D = SY_D
    ENDIF

    IF (.NOT. SWITCH) THEN
      SY_D=ARY_D(NR-IR+2,KPROF(J))
      SY_D=1/SY_D
      SB_D = SY_D
    ENDIF



    S4_D = 0.

    DO I=1,NP
      ARG = A_D + HINC_D*FLOAT(2*I - 1)
      IF (.NOT. SWITCH) THEN
        H_D=-(ARX(NR-IR+1,KPROF(J)) - ARX(NR-IR+2,KPROF(J)))/ &
        LOG(ARY_D(NR-IR+1,KPROF(J))/ARY_D(NR-IR+2,KPROF(J)))
        SY_D= ARY_D(NR-IR+2,KPROF(J))* &
        EXP(-(ARG-ARX(NR-IR+2,KPROF(J)))/H_D)
        SY_D=1/SY_D
        SS_D = SY_D
      ENDIF
      S4_D=S4_D+SS_D
    ENDDO


    S2_D=0
    DO I=1,NP-1
      ARG = A_D + HINC_D*FLOAT(2*I)
        IF (.NOT. SWITCH) THEN
         H_D=-(ARX(NR-IR+1,KPROF(J)) - ARX(NR-IR+2,KPROF(J)))/ &
         LOG(ARY_D(NR-IR+1,KPROF(J))/ARY_D(NR-IR+2,KPROF(J)))
         SY_D= ARY_D(NR-IR+2,KPROF(J))* &
         EXP(-(ARG-ARX(NR-IR+2,KPROF(J)))/H_D)
         SY_D=1/SY_D
         SS_D=SY_D
       ENDIF
      S2_D=S2_D+SS_D
    ENDDO


    XINT_D = HINC_D/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D)



    HINC=0.
    SA=0.
    SB=0.
    S4=0.
    S2=0.
    SS=0.
    SY=0.
    H=0.

    HINC=HINC+ZZ/3.0 * (SA_D + SB_D + 4.0*S4_D + 2.0*S2_D)
    SA=SA+HINC_D/3.0 *ZZ
    SB=SB+HINC_D/3.0 *ZZ
    S4=S4+HINC_D/3.0 * (4.0*ZZ)
    S2=S2+ HINC_D/3.0 * (2.0*ZZ)

    DO I=NP-1,1,-1
      ARG = A_D + HINC_D*FLOAT(2*I)     ! Argument for function values at
                                    ! odd integration points.
      SS=SS+S2
      IF (.NOT. SWITCH) THEN
        SY=SY+SS
        SS=0.

        H_D=-(ARX(NR-IR+1,KPROF(J)) - ARX(NR-IR+2,KPROF(J)))/ &
        LOG(ARY_D(NR-IR+1,KPROF(J))/ARY_D(NR-IR+2,KPROF(J)))
        SY2_D= ARY_D(NR-IR+2,KPROF(J))* &
        EXP(-(ARG-ARX(NR-IR+2,KPROF(J)))/H_D)

        SY= -SY/(SY2_D)**2
        ARY(NR-IR+2,J)=ARY(NR-IR+2,J)+SY
        H=H+SY*(ARG - ARX(NR-IR+2,KPROF(J)))
        ARY(NR-IR+1,J)=ARY(NR-IR+1,J)+H/(ARX(NR-IR+1,KPROF(J)) - &
        ARX(NR-IR+2,KPROF(J)))
        ARY(NR-IR+2,J)=ARY(NR-IR+2,J)-H/(ARX(NR-IR+1,KPROF(J)) - &
        ARX(NR-IR+2,KPROF(J)))
      ENDIF
      H=0
      SY=0
    ENDDO


    DO I=NP,1,-1
      ARG = A_D + HINC_D*FLOAT(2*I - 1) ! Argument for function values at even
                                    ! integration points.
      SS=SS+S4
      IF (.NOT. SWITCH) THEN
        SY=SY+SS
        SS=0.

        H_D=-(ARX(NR-IR+1,KPROF(J)) - ARX(NR-IR+2,KPROF(J)))/ &
        LOG(ARY_D(NR-IR+1,KPROF(J))/ARY_D(NR-IR+2,KPROF(J)))
        SY4_D= ARY_D(NR-IR+2,KPROF(J))* &
        EXP(-(ARG-ARX(NR-IR+2,KPROF(J)))/H_D)

        SY= -SY/(SY4_D)**2
        ARY(NR-IR+2,J)=ARY(NR-IR+2,J)+SY
        H=H+SY*(ARG - ARX(NR-IR+2,KPROF(J)))
        ARY(NR-IR+1,J)=ARY(NR-IR+1,J)+H/(ARX(NR-IR+1,KPROF(J)) - &
        ARX(NR-IR+2,KPROF(J)))
        ARY(NR-IR+2,J)=ARY(NR-IR+2,J)-H/(ARX(NR-IR+1,KPROF(J)) - &
        ARX(NR-IR+2,KPROF(J)))

        H=0.
        SY=0.

      ENDIF
    ENDDO


    IF (.NOT. SWITCH) THEN
      SY=SY+SB
      ARG=B_D
      SYB_D=ARY_D(NR-IR+2,KPROF(J))
      SY= -SY/(SYB_D)**2
      ARY(NR-IR+2,J)=ARY(NR-IR+2,J)+SY
      H=0.
      SY=0.
    ENDIF


    IF (.NOT. SWITCH) THEN
      SY=SY+SA
      ARG=A_D
      SYA_D=ARY_D(NR-IR+1,KPROF(J))
      SY= -SY/(SYA_D)**2
      ARY(NR-IR+1,J)=ARY(NR-IR+1,J)+SY
      H=0.
      SY=0.
    ENDIF


    PUS(IR-1,J)  =PUS(IR-1,J)+HINC/FLOAT(2*NP)
    PUS(IR,J)  =PUS(IR,J)-HINC/FLOAT(2*NP)

!----------------------------------------------------------------------

    ZZ=0.

  ENDDO

  HUS(1,J)  = 0.

  DO IR=NR,1,-1
    JJ = NR - IR + 1
    DUS(IR,J)=DUS(IR,J)+ARY(JJ,J)*(AIRMOL* &
    (1.0E6 - GASUS_D(IR,KPROF(J)))+ H2OMOL*GASUS_D(IR,KPROF(J)))/ &
    AVOG
    GASUS(IR,J)=GASUS(IR,J)+DUS_D(IR,KPROF(J))*(AIRMOL* &
    (-ARY(JJ,J))+H2OMOL*ARY(JJ,J))/AVOG
    ARY(JJ,J)=0.
  ENDDO
ENDDO

RETURN

END SUBROUTINE RTIASI_PRESHTK
!     Trace a ray  through a path
SUBROUTINE RTIASI_RAYTCE(NPTS,NP,A,B,ALT,TEMP,PRES,AIRNO,REARTH,SS,TFIN, &
                 MAXLEV)

!     Description:
!     Trace a ray through a path
!
!     Method:
!     A ray is traced through the sublayers of each layer
!     and the curtis-godson integrated quantites are
!     calculated. Simpsons rule is used to integrate over
!     two panels in each sub-layer. The refraction is not
!     acounted for.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            95/12/28    Original Code. 4.0   D.P. EDWARDS
!     2            12/11/99    Marco MAtricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER ,INTENT(IN)    :: MAXLEV
INTEGER ,INTENT(IN)    :: NPTS           ! Number of profile level points.
INTEGER ,INTENT(IN)    :: NP             ! Accuracy criterion for
                                         ! integration.
REAL    ,INTENT(IN)    :: REARTH         ! Earth radius for this profile
                                         ! [km].
REAL    ,INTENT(IN)    :: SS             ! Sine of path initial zenith
                                         ! angle.
REAL    ,INTENT(IN)    :: A              ! Path lower height [km].
REAL    ,INTENT(IN)    :: B              ! Path upper height [km].

!       Array arguments with intent in:
REAL    ,INTENT(IN)    :: ALT   (MAXLEV) ! Data point altitudes [km].
REAL    ,INTENT(IN)    :: TEMP  (MAXLEV) ! Data point temperatures [K].
REAL    ,INTENT(IN)    :: PRES  (MAXLEV) ! Data point pressures [mb].
REAL    ,INTENT(IN)    :: AIRNO (MAXLEV) ! Data point air densities [cm-3]

!       Scalar arguments with intent out:
REAL    ,INTENT(OUT)   :: TFIN           ! Path average temperature [K].

!     End of subroutine arguments


!       Local parameters:
INTEGER,PARAMETER      ::  NBYTES_DP = KIND (1.D0)

!       Local scalars:
INTEGER                :: IS
INTEGER                :: IGAS
INTEGER                :: NSL
REAL                   :: SZ1,SZ2,SZ3    ! Altitude of sub-layer boundaries.
REAL                   :: H1,H2,H3
REAL                   :: YA1,YA2,YA3
REAL                   :: YT1,YT2,YT3
REAL                   :: YP1,YP2,YP3
REAL                   :: ZMID
REAL (KIND=NBYTES_DP)  :: CC
REAL (KIND=NBYTES_DP)  :: FP1,FP2,FP3
REAL (KIND=NBYTES_DP)  :: FT1,FT2,FT3
REAL (KIND=NBYTES_DP)  :: FA1,FA2,FA3
REAL (KIND=NBYTES_DP)  :: FG1,FG2,FG3
REAL (KIND=NBYTES_DP)  :: DS12,DS23
REAL (KIND=NBYTES_DP)  :: DX12,DX23
REAL (KIND=NBYTES_DP)  :: X1,X2,X3
REAL (KIND=NBYTES_DP)  :: CS1,CS2,CS3
REAL (KIND=NBYTES_DP)  :: SN1,SN2,SN3
REAL (KIND=NBYTES_DP)  :: RR1,RR2,RR3
REAL (KIND=NBYTES_DP)  :: Z1,Z2,Z3
REAL (KIND=NBYTES_DP)  :: SINCR          ! Nominal path lenght.
REAL (KIND=NBYTES_DP)  :: DZMIN
REAL (KIND=NBYTES_DP)  :: DELZ
REAL (KIND=NBYTES_DP)  :: SRAY           ! Ray lenght.
REAL (KIND=NBYTES_DP)  :: ZUS
REAL (KIND=NBYTES_DP)  :: ZTS
REAL (KIND=NBYTES_DP)  :: ZPS
LOGICAL                :: FINAL


!-----End of header-----------------------------------------------------


IGAS=1

!-----Lower boundary altitude-------------------------------------------

Z1 = A
SZ1 = A


!-----Sine of initial zenith angle--------------------------------------

SN1 = SS


!-----Nominal path length set by user input-----------------------------

FINAL = .FALSE.
SINCR = (B - Z1)/NP
DZMIN = 0.5D0*SINCR**2/(Z1 + REARTH)
!     DZMIN = 0.001D0


!-----Distance from centre of earth-------------------------------------

RR1 = Z1 + REARTH


!-----Calculate snells law constant-------------------------------------

CC = SN1*RR1
NSL = 1



!-----INITIALIZE LOWER BOUNDARY VALUES----------------------------------

!-----Cos of initial zenith angle---------------------------------------
CS1 = DSQRT(1.D0 - SN1**2)

!-----Change to intermediate variable X----------------------------------
X1 = RR1*CS1

!-----Ray length is zero at this stage----------------------------------
SRAY = 0.D0



!-----INITIALISE LAYER LOWER BOUNDARY WEIGH. QUANTITIES FOR INTEGRATION-


!-----Path gas amount---------------------------------------------------

ZUS = 0.D0
CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ1,ALT,AIRNO,YA1,H1)
  IF (IGAS /= 91)THEN
    FG1 = YA1
  END IF
FA1 = YA1


!-----Path average temperature------------------------------------------

ZTS = 0.D0
CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ1,ALT,TEMP,YT1,H1)
FT1 = FG1*YT1


!-----Path average pressure---------------------------------------------

ZPS = 0.D0
CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ1,ALT,PRES,YP1,H1)
FP1 = FG1*YP1



!-----LOOP OVER SUB-LAYERS----------------------------------------------
!-----CALCULATE EXTENT OF THIS SUB-LAYER--------------------------------
!-----ALLOW FOR COS-->0 AT THE TANGENT HEIGHT---------------------------

DO IS=1,10000
DELZ = SINCR*CS1
  IF (DELZ <= DZMIN)THEN
    DELZ = DZMIN
  END IF
Z3 = Z1 + DELZ

!-----See if top of layer has been reached------------------------------

  IF (Z3 >= B) THEN
    Z3 = B
    DELZ = Z3 - Z1
    FINAL = .TRUE.
  ENDIF
Z2 = Z1 + 0.5D0*DELZ
SZ2 = Z2
SZ3 = Z3
RR2 = REARTH + Z2
RR3 = REARTH + Z3

SN2 = CC/RR2
SN3 = CC/RR3

!-----Need the cosine of the angle--------------------------------------

CS2 = DSQRT(1.D0 - SN2**2)
CS3 = DSQRT(1.D0 - SN3**2)

X2 = RR2*CS2
X3 = RR3*CS3
DX12 = X2 - X1
DX23 = X3 - X2


!-----Ray path length---------------------------------------------------

DS12 = DX12
DS23 = DX23
SRAY = SRAY + DS12 + DS23


!-----Path gas amount---------------------------------------------------

CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2,ALT,AIRNO,YA2,H2)
  IF (IGAS /= 91)THEN
    FG2 = YA2
  END IF
FA2 = YA2
CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3,ALT,AIRNO,YA3,H3)
  IF (IGAS /= 91)THEN
    FG3 = YA3
  END IF
FA3 = YA3
ZUS = ZUS + (DS12*(5.D0*FG1 + 8.D0*FG2 - FG3) &
           + DS23*(-FG1 + 8.D0*FG2 + 5.D0*FG3))/12.D0


!-----Path average temperature------------------------------------------

CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ2,ALT,TEMP,YT2,H2)
FT2 = FG2*YT2
CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ3,ALT,TEMP,YT3,H3)
FT3 = FG3*YT3
ZTS = ZTS + (DS12*(5.D0*FT1 + 8.D0*FT2 - FT3) &
           + DS23*(-FT1 + 8.D0*FT2 + 5.D0*FT3))/12.D0


!-----Path average pressure---------------------------------------------

CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2,ALT,PRES,YP2,H2)
FP2 = FG2*YP2
CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3,ALT,PRES,YP3,H3)
FP3 = FG3*YP3
ZPS = ZPS + (DS12*(5.D0*FP1 + 8.D0*FP2 - FP3) &
           + DS23*(-FP1 + 8.D0*FP2 + 5.D0*FP3))/12.D0


  IF (.NOT. FINAL) THEN

!---------Update quantites at new layer level 1 from old layer level 3--
    NSL = NSL + 1
    Z1 = Z3
    SZ1 = SZ3
    RR1 = RR3
    SN1 = SN3
    CS1 = CS3
    X1 = X3
    FA1=FA3
    FG1 = FG3
    FT1 = FT3
    FP1 = FP3
    CYCLE
  ELSE

!---------FORM FINAL INTEGRATED QUANTITES-------------------------------

      IF (ZUS /= 0.0) THEN
        TFIN = ZTS/ZUS
      ELSE

!-------------If gas amount is zero return averages for pressure and----
!             temperature.

        ZMID = 0.5*(A + B)
        CALL RTIASI_INTER(MAXLEV,NPTS,2,ZMID,ALT,TEMP,TFIN,H2)
      ENDIF
    EXIT
  ENDIF
  END DO

RETURN
END SUBROUTINE RTIASI_RAYTCE
!     Trace a ray through a path
SUBROUTINE RTIASI_RAYTCEK(NPTS,NP,A_D,A,B_D,B,ALT_D,TEMP_D,TEMP, &
                    PRES_D,PRES,AIRNO_D,AIRNO,REARTH,SS,TFIN, &
                    MAXLEV,NOLAY,KNCHPF,KNPF,KPROF)

!     Description:
!     K of subroutine RAYTCE
!
!     Trace a ray  through a path
!
!     Method:
!     A ray is traced through the sublayers of each layer
!     and the curtis-godson integrated quantites are
!     calculated. Simpsons rule is used to integrate over
!     two panels in each sub-layer. The refraction is not
!     acounted for.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            12/11/99    Marco MAtricardi. ECMWF.
!     1.01         04/04/00    Introduced DELTA_ALT_D and constrained to a
!                              minimum value to prevent divide by zeroes.
!                              Andrew Collard. 
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER,INTENT(IN)    :: NOLAY
INTEGER,INTENT(IN)    :: KNCHPF
INTEGER,INTENT(IN)    :: KNPF
INTEGER,INTENT(IN)    :: KPROF   (KNCHPF)
INTEGER,INTENT(IN)    :: MAXLEV
INTEGER,INTENT(IN)    :: NPTS                      ! Number of prof. level
                                                   ! points.
INTEGER,INTENT(IN)    :: NP                        ! Accuracy criterion r
                                                   ! for integration.
REAL,INTENT(IN)       :: REARTH                    ! Earth radius for this
                                                   ! profile [km].
REAL,INTENT(IN)       :: SS      (MAXLEV,KNPF)     ! Sine of path initial
                                                   ! zenith angle.
REAL,INTENT(IN)       :: A_D     (MAXLEV,KNPF)     ! Path lower height[km]
REAL,INTENT(IN)       :: B_D     (MAXLEV,KNPF)     ! Path upper height[km]

!       K scalar arguments with intent out:
REAL,INTENT(INOUT)    :: A       (MAXLEV,KNCHPF)   ! Path lower height[km]
REAL,INTENT(INOUT)    :: B       (MAXLEV,KNCHPF)   ! Path upper height[km]

!       Array arguments with intent in:
REAL,INTENT(IN)       :: ALT_D   (MAXLEV,KNPF)     ! Data point altitudes
                                                   ! [km].
REAL,INTENT(IN)       :: TEMP_D  (MAXLEV,KNPF)     ! Data point temp.
                                                   ! [K].
REAL,INTENT(IN)       :: PRES_D  (MAXLEV,KNPF)     ! Data point pressures
                                                   ! [mb].
REAL,INTENT(IN)       :: AIRNO_D (MAXLEV,KNPF)     ! Data point air
                                                   ! densities [cm-3]

!       K array arguments with intent out:
REAL,INTENT(INOUT)      :: TEMP    (MAXLEV,KNCHPF) ! Data point temp.
                                                   ! [K].
REAL,INTENT(INOUT)      :: PRES    (MAXLEV,KNCHPF) ! Data point pressures
                                                   ! [mb].
REAL,INTENT(INOUT)      :: AIRNO   (MAXLEV,KNCHPF) ! Data point air
                                                   ! densities [cm-3]

!       K scalar arguments with intent in:
REAL,INTENT(INOUT)      :: TFIN    (MAXLEV,KNCHPF) ! Path average
                                                   ! temperature [K].

!     End of subroutine arguments



!       Local parameters:
INTEGER,PARAMETER     :: NBYTES_DP = KIND (1.D0)


!       Local scalars:
INTEGER               :: IS,IPATH,J
INTEGER               :: IGAS
INTEGER               :: NLOOP
REAL                  :: DELTA_ALT_D   ! Height Interval
REAL                  :: SZ1_D   (NP,MAXLEV,KNPF)
REAL                  :: SZ2_D   (NP,MAXLEV,KNPF)
REAL                  :: SZ3_D   (NP,MAXLEV,KNPF)
REAL                  :: H1_D    (MAXLEV,KNPF)
REAL                  :: H2_D    (MAXLEV,KNPF)
REAL                  :: H3_D    (MAXLEV,KNPF)
REAL                  :: YA1_D   (NP,MAXLEV,KNPF)
REAL                  :: YA2_D   (NP,MAXLEV,KNPF)
REAL                  :: YA3_D   (NP,MAXLEV,KNPF)
REAL                  :: YT1_D   (NP,MAXLEV,KNPF)
REAL                  :: YT2_D   (NP,MAXLEV,KNPF)
REAL                  :: YT3_D   (NP,MAXLEV,KNPF)
REAL                  :: YP1_D   (NP,MAXLEV,KNPF)
REAL                  :: YP2_D   (NP,MAXLEV,KNPF)
REAL                  :: YP3_D   (NP,MAXLEV,KNPF)
REAL                  :: ZMID
REAL (KIND=NBYTES_DP) :: CC_D    (MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FP1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FP2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FP3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FT1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FT2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FT3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FA1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FA2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FA3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FG1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FG2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: FG3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DS12_D  (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DS23_D  (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DX12_D  (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DX23_D  (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: X1_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: X2_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: X3_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: CS1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: CS2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: CS3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: SN1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: SN2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: SN3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: RR1_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: RR2_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: RR3_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: Z1_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: Z2_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: Z3_D    (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: SINCR_D (MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DZMIN_D (MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: DELZ_D  (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: ZUS_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: ZTS_D   (NP,MAXLEV,KNPF)
REAL (KIND=NBYTES_DP) :: ZPS_D   (NP,MAXLEV,KNPF)
LOGICAL               :: FINAL

!       K local scalars:
REAL (KIND=NBYTES_DP) :: CC
REAL (KIND=NBYTES_DP) :: SINCR
REAL (KIND=NBYTES_DP) :: DZMIN
REAL (KIND=NBYTES_DP) :: ZUS
REAL                  :: YA1,YA2,YA3
REAL                  :: YT1,YT2,YT3
REAL                  :: YP1,YP2,YP3
REAL                  :: SZ1,SZ2,SZ3
REAL                  :: H1,H2,H3
REAL (KIND=NBYTES_DP) :: DELZ
REAL (KIND=NBYTES_DP) :: ZTS
REAL (KIND=NBYTES_DP) :: ZPS
REAL (KIND=NBYTES_DP) :: Z1,Z2,Z3
REAL (KIND=NBYTES_DP) :: RR1,RR2,RR3
REAL (KIND=NBYTES_DP) :: SN1,SN2,SN3
REAL (KIND=NBYTES_DP) :: CS1,CS2,CS3
REAL (KIND=NBYTES_DP) :: FP1,FP2,FP3
REAL (KIND=NBYTES_DP) :: FT1,FT2,FT3
REAL (KIND=NBYTES_DP) :: FA1,FA2,FA3
REAL (KIND=NBYTES_DP) :: FG1,FG2,FG3
REAL (KIND=NBYTES_DP) :: DS12,DS23
REAL (KIND=NBYTES_DP) :: DX12,DX23
REAL (KIND=NBYTES_DP) :: X1,X2,X3


!-----End of header-----------------------------------------------------


IGAS=1


!-----REPEAT DIRECT COMPUTATIONS----------------------------------------


!-----LOOP OVER SUB-LAYERS----------------------------------------------
!-----CALCULATE EXTENT OF THIS SUB-LAYER--------------------------------
!-----ALLOW FOR COS-->0 AT THE TANGENT HEIGHT---------------------------


DO J=1,KNPF

  DO IPATH=1,NOLAY

    DO IS=1,10000


      IF (IS == 1) THEN

        Z1_D   (IS,IPATH,J) = A_D(IPATH,J)
        SZ1_D  (IS,IPATH,J) = A_D(IPATH,J)
        SN1_D  (IS,IPATH,J) = SS(IPATH,J)
        FINAL               = .FALSE.
        SINCR_D(IPATH,J)    = (B_D(IPATH,J)-Z1_D(IS,IPATH,J))/NP
        DZMIN_D(IPATH,J)    = 0.5D0*SINCR_D(IPATH,J)**2/ &
                              (Z1_D(IS,IPATH,J) +REARTH)
        RR1_D  (IS,IPATH,J) = Z1_D(IS,IPATH,J)+REARTH
        CC_D   (IPATH,J)    = SN1_D(IS,IPATH,J)*(Z1_D(IS,IPATH,J)+ &
                              REARTH)
        CS1_D  (IS,IPATH,J) = DSQRT(1.D0 - SN1_D(IS,IPATH,J)**2)
        X1_D   (IS,IPATH,J) = RR1_D(IS,IPATH,J)*CS1_D(IS,IPATH,J)
        ZUS_D  (IS,IPATH,J) = 0.D0


        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ1_D(IS,IPATH,J),ALT_D(:,J), &
        AIRNO_D(:,J),YA1_D(IS,IPATH,J),H1_D(IPATH,J))

        IF (IGAS /= 91)THEN
          FG1_D(IS,IPATH,J) = YA1_D(IS,IPATH,J)
        END IF
        FA1_D(IS,IPATH,J) = YA1_D(IS,IPATH,J)


        ZTS_D(IS,IPATH,J) = 0.D0
        CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ1_D(IS,IPATH,J),ALT_D(:,J), &
               TEMP_D(:,J),YT1_D(IS,IPATH,J),H1_D(IPATH,J))
        FT1_D(IS,IPATH,J) = FG1_D(IS,IPATH,J)*YT1_D(IS,IPATH,J)


        ZPS_D(IS,IPATH,J) = 0.D0
        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ1_D(IS,IPATH,J),ALT_D(:,J), &
               PRES_D(:,J),YP1_D(IS,IPATH,J),H1_D(IPATH,J))
        FP1_D(IS,IPATH,J) = FG1_D(IS,IPATH,J)*YP1_D(IS,IPATH,J)


        DELZ_D(IS,IPATH,J)=SINCR_D(IPATH,J)*CS1_D(IS,IPATH,J)

        Z3_D(IS,IPATH,J) = Z1_D(IS,IPATH,J) + DELZ_D(IS,IPATH,J)

        IF (IS == 10) THEN
          Z3_D(IS,IPATH,J)=B_D(IPATH,J)
          DELZ_D(IS,IPATH,J)=Z3_D(IS,IPATH,J) - Z1_D(IS,IPATH,J)
          FINAL = .TRUE.
          NLOOP =IS
        ENDIF

       Z2_D(IS,IPATH,J)=Z1_D(IS,IPATH,J)+0.5D0*DELZ_D(IS,IPATH,J)

        SZ2_D(IS,IPATH,J) = Z2_D(IS,IPATH,J)

        SZ3_D(IS,IPATH,J)= Z3_D(IS,IPATH,J)

        RR2_D(IS,IPATH,J)=REARTH + Z2_D(IS,IPATH,J)
        RR3_D(IS,IPATH,J)=REARTH + Z3_D(IS,IPATH,J)


        SN2_D(IS,IPATH,J)=CC_D(IPATH,J)/RR2_D(IS,IPATH,J)
        SN3_D(IS,IPATH,J)=CC_D(IPATH,J)/RR3_D(IS,IPATH,J)


!-----Need the cosine of the angle--------------------------------------

        CS2_D(IS,IPATH,J) = DSQRT(1.D0 - SN2_D(IS,IPATH,J)**2)
        CS3_D(IS,IPATH,J) = DSQRT(1.D0 - SN3_D(IS,IPATH,J)**2)


        X2_D(IS,IPATH,J)=RR2_D(IS,IPATH,J)*CS2_D(IS,IPATH,J)
        X3_D(IS,IPATH,J)=RR3_D(IS,IPATH,J)*CS3_D(IS,IPATH,J)


        DX12_D(IS,IPATH,J)=X2_D(IS,IPATH,J)-X1_D(IS,IPATH,J)
        DX23_D(IS,IPATH,J)=X3_D(IS,IPATH,J)-X2_D(IS,IPATH,J)



!-----Ray path length---------------------------------------------------

        DS12_D(IS,IPATH,J)=DX12_D(IS,IPATH,J)
        DS23_D(IS,IPATH,J)=DX23_D(IS,IPATH,J)



!-----Path gas amount---------------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
               AIRNO_D(:,J),YA2_D(IS,IPATH,J),H2_D(IPATH,J))
        IF (IGAS /= 91)THEN
          FG2_D(IS,IPATH,J) = YA2_D(IS,IPATH,J)
        END IF
        FA2_D(IS,IPATH,J) = YA2_D(IS,IPATH,J)


        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               AIRNO_D(:,J),YA3_D(IS,IPATH,J),H3_D(IPATH,J))
        IF (IGAS /= 91)THEN
          FG3_D(IS,IPATH,J) = YA3_D(IS,IPATH,J)
        END IF
        FA3_D(IS,IPATH,J) = YA3_D(IS,IPATH,J)


        ZUS_D(IS,IPATH,J)= ZUS_D(IS,IPATH,J)+ (DS12_D(IS,IPATH,J)* &
        (5.D0*FG1_D(IS,IPATH,J) + 8.D0*FG2_D(IS,IPATH,J) - &
        FG3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FG1_D(IS,IPATH,J)+ &
        8.D0*FG2_D(IS,IPATH,J) +5.D0*FG3_D(IS,IPATH,J)))/12.D0



!-----Path average temperature------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
                TEMP_D(:,J),YT2_D(IS,IPATH,J),H2_D(IPATH,J))
        FT2_D(IS,IPATH,J) = FG2_D(IS,IPATH,J)*YT2_D(IS,IPATH,J)


        CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               TEMP_D(:,J),YT3_D(IS,IPATH,J),H3_D(IPATH,J))
        FT3_D(IS,IPATH,J) = FG3_D(IS,IPATH,J)*YT3_D(IS,IPATH,J)


        ZTS_D(IS,IPATH,J)= ZTS_D(IS,IPATH,J)+ (DS12_D(IS,IPATH,J)* &
        (5.D0*FT1_D(IS,IPATH,J) + 8.D0*FT2_D(IS,IPATH,J) - &
        FT3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FT1_D(IS,IPATH,J)+ &
        8.D0*FT2_D(IS,IPATH,J) +5.D0*FT3_D(IS,IPATH,J)))/12.D0



!-----Path average pressure---------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
                PRES_D(:,J),YP2_D(IS,IPATH,J),H2_D(IPATH,J))
        FP2_D(IS,IPATH,J)=FG2_D(IS,IPATH,J)*YP2_D(IS,IPATH,J)

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               PRES_D(:,J),YP3_D(IS,IPATH,J),H3_D(IPATH,J))
        FP3_D(IS,IPATH,J)=FG3_D(IS,IPATH,J)*YP3_D(IS,IPATH,J)

        ZPS_D(IS,IPATH,J)= ZPS_D(IS,IPATH,J)+ (DS12_D(IS,IPATH,J)* &
        (5.D0*FP1_D(IS,IPATH,J) + 8.D0*FP2_D(IS,IPATH,J) - &
        FP3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FP1_D(IS,IPATH,J)+ &
        8.D0*FP2_D(IS,IPATH,J) +5.D0*FP3_D(IS,IPATH,J)))/12.D0

      ELSE IF (IS /= 1 ) THEN

        Z1_D  (IS,IPATH,J) = Z3_D  (IS-1,IPATH,J)
        SZ1_D (IS,IPATH,J) = SZ3_D (IS-1,IPATH,J)
        RR1_D (IS,IPATH,J) = RR3_D (IS-1,IPATH,J)
        SN1_D (IS,IPATH,J) = SN3_D (IS-1,IPATH,J)
        CS1_D (IS,IPATH,J) = CS3_D (IS-1,IPATH,J)
        X1_D  (IS,IPATH,J) = X3_D  (IS-1,IPATH,J)
        FA1_D (IS,IPATH,J) = FA3_D (IS-1,IPATH,J)
        FG1_D (IS,IPATH,J) = FG3_D (IS-1,IPATH,J)
        FT1_D (IS,IPATH,J) = FT3_D (IS-1,IPATH,J)
        FP1_D (IS,IPATH,J) = FP3_D (IS-1,IPATH,J)

       DELZ_D(IS,IPATH,J)=SINCR_D(IPATH,J)*CS1_D(IS,IPATH,J)

        Z3_D(IS,IPATH,J) = Z1_D(IS,IPATH,J) + DELZ_D(IS,IPATH,J)

        IF (IS == 10) THEN
          Z3_D(IS,IPATH,J)=B_D(IPATH,J)
          DELZ_D(IS,IPATH,J)=Z3_D(IS,IPATH,J) - Z1_D(IS,IPATH,J)
          FINAL = .TRUE.
          NLOOP =IS
        ENDIF

       Z2_D(IS,IPATH,J)=Z1_D(IS,IPATH,J)+0.5D0*DELZ_D(IS,IPATH,J)

        SZ2_D(IS,IPATH,J) = Z2_D(IS,IPATH,J)

        SZ3_D(IS,IPATH,J)= Z3_D(IS,IPATH,J)

        RR2_D(IS,IPATH,J)=REARTH + Z2_D(IS,IPATH,J)
        RR3_D(IS,IPATH,J)=REARTH + Z3_D(IS,IPATH,J)


        SN2_D(IS,IPATH,J)=CC_D(IPATH,J)/RR2_D(IS,IPATH,J)
        SN3_D(IS,IPATH,J)=CC_D(IPATH,J)/RR3_D(IS,IPATH,J)


!-----Need the cosine of the angle--------------------------------------

        CS2_D(IS,IPATH,J) = DSQRT(1.D0 - SN2_D(IS,IPATH,J)**2)
        CS3_D(IS,IPATH,J) = DSQRT(1.D0 - SN3_D(IS,IPATH,J)**2)


        X2_D(IS,IPATH,J)=RR2_D(IS,IPATH,J)*CS2_D(IS,IPATH,J)
        X3_D(IS,IPATH,J)=RR3_D(IS,IPATH,J)*CS3_D(IS,IPATH,J)


        DX12_D(IS,IPATH,J)=X2_D(IS,IPATH,J)-X1_D(IS,IPATH,J)
        DX23_D(IS,IPATH,J)=X3_D(IS,IPATH,J)-X2_D(IS,IPATH,J)



!-----Ray path length---------------------------------------------------

        DS12_D(IS,IPATH,J)=DX12_D(IS,IPATH,J)
        DS23_D(IS,IPATH,J)=DX23_D(IS,IPATH,J)



!-----Path gas amount---------------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
               AIRNO_D(:,J),YA2_D(IS,IPATH,J),H2_D(IPATH,J))
        IF (IGAS /= 91)THEN
          FG2_D(IS,IPATH,J) = YA2_D(IS,IPATH,J)
        END IF
        FA2_D(IS,IPATH,J) = YA2_D(IS,IPATH,J)


        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               AIRNO_D(:,J),YA3_D(IS,IPATH,J),H3_D(IPATH,J))
        IF (IGAS /= 91)THEN
          FG3_D(IS,IPATH,J) = YA3_D(IS,IPATH,J)
        END IF
        FA3_D(IS,IPATH,J) = YA3_D(IS,IPATH,J)


        ZUS_D(IS,IPATH,J)=ZUS_D(IS-1,IPATH,J)+(DS12_D(IS,IPATH,J)* &
        (5.D0*FG1_D(IS,IPATH,J) + 8.D0*FG2_D(IS,IPATH,J) - &
        FG3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FG1_D(IS,IPATH,J)+ &
        8.D0*FG2_D(IS,IPATH,J) +5.D0*FG3_D(IS,IPATH,J)))/12.D0





!-----Path average Temperature------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
               TEMP_D(:,J),YT2_D(IS,IPATH,J),H2_D(IPATH,J))
        FT2_D(IS,IPATH,J) = FG2_D(IS,IPATH,J)*YT2_D(IS,IPATH,J)


        CALL RTIASI_INTER(MAXLEV,NPTS,2,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               TEMP_D(:,J),YT3_D(IS,IPATH,J),H3_D(IPATH,J))
        FT3_D(IS,IPATH,J) = FG3_D(IS,IPATH,J)*YT3_D(IS,IPATH,J)


        ZTS_D(IS,IPATH,J)=ZTS_D(IS-1,IPATH,J)+(DS12_D(IS,IPATH,J)* &
        (5.D0*FT1_D(IS,IPATH,J) + 8.D0*FT2_D(IS,IPATH,J) - &
        FT3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FT1_D(IS,IPATH,J)+ &
        8.D0*FT2_D(IS,IPATH,J) +5.D0*FT3_D(IS,IPATH,J)))/12.D0






!-----Path average pressure---------------------------------------------

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ2_D(IS,IPATH,J),ALT_D(:,J), &
               PRES_D(:,J),YP2_D(IS,IPATH,J),H2_D(IPATH,J))
        FP2_D(IS,IPATH,J)=FG2_D(IS,IPATH,J)*YP2_D(IS,IPATH,J)

        CALL RTIASI_INTER(MAXLEV,NPTS,1,SZ3_D(IS,IPATH,J),ALT_D(:,J), &
               PRES_D(:,J),YP3_D(IS,IPATH,J),H3_D(IPATH,J))
        FP3_D(IS,IPATH,J)=FG3_D(IS,IPATH,J)*YP3_D(IS,IPATH,J)

        ZPS_D(IS,IPATH,J)=ZPS_D(IS-1,IPATH,J)+(DS12_D(IS,IPATH,J)* &
        (5.D0*FP1_D(IS,IPATH,J) + 8.D0*FP2_D(IS,IPATH,J) - &
        FP3_D(IS,IPATH,J))+DS23_D(IS,IPATH,J)*(-FP1_D(IS,IPATH,J)+ &
        8.D0*FP2_D(IS,IPATH,J) +5.D0*FP3_D(IS,IPATH,J)))/12.D0

       IF(FINAL)THEN
          EXIT
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO


!-----END OF DIRECT COMPUTATIONS----------------------------------------





!-----DO K COMPUTATIONS-------------------------------------------------



DO J=1,KNCHPF
  DO IPATH=NOLAY,1,-1


!-----Initialize local K variables--------------------------------------

  FP3=0.
  FT3=0.
  FG3=0.
  FA3=0.
  FA1=0.
  X3=0.
  CS3=0.
  RR3=0.
  SZ3=0.
  SZ2=0.
  Z3=0.
  ZTS=0.
  ZUS=0.
  DS12=0.
  FP1=0.
  FP2=0.
  DS23=0.
  YP3=0.
  FG3=0.
  YP2=0.
  FG2=0.
  FT1=0.
  FT2=0.
  YT3=0.
  FG3=0.
  YT2=0.
  FG2=0.
  FG1=0.
  FG2=0.
  FG3=0.
  YA3=0.
  YA2=0.
  DX23=0.
  DX12=0.
  X2=0.
  RR2=0.
  CS2=0.
  SN3=0.
  SN2=0.
  CC=0.
  RR2=0.
  Z2=0.
  Z1=0.
  DELZ=0.
  YP1=0.
  FG1=0.
  YT1=0.
  FG1=0.
  YA1=0.
  CS1=0.
  RR1=0.
  SN1=0.
  SINCR=0.
  X1=0.
  SZ1=0.
  ZPS=0.
  H1=0.
  H2=0.
  H3=0.
  FA2=0.
  DZMIN=0.
!-----------------------------------------------------------------------


  IF (ZUS_D(NLOOP,IPATH,KPROF(J)) /= 0.0) THEN
    ZTS=ZTS+TFIN(IPATH,J)/ZUS_D(NLOOP,IPATH,KPROF(J))
    ZUS=ZUS-ZTS_D(NLOOP,IPATH,KPROF(J))*TFIN(IPATH,J)/ &
    ZUS_D(NLOOP,IPATH,KPROF(J))**2
  ELSE

!-----If gas amount is zero return averages for pressure and------------
!     temperature.

    ZMID = 0.5*(A_D(IPATH,KPROF(J)) + B_D(IPATH,KPROF(J)))
    CALL RTIASI_INTER(MAXLEV,NPTS,2,ZMID,ALT_D(:,KPROF(J)),TEMP(IPATH,J), &
    TFIN(IPATH,J),H2)
  ENDIF


!-----LOOP OVER SUB-LAYERS----------------------------------------------
!-----CALCULATE EXTENT OF THIS SUB-LAYER--------------------------------
!-----ALLOW FOR COS-->0 AT THE TANGENT HEIGHT---------------------------

    DO IS=NLOOP,1,-1

      IF (IS <= NLOOP-1) THEN
        FP3 =FP3+FP1
        FP1 =0.
        FT3 =FT3+FT1
        FT1 =0.
        FG3 =FG3+FG1
        FG1 =0.
        FA3 =FA3+FA1
        FA1 =0.
        X3  =X3+X1
        X1  = 0.
        CS3 =CS3+CS1
        CS1 =0.
        SN3 =SN3+SN1
        SN1 =0.
        RR3 =RR3+RR1
        RR1 =0.
        SZ3 =SZ3+SZ1
        SZ1 =0.
        Z3  =Z3+Z1
        Z1  =0.
      ENDIF


!-----Path average pressure---------------------------------------------

      DELTA_ALT_D = &    ! The smallest allowed layer depth is 1cm
           MAX(ALT_D(IPATH+1,KPROF(J))-ALT_D(IPATH,KPROF(J)),1.e-5)

      DS12=DS12+ZPS*(5.D0*FP1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FP2_D(IS,IPATH,KPROF(J)) - FP3_D(IS,IPATH,KPROF(J)))/ &
      12.D0
      FP1=FP1+DS12_D(IS,IPATH,KPROF(J))*5.D0*ZPS/12.D0
      FP2=FP2+DS12_D(IS,IPATH,KPROF(J))*8.D0*ZPS/12.D0
      FP3=FP3-DS12_D(IS,IPATH,KPROF(J))*ZPS/12.D0
      DS23=DS23+ZPS*(-FP1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FP2_D(IS,IPATH,KPROF(J))+ &
      5.D0*FP3_D(IS,IPATH,KPROF(J)))/12.D0
      FP1=FP1- DS23_D(IS,IPATH,KPROF(J))*ZPS/12.D0
      FP2=FP2+DS23_D(IS,IPATH,KPROF(J))*8.D0*ZPS/12.D0
      FP3=FP3+DS23_D(IS,IPATH,KPROF(J))*5.D0*ZPS/12.D0

      YP3=YP3+FP3*FG3_D(IS,IPATH,KPROF(J))
      FG3=FG3+FP3*YP3_D(IS,IPATH,KPROF(J))
      FP3=0.

      PRES(IPATH,J)=PRES(IPATH,J)+YP3
      PRES(IPATH+1,J)=PRES(IPATH+1,J)+YP3* &
      (SZ3_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      PRES(IPATH,J)=PRES(IPATH,J)-YP3* &
      (SZ3_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YP3=0


      YP2=YP2+FG2_D(IS,IPATH,KPROF(J))*FP2
      FG2=FG2+FP2*YP2_D(IS,IPATH,KPROF(J))
      FP2=0.


      PRES(IPATH,J)=PRES(IPATH,J)+YP2
      PRES(IPATH+1,J)=PRES(IPATH+1,J)+YP2* &
      (SZ2_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      PRES(IPATH,J)=PRES(IPATH,J)-YP2* &
      (SZ2_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YP2=0

!-----Path average temperature------------------------------------------

      DS12=DS12+ZTS*(5.D0*FT1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FT2_D(IS,IPATH,KPROF(J)) - FT3_D(IS,IPATH,KPROF(J)))/ &
      12.D0
      FT1=FT1+DS12_D(IS,IPATH,KPROF(J))*5.D0*ZTS/12.D0
      FT2=FT2+DS12_D(IS,IPATH,KPROF(J))*8.D0*ZTS/12.D0
      FT3=FT3-DS12_D(IS,IPATH,KPROF(J))*ZTS/12.D0
      DS23=DS23+ZTS*(-FT1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FT2_D(IS,IPATH,KPROF(J))+ &
      5.D0*FT3_D(IS,IPATH,KPROF(J)))/12.D0
      FT1=FT1-DS23_D(IS,IPATH,KPROF(J))*ZTS/12.D0
      FT2=FT2+DS23_D(IS,IPATH,KPROF(J))*8.D0*ZTS/12.D0
      FT3=FT3+DS23_D(IS,IPATH,KPROF(J))*5.D0*ZTS/12.D0

      YT3=YT3+FG3_D(IS,IPATH,KPROF(J))*FT3
      FG3=FG3+FT3*YT3_D(IS,IPATH,KPROF(J))
      FT3=0.

      TEMP(IPATH,J)=TEMP(IPATH,J)+YT3
      TEMP(IPATH+1,J)=TEMP(IPATH+1,J)+YT3* &
      (SZ3_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      TEMP(IPATH,J)=TEMP(IPATH,J)-YT3*(SZ3_D(IS,IPATH,KPROF(J))- &
      ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YT3=0


      YT2=YT2+FG2_D(IS,IPATH,KPROF(J))*FT2
      FG2=FG2+FT2*YT2_D(IS,IPATH,KPROF(J))
      FT2=0.

      TEMP(IPATH,J)=TEMP(IPATH,J)+YT2
      TEMP(IPATH+1,J)=TEMP(IPATH+1,J)+YT2* &
      (SZ2_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      TEMP(IPATH,J)=TEMP(IPATH,J)-YT2*(SZ2_D(IS,IPATH,KPROF(J))- &
      ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YT2=0



!-----Path gas amount---------------------------------------------------

      DS12=DS12+ZUS*(5.D0*FG1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FG2_D(IS,IPATH,KPROF(J)) - FG3_D(IS,IPATH,KPROF(J)))/ &
      12.D0
      FG1 =FG1+DS12_D(IS,IPATH,KPROF(J))*5.D0*ZUS/12.D0
      FG2=FG2+DS12_D(IS,IPATH,KPROF(J))*8.D0*ZUS/12.D0
      FG3=FG3-DS12_D(IS,IPATH,KPROF(J))*ZUS/12.D0
      DS23=DS23+ZUS*(-FG1_D(IS,IPATH,KPROF(J))+ &
      8.D0*FG2_D(IS,IPATH,KPROF(J))+ &
      5.D0*FG3_D(IS,IPATH,KPROF(J)))/12.D0
      FG1=FG1-DS23_D(IS,IPATH,KPROF(J))*ZUS/12.D0
      FG2=FG2+DS23_D(IS,IPATH,KPROF(J))*8.D0*ZUS/12.D0
      FG3=FG3+DS23_D(IS,IPATH,KPROF(J))*5.D0*ZUS/12.D0

      YA3=YA3+FA3
      FA3=0.
      IF (IGAS /= 91)THEN
        YA3=YA3+FG3
      FG3=0.
      END IF


      AIRNO(IPATH,J)=AIRNO(IPATH,J)+YA3
      AIRNO(IPATH+1,J)=AIRNO(IPATH+1,J)+YA3* &
      (SZ3_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      AIRNO(IPATH,J)=AIRNO(IPATH,J)-YA3*(SZ3_D(IS,IPATH,KPROF(J))- &
      ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YA3=0


      YA2=YA2+FA2
      FA2=0.

      IF (IGAS /= 91)THEN
        YA2=YA2+FG2
        FG2=0.
      END IF


      AIRNO(IPATH,J)=AIRNO(IPATH,J)+YA2
      AIRNO(IPATH+1,J)=AIRNO(IPATH+1,J)+YA2* &
      (SZ2_D(IS,IPATH,KPROF(J))-ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D
      AIRNO(IPATH,J)=AIRNO(IPATH,J)-YA2*(SZ2_D(IS,IPATH,KPROF(J))- &
      ALT_D(IPATH,KPROF(J)))/ &
      DELTA_ALT_D

      YA2=0


!-----Ray path length---------------------------------------------------

      DX23=DX23+DS23
      DS23=0.

      DX12=DX12+DS12
      DS12=0.

!-----Need the cosine of the angle--------------------------------------

      X3=X3+DX23
      X2=X2-DX23
      DX23=0.

      X2=X2+DX12
      X1=X1-DX12
      DX12=0.

      RR3=RR3+X3*CS3_D(IS,IPATH,KPROF(J))
      CS3=CS3+RR3_D(IS,IPATH,KPROF(J))*X3
      X3=0.

      RR2=RR2+X2*CS2_D(IS,IPATH,KPROF(J))
      CS2=CS2+RR2_D(IS,IPATH,KPROF(J))*X2
      X2=0.

      SN3=SN3+SN3_D(IS,IPATH,KPROF(J))*CS3/CS3_D(IS,IPATH,KPROF(J) &
      )
      CS3=0.

      SN2=SN2+SN2_D(IS,IPATH,KPROF(J))*CS2/CS2_D(IS,IPATH,KPROF(J) &
      )
      CS2=0.



!-----See if top of layer has been reached------------------------------

      CC=CC+SN3/RR3_D(IS,IPATH,KPROF(J))
      RR2=RR2-CC_D(IPATH,KPROF(J))*SN3/RR2_D(IS,IPATH,KPROF(J))**2
      SN3=0.

      CC=CC+SN2/RR2_D(IS,IPATH,KPROF(J))
      RR2=RR2-CC_D(IPATH,KPROF(J))*SN2/RR2_D(IS,IPATH,KPROF(J))**2
      SN2=0.

      Z3=Z3+RR3
      RR3=0.

      Z2=Z2+RR2
      RR2=0.

      Z3=Z3+SZ3
      SZ3=0.

      Z2=Z2+SZ2
      SZ2=0.

      Z1=Z1+Z2
      DELZ=DELZ+0.5D0*Z2
      Z2=0.

      IF (IS == NLOOP) THEN
        Z3=Z3+DELZ
        Z1=Z1-DELZ
        DELZ=0.

        B(IPATH,J)=B(IPATH,J)+Z3
        Z3=0.
      ENDIF

!-----------------------------------------------------------------------

      Z1=Z1+Z3
      DELZ=DELZ+Z3
      Z3=0.


      IF (DELZ_D(IS,IPATH,KPROF(J))<= DZMIN_D(IPATH,KPROF(J)))THEN
        DZMIN=DZMIN+DELZ
        DELZ=0.
      END IF

      CS1=CS1+SINCR_D(IPATH,KPROF(J))*DELZ
      SINCR=SINCR+DELZ*CS1_D(IS,IPATH,KPROF(J))
      DELZ=0.

   ENDDO

!-----Path average pressure---------------------------------------------

    YP1=YP1+FG1_D(1,IPATH,KPROF(J))*FP1
    FG1=FG1+FP1*YP1_D(1,IPATH,KPROF(J))

    PRES(IPATH,J)=PRES(IPATH,J)+YP1

    ZPS = 0.D0

!-----Path average temperature------------------------------------------

    YT1=YT1+FT1*FG1_D(1,IPATH,KPROF(J))
    FG1=FG1+FT1*YT1_D(1,IPATH,KPROF(J))

    TEMP(IPATH,J)=TEMP(IPATH,J)+YT1

    ZTS = 0.D0

!-----Path gas amount---------------------------------------------------

    YA1=YA1+FA1
    IF (IGAS /= 91)THEN
      YA1=YA1+FG1
    END IF
    AIRNO(IPATH,J)=AIRNO(IPATH,J)+YA1

    ZUS = 0.D0

!-----Change to intermediate variable X----------------------------------

    CS1=CS1+X1*RR1_D(1,IPATH,KPROF(J))
    RR1=RR1+X1*CS1_D(1,IPATH,KPROF(J))

!-----Cos of initial zenith angle---------------------------------------
    SN1=SN1-SN1_D(1,IPATH,KPROF(J))*CS1/CS1_D(1,IPATH,KPROF(J))


!-----Calculate snells law constant-------------------------------------

    Z1=Z1+SN1_D(1,IPATH,KPROF(J))*CC
    SN1=SN1+(Z1_D(1,IPATH,KPROF(J))+REARTH)*CC

!-----Distance from centre of earth-------------------------------------

    Z1=Z1+RR1

!-----Nominal path length set by user input-----------------------------

    SINCR=SINCR+1.D0*SINCR_D(IPATH,KPROF(J))*DZMIN/ &
    (Z1_D(1,IPATH,KPROF(J))+ REARTH)
    Z1=Z1-(0.5D0*SINCR_D(IPATH,KPROF(J))**2)*DZMIN/ &
    ((Z1_D(1,IPATH,KPROF(J)) + REARTH)**2)

    B(IPATH,J)=B(IPATH,J)+SINCR/NP

    Z1=Z1-SINCR/NP

!-----Sine of initial zenith angle--------------------------------------

    SN1 = 0.

!-----Lower boundary altitude-------------------------------------------

    A(IPATH,J)=A(IPATH,J)+SZ1

    A(IPATH,J)=A(IPATH,J)+Z1

  ENDDO
ENDDO


RETURN
END SUBROUTINE RTIASI_RAYTCEK
!     Atmospheric temps to Planck functions in many channels.
SUBROUTINE RTIASI_PLNCX(PRAD,KCHAN,PTEMP,KNCHPF)

!     Description:
!     Convert an array of atmospheric temperatures
!     to Planck functions in many channels
!
!     Method:
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code. J.R.Eyre. ECMWF.
!     2            12/8/92.    For version 2. Input satellite index and temps
!                              as vector of (chans*profs). J.R.Eyre. ECMWF.
!     3            18/1/1999.  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!
!     Module used:

USE RTIASI_IASCHN,ONLY : &
!     Imported scalars:
  ZPCON1,    & ! First planck function constant in mW/m**2/sr/cm-1
  ZPCON2,    & ! Second Planck function constant in K/cm**-1
!     Imported arrays:
  WVNUM      ! Wavenumber of IASI channel in cm**-1


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KNCHPF                ! Number of radiances.

!       Array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN(KNCHPF)         ! Array of channel indices.
REAL,INTENT(IN)     :: PTEMP(KNCHPF)         ! Temperatures.

!       Array arguments with intent out:
REAL,INTENT(OUT)    :: PRAD(KNCHPF)          ! Radiances.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,ICHAN
REAL :: BB
REAL :: F1


!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE PLANCK FUNCTIONS.
!-----------------------------------------------------------------------
DO J=1,KNCHPF
  ICHAN=KCHAN(J)
  BB=ZPCON1*WVNUM(ICHAN)**3
  F1=ZPCON2*WVNUM(ICHAN)/PTEMP(J)
  PRAD(J)=BB/(EXP(F1)-1.)
ENDDO
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_PLNCX
!     Atmospheric temps to Planck functions in many channels.
SUBROUTINE RTIASI_PLNCXK(PRAD_D,PRAD,KCHAN,PTEMP_D,PTEMP,KNCHPF)

!     Description:
!     K of subroutine PLNCX
!     Convert an array of atmospheric temperatures
!     to Planck functions in many channels
!
!     Method:
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code. J.R.Eyre. ECMWF.
!     2            12/8/92.    For version 2. Input satellite index and temps
!                              as vector of (chans*profs). J.R.Eyre. ECMWF.
!     3            18/1/1999.  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!
!     Module used:

USE RTIASI_IASCHN,ONLY : &
!     Imported scalars:
  ZPCON1,    & ! First planck function constant in mW/m**2/sr/cm-1
  ZPCON2,    & ! Second Planck function constant in K/cm**-1
!     Imported arrays:
  WVNUM      ! Wavenumber of IASI channel in cm**-1


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KNCHPF           ! Number of radiances.

!       Array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN(KNCHPF)    ! Array of channel indices.
REAL,INTENT(IN)     :: PRAD(KNCHPF)     ! Adjoint temperatures.
REAL,INTENT(IN)     :: PTEMP_D(KNCHPF)  ! Forward temperatures.
REAL,INTENT(IN)     :: PRAD_D(KNCHPF)   ! Forward radiances.

!       Array arguments with intent inout:
REAL,INTENT(INOUT)  :: PTEMP(KNCHPF)    ! K radiance.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,ICHAN
REAL    :: BB
REAL    :: F1
REAL    :: F1_D

!       Local arrays:



!-----End of header-----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE PLANCK FUNCTIONS.
!-----------------------------------------------------------------------
DO J=1,KNCHPF
  ICHAN=KCHAN(J)
  BB  =ZPCON1*WVNUM(ICHAN)**3
  F1_D=ZPCON2*WVNUM(ICHAN)/PTEMP_D(J)
  F1=-PRAD_D(J)**2*EXP(F1_D)*PRAD(J)/BB
  PTEMP(J)=PTEMP(J)+(-ZPCON2*WVNUM(ICHAN)/PTEMP_D(J)**2)*F1
ENDDO

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_PLNCXK
!     Radiance to brightness temperature in many channels.
SUBROUTINE RTIASI_BRIGV(PTB,PRAD,KCHAN,KNCHPF)

!    Description:
!    To convert an array of radiances in many channels
!    to brightness temperatures.
!
!     Method:
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code. J.R.Eyre. ECMWF.
!     2            12/8/92.    For version 2. Input satellite index and temps
!                              as vector of (chans*profs). J.R.Eyre. ECMWF.
!     3            18/1/1999.  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!
!     Module used:

USE RTIASI_IASCHN,ONLY : &
!     Imported scalars:
  ZPCON1,    & ! First planck function constant in mW/m**2/sr/cm-1
  ZPCON2,    & ! Second Planck function constant in K/cm**-1
!     Imported arrays:
  WVNUM      ! Wavenumber of IASI channel in cm**-1


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KNCHPF       ! Number of processed radiances.

!       Array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN(KNCHPF)! Array of channel indices.
REAL,INTENT(IN)     :: PRAD(KNCHPF) ! Radiances.

!       Array arguments with intent out:
REAL,INTENT(OUT)    :: PTB(KNCHPF)  ! Brightness temperatures.

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,ICHAN
REAL    :: T


!-----End of header----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE BRIGHTNESS TEMPERATURES
!-----------------------------------------------------------------------
DO J=1,KNCHPF
  ICHAN=KCHAN(J)
  T=ZPCON2*WVNUM(ICHAN)
    IF(PRAD(J) < 0)THEN
      PTB(J)=0.
    ELSE
      PTB(J)=T/ALOG(1.+ZPCON1*WVNUM(ICHAN)**3/PRAD(J))
    END IF
ENDDO
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_BRIGV
!     Radiance to brightness temperature in many channels.
SUBROUTINE RTIASI_BRIGVK(PTB_D,PTB,PRAD_D,PRAD,KCHAN,KNCHPF)

!     Description:
!     K of subroutine BRIGV.
!     To convert an array of radiances in many channels
!     to brightness temperatures.
!
!     Method:
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            6/12/90.    Original code. J.R.Eyre. ECMWF.
!     2            12/8/92.    For version 2. Input satellite index and temps
!                              as vector of (chans*profs). J.R.Eyre. ECMWF.
!     3            18/1/1999.  Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!
!     Module used:

USE RTIASI_IASCHN,ONLY : &
!     Imported scalars:
  ZPCON1,    & ! First planck function constant in mW/m**2/sr/cm-1
  ZPCON2,    & ! Second Planck function constant in K/cm**-1
!     Imported arrays:
  WVNUM      ! Wavenumber of IASI channel in cm**-1


IMPLICIT NONE


!     Subroutine arguments:

!       Scalar arguments with intent in:
INTEGER, INTENT(IN) :: KNCHPF        ! Number of processed radiances.

!       Direct array arguments with intent in:
INTEGER,INTENT(IN)  :: KCHAN(KNCHPF) ! Array of channel indices.
REAL,INTENT(IN)     :: PTB(KNCHPF)   ! Tangent linear Radiances.
REAL,INTENT(IN)     :: PRAD_D(KNCHPF)! Forward radiances.
REAL,INTENT(IN )    :: PTB_D(KNCHPF) ! Forward brightness temperatures.

!       Array arguments with intent out:
REAL,INTENT(INOUT)  :: PRAD(KNCHPF)  ! K brigh. temperature

!     End of subroutine arguments


!       Local scalars:
INTEGER :: J,ICHAN
REAL    :: T
REAL    :: T1
REAL    :: T2

!       Local arrays:
REAL    :: PRAD_AD(KNCHPF)


!-----End of header----------------------------------------------------


!-----------------------------------------------------------------------
!         1.   CALCULATE BRIGHTNESS TEMPERATURES
!-----------------------------------------------------------------------
DO J=1,KNCHPF
  ICHAN=KCHAN(J)
  T =ZPCON2*WVNUM(ICHAN)
  T1=ZPCON1*WVNUM(ICHAN)**3
  T2=1+T1/PRAD_D(J)
    IF(PRAD_D(J) > 0)THEN
      PRAD_AD(J)=PTB(J)*(PTB_D(J)**2/T)*T1/(T2*PRAD_D(J)**2)
    END IF
!        PRAD(J)=PRAD(J)+PRAD_AD(J)
   PRAD(J)=PRAD_AD(J)
ENDDO
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE RTIASI_BRIGVK
