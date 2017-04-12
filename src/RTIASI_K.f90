SUBROUTINE RTIASI_K( &
     PROF_D,   & ! IN
     KCHAN,    & ! IN
     PANGL,    & ! IN
     PANGS,    & ! IN
     KSURF,    & ! IN
     KNPF,     & ! IN 
     KNCHPRO,  & ! IN 
     KLENPF,   & ! IN
     SatIndex, & ! IN
     IPRINT,   & ! IN 
     PRAD_D,   & ! OUT
     TB_D,     & ! OUT
     PROF_K)     ! OUT


!     Description:
!     Calculate Jacobians for IASI selected channels.
!
!     IASRTM--->RTIAI--->IASCF             ------------------------------|
!           |                                                            |
!           |                                                            |
!           |-->SETPRO                                                   |
!           |                                                            |
!           |                                                            |
!           |-->RTIASI-->PRFIN---->PRSLEV                                |
!                     |       |                                          |
!                     |       |--->PRFTAU                                |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |----------->OPDEP                                 |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |RTIASI
!                     |                                                  |SUITE
!                     |                                                  |
!                     |----------->RTTAU                                 |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |                                                  |
!                     |--->RTINT---->PLNCX                               |
!                               |                                        |
!                               |--->EMISS----->INTEG                    |
!                               |                                        |
!                               |--->BRIGV                               |
!                               |                                        |
!                               |--->LAYERS---->PRESHT--->INTEG          |
!                                          |                             |
!                                          |                             |
!                                          |--->RAYTCE-->INTER  ---------|
!
!
!  RTIAI:    Initialising routine for RTIASI.
!
!
!  IASCF:    Initialise satellite-dependent data for IASI.
!
!
!  SETPRO:   Set up input profile to RTIASI. Read user supplied emissivities.
!
!
!  RTIASI:   Perform radiative transfer computations for IASI.
!
!
!  PRFIN:    Set up profile variables for radiative transfer calculations.
!
!
!  PRSLEV:   Locate given pressures on array of fixed levels.
!
!
!  PRFTAU:   Store profile variables for transmittance calculations.
!
!
!  OPDEP:    Calculate layer optical depths.
!
!
!  RTTAU:    Calculate transmittance on levels of radiative transfer model.
!
!
!  RTINT:    Performs integration of radiative transfer equation in RTIASI.
!
!
!  PLNCX:    Converts atmospheric temperatures to PLANCK functions.
!
!
!  EMISS:    Computes surface emissivity for RT calculation.
!
!
!  BRIGV:    Converts radiances to brightness temperatures: if negative
!            radiances are encountered, brightness temperature is set
!            equal to zero.
!
!
!  LAYERS:   Produces Curtis-Godson weighted average layer temperatures.
!
!
!  PRESHT:   Calculates the heights of pressure levels.
!
!
!  INTEG:    Calculates a  weighted integrated quantity using
!            SIMPSONS rule.
!
!
!  INTER:    Interpolation routine.
!
!
!  RAYTCE:   Traces a ray through each atmospheric layer and calculates
!            the Curtis-Godson integrated temperatures.
!
!
!     Method:
!     1) See M.Matricardi and R.Saunders:A fast radiative model for
!        simulation of IASI radiances, Applied Optics, 1999, submitted.
!     2) See ECMWF TECH MEM 176. 
!
!     Input files:
!       1) IASIMIXCOEFPOS.dat:   Regression coefficients to compute positive 
!                                mixed gases transmittances; read in subroutine  
!                                IASCF.
!       2) IASIWVCOEFPOSRG1.dat: Regression coefficients to compute positive  
!                                water vapour transmittances: weak absorption 
!                                regime; read in subroutine IASCF.
!       3) IASIWVCOEFPOSRG2.dat: Regression coefficients to compute positive 
!                                water vapour transmittances: strong absorption 
!                                regime; read in subroutin IASCF.
!       4) IASIOZCOEF.dat:       Regression coefficients to compute ozone 
!                                transmittances; read in subroutine IASCF.
!       5) IASIGAMMAS.dat:       Gamma values to tune transmittances; read in
!                                subroutine IASCF.
!       Either
!       6a)IASIWVPOPC.dat:       To be used to compute emissivities over fresh
!                                water; read in subroutine IASCF.
!       or
!       6b)IASIWVPOPCOW.dat:     To be used to compute emissivities over sea;
!                                read in subroutine IASCF.
!       7) IASIEMS.dat:          User defined emissivities; read in subroutine 
!                                SETPRO.
!       8) PROFIN' '.dat         Input profile to RTIASI; read in subroutine 
!                                SETPRO. KNPF files PROFIN1.dat, PROFIN2.dat..
!                                ....PROFIN'KNPF'.dat are read.
!                                 
!
!     Output files:
!       1)WARNMESS' '.dat        Warning message file; written in subroutine
!                                SETPRO. KNPF files WARNMESS1.dat, WARNMESS2.dat
!                                .....WARNMESS'KNPF'.dat are written.
!       2)BRIGHTTEM' '.out       Computed brightness temperatures; written in
!                                subroutine RTIASI. KNPF files BRIGHTTEM1.out,
!                                BRIGHTTEM2.out.....BRIGHTTEM'KNPF'.out are 
!                                written.
!       3)RADIANCES' '.out       Computed radiances; written in subroutine 
!                                RTIASI. KNPF files RADIANCES1.out,
!                                RADIANCES2.out....RADIANCES'KNPF'.out are
!                                written.
!       4)IASITRS' '.out         Computed transmittances for an user selected
!                                set of IASI channels; written in subroutine
!                                RTTAU. KNPF files IASITRS1.out,IASITRS2.out
!                                ......IASITRS'KNPF'.out are written.        
!       5)IASIEMS.out            Computed emissivities; written in subroutine
!                                EMISS.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            1/02/1999   Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code"
!
!     Module used:

  USE NWPSAFMod_Params, ONLY : &
       StatusFatal

  IMPLICIT NONE
  
  INCLUDE 'Gen_ErrorReport.interface'

  !     Program arguments:

  INTEGER, INTENT(IN) :: KNPF
  INTEGER, INTENT(IN) :: KNCHPRO
  INTEGER, INTENT(IN) :: KLENPF
  REAL, INTENT(IN)    :: PROF_D(KLENPF,KNPF)  
  INTEGER, INTENT(IN) :: KCHAN(KNCHPRO*KNPF)
  REAL, INTENT(IN)    :: PANGL(KNPF)
  REAL, INTENT(IN)    :: PANGS(KNPF)
  INTEGER, INTENT(IN) :: KSURF(KNPF)
  INTEGER, INTENT(IN) :: SatIndex
  INTEGER, INTENT(IN) :: IPRINT
  REAL, INTENT(OUT)   :: PRAD_D(KNCHPRO*KNPF)
  REAL, INTENT(OUT)   :: TB_D(KNCHPRO*KNPF)
  REAL, INTENT(OUT)   :: PROF_K(KLENPF,KNCHPRO*KNPF) 

  !       Local parameters:

  CHARACTER (LEN=80)  :: ErrorMessage(2)  ! Message for Gen_ErrorReport
  CHARACTER (LEN=*), PARAMETER :: RoutineName = "RTIASI_Direct"

  INTEGER, PARAMETER :: KNCHTRS=20     ! Number of channels for whom
                                       ! transmitances will be 
                                       ! written to output file.

  !       Local scalars:
  INTEGER :: KNCHPF  ! Number of processed radiances.
  INTEGER :: J,I,K
  
  INTEGER :: ErrStatRep       ! error status for Gen_ErrorReport

  !       Local arrays:

  INTEGER :: KPROF(KNCHPRO*KNPF)
  INTEGER :: KCHTRSID(KNCHTRS*KNPF)   ! Array id of channels for
                                      ! whom transmittances will
                                      ! be written to output file.
                                      ! temperatures.
  !     End of program arguments


  !-----End of header-----------------------------------------------------

  KNCHPF=(KNCHPRO*KNPF)
  
  IF (KNPF /= 1) THEN
     ErrStatRep = StatusFatal
     Errormessage(1) = 'KNPF is not 1'
     CALL Gen_ErrorReport( RoutineName,  &
          ErrorMessage, &
          ErrStatRep    )
  ENDIF

  !-----------------------------------------------------------------------
  !         1.    SETS UP SOME INPUT PARAMETERS AND ARRAYS
  !-----------------------------------------------------------------------

  !-----1.1 Profile dependent quantities----------------------------------
  
  DO J=1,KNPF
     DO I=1,KNCHPRO
        KPROF(KNCHPRO*(J-1)+I)=J
     END DO


     DO K=1,KNCHTRS
        KCHTRSID(KNCHTRS*(J-1)+K)=K
     END DO
  END DO

  !-----------------------------------------------------------------------
  !   SETPRO performs checks on input profiles and reads in emissivities.
  !-----------------------------------------------------------------------
  CALL RTIASI_SETPRO(PROF_D,KLENPF,KNPF)

  !-----------------------------------------------------------------------
  !          4.    PERFORM TANGENT LINEAR RADIATIVE TRANSFER CALCULATIONS
  !-----------------------------------------------------------------------




!/////PROF_K(KLENPF,KNPF): Tangent linear Input vector to RTIASI.
!                           Must be set as follows:
!
! Elements 1  TO  43 : Input temperature profile  : [K].
! Elements 44 TO  86 : Input water vapour profile : Volume mixing ratio [ppmv].
! Elements 87 TO  129: Input ozone profile        : Volume mixing ratio [ppmv].
! Element         130: Total column cloud liquid
!                      water (not used at present): [mm].
! Element         131: Surface air temperature    : [K].
! Element         132: Surface air humidity       : Volume mixing ratio [ppmv].
! Element         133: Surface pressure           : [mb].
! Element         134: U Surface wind             : [m/sec].
! Element         135: V Surface wind             : [m/sec].
! Element         136: Skin temperature           : [K].
! Element         137: Cloud top pressure         : [mb].
! Element         138: Effective cloud coverage
!///////////////////////////////////////////////////////////////////////

  CALL RTIASIK(PROF_K,PROF_D,KNPF,KLENPF,PANGL,PANGS,KSURF,KCHAN, &
       KPROF,KNCHPF,SatIndex,TB_D,PRAD_D,IPRINT)

END SUBROUTINE RTIASI_K
