MODULE RTIASI_ARRAYS

!
! Owner:
!     Andrew Collard
!
! Version Date     Comment
! ------- -------- -------
! 1.1     9/07/97 Original code.  ADC
! 2.0     2/3/01  Now a shadow of its former self.  ADC
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code"
!

  IMPLICIT NONE
 
CONTAINS

  !     Compute Lagrangian interpolation coefficients

  FUNCTION XL(MXDIM,KS,KF,K,ARG,ARR)

    !  Description:
    !  XL - To compute Lagrangian interpolation coefficients
    !  
    !  Method:
    !  
    !  Owner:
    !  Marco Matricardi.ECMWF
    !  
    !  History:
    !  Version      Date        Comment
    !  1            01/01/89    Original code. D.P. Edwards.
    !  2            01/02/1999  Marco Matricardi. ECMWF.
    !  
    !  Code description:
    !    Language:              Fortran 90.
    !    Software Standards:    "European Standards for Writing and Documenting
    !                           Exchangeable Fortran 90 code"


    IMPLICIT NONE

    !     Function arguments:

    !       Scalar arguments with intent in:
    INTEGER,INTENT(IN) :: MXDIM       ! Array dimension.
    INTEGER,INTENT(IN) :: KS          ! Lower limit of lagrangian sum.
    INTEGER,INTENT(IN) :: KF          ! Upper limit of lagrangian sum.
    REAL   ,INTENT(IN) :: ARG         ! Interpolation argument.
    
    !       Array arguments with intent in:
    REAL,INTENT(IN)    :: ARR(MXDIM)  ! Array to interpolate between.
    
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

  END FUNCTION XL

END MODULE RTIASI_ARRAYS


!     Parameters used in RTIASI suite
MODULE RTIASI_IASPAR

!     Description:
!     Set up parameters that define the maximum dimension of some arrays
!     used in the RTIASI suite.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".


IMPLICIT NONE


!       Global parameters:
INTEGER, PARAMETER :: JPCH=8461        ! Max. number of satellite coefficients
                                       ! stored.
INTEGER, PARAMETER :: JPMCH=8461       ! Max. number of IASI channels.
INTEGER, PARAMETER :: JPNSAT=1         ! Max. number of surface types.
INTEGER, PARAMETER :: JPLEV=43         ! Number of pressure levels.
INTEGER, PARAMETER :: JPCOFM=14        ! Max. number of mixed gas coefficients.
INTEGER, PARAMETER :: JPCOFW=14        ! Max. number of mixed gas coefficients.
INTEGER, PARAMETER :: JPCOFO=14        ! Max. number of ozone coefficients.
INTEGER, PARAMETER :: JPPF=1           ! Max. number of profiles to be
                                       ! processed.
INTEGER, PARAMETER :: JPCHPF=JPMCH*JPPF ! Max. number of processed radiances


!-----End of header----------------------------------------------------


END MODULE RTIASI_IASPAR
!     Constants for IASI channels
MODULE RTIASI_IASCHN

!     Description:
!     Set up constants for IASI channels.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPNSAT,    & ! Max. number of satellite coefficients stored.
  JPCH         ! Max. number of IASI channels.


IMPLICIT NONE


!     Module arguments:

!       Global scalars:
REAL, SAVE :: ZPCON1                ! First planck function constant in
                              ! mW/m^2/sr/cm-1.
REAL, SAVE :: ZPCON2                ! Second planck function constant in
                              ! K/cm^-1.

!       Global arrays:
REAL, SAVE    :: WVNUM(JPCH)        ! Wavenumber in cm^-1.
REAL, SAVE    :: GAMMA(JPCH,JPNSAT) ! "Gamma factor" transmittance
                              ! corrections.
REAL, SAVE    :: TRESMINM(JPCH)
REAL, SAVE    :: TRESMAXM(JPCH)
REAL, SAVE    :: TRESMINW(JPCH)
REAL, SAVE    :: TRESMAXW(JPCH)
!      COMPLEX :: WAOPC(JPCH)        ! Water optical constants.
REAL, SAVE    :: EMSCOEF(9,JPCH)

!     End of module arguments


!-----End of header-----------------------------------------------------


DATA ZPCON1/1.191044E-05/
DATA ZPCON2/1.438769/


END MODULE RTIASI_IASCHN
!     Transmittance coefficients
MODULE RTIASI_TAUCFN

!     Description:
!     Transmittance coefficients
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCH,      & ! Max. number of IASI channels.
  JPNSAT,    & ! Max. number of satellite coefficients stored.
  JPLEV,     & ! Number of pressure levels.
  JPCOFM,    & ! Max. number of mixed gas coefficients.
  JPCOFW,    & ! Max. number of water vapour coefficients.
  JPCOFO       ! Max. number of ozone coefficient.


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
!xxx REAL, SAVE :: CFMPOS(JPCOFM,JPLEV,JPCH,JPNSAT) ! Mixed gases coefficients
!xxx
!xxx REAL, SAVE :: CFWPOS(JPCOFW,JPLEV,JPCH,JPNSAT) ! Water vapour coefficients
!xxx                                              ! (positive transmittances).
!xxx REAL, SAVE :: CFWPOS1(JPCOFW,JPLEV,JPCH,JPNSAT) ! Water vapour coefficients
!xxx                                              ! (positive transmittances,
!xxx                                              ! weak absorption).
!xxx REAL, SAVE :: CFWPOS2(JPCOFW,JPLEV,JPCH,JPNSAT) ! Water vapour coefficients
!xxx                                              ! (positive transmittances,
!xxx                                              ! strong  absorption).
!xxx REAL, SAVE :: CFO(JPCOFO,JPLEV,JPCH,JPNSAT)  ! Ozone coefficients.

REAL, ALLOCATABLE :: CFMPOS(:,:,:,:)  ! Mixed gases coefficients
REAL, ALLOCATABLE :: CFWPOS(:,:,:,:)  ! Water vapour coefficients
                                      ! (positive transmittances).
REAL, ALLOCATABLE :: CFWPOS1(:,:,:,:) ! Water vapour coefficients
                                      ! (positive transmittances,
                                      ! weak absorption).
REAL, ALLOCATABLE :: CFWPOS2(:,:,:,:) ! Water vapour coefficients
                                      ! (positive transmittances,
                                      ! strong  absorption).
REAL, ALLOCATABLE :: CFO(:,:,:,:)     ! Ozone coefficients.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_TAUCFN
!     Satellite geometry constants
MODULE RTIASI_GEOCON

!     Description:
!     Set up Satellite geometry constants.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPNSAT     ! Max. number of satellite coefficients stored.


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
REAL, SAVE :: HTSAT(JPNSAT)  ! Satellite altitude in km
REAL, SAVE :: RATOE(JPNSAT)  ! Ratio satellite-orbit/earth radii

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_GEOCON
!     Satellite and solar geometry variables
MODULE RTIASI_GEOPTH

!     Description:
!     Satellite and solar geometry variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF,       & ! Max. number of profiles to be processed
  JPCHPF


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
REAL, SAVE :: XPATH(JPPF)  ! Secant of viewing path angle at surface
REAL, SAVE :: XPATH1(JPPF) ! XPATH - 1
REAL, SAVE :: SQTPTH(JPPF) ! Sqrt of XPATH
REAL, SAVE :: XPATHS(JPPF) ! Secant of solar zenith angle at surface

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_GEOPTH
!     Surface variables
MODULE RTIASI_SURF

!     Description:
!     Surface variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF       ! Max. number of profiles to be processed.


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
INTEGER, SAVE :: NLEVSF(JPPF)   ! Index of nearest std press level at/below
                          ! surface.
INTEGER, SAVE :: NSTYPE(JPPF)   ! Surface type index; 1=sea, 2=land.
REAL, SAVE   :: FRACPS(JPPF)   ! Fraction of std press level interval by which
                          ! surface is above level NLEVSF.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_SURF
!     Surface variables
MODULE RTIASI_SURFK

!     Description:
!     Tangent linear of module SURF
!     Surface variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF  ,      & ! Max. number of profiles to be processed.
  JPCHPF       ! Max. number of radiances to be processed.


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
INTEGER, SAVE :: NLEVSF(JPPF)   ! Index of nearest std press level at/below
                                ! surface.
INTEGER, SAVE :: NSTYPE(JPPF)   ! Surface type index; 1=sea, 2=land.
REAL, SAVE    :: FRACPS(JPCHPF) ! Fraction of std press level interval by
                                ! which surface is above level NLEVSF.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_SURFK
!     IR cloud variables
MODULE RTIASI_IRCLD

!     Description:
!     IR cloud variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF       ! Max. number of profiles to be processed


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
INTEGER, SAVE :: NLEVCD(JPPF)   ! Index of nearest std press level at/below
                                ! cloud top.
REAL, SAVE    :: FRACPC(JPPF)   ! Fraction of std press level interval by which
                                ! cloud is above level NLEVCD.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_IRCLD
!     IR cloud variables
MODULE RTIASI_IRCLDK

!     Description:
!     Tangent linear of module IRCLD
!     IR cloud variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCHPF       ! Max. number of radiances to be processed


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
INTEGER, SAVE :: NLEVCD(JPCHPF)   ! Index of nearest std press level at/below
                                  ! cloud top.
REAL, SAVE    :: FRACPC(JPCHPF)   ! Fraction of std press level interval by
                                  ! which cloud is above level NLEVCD.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_IRCLDK
!     Reference profile in transmittance calculations
MODULE RTIASI_PRFREF

!     Description:
!     Reference profile in transmittance calculations
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPLEV      ! Number of pressure levels


IMPLICIT NONE


!     Module arguments

!       Global arrays
REAL, SAVE :: TREF(JPLEV) ! Reference temperature profile in K.
REAL, SAVE :: WREF(JPLEV) ! Reference water vapour volume mixing ratio profile
                          ! in ppmv.
REAL, SAVE :: OREF(JPLEV) ! Reference ozone volume mixing ratio profile 
                          ! in ppmv.
REAL, SAVE :: TREFO(JPLEV)! Reference temperature profile in K (to be used for
                          ! the computation of ozone transmittances).
REAL, SAVE:: WREFO(JPLEV) ! Reference water vapour profile in ppmv (to be used
                          ! for the computation of ozone transmittances).

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_PRFREF
!     Atmospheric profile constants
MODULE RTIASI_PRFCON

!     Description:
!     Atmospheric profile constants
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPLEV      ! Number of pressure levels


IMPLICIT NONE


!     Module arguments:

!       Global scalars:
INTEGER, SAVE :: NLEVW   ! Upper level for water vapour transmittance
                         ! calculation (assumed transmittance=1. above
                         ! this).
REAL, SAVE :: UPPRES     ! Upper pressure level for transmittance
                         ! calculation.

!       Global arrays:
REAL, SAVE :: XPRES(JPLEV)  ! Standard pressure levels for transmittance
                            ! (and, currently, radiative transfer)
                            ! calculation; from top down; in mb=hpa.
REAL, SAVE :: XPRES2(JPLEV) ! XPRES**2; in mb^2.
REAL, SAVE :: DPRES(JPLEV)  ! Intervals between standard pressure levels;
                            ! in mb.
REAL, SAVE :: DPP(JPLEV)    ! DPRES*XPRES.

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_PRFCON
!     IR surface emissivities over sea
MODULE RTIASI_EMISIR

!     Description:
!     IR surface emissivities over sea
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCH,              & ! Max. number of radiances to be processed
  JPPF


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
REAL, SAVE :: EMSIR(JPCH*JPPF) !Infrared surface emissivities over land

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_EMISIR
!     IR surface emissivities over sea
MODULE RTIASI_EMISIRK

!     Description:
!     Tangent linear of module EMISIR
!     IR surface emissivities over sea
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCHPF              ! Max. no. of IASI channels


IMPLICIT NONE


!     Module arguments:

!       Global arrays:
REAL, SAVE :: EMSIR(JPCHPF) !Infrared surface emissivities over land

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_EMISIRK
!     Profile variables
MODULE RTIASI_PRFVAR

!     Description:
!     Profile variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF,      & ! Max. number of profiles to be processed
  JPCHPF, &
  JPLEV      ! Number of pressure levels


IMPLICIT NONE


!     Module arguments:

!       Global arrays
REAL, SAVE :: TEMPW(JPLEV,JPPF) ! Curtis-Godson Temperature profile in K.
REAL, SAVE :: TEMP(JPLEV,JPPF)  ! Temperature profile in K.
REAL, SAVE :: WMIX(JPLEV,JPPF)  ! Water vapour volume mixing ratio in ppmv.
REAL, SAVE :: TA(JPPF)          ! Surface air temperature in K.
REAL, SAVE :: WMIXS(JPPF)       ! Water vapour surface volume mixing ratio
                                ! in ppmv.
REAL, SAVE :: TS(JPPF)          ! Surface skin temperature in K.
REAL, SAVE :: SURFP(JPPF)       ! Surface pressure in mb=hPa.
REAL, SAVE :: CLDP(JPPF)        ! Cloud-top pressure in mb=hPa.
REAL, SAVE :: CLDF(JPPF)        ! Fractional (ir) cloud cover.
REAL, SAVE :: CLW(JPPF)         ! Total column cloud liquid water in mm.
REAL, SAVE :: SURFWU(JPPF)      ! U component of surface wind in m/sec.
REAL, SAVE :: SURFWV(JPPF)      ! V component of surface wind in m/sec.
REAL, SAVE :: OMIX(JPLEV,JPPF)  ! Ozone volume mixing ratio in ppmv

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_PRFVAR
!     Profile variables
MODULE RTIASI_PRFVARK

!     Description:
!     Tangent linear of module PRFVAR
!     Profile variables
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPCHPF,    & ! Max. number of profiles to be processed
  JPLEV      ! Number of pressure levels


IMPLICIT NONE


!     Module arguments:

!       Global arrays
REAL, SAVE :: TEMPW(JPLEV,JPCHPF) ! Curtis-Godson temperature profile.
REAL, SAVE :: TEMP(JPLEV,JPCHPF)  ! Temperature profile in K.
REAL, SAVE :: WMIX(JPLEV,JPCHPF)  ! Water vapour volume mixing ratio in ppmv.
REAL, SAVE :: TA(JPCHPF)          ! Surface air temperature in K.
REAL, SAVE :: WMIXS(JPCHPF)       ! Water vapour surface volume mixing ratio
                                  ! in ppmv.
REAL, SAVE :: TS(JPCHPF)          ! Surface skin temperature in K.
REAL, SAVE :: SURFP(JPCHPF)       ! Surface pressure in mb=hPa.
REAL, SAVE :: CLDP(JPCHPF)        ! Cloud-top pressure in mb=hPa.
REAL, SAVE :: CLDF(JPCHPF)        ! Fractional (ir) cloud cover.
REAL, SAVE :: CLW(JPCHPF)         ! Total column cloud liquid water in mm.
REAL, SAVE :: SURFWU(JPCHPF)      ! U component of surface wind in m/sec.
REAL, SAVE :: SURFWV(JPCHPF)      ! V component of surface wind in m/sec.
REAL, SAVE :: OMIX(JPLEV,JPCHPF)  ! Ozone volume mixing ratio in ppmv

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_PRFVARK
!     Profile dependent variables
MODULE RTIASI_PREDVAR

!     Description:
!     Profile dependent variables used to set up
!     IASI predictors.
!
!     Owner:
!     Marco Matricardi
!
!     History:
!     Version      Date        Comment
!     1            01/02/1999  Original code. Marco Matricardi. ECMWF.
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".

!     Module used

USE RTIASI_IASPAR,ONLY : &
!     Imported parameters:
  JPPF,      & ! Max. number of profiles to be processed
  JPLEV      ! Number of pressure levels


IMPLICIT NONE


!     Module arguments:

!       Global parameters:
REAL, SAVE :: RATT(JPLEV,JPPF)
REAL, SAVE :: RATO(JPLEV,JPPF)
REAL, SAVE :: RATW(JPLEV,JPPF)
REAL, SAVE :: PWTR(JPLEV,JPPF)
REAL, SAVE :: PWOR(JPLEV,JPPF)
REAL, SAVE :: PWORTR(JPLEV,JPPF)
REAL, SAVE :: DELT(JPLEV,JPPF)
REAL, SAVE :: DELTO(JPLEV,JPPF)
REAL, SAVE :: RATWO(JPLEV,JPPF)
REAL, SAVE :: RATTO(JPLEV,JPPF)
REAL, SAVE :: PUWOR(JPLEV,JPPF)

!     End of module arguments


!-----End of header-----------------------------------------------------


END MODULE RTIASI_PREDVAR
