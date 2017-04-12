!+ Physical and other constants for global access in NWPSAF_1DVar

Module NWPSAFMod_Constants

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
!---------------------------------------------------------------------------
! Description:
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     30/09/98 Original code developed from GLOSSMod_Constants.
! 1.1              Andrew J. Smith
! 1.8     12/11/98 Updated thresholds for cloud test. AJS
! 1.9     18/11/98 Threshold for scattering index increased from 10 to 25. RJR
! 1.10    19/11/98 Removed ozone profile. AJSmith
! 1.16    06/01/99 New AMSU scattering index threshold value. AJSmith
! 1.25    08/03/99 Added reference temperature profile. AJSmith
! NWPSAFMod_Constants:
! 1.0     12/06/00 Original version based on ATOVSMod_Constants.  A Collard 
! 1.1     12/06/00 Moved cloud cost thresholds to NWPSAFMod_Params.  A Collard 
! 1.2     21/11/00 Decreased min_q from 3.e-06 to 1.e-08 (was above largest
!                  value that RTIASI could cope with!).  A. Collard.
!                                                        Met Office.
! 1.3     01/03/01 Corrected orbital heights added NOAA-16 R. Saunders
! 1.4     01/03/01 Added Pi                               M.D.E. Szyndel
! 3.0.1   15/07/03 Add Aqua and Geostationary orbit heights.   A.D. Collard.
! 3.0.4   04/03/04 Remove surplus parameters. Add Gastropod limits. A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
! 28        22/02/12 Added parameters for cloud liquid water 
!                    retrieval.   TR Sreerekha
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

!1) Missing data
!---------------
!Assigned to OPS/UM/VAR values


INTEGER, PARAMETER :: &
         MissData_I = -9999     ! integer missing data indicator

REAL,    PARAMETER :: &  
         MissData_R = -9999., &  ! real missing data indicator
         Tolerance  = 0.01       ! Margin of error when testing for missing data


!2) Miscellaneous constants
!--------------------------

REAL, PARAMETER :: &
  MaxSurfaceP         =   1200.0,    &  ! ( hPa )
  MinSurfaceP         =    300.0,    &  ! ( hPa )
  Min_q               =      1.0E-8, &  ! ( kg / kg )
  MaxTemperature      =    340.0,    &  ! ( K ) 
 ! MinTemperature      =     90.0,    &  ! ( K )
  MinTemperature      =     10.0,    &  ! ( K ) !PM MODIFICATION TO DECREASE THE MIN BT ALLOWED
  Upper_Transition    =     10.0,    &  ! Define transition region for
  Lower_Transition    =     30.0,    &  ! stratospheric extrapolation (hPa)
!  Pi                  =      3.14159265358979323846,&   ! Pi
  ScaleFactorP        =     1.0,   &    ! Additional cost function coefficient
                                        ! for Cloud Top Pressure
  ScaleFactorF        =     1.e12       ! Additional cost function coefficient
                                        ! for Cloud Cover
REAL, PARAMETER              :: RH_1=0.95 ! below RH_1 no Cloud liquid water
REAL, PARAMETER              :: RH_2=1.05 ! above RH_2, water vapor content remains fixed
REAL, PARAMETER              :: SplitQtotal=0.50 ! split factor for water vapor and cloud liquid water


REAL, PARAMETER :: InverseTol = TINY( 0.0 )

! 3. Universal constants 
REAL, PARAMETER    :: gravity        = 9.80665        ! accel due to grav m/s^2
REAL, PARAMETER    :: epsilon_1000   = 621.98         ! 1000 x epsilon (= Mw / Ma, the
                                                      !   ratio of the molecular weight of
                                                      !   water to that of air )
! RTIASI pressure levels and profile limits
! -----------------------------------------

REAL, PARAMETER :: RTIASI_ProfMin_Temperature(43) = &
! Temperature
(/ 190.9450, 215.5050, 212.8999, 198.4644, 198.4496, 204.5798,  &
   189.3978, 186.9405, 185.8896, 186.1667, 183.1527, 183.6486, 184.5035,  &
   183.7715, 186.2423, 185.1949, 186.2353, 189.8716, 189.5201, 189.3538,  &
   191.6970, 196.3361, 201.8183, 204.3071, 207.7184, 208.7572, 211.7302,  &
   216.1444, 220.3337, 222.7711, 226.1762, 230.0194, 233.5735, 235.5786,  &
   236.9657, 237.2966, 235.5742, 234.810, 234.7339, 235.1377, 236.1328,  &
   233.8223, 232.2004 /)

REAL, PARAMETER :: RTIASI_ProfMin_Humidity(43) = &
! Humidity (log[kg/kg])
(/ -12.8175, -12.6020, -12.5612, -12.5985, -12.7113, -12.7591, &  
   -12.8375, -12.9131, -12.9708, -12.9845, -13.0233, -13.1935, &
   -13.2826, -13.6400, -13.4100, -13.7933, -13.7887, -13.3353, &
   -13.5590, -13.7434, -13.4927, -13.1266, -13.1236, -13.1341, &
   -13.9211, -15.4624, -13.5051, -13.2664, -12.6056, -12.2045, &
   -12.2045, -11.2378, -10.5643, -9.67729, -9.44552, -9.55928, &
   -9.56623, -9.32946, -9.38823, -9.37936, -9.29087, -9.95508, &
   -10.5446 /)

REAL, PARAMETER :: RTIASI_ProfMin_Ozone(43) = &
! Ozone
(/ 0.6800000, 1.249224, 1.409621, 1.413208, 1.414658, 1.415650,  &
   1.415996, 1.416218, 1.416394, 1.416506, 1.360321, 0.5584311, 0.3678274,  &
   0.1907501, 0.1429634, 0.1146927, 4.6684265E-02, 7.8972578E-03,  &
   6.9316402E-03, 5.8482438E-03, 4.9787536E-03, 4.7907233E-03,  &
   8.0226958E-03, 7.8008175E-03, 1.0673031E-02, 5.3272247E-03,  &
   7.5416565E-03, 4.4134855E-03, 4.1854382E-03, 3.3203363E-03,  &
   1.8514395E-03, 1.6137362E-03, 1.5238980E-03, 1.3898611E-03,  &
   4.8571825E-04, 4.8571825E-04, 4.8571825E-04, 0.0000000E+00,  &
   0.0000000E+0, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00,  &
   0.0000000E+00 /)

REAL, PARAMETER :: RTIASI_ProfMax_Temperature(43) = &
! Temperature
(/ 288.1028,296.6161,310.956,309.3185,300.2095,289.6117, &
   290.7231,285.7654,264.9886,261.3101,257.8293,253.1519,251.5889, &
   254.8312,250.4318,243.4454,242.3198,241.5646,240.8297,242.2247, &
   243.8277,245.8346,253.0286,253.6444,259.7492,266.1804,271.1689, &
   273.4873,274.3341,275.3654,279.0348,282.4702,284.9223,289.6340, &
   294.6081,299.0551,303.4185,306.0771,306.4893,308.8958,310.5137, &
   311.5381,311.9639 /)

REAL, PARAMETER :: RTIASI_ProfMax_Humidity(43) = &
! Humidity (log[kg/kg])
(/ -12.5211, -12.4605, -12.4699, -12.4910, -12.5428, -12.5495, &
   -12.5251, -12.5273, -12.6098, -12.6045, -12.5690, -12.7110, &
   -12.7658, -12.7962, -12.8241, -12.8978, -12.9322, -10.4216, &
   -9.54228, -8.89603, -8.36804, -7.90317, -7.35874, -6.97959, &
   -6.57029, -6.15705, -5.80294, -5.48029, -5.21796, -5.04204, &
   -4.89851, -4.73497, -4.57232, -4.44394, -4.32548, -4.21415, &
   -4.09663, -4.02783, -3.91096, -3.79562, -3.76198, -3.75066, &
   -3.74429 /)

REAL, PARAMETER :: RTIASI_ProfMax_Ozone(43) = &
! Ozone
(/ 8.950000,9.244218,9.340413,9.371043,9.385607,9.392838, &
   9.396525,9.403474,9.433378,8.851357,8.767143,6.276594,6.064847, &
   5.385192,3.544674,2.912089,2.250107,1.991318,1.680061,1.226486, &
   0.8982038,0.8089733,0.5956802,0.4170179,0.3237114,0.2100811, &
   0.1716347,0.1342554,0.1286316,0.1264586,0.1167932,0.1103659, &
   0.1134226,0.1211433,0.1026807,0.1048491,9.9005103E-02, &
   9.4605148E-02,9.2247009E-02,9.0234280E-02,8.9395046E-02, &
   8.8834248E-02,8.8595532E-02 /)

REAL, PARAMETER :: Gastro_ProfMin_Temperature(43) = &
! Temperature
(/    171.312, 223.875, 226.167, 227.865, 222.880, 218.267, 206.125, &
      197.691, 191.106, 190.530, 187.971, 184.769, 185.525, 187.073, &
      188.609, 190.084, 191.502, 193.022, 193.943, 194.038, 193.335, &
      199.443, 203.401, 207.939, 212.309, 214.271, 215.836, 217.543, &
      220.291, 222.910, 226.219, 229.909, 233.355, 235.384, 237.524, &
      238.791, 239.796, 240.300, 239.511, 238.869, 238.396, 238.124, &
      238.012 /)

REAL, PARAMETER :: Gastro_ProfMin_Humidity(43) = &
! Humidity (log[kg/kg])
(/   -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, &
     -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, &
     -14.5111, -14.5111, -14.5111, -14.5111, -14.5111, -14.5114, -14.5160, &
     -14.1118, -13.7358, -13.2441, -13.4027, -13.5505, -13.2644, -13.0077, &
     -12.5583, -11.2047, -10.9886, -10.7779, -10.6132, -10.4834, -10.3468, &
     -10.1791, -10.0334, -9.83407, -9.65981, -9.78593, -9.91676, -10.0146, &
     -10.0561 /)

REAL, PARAMETER :: Gastro_ProfMin_Ozone(43) = &
! Ozone
(/    0.409180, 0.548767,  1.04839,  1.38661,  1.49426,  1.68436, &
      2.314520, 2.996720,  3.39560,  3.06009,  2.55029,  1.38349, &
      0.458024, 0.320123, 0.242680, 0.170826, 0.133424, 0.107896, &
     0.0841100, 0.0613907, 0.0461618, 0.0381580, 0.0337791, 0.0317365, &
     0.0307591, 0.0303822, 0.0300289, 0.0288931, 0.0270191, 0.0252648, &
     0.0236316, 0.0220547, 0.0206544, 0.0192475, 0.0182956, 0.0175253, &
     0.0168025, 0.0161618, 0.0158886, 0.0158479, 0.0158086, 0.0157865, &
     0.0157800 /)

REAL, PARAMETER :: Gastro_ProfMax_Temperature(43) = &
! Temperature
(/    243.625, 269.859, 280.628, 281.293, 281.078, 277.768, 271.529, &
      251.682, 248.584, 246.837, 245.370, 243.983, 243.993, 242.288, &
      238.996, 237.648, 236.968, 235.678, 234.432, 235.555, 237.713, &
      242.347, 250.202, 254.862, 257.603, 258.715, 260.381, 263.515, &
      268.018, 271.890, 275.317, 278.455, 281.956, 284.880, 287.260, &
      291.199, 295.486, 299.354, 301.303, 303.421, 306.694, 308.515, &
      309.271 /)

REAL, PARAMETER :: Gastro_ProfMax_Humidity(43) = &
! Humidity (log[kg/kg])
(/   -12.3253, -12.2527, -12.2105, -12.2109, -12.2109, -12.2109, -12.2109, & 
     -12.2109, -12.2109, -12.2108, -12.1821, -12.1315, -12.1041, -12.0585, &
     -11.6700, -11.0277, -10.5697, -10.0203, -9.43068, -8.93175, -8.50495, &
     -8.10161, -7.71420, -7.34577, -7.01238, -6.46244, -6.02537, -5.75955, &
     -5.49375, -5.28480, -5.11991, -4.89784, -4.72106, -4.55607, -4.42214, &
     -4.34550, -4.20359, -4.02166, -3.86965, -3.71388, -3.60369, -3.53592, &
     -3.51030 /)

REAL, PARAMETER :: Gastro_ProfMax_Ozone(43) = &
! Ozone
(/  1.32568,   2.61090,   4.88074,   7.02962,   9.87535,   12.7107,   &
    13.8342,   13.3531,   11.4415,   10.8979,   9.52897,   8.69141,   &
    7.26722,   5.81307,   5.04974,   4.38379,   3.75770,   3.13931,   &
    2.57147,   2.16578,   1.86273,   1.38077,   1.07591,   0.807083,  &
    0.561650,  0.358937,  0.296759,  0.256264,  0.218057,  0.182080,  &
    0.147922,  0.115395,  0.0937204, 0.0843766, 0.0758238, 0.0743263, &
    0.0724164, 0.0705681, 0.0687623, 0.0669211, 0.0663099, 0.0658188, &
    0.0655494 /)

END MODULE NWPSAFMod_Constants
