# Define directories containing source 
SHELL = /bin/sh

# Define FORTRAN compiler
FC = gfortran
#FC = nagfor
#FC = pgf90
#FC = gfortran
#FC = xlf95_r

# Where to find the RTM. Replace this with your own path to the object files and modules
#RTM_DIR = /data/local/frpw/fcm/RTTOV/11.1
#RTM_INC = $(RTM_DIR)/include
#RTM_MOD = $(RTM_DIR)/mod 
#RTM_LIB = $(RTM_DIR)/lib
#RTMODEL = rttov11.1.0

#PM MODIF
#RTM_DIR = /home/francesco/Documenti/radiative_trasfer/rttov112_GB_K_nohdf
RTM_DIR=/home/martinet/rttov/rttov_visee_sol/rttov112_GB_K_v1.3
#RTM_DIR = /home/martinet/rttov/rttov_visee_sol/rttov102_r1094
#RTM_DIR = /home/martinet/rttov/rttov_10-2_r916/
RTM_INC = $(RTM_DIR)/include
RTM_MOD = $(RTM_DIR)/mod 
RTM_LIB = $(RTM_DIR)/lib
RTMODEL = rttov11.2.0

# Compiler and linker flags depend on RTTOV version
# Compiler and linker flags depend on RTTOV version
ifeq ($(findstring rttov1,$(RTMODEL)),rttov1)
  # RTTOV10 or later
  NCDF_LIB = /data/local/frpw/netcdf-4.1.1/build/ifort12/lib
  HDF5_LIB = /data/nwp1/frhg/hdf5-1.8.8/build/ifort12/lib
  # Preprocesssor flags
  CPPFLAGS = -I$(RTM_MOD) -I$(RTM_INC)
  # If reading in HDF5 format coefficients uncomment the line below
  #CPPFLAGS = -D_RTTOV_HDF -I$(RTM_MOD) -I$(RTM_INC)
  # Compiler flags depend on compiler
  ifeq ($(FC),ifort)
    CFLAGS = -cpp -g -r8 -warn all $(CPPFLAGS) 
  else ifeq ($(FC),nagfor)
    CFLAGS = -gline -C -fpp -kind=byte $(CPPFLAGS)
  else ifeq ($(FC),pgf90)
    CFLAGS = -C -Mfree -Mpreprocess $(CPPFLAGS)
  else ifeq ($(FC),gfortran)
#   CFLAGS= -g -O3 -fbounds-check $(CPPFLAGS)
   CFLAGS = -g -cpp -fdefault-real-8 -ffree-form $(CPPFLAGS)
#    CFLAGS= -g -cpp -fPIC -fdefault-real-8 -fbounds-check -finit-integer=-1 -finit-real=nan -finit-character=127 \
   #         -pedantic -Wall -fbacktrace -fcheck-array-temporaries -ffree-line-length-none -ffree-form $(CPPFLAGS)
  endif

  # Linker flags
  # Linker flags if RTTOV hasn't been compiled with HDF
  #LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_coef_io -l$(RTMODEL)_emis_atlas \
       # PM     -l$(RTMODEL)_main -L$(NCDF_LIB) -lnetcdf -lnetcdff
  # Linker flags if RTTOV has been compiled with HDF
  #LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_coef_io -l$(RTMODEL)_emis_atlas \
     # PM      -l$(RTMODEL)_hdf -l$(RTMODEL)_main -L$(NCDF_LIB) -lnetcdf -lnetcdff \
    # PM       -L$(HDF5_LIB) -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran

  # Comment out the above linker flags and uncomment the lines below when 
  # not using an emissivity atlas and netcdf library
  # Linker flags if RTTOV hasn't been compiled with HDF
  LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_other -l$(RTMODEL)_coef_io -l$(RTMODEL)_main 
  # Linker flags if RTTOV has been compiled with HDF
 # PM LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_coef_io -l$(RTMODEL)_hdf \
  # PM         -l$(RTMODEL)_main -L$(HDF5_LIB) -lhdf5 -lhdf5_fortran -lhdf5_hl \
  # PM         -lhdf5hl_fortran

#LFLAGS = -L$(RTM_LIB) -lrttov10.2.0_coef_io -l$(RTMODEL)_main -lrttov10.2.0_coef_io -lrttov10.2.0_other

else
  # RTTOV9 and before

  # Compiler flags depend on compiler
  ifeq ($(FC),ifort)
    CFLAGS = -g -r8 -warn all -I$(RTM_MOD) -I$(RTM_INC)
  else ifeq ($(FC),nagfor)
    CFLAGS = -gline -C -I$(RTM_MOD) -I$(RTM_INC)
  else ifeq ($(FC),pgf90)
    CFLAGS = -C -I$(RTM_MOD) -I$(RTM_INC)
  else ifeq ($(FC),gfortran)
    CFLAGS = -g -fdefault-real-8 -ffree-form -I$(RTM_MOD) -I$(RTM_INC)
  endif

  # Linker flags
  LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)

endif

# If compiling with the XL Fortran compiler on an AIX platform comment out 
# all of the if blocks above and uncomment the lines below

# Compiler flags
#CFLAGS = -pg -qsuffix=cpp=f90 -qextname -I$(RTM_MOD) -I$(RTM_INC)

# Linker flags 
#LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_coef_io -l$(RTMODEL)_emis_atlas \
   # PM      -l$(RTMODEL)_main -L$(NCDF_LIB) -lnetcdf

# Comment out the above linker flags and uncomment the lines below when 
# not using an emissivity atlas and netcdf library
#LFLAGS = -L$(RTM_LIB) -l$(RTMODEL)_coef_io -l$(RTMODEL)_main


SOURCE = NWPSAF_1DVar_Driver.f90

OBJS  =	NWPSAFMod_Channellist.o \
	NWPSAFMod_Constants.o \
	NWPSAFMod_CovarianceMatrices.o \
	NWPSAFMod_RTmodel.o \
	NWPSAFMod_RTTOV8.o \
	NWPSAFMod_RTTOV9.o \
	NWPSAFMod_RTTOV10.o \
	NWPSAFMod_RTTOV11.o \
	NWPSAFMod_Params.o \
	NWPSAFMod_ObsInfo.o \
	NWPSAFMod_SatVapor.o \
	Gen_ErrorReport.o \
	Gen_FreeUnit.o \
	Gen_GetUnit.o \
	Gen_MessageReport.o \
	NWPSAF_1DVar.o \
	NWPSAF_1DVar_Driver.o \
	NWPSAF_AdditionalCost_Cloud.o \
	NWPSAF_BandInverse.o \
	NWPSAF_BandMultiply.o \
	NWPSAF_Calculate_Cost_Function.o \
	NWPSAF_Channellist.o \
	NWPSAF_CheckIteration.o \
	NWPSAF_Check_Temperatures.o \
	NWPSAF_Cholesky.o \
        NWPSAF_CO2Slice.o \
	NWPSAF_CloudCost.o \
	NWPSAF_CloudyOrNot.o \
	NWPSAF_CloudStructure.o \
	NWPSAF_DeAllocate.o \
	NWPSAF_Fastmodel_Interface.o \
        NWPSAF_Gastro_Interface_Dummy.o \
	NWPSAF_RTTOV7_Interface_Dummy.o \
	NWPSAF_RTTOV8_Interface_Dummy.o \
	NWPSAF_RTTOV9_Interface_Dummy.o \
        NWPSAF_RTTOV11_Interface_nonetcdf.o \
	NWPSAF_RTTOV10_Interface_Dummy.o \
	NWPSAF_RTIASI_Interface_Dummy.o \
        NWPSAF_RTTOV11_Brute_Force_Jacobians.o \
	NWPSAF_InitBmatrix.o \
	NWPSAF_InitRmatrix.o \
	NWPSAF_Initialise.o \
	NWPSAF_Layers_to_LWP.o \
	NWPSAF_LWP_to_Layers.o \
	NWPSAF_LayerK_to_LWPK.o \
	NWPSAF_Qtot_to_q_ql.o \
	NWPSAF_Minimize.o \
	NWPSAF_Minimize_100.o \
	NWPSAF_Minimize_100ML.o \
	NWPSAF_Minimize_101.o \
	NWPSAF_OpenFile.o \
	NWPSAF_ProcessData.o \
	NWPSAF_RMatrix_ChanSelect.o \
	NWPSAF_Read_Background.o \
	NWPSAF_Read_ControlData.o \
	NWPSAF_ReadHeaders.o \
	NWPSAF_Read_Observations.o \
	NWPSAF_SatMatInv.o \
	NWPSAF_SetUpBackground.o \
	NWPSAF_SetUpRetrievals.o \
	NWPSAF_StratosExtrap.o \
	NWPSAF_TranslateDataIn.o \
	NWPSAF_TranslateDataOut.o \
	NWPSAF_IntegerSort.o \
	QSAT.o \
	NWPSAF_svp.o \
	NWPSAF_svp_deriv.o \
	RTIASI_Arrays.o \
	RTIASI_Direct.o \
	RTIASI_K.o \
	RTIASI_Subroutines.o \
	RTIASI_rttvi.o

MODULES = NWPSAFMod_Constants.o \
	NWPSAFMod_Params.o \
	NWPSAFMod_CovarianceMatrices.o \
	NWPSAFMod_Channellist.o \
	NWPSAFMod_RTmodel.o \
	NWPSAFMod_ObsInfo.o \
    RTIASI_Arrays.o \
	NWPSAFMod_RTTOV8.o \
	NWPSAFMod_RTTOV9.o \
    NWPSAFMod_RTTOV10.o \
    NWPSAFMod_RTTOV11.o


# Allowed SUFFIXES:

.SUFFIXES:.f90 .o .mod

# Object dependencies:

NWPSAF_1DVar : $(OBJS) $(MODULES) 
	$(FC) -o NWPSAF_1DVar $(OBJS) $(LFLAGS)

.f90.o : $<
	$(FC) -c $(CFLAGS) $<

# Here are the modules 

NWPSAF_1DVar_Driver.o : $(SOURCE) $(MODULES) 
	$(FC) -c $(CFLAGS) $(INCFLAG) $(SOURCE)

NWPSAFMod_Constants.o : NWPSAFMod_Constants.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_Constants.f90

NWPSAFMod_Params.o : NWPSAFMod_Params.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_Params.f90 

NWPSAFMod_CovarianceMatrices.o : NWPSAFMod_CovarianceMatrices.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_CovarianceMatrices.f90

NWPSAFMod_Channellist.o : NWPSAFMod_Channellist.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_Channellist.f90

NWPSAFMod_RTmodel.o : NWPSAFMod_RTmodel.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_RTmodel.f90

NWPSAFMod_ObsInfo.o : NWPSAFMod_ObsInfo.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_ObsInfo.f90

NWPSAFMod_SatVapor.o : NWPSAFMod_SatVapor.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_SatVapor.f90

RTIASI_Arrays.o : RTIASI_Arrays.f90
	$(FC) -c $(CFLAGS) RTIASI_Arrays.f90

NWPSAFMod_RTTOV8.o : NWPSAFMod_RTTOV8.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_RTTOV8.f90

NWPSAFMod_RTTOV9.o : NWPSAFMod_RTTOV9.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_RTTOV9.f90

NWPSAFMod_RTTOV10.o : NWPSAFMod_RTTOV10.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_RTTOV10.f90

NWPSAFMod_RTTOV11.o : NWPSAFMod_RTTOV11.f90
	$(FC) -c $(CFLAGS) NWPSAFMod_RTTOV11.f90


clean:
	rm -f *.o *.bak *~ NWPSAF_1DVar *.mod *genmod.f90
