AMREX_HOME ?= ../../amrex

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
PRECISION     = DOUBLE
BL_NO_FORT    = TRUE
USE_FFT       = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

VPATH_LOCATIONS   += ../src/analysis
INCLUDE_LOCATIONS += ../src/analysis

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
