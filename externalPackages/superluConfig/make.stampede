############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = 

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = tmglib$(PLAT).a
SUPERLULIB   = libsuperlu$(PLAT).a
BLASDEF = -DUSE_VENDOR_BLAS
BLASLIB = -lblas
PREFIX = ${PROTEUS_PREFIX}
#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = mpicc -fPIC
CFLAGS       = -O2
FORTRAN	     = mpif77 -fPIC
FFLAGS       = -O2
LOADER       = mpif77 -fPIC
LOADOPTS     = 

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DUpCase, or -DAdd__)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = 

