AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = PLUME-MoM
bindir=@bindir@

# list of sources for the binary
PLUMEMoM_SOURCES = \
		plume_model.f90

FCFLAGS         = -g -O3 -finit-real=snan
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lz -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

