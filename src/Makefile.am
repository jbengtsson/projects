SUBDIRS      = ptc

bin_PROGRAMS = main max_iv_1 vlad basic dynap leac dnu m4_1

#sc_f_LINK    = gfortran -g -O0 -Wall -fbounds-check -o sc_f

main_SOURCES     = main.cc
max_iv_1_SOURCES = max_iv_1.cc
vlad_SOURCES     = vlad.cc
basic_SOURCES    = basic.cc
dynap_SOURCES    = dynap.cc
dnu_SOURCES      = dnu.cc
leac_SOURCES     = leac.cc
m4_1_SOURCES     = m4_1.cc

#LIBS     = -L$(TRACY_LIB)/tracy -ltracy
# gfortran needs to be before stdc_++ for Fortran programs
LIBS	 = -L$(TRACY_LIB)/tracy/lib -ltracy \
	   -L$(NUM_REC)/lib -lnum_rec \
           -lstdc++ \
           -lm
#	   -lgfortran
#INCLUDES = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_FFLAGS   = -g -O0 -Wall -fbounds-check
AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend
#CXXFLAGS = -pg -O2 -Wall -Wno-non-template-friend
