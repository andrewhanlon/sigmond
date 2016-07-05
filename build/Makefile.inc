# Makefile for sigmond

.SUFFIXES:  .o .cc .h .c
.PRECIOUS:  .cc .h .c

SRC_DIR=/home/drew/git/research/lattice_qcd/sigmond/source

#GRACE_INCDIR=/usr/local/grace-5.1.25/grace/include
#GRACE_LIBDIR=/usr/local/lib64

#MINUIT2_INCDIR=/usr/local/include
#MINUIT2_LIBDIR=/usr/local/lib

CXX=g++
#CXX=/usr/local/gcc-4.9.2/bin/g++

  # uncomment one below to specify precision of data files
PREC=SINGLEPRECISION
#PREC=DOUBLEPRECISION

  # uncomment one below for real or complex type in data files
NUMTYPE=COMPLEXNUMBERS
#NUMTYPE=REALNUMBERS


#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR) -Wl,-rpath=/usr/local/gcc-4.9.2/lib64
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11   -g -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -DNO_CXX11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11  -O3

#LDFLAGS= -L$(GRACE_LIBDIR) -L$(MINUIT2_LIBDIR)
#LDFLAGS= -L$(GRACE_LIBDIR)
LDFLAGS=

LIBS=-lgrace_np -llapack -lMinuit2 -lgomp
#LIBS=-lgrace_np -llapack -lMinuit2 
#LIBS=-lgrace_np -llapack 


TASKS_SRCS = \
   scalar_defs.cc        \
   single_pivot.cc       \
   stopwatch.cc          \
   task_handler.cc       \
   task_check.cc         \
   task_fit.cc           \
   task_obsfunc.cc       \
   task_plot.cc          \
   task_print.cc         \
   task_rotate_corrs.cc  \
   task_utils.cc         \
   xml_handler.cc

TASKS_INCS = \
   args_handler.h   \
   log_helper.h     \
   multi_compare.h  \
   scalar_defs.h    \
   single_pivot.h   \
   stopwatch.h      \
   task_handler.h   \
   task_utils.h     \
   xml_handler.h

ANAL_SRCS = \
   bootstrapper.cc      \
   histogram.cc         \
   matrix.cc            \
   mc_estimate.cc       \
   mcobs_handler.cc

ANAL_INCS = \
   bootstrapper.h       \
   histogram.h          \
   matrix.h             \
   mc_estimate.h        \
   mcobs_handler.h

DATA_SRCS = \
   corr_data_handler.cc   \
   ensemble_info.cc       \
   filelist_info.cc       \
   io_handler.cc          \
   obs_get_handler.cc     \
   qdp_byteorder.cc       \
   qdp_crc32.cc           \
   vev_data_handler.cc

DATA_INCS = \
   array.h                \
   corr_data_handler.h    \
   data_io_handler.h      \
   ensemble_info.h        \
   filelist_info.h        \
   io_handler.h           \
   io_map.h               \
   obs_get_handler.h      \
   qdp_byteorder.h        \
   qdp_crc32.h            \
   vev_data_handler.h

OBS_SRCS = \
   diag_corr_set.cc          \
   correlator_info.cc        \
   correlator_matrix_info.cc \
   mcobs_info.cc             \
   momenta.cc                \
   operator_info.cc

OBS_INCS = \
   diag_corr_set.h  \
   correlator_info.h         \
   correlator_matrix_info.h  \
   mcobs_info.h              \
   momenta.h                 \
   operator_info.h

PLOT_SRCS = \
   create_plots.cc        \
   grace_plot.cc

PLOT_INCS = \
   create_plots.h         \
   grace_plot.h

FIT_SRCS = \
   chisq_base.cc           \
   chisq_fit.cc            \
   chisq_tcorr.cc          \
   lmder.cc                \
   minimizer.cc            \
   minpack.cc              \
   model_tcorr.cc          \
   nl2sol.cc   

FIT_INCS = \
   chisq_base.h            \
   chisq_fit.h             \
   chisq_tcorr.h           \
   minimizer.h             \
   minpack.h               \
   model_tcorr.h
