# Makefile for sigmond

.SUFFIXES:  .o .cc .h .c
.PRECIOUS:  .cc .h .c

SRC_DIR=../../source

GRACE_INCDIR=/usr/local/grace-5.1.25/grace/include
GRACE_LIBDIR=/usr/local/lib64

MINUIT2_INCDIR=/usr/local/include
MINUIT2_LIBDIR=/usr/local/lib

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
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11  -g -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR) -Wl,-rpath=/usr/local/gcc-4.9.2/lib64
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11   -g -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -DNO_CXX11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11  -O3
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address  -fsanitize=undefined  -fstack-protector -DNO_MINUIT  -I$(GRACE_INCDIR) 
#COMMONFLAGS= -Wall -D$(PREC) -D$(NUMTYPE) -std=c++11 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address  -fsanitize=undefined  -fstack-protector -DNO_MINUIT  -I$(GRACE_INCDIR)  -Wl,-rpath=/usr/local/gcc-4.9.2/lib64

#LDFLAGS= -L$(GRACE_LIBDIR) -L$(MINUIT2_LIBDIR)
#LDFLAGS= -L$(GRACE_LIBDIR)
LDFLAGS=

LIBS=-lgrace_np -llapack -lMinuit2 -lgomp
#LIBS=-lgrace_np -llapack -lMinuit2 
#LIBS=-lgrace_np -llapack 


TASKS_SRCS = \
   rolling_pivot.cc      \
   scalar_defs.cc        \
   single_pivot.cc       \
   stopwatch.cc          \
   task_handler.cc       \
   task_check.cc         \
   task_fit.cc           \
   task_get_from_pivot.cc\
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
   rolling_pivot.h  \
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
   mcobs_handler.cc     \
   sampling_info.cc

ANAL_INCS = \
   bootstrapper.h       \
   histogram.h          \
   matrix.h             \
   mc_estimate.h        \
   mcobs_handler.h      \
   sampling_info.h

DATA_SRCS = \
   bins_info.cc            \
   corr_data_handler.cc    \
   ensemble_info.cc        \
   filelist_info.cc        \
   io_handler.cc           \
   obs_get_handler.cc      \
   qdp_byteorder.cc        \
   qdp_crc32.cc            \
   reweightings_handler.cc \
   vev_data_handler.cc

DATA_INCS = \
   array.h                \
   bins_handler.h         \
   bins_info.h            \
   corr_data_handler.h    \
   data_io_handler.h      \
   Endian.h               \
   ensemble_info.h        \
   filelist_info.h        \
   io_handler.h           \
   io_map.h               \
   laph_data_io_handler.h \
   obs_get_handler.h      \
   qdp_byteorder.h        \
   qdp_crc32.h            \
   reweightings_handler.h \
   samplings_handler.h    \
   vev_data_handler.h

OBS_SRCS = \
   basic_laph_operator_info.cc \
   diag_corr_set.cc            \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   gen_irrep_operator_info.cc  \
   mcobs_info.cc               \
   momenta.cc                  \
   operator_info.cc

OBS_INCS = \
   basic_laph_operator_info.h  \
   diag_corr_set.h             \
   correlator_info.h           \
   correlator_matrix_info.h    \
   encoder.h                   \
   gen_irrep_operator_info.h   \
   mcobs_info.h                \
   momenta.h                   \
   operator_info.h

PLOT_SRCS = \
   create_plots.cc        \
   grace_plot.cc

PLOT_INCS = \
   create_plots.h         \
   grace_plot.h

FIT_SRCS = \
   chisq_anisotropy.cc     \
   chisq_base.cc           \
   chisq_fit.cc            \
   chisq_tcorr.cc          \
   lmder.cc                \
   minimizer.cc            \
   minpack.cc              \
   model_tcorr.cc          \
   nl2sol.cc   

FIT_INCS = \
   chisq_anisotropy.h      \
   chisq_base.h            \
   chisq_fit.h             \
   chisq_tcorr.h           \
   minimizer.h             \
   minpack.h               \
   model_tcorr.h

SIGIOLIB_SRCS = \
   basic_laph_operator_info.cc \
   bins_info.cc                \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   ensemble_info.cc            \
   filelist_info.cc            \
   gen_irrep_operator_info.cc  \
   io_handler.cc               \
   matrix.cc                   \
   mcobs_info.cc               \
   momenta.cc                  \
   operator_info.cc            \
   qdp_byteorder.cc            \
   qdp_crc32.cc                \
   scalar_defs.cc              \
   stopwatch.cc                \
   xml_handler.cc

SIGQUERY_SRCS =                \
   basic_laph_operator_info.cc \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   gen_irrep_operator_info.cc  \
   io_handler.cc               \
   matrix.cc                   \
   mcobs_info.cc               \
   momenta.cc                  \
   operator_info.cc            \
   qdp_byteorder.cc            \
   qdp_crc32.cc                \
   scalar_defs.cc              \
   sigmond_query.cc            \
   xml_handler.cc

SIGLIB_SRCS = \
   basic_laph_operator_info.cc \
   bins_info.cc                \
   bootstrapper.cc             \
   corr_data_handler.cc        \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   ensemble_info.cc            \
   filelist_info.cc            \
   gen_irrep_operator_info.cc  \
   io_handler.cc               \
   matrix.cc                   \
   mc_estimate.cc              \
   mcobs_handler.cc            \
   mcobs_info.cc               \
   momenta.cc                  \
   obs_get_handler.cc          \
   operator_info.cc            \
   qdp_byteorder.cc            \
   qdp_crc32.cc                \
   reweightings_handler.cc     \
   sampling_info.cc            \
   scalar_defs.cc              \
   stopwatch.cc                \
   vev_data_handler.cc         \
   xml_handler.cc
