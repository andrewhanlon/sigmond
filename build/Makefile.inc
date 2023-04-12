# Makefile for sigmond

.SUFFIXES:  .o .cc .h .c
.PRECIOUS:  .cc .h .c

SRC_DIR=../../source

#GRACE_INCDIR=/usr/local/include
#GRACE_LIBDIR=/usr/local/lib64
#GRACE_INCDIR=/usr/local/grace-5.1.25/grace/include
GRACE_INCDIR=/usr/include
#GRACE_LIBDIR=/usr/lib64

MINUIT2_INCDIR=/usr/local/include
#MINUIT2_INCDIR=/usr/include/Minuit2
MINUIT2_LIBDIR=/usr/local/lib
#MINUIT2_LIBDIR=/usr/lib

CXX=g++
#CXX=/usr/local/gcc-4.9.2/bin/g++
#CXX=/usr/local/apps/gcc/gcc-9.1.0/bin/g++

  # uncomment one below to specify precision of data files
  #PREC=SINGLEPRECISION
PREC=DOUBLEPRECISION

  # uncomment one below for real or complex type in data files
NUMTYPE=COMPLEXNUMBERS
#NUMTYPE=REALNUMBERS

  # specify location of default XML file containing known ensembles (do not use quotes)
KNOWN_ENSEMBLES_FILE=/home/sarahski/latticeQCD/D200_analysis/ensembles_simon.xml
#KNOWN_ENSEMBLES_FILE=/home/colin/research/sigmond/build/batch/ensembles.xml

  # uncomment one below for fstreams or HDF5 file as default file format
IOTYPE=DEFAULT_FSTREAM
#IOTYPE=DEFAULT_HDF5

COMMONFLAGS= -Wall -std=c++17  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -std=c++11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -std=c++11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR) -Wl,-rpath=/usr/local/apps/gcc/gcc-9.1.0/lib64
#COMMONFLAGS= -Wall -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -std=c++11 -DNO_MINUIT  -O3 -I$(GRACE_INCDIR)
#COMMONFLAGS= -Wall -std=c++11   -g -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -DNO_CXX11  -O3 -I$(GRACE_INCDIR) -I$(MINUIT2_INCDIR)
#COMMONFLAGS= -Wall -std=c++11 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address  -fsanitize=undefined  -fstack-protector -DNO_MINUIT  -I$(GRACE_INCDIR) 
#COMMONFLAGS= -Wall -std=c++11 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address  -fsanitize=undefined  -fstack-protector -DNO_MINUIT  -I$(GRACE_INCDIR)  -Wl,-rpath=/usr/local/gcc-4.9.2/lib64
COMMONFLAGS= -Wall -std=c++11 -O3 -I$(MINUIT2_INCDIR)

LDFLAGS= -L$(GRACE_LIBDIR) -L$(MINUIT2_LIBDIR)
#LDFLAGS= -L$(GRACE_LIBDIR)

#LIBS=-lgrace_np -llapack -lMinuit2 -lgomp -lMinuit2Math
#LIBS=-lgrace_np -llapack -lMinuit2 
#LIBS=-lgrace_np -llapack -lMinuit2 -lMinuit2Math -lhdf5 -lgomp
LIBS=-lgrace_np -llapack -lMinuit2 -lhdf5 -lgomp
#LIBS=-lgrace_np -llapack -lMinuit2 -lhdf5 
#LIBS=-lgrace_np -llapack 

#####################################################################

DCFLAGS = -DH5_USE_110_API -D$(PREC) -D$(NUMTYPE) -D$(IOTYPE) -DDEFAULTENSFILE="\""$(KNOWN_ENSEMBLES_FILE)"\""

DRIVER_SRC = sigmond.cc

TASKS_SRCS = \
   scalar_defs.cc        \
   single_pivot.cc       \
   rolling_pivot.cc       \
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
   scalar_defs.h    \
   pivoter.h   \
   single_pivot.h   \
   rolling_pivot.h   \
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
   bins_info.cc           \
   byte_handler.cc        \
   corr_data_handler.cc   \
   ensemble_info.cc       \
   filelist_info.cc       \
   io_handler_fstream.cc  \
   io_handler_hdf5.cc     \
   obs_get_handler.cc     \
   vev_data_handler.cc

DATA_INCS = \
   array.h                \
   bins_handler.h         \
   bins_info.h            \
   byte_handler.h         \
   corr_data_handler.h    \
   data_io_handler.h      \
   ensemble_info.h        \
   filelist_info.h        \
   io_handler_fstream.h   \
   io_handler_hdf5.h      \
   io_map.h               \
   io_map_base.h          \
   io_map_fstream.h       \
   io_map_hdf5.h          \
   laph_data_io_handler.h \
   obs_get_handler.h      \
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
   chisq_disp.cc     \
   chisq_fit.cc            \
   chisq_tcorr.cc          \
   chisq_logtcorr.cc       \
   lmder.cc                \
   minimizer.cc            \
   minpack.cc              \
   model_logtcorr.cc       \
   model_tcorr.cc          \
   nl2sol.cc   

FIT_INCS = \
   chisq_anisotropy.h      \
   chisq_base.h            \
   chisq_disp.h      \
   chisq_fit.h             \
   chisq_tcorr.h           \
   minimizer.h             \
   minpack.h               \
   model_logtcorr.h        \
   model_tcorr.h

SIGIOLIB_SRCS = \
   basic_laph_operator_info.cc \
   bins_info.cc                \
   byte_handler.cc             \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   ensemble_info.cc            \
   filelist_info.cc            \
   gen_irrep_operator_info.cc  \
   io_handler_fstream.cc       \
   io_handler_hdf5.cc          \
   matrix.cc                   \
   mcobs_info.cc               \
   momenta.cc                  \
   operator_info.cc            \
   scalar_defs.cc              \
   stopwatch.cc                \
   xml_handler.cc

SIGQUERY_SRCS =                \
   basic_laph_operator_info.cc \
   byte_handler.cc             \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   gen_irrep_operator_info.cc  \
   io_handler_fstream.cc       \
   io_handler_hdf5.cc          \
   matrix.cc                   \
   mcobs_info.cc               \
   momenta.cc                  \
   operator_info.cc            \
   scalar_defs.cc              \
   sigmond_query.cc            \
   xml_handler.cc

SIGCONV_SRCS =                 \
   basic_laph_operator_info.cc \
   bins_info.cc                \
   byte_handler.cc             \
   correlator_info.cc          \
   encoder.cc                  \
   ensemble_info.cc            \
   gen_irrep_operator_info.cc  \
   io_handler_fstream.cc       \
   io_handler_hdf5.cc          \
   mcobs_info.cc               \
   operator_info.cc            \
   scalar_defs.cc              \
   sigmond_convert.cc          \
   xml_handler.cc

SIGLIB_SRCS = \
   basic_laph_operator_info.cc \
   bins_info.cc                \
   bootstrapper.cc             \
   byte_handler.cc             \
   chisq_base.cc               \
   corr_data_handler.cc        \
   correlator_info.cc          \
   correlator_matrix_info.cc   \
   encoder.cc                  \
   ensemble_info.cc            \
   filelist_info.cc            \
   gen_irrep_operator_info.cc  \
   io_handler_fstream.cc       \
   io_handler_hdf5.cc          \
   lmder.cc                    \
   matrix.cc                   \
   mc_estimate.cc              \
   mcobs_handler.cc            \
   mcobs_info.cc               \
   minimizer.cc                \
   minpack.cc                  \
   momenta.cc                  \
   nl2sol.cc                   \
   obs_get_handler.cc          \
   operator_info.cc            \
   sampling_info.cc            \
   scalar_defs.cc              \
   stopwatch.cc                \
   task_utils.cc               \
   vev_data_handler.cc         \
   xml_handler.cc

