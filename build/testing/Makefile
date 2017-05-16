
include ../Makefile.inc

CXXFLAGS= $(COMMONFLAGS) -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks             \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/batch \
          -I$(SRC_DIR)/fitting -I$(SRC_DIR)/observables -I$(SRC_DIR)/testing

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/observables $(SRC_DIR)/batch   \
         $(SRC_DIR)/fitting $(SRC_DIR)/testing

TEST_SRCS = \
   test_args_handler.cc       \
   test_bins_info.cc          \
   test_bootstrap.cc          \
   test_chisq.cc              \
   test_corr_handler.cc       \
   test_correlator_info.cc    \
   test_corr_estimates.cc     \
   test_data_getput.cc        \
   test_graceplot.cc          \
   test_matrix.cc             \
   test_mcobs_info.cc         \
   test_mcobs_handler.cc      \
   test_multicompare.cc       \
   test_obs_get_handler.cc    \
   test_operator_info.cc      \
   test_reorder.cc            \
   test_rotate_corr.cc        \
   test_sampling_info.cc      \
   test_task_handler.cc       \
   test_tcorr.cc              \
   test_vev_handler.cc        \
   test_xml_handler.cc

TEST_INCS = \
   testing.h               \
   test_bootstrap.h        \
   test_mcobs_handler.h    \
   test_obs_get_handler.h  


BATCH_SRCS =  \
   user_interface.cc

BATCH_INCS =  \
   user_interface.h


SRCS = $(TASKS_SRCS) $(ANAL_SRCS) $(DATA_SRCS) $(PLOT_SRCS) $(FIT_SRCS)  $(OBS_SRCS) $(BATCH_SRCS)

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS)  $(OBS_INCS) $(BATCH_INCS)

OBJS = $(SRCS:.cc=.o) 

TESTOBJS = $(TEST_SRCS:.cc=.o) 


runtests: run_tests.o $(OBJS) $(TESTOBJS) $(INCS) ../Makefile.inc
	$(CXX) -o $@ run_tests.o $(OBJS) $(TESTOBJS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) -Wno-psabi $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o; rm -f runtests
