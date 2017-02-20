
include ../Makefile.inc

CXXFLAGS= $(COMMONFLAGS) -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks               \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/fitting \
          -I$(SRC_DIR)/batch -I$(SRC_DIR)/observables

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/fitting $(SRC_DIR)/observables \
         $(SRC_DIR)/batch 

BATCH_SRCS =  \
   sigmond_batch.cc   \
   user_interface.cc

BATCH_INCS =  \
   user_interface.h


SRCS = $(TASKS_SRCS) $(ANAL_SRCS) $(DATA_SRCS) $(PLOT_SRCS) $(FIT_SRCS) $(OBS_SRCS)  $(BATCH_SRCS)

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS) $(OBS_INCS)  $(BATCH_INCS)

OBJS = $(SRCS:.cc=.o)

QOBJS = $(SIGQUERY_SRCS:.cc=.o)

SIGIOLIBOBJS = $(SIGIOLIB_SRCS:.cc=.o) 

sigmond_batch: $(OBJS) $(INCS) ../Makefile.inc
	$(CXX) -o $@ $(OBJS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

sigmond_query: $(QOBJS) $(INCS) ../Makefile.inc
	$(CXX) -o $@ $(QOBJS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -c $< -o $@

lib: $(SIGIOLIBOBJS)
	mkdir -p ../../lib
	sed "s:fill_in_template:$(SRC_DIR):g" lib.template > ../../lib/sigmondiolib.inc
	ar rcs ../../lib/libsigmondio.a $(SIGIOLIBOBJS)

clean:
	rm -f *.o; rm -f sigmond_batch; rm -f ../../lib/*
