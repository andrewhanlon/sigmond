
include ../Makefile.inc

CXXFLAGS= $(COMMONFLAGS) -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/fitting \
          -I$(SRC_DIR)/batch -I$(SRC_DIR)/observables

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/fitting $(SRC_DIR)/observables \
         $(SRC_DIR)/batch $(SRC_DIR)/python_binding

BATCH_INCS =  \
   user_interface.h

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS) $(OBS_INCS)  $(BATCH_INCS)

SRCS = $(SIGLIB_SRCS)

OBJS = $(SRCS:.cc=.o)

sigmond_bind: sigmond_bind.cc lib $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -fPIC -shared `python -m pybind11 --includes` -L. -Wl,-R,`pwd` -lsigmond $< -o sigmondbind`python-config --extension-suffix`

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

lib: $(OBJS) $(INCS) ../Makefile.inc
	$(CXX) -o libsigmond.so $(CXXFLAGS) -Wl,-z,defs -shared $(OBJS)

install:
	cp sigmondbind`python-config --extension-suffix` `python -m site --user-site`

clean:
	rm -f *.o; rm -f libsigmond.so; rm -f sigmondbind.*.so
