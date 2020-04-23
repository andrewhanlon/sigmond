
include ../Makefile.inc

CXXFLAGS= $(COMMONFLAGS) $(DCFLAGS) -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks    \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/fitting \
          -I$(SRC_DIR)/observables

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/fitting $(SRC_DIR)/observables \
         $(SRC_DIR)/pysigmond

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS) $(OBS_INCS)

SRCS = $(SIGLIB_SRCS)

OBJS = $(SRCS:.cc=.o)

INSTALL_DIR=$(HOME)/.local/lib

lib: $(OBJS) $(INCS) ../Makefile.inc
	$(CXX) -o libsigmond.so $(CXXFLAGS) -shared $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

install: pysigmond.cc
	cp libsigmond.so $(INSTALL_DIR)
	$(CXX) $(CXXFLAGS) -fPIC -shared `python -m pybind11 --includes` -L$(INSTALL_DIR) -Wl,-R,$(INSTALL_DIR) -lsigmond $< -o sigmond`python-config --extension-suffix`
	cp sigmond`python-config --extension-suffix` `python -m site --user-site`

uninstall:
	rm -f $(INSTALL_DIR)/libsigmond.so
	rm -f `python -m site --user-site`/sigmond`python-config --extension-suffix`

clean:
	rm -f *.o; rm -f libsigmond.so; rm -f sigmond.*.so
