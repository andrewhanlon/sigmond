
include ../Makefile.inc

GRACE_LIB_DIR=../../grace_np

CXXFLAGS= $(COMMONFLAGS) $(DCFLAGS) -fPIC -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks    \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/fitting \
          -I$(SRC_DIR)/observables

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/fitting $(SRC_DIR)/observables \
         $(SRC_DIR)/pysigmond $(GRACE_LIB_DIR)

SRCS = $(DRIVER_SRC) $(TASKS_SRCS) $(ANAL_SRCS) $(DATA_SRCS) $(PLOT_SRCS) $(FIT_SRCS) $(OBS_SRCS)

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS) $(OBS_INCS)

OBJS = $(SRCS:.cc=.o)

#LD_LIBRARY_PATH=.
CC=gcc
CFLAGS=-O2 -fno-common -Wall -Wpointer-arith -Wnested-externs -I$(GRACE_LIB_DIR) -fPIC
GRACE_LIB_OJBS = grace_np.o gracef_np.o gracef_np_.o

PYTHON=python
PYTHON_CONFIG=python-config

INSTALL_DIR=$(HOME)/.local/lib

lib: $(OBJS) $(INCS) ../Makefile.inc libgrace_np.a
	$(CXX) -o libsigmond.so $(CXXFLAGS) -shared $(OBJS) $(LDFLAGS) $(LIBS) -L.

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -c $< -o $@

gracef_np_.o: gracef_np.c
	$(CC) $(CFLAGS) -DNEED_F77_UNDERSCORE -c -o $@ $(GRACE_LIB_DIR)/gracef_np.c

$(GRACE_LIB_OJBS): config.h grace_np.h

libgrace_np.a: $(GRACE_LIB_OJBS)
	rm -f $@
	ar cr $@ $(GRACE_LIB_OJBS)
	ranlib $@

install: pysigmond.cc
	mkdir -p $(INSTALL_DIR)
	cp libsigmond.so $(INSTALL_DIR)/
	$(CXX) $(CXXFLAGS) -shared `$(PYTHON) -m pybind11 --includes` -L$(INSTALL_DIR) -Wl,-R,$(INSTALL_DIR) -lsigmond $< -o sigmond`$(PYTHON_CONFIG) --extension-suffix`
	mkdir -p `$(PYTHON) -m site --user-site`
	cp sigmond`$(PYTHON_CONFIG) --extension-suffix` `$(PYTHON) -m site --user-site`

uninstall:
	rm -f $(INSTALL_DIR)/libsigmond.so
	rm -f `$(PYTHON) -m site --user-site`/sigmond`$(PYTHON_CONFIG) --extension-suffix`

clean:
	rm -f *.o; rm -f libsigmond.so; rm -f sigmond.*.so; rm -f libgrace_np.a
