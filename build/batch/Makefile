
include ../Makefile.inc

CXXFLAGS= $(COMMONFLAGS) $(DCFLAGS) -I$(SRC_DIR)/analysis -I$(SRC_DIR)/tasks    \
          -I$(SRC_DIR)/data_handling -I$(SRC_DIR)/plotting -I$(SRC_DIR)/fitting \
          -I$(SRC_DIR)/observables

VPATH =  $(SRC_DIR)/analysis $(SRC_DIR)/tasks $(SRC_DIR)/data_handling \
         $(SRC_DIR)/plotting $(SRC_DIR)/fitting $(SRC_DIR)/observables


SRCS = $(DRIVER_SRC) $(TASKS_SRCS) $(ANAL_SRCS) $(DATA_SRCS) $(PLOT_SRCS) $(FIT_SRCS) $(OBS_SRCS)

INCS = $(TASKS_INCS) $(ANAL_INCS) $(DATA_INCS) $(PLOT_INCS) $(FIT_INCS) $(OBS_INCS)

OBJS = $(SRCS:.cc=.o)

QOBJS = $(SIGQUERY_SRCS:.cc=.o)

SIGIOLIBOBJS = $(SIGIOLIB_SRCS:.cc=.o) 

INSTALL_DIR=$(HOME)/.local/bin

sigmond: $(OBJS) $(INCS) ../Makefile.inc
	$(CXX) -o $@ $(OBJS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

sigmond_query: $(QOBJS) $(INCS) ../Makefile.inc
	$(CXX) -o $@ $(QOBJS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(INCS) ../Makefile.inc
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: sigmond sigmond_query
	mkdir -p $(INSTALL_DIR)
	cp sigmond $(INSTALL_DIR)/
	cp sigmond_query $(INSTALL_DIR)/

uninstall:
	rm -f $(INSTALL_DIR)/sigmond
	rm -f $(INSTALL_DIR)/sigmond_query

lib: $(SIGIOLIBOBJS)
	mkdir -p ../../lib
	sed "s:fill_in_template:$(SRC_DIR):g" lib.template > ../../lib/sigmondiolib.inc
	ar rcs ../../lib/libsigmondio.a $(SIGIOLIBOBJS)
	$(CXX) -shared -o ../../lib/libsigmondio.so $(SIGIOLIBOBJS)

clean:
	rm -f *.o; rm -f sigmond sigmond_query; rm -f ../../lib/*
