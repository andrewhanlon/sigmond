
SIG_SRC_DIR=../../source

SIG_CXXFLAGS= -I$(SIG_SRC_DIR)/tasks -I$(SIG_SRC_DIR)/data_handling  \
              -I$(SIG_SRC_DIR)/observables -I$(SIG_SRC_DIR)/analysis

SIG_VPATH =  $(SIG_SRC_DIR)/tasks $(SIG_SRC_DIR)/data_handling \
             $(SIG_SRC_DIR)/observables -I$(SIG_SRC_DIR)/analysis
