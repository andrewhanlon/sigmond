
set(CMAKE_CXX_FLAGS "-DDEFAULTENSFILE='\"\"' -Wall -std=c++17 -O3 -lgrace_np -llapack -lMinuit2 -lhdf5 -lgomp ") # `python -m pybind11 --includes`
# add_compile_options( "`python -m pybind11 --includes`") #libsigmond BEFORE INTERFACE 
add_library(sigmond SHARED pysigmond.cc)
# target_compile_options(sigmond INTERFACE "`python -m pybind11 --includes`")

# set(CMAKE_CXX_FLAGS "-DDEFAULTENSFILE='\"\"' -Wall -std=c++17 -O3 -lgrace_np -llapack -lMinuit2 -lhdf5 -lgomp ")
