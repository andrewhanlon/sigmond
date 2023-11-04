project(sigmond) #CMAKE_CXX_COMPILER set this using python
cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})

# find_package(GRACE REQUIRED)
# find_package(Minuit2 REQUIRED)
find_package(Python COMPONENTS Interpreter Development)

#check this as a directory and also include other directories.
find_package(pybind11 REQUIRED HINTS "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")

# include_directories("${GRACE_INCLUDES}")
# include_directories("${MINUIT2_INCLUDES}")
include_directories("${Python_INCLUDE_DIRS}")


include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/tasks")
include_directories("${PROJECT_SOURCE_DIR}/observables")
include_directories("${PROJECT_SOURCE_DIR}/analysis")
include_directories("${PROJECT_SOURCE_DIR}/data_handling")
include_directories("${PROJECT_SOURCE_DIR}/fitting")
include_directories("${PROJECT_SOURCE_DIR}/plotting")
include_directories("${PROJECT_SOURCE_DIR}/pysigmond")

#set cxx compiler
add_compile_definitions(DOUBLEPRECISION COMPLEXNUMBERS DEFAULT_HDF5 NOGRACE NO_MINUIT)
set(CMAKE_CXX_FLAGS "-DDEFAULTENSFILE='\"\"' -Wall -std=c++17 -O3 -llapack -lhdf5 -lgomp")

# set the project name
project(Sigmond VERSION 1.0)

add_subdirectory(analysis)
add_subdirectory(data_handling)
add_subdirectory(fitting)
add_subdirectory(observables)
add_subdirectory(plotting)
add_subdirectory(tasks)
# add_subdirectory(pysigmond)

# add_subdirectory(pybind11)
pybind11_add_module(sigmond pysigmond/pysigmond.cc)
# target_include_directories(sigmond PUBLIC "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pybind11/include")

target_link_libraries(sigmond PUBLIC data_handling_libs )
target_link_libraries(sigmond PUBLIC analysis_libs )
target_link_libraries(sigmond PUBLIC tasks_libs)
target_link_libraries(sigmond PUBLIC fitting_libs)
target_link_libraries(sigmond PUBLIC observables_libs)
target_link_libraries(sigmond PUBLIC plotting_libs)


target_compile_definitions(sigmond PRIVATE VERSION_INFO=${VERSION_INFO})