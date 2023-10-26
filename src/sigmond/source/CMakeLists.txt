project(sigmond) #CMAKE_CXX_COMPILER set this using python
cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})


find_package(GRACE REQUIRED)
find_package(Minuit2 REQUIRED)
find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED HINTS "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")

message("${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pybind11/include")
message("${Python_INCLUDE_DIRS}")

include_directories(${GRACE_INCLUDES})
include_directories(${MINUIT2_INCLUDES})
include_directories("${Python_INCLUDE_DIRS}")

#check this as a directory and also include other directories.
include_directories("${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pybind11/include")

#set cxx compiler
add_compile_definitions(DOUBLEPRECISION COMPLEXNUMBERS DEFAULT_HDF5)
set(CMAKE_CXX_FLAGS "-DDEFAULTENSFILE='\"\"' -Wall -std=c++17 -O3 -lgrace_np -llapack -lMinuit2 -lhdf5 -lgomp ")

# set the project name
project(Sigmond VERSION 1.0)

add_subdirectory(analysis)
add_subdirectory(data_handling)
add_subdirectory(fitting)
add_subdirectory(observables)
add_subdirectory(plotting)
add_subdirectory(tasks)
add_subdirectory(pysigmond)

# add_subdirectory(pybind11)
# pybind11_add_module(sigmond SHARED pysigmond/pysigmond.cc)

target_include_directories(data_handling_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(data_handling_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(data_handling_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(data_handling_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/data_handling")
target_include_directories(analysis_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(analysis_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/data_handling")
target_include_directories(analysis_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(analysis_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/data_handling")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fitting")
target_include_directories(tasks_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/plotting")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/data_handling")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/plotting")
target_include_directories(fitting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fitting")
target_include_directories(observables_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(observables_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fitting")
target_include_directories(observables_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(observables_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/plotting")
target_include_directories(observables_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(plotting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(plotting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(plotting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(plotting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fitting")
target_include_directories(plotting_libs PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/plotting")


link_libraries(tasks_libs analysis_libs)
link_libraries(plotting_libs ${GRACE_LIBRARIES})
link_libraries(fitting_libs ${MINUIT2_LIBRARIES})

target_include_directories(sigmond PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/tasks")
target_include_directories(sigmond PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/data_handling")
target_include_directories(sigmond PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/analysis")
target_include_directories(sigmond PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/observables")
target_include_directories(sigmond PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fitting")
link_libraries(sigmond data_handling_libs)
link_libraries(sigmond analysis_libs)
link_libraries(sigmond tasks_libs)
link_libraries(sigmond fitting_libs)
link_libraries(sigmond observables_libs)
link_libraries(sigmond plotting_libs)


target_compile_definitions(sigmond PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

