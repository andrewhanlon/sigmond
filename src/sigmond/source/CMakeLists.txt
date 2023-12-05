set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
project(sigmond) #CMAKE_CXX_COMPILER set this using python
cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})

find_package(Python COMPONENTS Interpreter Development)

#check this as a directory and also include other directories.
set(WINDOWS_LIB_DIRS "${Python_LIBRARY_DIRS}")
string(REPLACE "miniconda3/libs" "miniconda3/Lib" WINDOWS_LIB_DIRS "${WINDOWS_LIB_DIRS}")

message("${Python_LIBRARY_DIRS}")
find_package(pybind11 REQUIRED 
                      HINTS "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
                            "${WINDOWS_LIB_DIRS}/site-packages"
			    "/usr/local/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/dist-packages"
			)

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
add_compile_definitions(DOUBLEPRECISION COMPLEXNUMBERS DEFAULT_HDF5 NOGRACE NO_MINUIT NOXML HDF5 LAPACK)

# set the project name
project(Sigmond VERSION 1.0)

add_subdirectory(analysis)
add_subdirectory(data_handling)
add_subdirectory(fitting)
add_subdirectory(observables)
add_subdirectory(plotting)
add_subdirectory(tasks)

#link_libraries(tasks_libs data_handling_lib analysis_libs fitting_libs observables_libs plotting_libs)

pybind11_add_module(sigmond pysigmond/pysigmond.cc)

# target_link_libraries(sigmond PUBLIC data_handling_libs analysis_libs tasks_libs fitting_libs observables_libs plotting_libs) 
target_link_libraries(sigmond PUBLIC data_handling_libs )
target_link_libraries(sigmond PUBLIC analysis_libs )
target_link_libraries(sigmond PUBLIC tasks_libs )
target_link_libraries(sigmond PUBLIC fitting_libs )
target_link_libraries(sigmond PUBLIC observables_libs )
target_link_libraries(sigmond PUBLIC plotting_libs)


target_compile_definitions(sigmond PRIVATE VERSION_INFO=${VERSION_INFO})
