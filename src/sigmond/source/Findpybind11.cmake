if (PYBIND11_INCLUDES)
  # Already in cache, be silent
  set (PYBIND11_FIND_QUIETLY TRUE)
endif (PYBIND11_INCLUDES) #do same for MINuit2

find_path (PYBIND11_INCLUDES site-packages/pybind11/include/pybind11/pybind11.h)

# find_library (PYBIND11_LIBRARIES NAMES libgrace_np.a)

# handle the QUIETLY and REQUIRED arguments and set GRACE_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (pybind11 DEFAULT_MSG PYBIND11_INCLUDES) #PYBIND11_LIBRARIES

mark_as_advanced ( PYBIND11_INCLUDES) #PYBIND11_LIBRARIES