# Create a static library from the GSL helper source files
add_library(gsl_helpers STATIC coupling_coefficients.cpp auxiliary.cpp)

# Find GSL in a portable way: try CONFIG, then pkg-config, then manual
find_package(GSL CONFIG QUIET)
if(TARGET GSL::gsl)
  message(STATUS "Found GSL via CMake config")
endif()

if(NOT TARGET GSL::gsl)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(GSL_PC QUIET gsl)
  endif()

  if(GSL_PC_FOUND)
    add_library(GSL::gsl INTERFACE IMPORTED)
    target_include_directories(GSL::gsl INTERFACE ${GSL_PC_INCLUDE_DIRS})
    target_link_libraries(GSL::gsl INTERFACE ${GSL_PC_LIBRARIES})
    message(STATUS "Found GSL via pkg-config")
  endif()
endif()

if(NOT TARGET GSL::gsl)
  find_path(GSL_INCLUDE_DIR
    NAMES gsl/gsl_math.h
    HINTS ENV GSL_ROOT
    PATH_SUFFIXES include
  )

  find_library(GSL_LIBRARY
    NAMES gsl
    HINTS ENV GSL_ROOT
    PATH_SUFFIXES lib lib64
  )

  if(GSL_INCLUDE_DIR AND GSL_LIBRARY)
    add_library(GSL::gsl UNKNOWN IMPORTED)
    set_target_properties(GSL::gsl PROPERTIES
      IMPORTED_LOCATION "${GSL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
    )
    message(STATUS "Found GSL via manual search")
  else()
    message(FATAL_ERROR "GSL not found")
  endif()
endif()

# Add include directories:
# - src/ for internal project headers like utils/constants.h
target_include_directories(gsl_helpers
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link against the GSL library using the portable target
if(TARGET HDF5::HDF5)
  target_link_libraries(gsl_helpers PRIVATE GSL::gsl HDF5::HDF5 ${HDF5_HL_TARGET})
elseif(TARGET hdf5::hdf5_cpp)
  target_link_libraries(gsl_helpers PRIVATE GSL::gsl hdf5::hdf5_cpp ${HDF5_HL_TARGET})
else()
  target_link_libraries(gsl_helpers PRIVATE GSL::gsl HDF5::HDF5_CXX)
endif()

# Link OpenMP flags and libraries via the interface target
target_link_libraries(gsl_helpers PRIVATE OpenMPHelper)