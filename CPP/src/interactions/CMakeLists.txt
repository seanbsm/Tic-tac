
# Collect all source files
file(GLOB_RECURSE INTERACTIONS_SOURCES CONFIGURE_DEPENDS *.cpp *.f90 *.f)

# Create object library
add_library(interactions OBJECT ${INTERACTIONS_SOURCES})

# Set compile flags for all sources
set_source_files_properties(${INTERACTIONS_SOURCES} PROPERTIES COMPILE_FLAGS "-O3 -fPIC")

# Special compile flags for Fortran sources
foreach(src ${INTERACTIONS_SOURCES})
    if(src MATCHES "\\.f90$")
        set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000")
    endif()
endforeach()

# Add private include path for internal headers
target_include_directories(interactions
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link HDF5 for includes
if(TARGET HDF5::HDF5)
  target_link_libraries(interactions PRIVATE HDF5::HDF5 ${HDF5_HL_TARGET})
elseif(TARGET hdf5::hdf5_cpp)
  target_link_libraries(interactions PRIVATE hdf5::hdf5_cpp ${HDF5_HL_TARGET})
else()
  target_link_libraries(interactions PRIVATE HDF5::HDF5_CXX)
endif()
