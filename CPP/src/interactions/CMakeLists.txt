
# List all source files explicitly
set(INTERACTIONS_SOURCES
    chiral_Idaho_N3LO.cpp
    chiral_LO_internal.cpp
    chiral_LO_simple.cpp
    chiral_N2LOopt.cpp
    chiral_twobody.cpp
    malfliet_tjon.cpp
    nijmegen.cpp
    potential_model.cpp
    OPEP.cpp
    chp/chiral-twobody-potentials.f90
    chp/chp-set.f90
    nijmegen/nijmegen_interface.f90
    nijmegen/pnijm.f
)

# Create object library
add_library(interactions OBJECT ${INTERACTIONS_SOURCES})

# Set compile options
target_compile_options(interactions PRIVATE -O3 -fPIC)
target_compile_options(interactions PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000>)

# Set fixed-form for specific Fortran files
set_source_files_properties(nijmegen/pnijm.f PROPERTIES Fortran_FORMAT FIXED)

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
