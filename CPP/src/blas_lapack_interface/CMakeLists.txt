# Create a static library from the source file
add_library(blas_lapack_interface STATIC blas_lapack_interface.cpp)

# Add include directories:
# - LAPACKE_INCLUDE_DIR is needed for LAPACK-related headers
# - src/ is added as a private include path for internal headers like utils/constants.h
target_include_directories(blas_lapack_interface
    PUBLIC
        ${LAPACKE_INCLUDE_DIR}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link against required libraries:
# - BLAS, LAPACKE, and OpenBLAS for numerical routines
target_link_libraries(blas_lapack_interface
    PRIVATE
        ${BLAS_LIBRARIES}
        ${LAPACKE_LIBRARY}
        ${OPENBLAS_LIB}
)

# Link OpenMP flags and libraries via the interface target
target_link_libraries(blas_lapack_interface PRIVATE OpenMPHelper)