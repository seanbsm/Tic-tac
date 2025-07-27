# Create a static library from the GSL helper source files
add_library(gsl_helpers STATIC coupling_coefficients.cpp auxiliary.cpp)

# Find GSL include and library paths (custom install location under $HOME/local)
find_path(GSL_INCLUDE_DIR gsl/gsl_math.h PATHS "$ENV{HOME}/local/gsl/include")
find_library(GSL_LIBRARY NAMES gsl PATHS "$ENV{HOME}/local/gsl/lib")

# Add include directories:
# - GSL_INCLUDE_DIR for GSL headers
# - src/ for internal project headers like utils/constants.h
target_include_directories(gsl_helpers
    PRIVATE
        ${GSL_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/src
)

# Link against the GSL library
target_link_libraries(gsl_helpers
    PRIVATE
        ${GSL_LIBRARY}
)

# Link OpenMP flags and libraries via the interface target
target_link_libraries(gsl_helpers PRIVATE OpenMPHelper)