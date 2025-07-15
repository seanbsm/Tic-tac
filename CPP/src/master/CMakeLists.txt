# Recursively collect all .cpp source files in this directory
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS *.cpp)

# Create an object library from the collected sources
add_library(master OBJECT ${SOURCES})

# Add private include path for internal headers (e.g., from utils/, core/, etc.)
target_include_directories(master
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
