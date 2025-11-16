cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Require Boost 1.89+ on all platforms for consistency
find_package(Boost 1.89 REQUIRED COMPONENTS
    filesystem
    iostreams
    program_options
    regex
)

message(STATUS "Boost version: ${Boost_VERSION}")
message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")

# Create a wrapper target for Boost that suppresses its deprecation warnings
add_library(boost_interface INTERFACE)
target_include_directories(boost_interface SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(boost_interface INTERFACE ${Boost_LIBRARIES})

# System headers don't trigger warnings
if(NOT MSVC)
    target_compile_options(boost_interface INTERFACE -Wno-deprecated-declarations)
endif()

# Global include dirs so old code builds without modification
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/muscle/libMUSCLE
    ${PROJECT_SOURCE_DIR}/sgEvolver
    ${PROJECT_SOURCE_DIR}/sgEvolver/src
    ${PROJECT_SOURCE_DIR}/mauveAligner
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)

# Platform-specific compiler flags (for your own code)
if(MSVC)
    add_compile_options(/W3 /permissive- /std:c++17)
else()
    add_compile_options(-Wall -Wextra -pedantic -Wdeprecated-declarations)
endif()

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)
add_subdirectory(mauveAligner)
