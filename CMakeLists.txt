cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Boost 1.89 REQUIRED COMPONENTS
    system
    filesystem
    iostreams
    program_options
    regex
)

message(STATUS "Boost version: ${Boost_VERSION}")
message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library dir: ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost 1.89 not found. Install with:")
    message(FATAL_ERROR "  macOS: brew install boost@1.89 or brew install boost")
    message(FATAL_ERROR "  Linux: sudo apt-get install libboost-all-dev")
endif()

add_library(BoostInterface INTERFACE IMPORTED)
set_target_properties(BoostInterface PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
    INTERFACE_LINK_DIRECTORIES "${Boost_LIBRARY_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Boost_LIBRARIES}"
)

if(NOT MSVC)
    target_compile_options(BoostInterface INTERFACE
        -Wno-deprecated-declarations
        -Wno-unused-parameter
    )
endif()

if(MSVC)
    add_compile_options(/W3 /permissive- /std:c++17)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(
        -Wall
        -Wextra
        -pedantic
    )
endif()

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
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

link_directories(${Boost_LIBRARY_DIRS})

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)
add_subdirectory(mauveAligner)
