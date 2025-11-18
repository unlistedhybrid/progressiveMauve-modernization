cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

find_package(Boost 1.89 REQUIRED 
    COMPONENTS
    filesystem
    iostreams
    program_options
    regex
)

include_directories(${Boost_INCLUDE_DIRS})

add_library(BoostInterface INTERFACE IMPORTED)
set_target_properties(BoostInterface PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Boost_LIBRARIES}"
)

if(NOT MSVC)
    target_compile_options(BoostInterface INTERFACE
        -Wno-deprecated-declarations
        -Wno-unused-parameter
    )
endif()

# ------------------------------------------------------------
# SHARED LIBRARIES
# ------------------------------------------------------------

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)

# sgEvolver is optional tools
add_subdirectory(sgEvolver)

# ------------------------------------------------------------
# EXECUTABLE: mauveAligner
# ------------------------------------------------------------
add_executable(mauveAligner
    mauveAligner/src/mauveAligner.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/alignmentProjector.cpp
)

target_include_directories(mauveAligner PRIVATE
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)

target_link_libraries(mauveAligner
    libGenome
    libMems
    libMUSCLE
    BoostInterface
)

# ------------------------------------------------------------
# EXECUTABLE: progressiveMauve
# ------------------------------------------------------------
add_executable(progressiveMauve
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/repeatoire.cpp
    mauveAligner/src/addUnalignedIntervals.cpp
    mauveAligner/src/backbone_global_to_local.cpp
)

target_include_directories(progressiveMauve PRIVATE
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)

target_link_libraries(progressiveMauve
    libGenome
    libMems
    libMUSCLE
    BoostInterface
)
