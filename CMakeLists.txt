cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ------------------------------------------------------------
# STATIC RUNTIME LINKING (WINDOWS FIX)
# ------------------------------------------------------------
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_options($<$<CONFIG:Release>:/MT>)
    add_compile_options($<$<CONFIG:Debug>:/MTd>)
    add_compile_options(/wd4834)
    add_compile_options(/wd4996)
endif()

# ------------------------------------------------------------
# WINDOWS COMPATIBILITY SETUP
# ------------------------------------------------------------
set(EXTRA_SRCS "")
set(EXTRA_LIBS "")
if(WIN32)
    list(APPEND EXTRA_SRCS mauveAligner/src/win_getopt.cpp)
    list(APPEND EXTRA_LIBS winmm psapi)
endif()

# ------------------------------------------------------------
# BOOST & COMPRESSION SETUP
# ------------------------------------------------------------
set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)

# macOS Hints
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew" "/usr/local")
endif()

find_package(Boost 1.89 REQUIRED 
    COMPONENTS
    filesystem
    iostreams
    program_options
    regex
)

find_package(ZLIB)
find_package(BZip2)
find_package(PkgConfig)

# ZSTD Logic (Required for Boost::iostreams)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(ZSTD libzstd)
endif()
if(NOT ZSTD_FOUND)
    find_library(ZSTD_LIB 
        NAMES zstd zstd_static libzstd
        HINTS "/opt/homebrew/lib" "/usr/local/lib"
    )
else()
    set(ZSTD_LIB ${ZSTD_LIBRARIES})
    link_directories(${ZSTD_LIBRARY_DIRS}) 
endif()

add_library(BoostInterface INTERFACE)
target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
)

# Link compression libs if found
if(ZLIB_FOUND)
    target_link_libraries(BoostInterface INTERFACE ZLIB::ZLIB)
endif()
if(BZIP2_FOUND)
    target_link_libraries(BoostInterface INTERFACE BZip2::BZip2)
endif()
if(ZSTD_LIB)
    target_link_libraries(BoostInterface INTERFACE ${ZSTD_LIB})
else()
    # Fallback if specific library file wasn't found
    target_link_libraries(BoostInterface INTERFACE zstd)
endif()

if(TARGET Boost::headers)
    target_link_libraries(BoostInterface INTERFACE Boost::headers)
elseif(TARGET Boost::boost)
    target_link_libraries(BoostInterface INTERFACE Boost::boost)
endif()

if(Boost_INCLUDE_DIRS)
    target_include_directories(BoostInterface INTERFACE ${Boost_INCLUDE_DIRS})
endif()

# ------------------------------------------------------------
# SUBDIRECTORIES
# ------------------------------------------------------------
add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)

# ------------------------------------------------------------
# EXECUTABLES
# ------------------------------------------------------------
add_executable(mauveAligner
    mauveAligner/src/mauveAligner.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/alignmentProjector.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
    ${EXTRA_SRCS}
)
target_include_directories(mauveAligner PRIVATE
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)
target_link_libraries(mauveAligner libGenome libMems libMUSCLE BoostInterface ${EXTRA_LIBS})

add_executable(progressiveMauve
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/repeatoire.cpp
    mauveAligner/src/addUnalignedIntervals.cpp
    mauveAligner/src/backbone_global_to_local.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/alignmentProjector.cpp
    mauveAligner/src/RepeatHashCat.cpp
    ${EXTRA_SRCS}
)
target_include_directories(progressiveMauve PRIVATE
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)
target_link_libraries(progressiveMauve libGenome libMems libMUSCLE BoostInterface ${EXTRA_LIBS})
