cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_options($<$<CONFIG:Release>:/MT>)
    add_compile_options($<$<CONFIG:Debug>:/MTd>)
    add_compile_options(/wd4834)
    add_compile_options(/wd4996)
endif()

# --- WINDOWS SETUP ---
set(EXTRA_SRCS "")
set(EXTRA_LIBS "")
if(WIN32)
    list(APPEND EXTRA_SRCS mauveAligner/src/win_getopt.cpp)
    list(APPEND EXTRA_LIBS winmm psapi)
endif()

# --- MACOS HOMEBREW HINTS ---
if(APPLE)
    # icu4c is keg-only on macOS (not linked to /usr/local), so we must hint explicitly
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/icu4c" "/usr/local/opt/icu4c")
endif()

# --- BOOST & DEPS SETUP ---
set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)

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

# --- NEW: Find ICU (Required for Boost::regex) ---
find_package(ICU COMPONENTS uc data i18n)

# ZSTD Logic
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

if(ZLIB_FOUND)
    target_link_libraries(BoostInterface INTERFACE ZLIB::ZLIB)
endif()
if(BZIP2_FOUND)
    target_link_libraries(BoostInterface INTERFACE BZip2::BZip2)
endif()
if(ZSTD_LIB)
    target_link_libraries(BoostInterface INTERFACE ${ZSTD_LIB})
else()
    target_link_libraries(BoostInterface INTERFACE zstd)
endif()

# --- NEW: Link ICU if found ---
if(ICU_FOUND)
    target_link_libraries(BoostInterface INTERFACE ICU::uc ICU::data ICU::i18n)
endif()

if(TARGET Boost::headers)
    target_link_libraries(BoostInterface INTERFACE Boost::headers)
elseif(TARGET Boost::boost)
    target_link_libraries(BoostInterface INTERFACE Boost::boost)
endif()

if(Boost_INCLUDE_DIRS)
    target_include_directories(BoostInterface INTERFACE ${Boost_INCLUDE_DIRS})
endif()

# --- SUBDIRECTORIES ---
add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)

# --- EXECUTABLES ---
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
