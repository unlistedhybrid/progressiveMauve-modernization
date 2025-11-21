cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ------------------ MSVC CONFIG ------------------
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_options($<$<CONFIG:Release>:/MT>)
    add_compile_options($<$<CONFIG:Debug>:/MTd>)
    add_compile_options(/wd4834)
    add_compile_options(/wd4996)
endif()

# ------------------ PLATFORM-SPECIFIC EXTRAS ------------------
set(EXTRA_SRCS "")
set(EXTRA_LIBS "")

# Windows needs getopt compatibility
if(WIN32)
    list(APPEND EXTRA_SRCS mauveAligner/src/win_getopt.cpp)
    list(APPEND EXTRA_LIBS winmm psapi)
endif()

# macOS ICU fixes
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/icu4c" "/usr/local/opt/icu4c")
    include_directories("/opt/homebrew/opt/icu4c/include" "/usr/local/opt/icu4c/include")
    link_directories("/opt/homebrew/opt/icu4c/lib" "/usr/local/opt/icu4c/lib")
endif()

# ------------------ BOOST + ZLIB + BZIP2 + ZSTD ------------------
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

# ICU (needed by Boost::regex)
find_package(ICU REQUIRED COMPONENTS uc data i18n)

# ZSTD detection
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

# Create interface for all Boost-related linking
add_library(BoostInterface INTERFACE)
target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
    ICU::uc ICU::data ICU::i18n
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

if(TARGET Boost::headers)
    target_link_libraries(BoostInterface INTERFACE Boost::headers)
elseif(TARGET Boost::boost)
    target_link_libraries(BoostInterface INTERFACE Boost::boost)
endif()

if(Boost_INCLUDE_DIRS)
    target_include_directories(BoostInterface INTERFACE ${Boost_INCLUDE_DIRS})
endif()

# ------------------ SUBDIRECTORIES ------------------
add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)

# ------------------ progressiveMauve REQUIRED SOURCES ------------------
set(PM_SRC
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/getopt.c
    mauveAligner/src/getopt1.c
    ${EXTRA_SRCS}  # win_getopt.cpp on Windows
)

# ------------------ EXECUTABLE: progressiveMauve ONLY ------------------
add_executable(progressiveMauve ${PM_SRC})

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
    ${EXTRA_LIBS}
)

# ------------------ OPTIONAL: OPENMP ACCELERATION ------------------
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(progressiveMauve PRIVATE OpenMP::OpenMP_CXX)
endif()
