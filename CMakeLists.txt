cmake_minimum_required(VERSION 3.20)
project(ProgressiveMauveModern LANGUAGES C CXX)

# ============================================================
# 1. GLOBAL BUILD SETTINGS
# ============================================================

# Modern Boost (1.89) requires modern C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ACCURACY CONFIGURATION:
# -O3: High optimization.
# -fno-fast-math: STRICT math. Prevents reordering. Ensures 1859 matches result.
# -DCOMMAND_LINE: Logic flag.
add_compile_options(-O3 -fno-fast-math -DCOMMAND_LINE -fno-strict-aliasing)

# Link standard libraries statically
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")

# ============================================================
# 2. OPTIONS
# ============================================================

option(ENABLE_SANITIZERS "Enable ASAN + UBSAN" OFF)
option(ENABLE_CCACHE "Enable CCache" OFF)
option(ENABLE_OPENMP "Enable OpenMP" ON)

# Disable ZLIB to match original logic
option(ENABLE_ZLIB "Enable ZLIB compression" OFF) 

# ============================================================
# 3. DEPENDENCIES (SYSTEM BOOST 1.89)
# ============================================================

# We use the System Boost (1.89)
set(Boost_USE_STATIC_LIBS ON)

# Allow CMake to search system paths (/usr/lib, etc.)
set(Boost_NO_SYSTEM_PATHS OFF)

# Find version 1.89 or higher
find_package(Boost 1.89 REQUIRED COMPONENTS filesystem iostreams program_options regex system)

if(ENABLE_ZLIB)
    find_package(ZLIB REQUIRED)
    find_package(BZip2 REQUIRED)
    add_definitions(-DHAVE_LIBZ)
endif()

find_package(ICU REQUIRED COMPONENTS uc data i18n)

# ============================================================
# 4. BOOST INTERFACE
# ============================================================

add_library(BoostInterface INTERFACE)
target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
    Boost::system
    ICU::uc ICU::data ICU::i18n
)

if(ENABLE_ZLIB)
    target_link_libraries(BoostInterface INTERFACE ZLIB::ZLIB BZip2::BZip2)
endif()

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(BoostInterface INTERFACE OpenMP::OpenMP_CXX)
endif()

# ============================================================
# 5. SUBDIRECTORIES
# ============================================================

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)

# ============================================================
# 6. EXECUTABLE
# ============================================================

set(PM_SRC
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
)

if(WIN32)
    list(APPEND PM_SRC mauveAligner/src/win_getopt.cpp)
endif()

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
)