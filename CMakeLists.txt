cmake_minimum_required(VERSION 3.20)

# Suppress the "FindBoost is removed" warning (CMP0167)
cmake_policy(SET CMP0167 NEW)

project(ProgressiveMauveModern LANGUAGES C CXX)

# ============================================================
# 1. GLOBAL BUILD SETTINGS
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --------------------------------------------------------
# STATIC LINKING & RUNTIME
# --------------------------------------------------------

# Force CMake to look for static Boost libraries (.a / .lib)
set(Boost_USE_STATIC_LIBS ON)

# Windows-Specific: Force Static Runtime (/MT)
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
    add_compile_options(/O2 /fp:precise /DCOMMAND_LINE)
else()
    # GCC / Clang Flags
    add_compile_options(-O3 -fno-fast-math -DCOMMAND_LINE -fno-strict-aliasing)
    
    # Linux: Link system libraries statically (Apple ignores this)
    if(NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    endif()
endif()

# ============================================================
# 2. OPTIONS
# ============================================================

option(ENABLE_SANITIZERS "Enable ASAN + UBSAN" OFF)
option(ENABLE_OPENMP "Enable OpenMP" ON)
option(ENABLE_ZLIB "Enable ZLIB compression" OFF) 

if(ENABLE_SANITIZERS AND NOT MSVC)
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# ============================================================
# 3. DEPENDENCIES
# ============================================================

# Use custom paths from build.yml (CMAKE_PREFIX_PATH)
set(Boost_NO_SYSTEM_PATHS OFF)

# REMOVED 'system' to fix the build error
find_package(Boost 1.89 REQUIRED COMPONENTS filesystem iostreams program_options regex)

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
    # REMOVED 'Boost::system' here as well
    ICU::uc ICU::data ICU::i18n
)

if(ENABLE_ZLIB)
    target_link_libraries(BoostInterface INTERFACE ZLIB::ZLIB BZip2::BZip2)
endif()

if(ENABLE_OPENMP)
    # This works on Mac/Linux/Windows because build.yml sets the paths
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
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mauveAligner/src/win_getopt.cpp")
        list(APPEND PM_SRC mauveAligner/src/win_getopt.cpp)
    endif()
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