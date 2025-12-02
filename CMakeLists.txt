cmake_minimum_required(VERSION 3.20)
project(ProgressiveMauveModern LANGUAGES C CXX)

# ============================================================
# GLOBAL SETTINGS
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable BoostConfig.cmake (b2 does NOT generate per-component config files)
set(Boost_NO_BOOST_CMAKE ON)

# Static libs everywhere
set(Boost_USE_STATIC_LIBS ON)

# Disable the new BoostConfig.cmake lookup (b2 does NOT create component configs)
cmake_policy(SET CMP0167 OLD)

# Force classic FindBoost mode
set(Boost_NO_BOOST_CMAKE ON)

# DO NOT search system Boost
set(Boost_NO_SYSTEM_PATHS ON)

find_package(Boost 1.89 REQUIRED COMPONENTS filesystem iostreams program_options regex system)


# ============================================================
# COMPILER FLAGS
# ============================================================

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
    add_compile_options(/O2 /fp:precise /DCOMMAND_LINE)
else()
    add_compile_options(-O3 -fno-fast-math -fno-strict-aliasing -DCOMMAND_LINE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

option(ENABLE_SANITIZERS "Enable ASAN/UBSAN" OFF)
option(ENABLE_OPENMP "Enable OpenMP" ON)
option(ENABLE_ZLIB "Enable ZLIB/BZip2 compression" OFF)

if(ENABLE_SANITIZERS AND NOT MSVC)
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# ============================================================
# DEPENDENCIES
# ============================================================

find_package(Boost 1.89 REQUIRED COMPONENTS
    filesystem iostreams program_options regex
)
find_package(ICU REQUIRED COMPONENTS uc data i18n)

if(ENABLE_ZLIB)
    find_package(ZLIB REQUIRED)
    find_package(BZip2 REQUIRED)
    add_definitions(-DHAVE_LIBZ)
endif()

# ============================================================
# BOOST + ICU INTERFACE
# ============================================================

add_library(BoostInterface INTERFACE)
target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
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
# SUBDIRECTORIES
# ============================================================

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)

# ============================================================
# EXECUTABLE
# ============================================================

set(PM_SRC
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
)

if(WIN32 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mauveAligner/src/win_getopt.cpp")
    list(APPEND PM_SRC mauveAligner/src/win_getopt.cpp)
endif()

add_executable(progressiveMauve ${PM_SRC})

target_include_directories(progressiveMauve PRIVATE
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
    ${PROJECT_SOURCE_DIR}/libMems/libMems
    ${PROJECT_SOURCE_DIR}/muscle/libMUSCLE
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)

target_link_libraries(progressiveMauve
    libGenome
    libMems
    libMUSCLE
    BoostInterface
)
