cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0167 NEW)

project(ProgressiveMauveModern LANGUAGES C CXX)

# ============================================================
# 1. OPENMP
# ============================================================
option(ENABLE_OPENMP "Enable OpenMP" ON)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    
    # Force flags into the global CXX flags.
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp") 
    elseif(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif()
    
    message(STATUS ">>> GLOBAL CXX FLAGS: ${CMAKE_CXX_FLAGS}")
endif()

# ============================================================
# 2. GLOBAL SETTINGS
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(Boost_USE_STATIC_LIBS ON)

option(ENABLE_SANITIZERS "Enable ASAN + UBSAN" OFF)
option(ENABLE_ZLIB "Enable ZLIB compression" OFF) 

# ============================================================
# 3. COMPILER FLAGS & RUNTIME
# ============================================================

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    set(Boost_USE_STATIC_RUNTIME ON)
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
    add_compile_options(/O2 /fp:precise /DCOMMAND_LINE)
else()
    add_compile_options(-O3 -fno-fast-math -DCOMMAND_LINE -fno-strict-aliasing)
    if(NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    endif()
endif()

# ============================================================
# 4. DEPENDENCIES
# ============================================================

set(Boost_NO_SYSTEM_PATHS OFF)
find_package(Boost 1.89 REQUIRED COMPONENTS filesystem iostreams program_options regex)

if(ENABLE_ZLIB)
    find_package(ZLIB REQUIRED)
    find_package(BZip2 REQUIRED)
    add_definitions(-DHAVE_LIBZ)
endif()

find_package(ICU REQUIRED COMPONENTS uc data i18n)

# ============================================================
# 5. BOOST INTERFACE
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

# ============================================================
# 6. SUBDIRECTORIES
# ============================================================

add_subdirectory(libGenome)
target_link_libraries(libGenome PUBLIC BoostInterface)

add_subdirectory(libMems)
target_link_libraries(libMems PUBLIC BoostInterface)

add_subdirectory(muscle)

# ============================================================
# 7. EXECUTABLE
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

if(ENABLE_OPENMP AND NOT MSVC)
    target_link_libraries(progressiveMauve PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
    # The functions timeGetTime/timeBeginPeriod are in winmm.lib
    # The function GetProcessMemoryInfo is in psapi.lib
    target_link_libraries(progressiveMauve
        PRIVATE
        winmm
        psapi
    )
endif()

target_link_libraries(progressiveMauve
    PUBLIC
    libGenome
    libMems
    libMUSCLE
    BoostInterface
)