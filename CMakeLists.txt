cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable specific MSVC warnings that are cluttering the output
if(MSVC)
    add_compile_options(/wd4834) # Disable [[nodiscard]] warning
    add_compile_options(/wd4996) # Disable deprecation warnings
endif()

# ------------------------------------------------------------
# WINDOWS COMPATIBILITY SETUP
# ------------------------------------------------------------
set(EXTRA_SRCS "")
set(EXTRA_LIBS "")

if(WIN32)
    # 1. Compile win_getopt.cpp for Windows (Fixes: unresolved external symbol getopt_long)
    # IMPORTANT: Ensure win_getopt.cpp is physically located at mauveAligner/src/win_getopt.cpp
    list(APPEND EXTRA_SRCS mauveAligner/src/win_getopt.cpp)
    
    # 2. Link System Libraries (Fixes: unresolved external symbol timeGetTime, GetProcessMemoryInfo)
    # winmm = Windows Multimedia API
    # psapi = Process Status API
    list(APPEND EXTRA_LIBS winmm psapi)
endif()

# ------------------------------------------------------------
# BOOST SETUP
# ------------------------------------------------------------
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

add_library(BoostInterface INTERFACE)

target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
)

if(TARGET Boost::headers)
    target_link_libraries(BoostInterface INTERFACE Boost::headers)
elseif(TARGET Boost::boost)
    target_link_libraries(BoostInterface INTERFACE Boost::boost)
endif()

if(Boost_INCLUDE_DIRS)
    target_include_directories(BoostInterface INTERFACE ${Boost_INCLUDE_DIRS})
endif()
if(Boost_INCLUDE_DIR)
    target_include_directories(BoostInterface INTERFACE ${Boost_INCLUDE_DIR})
endif()

# ------------------------------------------------------------
# SUBDIRECTORIES
# ------------------------------------------------------------

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)

# ------------------------------------------------------------
# EXECUTABLE: mauveAligner
# ------------------------------------------------------------
add_executable(mauveAligner
    mauveAligner/src/mauveAligner.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/alignmentProjector.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
    ${EXTRA_SRCS}  # <--- Adds win_getopt.cpp on Windows
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
    ${EXTRA_LIBS}  # <--- Links winmm.lib and psapi.lib on Windows
)

# ------------------------------------------------------------
# EXECUTABLE: progressiveMauve
# ------------------------------------------------------------
add_executable(progressiveMauve
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/repeatoire.cpp
    mauveAligner/src/addUnalignedIntervals.cpp
    mauveAligner/src/backbone_global_to_local.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
    mauveAligner/src/AlignmentTree.cpp
    mauveAligner/src/alignmentProjector.cpp
    mauveAligner/src/RepeatHashCat.cpp
    ${EXTRA_SRCS}  # <--- Adds win_getopt.cpp on Windows
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
    ${EXTRA_LIBS}  # <--- Links winmm.lib and psapi.lib on Windows
)
