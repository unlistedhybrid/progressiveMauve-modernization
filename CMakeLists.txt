cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ------------------ MSVC CONFIG ------------------
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
    add_compile_options($<$<CONFIG:Release>:/MT>)
    add_compile_options($<$<CONFIG:Debug>:/MTd>)
    add_compile_options(/wd4834) # C++17 attribute ignored
    add_compile_options(/wd4996) # Deprecation warnings
    add_compile_options(/bigobj) # Essential for large template projects
endif()

# ------------------ DEPENDENCIES ------------------
set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)

find_package(Boost 1.89 REQUIRED COMPONENTS filesystem iostreams program_options regex)
find_package(OpenMP)
find_package(ICU COMPONENTS uc data i18n) 

# Create an interface for Boost to share settings with sub-projects
add_library(BoostInterface INTERFACE)
target_link_libraries(BoostInterface INTERFACE
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
)
if(ICU_FOUND)
    target_link_libraries(BoostInterface INTERFACE ICU::uc ICU::data ICU::i18n)
endif()

# ------------------ SUBDIRECTORIES ------------------
# These match the folders in your filelist.txt
add_subdirectory(libGenome)
add_subdirectory(muscle)
add_subdirectory(libMems) 

# ------------------ EXECUTABLE ------------------
# Define sources based on mauveAligner/src content 
set(MAUVE_SRCS
    mauveAligner/src/progressiveMauve.cpp
    mauveAligner/src/UniqueMatchFinder.cpp
)

if(WIN32)
    list(APPEND MAUVE_SRCS mauveAligner/src/win_getopt.cpp)
endif()

add_executable(progressiveMauve ${MAUVE_SRCS})

# Include paths for the executable
target_include_directories(progressiveMauve PRIVATE
    "${CMAKE_SOURCE_DIR}/mauveAligner/src"
)

# Link everything together
target_link_libraries(progressiveMauve
    PRIVATE
    libMems
    libGenome
    libMUSCLE
    BoostInterface
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(progressiveMauve PRIVATE OpenMP::OpenMP_CXX)
endif()
