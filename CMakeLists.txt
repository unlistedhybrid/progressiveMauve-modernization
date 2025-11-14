cmake_minimum_required(VERSION 3.16)
project(ProgressiveMauveModern LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Boost 1.83 REQUIRED COMPONENTS
    filesystem
    iostreams
    program_options
    regex
)

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(libGenome)
add_subdirectory(libMems)
add_subdirectory(muscle)
add_subdirectory(sgEvolver)
add_subdirectory(mauveAligner)
