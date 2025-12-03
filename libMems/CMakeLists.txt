# ============================================================
# libMems/CMakeLists.txt
# ============================================================

find_package(OpenMP REQUIRED)

file(GLOB LIBMEMS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/HomologyHMM/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/dmSML/*.c"
)

add_library(libMems STATIC ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/libGenome"
    "${PROJECT_SOURCE_DIR}/muscle"
)

# 1. Link the library
target_link_libraries(libMems PUBLIC
    libGenome
    libMUSCLE
    BoostInterface
    OpenMP::OpenMP_CXX
)

# 2. FORCE the Compile Options (Crucial for ParallelMemHash.cpp)
if(OpenMP_CXX_FOUND)
    target_compile_options(libMems PRIVATE ${OpenMP_CXX_FLAGS})
endif()