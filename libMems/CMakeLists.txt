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

target_link_libraries(libMems PUBLIC
    libGenome
    libMUSCLE
    BoostInterface
    OpenMP::OpenMP_CXX
)

# CRITICAL FIX: Explicitly force the compilation flags so _OPENMP is defined
if(OpenMP_CXX_FOUND)
    if(MSVC)
        target_compile_options(libMems PRIVATE /openmp)
    elseif(APPLE)
        target_compile_options(libMems PRIVATE -Xpreprocessor -fopenmp)
    else()
        target_compile_options(libMems PRIVATE -fopenmp)
    endif()
endif()