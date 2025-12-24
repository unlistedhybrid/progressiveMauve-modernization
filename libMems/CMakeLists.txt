# ============================================================
# libMems/CMakeLists.txt
# ============================================================

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
)

# On Windows (MSVC), the /openmp flag handles everything.
# Explicitly linking OpenMP::OpenMP_CXX causes LNK1104 errors.
if(NOT MSVC)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(libMems PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()