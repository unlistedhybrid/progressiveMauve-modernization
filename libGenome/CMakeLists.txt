# ============================================================
# libGenome/CMakeLists.txt
# ============================================================

file(GLOB LIBGENOME_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.cpp"
)

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/"
)

target_link_libraries(libGenome PUBLIC
    BoostInterface
)

# On Windows (MSVC), the /openmp flag handles everything.
# Explicitly linking OpenMP::OpenMP_CXX causes LNK1104 'genomeomp.lib' errors.
if(NOT MSVC)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(libGenome PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()