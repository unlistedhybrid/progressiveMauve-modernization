# ============================================================
# libGenome/CMakeLists.txt
# ============================================================

find_package(OpenMP REQUIRED)

file(GLOB LIBGENOME_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.cpp"
)

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/"
)

target_link_libraries(libGenome PUBLIC
    BoostInterface
    OpenMP::OpenMP_CXX
)

# CRITICAL FIX: Explicitly force the compilation flags
if(OpenMP_CXX_FOUND)
    if(MSVC)
        target_compile_options(libGenome PRIVATE /openmp)
    elseif(APPLE)
        # Apple Clang specific flags
        target_compile_options(libGenome PRIVATE -Xpreprocessor -fopenmp)
    else()
        # GCC / Linux
        target_compile_options(libGenome PRIVATE -fopenmp)
    endif()
endif()