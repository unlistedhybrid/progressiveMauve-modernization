# ============================================================
# libGenome/CMakeLists.txt
# ============================================================

# 1. Find OpenMP to ensure flags (-fopenmp) are applied
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