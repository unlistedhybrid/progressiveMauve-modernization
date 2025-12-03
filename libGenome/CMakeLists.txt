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

# 1. Link the library (Handles linking)
target_link_libraries(libGenome PUBLIC
    BoostInterface
    OpenMP::OpenMP_CXX
)

# 2. FORCE the Compile Options (Handles the -fopenmp flag)
if(OpenMP_CXX_FOUND)
    target_compile_options(libGenome PRIVATE ${OpenMP_CXX_FLAGS})
endif()