file(GLOB LIBGENOME_SRC
    libGenome/*.cpp
)

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
)
