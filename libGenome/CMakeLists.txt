file(GLOB LIBGENOME_SRC
    libGenome/*.cpp
)

# Remove test/demo files
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "test")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "Test")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "SeqStringTest")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "SourceFactory")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "RevComp")

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
)
