file(GLOB LIBGENOME_SRC
    libGenome/*.cpp
)

# Remove test programs
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "test")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "Test")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "testgn")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "testSource")
list(FILTER LIBGENOME_SRC EXCLUDE REGEX "test-o-matic")

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
)
