file(GLOB LIBGENOME_SRC
    "libGenome/*.cpp"
    "libGenome/*.c"
)

add_library(libGenome ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
)
