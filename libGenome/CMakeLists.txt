# CMakeLists.txt for libGenome (cleaned version)

# Collect all .cpp and .c source files
file(GLOB LIBGENOME_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.c"
)

# Remove obsolete test/debug programs that should NOT be built
list(REMOVE_ITEM LIBGENOME_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testTests.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/test-o-matic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testSource.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testSourceSeq.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testgnSequence.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testgnSeqConverter.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/coordMapper.cpp"  # optional: remove coordMapper too
)

# list(REMOVE_ITEM LIBGENOME_SRC "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/coordMapper.cpp")

# Build libGenome as a static library
add_library(libGenome STATIC ${LIBGENOME_SRC})

# Include directories used across libGenome
target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link dependencies if needed
target_link_libraries(libGenome
    PUBLIC
        ${Boost_LIBRARIES}
)
