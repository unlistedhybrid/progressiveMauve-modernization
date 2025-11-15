# Gather all real library sources
file(GLOB LIBGENOME_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/*.c"
)

# Explicitly remove ALL legacy test/utility programs by exact name
foreach(f
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testTests.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/test-o-matic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testSource.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testSourceSeq.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testgnSequence.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/testgnSeqConverter.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/gnSeqStringTest.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/gnFastTranslator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libGenome/coordMapper.cpp"
)
    list(REMOVE_ITEM LIBGENOME_SRC "${f}")
endforeach()

# Build the core libGenome library
add_library(libGenome STATIC ${LIBGENOME_SRC})

# Include directories
target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link dependencies
target_link_libraries(libGenome PUBLIC ${Boost_LIBRARIES})
