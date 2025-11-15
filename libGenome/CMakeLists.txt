set(LIBGENOME_SRC
    libGenome/gnABISource.cpp
    libGenome/gnBaseFeature.cpp
    libGenome/gnCompare.cpp
    libGenome/gnContigSpec.cpp
    libGenome/gnDefs.cpp
    libGenome/gnException.cpp
    libGenome/gnFASSource.cpp
    libGenome/gnFastTranslator.cpp
    libGenome/gnFeature.cpp
    libGenome/gnFileContig.cpp
    libGenome/gnFileSource.cpp
    libGenome/gnFilter.cpp
    libGenome/gnFragmentSpec.cpp
    libGenome/gnGBKSource.cpp
    libGenome/gnGenomeSpec.cpp
    libGenome/gnLocation.cpp
    libGenome/gnPosSpecificTranslator.cpp
    libGenome/gnRAWSource.cpp
    libGenome/gnSEQSource.cpp
    libGenome/gnSequence.cpp
    libGenome/gnSourceFactory.cpp
    libGenome/gnSourceHeader.cpp
    libGenome/gnSourceQualifier.cpp
    libGenome/gnSourceSpec.cpp
    libGenome/gnStringSpec.cpp
    libGenome/gnStringTools.cpp
    libGenome/gnTranslator.cpp
)

# Build the static library
add_library(libGenome STATIC ${LIBGENOME_SRC})

# Include directories
target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link dependencies (Boost already found by top-level CMake)
target_link_libraries(libGenome PUBLIC ${Boost_LIBRARIES})
