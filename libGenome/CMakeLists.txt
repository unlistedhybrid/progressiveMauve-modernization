set(LIBGENOME_SRC
    libGenome/gnAlgebra.cpp
    libGenome/gnBaseContig.cpp
    libGenome/gnBaseFeature.cpp
    libGenome/gnBaseSource.cpp
    libGenome/gnBaseSpreadsheet.cpp
    libGenome/gnBaseSpec.cpp
    libGenome/gnBaseUtilities.cpp
    libGenome/gnCmdLine.cpp
    libGenome/gnContigSequence.cpp
    libGenome/gnContigSpec.cpp
    libGenome/gnDNA.cpp
    libGenome/gnDataRow.cpp
    libGenome/gnDebug.cpp
    libGenome/gnFeature.cpp
    libGenome/gnFeatureList.cpp
    libGenome/gnFeaturePosList.cpp
    libGenome/gnFileContig.cpp
    libGenome/gnFileSource.cpp
    libGenome/gnFileSpec.cpp
    libGenome/gnFilter.cpp
    libGenome/gnFragmentSequence.cpp
    libGenome/gnFragmentSpec.cpp
    libGenome/gnGenomeSpec.cpp
    libGenome/gnLocation.cpp
    libGenome/gnProtein.cpp
    libGenome/gnParser.cpp
    libGenome/gnSequence.cpp
    libGenome/gnSeqConverter.cpp
    libGenome/gnSpreadsheet.cpp
    libGenome/gnStringTools.cpp
    libGenome/gnTranslatedSequence.cpp
)

# Build the library
add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(libGenome
    PUBLIC
        ${Boost_LIBRARIES}
)
