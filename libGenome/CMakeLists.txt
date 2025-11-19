set(LIBGENOME_SRC
    libGenome/coordMapper.cpp
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
    libGenome/gnStringTools.cpp
    libGenome/gnStringSpec.cpp
    libGenome/gnTranslator.cpp
)

add_library(libGenome STATIC ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
)

# Link directly to Boost targets to ensure include paths are propagated on Windows
target_link_libraries(libGenome PUBLIC
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
)
