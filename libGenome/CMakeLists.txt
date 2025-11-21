# libGenome/CMakeLists.txt

# Explicitly list sources from filelist.txt [cite: 1-8]
set(LIBGENOME_SRC
    libGenome/gnBaseFeature.cpp
    libGenome/gnCompare.cpp
    libGenome/gnDefs.cpp
    libGenome/gnException.cpp
    libGenome/gnFASSource.cpp
    libGenome/gnFeature.cpp
    libGenome/gnFileContig.cpp
    libGenome/gnFileSource.cpp
    libGenome/gnFilter.cpp
    libGenome/gnFragmentSpec.cpp
    libGenome/gnGBKSource.cpp
    libGenome/gnGenomeSpec.cpp
    libGenome/gnLocation.cpp
    libGenome/gnRAWSource.cpp
    libGenome/gnSEQSource.cpp
    libGenome/gnSequence.cpp
    libGenome/gnSourceFactory.cpp
    libGenome/gnSourceHeader.cpp
    libGenome/gnSourceQualifier.cpp
    libGenome/gnSourceSpec.cpp
    libGenome/gnStringSpec.cpp
    libGenome/gnStringTools.cpp
)

add_library(libGenome STATIC ${LIBGENOME_SRC})

# Allow including files via #include "libGenome/gnDefs.h" OR #include "gnDefs.h"
target_include_directories(libGenome PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libGenome
)

target_link_libraries(libGenome PUBLIC BoostInterface)
