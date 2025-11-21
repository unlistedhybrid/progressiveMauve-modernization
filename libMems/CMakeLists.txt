# libMems/CMakeLists.txt

# Main libMems sources [cite: 9-15]
set(MEMS_MAIN_SRC
    libMems/Aligner.cpp
    libMems/Backbone.cpp
    libMems/DNAFileSML.cpp
    libMems/DNAMemorySML.cpp
    libMems/FileSML.cpp
    libMems/GappedAlignment.cpp
    libMems/gnAlignedSequences.cpp
    libMems/GreedyBreakpointElimination.cpp
    libMems/Interval.cpp
    libMems/IntervalList.cpp
    libMems/Islands.cpp
    libMems/MaskedMemHash.cpp
    libMems/MatchFinder.cpp
    libMems/MatchHashEntry.cpp
    libMems/MatchList.cpp
    libMems/MemHash.cpp
    libMems/MemorySML.cpp
    libMems/MuscleInterface.cpp
    libMems/PairwiseMatchFinder.cpp
    libMems/ProgressiveAligner.cpp
    libMems/SortedMerList.cpp
    libMems/SuperInterval.cpp
    libMems/twister.cpp
)

# dmSML subfolder sources [cite: 16-18]
set(MEMS_DMSML_SRC
    libMems/dmSML/alibc.cpp
    libMems/dmSML/aPOSIXaio.cpp
    libMems/dmSML/asyncio.cpp
    libMems/dmSML/awin32aio.cpp
    libMems/dmSML/buffer.cpp
    libMems/dmSML/dmsort.cpp
    libMems/dmSML/sml.cpp
    libMems/dmSML/sorting.cpp
    libMems/dmSML/timing.cpp 
    libMems/dmSML/util.cpp
)

# HomologyHMM subfolder sources [cite: 19]
set(MEMS_HMM_SRC
    libMems/HomologyHMM/algebras.cpp
    libMems/HomologyHMM/homology.cpp
)

add_library(libMems STATIC ${MEMS_MAIN_SRC} ${MEMS_DMSML_SRC} ${MEMS_HMM_SRC})

target_include_directories(libMems PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libMems
)

# SAFETY FIX: Add dmSML as a PRIVATE include.
# This ensures that if libMems code does #include "timing.h", 
# it finds its own version (dmSML/timing.h) before looking at libMUSCLE.
target_include_directories(libMems PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libMems/dmSML
)

target_link_libraries(libMems PUBLIC
    libGenome
    libMUSCLE
    BoostInterface
)
