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

# dmSML subfolder sources [cite: 16, 17]
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

# HomologyHMM subfolder sources [cite: 18, 19]
set(MEMS_HMM_SRC
    libMems/HomologyHMM/algebras.cpp
    libMems/HomologyHMM/homology.cpp
)

add_library(libMems STATIC ${MEMS_MAIN_SRC} ${MEMS_DMSML_SRC} ${MEMS_HMM_SRC})

# Include paths for libMems and its submodules
target_include_directories(libMems PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libMems
)

# libMems depends on Genome and MUSCLE
target_link_libraries(libMems PUBLIC
    libGenome
    libMUSCLE
    BoostInterface
)
