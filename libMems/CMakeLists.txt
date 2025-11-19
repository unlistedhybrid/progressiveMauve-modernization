set(LIBMEMS_SRC
    libMems/Aligner.cpp
    libMems/Backbone.cpp
    libMems/DNAFileSML.cpp
    libMems/DNAMemorySML.cpp
    libMems/FileSML.cpp
    libMems/GappedAlignment.cpp
    libMems/GreedyBreakpointElimination.cpp

    libMems/HomologyHMM/algebras.cpp
    libMems/HomologyHMM/homology.cpp
    libMems/HomologyHMM/homologymain.cpp

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
    libMems/PhyloTree.cpp
    libMems/ProgressiveAligner.cpp
    libMems/RepeatHash.cpp
    libMems/SlotAllocator.cpp
    libMems/SortedMerList.cpp
    libMems/SuperInterval.cpp
    libMems/gnAlignedSequences.cpp
    libMems/twister.cpp

    libMems/dmSML/asyncio.cpp
    libMems/dmSML/aPOSIXaio.cpp
    libMems/dmSML/alinuxaio.cpp
    libMems/dmSML/awin32aio.cpp
    libMems/dmSML/alibc.cpp
    libMems/dmSML/buffer.cpp
    libMems/dmSML/sml.cpp
    libMems/dmSML/sorting.cpp
    libMems/dmSML/timing.cpp
    libMems/dmSML/util.cpp
    libMems/dmSML/dmsort.cpp
)

add_library(libMems STATIC ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
    ${PROJECT_SOURCE_DIR}/muscle
)

# Link against our robust BoostInterface
target_link_libraries(libMems PUBLIC 
    libGenome
    BoostInterface
)
