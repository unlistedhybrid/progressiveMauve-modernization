set(LIBMEMS_SRC
    libMems/libMems/Aligner.cpp
    libMems/libMems/Backbone.cpp
    libMems/libMems/DNAFileSML.cpp
    libMems/libMems/DNAMemorySML.cpp
    libMems/libMems/FileSML.cpp
    libMems/libMems/GappedAlignment.cpp
    libMems/libMems/GreedyBreakpointElimination.cpp
    libMems/libMems/HomologyHMM/algebras.cpp
    libMems/libMems/HomologyHMM/homology.cpp
    libMems/libMems/HomologyHMM/homologymain.cpp
    libMems/libMems/Interval.cpp
    libMems/libMems/IntervalList.cpp
    libMems/libMems/Islands.cpp
    libMems/libMems/MaskedMemHash.cpp
    libMems/libMems/MatchFinder.cpp
    libMems/libMems/MatchHashEntry.cpp
    libMems/libMems/MatchList.cpp
    libMems/libMems/MemHash.cpp
    libMems/libMems/MuscleInterface.cpp
    libMems/libMems/PairwiseMatchFinder.cpp
    libMems/libMems/PhyloTree.cpp
    libMems/libMems/ProgressiveAligner.cpp
    libMems/libMems/RepeatHash.cpp
    libMems/libMems/SlotAllocator.cpp
    libMems/libMems/SortedMerList.cpp
    libMems/libMems/SuperInterval.cpp
    libMems/libMems/gnAlignedSequences.cpp
    libMems/libMems/twister.cpp
    libMems/libMems/dmSML/asyncio.cpp
    libMems/libMems/dmSML/aPOSIXaio.cpp
    libMems/libMems/dmSML/alinuxaio.cpp
    libMems/libMems/dmSML/awin32aio.cpp
    libMems/libMems/dmSML/alibc.cpp
    libMems/libMems/dmSML/buffer.cpp
    libMems/libMems/dmSML/sml.cpp
    libMems/libMems/dmSML/sorting.cpp
    libMems/libMems/dmSML/timing.cpp
    libMems/libMems/dmSML/util.cpp
    libMems/libMems/dmSML/dmsort.cpp
)

add_library(libMems STATIC ${LIBMEMS_SRC})
target_include_directories(libMems PUBLIC
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
)
target_link_libraries(libMems PUBLIC libGenome)
