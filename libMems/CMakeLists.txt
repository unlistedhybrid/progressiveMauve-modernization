set(LIBMEMS_SRC
    libMems/Aligner.cpp
    libMems/Backbone.cpp
    libMems/GappedAlignment.cpp
    libMems/GreedyBreakpointElimination.cpp
    libMems/Interval.cpp
    libMems/IntervalList.cpp
    libMems/Islands.cpp
    libMems/MaskedMemHash.cpp
    libMems/MatchFinder.cpp
    libMems/MatchHashEntry.cpp
    libMems/MatchList.cpp
    libMems/MemHash.cpp
    libMems/MuscleInterface.cpp
    libMems/PhyloTree.cpp
    libMems/ProgressiveAligner.cpp
    libMems/SlotAllocator.cpp
    libMems/SortedMerList.cpp
    libMems/SuperInterval.cpp
    libMems/gnAlignedSequences.cpp
    libMems/twister.cpp
)

add_library(libMems STATIC ${LIBMEMS_SRC})
target_include_directories(libMems PUBLIC
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
)
target_link_libraries(libMems PUBLIC libGenome)
