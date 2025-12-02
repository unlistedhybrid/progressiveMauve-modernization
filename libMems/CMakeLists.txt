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
    libMems/ProgressiveAligner.cpp
    libMems/RepeatHash.cpp
    libMems/SortedMerList.cpp
    libMems/SuperInterval.cpp
    libMems/gnAlignedSequences.cpp
    libMems/dmSML/asyncio.c
    libMems/dmSML/aPOSIXaio.c
    libMems/dmSML/alinuxaio.c
    libMems/dmSML/awin32aio.c
    libMems/dmSML/alibc.c
    libMems/dmSML/buffer.c
    libMems/dmSML/sml.c
    libMems/dmSML/sorting.c
    libMems/dmSML/timing.c
    libMems/dmSML/util.c
    libMems/dmSML/dmsort.c
)

add_library(libMems STATIC ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
    ${PROJECT_SOURCE_DIR}/muscle
)

target_link_libraries(libMems PUBLIC
    libGenome
    BoostInterface
)
