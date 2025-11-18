set(PROGRESSIVEMAUVE_SRC
    src/mauveAligner.cpp
    src/progressiveMauve.cpp
    src/AlignmentTree.cpp
    src/UniqueMatchFinder.cpp
    src/RepeatHashCat.cpp
    src/alignmentProjector.cpp
    src/backbone_global_to_local.cpp
    src/addUnalignedIntervals.cpp
    src/repeatoire.cpp
)

add_executable(progressiveMauve ${PROGRESSIVEMAUVE_SRC})
target_include_directories(progressiveMauve PUBLIC
    ${PROJECT_SOURCE_DIR}/mauveAligner
    ${PROJECT_SOURCE_DIR}/mauveAligner/src
)
target_link_libraries(progressiveMauve PUBLIC
    libGenome
    libMems
    muscle
    sgEvolver
    Boost::filesystem
    Boost::iostreams
    Boost::program_options
    Boost::regex
)
