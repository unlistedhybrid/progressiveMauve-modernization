set(PROGRESSIVEMAUVE_SRC
    src/mauveAligner.cpp
    src/progressiveMauve.cpp
    src/ProgressiveAligner.cpp
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

# Utility programs (optional - uncomment if needed)
# add_executable(xmfa2maf src/xmfa2maf.cpp)
# target_include_directories(xmfa2maf PUBLIC ${PROJECT_SOURCE_DIR}/mauveAligner/src)
# target_link_libraries(xmfa2maf PUBLIC libGenome libMems)
