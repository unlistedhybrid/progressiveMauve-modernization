file(GLOB MAUVEALIGNER_SRC
    src/*.cpp
)

add_executable(progressiveMauve ${MAUVEALIGNER_SRC})

target_include_directories(progressiveMauve PUBLIC
    ${PROJECT_SOURCE_DIR}/mauveAligner
    ${PROJECT_SOURCE_DIR}/mauveAligner/include
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
