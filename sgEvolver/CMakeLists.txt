set(SGEVOLVER_SRC
    src/Alignment.cpp
    src/Mutator.cpp
    src/gnGapSource.cpp
    src/gnStringSource.cpp
    src/twister.cpp
)

add_library(sgEvolver STATIC ${SGEVOLVER_SRC})
target_include_directories(sgEvolver PUBLIC
    ${PROJECT_SOURCE_DIR}/sgEvolver
    ${PROJECT_SOURCE_DIR}/sgEvolver/src
)
target_link_libraries(sgEvolver PUBLIC
    libGenome
    libMems
)
