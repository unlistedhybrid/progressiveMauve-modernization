file(GLOB SGEVOLVER_SRC
    src/*.cpp
)

add_library(sgEvolver STATIC ${SGEVOLVER_SRC})

target_include_directories(sgEvolver PUBLIC
    ${PROJECT_SOURCE_DIR}/sgEvolver
    ${PROJECT_SOURCE_DIR}/sgEvolver/include
    ${PROJECT_SOURCE_DIR}/sgEvolver/src
)

target_link_libraries(sgEvolver PUBLIC
    libGenome
    libMems
)
