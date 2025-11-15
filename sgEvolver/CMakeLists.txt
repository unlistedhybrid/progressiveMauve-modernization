file(GLOB SGEVOLVER_SRC
    "src/*.cpp"
    "src/*.c"
)

add_library(sgEvolver ${SGEVOLVER_SRC})

target_include_directories(sgEvolver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
