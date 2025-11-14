file(GLOB SGEVOLVER_SRC
    "*.cpp"
    "*.c"
)

add_library(sgEvolver ${SGEVOLVER_SRC})

target_include_directories(sgEvolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
