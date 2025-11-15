file(GLOB MUSCLE_SRC
    "libMUSCLE/*.cpp"
    "libMUSCLE/*.c"
)

add_library(muscle ${MUSCLE_SRC})

target_include_directories(muscle PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/libMUSCLE
)
