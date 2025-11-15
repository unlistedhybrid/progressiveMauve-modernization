file(GLOB MUSCLE_SRC
    libMUSCLE/*.cpp
)

add_library(muscle STATIC ${MUSCLE_SRC})

target_include_directories(muscle PUBLIC
    ${PROJECT_SOURCE_DIR}/muscle
    ${PROJECT_SOURCE_DIR}/muscle/libMUSCLE
)
