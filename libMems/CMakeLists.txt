file(GLOB LIBMEMS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/HomologyHMM/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMems/dmSML/*.c"
)

add_library(libMems STATIC ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/libGenome"
    "${PROJECT_SOURCE_DIR}/muscle"
)

target_link_libraries(libMems PUBLIC
    libGenome
    libMUSCLE
    BoostInterface
)