file(GLOB LIBMEMS_SRC
    "libMems/*.cpp"
    "libMems/*.c"
)

add_library(libMems ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/libMems
)
