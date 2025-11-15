file(GLOB LIBMEMS_SRC
    libMems/*.cpp
)

add_library(libMems STATIC ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
)

target_link_libraries(libMems PUBLIC libGenome)
