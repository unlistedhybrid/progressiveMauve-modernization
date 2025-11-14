file(GLOB LIBMEMS_SRC
    "*.cpp"
    "*.c"
)

add_library(libMems ${LIBMEMS_SRC})

target_include_directories(libMems PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libMems PUBLIC libGenome ${Boost_LIBRARIES})
