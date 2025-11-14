file(GLOB LIBGENOME_SRC
    "*.cpp"
    "*.c"
)

add_library(libGenome ${LIBGENOME_SRC})

target_include_directories(libGenome PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libGenome PUBLIC ${Boost_LIBRARIES})
