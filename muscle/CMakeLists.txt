file(GLOB MUSCLE_SRC
    "*.cpp"
    "*.c"
)

add_library(muscle ${MUSCLE_SRC})

target_include_directories(muscle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(muscle PUBLIC ${Boost_LIBRARIES})
