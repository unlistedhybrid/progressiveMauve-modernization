file(GLOB ALIGNER_SRC
    "*.cpp"
    "*.c"
)

add_library(mauveAligner ${ALIGNER_SRC})

target_include_directories(mauveAligner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mauveAligner
    PUBLIC libGenome libMems muscle sgEvolver ${Boost_LIBRARIES}
)
