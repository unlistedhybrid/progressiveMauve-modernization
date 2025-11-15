file(GLOB ALIGNER_SRC
    "src/*.cpp"
    "src/*.c"
)

add_library(mauveAligner STATIC ${ALIGNER_SRC})

target_include_directories(mauveAligner PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/libGenome
    ${PROJECT_SOURCE_DIR}/libGenome/libGenome
    ${PROJECT_SOURCE_DIR}/libMems
    ${PROJECT_SOURCE_DIR}/libMems/libMems
)

add_executable(progressiveMauve src/progressiveMauve.cpp)

target_link_libraries(progressiveMauve
    PRIVATE
        mauveAligner
        libMems
        libGenome
        muscle
        sgEvolver
        ${Boost_LIBRARIES}
)
