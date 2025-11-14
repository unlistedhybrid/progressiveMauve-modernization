file(GLOB ALIGNER_SRC
    "*.cpp"
    "*.c"
)

# Build mauveAligner as a library
add_library(mauveAligner STATIC ${ALIGNER_SRC})

target_include_directories(mauveAligner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mauveAligner
    PUBLIC libGenome libMems muscle sgEvolver ${Boost_LIBRARIES}
)

# --- Build the progressiveMauve CLI ---
# progressiveMauve.cpp is inside this folder

add_executable(progressiveMauve progressiveMauve.cpp)

target_link_libraries(progressiveMauve
    PRIVATE
        mauveAligner
        libMems
        libGenome
        muscle
        sgEvolver
        ${Boost_LIBRARIES}
)
