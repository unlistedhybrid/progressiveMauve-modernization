# ============================================================
# muscle/CMakeLists.txt
# ============================================================

file(GLOB MUSCLE_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/libMUSCLE/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libMUSCLE/*.c"
)

add_library(libMUSCLE STATIC ${MUSCLE_SRC})

target_include_directories(libMUSCLE PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(libMUSCLE PUBLIC
    BoostInterface
)