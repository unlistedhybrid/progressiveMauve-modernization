# muscle/CMakeLists.txt

# Sources in libMUSCLE subfolder [cite: 21-25]
set(MUSCLE_SRC
    libMUSCLE/alpha.cpp
    libMUSCLE/clust.cpp
    libMUSCLE/distfunc.cpp
    libMUSCLE/enumopts.cpp
    libMUSCLE/intmath.cpp
    libMUSCLE/msa.cpp
    libMUSCLE/muscle.cpp
    libMUSCLE/objscore.cpp
    libMUSCLE/params.cpp
    libMUSCLE/profile.cpp
    libMUSCLE/pwpath.cpp
    libMUSCLE/seq.cpp
    libMUSCLE/seqvect.cpp
    libMUSCLE/textfile.cpp
)

add_library(libMUSCLE STATIC ${MUSCLE_SRC})

target_include_directories(libMUSCLE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/libMUSCLE
)
