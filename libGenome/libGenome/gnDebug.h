/////////////////////////////////////////////////////////////////////////////
// File:            libGenome/gnDebug.h
// Purpose:         Debug header used for libGenome (C++17 modernized)
// Author:          Aaron Darling (modernized)
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"
#include <string>
#include <iostream>

#if defined(_WIN32)
#include <intrin.h>
#endif

namespace genome {

/** Traps into the MSVC debugger */
inline void breakHere() {
#if defined(_WIN32)
    __debugbreak();
#endif
}

// Determine command-line/GUI mode
#if defined(COMMAND_LINE) || defined(_CONSOLE)
constexpr bool USE_COMMAND_LINE = true;
constexpr bool USE_GUI = false;
#elif defined(GN_GUI)
constexpr bool USE_COMMAND_LINE = false;
constexpr bool USE_GUI = true;
#else
constexpr bool USE_COMMAND_LINE = false;
constexpr bool USE_GUI = false;
#endif

// Debug message functions
inline void DebugMsg(const std::string& msg) {
    if (USE_COMMAND_LINE) {
        std::cout << msg;
    } else if (USE_GUI) {
        // GUI message logic can be added here
    }
}

inline void ErrorMsg(const std::string& msg) {
    if (USE_COMMAND_LINE) {
        std::cerr << msg;
    } else if (USE_GUI) {
        // GUI error logic can be added here
    }
}

} // namespace genome
