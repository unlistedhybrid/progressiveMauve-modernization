/////////////////////////////////////////////////////////////////////////////
// File:            libGenome/gnDefs.h
// Purpose:         Defines common constants in libGenome.
// Description:     Modernized for C++17: fixed illegal typedefs, replaced
//                  legacy integer types with <cstdint>, removed BOOL hacks,
//                  preserved public API constants and enums.
// Author:          Aaron Darling (modernized version)
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdint>
#include <cfloat>
#include <climits>
#include <string>

#include "libGenome/gnSetup.h"

namespace genome {

// -----------------------------
// Modern fixed‑width typedefs
// -----------------------------

using int1   = std::uint8_t;        // kept for legacy compatibility
using int8   = std::int8_t;
using int16  = std::int16_t;
using int32  = std::int32_t;
using int64  = std::int64_t;

using sint8  = std::int8_t;
using sint16 = std::int16_t;
using sint32 = std::int32_t;
using sint64 = std::int64_t;

using uint8  = std::uint8_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

using uint = unsigned int;

// floats
using float32 = float;
using float64 = double;

// -----------------------------
// Sequence types
// -----------------------------

using gnSeqC = char;        // sequence character
using gnSeqI = uint64;      // sequence index

constexpr uint32 GNSEQI_ERROR = UINT32_MAX;
constexpr uint32 GNSEQI_END   = UINT32_MAX;
constexpr uint64 GNSEQI_BEGIN = 0;
constexpr int8   GNSEQC_NULL  = 0;
constexpr int8   GNSEQC_MIN   = INT8_MIN;
constexpr int8   GNSEQC_MAX   = INT8_MAX;

// -----------------------------
// Constants
// -----------------------------

constexpr double PI = 3.14159265358979323846;

constexpr uint32 CONTIG_SECTION_SIZE = 3;

constexpr uint32 ALL_CONTIGS = UINT32_MAX;
constexpr uint32 BUFFER_SIZE = 100000;

// -----------------------------
// Enums
// -----------------------------

enum class gnContigSection : uint8 {
    Header      = 0,
    Annotation  = 1,
    Sequence    = 2
};

enum class gnNewlineType : uint8 {
    Unix     = 0,
    Windows  = 1,
    Mac      = 2
};

// -----------------------------
// Utility templates
// -----------------------------

template <typename T>
inline T absolut(const T& t) {
    return (t < 0) ? -t : t;
}

// -----------------------------
// Array RAII helper (legacy API)
// Modern users should prefer std::unique_ptr<T[]>
// But Muscle / progressiveMauve still uses this class.
// -----------------------------

template<class T>
class Array {
public:
    explicit Array(uint64 bufsize)
        : data(new T[bufsize]) {}

    ~Array() {
        delete[] data;
    }

    T* data;

private:
    Array(const Array&) = delete;
    Array& operator=(const Array&) = delete;
    Array() = delete;
};

// -----------------------------

#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ "Unknown()"
#endif

} // namespace genome
