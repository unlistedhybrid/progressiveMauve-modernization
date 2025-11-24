#ifndef _sml_h_
#define _sml_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"
#include <cstring>
#include <cstdint>      // For standard integer types
#include <memory>       // For std::unique_ptr
#include <string>
#include "libGenome/gnDefs.h" // For boolean, etc.

namespace sml {

// --- Type aliases for clarity ---
using uint8     = std::uint8_t;
using uint32    = std::uint32_t;
using uint64    = std::uint64_t;
using position_t = uint32;
using mask_t    = uint64;
using sarID_t   = std::int16_t;

// --- Constants ---
constexpr std::size_t UINT8_MAX_VALUE = 256;
constexpr std::size_t SML_DESCRIPTION_SIZE = 2048;
constexpr int MASK_T_BYTES = 8;

// Original default mask values
constexpr mask_t DEFAULT_SEED_MASK = 0x7FFFFFFF;
constexpr int DEFAULT_MASK_LENGTH = 31;
constexpr int DEFAULT_MASK_WEIGHT = 31;

// --- Masking parameters (Inline variables allow definition in the header) ---
// These replace the old static global variables
inline mask_t seed_mask = DEFAULT_SEED_MASK;
inline int mask_length = DEFAULT_MASK_LENGTH;
inline int mask_weight = DEFAULT_MASK_WEIGHT;

inline uint8* CreateBasicDNATable() {
    // Note: This function now uses the standard C++ library for memory allocation, 
    // but returns a raw pointer which is prone to leaks if not manually managed 
    // (the original code relies on this pattern).
    auto bdt = std::make_unique<uint8[]>(UINT8_MAX_VALUE);
    std::memset(bdt.get(), 0, UINT8_MAX_VALUE);

    bdt['c'] = 1;  bdt['C'] = 1;  bdt['b'] = 1;  bdt['B'] = 1;  bdt['y'] = 1;  bdt['Y'] = 1;
    bdt['g'] = 2;  bdt['G'] = 2;  bdt['s'] = 2;  bdt['S'] = 2;  bdt['k'] = 2;  bdt['K'] = 2;
    bdt['t'] = 3;  bdt['T'] = 3;
    // We must return the raw pointer to satisfy the assignment in dmsort.cpp:
    // DNA_TABLE = CreateBasicDNATable();
    return bdt.release(); // release ownership to the caller (dmsort.cpp)
}

// --- SML Header struct ---
struct SMLHeader_t {
    uint32 version;                         // Format version - 4 bytes
    uint32 alphabet_bits;                   // Bits per character in the alphabet - 4 bytes
    uint64 seed;                            // The pattern used in each seed
    uint32 seed_length;                     // The length of the seed mask
    uint32 seed_weight;                     // The weight of the seed mask
    uint64 length;                          // Length of the sequence before circularity - 8 bytes
    uint32 unique_mers;                     // Number of unique mers in the sequence 4 bytes
    uint32 word_size;                       // Word size on the machine the sequence was translated
    bool little_endian;                     // Is the byte order little endian?
    sarID_t id;                             // Obsolete ID value
    bool circular;                          // Circularity of sequence
    uint8 translation_table[UINT8_MAX_VALUE];    // Translation table for ascii chars to binary values -- 256 bytes
    char description[SML_DESCRIPTION_SIZE];      // Freeform text description of sequence data -- 2048 bytes
};

// --- SML record struct ---
struct sml_t {
    char key[8];
    position_t pos;
};

// --- SMLHeader init prototype (should be defined elsewhere) ---
SMLHeader_t InitSML(aFILE* file, uint64 file_size, uint64 seed);

} // namespace sml

#endif // _sml_h_
