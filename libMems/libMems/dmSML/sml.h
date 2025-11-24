#ifndef _sml_h_
#define _sml_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"
#include <cstring>
#include <cstdint> // For standard integer types
#include <memory>  // For std::unique_ptr
#include "libGenome/gnDefs.h" // For boolean, etc.

namespace sml {

// --- Type aliases using standard C++ types ---
using uint8     = std::uint8_t;
using uint32    = std::uint32_t;
using uint64    = std::uint64_t;
using position_t = uint32;
using mask_t    = uint64;
using sarID_t   = std::int16_t;

// --- Constants (Renamed to avoid macro collisions) ---
// UINT8_MAX
constexpr std::size_t UINT8_MAX_VALUE = 256;
// DESCRIPTION_SIZE
constexpr std::size_t SML_DESCRIPTION_SIZE = 2048; 
// DNA_ALPHA_BITS value (used in MemorySML.h)
constexpr std::uint32_t SML_DNA_ALPHA_BITS = 2;

// --- Masking parameters (Made inline for header definition) ---
constexpr int MASK_T_BYTES = 8;
inline mask_t seed_mask = 0x7FFFFFFF;
inline int mask_length = 31;
inline int mask_weight = 31;

// --- Create the basic DNA translation table ---
// Returns a raw pointer whose ownership is released to the caller (dmsort.cpp)
inline uint8* CreateBasicDNATable() {
    auto bdt = std::make_unique<uint8[]>(UINT8_MAX_VALUE);
    std::memset(bdt.get(), 0, UINT8_MAX_VALUE);

    bdt['c'] = 1;  bdt['C'] = 1;  bdt['b'] = 1;  bdt['B'] = 1;  bdt['y'] = 1;  bdt['Y'] = 1;
    bdt['g'] = 2;  bdt['G'] = 2;  bdt['s'] = 2;  bdt['S'] = 2;  bdt['k'] = 2;  bdt['K'] = 2;
    bdt['t'] = 3;  bdt['T'] = 3;
    return bdt.release();
}

// FIX: DNA_TABLE is no longer static in the header; dmsort.cpp defines it globally.
// static uint8* DNA_TABLE; 

typedef struct SMLHeader_s{
	uint32 version;
	uint32 alphabet_bits;
	uint64 seed;
	uint32 seed_length;
	uint32 seed_weight;
	uint64 length;
	uint32 unique_mers;
	uint32 word_size;
	boolean little_endian;
	signed short id;
	boolean circular;
	uint8 translation_table[UINT8_MAX_VALUE]; // Use C++ constant
	char description[SML_DESCRIPTION_SIZE];   // Use renamed constant
} SMLHeader_t;


typedef struct sml_s {
	char key[8];
	position_t pos;
} sml_t;

SMLHeader_t InitSML( aFILE* file, uint64 file_size, uint64 seed );

} // namespace sml

#endif /* _sml_h_ */
