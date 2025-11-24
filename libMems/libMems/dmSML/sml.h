#ifndef _sml_h_
#define _sml_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"
#include <cstring>
#include <cstdint>
#include <memory>
#include <string>
#include "libGenome/gnDefs.h"

namespace sml {

using uint8     = std::uint8_t;
using uint32    = std::uint32_t;
using uint64    = std::uint64_t;
using position_t = uint32;
using mask_t    = uint64;
using sarID_t   = std::int16_t;

constexpr std::size_t UINT8_MAX_VALUE = 256;
constexpr std::size_t SML_DESCRIPTION_SIZE = 2048; 
constexpr std::uint32_t SML_DNA_ALPHA_BITS = 2;
constexpr mask_t SEED_MASK = 0x7FFFFFFF;
constexpr int MASK_T_BYTES = 8;

inline mask_t seed_mask = SEED_MASK;
inline int mask_length = 31;
inline int mask_weight = 31;

inline uint8* CreateBasicDNATable() {
    auto bdt = std::make_unique<uint8[]>(UINT8_MAX_VALUE);
    std::memset(bdt.get(), 0, UINT8_MAX_VALUE);

    bdt['c'] = 1;  bdt['C'] = 1;  bdt['b'] = 1;  bdt['B'] = 1;  bdt['y'] = 1;  bdt['Y'] = 1;
    bdt['g'] = 2;  bdt['G'] = 2;  bdt['s'] = 2;  bdt['S'] = 2;  bdt['k'] = 2;  bdt['K'] = 2;
    bdt['t'] = 3;  bdt['T'] = 3;
    return bdt.release();
}

struct SMLHeader_t {
    uint32 version;
    uint32 alphabet_bits;
    uint64 seed;
    uint32 seed_length;
    uint32 seed_weight;
    uint64 length;
    uint32 unique_mers;
    uint32 word_size;
    bool little_endian;
    sarID_t id;
    bool circular;
    uint8 translation_table[UINT8_MAX_VALUE];
    char description[SML_DESCRIPTION_SIZE];
};

struct sml_t {
    char key[8];
    position_t pos;
};

SMLHeader_t InitSML(aFILE* file, uint64 file_size, uint64 seed);

}

#endif
