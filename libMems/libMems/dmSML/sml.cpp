#include <iostream> // for std::cerr
#include <cstring>
#include <limits>

#include "libMems/dmSML/sml.h"
#include "libMems/SeedMasks.h"

namespace sml {

SMLHeader_t InitSML(aFILE* file, uint64 file_size, uint64 seed) {
    SMLHeader_t header{};
    header.version = 5;
    header.alphabet_bits = 2;
    header.seed = seed;
    header.seed_length = getSeedLength(seed);
    header.seed_weight = getSeedWeight(seed);
    header.length = file_size;
    header.unique_mers = std::numeric_limits<uint32_t>::max(); // "unset"
    header.word_size = 32;
    header.little_endian = true;
    header.id = 0;
    header.circular = false;
    std::memset(header.translation_table, 0, UINT8_MAX_VALUE);

    // Fill DNA translation table
    {
        auto bdt = CreateBasicDNATable();
        std::memcpy(header.translation_table, bdt.get(), UINT8_MAX_VALUE);
    }
    header.description[0] = '\0';

    int retcode = aWrite(static_cast<void*>(&header), sizeof(header), 1, file, 0);
    if (retcode == 0) {
        std::cerr << "Error writing to SML\n";
    }
    aWaitComplete(file, retcode);
    return header;
}

} // namespace sml
