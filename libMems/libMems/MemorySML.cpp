/*******************************************************************************
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing details.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemorySML.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include "libMems/dmSML/sml.h"

namespace mems {

MemorySML::MemorySML(const std::uint8_t* table, std::uint32_t alpha_bits) {
    header.alphabet_bits = alpha_bits;
    std::memcpy(header.translation_table, table, sml::UINT8_MAX_VALUE);
    header.version = 0;
}

MemorySML::MemorySML(const MemorySML& msa) : SortedMerList(msa) {
    positions = msa.positions;
}

MemorySML& MemorySML::operator=(const MemorySML& msa) {
    if (this != &msa) {
        SortedMerList::operator=(msa);
        positions = msa.positions;
    }
    return *this;
}

MemorySML* MemorySML::Clone() const {
    return new MemorySML(*this);
}

void MemorySML::Clear() {
    SortedMerList::Clear();
    positions.clear();
}

void MemorySML::Create(const genome::gnSequence& seq, uint64_t seed) {
    SortedMerList::Create(seq, seed);

    std::vector<bmer> sml_array;
    bool is_spaced_seed = header.seed_length != header.seed_weight;
    if (is_spaced_seed)
        FillDnaSeedSML(seq, sml_array);
    else
        FillSML(seq, sml_array);

    std::sort(sml_array.begin(), sml_array.end(), &bmer_lessthan);

    positions.clear();
    positions.reserve(sml_array.size());
    for (gnSeqI merI = 0; merI < static_cast<gnSeqI>(sml_array.size()); ++merI) {
        positions.push_back(sml_array[merI].position);
    }
}

bool MemorySML::Read(std::vector<bmer>& readVector, gnSeqI size, gnSeqI offset) {
    readVector.clear();
    if (offset > positions.size())
        return false;

    gnSeqI last_mer = offset + size;
    bool success = true;
    if (last_mer > positions.size()) {
        last_mer = positions.size();
        success = false;
    }

    bmer cur_mer;
    for (gnSeqI merI = offset; merI < last_mer; ++merI) {
        cur_mer.position = positions[merI];
        cur_mer.mer = GetSeedMer(cur_mer.position);
        readVector.push_back(cur_mer);
    }
    return success;
}

void MemorySML::Merge(SortedMerList& /*sa*/, SortedMerList& /*sa2*/) {
    // Intentionally left empty or TODO
}

bmer MemorySML::operator[](gnSeqI index) {
    bmer cur_mer;
    cur_mer.position = positions[index];
    cur_mer.mer = GetSeedMer(cur_mer.position);
    return cur_mer;
}

} // namespace mems
