/******************************************************************************
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 ******************************************************************************/

#ifndef _MemorySML_h_
#define _MemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SortedMerList.h"
#include <vector>
#include <cstdint>

namespace mems {

/**
 * MemorySML: A SortedMerList stored entirely in memory.
 * Consumes roughly 32 + alpha_bits bits per character in the sequences.
 * For unambiguous DNA sequences, requires ~4.25 bytes per base.
 */
class MemorySML : public SortedMerList {
public:
    /**
     * Create an empty MemorySML.
     * @param table The translation table for characters to binary code.
     * @param alpha_bits The number of bits per character (default: DNA settings).
     */
    MemorySML(const uint8_t* table = SortedMerList::BasicDNATable(), uint32_t alpha_bits = DNA_ALPHA_BITS);
    MemorySML(const MemorySML& msa);
    MemorySML& operator=(const MemorySML& msa);
    MemorySML* Clone() const;

    virtual void Clear();

    virtual void Create(const genome::gnSequence& seq, uint64_t seed);
    virtual bool Read(std::vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0);
    virtual void Merge(SortedMerList& sa, SortedMerList& sa2);

    virtual bmer operator[](gnSeqI index);

protected:
    std::vector<smlSeqI_t> positions;
};

} // namespace mems

#endif // _MemorySML_h_
