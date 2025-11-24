#ifndef _MemorySML_h_
#define _MemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SortedMerList.h"
#include "libMems/dmSML/sml.h" // Needed for namespace constants and types
#include <vector>
#include <cstdint> // Needed for std::uint*_t

namespace mems {

class MemorySML : public SortedMerList
{
public:
    // FIX: Use std::uint8_t, std::uint32_t, and the new non-colliding constant.
	MemorySML(const std::uint8_t* table = SortedMerList::BasicDNATable(), const std::uint32_t alpha_bits = sml::SML_DNA_ALPHA_BITS);
	MemorySML(const MemorySML& msa);
	MemorySML& operator=(const MemorySML& msa );
	MemorySML* Clone() const;
	
	virtual void Clear();

    // FIX: Use std::uint64_t
	virtual void Create(const genome::gnSequence& seq, const std::uint64_t seed);
    // FIX: Use bool instead of boolean (assuming gnAlignedSequences fixed this)
	virtual bool Read(std::vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0); 
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2);
	
	virtual bmer operator[](gnSeqI index);
	
protected:

	std::vector<sml::position_t> positions; // FIX: Use sml::position_t

};

}

#endif //_MemorySML_h_
