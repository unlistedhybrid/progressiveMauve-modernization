#ifndef _DNAMemorySML_h_
#define _DNAMemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// FIX: Must include the base class of MemorySML (SortedMerList) first,
// and then MemorySML itself, to resolve the incomplete type dependency.
#include "libMems/SortedMerList.h" 
#include "libMems/MemorySML.h"    

#include "libGenome/gnSequence.h"
#include "libMems/dmSML/sml.h"
#include <vector>
#include <cstdint>

namespace mems {

// The previous namespace qualification 'mems::MemorySML' is correct, 
// and the error is resolved by the updated include order above.
class DNAMemorySML : public mems::MemorySML
{
public:
    // This signature uses the new, non-colliding constant sml::SML_DNA_ALPHA_BITS
	DNAMemorySML(const std::uint8_t* table = SortedMerList::BasicDNATable(), const std::uint32_t alpha_bits = sml::SML_DNA_ALPHA_BITS);
	DNAMemorySML(const DNAMemorySML& msa);
	DNAMemorySML(const SortedMerList& sa);
	DNAMemorySML& operator=(const DNAMemorySML& msa );
	DNAMemorySML* Clone() const;
	
	
	virtual std::uint64_t GetMer(gnSeqI offset) const;
	virtual std::uint64_t GetSeedMer( gnSeqI offset ) const;
	
protected:

	virtual void FillSML(const genome::gnSequence& seq, std::vector<bmer>& sml_array);

};

}

#endif
