#ifndef _DNAMemorySML_h_
#define _DNAMemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SortedMerList.h"
#include "libMems/MemorySML.h" 

#include "libGenome/gnSequence.h"
#include "libMems/dmSML/sml.h"
#include <vector>
#include <cstdint>

namespace mems {

class DNAMemorySML : public mems::MemorySML
{
public:
	DNAMemorySML(const std::uint8_t* table = SortedMerList::BasicDNATable(), const std::uint32_t alpha_bits = sml::SML_DNA_ALPHA_BITS);
	DNAMemorySML(const DNAMemorySML& msa);
	DNAMemorySML(const SortedMerList& sa);
	DNAMemorySML& operator=(const DNAMemorySML& msa );
	DNAMemorySML* Clone() const;

	virtual uint64 GetMer(gnSeqI offset) const;
	virtual uint64 GetSeedMer( gnSeqI offset ) const;
	
protected:

	virtual void FillSML(const genome::gnSequence& seq, std::vector<bmer>& sml_array);

};

}

#endif
