/*******************************************************************************
 * $Id: DNAMemorySML.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _DNAMemorySML_h_
#define _DNAMemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SortedMerList.h"
#include "libMems/dmSML/sml.h"
#include <vector>

namespace mems {

class DNAMemorySML : public mems::MemorySML
{
public:
	DNAMemorySML(const std::uint8_t* table = SortedMerList::BasicDNATable(), const std::uint32_t alpha_bits = sml::DNA_ALPHA_BITS);
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

