/*******************************************************************************
 * $Id: RepeatHashCat.h,v 1.2 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _RepeatHashCat_h_
#define _RepeatHashCat_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatHash.h"
#include <vector>
#include <cstdint>

namespace mems {

/**
 * RepeatHashCat extends RepeatHash to handle repeat finding in concatenated sequences.
 * This class tracks contig boundaries across concatenated genomic sequences and
 * properly adjusts match positions accordingly.
 */
class RepeatHashCat : public RepeatHash
{
public:
	RepeatHashCat();
	virtual ~RepeatHashCat();
	
	RepeatHashCat(const RepeatHashCat& rhc);
	RepeatHashCat& operator=(const RepeatHashCat& rhc);
	
	[[nodiscard]] virtual RepeatHashCat* Clone() const;
	
	/**
	 * Tracks the starting position of each contig in the concatenated sequence.
	 * This is used to map hash match positions back to individual contigs.
	 * @return The concatenated contig start positions
	 */
	const std::vector<uint32>& GetConcatContigStart() const { return concat_contig_start; }
	
	/**
	 * Sets the concatenated contig start positions.
	 * @param starts Vector of contig start positions in the concatenated sequence
	 */
	void SetConcatContigStart(const std::vector<uint32>& starts) { concat_contig_start = starts; }

protected:
	// Tracks where each contig starts in the concatenated sequence
	// Used to properly map repeat positions back to individual contigs
	std::vector<uint32> concat_contig_start;
};

} // namespace mems

#endif // _RepeatHashCat_h_
