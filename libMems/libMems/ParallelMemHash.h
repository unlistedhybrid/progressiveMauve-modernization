/*******************************************************************************
 * $Id: ParallelMemHash.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _ParallelMemHash_h_
#define _ParallelMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP

#include "libMUSCLE/threadstorage.h"
#include <omp.h>
#include "libMems/MemHash.h"

namespace mems {

/**
 * ParallelMemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences.
 */
class ParallelMemHash : public MemHash {
public:
	ParallelMemHash() : MemHash() {}

	ParallelMemHash(const ParallelMemHash& mh)
		: MemHash(mh),
		  thread_mem_table(mh.thread_mem_table) {}

	ParallelMemHash& operator=(const ParallelMemHash& mh)
	{
		if (this != &mh) {
			MemHash::operator=(mh);
			thread_mem_table = mh.thread_mem_table;
		}
		return *this;
	}

	virtual ParallelMemHash* Clone() const override {
		return new ParallelMemHash(*this);
	}

	/**
	 * Finds (in parallel) all matches in the sequences contained by "match_list"
	 * The resulting list of matches is stored within "match_list"
	 */
	virtual void FindMatches(MatchList& match_list) override;

protected:
	virtual MatchHashEntry* AddHashEntry(MatchHashEntry& mhe) override;
	virtual void MergeTable() override;

	TLS<std::vector<std::vector<MatchHashEntry*>>> thread_mem_table;
};

} // namespace mems

#else // _OPENMP  (No OpenMP)

namespace mems {

/**
 * When built without OpenMP, ParallelMemHash is a thin wrapper around MemHash.
 */
class ParallelMemHash : public MemHash {
public:
	ParallelMemHash() : MemHash() {}

	ParallelMemHash(const ParallelMemHash& mh)
		: MemHash(mh) {}

	ParallelMemHash& operator=(const ParallelMemHash& mh)
	{
		if (this != &mh)
			MemHash::operator=(mh);
		return *this;
	}

	virtual ParallelMemHash* Clone() const override {
		return new ParallelMemHash(*this);
	}
};

} // namespace mems

#endif // _OPENMP

#endif //_ParallelMemHash_h_
