/*******************************************************************************
 * $Id: ParallelMemHash.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007
 * Aaron Darling and authors listed in the AUTHORS file.
 * Licensed under GPL.
 ******************************************************************************/

#ifndef _ParallelMemHash_h_
#define _ParallelMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"
#include "libMUSCLE/threadstorage.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

namespace mems {

#ifdef _OPENMP

/**
 * ParallelMemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences using OpenMP parallelism.
 *
 * Uses thread-local hash tables for parallel match discovery.
 */
class ParallelMemHash : public MemHash {
public:

	ParallelMemHash() : MemHash() {}

	ParallelMemHash(const ParallelMemHash& mh)
		: MemHash(mh),
		  thread_mem_table(mh.thread_mem_table)
	{}

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

	virtual void FindMatches(MatchList& match_list) override;

protected:

	virtual MatchHashEntry* AddHashEntry(MatchHashEntry& mhe) override;
	virtual void MergeTable() override;

	TLS<std::vector<std::vector<MatchHashEntry*>>> thread_mem_table;
};

#else

/**
 * When built without OpenMP, ParallelMemHash behaves exactly like MemHash,
 * but maintains API compatibility.
 */
class ParallelMemHash : public MemHash {
public:

	ParallelMemHash() : MemHash() {}

	ParallelMemHash(const ParallelMemHash& mh)
		: MemHash(mh)
	{}

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

#endif

}

#endif
