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

#include "libMems/MemHash.h"        // MUST be included before using MemHash
#include "libMUSCLE/threadstorage.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

namespace mems {

/**
 * ParallelMemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences.
 *
 * The OpenMP version uses thread‑local hash tables for parallel match discovery.
 */
class ParallelMemHash : public MemHash {
public:

    // ────────────────────────────────────────────────────────────────
    // Constructors / Assignment
    // ────────────────────────────────────────────────────────────────

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

    // Clone
    virtual ParallelMemHash* Clone() const override {
        return new ParallelMemHash(*this);
    }

    // ────────────────────────────────────────────────────────────────
    // Main API
    // ────────────────────────────────────────────────────────────────
    virtual void FindMatches(MatchList& match_list) override;

protected:

    // Thread-local override versions of AddHashEntry / MergeTable
    virtual MatchHashEntry* AddHashEntry(MatchHashEntry& mhe) override;
    virtual void MergeTable() override;

    // Each thread gets its own temporary hash table
    TLS<std::vector<std::vector<MatchHashEntry*>>> thread_mem_table;
};

} // namespace mems


// ======================================================================
// Non-OpenMP fallback
// ======================================================================
#ifndef _OPENMP

namespace mems {

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

    // Uses base class implementation
};

} // namespace mems

#endif // !_OPENMP

#endif // _ParallelMemHash_h_
