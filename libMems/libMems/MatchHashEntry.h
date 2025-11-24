/*******************************************************************************
 * $Id: Match.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * Modernized C++17 version — behavior preserved.
 ******************************************************************************/

#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libMems/Match.h"

#include <set>
#include <cstdint>
#include <iostream>

namespace mems {

/**
 * MatchHashEntry stores the location of a fixed-size matching block ("MEM")
 * that occurs across multiple sequences.
 *
 * It extends Match and adds bookkeeping for offset computation and extension.
 */
class MatchHashEntry : public Match {
public:
    enum class MemType {
        seed,
        extended
    };

public:
    MatchHashEntry();

    MatchHashEntry(uint seq_count,
                   gnSeqI mersize,
                   MemType m_type = MemType::seed);

    MatchHashEntry(const MatchHashEntry& mhe) : Match(mhe) { }
    MatchHashEntry& operator=(const MatchHashEntry& mhe);

    MatchHashEntry* Clone() const { return new MatchHashEntry(*this); }

    // SlotAllocator-based clone/free interface
    MatchHashEntry* Copy() const { return m_allocateAndCopy(*this); }
    void Free() { m_free(this); }

    // Whether this seed has been extended to a full MEM
    bool Extended() const { return m_extended; }
    void SetExtended(bool extended) { m_extended = extended; }

    // Mer-size used to detect this match
    uint MerSize() const { return m_mersize; }

    // Required: recompute offset after changing coordinates
    void CalculateOffset();

    int64_t Offset() const { return m_offset; }
    void SetOffset(int64_t offset) { m_offset = offset; }

    // Comparisons
    bool operator==(const MatchHashEntry& mhe) const;

    static bool offset_lessthan(const MatchHashEntry& a,
                                const MatchHashEntry& b);

    static bool start_lessthan(const MatchHashEntry& a,
                               const MatchHashEntry& b);

    static bool start_lessthan_ptr(const MatchHashEntry* a,
                                   const MatchHashEntry* b);

    static bool strict_start_lessthan_ptr(const MatchHashEntry* a,
                                          const MatchHashEntry* b);

    static int64_t end_to_start_compare(const MatchHashEntry& a,
                                        const MatchHashEntry& b);

    static int64_t start_compare(const MatchHashEntry& a,
                                 const MatchHashEntry& b);

    // Containment test ("a contains b")
    bool Contains(const MatchHashEntry& mhe) const;

private:
    bool     m_extended = false;
    gnSeqI   m_mersize  = 0;
    int64_t  m_offset   = 0;
};


// Convenience non-pointer comparator
inline bool MatchHashEntry::start_lessthan(const MatchHashEntry& a,
                                           const MatchHashEntry& b)
{
    return start_lessthan_ptr(&a, &b);
}


// Sorting functor for pointer-to-MatchHashEntry
class MheCompare {
public:
    bool operator()(const MatchHashEntry* a,
                    const MatchHashEntry* b) const
    {
        // Primary: smaller "FirstStart()" sorts first
        if (a->FirstStart() != b->FirstStart())
            return a->FirstStart() < b->FirstStart();

        // Secondary: presence/absence of coordinates per sequence
        uint n = std::min(a->SeqCount(), b->SeqCount());
        for (uint i = 0; i < n; i++) {
            bool a_missing = (a->LeftEnd(i) == NO_MATCH);
            bool b_missing = (b->LeftEnd(i) == NO_MATCH);

            if (a_missing != b_missing)
                return b_missing; // defined > undefined
        }

        // Tertiary: containment
        if (a->Contains(*b) || b->Contains(*a))
            return false;

        // Last resort: strict start comparison
        return MatchHashEntry::strict_start_lessthan_ptr(a, b);
    }
};

} // namespace mems

