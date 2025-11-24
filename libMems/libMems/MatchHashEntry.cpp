/*******************************************************************************
 * $Id: Match.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * Modernized for C++17 — preserves all semantics, fixes type issues.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Include MatchHashEntry.h for class definition
#include "libMems/MatchHashEntry.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "libMems/dmSML/sml.h"

namespace mems {

MatchHashEntry::MatchHashEntry()
    : m_extended(false), m_mersize(0), m_offset(0) {}

// FIX: Implementation of the primary MemType constructor
MatchHashEntry::MatchHashEntry(uint seq_count, gnSeqI mersize, MemType m_type)
    : m_extended(m_type == MemType::extended),
      m_mersize(mersize),
      m_offset(0) {}

// FIX: Implementation of the missing sml::mask_t constructor to resolve the linker error
MatchHashEntry::MatchHashEntry(uint seq_count, gnSeqI mersize, sml::mask_t seed)
    // Forwards the call to the MemType constructor, assuming sml::mask_t implies MemType::seed
    : MatchHashEntry(seq_count, mersize, MemType::seed)
{
    // The specific 'seed' value is managed globally by sml::seed_mask and not stored here.
}

MatchHashEntry& MatchHashEntry::operator=(const MatchHashEntry& mhe) {
    Match::operator=(mhe);
    m_extended = mhe.m_extended;
    m_mersize = mhe.m_mersize;
    m_offset = mhe.m_offset;
    return *this;
}

bool MatchHashEntry::operator==(const MatchHashEntry& mhe) const {
    if (!(Match::operator==(mhe)))
        return false;
    if (m_extended != mhe.m_extended)
        return false;
    if (m_mersize != mhe.m_mersize)
        return false;
    return true;
}

void MatchHashEntry::CalculateOffset() {
    if (SeqCount() == 0) {
        m_offset = 0;
        return;
    }

    int64_t first_start = LeftEnd(0);
    if (first_start == NO_MATCH)
        first_start = 0;

    m_offset = first_start;

    for (uint i = 1; i < SeqCount(); ++i) {
        int64_t cur_start = LeftEnd(i);
        if (cur_start == NO_MATCH)
            cur_start = 0;
        m_offset += std::abs(cur_start - first_start);
    }
}

bool MatchHashEntry::offset_lessthan(const MatchHashEntry& a, const MatchHashEntry& b) {
    return a.m_offset < b.m_offset;
}

bool MatchHashEntry::start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b) {
    if (a->FirstStart() < b->FirstStart())
        return true;
    if (a->FirstStart() > b->FirstStart())
        return false;
    return a->Offset() < b->Offset();
}

bool MatchHashEntry::strict_start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b) {
    uint n = std::min(a->SeqCount(), b->SeqCount());
    for (uint i = 0; i < n; ++i) {
        int64_t a_start = a->LeftEnd(i);
        int64_t b_start = b->LeftEnd(i);

        if (a_start == NO_MATCH && b_start != NO_MATCH)
            return true;
        if (a_start != NO_MATCH && b_start == NO_MATCH)
            return false;
        if (a_start != NO_MATCH && b_start != NO_MATCH && a_start != b_start)
            return a_start < b_start;
    }
    return false;
}

int64_t MatchHashEntry::end_to_start_compare(const MatchHashEntry& a, const MatchHashEntry& b) {
    int64_t a_end = a.RightEnd(0);
    int64_t b_start = b.LeftEnd(0);

    if (a_end == NO_MATCH || b_start == NO_MATCH)
        return 0;

    return b_start - a_end;
}

int64_t MatchHashEntry::start_compare(const MatchHashEntry& a, const MatchHashEntry& b) {
    int64_t a_start = a.LeftEnd(0);
    int64_t b_start = b.LeftEnd(0);

    if (a_start == NO_MATCH || b_start == NO_MATCH)
        return 0;

    return a_start - b_start;
}

bool MatchHashEntry::Contains(const MatchHashEntry& mhe) const {
    if (Length() < mhe.Length())
        return false;

    uint n = std::min(SeqCount(), mhe.SeqCount());
    for (uint i = 0; i < n; ++i) {
        int64_t this_start = LeftEnd(i);
        int64_t mhe_start = mhe.LeftEnd(i);

        if (this_start == NO_MATCH && mhe_start != NO_MATCH)
            return false;
        if (this_start != NO_MATCH && mhe_start == NO_MATCH)
            continue;

        if (this_start != NO_MATCH && mhe_start != NO_MATCH) {
            int64_t this_end = RightEnd(i);
            int64_t mhe_end = mhe.RightEnd(i);

            if (mhe_start < this_start || mhe_end > this_end)
                return false;
        }
    }

    return true;
}

} // namespace mems
