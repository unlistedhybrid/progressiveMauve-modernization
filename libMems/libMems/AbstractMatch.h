/*******************************************************************************
 * $Id: AbstractMatch.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 ******************************************************************************/

#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <vector>
#include <algorithm>
#include <type_traits>
#include <boost/dynamic_bitset.hpp>
#include <libMems/SlotAllocator.h>
#include <libMems/configuration.h>

namespace mems {

static const gnSeqI NO_MATCH = 0;

using bitset_t = boost::dynamic_bitset<>;

// SlotAllocator helpers (still needed for memory pool logic)
template<typename T>
T* m_allocateAndCopy(const T& t) {
    SlotAllocator<T>& sat = SlotAllocator<T>::GetSlotAllocator();
    T* newt = sat.Allocate();
    new (newt) T(t);
    return newt;
}

template<typename T>
void m_free(T* t) {
    SlotAllocator<T>& sat = SlotAllocator<T>::GetSlotAllocator();
    sat.Free(t);
}

/**
 * AbstractMatch is a pure virtual base class that defines an interface for 
 * both gapped and ungapped alignments among several sequences or several regions
 * of the same sequence 
 */
class AbstractMatch : public genome::gnClone {
public:
    enum orientation {
        forward,    // alignment on the forward strand
        reverse,    // alignment on the reverse strand
        undefined   // no alignment on either strand
    };

    virtual ~AbstractMatch() = default;

    // Memory pool API
    virtual AbstractMatch* Copy() const = 0;
    virtual void Free() = 0;

    // Coordinate and length queries (pure virtual)
    virtual gnSeqI Length(uint seqI) const = 0;
    virtual void SetLength(gnSeqI len, uint seqI) = 0;
    virtual int64 Start(uint startI) const = 0;
    virtual void SetStart(uint seqI, int64 start) = 0;

    // Convenient synonym for Start()
    int64 operator[](uint seqI) const { return Start(seqI); }

    // Deprecated: Use LeftEnd() instead
    virtual int64 End(uint seqI) const;

    // Alignment endpoints and orientation
    virtual gnSeqI LeftEnd(uint seqI) const = 0;
    virtual gnSeqI RightEnd(uint seqI) const { return LeftEnd(seqI) + Length(seqI) - 1; }
    virtual orientation Orientation(uint seqI) const = 0;

    virtual void SetLeftEnd(uint seqI, gnSeqI start) = 0;
    virtual void SetOrientation(uint seqI, orientation o) = 0;

    virtual void MoveStart(int64 move_amount) = 0;
    virtual void MoveEnd(int64 move_amount) = 0;

    // Sequence coverage/multiplicity
    virtual uint Multiplicity() const = 0;
    virtual uint SeqCount() const = 0;
    virtual uint FirstStart() const = 0;
    virtual gnSeqI AlignmentLength() const = 0;

    virtual void Invert() = 0;

    // Cropping APIs (some deprecated)
    virtual void CropStart(gnSeqI crop_amount) = 0;
    virtual void CropEnd(gnSeqI crop_amount) = 0;
    virtual void CropLeft(gnSeqI crop_amount, uint seqI) = 0;
    virtual void CropRight(gnSeqI crop_amount, uint seqI) = 0;

    // Alignment structure
    virtual void GetAlignment(std::vector<bitset_t>& align_matrix) const = 0;
    virtual void GetColumn(gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column) const = 0;

    virtual bool IsGap(uint seq, gnSeqI col) const = 0;
    virtual uint UsedSeq(uint seqI) const = 0;
};

// Inline for End()
inline int64 AbstractMatch::End(uint endI) const {
    if (Start(endI) > 0)
        return Start(endI) + Length(endI) - 1;
    return Start(endI);
}

// Comparator functors (modernized, bool instead of boolean)

template<typename MatchType>
class AbstractMatchStartComparator {
public:
    explicit AbstractMatchStartComparator(unsigned seq = 0) : m_seq(seq) {}
    AbstractMatchStartComparator(const AbstractMatchStartComparator& msc) = default;
    AbstractMatchStartComparator& operator=(const AbstractMatchStartComparator& msc) = default;

    bool operator()(const MatchType& a, const MatchType& b) const {
        int start_diff = std::max(a.FirstStart(), m_seq) - std::max(b.FirstStart(), m_seq);
        if (start_diff == 0) {
            uint m_count = std::min(a.SeqCount(), b.SeqCount());
            for (uint seqI = m_seq; seqI < m_count; seqI++) {
                gnSeqI a_start = a.Orientation(seqI) == AbstractMatch::forward ? a.LeftEnd(seqI) : a.RightEnd(seqI);
                gnSeqI b_start = b.Orientation(seqI) == AbstractMatch::forward ? b.LeftEnd(seqI) : b.RightEnd(seqI);
                if (a_start == NO_MATCH || b_start == NO_MATCH)
                    continue;
                else if (a_start == b_start)
                    continue;
                else
                    return a_start < b_start;
            }
        }
        return start_diff < 0;
    }
private:
    unsigned m_seq;
};

template<typename MatchType>
class AbstractMatchSingleStartComparator {
public:
    explicit AbstractMatchSingleStartComparator(unsigned seq = 0) : m_seq(seq) {}
    AbstractMatchSingleStartComparator(const AbstractMatchSingleStartComparator& msc) = default;
    AbstractMatchSingleStartComparator& operator=(const AbstractMatchSingleStartComparator& msc) = default;

    bool operator()(const MatchType& a, const MatchType& b) const {
        int64 a_start = a.LeftEnd(m_seq), b_start = b.LeftEnd(m_seq);
        if (a_start == NO_MATCH || b_start == NO_MATCH)
            return b_start != NO_MATCH;
        return a_start < b_start;
    }
private:
    unsigned m_seq;
};

template<typename MatchType>
class MatchStartComparator {
public:
    explicit MatchStartComparator(unsigned seq = 0) : m_seq(seq) {}
    MatchStartComparator(const MatchStartComparator& msc) = default;
    MatchStartComparator& operator=(const MatchStartComparator& msc) = default;

    bool operator()(const MatchType* a, const MatchType* b) const {
        int start_diff = std::max(a->FirstStart(), m_seq) - std::max(b->FirstStart(), m_seq);
        if (start_diff == 0) {
            uint m_count = std::min(a->SeqCount(), b->SeqCount());
            for (uint seqI = m_seq; seqI < m_count; seqI++) {
                gnSeqI a_start = a->Orientation(seqI) == AbstractMatch::forward ? a->LeftEnd(seqI) : a->RightEnd(seqI);
                gnSeqI b_start = b->Orientation(seqI) == AbstractMatch::forward ? b->LeftEnd(seqI) : b->RightEnd(seqI);
                if (a_start == NO_MATCH || b_start == NO_MATCH)
                    continue;
                else if (a_start == b_start)
                    continue;
                else
                    return a_start < b_start;
            }
        }
        return start_diff < 0;
    }
private:
    unsigned m_seq;
};

template<typename MatchType>
class SingleStartComparator {
public:
    explicit SingleStartComparator(unsigned seq = 0) : m_seq(seq) {}
    SingleStartComparator(const SingleStartComparator& msc) = default;
    SingleStartComparator& operator=(const SingleStartComparator& msc) = default;

    bool operator()(const MatchType* a, const MatchType* b) const {
        int64 a_start = a->LeftEnd(m_seq), b_start = b->LeftEnd(m_seq);
        if (a_start == NO_MATCH || b_start == NO_MATCH)
            return b_start != NO_MATCH;
        return a_start < b_start;
    }
private:
    unsigned m_seq;
};

template<typename MatchType>
class SSC {
public:
    explicit SSC(unsigned seq = 0) : m_seq(seq) {}
    SSC(const SSC<MatchType>& msc) = default;
    SSC<MatchType>& operator=(const SSC<MatchType>& msc) = default;

    bool operator()(const typename std::add_pointer<MatchType>::type& a,
                    const typename std::add_pointer<MatchType>::type& b) const {
        return operator()(*a, *b);
    }
    bool operator()(const typename std::remove_pointer<MatchType>::type& a,
                    const typename std::remove_pointer<MatchType>::type& b) const {
        int64 a_start = a.LeftEnd(m_seq), b_start = b.LeftEnd(m_seq);
        if (a_start == NO_MATCH || b_start == NO_MATCH)
            return b_start != NO_MATCH;
        return a_start < b_start;
    }
private:
    unsigned m_seq;
};

} // namespace mems
