/////////////////////////////////////////////////////////////////////////////
// File:            gnFileContig.cpp
// Purpose:         File Position holder.
// Description:     Implements gnFileContig methods.
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFileContig.h"

namespace genome {

gnFileContig::gnFileContig() noexcept {
    Clear();
}

gnFileContig::gnFileContig(std::string nameStr, const uint64 s, const uint64 e) noexcept
    : m_name(std::move(nameStr)), m_seqLength(0), m_fileStartEnd{s, e},
      m_repeatSeqGap(false), m_repeatSeqGapSize{0, 0}
{
    for (uint32 i = 0; i < CONTIG_SECTION_SIZE; ++i)
        m_startEndArray[i] = {0, 0};
}

void gnFileContig::Clear() noexcept {
    m_name.clear();
    m_seqLength = 0;
    m_fileStartEnd = {0, 0};
    for (uint32 i = 0; i < CONTIG_SECTION_SIZE; ++i)
        m_startEndArray[i] = {0, 0};
    m_repeatSeqGap = false;
    m_repeatSeqGapSize = {0, 0};
}

boolean gnFileContig::SetRepeatSeqSize(const uint64 seqSize) noexcept {
    if (!m_repeatSeqGap)
        return false;
    if (m_repeatSeqGapSize.first == seqSize)
        return true;
    if (m_repeatSeqGapSize.first == 0) {
        m_repeatSeqGapSize.first = seqSize;
        return true;
    }
    m_repeatSeqGap = false;
    return false;
}

boolean gnFileContig::SetRepeatGapSize(const uint64 gapSize) noexcept {
    if (!m_repeatSeqGap)
        return false;
    if (m_repeatSeqGapSize.second == gapSize)
        return true;
    if (m_repeatSeqGapSize.second == 0) {
        m_repeatSeqGapSize.second = gapSize;
        return true;
    }
    m_repeatSeqGap = false;
    return false;
}

} // namespace genome
