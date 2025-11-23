/////////////////////////////////////////////////////////////////////////////
// File:            gnStringSpec.cpp
// Purpose:         implements gnContigSpec for strings
// Description:     stores sequence in memory
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnStringSpec.h"
#include <string>
#include <cstring>
#include <algorithm>

namespace genome {

gnStringSpec::gnStringSpec() {
    gnContigSpec::Clear();
}

gnStringSpec::gnStringSpec(const std::string& m_string, gnSeqI start, gnSeqI endI, boolean revComp) {
    m_seqString = m_string;
    m_start = start;
    gnSeqI actual_len = static_cast<gnSeqI>(m_seqString.length());

    // reverse comp: end bp first, then start; switch them
    m_start = revComp ? endI : start;
    gnSeqI actual_end = revComp ? start : endI;

    // trim start and end if too big
    actual_end = std::min(actual_end, actual_len ? actual_len - 1 : 0);
    m_start   = std::min(m_start, actual_len ? actual_len - 1 : 0);
    if (actual_len == 0)
        m_start = 0;

    // if start is after end, use as circular
    m_circular = (m_start > actual_end);

    // if circular, length is different
    m_length = m_circular ? (actual_len - m_start) + actual_end : actual_end - m_start + 1;

    m_reverseComplement = revComp;
    m_SourceContigIndex = ALL_CONTIGS;
}

gnStringSpec::gnStringSpec(const gnStringSpec& s)
    : gnContigSpec(s)
{
    m_seqString          = s.m_seqString;
    m_sourceName         = s.m_sourceName;
    m_name               = s.m_name;
    m_start              = s.m_start;
    m_length             = s.m_length;
    m_reverseComplement  = s.m_reverseComplement;
    m_circular           = s.m_circular;
    m_SourceContigIndex  = s.m_SourceContigIndex;
}

gnStringSpec::~gnStringSpec() {
    Clear();
}

void gnStringSpec::Clear() {
    gnContigSpec::Clear();
    m_seqString.clear();
}

gnStringSpec* gnStringSpec::CloneRange(const gnSeqI startI, const gnSeqI len) const {
    auto* destSpec = new gnStringSpec();
    destSpec->m_seqString = m_seqString.substr(m_start + startI, len);
    destSpec->m_sourceName = m_sourceName;
    destSpec->m_name = m_name;
    destSpec->m_start = 0;
    destSpec->m_length = static_cast<gnSeqI>(destSpec->m_seqString.length());
    destSpec->m_reverseComplement = m_reverseComplement;
    destSpec->m_circular = m_circular;
    destSpec->m_SourceContigIndex = m_SourceContigIndex;
    return destSpec;
}

} // end namespace genome
