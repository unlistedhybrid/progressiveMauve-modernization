/////////////////////////////////////////////////////////////////////////////
// File:            gnFilter.h
// Purpose:         Filter for all Sequences
// Description:     Filters sequences, translates, reverse complement, converts
//                   additions, etc.
// Author:          Aaron Darling
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnDefs.h"
#include "libGenome/gnBaseFilter.h"

namespace genome {

constexpr gnSeqC NO_REVCOMP_CHAR = 0;

class GNDLLEXPORT gnFilter : public gnBaseFilter
{
public:
    static const gnFilter* alphabetCharacterFilter();
    static const gnFilter* numberCharacterFilter();

    static const gnFilter* proteinSeqFilter();
    static const gnFilter* basicDNASeqFilter();
    static const gnFilter* fullDNASeqFilter();
    static const gnFilter* basicRNASeqFilter();
    static const gnFilter* fullRNASeqFilter();

    static const gnFilter* DNAtoRNAFilter();
    static const gnFilter* RNAtoDNAFilter();
    static const gnFilter* DNAComplementFilter();
    static const gnFilter* RNAComplementFilter();

    enum gnFilterType {
        alphabetCharacterFilterType,
        numberCharacterFilterType,
        proteinSeqFilterType,
        basicDNASeqFilterType,
        fullDNASeqFilterType,
        basicRNASeqFilterType,
        fullRNASeqFilterType,
        DNAtoRNAFilterType,
        RNAtoDNAFilterType,
        DNAComplementFilterType,
        RNAComplementFilterType
    };

public:
    gnFilter();
    /**
     * Creates a gnFilter for a predefined filter type.
     * Used by the static sequence constructors to avoid the "static initialization order fiasco"
     * @param f_type The type of filter to create.
     * @see gnFilterType
     */
    gnFilter(const gnFilterType f_type);
    gnFilter(const gnSeqC defaultChar, const gnSeqC rdefaultChar);
    gnFilter(const gnFilter& sf);
    ~gnFilter();

    gnFilter* Clone() const;

    // gnSeqC 
    bool IsValid(const gnSeqC ch) const;
    gnSeqC MakeValid(const gnSeqC ch) const;
    gnSeqC Filter(const gnSeqC ch) const;

    // gnSeqC[]
    /**
     * IsValid() scans the supplied character array for invalid characters.
     * @param seq The sequence to scan.  This is a generic character array.
     * @param len The length of the sequence to scan.
     * @return The index of the first invalid character, or len if none exists
     */
    uint32 IsValid(const gnSeqC* seq, uint32 len) const;
    void MakeValid(gnSeqC* seq, uint32 len) const;
    void Filter(gnSeqC** seq, gnSeqI& len) const;
    void ReverseFilter(gnSeqC** seq, gnSeqI& len) const;

    // std::string
    uint32 IsValid(const std::string& seq) const;
    void MakeValid(std::string& seq) const;
    void Filter(std::string& seq) const;
    void ReverseFilter(std::string& seq) const;

    // Default gnSeqC
    void SetDefaultChar(const gnSeqC ch1, const gnSeqC ch2);
    gnSeqC GetDefaultChar() const;
    gnSeqC GetRDefaultChar() const;

    // fill map
    void SetSingle(const gnSeqC ch);
    void SetPair(const gnSeqC ch1, const gnSeqC ch2);
    bool RemovePair(const gnSeqC ch);
    bool RemoveSingle(const gnSeqC ch);

private:
    void CreateAlphabetCharacterFilter();
    void CreateNumberCharacterFilter();

    void CreateProteinFilter();

    void CreateBasicDNAFilter();
    void CreateFullDNAFilter();

    void CreateBasicRNAFilter();
    void CreateFullRNAFilter();

    void CreateDNAtoRNAFilter();
    void CreateRNAtoDNAFilter();
    void CreateDNAComplementFilter();
    void CreateRNAComplementFilter();

    gnSeqC m_pairArray[GNSEQC_MAX];
    gnSeqC m_defaultChar;
    gnSeqC m_rDefaultChar;
};

// === Inline method implementations ===

inline gnFilter* gnFilter::Clone() const {
    return new gnFilter(*this);
}

// gnSeqC
inline bool gnFilter::IsValid(const gnSeqC ch) const {
    return m_pairArray[static_cast<unsigned char>(ch)] != NO_REVCOMP_CHAR;
}
inline gnSeqC gnFilter::MakeValid(const gnSeqC ch) const {
    return (m_pairArray[static_cast<unsigned char>(ch)] != NO_REVCOMP_CHAR ? ch : m_defaultChar);
}
inline gnSeqC gnFilter::Filter(const gnSeqC ch) const {
    return m_pairArray[static_cast<unsigned char>(ch)] != NO_REVCOMP_CHAR
        ? m_pairArray[static_cast<unsigned char>(ch)]
        : m_defaultChar;
}

// gnSeqC[]
inline uint32 gnFilter::IsValid(const gnSeqC* seq, uint32 len) const {
    for (uint32 i = 0; i < len; ++i) {
        if (!IsValid(seq[i]))
            return i;
    }
    return len;
}
inline void gnFilter::MakeValid(gnSeqC* seq, uint32 len) const {
    for (uint32 i = 0; i < len; ++i) {
        seq[i] = MakeValid(seq[i]);
    }
}

// std::string
inline uint32 gnFilter::IsValid(const std::string& seq) const {
    return IsValid(reinterpret_cast<const gnSeqC*>(seq.data()), static_cast<uint32>(seq.length()));
}
inline void gnFilter::MakeValid(std::string& seq) const {
    MakeValid(reinterpret_cast<gnSeqC*>(const_cast<char*>(seq.data())), static_cast<uint32>(seq.length()));
}

// Default gnSeqC
inline void gnFilter::SetDefaultChar(const gnSeqC ch1, const gnSeqC ch2) {
    m_defaultChar = ch1;
    m_rDefaultChar = ch2;
}
inline gnSeqC gnFilter::GetDefaultChar() const {
    return m_defaultChar;
}
inline gnSeqC gnFilter::GetRDefaultChar() const {
    return m_rDefaultChar;
}
// fill map
inline void gnFilter::SetSingle(const gnSeqC ch) {
    m_pairArray[static_cast<unsigned char>(ch)] = ch;
}
inline void gnFilter::SetPair(const gnSeqC ch1, const gnSeqC ch2) {
    m_pairArray[static_cast<unsigned char>(ch1)] = ch2;
}
inline bool gnFilter::RemovePair(const gnSeqC ch) {
    m_pairArray[static_cast<unsigned char>(ch)] = NO_REVCOMP_CHAR;
    return true;
}
inline bool gnFilter::RemoveSingle(const gnSeqC ch) {
    m_pairArray[static_cast<unsigned char>(ch)] = NO_REVCOMP_CHAR;
    return true;
}

} // end namespace genome
