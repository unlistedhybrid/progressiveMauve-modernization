/////////////////////////////////////////////////////////////////////////////
// File:            gnFilter.cpp
// Purpose:         Filter for all Sequences
// Description:     Filters sequences, translates, reverse complement, converts
//                   additions, etc.
// Author:          Aaron Darling
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFilter.h"
#include "libGenome/gnDebug.h"
#include <cstring>

namespace genome {

// Static filter singletons
const gnFilter* gnFilter::alphabetCharacterFilter() {
    static const gnFilter* t_filt = new gnFilter(alphabetCharacterFilterType);
    return t_filt;
}

const gnFilter* gnFilter::numberCharacterFilter() {
    static const gnFilter* t_filt = new gnFilter(numberCharacterFilterType);
    return t_filt;
}

const gnFilter* gnFilter::proteinSeqFilter() {
    static const gnFilter* t_filt = new gnFilter(proteinSeqFilterType);
    return t_filt;
}

const gnFilter* gnFilter::basicDNASeqFilter() {
    static const gnFilter* t_filt = new gnFilter(basicDNASeqFilterType);
    return t_filt;
}

const gnFilter* gnFilter::fullDNASeqFilter() {
    static const gnFilter* t_filt = new gnFilter(fullDNASeqFilterType);
    return t_filt;
}

const gnFilter* gnFilter::basicRNASeqFilter() {
    static const gnFilter* t_filt = new gnFilter(basicRNASeqFilterType);
    return t_filt;
}

const gnFilter* gnFilter::fullRNASeqFilter() {
    static const gnFilter* t_filt = new gnFilter(fullRNASeqFilterType);
    return t_filt;
}

const gnFilter* gnFilter::DNAtoRNAFilter() {
    static const gnFilter* t_filt = new gnFilter(DNAtoRNAFilterType);
    return t_filt;
}

const gnFilter* gnFilter::RNAtoDNAFilter() {
    static const gnFilter* t_filt = new gnFilter(RNAtoDNAFilterType);
    return t_filt;
}

const gnFilter* gnFilter::DNAComplementFilter() {
    static const gnFilter* t_filt = new gnFilter(DNAComplementFilterType);
    return t_filt;
}

const gnFilter* gnFilter::RNAComplementFilter() {
    static const gnFilter* t_filt = new gnFilter(RNAComplementFilterType);
    return t_filt;
}

// Constructors & destructor
gnFilter::gnFilter() {
    m_defaultChar = 'n';
    m_rDefaultChar = 'n';
    for (int i = 0; i < GNSEQC_MAX; ++i)
        m_pairArray[static_cast<unsigned char>(i)] = NO_REVCOMP_CHAR;
}

gnFilter::gnFilter(const gnSeqC defaultChar, const gnSeqC rdefaultChar) {
    m_defaultChar = defaultChar;
    m_rDefaultChar = rdefaultChar;
    for (int i = 0; i < GNSEQC_MAX; ++i)
        m_pairArray[static_cast<unsigned char>(i)] = NO_REVCOMP_CHAR;
}

gnFilter::gnFilter(const gnFilter& sf) {
    m_name = sf.m_name;
    for (int i = 0; i < GNSEQC_MAX; ++i)
        m_pairArray[static_cast<unsigned char>(i)] = sf.m_pairArray[static_cast<unsigned char>(i)];
    m_defaultChar = sf.m_defaultChar;
    m_rDefaultChar = sf.m_rDefaultChar;
}

gnFilter::gnFilter(const gnFilterType f_type) {
    for (int i = 0; i < GNSEQC_MAX; ++i)
        m_pairArray[static_cast<unsigned char>(i)] = NO_REVCOMP_CHAR;
    switch (f_type) {
        case alphabetCharacterFilterType:  CreateAlphabetCharacterFilter(); break;
        case numberCharacterFilterType:    CreateNumberCharacterFilter(); break;
        case proteinSeqFilterType:         CreateProteinFilter(); break;
        case basicDNASeqFilterType:        CreateBasicDNAFilter(); break;
        case fullDNASeqFilterType:         CreateFullDNAFilter(); break;
        case basicRNASeqFilterType:        CreateBasicRNAFilter(); break;
        case fullRNASeqFilterType:         CreateFullRNAFilter(); break;
        case DNAtoRNAFilterType:           CreateDNAtoRNAFilter(); break;
        case RNAtoDNAFilterType:           CreateRNAtoDNAFilter(); break;
        case DNAComplementFilterType:      CreateDNAComplementFilter(); break;
        case RNAComplementFilterType:      CreateRNAComplementFilter(); break;
    }
}

gnFilter::~gnFilter() = default;

// Filtering and reverse filtering on raw sequence arrays
void gnFilter::Filter(gnSeqC** seq, gnSeqI& len) const {
    Array<gnSeqC> array_buf(len);
    gnSeqC* tmp = array_buf.data;
    gnSeqI c = 0;
    for (uint32 i = 0; i < len; ++i)
        if (IsValid((*seq)[i]))
            tmp[c++] = m_pairArray[static_cast<unsigned char>((*seq)[i])];
    len = c;
    std::memcpy(*seq, tmp, len);
}

void gnFilter::ReverseFilter(gnSeqC** seq, gnSeqI& len) const {
    gnSeqC tmp, dum;
    uint32 halfLen = len / 2;
    uint32 end = len - 1;
    uint32 curB = 0;
    uint32 curE = end;
    for (uint32 i = 0; i < halfLen; ++i) {
        tmp = m_pairArray[static_cast<unsigned char>((*seq)[i])];
        dum = m_pairArray[static_cast<unsigned char>((*seq)[end - i])];
        if (dum != NO_REVCOMP_CHAR)
            (*seq)[curB++] = dum;
        if (tmp != NO_REVCOMP_CHAR)
            (*seq)[curE--] = tmp;
    }
    if (len & 0x1) {
        tmp = m_pairArray[static_cast<unsigned char>((*seq)[halfLen])];
        if (tmp != NO_REVCOMP_CHAR)
            (*seq)[curB++] = tmp;
    }
    if (curE >= curB) {
        std::memmove(*seq + curB, *seq + curE + 1, end - curE);
        len = end - curE + curB;
    }
}

// Filtering and reverse filtering on std::string
void gnFilter::Filter(std::string& seq) const {
    gnSeqI c = 0;
    for (uint32 i = 0; i < seq.length(); ++i)
        if (IsValid(seq[i]))
            seq[c++] = m_pairArray[static_cast<unsigned char>(seq[i])];
    seq.resize(c);
}

void gnFilter::ReverseFilter(std::string& seq) const {
    gnSeqC tmp, dum;
    uint32 halfLen = seq.length() / 2;
    uint32 end = seq.length() - 1;
    uint32 curB = 0;
    uint32 curE = end;
    for (uint32 i = 0; i < halfLen; ++i) {
        tmp = m_pairArray[static_cast<unsigned char>(seq[i])];
        dum = m_pairArray[static_cast<unsigned char>(seq[end - i])];
        if (dum != NO_REVCOMP_CHAR)
            seq[curB++] = dum;
        if (tmp != NO_REVCOMP_CHAR)
            seq[curE--] = tmp;
    }
    if (seq.length() & 0x1) {
        tmp = m_pairArray[static_cast<unsigned char>(seq[halfLen])];
        if (tmp != NO_REVCOMP_CHAR)
            seq[curB++] = tmp;
    }
    if (curE >= curB) {
        seq.erase(curB, curE - curB);
    }
}

// Filter creator implementations
void gnFilter::CreateAlphabetCharacterFilter() {
    SetDefaultChar(0, 0);
    SetName("Alphabet Character Filter");
    for (char c = 'A'; c <= 'Z'; ++c) SetPair(c, std::tolower(c));
    for (char c = 'a'; c <= 'z'; ++c) SetPair(c, c);
}

void gnFilter::CreateNumberCharacterFilter() {
    SetDefaultChar(0, 0);
    SetName("Number Character Filter");
    for (char c = '0'; c <= '9'; ++c) SetSingle(c);
}

void gnFilter::CreateProteinFilter() {
    SetDefaultChar('u', 'u');
    SetName("Protein Filter");
    const char* aminoAcids = "ARNDCQEGHILKMFPSTWYV";
    for (const char* p = aminoAcids; *p; ++p) SetSingle(*p);
    for (const char* p = aminoAcids; *p; ++p) SetSingle(std::tolower(*p));
}

void gnFilter::CreateBasicDNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Basic DNA Filter");
    const char chars[] = "acgtACGTnNxX-";
    for (const char* p = chars; *p; ++p) SetSingle(*p);
}

void gnFilter::CreateFullDNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full DNA Filter");
    const char chars[] = "acgtACGTrykmbvdhRYKMBVDHsSwWnNxX-";
    for (const char* p = chars; *p; ++p) SetSingle(*p);
}

void gnFilter::CreateBasicRNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Basic RNA Filter");
    const char chars[] = "acguACGUnN-";
    for (const char* p = chars; *p; ++p) SetSingle(*p);
}

void gnFilter::CreateFullRNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full RNA Filter");
    const char chars[] = "acguACGUrykmbvdhRYKMBVDHsSwWnN-";
    for (const char* p = chars; *p; ++p) SetSingle(*p);
}

void gnFilter::CreateDNAtoRNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full DNA to RNA Filter");
    const char unchanged[] = "acgACGrkybvdhRKYBVDHsSwWnN-";
    for (const char* p = unchanged; *p; ++p) SetSingle(*p);
    SetPair('t', 'u'); SetPair('T', 'U');
}

void gnFilter::CreateRNAtoDNAFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full RNA to DNA Filter");
    const char unchanged[] = "acgACGrkybvdhRKYBVDHsSwWnN-";
    for (const char* p = unchanged; *p; ++p) SetSingle(*p);
    SetPair('u', 't'); SetPair('U', 'T');
}

void gnFilter::CreateDNAComplementFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full DNA Complement Filter");
    SetPair('a', 't'); SetPair('A', 'T');
    SetPair('t', 'a'); SetPair('T', 'A');
    SetPair('c', 'g'); SetPair('C', 'G');
    SetPair('g', 'c'); SetPair('G', 'C');
    SetPair('r', 'y'); SetPair('R', 'Y');
    SetPair('y', 'r'); SetPair('Y', 'R');
    SetPair('k', 'm'); SetPair('K', 'M');
    SetPair('m', 'k'); SetPair('M', 'K');
    SetSingle('s'); SetSingle('S');
    SetSingle('w'); SetSingle('W');
    SetPair('b', 'v'); SetPair('B', 'V');
    SetPair('v', 'b'); SetPair('V', 'B');
    SetPair('d', 'h'); SetPair('D', 'H');
    SetPair('h', 'd'); SetPair('H', 'D');
    SetSingle('n'); SetSingle('N'); SetSingle('x'); SetSingle('X'); SetSingle('-');
}

void gnFilter::CreateRNAComplementFilter() {
    SetDefaultChar('n', 'n');
    SetName("Full RNA Complement Filter");
    SetPair('a', 'u'); SetPair('A', 'U');
    SetPair('u', 'a'); SetPair('U', 'A');
    SetPair('c', 'g'); SetPair('C', 'G');
    SetPair('g', 'c'); SetPair('G', 'C');
    SetPair('r', 'y'); SetPair('R', 'Y');
    SetPair('y', 'r'); SetPair('Y', 'R');
    SetPair('k', 'm'); SetPair('K', 'M');
    SetPair('m', 'k'); SetPair('M', 'K');
    SetSingle('s'); SetSingle('S');
    SetSingle('w'); SetSingle('W');
    SetPair('b', 'v'); SetPair('B', 'V');
    SetPair('v', 'b'); SetPair('V', 'B');
    SetPair('d', 'h'); SetPair('D', 'H');
    SetPair('h', 'd'); SetPair('H', 'D');
    SetSingle('n'); SetSingle('N'); SetSingle('-');
}

} // namespace genome
