/////////////////////////////////////////////////////////////////////////////
// File:            gnCompare.h
// Purpose:         Compares all sequences
// Description:     Compares sequences (protein, DNA, RNA)
// Author:          Aaron Darling 
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"

#include <string>
#include <cstring>
#include <cctype>
#include "libGenome/gnClone.h"

namespace genome {

class GNDLLEXPORT gnCompare : public gnClone
{
public:
    // Static accessors (singleton-style)
    static const gnCompare* ProteinSeqCompare();
    static const gnCompare* DNASeqCompare();
    static const gnCompare* RNASeqCompare();

    enum gnCompareType {
        ProteinSeqCompareType,
        DNASeqCompareType,
        RNASeqCompareType
    };

    gnCompare();
    gnCompare(gnCompareType c_type);
    ~gnCompare();

    gnCompare* Clone() const;

    const std::string& GetName() const;
    void SetName(const std::string& name);

    // Comparison operations
    bool LessThan(gnSeqC ch, gnSeqC ch2, bool case_sensitive = false) const;
    bool Contains(gnSeqC ch, gnSeqC ch2, bool case_sensitive = false) const;

    bool LessThan(const gnSeqC* seq, const gnSeqC* seq2,
                  uint32 len, bool case_sensitive = false) const;

    bool Contains(const gnSeqC* seq, const gnSeqC* seq2,
                  uint32 len, bool case_sensitive = false) const;

    bool LessThan(const std::string& seq, const std::string& seq2,
                  bool case_sensitive = false) const;

    bool Contains(const std::string& seq, const std::string& seq2,
                  bool case_sensitive = false) const;

    // Fill map
    void SetSingle(gnSeqC ch);
    void SetPair(gnSeqC ch, gnSeqC ch2);
    void SetContained(gnSeqC ch, gnSeqC ch2);

    void RemoveSingle(gnSeqC ch);
    void RemovePair(gnSeqC ch, gnSeqC ch2);
    void RemoveContained(gnSeqC ch, gnSeqC ch2);

private:
    // Private copy constructor to control cloning
    gnCompare(const gnCompare&) {}

    void CreateProteinComparator();
    void CreateDNAComparator();
    void CreateRNAComparator();

    void AddArrayEntry(gnSeqC* array[GNSEQC_MAX], gnSeqC ch, gnSeqC ch2);
    void DelArrayEntry(gnSeqC* array[GNSEQC_MAX], gnSeqC ch, gnSeqC ch2);

    std::string m_name;
    bool m_ignoreCase = false;

    // Arrays of dynamically-built char arrays
    gnSeqC* m_pairArray[GNSEQC_MAX];
    gnSeqC* m_containArray[GNSEQC_MAX];
};

// -----------------------------------------------------------------------------
// Inline implementations
// -----------------------------------------------------------------------------

inline gnCompare* gnCompare::Clone() const {
    return new gnCompare(*this);
}

inline const std::string& gnCompare::GetName() const {
    return m_name;
}

inline void gnCompare::SetName(const std::string& name) {
    m_name = name;
}

inline bool gnCompare::LessThan(gnSeqC ch, gnSeqC ch2, bool case_sensitive) const {
    if (!case_sensitive) {
        ch = static_cast<gnSeqC>(std::toupper(static_cast<unsigned char>(ch)));
        ch2 = static_cast<gnSeqC>(std::toupper(static_cast<unsigned char>(ch2)));
    }

    const auto idx = static_cast<unsigned char>(ch);
    const gnSeqC* arr = m_pairArray[idx];

    if (arr && std::strchr(arr, ch2) == nullptr) {
        return ch < ch2;
    }

    return false;
}

inline bool gnCompare::LessThan(const gnSeqC* seq, const gnSeqC* seq2,
                                uint32 len, bool case_sensitive) const {
    for (uint32 i = 0; i < len; ++i) {
        if (LessThan(seq[i], seq2[i], case_sensitive))
            return true;
    }
    return false;
}

inline bool gnCompare::LessThan(const std::string& seq, const std::string& seq2,
                                bool case_sensitive) const {
    const uint32 shorter = static_cast<uint32>(std::min(seq.size(), seq2.size()));
    return LessThan(reinterpret_cast<const gnSeqC*>(seq.data()),
                    reinterpret_cast<const gnSeqC*>(seq2.data()),
                    shorter, case_sensitive);
}

inline void gnCompare::SetSingle(gnSeqC ch) {
    AddArrayEntry(m_pairArray, ch, ch);
    AddArrayEntry(m_containArray, ch, ch);
}

inline void gnCompare::SetPair(gnSeqC ch, gnSeqC ch2) {
    AddArrayEntry(m_pairArray, ch, ch2);
    AddArrayEntry(m_pairArray, ch2, ch);
}

inline void gnCompare::SetContained(gnSeqC ch, gnSeqC ch2) {
    AddArrayEntry(m_containArray, ch2, ch);
}

inline void gnCompare::RemoveSingle(gnSeqC ch) {
    DelArrayEntry(m_pairArray, ch, ch);
    DelArrayEntry(m_containArray, ch, ch);
}

inline void gnCompare::RemovePair(gnSeqC ch, gnSeqC ch2) {
    DelArrayEntry(m_pairArray, ch, ch2);
    DelArrayEntry(m_pairArray, ch2, ch);
}

inline void gnCompare::RemoveContained(gnSeqC ch, gnSeqC ch2) {
    DelArrayEntry(m_containArray, ch2, ch);
}

} // namespace genome
