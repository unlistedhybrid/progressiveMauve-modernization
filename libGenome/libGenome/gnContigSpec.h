/////////////////////////////////////////////////////////////////////////////
// File:            gnContigSpec.h
// Purpose:         Abstract Contig Spec class
// Description:     Defines an interface for contig specs
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnContigSpec_h_
#define _gnContigSpec_h_

#include "libGenome/gnDefs.h"

#include <vector>
#include <string>
#include <utility>  // For std::move if needed

#include "libGenome/gnBaseSpec.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnDebug.h"

namespace genome {

/**
 * gnContigSpec is an interface for classes which store contigs, or reads,
 * of DNA or protein sequence
 */
class GNDLLEXPORT gnContigSpec : public gnBaseSpec {
public:
    gnContigSpec() = default;
    virtual ~gnContigSpec() = default;

    // Rule of 5: Default copy/move, can be deleted or defaulted as needed
    gnContigSpec(const gnContigSpec&) = default;
    gnContigSpec(gnContigSpec&&) noexcept = default;
    gnContigSpec& operator=(const gnContigSpec&) = default;
    gnContigSpec& operator=(gnContigSpec&&) noexcept = default;

    // Polymorphic interface
    [[nodiscard]] virtual gnContigSpec* Clone() const = 0;
    [[nodiscard]] virtual gnContigSpec* CloneRange(gnSeqI startI, gnSeqI len) const = 0;

    [[nodiscard]] virtual std::string GetSourceName() const;
    [[nodiscard]] virtual gnSeqI GetStart() const noexcept;
    [[nodiscard]] virtual gnSeqI GetLength() const noexcept;
    [[nodiscard]] virtual gnSeqI GetSourceLength() const = 0;
    [[nodiscard]] virtual uint32 GetSourceContigIndex() const noexcept;

    virtual void SetSourceName(const std::string& sourceName);
    virtual void SetStart(gnSeqI start) noexcept;
    virtual void SetLength(gnSeqI len) noexcept;
    virtual void SetSourceContigIndex(uint32 contigI) noexcept;
    virtual void SetReverseComplement(boolean value);

    virtual void CropStart(gnSeqI cropLen);
    virtual void CropEnd(gnSeqI cropLen);

    virtual boolean SeqRead(gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, uint32 contigI) const;
    virtual void Clear();

protected:
    gnSeqI m_start{};                // start within the genome.
    gnSeqI m_length{};
    uint32 m_SourceContigIndex{};
    std::string m_sourceName;        // Source name for this contig.

    /**
     * All derived classes must implement this!
     * Reads the specified bases into buf, disregarding circularity and reverse complement.
     */
    virtual boolean Read(gnSeqI start, gnSeqC* buf, gnSeqI& bufLen) const = 0;

private:
    // (No private fields or methods)
};

// Inline methods

inline std::string gnContigSpec::GetSourceName() const {
    return m_sourceName;
}
inline void gnContigSpec::SetSourceName(const std::string& sourceName) {
    m_sourceName = sourceName;
}

inline gnSeqI gnContigSpec::GetStart() const noexcept {
    return m_start;
}
inline gnSeqI gnContigSpec::GetLength() const noexcept {
    return m_length;
}
inline void gnContigSpec::SetStart(gnSeqI start) noexcept {
    m_start = start;
}
inline void gnContigSpec::SetLength(gnSeqI len) noexcept {
    m_length = len;
}
inline uint32 gnContigSpec::GetSourceContigIndex() const noexcept {
    return m_SourceContigIndex;
}
inline void gnContigSpec::SetSourceContigIndex(uint32 contigI) noexcept {
    m_SourceContigIndex = contigI;
}

}   // namespace genome

#endif  // _gnContigSpec_h_
