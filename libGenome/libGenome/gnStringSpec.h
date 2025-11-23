/////////////////////////////////////////////////////////////////////////////
// File:            gnStringSpec.h
// Purpose:         implements gnContigSpec for strings
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnStringSpec_h_
#define _gnStringSpec_h_

#include "libGenome/gnDefs.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnBaseSource.h"

namespace genome {

/**
 * gnStringSpec stores a sequence and annotation data in memory.
 * For a more complete description see the gnBaseSpec documentation.
 */
class GNDLLEXPORT gnStringSpec : public gnContigSpec
{
public:
    gnStringSpec();
    gnStringSpec(const std::string& m_string, gnSeqI startI = 0, gnSeqI endI = GNSEQI_END, boolean revComp = false);
    gnStringSpec(const gnStringSpec& s);
    ~gnStringSpec() override;

    gnStringSpec* Clone() const override;
    void Clear() override;

    [[nodiscard]] gnSeqI GetSourceLength() const noexcept override;
    [[nodiscard]] gnBaseSource* GetSource() const override;

    /**
     * Copies a specified range of bases and returns a pointer to
     * the resulting gnStringSpec. You must delete the copy when done.
     */
    gnStringSpec* CloneRange(gnSeqI startI, gnSeqI len) const override;

protected:
    boolean Read(gnSeqI start, gnSeqC* buf, gnSeqI& bufLen) const override;

    std::string m_seqString;
}; // class gnStringSpec

// ---- Inline functions ----

inline gnStringSpec* gnStringSpec::Clone() const {
    return new gnStringSpec(*this);
}
inline gnSeqI gnStringSpec::GetSourceLength() const noexcept {
    return static_cast<gnSeqI>(m_seqString.length());
}
inline gnBaseSource* gnStringSpec::GetSource() const {
    return nullptr;
}
inline boolean gnStringSpec::Read(gnSeqI start, gnSeqC* buf, gnSeqI& bufLen) const {
    std::memcpy(buf, m_seqString.data() + start, bufLen);
    return true;
}

} // end namespace genome

#endif // _gnStringSpec_h_
