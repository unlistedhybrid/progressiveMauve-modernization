/////////////////////////////////////////////////////////////////////////////
// File:            gnFragmentSpec.h
// Purpose:         abstract Spec class (Genome-level)
// Description:     Contains a list of specs/features/headers for a genome fragment.
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnFragmentSpec_h_
#define _gnFragmentSpec_h_

#include "libGenome/gnDefs.h"
#include <vector>
#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnMultiSpec.h"
#include "libGenome/gnException.h"

namespace genome {

/**
 * gnFragmentSpec: manages a list of gnContigSpec fragments and associated features/headers.
 * Used by file reader classes (like gnGBKSource).
 */
class GNDLLEXPORT gnFragmentSpec : public gnMultiSpec<gnContigSpec>
{
public:
    gnFragmentSpec();
    ~gnFragmentSpec() override;
    gnFragmentSpec(const gnFragmentSpec& s);

    gnFragmentSpec* Clone() const override;
    void Clear() override;
    void SetReverseComplement(const boolean value) override;

    void CropStart(gnSeqI cropLen) override;
    void CropEnd(gnSeqI cropLen) override;

    uint32 AddFeature(gnBaseFeature* feat);
    [[nodiscard]] uint32 GetFeatureListLength() const;
    [[nodiscard]] gnBaseFeature* GetFeature(uint32 i) const;
    void GetContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
    void GetIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
    void GetBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const;
    void RemoveFeature(uint32 i);

    /**
     * Copies a specified range of bases and returns a pointer to
     * the resulting gnFragmentSpec.  You must delete the copy when done.
     */
    gnFragmentSpec* CloneRange(gnSeqI startI, gnSeqI len) const override;

protected:
    std::vector<gnBaseFeature*> m_featureList;
}; // class gnFragmentSpec

// ---- Inline ----

inline gnFragmentSpec* gnFragmentSpec::Clone() const {
    return new gnFragmentSpec(*this);
}
inline uint32 gnFragmentSpec::AddFeature(gnBaseFeature* feat) {
    m_featureList.push_back(feat);
    feat->SetSpec(this);
    return static_cast<uint32>(m_featureList.size() - 1);
}
inline uint32 gnFragmentSpec::GetFeatureListLength() const {
    return static_cast<uint32>(m_featureList.size());
}
inline gnBaseFeature* gnFragmentSpec::GetFeature(uint32 i) const {
    if (i >= m_featureList.size())
        Throw_gnEx(FeatureIndexOutOfBounds());
    return m_featureList[i]->Clone();
}
inline void gnFragmentSpec::RemoveFeature(uint32 i) {
    if (i >= m_featureList.size())
        Throw_gnEx(FeatureIndexOutOfBounds());
    delete m_featureList[i];
    m_featureList.erase(m_featureList.begin() + i);
}

} // end namespace genome

#endif // _gnFragmentSpec_h_
