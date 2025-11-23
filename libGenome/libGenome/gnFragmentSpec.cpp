/////////////////////////////////////////////////////////////////////////////
// File:            gnFragmentSpec.cpp
// Purpose:         implements gnMultiSpec< gnContigSpec > for sequence fragments
// Description:     
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFragmentSpec.h"
#include <string>
#include <vector>

namespace genome {

gnFragmentSpec::gnFragmentSpec() {
    gnBaseSpec::Clear();
}

gnFragmentSpec::gnFragmentSpec(const gnFragmentSpec& s)
    : gnMultiSpec<gnContigSpec>(s)
{
    m_sourceName = s.m_sourceName;
    m_name = s.m_name;
    m_reverseComplement = s.m_reverseComplement;
    m_circular = s.m_circular;

    // Copy the header list.
    uint32 list_size = static_cast<uint32>(s.m_headerList.size());
    m_headerList.reserve(list_size);
    for (uint32 i = 0; i < list_size; i++)
        m_headerList.push_back(s.m_headerList[i]->Clone());

    // Copy the contig list.
    list_size = static_cast<uint32>(s.m_SpecList.size());
    m_SpecList.reserve(list_size);
    for (uint32 i = 0; i < list_size; i++)
        m_SpecList.push_back(s.m_SpecList[i]->Clone());

    // Copy the feature list.
    list_size = static_cast<uint32>(s.m_featureList.size());
    m_featureList.reserve(list_size);
    for (uint32 i = 0; i < list_size; i++)
        m_featureList.push_back(s.m_featureList[i]->Clone());
}

gnFragmentSpec::~gnFragmentSpec() {
    Clear();
}

void gnFragmentSpec::Clear() {
    for (auto* s : m_SpecList)
        delete s;
    m_SpecList.clear();
    for (auto* f : m_featureList)
        delete f;
    m_featureList.clear();
    gnMultiSpec<gnContigSpec>::Clear();
}

void gnFragmentSpec::GetContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const {
    for (uint32 i = 0; i < m_featureList.size(); i++) {
        if (m_featureList[i]->IsContainedBy(lt)) {
            feature_vector.push_back(m_featureList[i]->Clone());
            index_vector.push_back(i);
        }
    }
}
void gnFragmentSpec::GetIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const {
    for (uint32 i = 0; i < m_featureList.size(); i++) {
        if (m_featureList[i]->Intersects(lt)) {
            feature_vector.push_back(m_featureList[i]->Clone());
            index_vector.push_back(i);
        }
    }
}
void gnFragmentSpec::GetBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const {
    for (uint32 i = 0; i < m_featureList.size(); i++)
        if (m_featureList[i]->IsBroken() && m_featureList[i]->IsContainedBy(lt))
            feature_vector.push_back(m_featureList[i]->Clone());
}

void gnFragmentSpec::CropStart(gnSeqI cropLen) {
    for (auto* feature : m_featureList)
        feature->CropStart(cropLen);
    gnMultiSpec<gnContigSpec>::CropStart(cropLen);
}

void gnFragmentSpec::CropEnd(gnSeqI cropLen) {
    for (auto* feature : m_featureList)
        feature->CropEnd(cropLen);
    gnMultiSpec<gnContigSpec>::CropEnd(cropLen);
}

gnFragmentSpec* gnFragmentSpec::CloneRange(const gnSeqI startI, const gnSeqI len) const {
    if (len == 0)
        return new gnFragmentSpec();

    // Find the valid range of specs to copy
    uint32 firstSpec = GetSpecIndexByBase(startI);
    gnSeqI total_copylen = len;
    uint32 endSpec;
    if (len != GNSEQI_END) {
        endSpec = GetSpecIndexByBase(startI + len - 1);
    } else {
        endSpec = GetSpecListLength() - 1;
        total_copylen = GetLength() - startI;
    }

    // Find their starting and ending bases
    gnSeqI firstBase = startI - GetSpecStartBase(firstSpec);
    gnSeqI firstSpecLen = GetSpec(firstSpec)->GetLength();
    boolean spans_specs = true;
    gnSeqI firstCopyLen = firstSpecLen - firstBase;
    if (firstCopyLen >= total_copylen) {
        spans_specs = false;
        firstCopyLen = total_copylen;
    }
    gnFragmentSpec* destSpec = new gnFragmentSpec();
    gnContigSpec* newSpec = m_SpecList[firstSpec]->CloneRange(firstBase, firstCopyLen);
    destSpec->AddSpec(newSpec);

    gnSeqI cur_copylen = firstCopyLen;
    // Add all the completely covered specs in the middle
    for (uint32 specI = firstSpec + 2; specI <= endSpec; specI++) {
        destSpec->AddSpec(GetSpec(specI - 1)->Clone());
        cur_copylen += GetSpec(specI - 1)->GetLength();
    }
    // Add the last spec if necessary
    if (spans_specs) {
        newSpec = m_SpecList[endSpec]->CloneRange(0, total_copylen - cur_copylen);
        destSpec->AddSpec(newSpec);
    }

    // Now clone all the appropriate features
    gnLocation lt;
    std::vector<gnBaseFeature*> feature_vector;
    std::vector<uint32> index_vector;
    lt.SetStart(startI);
    lt.SetEnd(startI + total_copylen);
    GetIntersectingFeatures(lt, destSpec->m_featureList, index_vector);

    return destSpec;
}

void gnFragmentSpec::SetReverseComplement(const boolean value) {
    if (value == m_reverseComplement)
        return;
    // Reverse the spec list entries
    std::vector<gnContigSpec*> tmp_spec_list;
    for (uint32 i = 0; i < GetSpecListLength(); i++) {
        // Transmit rev_comp down the tree
        GetSpec(i)->SetReverseComplement(!GetSpec(i)->IsReverseComplement());
        tmp_spec_list.insert(tmp_spec_list.begin(), GetSpec(i));
    }
    m_SpecList = std::move(tmp_spec_list);
    m_reverseComplement = value;
}

} // end namespace genome
