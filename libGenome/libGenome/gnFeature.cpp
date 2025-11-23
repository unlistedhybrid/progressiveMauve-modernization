// gnFeature.cpp
/////////////////////////////////////////////////////////////////////////////
// File:            gnFeature.cpp
// Purpose:         implements the gnBaseFeature for generic features
// Author:          Aaron Darling 
// License:         See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFeature.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringQualifier.h"
#include "libGenome/gnDebug.h"
#include <list>

namespace genome {

gnFeature::gnFeature()
    : gnBaseFeature()  // If gnBaseFeature has a default constructor
{
    // m_id = 0;
    // m_name = "";
    // m_locationType = gnLocation::LT_Nothing;
    // m_broken = false;
}

gnFeature::gnFeature(const std::string& name, uint32 id,
                     gnLocation::gnLocationType lt, bool broken)
    : gnBaseFeature(name, id, nullptr, lt, broken)
{
}

gnFeature::gnFeature(const gnFeature& s)
    : gnBaseFeature(s)
{
    m_id = s.m_id;
    m_name = s.m_name;
    m_locationType = s.m_locationType;
    m_broken = s.m_broken;
    m_spec = s.m_spec;
    for (const auto& loc : s.m_locationList)
        m_locationList.push_back(loc);
    for (const auto& qual : s.m_qualifierList)
        m_qualifierList.push_back(qual->Clone());
}

gnFeature::~gnFeature()
{
    // If m_qualifierList stores owning pointers, consider deleting them here.
    // If using smart pointers, this is not needed.
}

} // end namespace genome
