// gnABISource.cpp
#include "libGenome/gnABISource.h"
#include <memory>

namespace genome {

gnABISource::gnABISource()
    : gnFileSource(), m_spec(nullptr)
{
}

gnABISource::gnABISource(const gnABISource& s)
    : gnFileSource(s), m_spec(nullptr)
{
    if (s.m_spec)
        m_spec = s.m_spec->Clone();
    m_contigList.reserve(s.m_contigList.size());
    for (auto* contig : s.m_contigList)
        m_contigList.push_back(contig ? contig->Clone() : nullptr);
}

gnABISource::~gnABISource()
{
    delete m_spec;
    for (auto* contig : m_contigList)
        delete contig;
}

uint32 gnABISource::GetContigListLength() const
{
    return static_cast<uint32>(m_contigList.size());
}

bool gnABISource::HasContig(const std::string& name) const
{
    for (const auto* contig : m_contigList) {
        if (contig && contig->GetName() == name)
            return true;
    }
    return false;
}

uint32 gnABISource::GetContigID(const std::string& name) const
{
    for (uint32 i = 0; i < m_contigList.size(); ++i) {
        if (m_contigList[i] && m_contigList[i]->GetName() == name)
            return i;
    }
    return static_cast<uint32>(-1);
}

std::string gnABISource::GetContigName(const uint32 i) const
{
    if (i < m_contigList.size() && m_contigList[i])
        return m_contigList[i]->GetName();
    return {};
}

gnSeqI gnABISource::GetContigSeqLength(const uint32 i) const
{
    if (i < m_contigList.size() && m_contigList[i])
        return m_contigList[i]->GetSeqLength();
    return 0;
}

gnFileContig* gnABISource::GetContig(const uint32 i) const
{
    if (i < m_contigList.size())
        return m_contigList[i];
    return nullptr;
}

bool gnABISource::SeqRead(const gnSeqI, char*, gnSeqI&, const uint32)
{
    // Not implemented
    return false;
}

bool gnABISource::Write(gnSequence&, const std::string&)
{
    // Not implemented
    return false;
}

gnGenomeSpec* gnABISource::GetSpec() const
{
    return m_spec ? m_spec->Clone() : nullptr;
}

gnFileContig* gnABISource::GetFileContig(const uint32 contigI) const
{
    return GetContig(contigI);
}

// Private methods - not implemented
bool gnABISource::SeqSeek(const gnSeqI, const uint32&, uint64&, uint64&)
{
    return false;
}

bool gnABISource::SeqStartPos(const gnSeqI, gnFileContig&, uint64&, uint64&)
{
    return false;
}

bool gnABISource::ParseStream(std::istream&)
{
    return false;
}

} // namespace genome
