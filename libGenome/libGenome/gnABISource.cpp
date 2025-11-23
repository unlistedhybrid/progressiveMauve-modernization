#include "libGenome/gnABISource.h"
#include <algorithm>

namespace genome {

gnABISource::gnABISource()
    : gnFileSource(), m_spec(nullptr), m_openString("")
{
    m_pFilter = gnFilter::fullDNASeqFilter();
}

gnABISource::gnABISource(const gnABISource& s)
    : gnFileSource(s), m_spec(nullptr), m_openString(s.m_openString), m_pFilter(s.m_pFilter)
{
    if (s.m_spec)
        m_spec = s.m_spec->Clone();
    m_contigList.reserve(s.m_contigList.size());
    for (auto* contig : s.m_contigList) {
        m_contigList.push_back(contig ? contig->Clone() : nullptr);
    }
}

gnABISource::~gnABISource()
{
    for (auto* contig : m_contigList) {
        delete contig;
    }
    m_contigList.clear();
    delete m_spec;
}

bool gnABISource::HasContig(const std::string& name) const
{
    return std::any_of(m_contigList.begin(), m_contigList.end(),
        [&name](const gnFileContig* contig) {
            return contig && contig->GetName() == name;
        });
}

uint32 gnABISource::GetContigID(const std::string& name) const
{
    for (uint32 i = 0; i < m_contigList.size(); ++i) {
        if (m_contigList[i] && m_contigList[i]->GetName() == name)
            return i;
    }
    return ALL_CONTIGS;
}

std::string gnABISource::GetContigName(uint32 i) const
{
    if (i < m_contigList.size() && m_contigList[i])
        return m_contigList[i]->GetName();
    return "";
}

gnSeqI gnABISource::GetContigSeqLength(uint32 i) const
{
    if (i < m_contigList.size() && m_contigList[i])
        return m_contigList[i]->GetSeqLength();
    else if (i == ALL_CONTIGS) {
        gnSeqI seqlen = 0;
        for (const auto* contig : m_contigList)
            if (contig) seqlen += contig->GetSeqLength();
        return seqlen;
    }
    return GNSEQI_ERROR;
}

bool gnABISource::SeqRead(const gnSeqI /*start*/, char* /*buf*/, gnSeqI& /*bufLen*/, const uint32 /*contigI*/)
{
    return false;
}

bool gnABISource::SeqSeek(const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes)
{
    if (contigI < m_contigList.size() && m_contigList[contigI])
        return SeqStartPos(start, *m_contigList[contigI], startPos, readableBytes);
    return false;
}

bool gnABISource::SeqStartPos(const gnSeqI /*start*/, gnFileContig& /*contig*/, uint64& /*startPos*/, uint64& /*readableBytes*/)
{
    return false;
}

gnFileContig* gnABISource::GetFileContig(const uint32 contigI) const
{
    if (contigI < m_contigList.size())
        return m_contigList[contigI];
    return nullptr;
}

gnGenomeSpec* gnABISource::GetSpec() const
{
    return m_spec ? m_spec->Clone() : nullptr;
}

bool gnABISource::ParseStream(std::istream& /*fin*/)
{
    return false;
}

bool gnABISource::Write(gnSequence& /*seq*/, const std::string& /*filename*/)
{
    return false;
}

} // namespace genome
