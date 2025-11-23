#ifndef GN_ABI_SOURCE_H
#define GN_ABI_SOURCE_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <iosfwd>

#include "libGenome/gnFileSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnFilter.h"

namespace genome {

using uint32 = std::uint32_t;
using gnSeqI = std::int64_t;
using uint64 = std::uint64_t;

constexpr uint32 ALL_CONTIGS = static_cast<uint32>(-1);
constexpr gnSeqI GNSEQI_ERROR = -1;

class gnABISource : public gnFileSource {
public:
    gnABISource();
    gnABISource(const gnABISource& s);
    ~gnABISource();

    // Contig access methods
    bool HasContig(const std::string& name) const;
    uint32 GetContigID(const std::string& name) const;
    std::string GetContigName(uint32 i) const;
    gnSeqI GetContigSeqLength(uint32 i) const;
    
    gnFileContig* GetFileContig(const uint32 contigI) const;
    gnGenomeSpec* GetSpec() const;

    // Sequence reading
    bool SeqRead(const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI);
    bool SeqSeek(const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes);
    bool SeqStartPos(const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes);

    // Stream parsing
    bool ParseStream(std::istream& fin);

protected:
    std::vector<gnFileContig*> m_contigList;
    gnGenomeSpec* m_spec;
    std::string m_openString;
    std::shared_ptr<gnFilter> m_pFilter;

private:
    // Not implemented
    bool Write(gnSequence& seq, const std::string& filename);
};

} // namespace genome

#endif // GN_ABI_SOURCE_H
