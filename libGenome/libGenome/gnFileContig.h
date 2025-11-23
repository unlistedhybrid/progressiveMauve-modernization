/////////////////////////////////////////////////////////////////////////////
// File:            gnFileContig.h
// Purpose:         File Position holder.
// Description:     Holds start/end offsets for sequence data on disk,
//                  repeat gap info, and validates section layout.
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnFileContig_h_
#define _gnFileContig_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <utility>
#include "libGenome/gnClone.h"

namespace genome {

/**
 * Tracks on-disk location and repeat-gaps of contig sequence data.
 */
class GNDLLEXPORT gnFileContig : public gnClone
{
public:
    gnFileContig() noexcept;
    gnFileContig(std::string nameStr, uint64 pos, uint64 len) noexcept;
    gnFileContig(const gnFileContig& fc) = default;
    gnFileContig& operator=(const gnFileContig& fc) = default;
    ~gnFileContig() override = default;

    gnFileContig* Clone() const override { return new gnFileContig(*this); }
    void Clear() noexcept;

    [[nodiscard]] std::string GetName() const noexcept;
    [[nodiscard]] gnSeqI GetSeqLength() const noexcept;
    [[nodiscard]] std::pair<uint64, uint64> GetFileStartEnd() const noexcept;
    [[nodiscard]] uint64 GetFileLength() const noexcept;
    [[nodiscard]] std::pair<uint64, uint64> GetSectStartEnd(gnContigSection i) const noexcept;
    [[nodiscard]] uint64 GetSectLength(gnContigSection i) const noexcept;
    [[nodiscard]] bool HasRepeatSeqGap() const noexcept;
    [[nodiscard]] std::pair<uint64, uint64> GetRepeatSeqGapSize() const noexcept;

    bool SetName(const std::string& nameStr) noexcept;
    bool SetSeqLength(gnSeqI len) noexcept;
    bool AddToSeqLength(gnSeqI len) noexcept;
    bool SetFileStart(uint64 s) noexcept;
    bool SetFileEnd(uint64 e) noexcept;
    bool SetFileStartEnd(const std::pair<uint64, uint64>& se) noexcept;
    bool SetSectStart(gnContigSection i, uint64 s) noexcept;
    bool SetSectEnd(gnContigSection i, uint64 e) noexcept;
    bool SetSectStartEnd(gnContigSection i, const std::pair<uint64, uint64>& se) noexcept;
    bool SetRepeatSeqGap(bool rsg) noexcept;
    bool SetRepeatSeqGapSize(const std::pair<uint64, uint64>& rsgSize) noexcept;
    bool SetRepeatSeqSize(uint64 seqSize) noexcept;
    bool SetRepeatGapSize(uint64 gapSize) noexcept;

private:
    std::string m_name;
    gnSeqI m_seqLength{0};
    std::pair<uint64, uint64> m_fileStartEnd{0, 0};
    std::pair<uint64, uint64> m_startEndArray[CONTIG_SECTION_SIZE]{};
    bool m_repeatSeqGap{false};
    std::pair<uint64, uint64> m_repeatSeqGapSize{0, 0};
};

// ---- Inline implementation ----

inline std::string gnFileContig::GetName() const noexcept {
    return m_name;
}
inline gnSeqI gnFileContig::GetSeqLength() const noexcept {
    return m_seqLength;
}
inline std::pair<uint64, uint64> gnFileContig::GetFileStartEnd() const noexcept {
    return m_fileStartEnd;
}
inline uint64 gnFileContig::GetFileLength() const noexcept {
    return m_fileStartEnd.second - m_fileStartEnd.first + 1;
}
inline std::pair<uint64, uint64> gnFileContig::GetSectStartEnd(gnContigSection i) const noexcept {
    if (static_cast<uint32>(i) < CONTIG_SECTION_SIZE)
        return m_startEndArray[static_cast<uint32>(i)];
    return {0, 0};
}
inline uint64 gnFileContig::GetSectLength(gnContigSection i) const noexcept {
    if (static_cast<uint32>(i) < CONTIG_SECTION_SIZE)
        return m_startEndArray[static_cast<uint32>(i)].second - m_startEndArray[static_cast<uint32>(i)].first + 1;
    return 0;
}
inline bool gnFileContig::HasRepeatSeqGap() const noexcept {
    return m_repeatSeqGap;
}
inline std::pair<uint64, uint64> gnFileContig::GetRepeatSeqGapSize() const noexcept {
    return m_repeatSeqGapSize;
}

// Setters
inline bool gnFileContig::SetName(const std::string& nameStr) noexcept {
    m_name = nameStr;
    return true;
}
inline bool gnFileContig::SetSeqLength(gnSeqI len) noexcept {
    m_seqLength = len;
    return true;
}
inline bool gnFileContig::AddToSeqLength(gnSeqI len) noexcept {
    m_seqLength += len;
    return true;
}
inline bool gnFileContig::SetFileStart(uint64 s) noexcept {
    m_fileStartEnd.first = s;
    return true;
}
inline bool gnFileContig::SetFileEnd(uint64 e) noexcept {
    m_fileStartEnd.second = e;
    return true;
}
inline bool gnFileContig::SetFileStartEnd(const std::pair<uint64, uint64>& se) noexcept {
    m_fileStartEnd = se;
    return true;
}
inline bool gnFileContig::SetSectStart(gnContigSection i, uint64 s) noexcept {
    if (static_cast<uint32>(i) < CONTIG_SECTION_SIZE) {
        m_startEndArray[static_cast<uint32>(i)].first = s;
        return true;
    }
    return false;
}
inline bool gnFileContig::SetSectEnd(gnContigSection i, uint64 e) noexcept {
    if (static_cast<uint32>(i) < CONTIG_SECTION_SIZE) {
        m_startEndArray[static_cast<uint32>(i)].second = e;
        return true;
    }
    return false;
}
inline bool gnFileContig::SetSectStartEnd(gnContigSection i, const std::pair<uint64, uint64>& se) noexcept {
    if (static_cast<uint32>(i) < CONTIG_SECTION_SIZE) {
        m_startEndArray[static_cast<uint32>(i)] = se;
        return true;
    }
    return false;
}
inline bool gnFileContig::SetRepeatSeqGap(bool rsg) noexcept {
    m_repeatSeqGap = rsg;
    return true;
}
inline bool gnFileContig::SetRepeatSeqGapSize(const std::pair<uint64, uint64>& rsgSize) noexcept {
    return SetRepeatSeqSize(rsgSize.first) && SetRepeatGapSize(rsgSize.second);
}

} // end namespace genome

#endif // _gnFileContig_h_
