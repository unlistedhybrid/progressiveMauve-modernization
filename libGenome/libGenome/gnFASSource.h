/////////////////////////////////////////////////////////////////////////////
// File:            gnFASSource.h
// Purpose:         Implements gnBaseSource for .FAS files (FastA)
// Description:     FastA file read/write interface
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnFASSource_h_
#define _gnFASSource_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <fstream>
#include <vector>
#include "libGenome/gnFileSource.h"
#include "libGenome/gnSequence.h"

namespace genome {

#define FAS_LINE_WIDTH 80

/**
 * gnFASSource reads/writes FastA files for gnSourceFactory.
 */
class GNDLLEXPORT gnFASSource : public gnFileSource
{
public:
    gnFASSource();
    gnFASSource(const gnFASSource& s);
    ~gnFASSource() override;
    gnFASSource* Clone() const override;

    [[nodiscard]] uint32 GetContigListLength() const noexcept;
    [[nodiscard]] boolean HasContig(const std::string& name) const;
    [[nodiscard]] uint32 GetContigID(const std::string& name) const;
    [[nodiscard]] std::string GetContigName(uint32 i) const;
    [[nodiscard]] gnSeqI GetContigSeqLength(uint32 i) const;
    [[nodiscard]] gnFileContig* GetContig(uint32 i) const;

    boolean SeqRead(gnSeqI start, char* buf, gnSeqI& bufLen, uint32 contigI = ALL_CONTIGS);

    // Write sequence to FastA file
    static void Write(gnSequence& sequence, const std::string& filename,
                      boolean write_coords = true, boolean enforce_unique_names = true);

    // Write sequence to an output stream
    static void Write(gnSequence& sequence, std::ostream& m_ostream,
                      boolean write_coords = true, boolean enforce_unique_names = true);

    // Deprecated: Write a gnBaseSource to FastA file
    [[deprecated("Use Write(gnSequence&, ...) instead.")]]
    static boolean Write(gnBaseSource* source, const std::string& filename);

    [[nodiscard]] gnGenomeSpec* GetSpec() const;
    gnFileContig* GetFileContig(uint32 contigI) const override;

private:
    boolean SeqReadImpl(gnSeqI start, char* buf, gnSeqI& bufLen, uint32 contigI = ALL_CONTIGS);
    boolean SeqSeek(gnSeqI start, uint32 contigI, uint64& startPos, uint64& readableBytes);
    boolean SeqStartPos(gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes);
    boolean ParseStream(std::istream& fin) override;

    std::vector<gnFileContig*> m_contigList;
}; // class gnFASSource

// ---- Inline methods ----
inline gnFASSource* gnFASSource::Clone() const {
    return new gnFASSource(*this);
}
inline uint32 gnFASSource::GetContigListLength() const noexcept {
    return static_cast<uint32>(m_contigList.size());
}

} // end namespace genome

#endif // _gnFASSource_h_
