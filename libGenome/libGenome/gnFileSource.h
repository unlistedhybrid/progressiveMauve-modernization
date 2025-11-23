/////////////////////////////////////////////////////////////////////////////
// File:            gnFileSource.h
// Purpose:         Implements gnBaseSource for .File files
// Description:     Standard interface for file-based genetic sources
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnFileSource_h_
#define _gnFileSource_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <fstream>
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnException.h"

namespace genome {

/**
 * Interface for all file-based sources of genetic information.
 */
class GNDLLEXPORT gnFileSource : public gnBaseSource {
public:
    gnFileSource();
    gnFileSource(const gnFileSource& gnfs);
    virtual ~gnFileSource() = default;

    [[nodiscard]] virtual gnFileSource* Clone() const = 0;

    // Open, Close
    virtual void Open(const std::string& openString);
    virtual void Open();
    virtual void Close();

    [[nodiscard]] virtual std::string GetOpenString() const noexcept;

    // Filter
    [[nodiscard]] virtual const gnFilter* GetFilter() const noexcept;
    virtual void SetFilter(gnFilter* filter);

    virtual boolean Read(uint64 pos, char* buf, gnSeqI& bufLen);

    /**
     * Returns a pointer to the file contig for contigI, or nullptr if not found.
     */
    virtual gnFileContig* GetFileContig(uint32 contigI) const = 0;

protected:
    void DetermineNewlineType();

    std::string m_openString;
    std::ifstream m_ifstream;
    const gnFilter* m_pFilter{nullptr};
    gnNewlineType m_newlineType{GN_NEWLINE_UNKNOWN};
    uint32 m_newlineSize{0};

private:
    virtual boolean ParseStream(std::istream& fin) = 0;
}; // class gnFileSource

// ---- Inline implementations ----

inline std::string gnFileSource::GetOpenString() const noexcept {
    return m_openString;
}
inline const gnFilter* gnFileSource::GetFilter() const noexcept {
    return m_pFilter;
}
inline void gnFileSource::SetFilter(gnFilter* filter) {
    if (filter == nullptr) {
        Throw_gnEx(NullPointer());
    }
    m_pFilter = filter;
}

} // end namespace genome

#endif // _gnFileSource_h_
