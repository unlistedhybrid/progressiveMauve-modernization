/////////////////////////////////////////////////////////////////////////////
// File:            gnFASSource.cpp
// Purpose:         Implements gnBaseSource for .FAS files
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFASSource.h"
#include "libGenome/gnBaseSpec.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnDebug.h"
#include <fstream>

namespace genome {

gnFASSource::gnFASSource() {
    m_openString.clear();
    m_pFilter = gnFilter::fullDNASeqFilter();
    if (m_pFilter == nullptr) {
        DebugMsg("Error using static sequence filters.");
    }
}

gnFASSource::gnFASSource(const gnFASSource& s) : gnFileSource(s) {
    for (const auto* contig : s.m_contigList) {
        m_contigList.push_back(contig->Clone());
    }
}

gnFASSource::~gnFASSource() {
    m_ifstream.close();
    for (auto* fg : m_contigList) {
        delete fg;
    }
    m_contigList.clear();
}

bool gnFASSource::HasContig(const std::string& nameStr) const {
    for (const auto* contig : m_contigList) {
        if (nameStr == contig->GetName())
            return true;
    }
    return false;
}

uint32 gnFASSource::GetContigID(const std::string& name) const {
    for (uint32 i = 0; i < m_contigList.size(); ++i) {
        if (name == m_contigList[i]->GetName())
            return i;
    }
    return ALL_CONTIGS;
}

std::string gnFASSource::GetContigName(const uint32 i) const {
    if (i < m_contigList.size()) {
        return m_contigList[i]->GetName();
    }
    return "";
}

gnSeqI gnFASSource::GetContigSeqLength(const uint32 i) const {
    if (i < m_contigList.size()) {
        return m_contigList[i]->GetSeqLength();
    } else if (i == ALL_CONTIGS) {
        gnSeqI seqlen = 0;
        for (const auto* contig : m_contigList)
            seqlen += contig->GetSeqLength();
        return seqlen;
    }
    return GNSEQI_ERROR;
}

gnFileContig* gnFASSource::GetContig(const uint32 i) const {
    if (i < m_contigList.size()) {
        return m_contigList[i];
    }
    return nullptr;
}

bool gnFASSource::SeqRead(const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI) {
    return SeqReadImpl(start, buf, bufLen, contigI);
}

bool gnFASSource::SeqReadImpl(const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI) {
    m_ifstream.clear();
    uint32 contigIndex = contigI;
    uint64 startPos = 0;
    uint64 readableBytes = 0;
    if (!SeqSeek(start, contigIndex, startPos, readableBytes)) {
        bufLen = 0;
        return false;
    }

    if (contigI == ALL_CONTIGS) {
        uint32 curLen = 0;
        uint64 bytesRead = 0;
        while (curLen < bufLen) {
            if (readableBytes <= 0) {
                if (!SeqSeek(start + curLen, contigIndex, startPos, readableBytes)) {
                    bufLen = curLen;
                    return true;
                }
            }
            uint64 readLen = std::min<uint64>(bufLen - curLen, readableBytes);
            Array<gnSeqC> array_buf(readLen);
            gnSeqC* tmpBuf = array_buf.data;

            m_ifstream.read(tmpBuf, readLen);
            uint64 gc = m_ifstream.gcount();
            bytesRead += gc;
            readableBytes -= gc;

            for (uint32 i = 0; i < gc; i++) {
                if (m_pFilter->IsValid(tmpBuf[i])) {
                    buf[curLen] = tmpBuf[i];
                    curLen++;
                }
            }
            if (m_ifstream.eof()) {
                m_ifstream.clear();
                bufLen = curLen;
                return true;
            }
        }
        bufLen = curLen;
    } else if (contigI < m_contigList.size()) {
        uint32 curLen = 0;
        gnSeqI contigSize = m_contigList[contigI]->GetSeqLength();
        bufLen = std::min(bufLen, contigSize);
        while (curLen < bufLen) {
            uint64 readLen = bufLen - curLen;
            Array<gnSeqC> array_buf(readLen);
            gnSeqC* tmpBuf = array_buf.data;

            m_ifstream.read(tmpBuf, readLen);
            uint64 gc = m_ifstream.gcount();

            for (uint32 i = 0; i < gc; i++) {
                if (m_pFilter->IsValid(tmpBuf[i])) {
                    buf[curLen] = tmpBuf[i];
                    curLen++;
                }
            }
            if (m_ifstream.eof()) {
                m_ifstream.clear();
                bufLen = curLen;
                return true;
            }
        }
        bufLen = curLen;
    }
    return true;
}

bool gnFASSource::SeqSeek(const gnSeqI start, const uint32 contigI, uint64& startPos, uint64& readableBytes) {
    if (contigI == ALL_CONTIGS) {
        gnSeqI curIndex = 0;
        auto iter = m_contigList.begin();
        for (; iter != m_contigList.end(); ++iter) {
            uint64 len = (*iter)->GetSeqLength();
            if ((curIndex + len) > start)
                break;
            curIndex += len;
        }
        if (iter == m_contigList.end())
            return false;
        gnSeqI startIndex = start - curIndex;
        return SeqStartPos(startIndex, *(*iter), startPos, readableBytes);
    } else if (contigI < m_contigList.size()) {
        return SeqStartPos(start, *(m_contigList[contigI]), startPos, readableBytes);
    }
    return false;
}

bool gnFASSource::SeqStartPos(const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes) {
    readableBytes = 0;
    uint32 curLen = 0;
    startPos = contig.GetSectStartEnd(gnContigSequence).first;

    if (contig.HasRepeatSeqGap())
        if (contig.GetRepeatSeqGapSize().first > 0)
            if (contig.GetRepeatSeqGapSize().second > 0) {
                startPos += start + (start / contig.GetRepeatSeqGapSize().first) * contig.GetRepeatSeqGapSize().second;
                readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
                m_ifstream.seekg(startPos, std::ios::beg);
                return true;
            }

    m_ifstream.seekg(startPos, std::ios::beg);
    if (m_ifstream.eof()) {
        ErrorMsg("ERROR in gnFASSource::Incorrect contig start position, End of file reached!\n");
        return false;
    }
    while (true) {
        uint32 tmpbufsize = contig.GetSectStartEnd(gnContigSequence).second - startPos;
        if (tmpbufsize == 0) {
            ErrorMsg("ERROR in gnFASSource: stored contig size is incorrect.\n");
            return false;
        }
        tmpbufsize = std::min<uint32>(tmpbufsize, BUFFER_SIZE);
        Array<char> array_buf(tmpbufsize);
        char* tmpbuf = array_buf.data;
        m_ifstream.read(tmpbuf, tmpbufsize);
        if (m_ifstream.eof()) {
            ErrorMsg("ERROR in gnFASSource::Read End of file reached!\n");
            return false;
        }
        for (uint32 i = 0; i < tmpbufsize; ++i) {
            if (m_pFilter->IsValid(tmpbuf[i])) {
                if (curLen >= start) {
                    startPos += i;
                    m_ifstream.seekg(startPos, std::ios::beg);
                    readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
                    return true;
                }
                ++curLen;
            }
        }
        startPos += tmpbufsize;
    }
    return true;
}

// Deprecated write
bool gnFASSource::Write(gnBaseSource* source, const std::string& filename) {
    std::ofstream m_ofstream(filename, std::ios::out | std::ios::binary);
    if (!m_ofstream.is_open())
        return false;
    uint32 contigCount = source->GetContigListLength();
    for (uint32 contigI = 0; contigI < contigCount; contigI++) {
        std::string contigName = source->GetContigName(contigI);
        m_ofstream << ">" << contigName << ";\n";
        gnSeqI seqLength = source->GetContigSeqLength(contigI);
        while (seqLength > 0) {
            gnSeqI writeLen = std::min(seqLength, static_cast<gnSeqI>(BUFFER_SIZE));
            Array<gnSeqC> array_buf(writeLen + 1);
            gnSeqC* bases = array_buf.data;
            bool success = source->SeqRead(0, bases, writeLen, contigI);
            if (!success)
                return false;
            bases[writeLen] = 0;
            m_ofstream << bases << "\n";
            seqLength -= writeLen;
        }
    }
    m_ofstream.close();
    return true;
}

void gnFASSource::Write(gnSequence& seq, const std::string& filename, bool write_coords, bool enforce_unique_names) {
    std::ofstream m_ofstream(filename, std::ios::out | std::ios::binary);
    if (!m_ofstream.is_open())
        Throw_gnEx(FileNotOpened());
    Write(seq, m_ofstream, write_coords, enforce_unique_names);
    m_ofstream.close();
}

void gnFASSource::Write(gnSequence& seq, std::ostream& m_ostream, bool write_coords, bool enforce_unique_names) {
    std::vector<std::string> contigNameList;
    Array<gnSeqC> array_buf(BUFFER_SIZE);
    gnSeqC* bases = array_buf.data;
    gnGenomeSpec* spec = seq.GetSpec();

    std::string newline = "\r\n";
    gnSeqI readOffset = 1;

    for (uint32 fragI = 0; fragI < seq.contigListLength(); fragI++) {
        std::string contigName = seq.contigName(fragI);

        if (enforce_unique_names) {
            uint32 name_count = 0;
            for (const auto& s : contigNameList)
                if (s == contigName)
                    name_count++;
            contigNameList.push_back(contigName);
            if (name_count > 0)
                contigName += "_" + uintToString(name_count);
        }

        gnFragmentSpec* subSpec = spec->GetSpec(fragI);
        gnBaseHeader* gpbh = subSpec->GetHeader(0);
        std::string header;
        if (gpbh != nullptr) {
            header = gpbh->GetHeader();
            auto newlinepos = header.find_first_of('\n', 0);
            if (newlinepos != std::string::npos) {
                if (newlinepos > 0 && header[newlinepos - 1] == '\r')
                    newlinepos--;
                header = header.substr(0, newlinepos);
            }
        }

        gnSeqI readLength = seq.contigLength(fragI);
        m_ostream << ">" << contigName;
        if (write_coords)
            m_ostream << " (" << readOffset << ", " << readOffset + readLength - 1 << ")";
        m_ostream << "; " << header << std::endl;

        gnSeqI linePos = 0;
        while (readLength > 0) {
            gnSeqI read_chunk_size = (BUFFER_SIZE / FAS_LINE_WIDTH) * FAS_LINE_WIDTH;
            gnSeqI writeLen = std::min(readLength, read_chunk_size);

            if (!seq.ToArray(bases, writeLen, readOffset))
                return;
            for (gnSeqI curbase = 0; curbase < writeLen; curbase += FAS_LINE_WIDTH) {
                gnSeqI writeout_size = std::min(writeLen - curbase, static_cast<gnSeqI>(FAS_LINE_WIDTH));
                m_ostream << std::string(bases + curbase, writeout_size) << std::endl;
            }
            readLength -= writeLen;
            readOffset += writeLen;
        }
        if (linePos != 0)
            m_ostream << std::endl;
    }
}

gnGenomeSpec* gnFASSource::GetSpec() const {
    auto* spec = new gnGenomeSpec();
    for (uint32 i = 0; i < m_contigList.size(); i++) {
        gnFragmentSpec* fragmentSpec = new gnFragmentSpec();
        gnSourceSpec* contigSpec = new gnSourceSpec(static_cast<gnBaseSource*>(const_cast<gnFASSource*>(this)), i);
        spec->AddSpec(fragmentSpec, i);
        fragmentSpec->AddSpec(contigSpec);

        fragmentSpec->SetName(m_contigList[i]->GetName());
        fragmentSpec->SetSourceName(m_openString);
        contigSpec->SetName(m_contigList[i]->GetName());
        contigSpec->SetSourceName(m_openString);

        auto headsect = m_contigList[i]->GetSectStartEnd(gnContigHeader);
        if (headsect.first != headsect.second) {
            gnSourceHeader* gpsh = new gnSourceHeader(static_cast<gnBaseSource*>(const_cast<gnFASSource*>(this)), std::string(""), headsect.first, headsect.second - headsect.first);
            fragmentSpec->AddHeader(gpsh, 0);
        }
    }
    return spec;
}

gnFileContig* gnFASSource::GetFileContig(const uint32 contigI) const {
    if (m_contigList.size() > contigI)
        return m_contigList[contigI];
    return nullptr;
}

bool gnFASSource::ParseStream(std::istream& fin) {
    uint32 readState = 0;
    gnFileContig* currentContig = nullptr;
    std::string nameFStr;
    uint64 seqLength = 0, gapLength = 0;
    uint64 streamPos = 0;
    uint64 bufReadLen = 0;
    Array<gnSeqC> array_buf(BUFFER_SIZE);
    char* buf = array_buf.data;
    bool paren_hit = false;
    uint32 repeatSeqSize = 0;
    bool corrupt_msg = false;

    DetermineNewlineType();

    while (!fin.eof()) {
        fin.read(buf, BUFFER_SIZE);
        streamPos += bufReadLen;
        bufReadLen = fin.gcount();
        for (uint32 i = 0; i < bufReadLen; i++) {
            char ch = buf[i];
            switch (readState) {
                case 0:
                    if (!((buf[0] == '>') || !m_pFilter->IsValid(buf[0]))) {
                        return false;
                    }
                    readState = 1;
                    [[fallthrough]];
                case 1:
                    if (ch == '>') {
                        seqLength = 0; gapLength = 0;
                        currentContig = new gnFileContig();
                        currentContig->SetFileStart(streamPos + i);
                        currentContig->SetRepeatSeqGap(true);
                        currentContig->SetRepeatSeqSize(repeatSeqSize);
                        currentContig->SetRepeatGapSize(m_newlineSize);
                        readState = 2;
                        paren_hit = false;
                        nameFStr.clear();
                    } else {
                        ++gapLength;
                    }
                    break;
                case 2:
                    if (isNewLine(ch) || ch == ';') {
                        currentContig->SetName(nameFStr);
                        currentContig->SetSectStart(gnContigHeader, streamPos + i + 1);
                        if (ch == ';')
                            readState = 3;
                        else {
                            readState = 4;
                            if (ch == '\r')
                                currentContig->SetSectStart(gnContigHeader, streamPos + i + 2);
                        }
                    } else if (ch == '(') {
                        if (i > 0 && isSpace(buf[i - 1])) {
                            nameFStr = nameFStr.substr(0, nameFStr.length() - 1);
                        }
                        paren_hit = true;
                    } else if ((!isSpace(ch) || !nameFStr.empty()) && !paren_hit) {
                        nameFStr += ch;
                    }
                    break;
                case 3:
                    if (isNewLine(ch))
                        readState = 4;
                    break;
                case 4:
                    if (ch == '>') {
                        readState = 3;
                    } else if (m_pFilter->IsValid(ch)) {
                        currentContig->SetSectEnd(gnContigHeader, streamPos + i);
                        currentContig->SetSectStart(gnContigSequence, streamPos + i);
                        seqLength = 1; gapLength = 0;
                        readState = 5;
                    }
                    break;
                case 5:
                    while (i < bufReadLen) {
                        ch = buf[i];
                        if (m_pFilter->IsValid(ch)) {
                            if (gapLength > 0) {
                                if (seqLength != repeatSeqSize) {
                                    if (!corrupt_msg) {
                                        corrupt_msg = true;
                                        ErrorMsg("Sequence file appears corrupt, proceeding with caution\n");
                                    }
                                    currentContig->SetRepeatSeqGap(false);
                                }
                                if (gapLength != m_newlineSize) {
                                    if (!corrupt_msg) {
                                        corrupt_msg = true;
                                        ErrorMsg("Sequence file appears corrupt, proceeding with caution\n");
                                    }
                                    currentContig->SetRepeatSeqGap(false);
                                }
                                currentContig->AddToSeqLength(seqLength);
                                seqLength = 0;
                                gapLength = 0;
                            }
                            seqLength++;
                        } else if (ch == '>') {
                            currentContig->AddToSeqLength(seqLength);
                            currentContig->SetSectEnd(gnContigSequence, streamPos + i - 1);
                            currentContig->SetFileEnd(streamPos + i - 1);
                            m_contigList.push_back(currentContig);
                            readState = 1;
                            i--;
                            break;
                        } else if (isNewLine(ch)) {
                            if (repeatSeqSize == 0) {
                                repeatSeqSize = seqLength;
                                currentContig->SetRepeatSeqSize(repeatSeqSize);
                            }
                            gapLength++;
                        } else {
                            currentContig->SetRepeatSeqGap(false);
                            if (!corrupt_msg) {
                                corrupt_msg = true;
                                ErrorMsg("Sequence file appears corrupt, proceeding with caution\n");
                            }
                        }
                        i++;
                    }
                    break;
                default:
                    ErrorMsg("ERROR");
                    return false;
            }
        }
    }
    if (currentContig != nullptr) {
        if (readState == 2) {
            currentContig->SetName(nameFStr);
        }
        if ((readState >= 2) && (readState < 5)) {
            currentContig->SetSectEnd(gnContigHeader, streamPos + bufReadLen);
        } else if (readState == 5) {
            currentContig->AddToSeqLength(seqLength);
            currentContig->SetSectEnd(gnContigSequence, streamPos + bufReadLen);
        }
        currentContig->SetFileEnd(streamPos + bufReadLen);
        m_contigList.push_back(currentContig);
    }
    m_ifstream.clear();
    return true;
}

} // end namespace genome
