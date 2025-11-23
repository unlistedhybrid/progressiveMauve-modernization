/////////////////////////////////////////////////////////////////////////////
// File:            gnSEQSource.h
// Purpose:         Implements gnBaseSource for .SEQ files
// Description:     
// Changes:        
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFilter.h"
#include "libGenome/gnFeature.h"
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnSourceQualifier.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnDebug.h"
#include <cstring>


using namespace std;
namespace genome {


gnSEQSource::gnSEQSource()
{
	m_openString = "";
	m_pFilter = gnFilter::fullDNASeqFilter();
	if(m_pFilter == nullptr){
		DebugMsg("Error using static sequence filters.\n");
	}
}
gnSEQSource::gnSEQSource( const gnSEQSource& s ) : gnFileSource(s)
{
	vector< gnFileContig* >::const_iterator iter = s.m_contigList.begin();
	for( ; iter != s.m_contigList.end(); ++iter )
	{
		m_contigList.push_back( (*iter)->Clone() );
	}
}
gnSEQSource::~gnSEQSource()
{
	m_ifstream.close();
	vector< gnFileContig* >::iterator iter = m_contigList.begin();
	for( ; iter != m_contigList.end(); ++iter )
	{
		gnFileContig* fg = *iter;
		*iter = nullptr;
		delete fg;
	}
}
bool gnSEQSource::HasContig( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return true;
	}
	return false;
}
uint32 gnSEQSource::GetContigID( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return i;
	}
	return ALL_CONTIGS;
}
string gnSEQSource::GetContigName( const uint32 i ) const
{
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetName();
	}
	return "";
}
gnSeqI gnSEQSource::GetContigSeqLength( const uint32 i ) const
{
	if( i == ALL_CONTIGS)
		return m_spec->GetLength();
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetSeqLength();
	}
	return GNSEQI_ERROR;
}

bool gnSEQSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ){
	uint64 startPos = 0;
	uint64 readableBytes = 0;
	if( !SeqSeek( start, contigI, startPos, readableBytes ) )
	{
		bufLen = 0;
		return false;
	}
	
	if( contigI == ALL_CONTIGS )
	{
		uint32 curLen = 0;
		uint64 bytesRead = 0;
		while (curLen < bufLen)
		{
//SeqSeek to start, Figure out how much can be read before SeqSeeking again.
			if(readableBytes <= 0)	//Look out for zero length contigs!  IMPLEMENT ME
				if( !SeqSeek( start + curLen, contigI, startPos, readableBytes ) ){
					bufLen = curLen;
					return true;
				}
			//readLen is the amount to read on this pass
			uint64 readLen = (bufLen - curLen) < readableBytes ? (bufLen - curLen) : readableBytes;	
			gnSeqC* tmpBuf = new gnSeqC[readLen];	//read into tmpBuf, then filter tmpBuf into curBuf

			// read chars and filter
			m_ifstream.read(tmpBuf, readLen);
			uint64 gc = m_ifstream.gcount();
			bytesRead += gc;
			readableBytes -= gc;
			for(uint32 i=0; i < gc; i++){
				if( m_pFilter->IsValid(tmpBuf[i]) ){
					buf[curLen] = tmpBuf[i];
					curLen++;
				}
			}
			delete[] tmpBuf;
			if(m_ifstream.eof()){	//we hit the end of the file.  bail out.
				m_ifstream.clear();
				bufLen = curLen;
				return true;
			}
		}
		bufLen = curLen;
	}
	else if( contigI < m_contigList.size() )
	{
		uint32 curLen = 0;
		//check to see if the buffer is bigger than the contig.  if so truncate it.
		gnSeqI contigSize = m_contigList[contigI]->GetSeqLength();
		bufLen = bufLen < contigSize ? bufLen : contigSize;
		while (curLen < bufLen)
		{
			uint64 readLen = bufLen - curLen;	//the amount to read on this pass
			gnSeqC* tmpBuf = new gnSeqC[readLen];	//read into tmpBuf, then filter tmpBuf into curBuf

			// read chars and filter
			m_ifstream.read(tmpBuf, readLen);
			uint64 gc = m_ifstream.gcount();
//			cout << "Read " << gc << " chars from " << m_openString << "\n";
//			cout << "Checking character validity on: " << tmpBuf << "\n";
			for(uint32 i=0; i < gc; i++){
				if( m_pFilter->IsValid(tmpBuf[i]) ){
					buf[curLen] = tmpBuf[i];
					curLen++;
				}
			}
			if(m_ifstream.eof()){	//we hit the end of the file.  bail out.
				m_ifstream.clear();
				bufLen = curLen;
				return true;
			}
			delete[] tmpBuf;
		}
		bufLen = curLen;
	}
	return true;

}
// private:
// figures out which contig the sequence starts at then calls SeqStartPos to get the offset within that contig
// returns startPos, the file offset where the sequence starts
// returns true if successful, false otherwise
bool gnSEQSource::SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes )
{
	if( contigI == ALL_CONTIGS )
	{
		// find first contig
		gnSeqI curIndex = 0;
		vector< gnFileContig* >::iterator iter = m_contigList.begin();
		for( ; iter != m_contigList.end(); ++iter )
		{
			uint64 len = (*iter)->GetSeqLength();
			if( (curIndex + len) > start )
				break;
			curIndex += len;
		}
		if( iter == m_contigList.end() )
			return false;
		// seek to start
		gnSeqI startIndex = start - curIndex;  //startIndex is starting pos. within the contig
		return SeqStartPos( startIndex, *(*iter), startPos, readableBytes );
	}
	else if( contigI < m_contigList.size() )
	{
		return SeqStartPos( start, *(m_contigList[contigI]), startPos, readableBytes );
	}
	return false;
}
//Returns startPos, the file offset where the sequence starts.
bool gnSEQSource::SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes )
{
	readableBytes = 0;
	uint32 curLen = 0;
	//seek to the file offset where the contig starts
	startPos = contig.GetSectStartEnd(gnContigSequence).first;	//set startPos to start where the contig starts
	m_ifstream.seekg( startPos, ios::beg );
	if( m_ifstream.eof() ){
		DebugMsg("ERROR in gnSEQSource::Incorrect contig start position, End of file reached!\n");
		return false;
	}
	while( true )
	{
		  // READ the rest of the contig skipping over invalid characters until we get to the starting base pair.
		  // startPos will contain the file offset with the starting base pair
		uint32 tmpbufsize = contig.GetSectStartEnd(gnContigSequence).second - startPos;
		if(tmpbufsize == 0){
			DebugMsg("ERROR in gnSEQSource: stored contig size is incorrect.");
			return false;
		}
		uint64 startOffset = start;
		if(contig.HasRepeatSeqGap()){
			if(contig.GetRepeatSeqGapSize().first > 0){
				if(contig.GetRepeatSeqGapSize().second > 0){
					startOffset += (start*contig.GetRepeatSeqGapSize().second)/contig.GetRepeatSeqGapSize().first;
					startPos+=startOffset;
					m_ifstream.seekg(startPos , ios::beg);
					readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
					return true;
				}
			}else{
				startPos+=start;
				m_ifstream.seekg(startPos , ios::beg);
				readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
				return true;
			}
		}
		tmpbufsize = tmpbufsize < BUFFER_SIZE ? tmpbufsize : BUFFER_SIZE;  //read in the smaller of the two.
		char *tmpbuf = new char[tmpbufsize];
		m_ifstream.read( tmpbuf, tmpbufsize );
		if( m_ifstream.eof() ){
			ErrorMsg("ERROR in gnSEQSource::Read End of file reached!\n");
			delete[] tmpbuf;
			return false;
		}
		for( uint32 i=0; i < tmpbufsize; ++i ){
			if( m_pFilter->IsValid(tmpbuf[i]) ){
				if(curLen == start){
					startPos += i;
					readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
					delete[] tmpbuf;
					return true;
				}
				curLen++;
			}
		}
		startPos += tmpbufsize;
		delete[] tmpbuf;
	}
	return false;
}

bool gnSEQSource::ParseStream( std::istream& fin )
{
	if(m_spec == nullptr)
		m_spec = new gnGenomeSpec();
	char *buf = new char[BUFFER_SIZE];	//DNA char buffer
	uint64 bufReadLen = 0;
	std::string::size_type loc = 0;
	uint64 streamPos = 0;
	const char *pRead = NULL;
	char ch, cht;

	gnFileContig *curContig = nullptr;
	gnSourceSpec *curSpec = nullptr;
	gnFeature *curFeature = nullptr;
	gnFragment *curFrag = nullptr;
	boolean readstate = true;

	// Parse the stream...
	while( readstate && (fin.read(buf, BUFFER_SIZE) || fin.gcount() > 0))
	{
		bufReadLen = fin.gcount();

		//	cout << "bufReadLen: " << bufReadLen << "\n";	
		string strBuf(buf, bufReadLen);
		uint32 readState = 0;
		uint32 sectionStart = 0;
		uint32 lineStart = 0;
		gnSeqI seqLength = 0;
		gnSeqI curLocationStart = 0;
		int32 curStartLength = 0;
		int32 curEndLength = 0;
		gnLocation::gnLocationType curBaseLocationType = gnLocation::LT_Standard;
		std::string curLocContig = "";
		std::string curQualifierName = "";
		uint64 curQualifierStart = 0;
		int32 gapstart = -1;
		uint32 lineSeqSize = 0;

		for(uint32 i = 0; i < bufReadLen; i++){
			ch = buf[i];
			switch(readState)
			{
				case 0:  //look for lines that don't start with whitespace or "//", treat them as a tag.  skip the tag and get
					if((ch == '\n') || (ch == '\r')){
						lineStart = i + 1;
						break;
					}
					if((ch == '/') && (i < (bufReadLen - 1)) && (buf[i+1] == '/')){
						readState = 0;
						i++;
						lineStart = i + 1;
						break;
					}
					if((ch != ' ') && (ch != '	')){
						if((i - lineStart) < SEQ_SUBTAG_COLUMN){
							if((i - lineStart) == 0){
								sectionStart = i;
								readState = 1;
							}
						}
					}
					break;
				case 1:	//look for feature tag in column six.  ignore whitespace before feature.
					if((ch == ' ')||(ch == '	')){
						break;
					}else if(ch == '\n'){
						lineStart = i + 1;
						sectionStart = i + 1;
						break;
					}else if(sectionStart == i){ //there was no whitespace, we hit a TAG instead
						i--;
						readState = 0; //Deal with a Header TAG
						sectionStart = i + 1;
						break;
					}else if((i - lineStart == SEQ_SUBTAG_COLUMN)||((buf[lineStart]=='	')&&(i==lineStart+1))){
						sectionStart = i;
						readState = 2;
					} //
				case 2:  //Get the feature name.  stop on whitespace
					if((ch == ' ')||(ch == '	')){
						string featureName(buf+sectionStart, i - sectionStart);
						curFeature = new gnFeature(featureName);
						curFrag->AddFeature(curFeature);
						sectionStart = i + 1;
						readState = 3;
					}
					break;
				case 3:   //Ignore whitespace before feature location
					if((ch == ' ')||(ch == '	')){
						break;
					}else if((ch == '\r')||(ch == '\n')){
						lineStart = i+1;
						break;
					}
					sectionStart = i;
					readState = 4;
				case 4:		//Read a location start.  stop on (<.:^ and whitespace
					if((ch == ' ')||(ch == '	')||(ch == '(')||(ch == '.')||(ch=='^')||(ch==':')){
						string starter(buf+sectionStart, i - sectionStart);
						if(ch == '('){
							if(starter == "complement")
								curFeature->SetLocationType(gnLocation::LT_Complement);
							else if(starter == "order")
								curFeature->SetLocationType(gnLocation::LT_Order);
							else if(starter == "group")
								curFeature->SetLocationType(gnLocation::LT_Group);
							else if(starter == "one-of")
								curFeature->SetLocationType(gnLocation::LT_OneOf);
							sectionStart = i + 1;	//ignore join since it is default.
							break;
						}else if(ch == ':'){
							curLocContig = starter;
							sectionStart = i + 1;
							break;
						}
						curLocationStart = atoi(starter.c_str());
						readState = 6;	//read in end base by default.
						if(ch == '.'){
							//go to special state to look for another one.
							readState = 5;
							sectionStart = i + 1;
							break;
						}else if(ch == '^'){
							curBaseLocationType = gnLocation::LT_BetweenBases;
						}else if((ch == ' ')||(ch == '	')){
							//no end location go to qualifier
							gnLocation curLocation(curLocationStart, curLocationStart);
							curFeature->AddLocation(curLocation, curFeature->GetLocationListLength());
							readState = 7;
						}
						sectionStart = i + 1;

					}else if(ch == '<'){
						curStartLength = -1;
						sectionStart = i + 1;
					}else if(ch == '>'){
						curStartLength = 1;
						sectionStart = i + 1;
					}
					break;
				case 5: //look for another period or location start.
					if(ch == '.'){
						curBaseLocationType = gnLocation::LT_Standard;
						readState = 6;
						sectionStart = i + 1;
						break;
					}
					curBaseLocationType = gnLocation::LT_OneOf;
				case 6:	//see if there's a second location value.  stop on >, and whitespace
					if(ch == '>'){
						curEndLength = 1;
						sectionStart = i + 1;
					}else if(ch == '<'){
						curEndLength = -1;
						sectionStart = i + 1;
					}else if((ch == ' ')||(ch == '	')||(ch == ',')){
						//read end location
						string ender(buf+sectionStart, i - sectionStart);
						gnSeqI curLocationEnd = atoi(ender.c_str());
						gnLocation curLocation(curLocationStart, curStartLength, curLocationEnd, curEndLength, curBaseLocationType);
						curEndLength = 0;
						curStartLength = 0;
						curFeature->AddLocation(curLocation, curFeature->GetLocationListLength());
						readState = ch == ',' ? 3 : 7;  //read another loc if we need to.
						sectionStart = i+1;
					}
					break;
				case 7:  //skip to start of qualifier
					if((ch != ' ')&&(ch != '	')&&(lineStart == i)){
						sectionStart = i;	// Hit a tag.  go deal with it.
						readState = 0;
						i--;
					}else if((ch != ' ')&&(ch != '	')&&((lineStart == i - SEQ_SUBTAG_COLUMN)||((buf[lineStart]=='	')&&(i==lineStart+1)))){
						sectionStart = i;	// Hit a feature.  go deal with it.
						readState = 2;
						i--;
					}else if(ch == ','){  //oops!  another location to read!
						sectionStart = i+1;
						readState = 3;
					}else if(ch == '/'){  //finally, a qualifier.
						sectionStart = i+1;
						readState = 8;
					}else if(ch == '\n')
						lineStart = i + 1;
					break;
				case 8:		//get a qualifier, stop on =
					if(ch == '='){
						curQualifierName = string(buf+sectionStart, i - sectionStart);
						readState = 9;
						sectionStart = i+1;
					}
					break;
				case 9:		//are we getting a string? look for " or [
					if(ch == '"'){
						readState = 10;
						sectionStart = i;
						curQualifierStart = i + streamPos;
					}else if(ch == '['){
						readState = 11;
						sectionStart = i;
					}else if((ch == '\r')||(ch == '\n')){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, sectionStart + streamPos, i - sectionStart));
						sectionStart = i+1;
						readState = 7; //look for another qualifier
					}
					break;
				case 10:		//read until the end of the quotation. look out for escaped quotes
					if(ch == '"')
						readState = 11;
					if(ch == '\n'){
						lineStart = i + 1;
					}
					break;
				case 11:
					if(ch != '"'){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, curQualifierStart, i - sectionStart));
						sectionStart = i+1;
						readState = 7;	//look for another qualifier.
						if(ch == '\n')
							lineStart = i + 1;
					}else
						readState = 10;  //quote was escaped.  look for another.
					break;
				case 12:
					if(ch == ']'){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, sectionStart + streamPos, i - sectionStart));
						sectionStart = i+1;
						readState = 7;	//look for another qualifier.
					}
					break;
				case 13:	//start the sequence read.
					if(ch == '^')	//stupid blattlab file format.
						readState = 14;
					else{
						curContig->SetSectStart(gnContigSequence, i + 1 + streamPos);
						readState = 16;
					}
					break;
				case 14:	//wait for newline before sequence starts.
					if(ch == '\n'){
						curContig->SetRepeatSeqGap(true);
						lineStart = i + 1;
						sectionStart = i + 1;
						curContig->SetSectStart(gnContigSequence, i + 1 + streamPos);
						readState = 15;
					}
					break;
				case 15:
					if(m_pFilter->IsValid(ch))
						seqLength++;
					else
						curContig->SetRepeatSeqGap(false);
					break;
				case 16:
					if((ch == '/')&&(i==lineStart)){
						readState = 17;
					}else if(m_pFilter->IsValid(ch)){
						seqLength++;
						lineSeqSize++;
						if(gapstart >= 0){
							curContig->SetRepeatGapSize(i - gapstart);
							gapstart = -1;
						}
					}else if(ch == '\n'){	//IMPLEMENT ME! Needs consistent gap size checking
						if(sectionStart == lineStart){
							curContig->SetRepeatSeqGap(true);
							curContig->SetRepeatSeqSize(seqLength);
							gapstart = i;
							for(; gapstart >= lineStart; gapstart--)
								if(m_pFilter->IsValid(buf[gapstart]))
									break;
							gapstart++;
						}else if(lineSeqSize != curContig->GetRepeatSeqGapSize().first)
							curContig->SetRepeatSeqGap(false);
						lineSeqSize = 0;
						lineStart = i + 1;
					}
					break;
				case 17:
					if((ch == '\n')&&(buf[lineStart+1] == '/')){
						curContig->SetSectEnd(gnContigSequence, lineStart - 2 + streamPos);
						curContig->SetSeqLength(seqLength);
						m_contigList.push_back(curContig);
						curContig = nullptr;
						curSpec->SetLength(seqLength);
						curSpec = nullptr;
						seqLength = 0;
						lineStart = i + 1;
						sectionStart = i + 1;
						readState = 0;
					}
					break;
			}
		}
		streamPos += bufReadLen;
	}
	if(curContig != nullptr){
		curContig->SetSectEnd(gnContigSequence, streamPos - 1);
		curContig->SetSeqLength(seqLength);
		m_contigList.push_back(curContig);
		curSpec->SetLength(seqLength);
	}
	if(curSpec != nullptr)
		if((curFrag->GetFeatureListLength() == 0) && (curFrag->GetHeaderListLength() == 0)
			&&(curSpec->GetLength() == 0)){
			m_spec->RemoveSpec(m_spec->GetSpecListLength() - 1);
			delete curFrag;
		}
	m_ifstream.clear();
	delete[] buf;
	return true;
}

}	// end namespace genome
