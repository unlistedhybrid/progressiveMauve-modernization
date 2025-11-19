/////////////////////////////////////////////////////////////////////////////
// File:            gnRAWSequence.h
// Purpose:         Optimized Sequence class for RAW sequence files
// Description:     Provides a high level sequence interface to all types of
//					sequence data.
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

#ifndef _gnRAWSequence_h_
#define _gnRAWSequence_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <iostream>
#include <list>
#include "libGenome/gnSequence.h"

namespace genome {


/**
 * gnRAWSequence is a bastardization of gnSequence that creates a lightweight wrapper
 * around a memory-mapped file of raw sequence data
 */
class GNDLLEXPORT gnRAWSequence : public gnSequence
{
public:
	/**
	 * Empty Constructor, creates an empty gnRAWSequence.
	 */
	gnRAWSequence();
	/**
	 * Creates a gnRAWSequence based on the file specified by filename
	 */
	gnRAWSequence( const std::string& filename )
	{
		this->filename = filename;
		data.open( filename );
	}

 virtual gnRAWSequence* Clone() const override {return new gnRAWSequence(*this);}

	virtual gnSeqI contigListSize() const {return 1;}
	virtual gnSeqI contigListLength() const {return 1;}
	virtual uint32 contigIndexByBase( const gnSeqI baseI) const {
		if(baseI >= data.size()) Throw_gnEx(SeqIndexOutOfBounds());
		return 0;
	}
	gnRAWSequence contig( const uint32 contigI) const { 
		if(contigI>0) Throw_gnEx(FragmentIndexOutOfBounds()); 
		return *this;
	}
	gnRAWSequence contigByBase( const gnSeqI baseI) const {
		if(baseI >= data.size()) Throw_gnEx(SeqIndexOutOfBounds());
		return *this;
	}
 virtual gnSeqI contigStart( const uint32 contigI) const override {
		if(contigI>0) Throw_gnEx(FragmentIndexOutOfBounds()); 
		return 0;
	}
 virtual gnSeqI contigLength( const uint32 contigI) const override {
		if(contigI>0) Throw_gnEx(FragmentIndexOutOfBounds()); 
		return data.size();
	}
 virtual uint32 contigIndexByName( std::string& contigName) const override {
		return 0;
	}
 virtual std::string contigName( const uint32 contigI) const override {
		if(contigI>0) Throw_gnEx(FragmentIndexOutOfBounds()); 
		return "";
	}
 virtual gnSequence contigByName( std::string& contigName) const override {
		Throw_gnEx(FragmentIndexOutOfBounds());
	}
 virtual void merge([[maybe_unused]] const gnSeqI startI, [[maybe_unused]] const gnSeqI endI) override { throw; }
 virtual void mergeContigs([[maybe_unused]] const uint32 startC, [[maybe_unused]] const uint32 endC) override { throw; }
 virtual void splitContig([[maybe_unused]] const gnSeqI splitI, [[maybe_unused]] const uint32 contigI=ALL_CONTIGS) override { throw; }

 virtual void setContigName( [[maybe_unused]] const uint32 contigI, [[maybe_unused]] const std::string& contig_name) override { throw; }

 virtual uint32 getFeatureListLength() const override {
		return 0;
	}
 virtual gnBaseFeature* getFeature([[maybe_unused]] const uint32 featureI) const override { Throw_gnEx(FeatureIndexOutOfBounds()); }
 virtual void getContainedFeatures([[maybe_unused]] const gnLocation& lt, [[maybe_unused]] std::vector<gnBaseFeature*>& feature_vector, [[maybe_unused]] std::vector<uint32>& index_vector) const override {}
 virtual void getIntersectingFeatures([[maybe_unused]] const gnLocation& lt, [[maybe_unused]] std::vector<gnBaseFeature*>& feature_vector, [[maybe_unused]] std::vector<uint32>& index_vector) const override {}
 virtual uint32 addFeature([[maybe_unused]] gnBaseFeature* feature) override { throw; }
 virtual void removeFeature([[maybe_unused]] const uint32 featureI) override { Throw_gnEx(FeatureIndexOutOfBounds()); }
 virtual void getBrokenFeatures([[maybe_unused]] const gnLocation& lt, [[maybe_unused]] std::vector<gnBaseFeature*>& feature_vector) const override {};
 virtual uint32 getHeaderListLength([[maybe_unused]] const uint32 contigI) const override { return 0; }
 virtual gnBaseHeader* getHeader([[maybe_unused]] const uint32 contigI, [[maybe_unused]] const uint32 headerI) const override {Throw_gnEx(HeaderIndexOutOfBounds());};
 virtual void addHeader([[maybe_unused]] const uint32 contigI, [[maybe_unused]] gnBaseHeader* header, [[maybe_unused]] const uint32 headerI) override {Throw_gnEx(FragmentIndexOutOfBounds());}
 virtual void removeHeader([[maybe_unused]] const uint32 contigI, [[maybe_unused]] const uint32 headerI) override { Throw_gnEx(HeaderIndexOutOfBounds()); }
 virtual void setReverseComplement( const bool revComp, const uint32 contigI=ALL_CONTIGS) override {Throw_gnEx(FragmentIndexOutOfBounds());};
 virtual bool isReverseComplement( const uint32 contigI=ALL_CONTIGS ) override {return false;}
 virtual bool isCircular() const override { return false; }
 virtual void setCircular( const bool value ) override {}
	
 virtual void globalToLocal(uint32& contigI, gnSeqI& baseI) const override {};
 virtual void localToGlobal(const uint32 contigI, gnSeqI& baseI) const override {};
 virtual void globalToSource(uint32& contigI, gnSeqI& baseI) const override {};
 virtual void localToSource(uint32& contigI, [[maybe_unused]] gnSeqI& baseI) const override {};
 [[nodiscard]] virtual bool LoadSource(const std::string sourcename) override {
		data.open( sourcename );
		filename = sourcename;
		return true;
	}

	/**
	 * Appends the bases in "seq" to this sequence.
	 */
	gnRAWSequence& operator+=(const gnRAWSequence& seq);

	/**
	 * Compares the bases in "seq" to this sequence.
	 * @param seq The sequence to compare this sequence to.
	 * @return Negative if this sequence is lesser, 0 if the two sequences are
	 * equal, and positive if this sequence is greater.
	 */
/*	virtual int compare(const gnRAWSequence& seq) const;
	virtual int compare(const std::string& str) const;

	virtual void append( const gnRAWSequence& seq);
	virtual void insert( const gnSeqI offset, const gnSeqC *bases, const gnSeqI length);
	virtual void insert( const gnSeqI offset, const gnRAWSequence& seq);
	virtual void insert( const gnSeqI offset, const gnGenomeSpec& gnbs);
	gnRAWSequence const operator+(const gnRAWSequence& seq) const;
	virtual void erase( const gnSeqI offset=0, const gnSeqI length=GNSEQI_END );
*/
	gnRAWSequence subseq(const gnSeqI offset, const gnSeqI length) const
	{
		gnRAWSequence gnrs;
		gnrs.data.open(filename, length, offset - 1);
		return gnrs;
	}
//	friend std::istream& operator>>(std::istream& is, gnRAWSequence& gns);	//read from source.
	/**
	 * Writes the bases in this sequence to the specified output stream (e.g. cout).
	 */
//	friend std::ostream& operator<<(std::ostream& os, const gnRAWSequence& gns); //write to source.

 virtual gnSeqI length() const override { return data.size(); }
 virtual gnSeqI size() const override { return data.size(); }

 virtual std::string ToString( const gnSeqI length=GNSEQI_END, const gnSeqI offset=1 ) const override
	{
		gnSeqI len = length == GNSEQI_END ?  data.size() - offset - 1 : length;
		std::string asdf(data.data()+offset-1, len);
		return asdf;
	}

 [[nodiscard]] virtual bool ToString( std::string& str, const gnSeqI length=GNSEQI_END, const gnSeqI offset=1 ) const override
	{
		gnSeqI len = length == GNSEQI_END ? data.size() - offset - 1 : length;
		str.assign(data.data()+offset-1,len);
		return true;
	}
 [[nodiscard]] virtual bool ToArray( gnSeqC* pSeqC, gnSeqI length, const gnSeqI offset=1 ) const override
	{
		gnSeqI len = length == GNSEQI_END ? data.size() - offset - 1 : length;
		memcpy(pSeqC, data.data()+offset-1, len);
		return true;
	}
 virtual gnSeqC GetSeqC( const gnSeqI offset ) const override
	{
		return *(data.data()+(offset-1));
	}
	gnSeqC operator[]( const gnSeqI offset ) const
	{
		return *(data.data()+(offset-1));
	}

	virtual gnSeqI find([[maybe_unused]] const gnSequence& search, [[maybe_unused]] const gnSeqI offset=0) const override {return GNSEQI_ERROR;}
	virtual gnSeqI find([[maybe_unused]] const gnRAWSequence& search, [[maybe_unused]] const gnSeqI offset=0) const {return GNSEQI_ERROR;}
	
private:
	boost::iostreams::mapped_file_source data;
	std::string filename;
}; // class gnRAWSequence

/*
GNDLLEXPORT
std::istream& operator>>(std::istream& is, gnRAWSequence& gns);	//read from source.
GNDLLEXPORT
std::ostream& operator<<(std::ostream& os, const gnRAWSequence& gns); //write to source.
*/



}	// end namespace genome

#endif
	// _gnRAWSequence_h_
