/*******************************************************************************
 * $Id: MatchFinder.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MatchFinder_h_
#define _MatchFinder_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SortedMerList.h"
#include "libMems/Match.h"
#include "libMems/MatchList.h"
#include <list>
#include <iostream>
#include <cstdint>
#include <boost/pool/pool_alloc.hpp>

namespace mems {

struct idmer{
	gnSeqI	position;	//starting position of this mer in the genome
	uint64_t 	mer; 		//the actual sequence
	sarID_t	id;			//the sequence identifier.
};

// typedef std::list<idmer, boost::fast_pool_allocator<idmer> > IdmerList;
// using boost::fast_pool_allocator<idmer> results in a significant speedup
// over std::allocator.  testing on a Salmonella vs. Y. pestis comparison shows
// a 30% speedup
typedef std::list<idmer> IdmerList;

const unsigned int PROGRESS_GRANULARITY = 100;

/**
 * This pure virtual class implements a general framework for finding
 * exactly matching mers.  It is extended by the MemHash and MemScorer
 * classes.
 * @see MemHash
 * @see MemScorer
 */
class MatchFinder : public genome::gnClone{
public:
	MatchFinder();
	~MatchFinder();
	MatchFinder(const MatchFinder& mf);
	virtual void Clear();
	/**
	 * Adds a sequence to use when searching for exact matches.
	 * @param sar A pointer to the sorted mer list for the new sequence
	 * @param seq A pointer to the genome::gnSequence corresponding to the new sequence.
	 */
	virtual bool AddSequence( SortedMerList* sar, genome::gnSequence* seq = NULL );
	/**
	 * Given the index of a sequence and an index into the sorted mer list, this function
	 * will search the other sorted mer lists for the same mer.  This function returns the
	 * position of the mer in each sequence in the breakpoints vector.
	 */
	virtual void GetBreakpoint( uint32_t sarI, gnSeqI startI, std::vector<gnSeqI>& breakpoints ) const;
	virtual uint32_t Multiplicity(void){return seq_count;};
	/** NOT IMPLEMENTED: Sets the number of ambiguities allowed in a mer match*/
	virtual void SetAmbiguityTolerance(uint32_t ambiguity_tol){ambiguity_tolerance = ambiguity_tol;}
	/** @return the number of ambiguities allowed in a mer match */
	virtual uint32_t AmbiguityTolerance(){return ambiguity_tolerance;}
	/** @return The progress of the current operation.  Ranges from 0 to 100.  -1 indicates no computation is being performed */
	virtual float GetProgress() const {return m_progress;}

	/** Finds all the matches between the sequences */
	virtual void FindMatchSeeds();
	/** Finds all the matches between the sequences, starting at a particular offset */
	virtual void FindMatchSeeds( const std::vector<gnSeqI>& start_offsets );

	/**
	 * Logs progress to the designated ostream.  Set to null to skip progress logging.
	 */
	virtual void LogProgress( std::ostream* os );
	void SetOffsetLog( std::ostream* offset_stream ){ this->offset_stream = offset_stream; }
protected:
	/** 
	 * Searches for mer matches in a designated range of the sequence's sorted mer lists 
	 * @throws InvalidData thrown if the start_points are bad or if the sorted mer lists were sorted on different mer sizes
	 * @return true if completed searching, false if repetitive mers were encountered and FindMatches must be called again.
	 */
	[[nodiscard]] virtual bool SearchRange(std::vector<gnSeqI>& start_points, std::vector<gnSeqI>& search_len);
	/** Called whenever a mer match is found */
	[[nodiscard]] virtual bool HashMatch(IdmerList& match_list) = 0;
	[[nodiscard]] virtual bool EnumerateMatches(IdmerList& match_list);

	template< class MatchType >
	void FindSubsets(const MatchType& mhe, std::vector<MatchType>& subset_matches);

	template< class UngappedMatchType >
	void ExtendMatch(UngappedMatchType& mhe, std::vector<UngappedMatchType>& subset_matches, gnSeqI max_backward = GNSEQI_END, gnSeqI max_forward = GNSEQI_END);

	virtual SortedMerList* GetSar(uint32_t sarI) const;
	std::vector<SortedMerList*> sar_table;
	std::vector<genome::gnSequence*> seq_table;
	
	uint32_t mer_size;
	uint32_t seq_count;
	uint32_t ambiguity_tolerance;
	
	// for subset matches
	std::vector< std::vector< uint32_t > > alpha_map;
	uint32_t alpha_map_size;
	uint32_t alphabet_bits;
	
	float m_progress;
	std::ostream* log_stream;

	uint64_t mers_processed;	/**< The number of mers processed thus far */
	uint64_t total_mers;	/**< The total number of mers to search */
	std::ostream* offset_stream;	/**< log for the current offset in each SML */
};

/** 
 * InvalidData exceptions are thrown when the input to an algorithm is invalid
 */
CREATE_EXCEPTION( InvalidData );

inline
SortedMerList* MatchFinder::GetSar(uint32_t sarI) const{
	if(sarI >= sar_table.size()){
		std::cerr << "ERROR: GetSar() out of bounds! sarI=" << sarI << " sar_table.size()=" << sar_table.size() << "\n";
		return NULL;
	}
	return sar_table[sarI];
}

inline
bool idmer_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

//id less than function for STL sort functions
inline
bool idmer_id_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.id < m_v.id);// ? true : false;
};



// takes as input a fully extended mem and returns the subset matches on the lower side
template< class MatchType >
void MatchFinder::FindSubsets(const MatchType& mhe, std::vector<MatchType>& subset_matches){

	SMLHeader head = GetSar( 0 )->GetHeader();
	uint32_t shift_amt = 64 - head.alphabet_bits;
	uint32_t rshift_amt = head.alphabet_bits * ( GetSar(0)->SeedLength() - 1 );

	uint32_t seqI, alphaI;

	// initialize subset match data structures
	alpha_map_size = 1;
	alpha_map_size <<= alphabet_bits;
	if( alpha_map.size() != alpha_map_size ){
		alpha_map.clear();
		alpha_map.reserve( alpha_map_size );
		std::vector< uint32_t > tmp_list;
		tmp_list.reserve( seq_count );
		for( uint32_t alphaI = 0; alphaI < alpha_map_size; ++alphaI )
			alpha_map.push_back( tmp_list );
	}else{
		for( uint32_t alphaI = 0; alphaI < alpha_map_size; ++alphaI )
			alpha_map[ alphaI ].clear();
	}
	
	
	for( seqI = 0; seqI < seq_count; ++seqI ){
		//check that all mers at the new position match
		int64_t mer_to_get = mhe[ seqI ];
		if( mer_to_get == NO_MATCH )
			continue;
		if(mer_to_get < 0){
			mer_to_get *= -1;
			mer_to_get += mhe.Length() - GetSar(0)->SeedLength();
		}

		uint64_t cur_mer = GetSar( seqI )->GetMer( mer_to_get - 1 );

		bool parity;
		if( mhe[ seqI ] < 0 )
			parity = cur_mer & 0x1;
		else
			parity = !(cur_mer & 0x1);

		if( parity ){
			cur_mer >>= shift_amt;
		}else{
			cur_mer <<= rshift_amt;
			cur_mer = ~cur_mer;
			cur_mer >>= shift_amt;
		}

		alpha_map[ cur_mer ].push_back( seqI );

	}
	
	for( alphaI = 0; alphaI < alpha_map_size; ++alphaI ){
		if( alpha_map[ alphaI ].size() < 2 ){
			alpha_map[ alphaI ].clear();
			continue;
		}
		// this is a subset
		MatchType cur_subset = mhe;
		cur_subset.SetLength( mhe.Length() );
		for( uint32_t sqI = 0; sqI < mhe.SeqCount(); ++sqI )
			cur_subset.SetStart( sqI, NO_MATCH );	// init everything to NO_MATCH
		for( uint32_t subI = 0; subI < alpha_map[ alphaI ].size(); ++subI )
			cur_subset.SetStart( alpha_map[ alphaI ][ subI ], mhe[ alpha_map[ alphaI ][ subI ] ] );
		subset_matches.push_back( cur_subset );
		alpha_map[ alphaI ].clear();
	}
}

// BUGS:
// matches which span the end-start of a circular sequence will be hashed a second time
template< class UngappedMatchType >
void MatchFinder::ExtendMatch(UngappedMatchType& mhe, std::vector<UngappedMatchType>& subset_matches, gnSeqI max_backward, gnSeqI max_forward){
	uint64_t cur_mer;
	uint64_t mer_mask = GetSar(0)->GetSeedMask();

	//which sequences are used in this match?
	if(GetSar(0) == NULL){
		std::cerr << "ERROR: GetSar(0) is NULL in ExtendMatch!\n";
		return;
	}
	
	uint32_t* cur_seqs = new uint32_t[mhe.SeqCount()];
	uint32_t used_seqs = 0;
	for(uint32_t seqI = 0; seqI < mhe.SeqCount(); ++seqI){
		if(seqI >= seq_count){
			std::cerr << "ERROR: seqI " << seqI << " >= seq_count " << seq_count << "\n";
			delete[] cur_seqs;
			return;
		}
		if(mhe[seqI] != NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			++used_seqs;
		}
	}
	//First extend backwards then extend forwards.  The following loop does them both.
	int jump_size = GetSar(0)->SeedLength();
	uint32_t extend_limit = 0;	/**< Tracks the distance to the most distant overlapping matching seed */
	uint32_t extend_attempts = 0;	/**< Counts the total number of overlapping seeds checked */
	bool extend_again = false;	/**< Set to true if any overlapping seeds matched, the search will be restarted from that point */
	for(uint32_t directionI = 0; directionI < 4; ++directionI){
		//how far can we go?	
		//first calculate the maximum amount of traversal
		//then do fewer comparisons.
		int64_t maxlen = GNSEQI_END;
		if(directionI == 0)
			maxlen = static_cast<int64_t>(max_backward);
		else if(directionI == 1)
			maxlen = static_cast<int64_t>(max_forward);
		else
			maxlen = static_cast<int64_t>(GetSar(0)->SeedLength());
		for(uint32_t maxI = 0; maxI < used_seqs; ++maxI)
			if(GetSar(cur_seqs[maxI])->IsCircular()){
				if(static_cast<int64_t>(GetSar(cur_seqs[maxI])->Length()) < maxlen)
					maxlen = static_cast<int64_t>(GetSar(cur_seqs[maxI])->Length());
			}else if(mhe[cur_seqs[maxI]] < 0){
				int64_t rc_len = static_cast<int64_t>(GetSar(cur_seqs[maxI])->Length()) - static_cast<int64_t>(mhe.Length()) + mhe[cur_seqs[maxI]] + 1;
				if( rc_len < maxlen)
					maxlen = rc_len;
			}else if(mhe[cur_seqs[maxI]] - 1 < maxlen)
				maxlen = mhe[cur_seqs[maxI]] - 1;
		uint32_t j=0;
		uint32_t i = used_seqs;	// set to used_seqs in case maxlen is already less than jump size.

		extend_limit = 0;
		extend_attempts = 0;

		while(maxlen >= static_cast<int64_t>(jump_size)){
			mhe.SetLength(mhe.Length() + jump_size);
			maxlen -= jump_size;
			for(j=0; j < used_seqs; ++j){
				if(cur_seqs[j] >= seq_count){
					std::cerr << "ERROR: cur_seqs[" << j << "]=" << cur_seqs[j] << " >= seq_count=" << seq_count << "\n";
					delete[] cur_seqs;
					return;
				}
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - jump_size);
					if(mhe[cur_seqs[j]] <= 0)
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + static_cast<int64_t>(GetSar(cur_seqs[j])->Length()));
				}
			}
			//check that all mers at the new position match
			if(used_seqs == 0){
				std::cerr << "ERROR: used_seqs is 0 in ExtendMatch loop\n";
				delete[] cur_seqs;
				return;
			}
			if(cur_seqs[0] >= seq_count){
				std::cerr << "ERROR: cur_seqs[0]=" << cur_seqs[0] << " >= seq_count=" << seq_count << "\n";
				delete[] cur_seqs;
				return;
			}
			int64_t mer_to_get = mhe[cur_seqs[0]];
			if(mer_to_get < 0){
				mer_to_get *= -1;
				mer_to_get += static_cast<int64_t>(mhe.Length()) - static_cast<int64_t>(GetSar(0)->SeedLength());
			}
			cur_mer = GetSar(cur_seqs[0])->GetSeedMer(mer_to_get - 1);
			bool parity;
			if( mhe[cur_seqs[0]] < 0 )
				parity = cur_mer & 0x1;
			else
				parity = !(cur_mer & 0x1);
			cur_mer &= mer_mask;

			for(i=1; i < used_seqs; ++i){
				if(cur_seqs[i] >= seq_count){
					std::cerr << "ERROR: cur_seqs[" << i << "]=" << cur_seqs[i] << " >= seq_count=" << seq_count << "\n";
					delete[] cur_seqs;
					return;
				}
				mer_to_get = mhe[cur_seqs[i]];
				if(mer_to_get < 0){
					//Convert the cur_seqs[i] entry since negative implies reverse complement
					mer_to_get *= -1;
					mer_to_get += static_cast<int64_t>(mhe.Length()) - static_cast<int64_t>(GetSar(0)->SeedLength());
				}
				uint64_t comp_mer = GetSar(cur_seqs[i])->GetSeedMer(mer_to_get - 1);
				bool comp_parity;				
				if( mhe[cur_seqs[i]] < 0 )
					comp_parity = comp_mer & 0x1;
				else
					comp_parity = !(comp_mer & 0x1);
				comp_mer &= mer_mask;
				
				if(cur_mer != comp_mer || parity != comp_parity ){
					if( directionI < 2 )
						maxlen = 0;
					break;
				}
			}
			extend_attempts += jump_size;
			if( i == used_seqs )
				extend_limit = extend_attempts;
			if( directionI > 1 && extend_attempts == static_cast<uint32_t>(GetSar(0)->SeedLength()) )
				break;
		}
		//this stuff cleans up if there was a mismatch
		if(i < used_seqs){
			mhe.SetLength(mhe.Length() - jump_size);
			for(;j > 0; j--){
				if(cur_seqs[j-1] >= seq_count){
					std::cerr << "ERROR: cur_seqs[" << (j-1) << "]=" << cur_seqs[j-1] << " >= seq_count=" << seq_count << "\n";
					delete[] cur_seqs;
					return;
				}
				if(mhe[cur_seqs[j - 1]] >= 0)
					mhe.SetStart(cur_seqs[j - 1], mhe[cur_seqs[j - 1]] + jump_size);
			}
		}
		// check whether any of the overlapping seeds matched.
		// if so, set the match to that length and set the flag to start the search again
		if( directionI > 1 && extend_attempts > 0 ){
			if( extend_limit > 0 )
				extend_again = true;
			// minus jump_size because the cleanup above already moved the length back a little
			int unmatched_diff = static_cast<int>(extend_attempts) - static_cast<int>(extend_limit);
			if( i < used_seqs )
				unmatched_diff -= jump_size;
			if( (unmatched_diff > static_cast<int>(mhe.Length())) && unmatched_diff >= 0 )
				std::cerr << "oh sheethockey mushrooms\n";
			mhe.SetLength(mhe.Length() - unmatched_diff);
			for(j=0; j < used_seqs; ++j){
				if(cur_seqs[j] >= seq_count){
					std::cerr << "ERROR: cur_seqs[" << j << "]=" << cur_seqs[j] << " >= seq_count=" << seq_count << "\n";
					delete[] cur_seqs;
					return;
				}
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + unmatched_diff);
					if(mhe[cur_seqs[j]] > static_cast<int64_t>(GetSar(cur_seqs[j])->Length()) )
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - static_cast<int64_t>(GetSar(cur_seqs[j])->Length()) );
				}
			}
		}
		//Invert the sequence directions so that we extend in the other direction
		//next time through the loop.  The second time we do this we are setting
		//sequence directions back to normal.
		mhe.Invert();

		//if we've already been through twice then decrease the jump size
		if(directionI >= 1)
			jump_size = 1;
		if( directionI == 3 && extend_again ){
			directionI = -1;	// will become 0 on next iteration
			jump_size = GetSar(0)->SeedLength();
			extend_again = false;
		}
	}
	// after the match has been fully extended, search for subset matches
	// this code only works when using SOLID seeds-- so it's been disabled
/*	if( used_seqs > 2 ){
		FindSubsets( mhe, subset_matches );
		mhe.Invert();
		FindSubsets( mhe, subset_matches );
		mhe.Invert();
	}
*/
	// set the subsets so their reference sequence is always positive
	for(uint32_t subsetI = 0; subsetI < subset_matches.size(); ++subsetI){
		if( subset_matches[subsetI][subset_matches[subsetI].FirstStart()] < 0 )
			subset_matches[subsetI].Invert();
		subset_matches[subsetI].CalculateOffset();
	}

	delete[] cur_seqs;
}



}

#endif	//_MatchFinder_h_
