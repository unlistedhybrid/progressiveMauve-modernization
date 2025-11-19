/*******************************************************************************
 * $Id: progressiveAligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * BEWARE!!
 * This code was created in the likeness of the flying spaghetti monster
 *
 * dedicated to Loren...
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libMems/ProgressiveAligner.h"
#include "libMems/GreedyBreakpointElimination.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libMems/gnAlignedSequences.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/PairwiseMatchFinder.h"
#include "libMems/TreeUtilities.h"
#include "libMems/PairwiseMatchAdapter.h"
#include "libMems/DistanceMatrix.h"

#include <boost/dynamic_bitset.hpp>
#include <tuple>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>
#include <iomanip>

#include "stdlib.h"
#include <cstdint>

using namespace std;
using namespace genome;

namespace mems {


bool progress_msgs = false;

bool debug_me = false;
static int dbg_count = 0; 	 


double min_window_size = 200;
double max_window_size = 20000;  // don't feed MUSCLE anything bigger than this
double min_density = .5;
double max_density = .9;
size_t max_gap_length = 5000;
size_t lcb_hangover = 300;


void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list );

/**
 * Test code to ensure that an individual LCB is truly collinear
 * @return	true if the LCB is good
 */
bool my_validateLCB( MatchList& lcb ){
	vector< Match* >::iterator lcb_iter = lcb.begin();
	if( lcb.size() == 0 )
		return true;
	uint seq_count = (*lcb_iter)->SeqCount();
	uint seqI = 0;
	bool complain = false;
	for(; seqI < seq_count; seqI++ ){
		lcb_iter = lcb.begin();
		int64 prev_coord = 0;
		for(; lcb_iter != lcb.end(); ++lcb_iter ){
			if( (*lcb_iter)->Start( seqI ) == NO_MATCH )
				continue;
			else if( prev_coord != 0 && (*lcb_iter)->Start( seqI ) < prev_coord ){
				complain = true;
			}
			prev_coord = (*lcb_iter)->Start( seqI );
		}
	}
	return !complain;
}

template< class BoostMatType >
void print2d_matrix( BoostMatType& mat, std::ostream& os )
{
	for( size_t i = 0; i < mat.shape()[0]; ++i )
	{
		for( size_t j = 0; j < mat.shape()[1]; ++j )
		{
			if( j > 0 )
				os << "\t";
			os << mat[i][j];
		}
		os << endl;
	}
}

double getDefaultBreakpointPenalty( std::vector< gnSequence* >& sequences )
{
	uint default_mer_size = MatchList::GetDefaultMerSize( sequences );
	double avg_seq_len = 0;
	for( size_t seqI = 0; seqI < sequences.size(); ++seqI )
		avg_seq_len += (double)sequences[seqI]->length();
	avg_seq_len /= (double)sequences.size();
	avg_seq_len = log( avg_seq_len ) / log( 2.0 );
	return avg_seq_len * 7000;	  // seems to work reasonably well?
}


double getDefaultBpDistEstimateMinScore( std::vector< gnSequence* >& sequences )
{
	// this value was empirically derived by a process that involved burning incense
	// and uttering arcane words
	return 3.0 * getDefaultBreakpointPenalty(sequences);
}



/*
 * A progressive alignment algorithm for genomes with rearrangements.
 * Start simple, add complexity later.
 * TODO: rewrite the algorithm outline
 */

ProgressiveAligner::ProgressiveAligner( uint seq_count ) :
Aligner( seq_count ),
breakpoint_penalty( -1 ),
min_breakpoint_penalty( 4000 ),
debug(false),
refine(true),
scoring_scheme(ExtantSumOfPairsScoring),
use_weight_scaling(true),
conservation_dist_scale(1),
bp_dist_scale(.9),
max_gapped_alignment_length(20000),
bp_dist_estimate_score(-1),
use_seed_families(false),
using_cache_db(true)
{
	gapped_alignment = true;
	max_window_size = max_gapped_alignment_length;
}

void ProgressiveAligner::SetMaxGappedAlignmentLength( size_t len )
{ 
	max_gapped_alignment_length = len; 
	max_window_size = max_gapped_alignment_length;
}

/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getAlignedChildren( node_id_t node, vector< node_id_t >& descendants )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( node );
	vector< bool > visited( alignment_tree.size(), false );
	descendants.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		if(progress_msgs) cout << "Evaluating aligned nodes linked to node " << cur_node << endl;
		node_stack.pop();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].children_aligned[childI] && !visited[child_id])
				node_stack.push( child_id );
		}
		if( alignment_tree[ cur_node ].sequence != NULL )
			descendants.push_back( cur_node );
	}
}


/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getPath( node_id_t first_n, node_id_t last_n, vector< node_id_t >& path )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( last_n );
	vector< bool > visited( alignment_tree.size(), false );
	while( node_stack.top() != first_n )
	{
		node_id_t cur_node = node_stack.top();
		size_t pre_size = node_stack.size();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if(!visited[child_id])
			{
				node_stack.push( child_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		for( uint parentI = 0; parentI < alignment_tree[cur_node].parents.size(); parentI++ )
		{
			node_id_t parent_id = alignment_tree[cur_node].parents[parentI];
			if(!visited[parent_id])
			{
				node_stack.push( parent_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		node_stack.pop();	// didn't make any progress
	}
	path = vector< node_id_t >( node_stack.size() );
	for( size_t pI = 0; pI < path.size(); pI++ )
	{
		path[pI] = node_stack.top();
		node_stack.pop();
	}
}






template<class MatchType>
void ProgressiveAligner::propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector<MatchType*>& iv_list )
{
	SSC<MatchType> ilc(seqI);
	sort( iv_list.begin(), iv_list.end(), ilc );
	vector< SuperInterval >& ord = alignment_tree[ node1 ].ordering;
	vector<gnSeqI> bp_list;
	for( size_t sI = 0; sI < ord.size(); sI++ )
		bp_list.push_back( ord[sI].LeftEnd() );

	GenericMatchSeqManipulator<MatchType> ism( seqI );
	applyBreakpoints( bp_list, iv_list, ism );
}

// T should be a pointer type
template<class T, class Manipulator>
void applyAncestralBreakpoints( const vector< SuperInterval >& siv_list, vector<T>& ord, uint seqI, Manipulator& m )
{
	// make bp list
	vector<gnSeqI> bp_list(siv_list.size()*2, 0);
	size_t cur = 0;
	for( size_t i = 0; i < siv_list.size(); i++ )
	{
		if( siv_list[i].reference_iv.Start(seqI) == NO_MATCH )
			continue;
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI);
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI) + siv_list[i].reference_iv.Length(seqI);
	}
	bp_list.resize(cur);
	// sort the breakpoints and apply...
	sort( bp_list.begin(), bp_list.end() );
	applyBreakpoints( bp_list, ord, m );
}


// assuming breakpoints have been propagated in both directions
// there should now be a 1-to-1 correspondence between superintervals
// in the ancestor and descendants.
void ProgressiveAligner::linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor )
{
	// TODO: speed this up by implementing O(N) instead of O(N^2)
	vector<SuperInterval>& a_ord = alignment_tree[ancestor].ordering;
	vector<SuperInterval>& c_ord = alignment_tree[node1].ordering;
	// initialize all linkages to nothing
	for( size_t aI = 0; aI < a_ord.size(); aI++ )
		if( seqI == 0 )
			a_ord[aI].c1_siv = (std::numeric_limits<size_t>::max)();
		else
			a_ord[aI].c2_siv = (std::numeric_limits<size_t>::max)();
	for( size_t cI = 0; cI < c_ord.size(); cI++ )
		c_ord[cI].parent_siv = (std::numeric_limits<size_t>::max)();

	for( size_t aI = 0; aI < a_ord.size(); aI++ )
	{
		if( a_ord[aI].reference_iv.LeftEnd(seqI) == NO_MATCH )
			continue;
		size_t cI = 0;
		for( ; cI < c_ord.size(); cI++ )
		{
			if( absolut(a_ord[aI].reference_iv.Start(seqI)) != c_ord[cI].LeftEnd() )
				continue;
			if( a_ord[aI].reference_iv.Length(seqI) != c_ord[cI].Length() )
			{
				breakHere();
				cerr << "mapping length mismatch\n";
				cerr << "ancestor: " << ancestor << "\t node1: " << node1 << endl;
				cerr << "a_ord[" << aI << "].reference_iv.Length(" << seqI << "): " << a_ord[aI].reference_iv.Length(seqI) << endl;
				cerr << "a_ord[" << aI << "].reference_iv.LeftEnd(" << seqI << "): " << a_ord[aI].reference_iv.LeftEnd(seqI) << endl;
				cerr << "c_ord[" << cI << "].Length(): " << c_ord[cI].Length() << endl;
				cerr << "c_ord[" << cI << "].LeftEnd(): " << c_ord[cI].LeftEnd() << endl;
				cerr << "";
				cerr << "";
			}
			// link these
			if( seqI == 0 )
				a_ord[aI].c1_siv = cI;
			else
				a_ord[aI].c2_siv = cI;
			c_ord[cI].parent_siv = aI;
			break;
		}
		if( cI == c_ord.size() )
		{
			breakHere();
			cerr << "error no mapping\n";
		}
	}
}


void ProgressiveAligner::translateGappedCoordinates( vector<AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor )
{
	// determine the path that must be traversed
	vector< node_id_t > trans_path;
	getPath( extant, ancestor, trans_path );

	// set seqI to forward orientation 
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();

	// for each node on the path, construct a complete coordinate translation
	for( size_t nI = 1; nI < trans_path.size(); nI++ )
	{
		// first sort matches on start pos and make them all forward oriented
		// then split them on superinterval boundaries and assign each to a superinterval
		// then convert each match's coordinates to be superinterval-local
		// then apply the coordinate translation with transposeCoordinates
		// then shift each match's coordinates to the global ancestral coordinate space
		SSC<AbstractMatch> ssc(seqI);
		sort(ml.begin(), ml.end(), ssc);

		// split on superinterval boundaries
		vector< SuperInterval >& siv_list = alignment_tree[trans_path[nI]].ordering;
		vector< vector< AbstractMatch* > > siv_matches = vector< vector< AbstractMatch* > >(siv_list.size());
		size_t cur_child = 0;
		if( alignment_tree[trans_path[nI]].children[0] == trans_path[nI-1] )
			cur_child = 0;
		else if( alignment_tree[trans_path[nI]].children[1] == trans_path[nI-1] )
			cur_child = 1;
		else 
		{
			breakHere();
			cerr << "forest fire\n";
		}

		AbstractMatchSeqManipulator amsm( seqI );
		applyAncestralBreakpoints(siv_list, ml, cur_child, amsm );
		
		// sort matches again because new ones were added at the end
		sort(ml.begin(), ml.end(), ssc);

		// assign each match to a siv, and convert coords to siv-local
		for( size_t mI = 0; mI < ml.size(); mI++ )
		{
			if( ml[mI]->LeftEnd(seqI) == 0 )
			{
				breakHere();
				cerr << "fefefe";
			}
			size_t sivI = 0;
			for( ; sivI < siv_list.size(); sivI++ )
			{
				if( siv_list[sivI].reference_iv.LeftEnd(cur_child) == NO_MATCH )
					continue;
				if( ml[mI]->LeftEnd(seqI) >= siv_list[sivI].reference_iv.LeftEnd(cur_child) &&
					ml[mI]->LeftEnd(seqI) < siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
					break;
			}
			if( sivI == siv_list.size() )
			{
				cerr << "nI is: "<< nI << endl;
				cerr << "trans_path: ";
				for( size_t ttI = 0; ttI < trans_path.size(); ttI++ )
					cerr << "  " << trans_path[ttI];
				cerr << endl;
				cerr << "problem seq: " << seqI << std::endl;
				cerr << "ml[" << mI << "]->Start(0) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(0) == " << ml[mI]->Length(1) << endl;
				cerr << "ml[" << mI << "]->Start(1) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(1) == " << ml[mI]->Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				for( sivI = 0; sivI < siv_list.size(); sivI++ )
				{
					cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
						cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
					cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				}
				breakHere();
			}
			if( ml[mI]->LeftEnd(seqI) + ml[mI]->Length(seqI) > 
				siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
			{
				cerr << "doesn't fit\n";
				cerr << "ml[" << mI << "]->LeftEnd(" << seqI << "): " << ml[mI]->LeftEnd(seqI) << endl;
				cerr << "ml[" << mI << "]->RightEnd(" << seqI << "): " << ml[mI]->RightEnd(seqI) << endl;
				cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
				if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
					cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
				cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				cerr << "siv_list.size(): " << siv_list.size() << endl;
				cerr << "trans_path:";
				for( size_t tI = 0; tI < trans_path.size(); tI++ )
					cerr << " " << trans_path[tI];
				cerr << endl;
				cerr << "trans_path[" << nI << "]: " << trans_path[nI] << endl;
				breakHere();
			}

			ml[mI]->SetLeftEnd( seqI, ml[mI]->LeftEnd(seqI) - siv_list[sivI].reference_iv.LeftEnd(cur_child) + 1 );
			// if this interval matches the reverse strand then we should effectively invert all matches
			if( siv_list[sivI].reference_iv.Start(cur_child) < 0 )
			{
				int64 new_lend = siv_list[sivI].reference_iv.Length(cur_child) - ml[mI]->LeftEnd(seqI);
				new_lend -= ml[mI]->Length( seqI ) - 2;
				new_lend *= ml[mI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				ml[mI]->Invert();
				ml[mI]->SetStart( seqI, new_lend ); 
			}
			siv_matches[sivI].push_back( ml[mI] );
		}

		// apply the coordinate translation
		ml.clear();
		for( size_t sivI = 0; sivI < siv_matches.size(); sivI++ )
		{
			if( siv_matches[sivI].size() == 0 )
				continue;
			
			// get a CompactGappedAlignment<> for this interval
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			bool alloc_new_siv = false;
			CompactGappedAlignment<> tmp_cga;
			if( siv_cga == NULL )
			{
				alloc_new_siv = true;
				siv_cga = tmp_cga.Copy();
				CompactGappedAlignment<> dorkas(siv_list[sivI].reference_iv);
				*siv_cga = dorkas;
			}

			// now translate each match...
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				CompactGappedAlignment<>* match_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_matches[sivI][mI]);
				bool alloc_new = false;
				if( match_cga == NULL )
				{
					match_cga = tmp_cga.Copy();
					*match_cga = CompactGappedAlignment<>(*(siv_matches[sivI][mI]));
					alloc_new = true;
				}
				siv_cga->translate( *match_cga, seqI, cur_child );

				if( alloc_new )
				{
					siv_matches[sivI][mI]->Free();
					siv_matches[sivI][mI] = match_cga;
				}
			}

			// shift coordinates back to global space
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				int64 cur_start = siv_matches[sivI][mI]->Start(seqI);
				if( cur_start > 0 )
					siv_matches[sivI][mI]->SetStart( seqI, cur_start + siv_list[sivI].LeftEnd() - 1 );
				else
					siv_matches[sivI][mI]->SetStart( seqI, cur_start - siv_list[sivI].LeftEnd() + 1);
				if( (siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) > siv_list.back().LeftEnd() + siv_list.back().Length() )
					 )
				{
					// is there something wrong with the translation table?
					cerr << "siv left is: " << siv_list[sivI].LeftEnd() << endl;
					cerr << "siv right is: " << siv_list[sivI].LeftEnd() + siv_list[sivI].Length() << endl;
					cerr << "match right is: " << siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) << endl;
					cerr << "superseq right is: " << siv_list.back().LeftEnd() + siv_list.back().Length() << endl;
					cerr << "";
					breakHere();
				}
				if( debug_aligner && siv_matches[sivI][mI]->Start(seqI) == 0 )
				{
					breakHere();
				}
			}
			if(alloc_new_siv)
				siv_cga->Free();
			ml.insert( ml.end(), siv_matches[sivI].begin(), siv_matches[sivI].end() );
		}
	}
	// restore forward orientation seqI
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();
}

class SuperIntervalPtrComp
{
public:
	bool operator()( const SuperInterval* a, const SuperInterval* b )
	{
		return (*a) < (*b);
	}
};

void ProgressiveAligner::recursiveApplyAncestralBreakpoints( node_id_t ancestor )
{
	stack<node_id_t> node_stack;
	node_stack.push(ancestor);
	while( node_stack.size() > 0 )
	{
		// pop the current node, apply ancestral breakpoints, recurse on children
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		SuperIntervalManipulator sim;
		if( progress_msgs ) cout << "cur node: " << cur_node << endl;
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			AlignmentTreeNode& atn = alignment_tree[alignment_tree[cur_node].children[childI]];
			if( progress_msgs ) cout << "childI " << childI << " aab\n";
			applyAncestralBreakpoints( alignment_tree[cur_node].ordering, atn.ordering, childI, sim );
			if( progress_msgs ) cout << "sort childI " << childI << "\n";
			vector<SuperInterval*> siv_ptr_list(atn.ordering.size());
			for( size_t sivI = 0; sivI < atn.ordering.size(); ++sivI )
				siv_ptr_list[sivI] = &(atn.ordering[sivI]);
			SuperIntervalPtrComp sipc;
			sort( siv_ptr_list.begin(), siv_ptr_list.end(), sipc );
			vector< SuperInterval > siv_list;
			for( size_t sivI = 0; sivI < siv_ptr_list.size(); ++sivI )
				siv_list.push_back(*siv_ptr_list[sivI]);
			swap(siv_list, atn.ordering);
			node_stack.push( alignment_tree[cur_node].children[childI] );
		}
		if( debug_aligner && alignment_tree[cur_node].children.size() > 0 )
			validateSuperIntervals(alignment_tree[cur_node].children[0], alignment_tree[cur_node].children[1], cur_node);
		if( progress_msgs ) cout << "linking node " << cur_node << "'s" << alignment_tree[cur_node].ordering.size() << " superintervals\n"; 
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
			linkSuperIntervals( alignment_tree[cur_node].children[childI], childI, cur_node );
	}
}


bool getInterveningCoordinates( const AbstractMatch* iv, uint oseqI, Match* r_begin, Match* r_end, uint seqI, int64& gap_lend, int64& gap_rend ){
	// skip this sequence if it's undefined
	if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH) ||
		(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
		gap_lend = 0;
		gap_rend = 0;
		return true;
	}
			
	// determine the size of the gap
	gap_rend = r_end != NULL ? r_end->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
	gap_lend = r_begin != NULL ? r_begin->End( seqI ) + 1 : iv->LeftEnd( oseqI );
	if( gap_rend < 0 || gap_lend < 0 ){
		gap_rend = r_begin != NULL ? -r_begin->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
		gap_lend = r_end != NULL ? -r_end->Start( seqI ) + r_end->Length() : 1;
	}
	if( gap_rend <= 0 || gap_lend <= 0 ){
		// if either is still < 0 then there's a problem...
		genome::ErrorMsg( "Error constructing intervening coordinates" );
	}
	return true;
}


void ProgressiveAligner::pairwiseAnchorSearch( MatchList& r_list, Match* r_begin, Match* r_end, const AbstractMatch* iv, uint oseqI, uint oseqJ )
{
	uint seqI = 0;
	MatchList gap_list;
	vector< int64 > starts;
// 
//	Get the sequence in the intervening gaps between these two matches
//
	for( seqI = 0; seqI < 2; seqI++ )
	{
		int64 gap_end = 0;
		int64 gap_start = 0;
		getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_start, gap_end);
		int64 diff = gap_end - gap_start;
		diff = diff > 0 ? diff - 1 : 0;

		starts.push_back( gap_start );
		gnSequence* new_seq = NULL;
		if(diff > 0 && gap_start + diff - 1 <= r_list.seq_table[ seqI ]->length())
			new_seq = new gnSequence( r_list.seq_table[ seqI ]->ToString( diff, gap_start ) );
		else
			new_seq = new gnSequence();
		gap_list.seq_table.push_back( new_seq );
		gap_list.sml_table.push_back( new DNAMemorySML() );
	}

	gnSeqI avg_len = (gap_list.seq_table[0]->length() + gap_list.seq_table[1]->length())/2;
	uint search_seed_size = getDefaultSeedWeight( avg_len );
	gap_mh.get().Clear();
	
	uint seed_count = use_seed_families ? 3 : 1;
	for( size_t seedI = 0; seedI < seed_count; seedI++ )
	{
		//
		//	Create sorted mer lists for the intervening gap region
		//
		uint64 default_seed = getSeed( search_seed_size, seedI );
		if( search_seed_size < MIN_DNA_SEED_WEIGHT )
		{
			for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
				delete gap_list.seq_table[ seqI ];
			for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
				delete gap_list.sml_table[ seqI ];
			return;
		}
		for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ ){
			gap_list.sml_table[ seqI ]->Clear();
			gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
		}

		//
		//	Find all matches in the gap region
		//
		gap_mh.get().ClearSequences();
		if(seed_count>1)
		{
			MatchList cur_list = gap_list;
			gap_mh.get().FindMatches( cur_list );
			for( size_t mI = 0; mI < cur_list.size(); mI++ )
				cur_list[mI]->Free();
		}else
			gap_mh.get().FindMatches( gap_list );
	}
	if(seed_count>1)
		gap_mh.get().GetMatchList(gap_list);

	EliminateOverlaps_v2( gap_list );

	// for anchor accuracy, throw out any anchors that are shorter than the minimum
	// anchor length after EliminateOverlaps()
	gap_list.LengthFilter( MIN_ANCHOR_LENGTH + 3 );

	for( size_t gI = 0; gI < gap_list.size(); gI++ )
	{
		for( seqI = 0; seqI < 2; seqI++ )
		{
			int64 gap_rend = 0;
			int64 gap_lend = 0;
			getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_lend, gap_rend);
			gap_list[gI]->SetLeftEnd(seqI, gap_list[gI]->LeftEnd(seqI) + gap_lend - 1);
		}
	}
	r_list.insert(r_list.end(), gap_list.begin(), gap_list.end());

	// delete sequences and smls
	for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
		delete gap_list.seq_table[ seqI ];
	for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
		delete gap_list.sml_table[ seqI ];
}

template<class GappedAlignmentType>
void ProgressiveAligner::recurseOnPairs( const vector<node_id_t>& node1_seqs, const vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, Matrix<MatchList>& matches, Matrix< std::vector< search_cache_t > >& search_cache_db, Matrix< std::vector< search_cache_t > >& new_cache_db, boost::multi_array< vector< vector< int64 > >, 2 >& iv_regions )
{
	matches = Matrix<MatchList>(node1_seqs.size(),node2_seqs.size());

	std::vector< bitset_t > aln_matrix;
	iv.GetAlignment(aln_matrix);
	Match tmp(2);
	const size_t sizer = node1_seqs.size() * node2_seqs.size();
	std::vector< std::pair<size_t,size_t> > node_pairs(sizer);
	int nni = 0;
	for( size_t n1 = 0; n1 < node1_seqs.size(); n1++ )
		for( size_t n2 = 0; n2 < node2_seqs.size(); n2++ )
			node_pairs[nni++] = make_pair(n1,n2);

#pragma omp parallel for
	for(int ni = 0; ni < node_pairs.size(); ni++)
	{
		size_t n1 = node_pairs[ni].first;
		size_t n2 = node_pairs[ni].second;
		vector<node_id_t>::const_iterator n1_iter = node1_seqs.begin() + n1;
		vector<node_id_t>::const_iterator n2_iter = node2_seqs.begin() + n2;
		
		uint seqI = node_sequence_map[*n1_iter];
		uint seqJ = node_sequence_map[*n2_iter];
		MatchList& mlist = matches(n1, n2);
		std::vector< search_cache_t >& cache = search_cache_db(n1, n2);
		std::vector< search_cache_t >& new_cache = new_cache_db(n1, n2);
		mlist.seq_table.push_back( alignment_tree[*n1_iter].sequence );
		mlist.seq_table.push_back( alignment_tree[*n2_iter].sequence );

		if( iv.LeftEnd(seqI) == NO_MATCH )
		{
			if( iv.LeftEnd(seqJ) != NO_MATCH )
			{
				iv_regions[n1][n2][1].push_back(iv.LeftEnd(seqJ));
				iv_regions[n1][n2][1].push_back(iv.RightEnd(seqJ));
			}
			continue;	// no sense searching one isn't defined!
		}
		if(iv.LeftEnd(seqJ) == NO_MATCH )
		{
			if( iv.LeftEnd(seqI) != NO_MATCH )
			{
				iv_regions[n1][n2][0].push_back(iv.LeftEnd(seqI));
				iv_regions[n1][n2][0].push_back(iv.RightEnd(seqI));
			}
			continue;	// no sense searching one isn't defined!
		}

		gnSeqI charI = 0;
		gnSeqI charJ = 0;
		const size_t iv_aln_length = iv.AlignmentLength();

// first determine the outer aligned boundaries of the LCB and record them for
// later use
		pair< int64, int64 > pair_1l(0,0);
		pair< int64, int64 > pair_1r(0,0);
		pair< int64, int64 > pair_2l(0,0);
		pair< int64, int64 > pair_2r(0,0);
		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || (aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( colI == 0 )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				else
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				else
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				break;
			}
			if( colI < iv_aln_length && aln_matrix[seqI].test(colI) )
				++charI;
			if( colI < iv_aln_length && aln_matrix[seqJ].test(colI) )
				++charJ;
		}

		charI = 0;
		charJ = 0;
		for( uint colI = iv_aln_length; colI > 0 ; colI-- )
		{
			if( (aln_matrix[seqI].test(colI-1) && aln_matrix[seqJ].test(colI-1)) )
			{
				if( colI == iv_aln_length )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				else
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				else
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				break;
			}
			if( aln_matrix[seqI].test(colI-1) )
				++charI;
			if( aln_matrix[seqJ].test(colI-1) )
				++charJ;
		}
		if( pair_1l.first < pair_1l.second )
		{
			iv_regions[n1][n2][0].push_back(pair_1l.first);
			iv_regions[n1][n2][0].push_back(pair_1l.second);
		}
		if( pair_1r.first < pair_1r.second )
		{
			if( pair_1l.first < pair_1l.second && pair_1r.first == pair_1l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][0].back() = pair_1r.second;
			}else{
				iv_regions[n1][n2][0].push_back(pair_1r.first);
				iv_regions[n1][n2][0].push_back(pair_1r.second);
				if( pair_1r.first <= pair_1l.second && pair_1r.second >= pair_1l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_1l.first << '\t' << pair_1l.second << " right:  " << pair_1r.first << '\t' << pair_1r.second << endl;
					cout << "0 iv.Start(" << seqI << "): " << iv.Start(seqI) << '\t' << "iv.RightEnd(" << seqI << "): " << iv.RightEnd(seqI) << endl;
					if( pair_1l.first == 0 )
						genome::breakHere();
				}
			}
		}

		if( pair_2l.first < pair_2l.second )
		{
			iv_regions[n1][n2][1].push_back(pair_2l.first);
			iv_regions[n1][n2][1].push_back(pair_2l.second);
		}
		if( pair_2r.first < pair_2r.second )
		{
			if( pair_2l.first < pair_2l.second && pair_2r.first == pair_2l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][1].back() = pair_2r.second;
			}else{
				iv_regions[n1][n2][1].push_back(pair_2r.first);
				iv_regions[n1][n2][1].push_back(pair_2r.second);
				if( pair_2r.first <= pair_2l.second && pair_2r.second >= pair_2l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_2l.first << '\t' << pair_2l.second << " right:  " << pair_2r.first << '\t' << pair_2r.second << endl;
					cout << "1 iv.Start(" << seqJ << "): " << iv.Start(seqJ) << '\t' << "iv.RightEnd(" << seqJ << "): " << iv.RightEnd(seqJ) << endl;
					cout << "charI " << charI << "\tcharJ" << charJ << endl;
					if( pair_2l.first == 0 )
						genome::breakHere();
				}
			}
		}

		charI = 0;
		charJ = 0;
		gnSeqI prev_charI = 0;
		gnSeqI prev_charJ = 0;
		bool in_gap = false;

		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || 
				(aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( in_gap && 
					charI - prev_charI > min_recursive_gap_length &&
					charJ - prev_charJ > min_recursive_gap_length )
				{

					Match* l_match = NULL;
					l_match = tmp.Copy();
					if(iv.Orientation(seqI) == AbstractMatch::forward)
						l_match->SetLeftEnd(0, iv.LeftEnd(seqI)+prev_charI);
					else
					{
						l_match->SetLeftEnd(0, iv.RightEnd(seqI)-prev_charI);
						l_match->SetOrientation(0, AbstractMatch::reverse );
				/*******************************************************************************
 * $Id: progressiveAligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * BEWARE!!
 * This code was created in the likeness of the flying spaghetti monster
 *
 * dedicated to Loren...
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libMems/ProgressiveAligner.h"
#include "libMems/GreedyBreakpointElimination.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libMems/gnAlignedSequences.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/PairwiseMatchFinder.h"
#include "libMems/TreeUtilities.h"
#include "libMems/PairwiseMatchAdapter.h"
#include "libMems/DistanceMatrix.h"

#include <boost/dynamic_bitset.hpp>
#include <tuple>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>
#include <iomanip>

#include "stdlib.h"
#include <cstdint>

using namespace std;
using namespace genome;

namespace mems {


bool progress_msgs = false;

bool debug_me = false;
static int dbg_count = 0; 	 


double min_window_size = 200;
double max_window_size = 20000;  // don't feed MUSCLE anything bigger than this
double min_density = .5;
double max_density = .9;
size_t max_gap_length = 5000;
size_t lcb_hangover = 300;


void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list );

/**
 * Test code to ensure that an individual LCB is truly collinear
 * @return	true if the LCB is good
 */
bool my_validateLCB( MatchList& lcb ){
	vector< Match* >::iterator lcb_iter = lcb.begin();
	if( lcb.size() == 0 )
		return true;
	uint seq_count = (*lcb_iter)->SeqCount();
	uint seqI = 0;
	bool complain = false;
	for(; seqI < seq_count; seqI++ ){
		lcb_iter = lcb.begin();
		int64 prev_coord = 0;
		for(; lcb_iter != lcb.end(); ++lcb_iter ){
			if( (*lcb_iter)->Start( seqI ) == NO_MATCH )
				continue;
			else if( prev_coord != 0 && (*lcb_iter)->Start( seqI ) < prev_coord ){
				complain = true;
			}
			prev_coord = (*lcb_iter)->Start( seqI );
		}
	}
	return !complain;
}

template< class BoostMatType >
void print2d_matrix( BoostMatType& mat, std::ostream& os )
{
	for( size_t i = 0; i < mat.shape()[0]; ++i )
	{
		for( size_t j = 0; j < mat.shape()[1]; ++j )
		{
			if( j > 0 )
				os << "\t";
			os << mat[i][j];
		}
		os << endl;
	}
}

double getDefaultBreakpointPenalty( std::vector< gnSequence* >& sequences )
{
	uint default_mer_size = MatchList::GetDefaultMerSize( sequences );
	double avg_seq_len = 0;
	for( size_t seqI = 0; seqI < sequences.size(); ++seqI )
		avg_seq_len += (double)sequences[seqI]->length();
	avg_seq_len /= (double)sequences.size();
	avg_seq_len = log( avg_seq_len ) / log( 2.0 );
	return avg_seq_len * 7000;	  // seems to work reasonably well?
}


double getDefaultBpDistEstimateMinScore( std::vector< gnSequence* >& sequences )
{
	// this value was empirically derived by a process that involved burning incense
	// and uttering arcane words
	return 3.0 * getDefaultBreakpointPenalty(sequences);
}



/*
 * A progressive alignment algorithm for genomes with rearrangements.
 * Start simple, add complexity later.
 * TODO: rewrite the algorithm outline
 */

ProgressiveAligner::ProgressiveAligner( uint seq_count ) :
Aligner( seq_count ),
breakpoint_penalty( -1 ),
min_breakpoint_penalty( 4000 ),
debug(false),
refine(true),
scoring_scheme(ExtantSumOfPairsScoring),
use_weight_scaling(true),
conservation_dist_scale(1),
bp_dist_scale(.9),
max_gapped_alignment_length(20000),
bp_dist_estimate_score(-1),
use_seed_families(false),
using_cache_db(true)
{
	gapped_alignment = true;
	max_window_size = max_gapped_alignment_length;
}

void ProgressiveAligner::SetMaxGappedAlignmentLength( size_t len )
{ 
	max_gapped_alignment_length = len; 
	max_window_size = max_gapped_alignment_length;
}

/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getAlignedChildren( node_id_t node, vector< node_id_t >& descendants )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( node );
	vector< bool > visited( alignment_tree.size(), false );
	descendants.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		if(progress_msgs) cout << "Evaluating aligned nodes linked to node " << cur_node << endl;
		node_stack.pop();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].children_aligned[childI] && !visited[child_id])
				node_stack.push( child_id );
		}
		if( alignment_tree[ cur_node ].sequence != NULL )
			descendants.push_back( cur_node );
	}
}


/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getPath( node_id_t first_n, node_id_t last_n, vector< node_id_t >& path )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( last_n );
	vector< bool > visited( alignment_tree.size(), false );
	while( node_stack.top() != first_n )
	{
		node_id_t cur_node = node_stack.top();
		size_t pre_size = node_stack.size();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if(!visited[child_id])
			{
				node_stack.push( child_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		for( uint parentI = 0; parentI < alignment_tree[cur_node].parents.size(); parentI++ )
		{
			node_id_t parent_id = alignment_tree[cur_node].parents[parentI];
			if(!visited[parent_id])
			{
				node_stack.push( parent_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		node_stack.pop();	// didn't make any progress
	}
	path = vector< node_id_t >( node_stack.size() );
	for( size_t pI = 0; pI < path.size(); pI++ )
	{
		path[pI] = node_stack.top();
		node_stack.pop();
	}
}






template<class MatchType>
void ProgressiveAligner::propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector<MatchType*>& iv_list )
{
	SSC<MatchType> ilc(seqI);
	sort( iv_list.begin(), iv_list.end(), ilc );
	vector< SuperInterval >& ord = alignment_tree[ node1 ].ordering;
	vector<gnSeqI> bp_list;
	for( size_t sI = 0; sI < ord.size(); sI++ )
		bp_list.push_back( ord[sI].LeftEnd() );

	GenericMatchSeqManipulator<MatchType> ism( seqI );
	applyBreakpoints( bp_list, iv_list, ism );
}

// T should be a pointer type
template<class T, class Manipulator>
void applyAncestralBreakpoints( const vector< SuperInterval >& siv_list, vector<T>& ord, uint seqI, Manipulator& m )
{
	// make bp list
	vector<gnSeqI> bp_list(siv_list.size()*2, 0);
	size_t cur = 0;
	for( size_t i = 0; i < siv_list.size(); i++ )
	{
		if( siv_list[i].reference_iv.Start(seqI) == NO_MATCH )
			continue;
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI);
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI) + siv_list[i].reference_iv.Length(seqI);
	}
	bp_list.resize(cur);
	// sort the breakpoints and apply...
	sort( bp_list.begin(), bp_list.end() );
	applyBreakpoints( bp_list, ord, m );
}


// assuming breakpoints have been propagated in both directions
// there should now be a 1-to-1 correspondence between superintervals
// in the ancestor and descendants.
void ProgressiveAligner::linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor )
{
	// TODO: speed this up by implementing O(N) instead of O(N^2)
	vector<SuperInterval>& a_ord = alignment_tree[ancestor].ordering;
	vector<SuperInterval>& c_ord = alignment_tree[node1].ordering;
	// initialize all linkages to nothing
	for( size_t aI = 0; aI < a_ord.size(); aI++ )
		if( seqI == 0 )
			a_ord[aI].c1_siv = (std::numeric_limits<size_t>::max)();
		else
			a_ord[aI].c2_siv = (std::numeric_limits<size_t>::max)();
	for( size_t cI = 0; cI < c_ord.size(); cI++ )
		c_ord[cI].parent_siv = (std::numeric_limits<size_t>::max)();

	for( size_t aI = 0; aI < a_ord.size(); aI++ )
	{
		if( a_ord[aI].reference_iv.LeftEnd(seqI) == NO_MATCH )
			continue;
		size_t cI = 0;
		for( ; cI < c_ord.size(); cI++ )
		{
			if( absolut(a_ord[aI].reference_iv.Start(seqI)) != c_ord[cI].LeftEnd() )
				continue;
			if( a_ord[aI].reference_iv.Length(seqI) != c_ord[cI].Length() )
			{
				breakHere();
				cerr << "mapping length mismatch\n";
				cerr << "ancestor: " << ancestor << "\t node1: " << node1 << endl;
				cerr << "a_ord[" << aI << "].reference_iv.Length(" << seqI << "): " << a_ord[aI].reference_iv.Length(seqI) << endl;
				cerr << "a_ord[" << aI << "].reference_iv.LeftEnd(" << seqI << "): " << a_ord[aI].reference_iv.LeftEnd(seqI) << endl;
				cerr << "c_ord[" << cI << "].Length(): " << c_ord[cI].Length() << endl;
				cerr << "c_ord[" << cI << "].LeftEnd(): " << c_ord[cI].LeftEnd() << endl;
				cerr << "";
				cerr << "";
			}
			// link these
			if( seqI == 0 )
				a_ord[aI].c1_siv = cI;
			else
				a_ord[aI].c2_siv = cI;
			c_ord[cI].parent_siv = aI;
			break;
		}
		if( cI == c_ord.size() )
		{
			breakHere();
			cerr << "error no mapping\n";
		}
	}
}


void ProgressiveAligner::translateGappedCoordinates( vector<AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor )
{
	// determine the path that must be traversed
	vector< node_id_t > trans_path;
	getPath( extant, ancestor, trans_path );

	// set seqI to forward orientation 
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();

	// for each node on the path, construct a complete coordinate translation
	for( size_t nI = 1; nI < trans_path.size(); nI++ )
	{
		// first sort matches on start pos and make them all forward oriented
		// then split them on superinterval boundaries and assign each to a superinterval
		// then convert each match's coordinates to be superinterval-local
		// then apply the coordinate translation with transposeCoordinates
		// then shift each match's coordinates to the global ancestral coordinate space
		SSC<AbstractMatch> ssc(seqI);
		sort(ml.begin(), ml.end(), ssc);

		// split on superinterval boundaries
		vector< SuperInterval >& siv_list = alignment_tree[trans_path[nI]].ordering;
		vector< vector< AbstractMatch* > > siv_matches = vector< vector< AbstractMatch* > >(siv_list.size());
		size_t cur_child = 0;
		if( alignment_tree[trans_path[nI]].children[0] == trans_path[nI-1] )
			cur_child = 0;
		else if( alignment_tree[trans_path[nI]].children[1] == trans_path[nI-1] )
			cur_child = 1;
		else 
		{
			breakHere();
			cerr << "forest fire\n";
		}

		AbstractMatchSeqManipulator amsm( seqI );
		applyAncestralBreakpoints(siv_list, ml, cur_child, amsm );
		
		// sort matches again because new ones were added at the end
		sort(ml.begin(), ml.end(), ssc);

		// assign each match to a siv, and convert coords to siv-local
		for( size_t mI = 0; mI < ml.size(); mI++ )
		{
			if( ml[mI]->LeftEnd(seqI) == 0 )
			{
				breakHere();
				cerr << "fefefe";
			}
			size_t sivI = 0;
			for( ; sivI < siv_list.size(); sivI++ )
			{
				if( siv_list[sivI].reference_iv.LeftEnd(cur_child) == NO_MATCH )
					continue;
				if( ml[mI]->LeftEnd(seqI) >= siv_list[sivI].reference_iv.LeftEnd(cur_child) &&
					ml[mI]->LeftEnd(seqI) < siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
					break;
			}
			if( sivI == siv_list.size() )
			{
				cerr << "nI is: "<< nI << endl;
				cerr << "trans_path: ";
				for( size_t ttI = 0; ttI < trans_path.size(); ttI++ )
					cerr << "  " << trans_path[ttI];
				cerr << endl;
				cerr << "problem seq: " << seqI << std::endl;
				cerr << "ml[" << mI << "]->Start(0) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(0) == " << ml[mI]->Length(1) << endl;
				cerr << "ml[" << mI << "]->Start(1) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(1) == " << ml[mI]->Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				for( sivI = 0; sivI < siv_list.size(); sivI++ )
				{
					cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
						cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
					cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				}
				breakHere();
			}
			if( ml[mI]->LeftEnd(seqI) + ml[mI]->Length(seqI) > 
				siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
			{
				cerr << "doesn't fit\n";
				cerr << "ml[" << mI << "]->LeftEnd(" << seqI << "): " << ml[mI]->LeftEnd(seqI) << endl;
				cerr << "ml[" << mI << "]->RightEnd(" << seqI << "): " << ml[mI]->RightEnd(seqI) << endl;
				cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
				if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
					cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
				cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				cerr << "siv_list.size(): " << siv_list.size() << endl;
				cerr << "trans_path:";
				for( size_t tI = 0; tI < trans_path.size(); tI++ )
					cerr << " " << trans_path[tI];
				cerr << endl;
				cerr << "trans_path[" << nI << "]: " << trans_path[nI] << endl;
				breakHere();
			}

			ml[mI]->SetLeftEnd( seqI, ml[mI]->LeftEnd(seqI) - siv_list[sivI].reference_iv.LeftEnd(cur_child) + 1 );
			// if this interval matches the reverse strand then we should effectively invert all matches
			if( siv_list[sivI].reference_iv.Start(cur_child) < 0 )
			{
				int64 new_lend = siv_list[sivI].reference_iv.Length(cur_child) - ml[mI]->LeftEnd(seqI);
				new_lend -= ml[mI]->Length( seqI ) - 2;
				new_lend *= ml[mI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				ml[mI]->Invert();
				ml[mI]->SetStart( seqI, new_lend ); 
			}
			siv_matches[sivI].push_back( ml[mI] );
		}

		// apply the coordinate translation
		ml.clear();
		for( size_t sivI = 0; sivI < siv_matches.size(); sivI++ )
		{
			if( siv_matches[sivI].size() == 0 )
				continue;
			
			// get a CompactGappedAlignment<> for this interval
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			bool alloc_new_siv = false;
			CompactGappedAlignment<> tmp_cga;
			if( siv_cga == NULL )
			{
				alloc_new_siv = true;
				siv_cga = tmp_cga.Copy();
				CompactGappedAlignment<> dorkas(siv_list[sivI].reference_iv);
				*siv_cga = dorkas;
			}

			// now translate each match...
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				CompactGappedAlignment<>* match_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_matches[sivI][mI]);
				bool alloc_new = false;
				if( match_cga == NULL )
				{
					match_cga = tmp_cga.Copy();
					*match_cga = CompactGappedAlignment<>(*(siv_matches[sivI][mI]));
					alloc_new = true;
				}
				siv_cga->translate( *match_cga, seqI, cur_child );

				if( alloc_new )
				{
					siv_matches[sivI][mI]->Free();
					siv_matches[sivI][mI] = match_cga;
				}
			}

			// shift coordinates back to global space
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				int64 cur_start = siv_matches[sivI][mI]->Start(seqI);
				if( cur_start > 0 )
					siv_matches[sivI][mI]->SetStart( seqI, cur_start + siv_list[sivI].LeftEnd() - 1 );
				else
					siv_matches[sivI][mI]->SetStart( seqI, cur_start - siv_list[sivI].LeftEnd() + 1);
				if( (siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) > siv_list.back().LeftEnd() + siv_list.back().Length() )
					 )
				{
					// is there something wrong with the translation table?
					cerr << "siv left is: " << siv_list[sivI].LeftEnd() << endl;
					cerr << "siv right is: " << siv_list[sivI].LeftEnd() + siv_list[sivI].Length() << endl;
					cerr << "match right is: " << siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) << endl;
					cerr << "superseq right is: " << siv_list.back().LeftEnd() + siv_list.back().Length() << endl;
					cerr << "";
					breakHere();
				}
				if( debug_aligner && siv_matches[sivI][mI]->Start(seqI) == 0 )
				{
					breakHere();
				}
			}
			if(alloc_new_siv)
				siv_cga->Free();
			ml.insert( ml.end(), siv_matches[sivI].begin(), siv_matches[sivI].end() );
		}
	}
	// restore forward orientation seqI
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();
}

class SuperIntervalPtrComp
{
public:
	bool operator()( const SuperInterval* a, const SuperInterval* b )
	{
		return (*a) < (*b);
	}
};

void ProgressiveAligner::recursiveApplyAncestralBreakpoints( node_id_t ancestor )
{
	stack<node_id_t> node_stack;
	node_stack.push(ancestor);
	while( node_stack.size() > 0 )
	{
		// pop the current node, apply ancestral breakpoints, recurse on children
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		SuperIntervalManipulator sim;
		if( progress_msgs ) cout << "cur node: " << cur_node << endl;
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			AlignmentTreeNode& atn = alignment_tree[alignment_tree[cur_node].children[childI]];
			if( progress_msgs ) cout << "childI " << childI << " aab\n";
			applyAncestralBreakpoints( alignment_tree[cur_node].ordering, atn.ordering, childI, sim );
			if( progress_msgs ) cout << "sort childI " << childI << "\n";
			vector<SuperInterval*> siv_ptr_list(atn.ordering.size());
			for( size_t sivI = 0; sivI < atn.ordering.size(); ++sivI )
				siv_ptr_list[sivI] = &(atn.ordering[sivI]);
			SuperIntervalPtrComp sipc;
			sort( siv_ptr_list.begin(), siv_ptr_list.end(), sipc );
			vector< SuperInterval > siv_list;
			for( size_t sivI = 0; sivI < siv_ptr_list.size(); ++sivI )
				siv_list.push_back(*siv_ptr_list[sivI]);
			swap(siv_list, atn.ordering);
			node_stack.push( alignment_tree[cur_node].children[childI] );
		}
		if( debug_aligner && alignment_tree[cur_node].children.size() > 0 )
			validateSuperIntervals(alignment_tree[cur_node].children[0], alignment_tree[cur_node].children[1], cur_node);
		if( progress_msgs ) cout << "linking node " << cur_node << "'s" << alignment_tree[cur_node].ordering.size() << " superintervals\n"; 
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
			linkSuperIntervals( alignment_tree[cur_node].children[childI], childI, cur_node );
	}
}


bool getInterveningCoordinates( const AbstractMatch* iv, uint oseqI, Match* r_begin, Match* r_end, uint seqI, int64& gap_lend, int64& gap_rend ){
	// skip this sequence if it's undefined
	if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH) ||
		(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
		gap_lend = 0;
		gap_rend = 0;
		return true;
	}
			
	// determine the size of the gap
	gap_rend = r_end != NULL ? r_end->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
	gap_lend = r_begin != NULL ? r_begin->End( seqI ) + 1 : iv->LeftEnd( oseqI );
	if( gap_rend < 0 || gap_lend < 0 ){
		gap_rend = r_begin != NULL ? -r_begin->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
		gap_lend = r_end != NULL ? -r_end->Start( seqI ) + r_end->Length() : 1;
	}
	if( gap_rend <= 0 || gap_lend <= 0 ){
		// if either is still < 0 then there's a problem...
		genome::ErrorMsg( "Error constructing intervening coordinates" );
	}
	return true;
}


void ProgressiveAligner::pairwiseAnchorSearch( MatchList& r_list, Match* r_begin, Match* r_end, const AbstractMatch* iv, uint oseqI, uint oseqJ )
{
	uint seqI = 0;
	MatchList gap_list;
	vector< int64 > starts;
// 
//	Get the sequence in the intervening gaps between these two matches
//
	for( seqI = 0; seqI < 2; seqI++ )
	{
		int64 gap_end = 0;
		int64 gap_start = 0;
		getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_start, gap_end);
		int64 diff = gap_end - gap_start;
		diff = diff > 0 ? diff - 1 : 0;

		starts.push_back( gap_start );
		gnSequence* new_seq = NULL;
		if(diff > 0 && gap_start + diff - 1 <= r_list.seq_table[ seqI ]->length())
			new_seq = new gnSequence( r_list.seq_table[ seqI ]->ToString( diff, gap_start ) );
		else
			new_seq = new gnSequence();
		gap_list.seq_table.push_back( new_seq );
		gap_list.sml_table.push_back( new DNAMemorySML() );
	}

	gnSeqI avg_len = (gap_list.seq_table[0]->length() + gap_list.seq_table[1]->length())/2;
	uint search_seed_size = getDefaultSeedWeight( avg_len );
	gap_mh.get().Clear();
	
	uint seed_count = use_seed_families ? 3 : 1;
	for( size_t seedI = 0; seedI < seed_count; seedI++ )
	{
		//
		//	Create sorted mer lists for the intervening gap region
		//
		uint64 default_seed = getSeed( search_seed_size, seedI );
		if( search_seed_size < MIN_DNA_SEED_WEIGHT )
		{
			for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
				delete gap_list.seq_table[ seqI ];
			for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
				delete gap_list.sml_table[ seqI ];
			return;
		}
		for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ ){
			gap_list.sml_table[ seqI ]->Clear();
			gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
		}

		//
		//	Find all matches in the gap region
		//
		gap_mh.get().ClearSequences();
		if(seed_count>1)
		{
			MatchList cur_list = gap_list;
			gap_mh.get().FindMatches( cur_list );
			for( size_t mI = 0; mI < cur_list.size(); mI++ )
				cur_list[mI]->Free();
		}else
			gap_mh.get().FindMatches( gap_list );
	}
	if(seed_count>1)
		gap_mh.get().GetMatchList(gap_list);

	EliminateOverlaps_v2( gap_list );

	// for anchor accuracy, throw out any anchors that are shorter than the minimum
	// anchor length after EliminateOverlaps()
	gap_list.LengthFilter( MIN_ANCHOR_LENGTH + 3 );

	for( size_t gI = 0; gI < gap_list.size(); gI++ )
	{
		for( seqI = 0; seqI < 2; seqI++ )
		{
			int64 gap_rend = 0;
			int64 gap_lend = 0;
			getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_lend, gap_rend);
			gap_list[gI]->SetLeftEnd(seqI, gap_list[gI]->LeftEnd(seqI) + gap_lend - 1);
		}
	}
	r_list.insert(r_list.end(), gap_list.begin(), gap_list.end());

	// delete sequences and smls
	for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
		delete gap_list.seq_table[ seqI ];
	for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
		delete gap_list.sml_table[ seqI ];
}

template<class GappedAlignmentType>
void ProgressiveAligner::recurseOnPairs( const vector<node_id_t>& node1_seqs, const vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, Matrix<MatchList>& matches, Matrix< std::vector< search_cache_t > >& search_cache_db, Matrix< std::vector< search_cache_t > >& new_cache_db, boost::multi_array< vector< vector< int64 > >, 2 >& iv_regions )
{
	matches = Matrix<MatchList>(node1_seqs.size(),node2_seqs.size());

	std::vector< bitset_t > aln_matrix;
	iv.GetAlignment(aln_matrix);
	Match tmp(2);
	const size_t sizer = node1_seqs.size() * node2_seqs.size();
	std::vector< std::pair<size_t,size_t> > node_pairs(sizer);
	int nni = 0;
	for( size_t n1 = 0; n1 < node1_seqs.size(); n1++ )
		for( size_t n2 = 0; n2 < node2_seqs.size(); n2++ )
			node_pairs[nni++] = make_pair(n1,n2);

#pragma omp parallel for
	for(int ni = 0; ni < node_pairs.size(); ni++)
	{
		size_t n1 = node_pairs[ni].first;
		size_t n2 = node_pairs[ni].second;
		vector<node_id_t>::const_iterator n1_iter = node1_seqs.begin() + n1;
		vector<node_id_t>::const_iterator n2_iter = node2_seqs.begin() + n2;
		
		uint seqI = node_sequence_map[*n1_iter];
		uint seqJ = node_sequence_map[*n2_iter];
		MatchList& mlist = matches(n1, n2);
		std::vector< search_cache_t >& cache = search_cache_db(n1, n2);
		std::vector< search_cache_t >& new_cache = new_cache_db(n1, n2);
		mlist.seq_table.push_back( alignment_tree[*n1_iter].sequence );
		mlist.seq_table.push_back( alignment_tree[*n2_iter].sequence );

		if( iv.LeftEnd(seqI) == NO_MATCH )
		{
			if( iv.LeftEnd(seqJ) != NO_MATCH )
			{
				iv_regions[n1][n2][1].push_back(iv.LeftEnd(seqJ));
				iv_regions[n1][n2][1].push_back(iv.RightEnd(seqJ));
			}
			continue;	// no sense searching one isn't defined!
		}
		if(iv.LeftEnd(seqJ) == NO_MATCH )
		{
			if( iv.LeftEnd(seqI) != NO_MATCH )
			{
				iv_regions[n1][n2][0].push_back(iv.LeftEnd(seqI));
				iv_regions[n1][n2][0].push_back(iv.RightEnd(seqI));
			}
			continue;	// no sense searching one isn't defined!
		}

		gnSeqI charI = 0;
		gnSeqI charJ = 0;
		const size_t iv_aln_length = iv.AlignmentLength();

// first determine the outer aligned boundaries of the LCB and record them for
// later use
		pair< int64, int64 > pair_1l(0,0);
		pair< int64, int64 > pair_1r(0,0);
		pair< int64, int64 > pair_2l(0,0);
		pair< int64, int64 > pair_2r(0,0);
		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || (aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( colI == 0 )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				else
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				else
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				break;
			}
			if( colI < iv_aln_length && aln_matrix[seqI].test(colI) )
				++charI;
			if( colI < iv_aln_length && aln_matrix[seqJ].test(colI) )
				++charJ;
		}

		charI = 0;
		charJ = 0;
		for( uint colI = iv_aln_length; colI > 0 ; colI-- )
		{
			if( (aln_matrix[seqI].test(colI-1) && aln_matrix[seqJ].test(colI-1)) )
			{
				if( colI == iv_aln_length )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				else
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				else
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				break;
			}
			if( aln_matrix[seqI].test(colI-1) )
				++charI;
			if( aln_matrix[seqJ].test(colI-1) )
				++charJ;
		}
		if( pair_1l.first < pair_1l.second )
		{
			iv_regions[n1][n2][0].push_back(pair_1l.first);
			iv_regions[n1][n2][0].push_back(pair_1l.second);
		}
		if( pair_1r.first < pair_1r.second )
		{
			if( pair_1l.first < pair_1l.second && pair_1r.first == pair_1l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][0].back() = pair_1r.second;
			}else{
				iv_regions[n1][n2][0].push_back(pair_1r.first);
				iv_regions[n1][n2][0].push_back(pair_1r.second);
				if( pair_1r.first <= pair_1l.second && pair_1r.second >= pair_1l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_1l.first << '\t' << pair_1l.second << " right:  " << pair_1r.first << '\t' << pair_1r.second << endl;
					cout << "0 iv.Start(" << seqI << "): " << iv.Start(seqI) << '\t' << "iv.RightEnd(" << seqI << "): " << iv.RightEnd(seqI) << endl;
					if( pair_1l.first == 0 )
						genome::breakHere();
				}
			}
		}

		if( pair_2l.first < pair_2l.second )
		{
			iv_regions[n1][n2][1].push_back(pair_2l.first);
			iv_regions[n1][n2][1].push_back(pair_2l.second);
		}
		if( pair_2r.first < pair_2r.second )
		{
			if( pair_2l.first < pair_2l.second && pair_2r.first == pair_2l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][1].back() = pair_2r.second;
			}else{
				iv_regions[n1][n2][1].push_back(pair_2r.first);
				iv_regions[n1][n2][1].push_back(pair_2r.second);
				if( pair_2r.first <= pair_2l.second && pair_2r.second >= pair_2l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_2l.first << '\t' << pair_2l.second << " right:  " << pair_2r.first << '\t' << pair_2r.second << endl;
					cout << "1 iv.Start(" << seqJ << "): " << iv.Start(seqJ) << '\t' << "iv.RightEnd(" << seqJ << "): " << iv.RightEnd(seqJ) << endl;
					cout << "charI " << charI << "\tcharJ" << charJ << endl;
					if( pair_2l.first == 0 )
						genome::breakHere();
				}
			}
		}

		charI = 0;
		charJ = 0;
		gnSeqI prev_charI = 0;
		gnSeqI prev_charJ = 0;
		bool in_gap = false;

		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || 
				(aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( in_gap && 
					charI - prev_charI > min_recursive_gap_length &&
					charJ - prev_charJ > min_recursive_gap_length )
				{

					Match* l_match = NULL;
					l_match = tmp.Copy();
					if(iv.Orientation(seqI) == AbstractMatch::forward)
						l_match->SetLeftEnd(0, iv.LeftEnd(seqI)+prev_charI);
					else
					{
						l_match->SetLeftEnd(0, iv.RightEnd(seqI)-prev_charI);
						l_match->SetOrientation(0, AbstractMatch::reverse );
					}
					if(iv.Orientation(seqJ) == AbstractMatch::forward)
						l_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+prev_charJ);
					else
					{
						l_match->SetLeftEnd(1, iv.RightEnd(seqJ)-prev_charJ);
						l_match->SetOrientation(1, AbstractMatch::reverse );
					}
					l_match->SetLength(0);
					Match* r_match = NULL;
					if( charJ != iv.RightEnd(seqJ) && charI != iv.RightEnd(seqI) )
					{
						r_match = tmp.Copy();
						if(iv.Orientation(seqI) == AbstractMatch::forward)
							r_match->SetLeftEnd(0, iv.LeftEnd(seqI)+charI);
						else
						{
							r_match->SetLeftEnd(0, iv.RightEnd(seqI)-charI);
							r_match->SetOrientation(0, AbstractMatch::reverse );
						}
						if(iv.Orientation(seqJ) == AbstractMatch::forward)
							r_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+charJ);
						else
						{
							r_match->SetLeftEnd(1, iv.RightEnd(seqJ)-charJ);
							r_match->SetOrientation(1, AbstractMatch::reverse );
						}
						r_match->SetLength(0);
					}


					if( iv.Orientation(seqI) == AbstractMatch::reverse )
					{
						swap(l_match,r_match);
						if( l_match != NULL ) l_match->Invert();
						if( r_match != NULL ) r_match->Invert();
					}
					// check whether the current cache already has the searched region
					search_cache_t cacheval = make_pair( l_match, r_match );
					std::vector< search_cache_t >::iterator cache_entry = std::upper_bound( cache.begin(), cache.end(), cacheval, mems::cache_comparator );
					if( cache_entry == cache.end() || 
						(mems::cache_comparator( cacheval, *cache_entry ) || mems::cache_comparator( *cache_entry, cacheval )) )
					{
						// search this region
							pairwiseAnchorSearch(mlist, l_match, r_match, &iv, seqI, seqJ);
					}
					if(using_cache_db)
						new_cache.push_back( cacheval );
				}
				prev_charI = charI;
				prev_charJ = charJ;
				in_gap = false;
			}
			else
				in_gap = true;
			if( colI < iv.AlignmentLength() )
			{
				if( aln_matrix[seqI].test(colI) )
					++charI;
				if( aln_matrix[seqJ].test(colI) )
					++charJ;
			}
		}
	}
}
