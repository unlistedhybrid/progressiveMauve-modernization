/*******************************************************************************
 * $Id: PairwiseMatchFinder.cpp,v 1.13 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/PairwiseMatchFinder.h"
#include <list>

using namespace std;
using namespace genome;

namespace mems {

PairwiseMatchFinder::PairwiseMatchFinder(){
}

PairwiseMatchFinder::~PairwiseMatchFinder(){
}

PairwiseMatchFinder::PairwiseMatchFinder(const PairwiseMatchFinder& mh) : MemHash(mh){

}

PairwiseMatchFinder* PairwiseMatchFinder::Clone() const{
	return new PairwiseMatchFinder(*this);
}


// enumerate out every pairwise match
bool PairwiseMatchFinder::EnumerateMatches( IdmerList& match_list ){

	std::cerr << "DEBUG: PairwiseMatchFinder::EnumerateMatches called with " << match_list.size() << " items\n";
	
	if(match_list.empty()){
		std::cerr << "DEBUG: match_list is empty, returning true\n";
		return true;
	}

	match_list.sort(&idmer_id_lessthan);
	std::cerr << "DEBUG: List sorted\n";
	
	IdmerList::iterator iter = match_list.begin();
	IdmerList::iterator iter_end = match_list.end();
	unsigned int cur_id_count = 1;
	IdmerList unique_list;
	
	// identify all of the unique seeds and add them to unique_list
	while(iter != iter_end){
		IdmerList::iterator iter_next = iter;
		++iter_next;
		
		if(iter_next == iter_end || iter->id != iter_next->id){
			if(cur_id_count == 1){
				std::cerr << "DEBUG: Adding unique seed with id " << iter->id << " to unique_list\n";
				unique_list.push_back(*iter);
			}
			cur_id_count = 1;
		}else{
			cur_id_count++;
		}
		++iter;
	}
	
	std::cerr << "DEBUG: Found " << unique_list.size() << " unique seeds\n";
	
	// hash each pair of unique seeds
	bool success = true;
	IdmerList::iterator iter_a = unique_list.begin();
	IdmerList::iterator iter_a_end = unique_list.end();
	
	for(; iter_a != iter_a_end; ++iter_a){
		IdmerList::iterator iter_b = iter_a;
		++iter_b;
		
		for(; iter_b != iter_a_end; ++iter_b){
			std::cerr << "DEBUG: Hashing pair of seeds\n";
			IdmerList hash_list;
			hash_list.push_back(*iter_a);
			hash_list.push_back(*iter_b);
			success = success && HashMatch(hash_list);
		}
	}
	
	std::cerr << "DEBUG: PairwiseMatchFinder::EnumerateMatches returning " << (success ? "true" : "false") << "\n";
	return success;
}

}  // namespace mems
