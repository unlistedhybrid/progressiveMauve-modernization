/*******************************************************************************
 * $Id: UniqueMatchFinder.cpp,v 1.13 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "UniqueMatchFinder.h"
#include <list>
#include <algorithm>

using namespace std;
using namespace genome;
using namespace mems;

// Explicitly define constructors/destructors in the mems namespace

mems::UniqueMatchFinder::UniqueMatchFinder(){
}

mems::UniqueMatchFinder::~UniqueMatchFinder(){
}

mems::UniqueMatchFinder::UniqueMatchFinder(const UniqueMatchFinder& mh) : MemHash(mh){

}

mems::UniqueMatchFinder* mems::UniqueMatchFinder::Clone() const{
	return new UniqueMatchFinder(*this);
}


// enumerate out every match with unique seeds
bool mems::UniqueMatchFinder::EnumerateMatches( IdmerList& match_list ){

	std::cerr << "DEBUG: UniqueMatchFinder::EnumerateMatches called with " << match_list.size() << " items\n";
	
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
	
	// hash all unique seeds
	bool success = true;
	if(unique_list.size() >= 2){
		std::cerr << "DEBUG: Calling HashMatch with " << unique_list.size() << " unique seeds\n";
		success = HashMatch(unique_list);
	}
	
	std::cerr << "DEBUG: UniqueMatchFinder::EnumerateMatches returning " << (success ? "true" : "false") << "\n";
	return success;
}
