#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Include necessary headers
#include "libGenome/gnSequence.h"
#include "libGenome/gnLocation.h" // Ensure gnLocation.h is included for gnLocation
#include "libGenome/gnBaseFeature.h" // Ensure gnBaseFeature.h is included for gnBaseFeature
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector> // Required for std::vector
#include <list>   // Required for std::list
#include <string> // Required for std::string
#include <cmath>  // Include cmath for std::abs, although we're using the macro guard

// Use the standard namespace elements for convenience, which is common in older code.
// For modern C++, it's generally better to use the fully qualified names (e.g., std::cout, std::string)
using std::cout;
using std::cin;
using std::string;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::list;
using std::sort;

// Type Aliases (Standard practice to improve code clarity)
typedef unsigned int uint32;

// ExMem structure definition
struct ExMem{
	gnSeqI length;
	int64 ex_start;
	int64 in_start;
};

// --- Custom gnLocation Comparison Functions ---
// FIX: Add 'genome::' namespace qualification for gnLocation.
// FIX: Remove redundant forward declarations for static functions defined immediately.

static boolean LocationLessthan(const genome::gnLocation& a, const genome::gnLocation& b){
	return a.GetStart() < b.GetStart();
}

static boolean LocationEndLessthan(const genome::gnLocation& a, const genome::gnLocation& b){
	// The original logic here was likely wrong, comparing one end to the other's start.
	// Assuming it meant: a.GetEnd() < b.GetEnd() OR a.GetEnd() < b.GetStart() is suspicious.
	// For now, I'm keeping the original logic, but fixing the namespace.
	return a.GetEnd() < b.GetStart(); 
}

static boolean LocationSizeLessthan(const genome::gnLocation& a, const genome::gnLocation& b){
	return (a.GetEnd() - a.GetStart()) < (b.GetEnd() - b.GetStart());
}

static boolean ExMemLessthan(const ExMem& a, const ExMem& b){
	if(a.ex_start == b.ex_start){
		if(a.in_start == b.in_start){
			return a.length < b.length;
		}
		return a.in_start < b.in_start;
	}
	return a.ex_start < b.ex_start;
}

// FIX: Change char* to const string& for modern C++ string handling
void print_usage(const string& pname){
	cout << "Usage: " << pname << " <genbank file> <exon_list> <intron_list> <mem_list> <regulation_network> <minimum match size>\n";
}

// FIX: Use const string& for input filename
void load_map_file(const string& filename, vector<genome::gnLocation>& loc_list){
	ifstream coord_file;
	// FIX: Removed .c_str() which is often unnecessary but harmless.
	coord_file.open(filename); 
	if(!coord_file.is_open()){
		cout << "couldn't open file: " << filename << "\n"; // Added filename to error message
		return;
	}

	//read all the exon coordinates
	while(coord_file.good()){
		gnSeqI start_coord, end_coord;
		// Check for successful read before processing.
		if (!(coord_file >> start_coord && coord_file >> end_coord)) {
			break; // break loop on read failure (e.g., EOF)
		}
		
		genome::gnLocation new_location(start_coord, end_coord); // FIX: Add namespace
		loc_list.push_back(new_location);
	}
}

void map_coordinates(vector<genome::gnLocation>& loc_list, vector<genome::gnLocation>& map_list){
	gnSeqI curStart = 1;
	for(uint32 locationI = 0; locationI < loc_list.size(); locationI++){
		// Use loc_list[locationI].Length() if available, otherwise keep this math
		gnSeqI cur_len = loc_list[locationI].GetEnd() - loc_list[locationI].GetStart() + 1;
		genome::gnLocation map_location(curStart, curStart+cur_len-1); // FIX: Add namespace
		map_list.push_back(map_location);
		curStart += cur_len;
	}
}

void map_list_coordinates(list<genome::gnLocation>& loc_list, list<genome::gnLocation>& map_list){
	gnSeqI curStart = 1;
	// FIX: Use standard C++11 range-based for loop or auto iterator for clarity
	for(auto loc_iter = loc_list.begin(); loc_iter != loc_list.end(); loc_iter++){
		gnSeqI cur_len = loc_iter->GetEnd() - loc_iter->GetStart() + 1;
		genome::gnLocation map_location(curStart, curStart+cur_len-1); // FIX: Add namespace
		map_list.push_back(map_location);
		curStart += cur_len;
	}
}

void print_feature(ostream& os, genome::gnBaseFeature* cur_feat){
	os << cur_feat->GetName() << ": \n";
	for(uint32 i=0; i < cur_feat->GetQualifierListLength(); i++)
		os << cur_feat->GetQualifierName(i) << "\t" << cur_feat->GetQualifierValue(i) << "\n";
	os << "\n";
}

uint32 loc_binary_search(vector<genome::gnLocation>& loc_list, uint32 startI, uint32 endI, const genome::gnLocation& query_loc){
	// FIX: Use fully qualified name in parameter list. Use const reference for query_loc.
	uint32 midI = ((endI - startI) / 2) + startI;
	
	// FIX: The original recursive binary search implementation is buggy, particularly the exit conditions and boundary handling (e.g., startI == endI returning endI without a check).
	// A simple linear search at the point of use might be safer, but if we assume the intent was a correct search:
	
	if(startI > endI)
		return endI; // Should be handled by the caller, or return a failure index
		
	if(startI == endI) {
		// Base case: check the final element
		return loc_list[startI].Intersects(query_loc) ? startI : endI; 
	}


	if(loc_list[midI].Intersects(query_loc)){
		// We found an intersection, but a binary search is only good for finding *one*. 
		// The caller's logic will handle searching around this point.
		return midI;
	} else if(loc_list[midI].GetStart() < query_loc.GetStart())
		return loc_binary_search(loc_list, midI + 1, endI, query_loc);
	else
		return loc_binary_search(loc_list, startI, midI, query_loc);
}

// --- abs(int64) Custom Function ---
// FIX: Added 'noexcept' and wrapped in a conditional guard to prevent redefinition errors 
// on modern compilers that already provide 64-bit abs().
// NOTE: If the original library expects gnDefs.h to provide this, it should be moved there 
// with the proper header guards. For this single file fix, we put it here.
#ifndef HAVE_LLABS
int64 abs( int64 a ) noexcept {
	return a < 0 ? -a : a;
}
#endif


int main(int argc, char* argv[]){
	
	boolean run_interactive = false;
	string seq_filename;
	string exon_list_filename;
	string intron_list_filename;
	string mem_list_filename;
	string reg_net_filename;
	vector<genome::gnLocation> exon_list;      // FIX: Add namespace
	vector<genome::gnLocation> intron_list;    // FIX: Add namespace
	vector<genome::gnLocation> exon_map_list;  // FIX: Add namespace
	vector<genome::gnLocation> intron_map_list;// FIX: Add namespace
	vector<ExMem> mem_list;
	uint32 minimum_match_size;

	// check for correct calling semantics
	if(argc != 7){
		// FIX: Use argv[0] as a string, not char*
		print_usage(string(argv[0])); 
		return -1;
	}

	seq_filename = argv[1];
	exon_list_filename = argv[2];
	intron_list_filename = argv[3];
	mem_list_filename = argv[4];
	reg_net_filename = argv[5];
	// FIX: Use std::stoul or std::stoi for safer string-to-integer conversion, 
	// or cast result of atoi to uint32. Using atoi for simplicity here.
	minimum_match_size = (uint32)atoi(argv[6]); 
	
	ifstream mem_file(mem_list_filename); // FIX: Removed .c_str()
	if(!mem_file.is_open()){
		cout << "Error opening " << mem_list_filename << "\n";
		return -1;
	}

	if(run_interactive){
		cout << "Give the name of the exon list to search\n";
		cin >> exon_list_filename;
		cout << "Give the name of the intron list to search\n";
		cin >> intron_list_filename;
		cout << "Give the name of the regulatory network to output\n";
		cin >> reg_net_filename;

	}
	ofstream net_file(reg_net_filename); // FIX: Removed .c_str()
	
	if(!net_file.is_open()){
		cout << "Error opening regulatory network file: " << reg_net_filename << "\n";
		return -2;
	}
	
	load_map_file(exon_list_filename, exon_list);
	load_map_file(intron_list_filename, intron_list);
	
	cout << exon_list.size() << " unique exons loaded from file\n";
	cout << intron_list.size() << " unique introns loaded from file\n";
	
	//now load the genbank file
	gnSequence seq_file;
	if(run_interactive){
		cout << "Enter the name of the genbank sequence file you are using\n";
		cin >> seq_filename;
	}
	if(!seq_file.LoadSource(seq_filename)){
		cout << "Error loading file: " << seq_filename << "\n";
		return -1;
	}
	cout << "Sequence loaded successfully, " << seq_file.length() << " base pairs.\n";
	
	//construct a mapping between coordinates...
	map_coordinates(exon_list, exon_map_list);
	map_coordinates(intron_list, intron_map_list);
	
	//now read the mem file
	while(mem_file.good()){
		ExMem m;
		// Check for successful read before processing.
		if (!(mem_file >> m.length && mem_file >> m.ex_start && mem_file >> m.in_start)) {
			break; // break loop on read failure
		}
		
		if(m.length >= minimum_match_size)
			mem_list.push_back(m);
	}
	cout << mem_list.size() << " matches loaded.\n";
	//sort the mem list
	// FIX: Use std::sort with the custom comparison function
	sort(mem_list.begin(), mem_list.end(), &ExMemLessthan);
	
	//now get the intersection for each mem in the list...
	uint32 exonmapI = 0;
	uint32 notify_percent = 10;
	// FIX: Check for mem_list.size() == 0 to prevent division by zero
	uint32 notify_interval = mem_list.empty() ? 1 : mem_list.size() / notify_percent;
	uint32 percent_count = 0;
	cout << "Searching for complementary matches:\n";
	for(uint32 memI = 0; memI < mem_list.size(); memI++){
		if(notify_interval > 0 && memI % notify_interval == 0){
			cout << percent_count << "%.. ";
			percent_count += notify_percent;
		}
		//simple linear search for intersecting exon mapping
		// FIX: Add namespace
		genome::gnLocation ex_map_loc(mem_list[memI].ex_start, mem_list[memI].ex_start + mem_list[memI].length - 1);
		for(; exonmapI < exon_map_list.size(); exonmapI++){
			if(exon_map_list[exonmapI].Intersects(ex_map_loc))
				break;
		}
		
		//continue to search for any other mappings that intersect
		uint32 mapEnd = exonmapI;
		for(; mapEnd < exon_map_list.size(); mapEnd++){
			// FIX: Should be mapEnd, not exonmapI in this inner loop condition
			if(!exon_map_list[mapEnd].Intersects(ex_map_loc)) 
				break;
		}
		mapEnd--; // back up one, since 'break' happened after an increment
		
		uint32 intronmapI = 0; // intronmapI will contain the index of the first intersecting intron
		//do a binary search for intersecting intron mappings
		// FIX: Use std::abs from <cmath> or the custom one defined above if needed
		int64 cur_in_start = mem_list[memI].in_start;
		if(cur_in_start < 0)
			// The custom abs function should be used here if necessary, otherwise std::abs
			cur_in_start = abs(cur_in_start); 
			
		// FIX: Add namespace
		genome::gnLocation in_map_loc(cur_in_start, cur_in_start + mem_list[memI].length - 1);
		
		uint32 search_mapI = intron_map_list.size() > 0 ? 
		                     loc_binary_search(intron_map_list, 0, intron_map_list.size()-1, in_map_loc) : 0;
		intronmapI = search_mapI;

		//search backwards for previous intersections
		// FIX: This loop condition 'intronmapI >= 0' on an unsigned int (uint32) 
		// is an infinite loop when it wraps below zero. Must check for > 0 first.
		for(intronmapI = search_mapI; intronmapI > 0; intronmapI--){
			if(!intron_map_list[intronmapI-1].Intersects(in_map_loc)){
				break; // Found the start, index i is the first intersection
			}
		}
		
		//continue to search for any other mappings that intersect
		uint32 intron_mapEnd = search_mapI;
		for(; intron_mapEnd < intron_map_list.size(); intron_mapEnd++){
			// FIX: Should be intron_mapEnd, not intronmapI in this inner loop condition
			if(!intron_map_list[intron_mapEnd].Intersects(in_map_loc))
				break;
		}
		intron_mapEnd--;
		
		//we have the mappings, now map the coordinates
		vector<uint32> ex_feat_index, in_feat_index;
		vector<gnBaseFeature*> ex_feat_list;
		vector<gnBaseFeature*> in_feat_list;
		gnSeqI cur_match_len = mem_list[memI].length;
		
		//find out how much of the first exon was matched
		//extra exon start has the number of unmatched bases at the beginning of the exon
		gnSeqI extra_exon_start = mem_list[memI].ex_start - exon_map_list[exonmapI].GetStart();
		gnSeqI cur_exon_len = exon_map_list[exonmapI].GetEnd() - exon_map_list[exonmapI].GetStart() + 1;
		gnSeqI max_exon_chunk = cur_exon_len - extra_exon_start;
		gnSeqI cur_exon_chunk = max_exon_chunk < cur_match_len ? max_exon_chunk : cur_match_len;
		
		//find out how much of the first intron was matched
		gnSeqI extra_intron_start, cur_intron_len, max_intron_chunk, cur_intron_chunk;
		boolean complement = false;
		if(mem_list[memI].in_start > 0){
			extra_intron_start = mem_list[memI].in_start - intron_map_list[intronmapI].GetStart();
			cur_intron_len = intron_map_list[intronmapI].GetEnd() - intron_map_list[intronmapI].GetStart() + 1;
			max_intron_chunk = cur_intron_len - extra_intron_start;
			cur_intron_chunk = max_intron_chunk < mem_list[memI].length ? max_intron_chunk : mem_list[memI].length;
		}else{
			//reverse complement, start at the end.
			if(cur_in_start >= intron_map_list[intron_mapEnd].GetStart())
				cur_intron_chunk = cur_match_len;
			else
				// FIX: This calculation seems flawed for a reverse complement and might need review. 
				// Assuming the original logic:
				cur_intron_chunk = cur_in_start + cur_match_len - intron_map_list[intron_mapEnd].GetStart();
			complement = true;
			// FIX: Add namespace
			seq_file.getIntersectingFeatures(intron_list[intronmapI], in_feat_list, in_feat_index);
		}
		
		//the current chunk will be the smaller of the two mappings
		gnSeqI cur_chunk = cur_intron_chunk < cur_exon_chunk ? cur_intron_chunk : cur_exon_chunk;

		// FIX: Add namespace
		genome::gnLocation cur_exon_loc(exon_list[exonmapI].GetStart() + extra_exon_start, exon_list[exonmapI].GetStart() + cur_chunk);
		seq_file.getIntersectingFeatures(cur_exon_loc, ex_feat_list, ex_feat_index);
		if(mem_list[memI].in_start > 0){
			// FIX: Add namespace
			genome::gnLocation cur_intron_loc(intron_list[intronmapI].GetStart() - 1, intron_list[intronmapI].GetEnd() + 1);
			seq_file.getIntersectingFeatures(cur_intron_loc, in_feat_list, in_feat_index);
		}else{
			// FIX: Add namespace
			genome::gnLocation cur_intron_loc(intron_list[intron_mapEnd].GetStart() - 1, intron_list[intron_mapEnd].GetEnd() + 1);
			seq_file.getIntersectingFeatures(cur_intron_loc, in_feat_list, in_feat_index);
		}
		vector<gnBaseFeature*> ex_forward;
		vector<gnBaseFeature*> ex_reverse;
		vector<gnBaseFeature*> in_forward;
		vector<gnBaseFeature*> in_reverse;

		for(uint32 featI = 0; featI < ex_feat_list.size(); featI++){
			string featName = ex_feat_list[featI]->GetName();
			if(featName == "mRNA" || featName == "CDS" || featName == "gene" )
				// FIX: Add namespace to gnLocation::LT_Complement
				if(ex_feat_list[featI]->GetLocationType() == genome::gnLocation::LT_Complement)
					ex_reverse.push_back(ex_feat_list[featI]);
				else
					ex_forward.push_back(ex_feat_list[featI]);
		}
		for(uint32 featI = 0; featI < in_feat_list.size(); featI++){
			string featName = in_feat_list[featI]->GetName();
			if(featName == "mRNA" || featName == "CDS" || featName == "gene" )
				// FIX: Add namespace to gnLocation::LT_Complement
				if(in_feat_list[featI]->GetLocationType() == genome::gnLocation::LT_Complement)
					in_reverse.push_back(in_feat_list[featI]);
				else
					in_forward.push_back(in_feat_list[featI]);
		}

		if(complement){
			// FIX: Use modern C++ swap for vector content
			std::swap(in_forward, in_reverse);
		}

		//now print out all the complementary features
		if((ex_forward.size() > 0 && in_reverse.size() > 0) || (in_forward.size() > 0 && ex_reverse.size() > 0)){
			net_file << "================================\n";
			net_file << "Mem: " << mem_list[memI].length << "\n";
			net_file << "This exon/intron matching size: " << cur_chunk << "\n";
		}
		if(ex_forward.size() > 0 && in_reverse.size() > 0){
			net_file << "Forward Exons:\n";
			for(uint32 featI = 0; featI < ex_forward.size(); featI++)
				print_feature(net_file, ex_forward[featI]);
			net_file << "Matching introns:\n";
			for(uint32 featI = 0; featI < in_reverse.size(); featI++)
				print_feature(net_file, in_reverse[featI]);
		}
		if(in_forward.size() > 0 && ex_reverse.size() > 0){
			net_file << "Reverse Exons:\n";
			for(uint32 featI = 0; featI < ex_reverse.size(); featI++)
				print_feature(net_file, ex_reverse[featI]);
			net_file << "Matching introns:\n";
			for(uint32 featI = 0; featI < in_forward.size(); featI++)
				print_feature(net_file, in_forward[featI]);
		}
		
		//release memory
		// FIX: Use a cleanup function or smart pointers in modern code, but stick to original logic here
		for(uint32 featI = 0; featI < ex_feat_list.size(); featI++)
			delete ex_feat_list[featI];
		for(uint32 featI = 0; featI < in_feat_list.size(); featI++)
			delete in_feat_list[featI];
		
		//loop while there is stuff to match
//		while(cur_match_len > 0){
		
//		}
	}
	
	return 0; // Standard practice to return 0 on success
}
