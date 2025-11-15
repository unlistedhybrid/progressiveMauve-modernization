#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>

#include "gnBaseFeature.h"
#include "gnLocation.h"
#include "gnSequence.h"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::list;

struct ExMem{
    gnSeqI length;
    int64 ex_start;
    int64 in_start;
};

static boolean LocationLessthan(const gnLocation& a, const gnLocation& b){
    return a.GetStart() < b.GetStart();
}

static boolean LocationEndLessthan(const gnLocation& a, const gnLocation& b){
    return a.GetEnd() < b.GetStart();
}

static boolean LocationSizeLessthan(const gnLocation& a, const gnLocation& b){
    return (a.GetEnd() - a.GetStart()) < (b.GetEnd() - b.GetStart());
}

static boolean ExMemLessthan(const ExMem& a, const ExMem& b){
    if(a.ex_start == b.ex_start){
        if(a.in_start == b.in_start)
            return a.length < b.length;
        return a.in_start < b.in_start;
    }
    return a.ex_start < b.ex_start;
}

void print_usage(char* pname){
    cout << "Usage: " << pname
        << " <genbank file> <exon_list> <intron_list> <mem_list>"
        << " <regulation_network> <minimum match size>\n";
}

void load_map_file(string& filename, vector<gnLocation>& loc_list){
    ifstream coord_file(filename.c_str());
    if(!coord_file.is_open()){
        cout << "couldn't open file: " << filename << "\n";
        return;
    }

    while(coord_file.good()){
        gnSeqI start_coord, end_coord;
        coord_file >> start_coord;
        coord_file >> end_coord;

        if(coord_file.good()){
            loc_list.push_back(gnLocation(start_coord, end_coord));
        }
    }
}

void map_coordinates(vector<gnLocation>& loc_list, vector<gnLocation>& map_list){
    gnSeqI curStart = 1;
    for(uint32 i = 0; i < loc_list.size(); i++){
        gnSeqI cur_len = loc_list[i].GetEnd() - loc_list[i].GetStart() + 1;
        map_list.push_back(gnLocation(curStart, curStart + cur_len - 1));
        curStart += cur_len;
    }
}

void map_list_coordinates(list<gnLocation>& loc_list, list<gnLocation>& map_list){
    gnSeqI curStart = 1;
    for(list<gnLocation>::iterator it = loc_list.begin();
        it != loc_list.end(); ++it){
        gnSeqI cur_len = it->GetEnd() - it->GetStart() + 1;
        map_list.push_back(gnLocation(curStart, curStart + cur_len - 1));
        curStart += cur_len;
    }
}

uint32 loc_binary_search(vector<gnLocation>& loc_list,
                         uint32 startI, uint32 endI,
                         gnLocation& query_loc)
{
    uint32 midI = ((endI - startI) / 2) + startI;
    if(startI == endI)
        return endI;

    if(loc_list[midI].Intersects(query_loc)){
        return midI;
    }else if(loc_list[midI].GetStart() < query_loc.GetStart()){
        return loc_binary_search(loc_list, midI + 1, endI, query_loc);
    }else{
        return loc_binary_search(loc_list, startI, midI, query_loc);
    }
}

int main(int argc, char* argv[]){

    boolean run_interactive = false;
    string seq_filename;
    string exon_list_filename;
    string intron_list_filename;
    string mem_list_filename;
    string reg_net_filename;

    vector<gnLocation> exon_list;
    vector<gnLocation> intron_list;
    vector<gnLocation> exon_map_list;
    vector<gnLocation> intron_map_list;

    vector<ExMem> mem_list;
    uint32 minimum_match_size;

    if(argc != 7){
        print_usage(argv[0]);
        return -1;
    }

    seq_filename = argv[1];
    exon_list_filename = argv[2];
    intron_list_filename = argv[3];
    mem_list_filename = argv[4];
    reg_net_filename = argv[5];
    minimum_match_size = atoi(argv[6]);

    ifstream mem_file(mem_list_filename.c_str());
    if(!mem_file.is_open()){
        cout << "Error opening " << mem_list_filename << "\n";
        return -1;
    }

    ofstream net_file(reg_net_filename.c_str());
    if(!net_file.is_open()){
        cout << "Error opening regulatory network file: "
             << reg_net_filename << "\n";
        return -2;
    }

    load_map_file(exon_list_filename, exon_list);
    load_map_file(intron_list_filename, intron_list);

    cout << exon_list.size() << " unique exons loaded\n";
    cout << intron_list.size() << " unique introns loaded\n";

    gnSequence seq_file;
    if(!seq_file.LoadSource(seq_filename)){
        cout << "Error loading genbank file\n";
        return -1;
    }

    cout << "Sequence loaded, " << seq_file.length() << " bp\n";

    map_coordinates(exon_list, exon_map_list);
    map_coordinates(intron_list, intron_map_list);

    while(mem_file.good()){
        ExMem m;
        mem_file >> m.length;
        mem_file >> m.ex_start;
        mem_file >> m.in_start;

        if(mem_file.good() && m.length >= minimum_match_size)
            mem_list.push_back(m);
    }

    cout << mem_list.size() << " MEMs loaded\n";

    sort(mem_list.begin(), mem_list.end(), &ExMemLessthan);

    uint32 exonmapI = 0;
    uint32 notify_percent = 10;
    uint32 notify_interval =
        (mem_list.size() > notify_percent ? mem_list.size() / notify_percent : 1);
    uint32 percent_count = 0;

    cout << "Searching for complementary matches:\n";

    // Reintroduced variables (previously inside removed block)
    gnLocation ex_map_loc(0, 0);
    gnLocation in_map_loc(0, 0);

    for(uint32 memI = 0; memI < mem_list.size(); memI++){

        if(memI % notify_interval == 0){
            cout << percent_count << "%.. ";
            percent_count += notify_percent;
        }

        // Compute exon mapping location
        ex_map_loc = gnLocation(mem_list[memI].ex_start,
                                mem_list[memI].ex_start +
                                mem_list[memI].length - 1);

        for(; exonmapI < exon_map_list.size(); exonmapI++){
            if(exon_map_list[exonmapI].Intersects(ex_map_loc))
                break;
        }

        uint32 mapEnd = exonmapI;
        for(; mapEnd < exon_map_list.size(); mapEnd++){
            if(!exon_map_list[exonmapI].Intersects(ex_map_loc))
                break;
        }
        mapEnd--;

        int64 cur_in_start = mem_list[memI].in_start;
        if(cur_in_start < 0)
            cur_in_start = -cur_in_start;

        in_map_loc = gnLocation(cur_in_start,
                                cur_in_start + mem_list[memI].length - 1);

        uint32 search_mapI =
            loc_binary_search(intron_map_list, 0,
                              intron_map_list.size() - 1, in_map_loc);

        uint32 intronmapI = search_mapI;

        for(; intronmapI >= 0; intronmapI--){
            if(!intron_map_list[intronmapI].Intersects(in_map_loc)){
                intronmapI++;
                break;
            }
            if(intronmapI == 0)
                break;
        }

        uint32 intron_mapEnd = search_mapI;
        for(; intron_mapEnd < intron_map_list.size(); intron_mapEnd++){
            if(!intron_map_list[intronmapI].Intersects(in_map_loc))
                break;
        }
        intron_mapEnd--;
    }
    return 0;
}
