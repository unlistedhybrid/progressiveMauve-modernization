#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/ParallelMemHash.h"
#include <vector>

#ifdef _OPENMP

using namespace std;
using namespace genome;
namespace mems {

ParallelMemHash::ParallelMemHash() : MemHash() {}

ParallelMemHash::ParallelMemHash(const ParallelMemHash& mh) : MemHash(mh) {
    *this = mh;
}

ParallelMemHash& ParallelMemHash::operator=(const ParallelMemHash& mh) {
    thread_mem_table = mh.thread_mem_table;
    return *this;
}

ParallelMemHash* ParallelMemHash::Clone() const {
    return new ParallelMemHash(*this);
}

void ParallelMemHash::FindMatches(MatchList& ml) {
    for (uint32 seqI = 0; seqI < ml.seq_table.size(); ++seqI) {
        if (!AddSequence(ml.sml_table[seqI], ml.seq_table[seqI])) {
            ErrorMsg("Error adding " + ml.seq_filename[seqI] + "\n");
            return;
        }
    }

    size_t CHUNK_SIZE = 200000;
    vector<vector<gnSeqI>> chunk_starts;
    vector<gnSeqI> chunk_lengths;

    mers_processed = 0;
    total_mers = 0;
    m_progress = -1;
    for (size_t i = 0; i < ml.sml_table.size(); i++)
        total_mers += ml.sml_table[i]->Length();

    int max_length_sml = -1;
    size_t maxlen = 0;
    for (size_t i = 0; i < ml.sml_table.size(); i++) {
        if (ml.sml_table[i]->Length() > maxlen) {
            maxlen = ml.sml_table[i]->Length();
            max_length_sml = i;
        }
    }

    chunk_starts.push_back(vector<gnSeqI>(seq_count, 0));

    while (chunk_starts.back()[max_length_sml] + CHUNK_SIZE < ml.sml_table[max_length_sml]->Length()) {
        vector<gnSeqI> tmp(seq_count, 0);
        GetBreakpoint(max_length_sml, chunk_starts.back()[max_length_sml] + CHUNK_SIZE, tmp);
        chunk_starts.push_back(tmp);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(chunk_starts.size()); i++) {
        if (thread_mem_table.get().size() != mem_table.size())
            thread_mem_table.get().resize(mem_table.size());

        vector<gnSeqI> chunk_lens(seq_count);
        if (i + 1 < static_cast<int>(chunk_starts.size())) {
            for (size_t j = 0; j < seq_count; j++)
                chunk_lens[j] = chunk_starts[i + 1][j] - chunk_starts[i][j];
        } else
            chunk_lens = vector<gnSeqI>(seq_count, GNSEQI_END);
        SearchRange(chunk_starts[i], chunk_lens);
        MergeTable();
    }
    GetMatchList(ml);
}

void ParallelMemHash::MergeTable() {
#pragma omp critical
    {
        size_t buckets = thread_mem_table.get().size();
        for (size_t bI = 0; bI < buckets; bI++) {
            vector<MatchHashEntry*>& bucket = thread_mem_table.get()[bI];
            for (size_t mI = 0; mI < bucket.size(); mI++) {
                MemHash::AddHashEntry((*(bucket[mI])), mem_table);
            }
        }
        thread_mem_table.get() = mem_table;
    }
}

MatchHashEntry* ParallelMemHash::AddHashEntry(MatchHashEntry& mhe) {
    return MemHash::AddHashEntry(mhe, thread_mem_table.get());
}

} // namespace mems

#endif // _OPENMP
