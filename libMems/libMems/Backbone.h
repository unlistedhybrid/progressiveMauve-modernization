/*******************************************************************************
 * $Id: Backbone.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * Copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Licensed under the GPL. See COPYING for details.
 ******************************************************************************/

#ifndef __Backbone_h__
#define __Backbone_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/IntervalList.h"
#include "libMems/NumericMatrix.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include <boost/multi_array.hpp>
#include <sstream>
#include <vector>
#include <string>

namespace mems {

using ULA = mems::UngappedLocalAlignment<mems::HybridAbstractMatch<>>;
using backbone_list_t = std::vector<std::vector<ULA*>>;
using pairwise_genome_hss_t = boost::multi_array<std::vector<std::pair<size_t, size_t>>, 3>;

class HssDetector;

/** Compute the GC content of a set of sequences */
double computeGC(std::vector<genome::gnSequence*>& seq_table);

/** Collapse Intervals that are trivially collinear with each other */
void collapseCollinear(IntervalList& iv_list);

/** Sanity checks for alignment columns that contain only gaps */
void checkForAllGapColumns(IntervalList& iv_list);

/**
 * Applies pairwise transitive homology statistics to detect backbone in a single collinear alignment.
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list.
 */
void detectAndApplyBackbone(
    AbstractMatch* m,
    std::vector<genome::gnSequence*>& seq_table,
    CompactGappedAlignment<>*& result,
    backbone_list_t& bb_list,
    const Params& hmm_params,
    boolean left_homologous = false,
    boolean right_homologous = false
);

void detectAndApplyBackbone(
    IntervalList& iv_list,
    backbone_list_t& bb_list,
    const Params& hmm_params
);

void detectBackbone(
    IntervalList& iv_list,
    backbone_list_t& bb_list,
    const HssDetector* detector
);

void writeBackboneColumns(std::ostream& bb_out, backbone_list_t& bb_list);
void writeBackboneSeqCoordinates(backbone_list_t& bb_list, IntervalList& iv_list, std::ostream& bb_out);

class HssDetector {
public:
    using MatchListType = std::vector<CompactGappedAlignment<>*>;
    virtual ~HssDetector() = default;
    virtual void operator()(
        const MatchListType& iv_list,
        std::vector<genome::gnSequence*>& seq_table,
        hss_array_t& hss_array
    ) const = 0;
};

class HomologyHmmDetector : public HssDetector {
public:
    HomologyHmmDetector(const Params& hmm_params, bool left_homologous, bool right_homologous)
        : p(hmm_params), left(left_homologous), right(right_homologous) {}
    void operator()(
        const MatchListType& iv_list,
        std::vector<genome::gnSequence*>& seq_table,
        hss_array_t& hss_array
    ) const override
    {
        findHssHomologyHMM(iv_list, seq_table, hss_array, p, left, right);
    }
private:
    const Params& p;
    bool left;
    bool right;
};

class BigGapsDetector : public HssDetector {
public:
    BigGapsDetector(size_t big_gap_size) : big(big_gap_size) {}
    void operator()(
        const MatchListType& iv_list,
        std::vector<genome::gnSequence*>& seq_table,
        hss_array_t& hss_array
    ) const override
    {
        hss_array_t gap_array;
        findBigGaps(iv_list, seq_table, gap_array, big);
        HssColsToIslandCols(iv_list, seq_table, gap_array, hss_array);
    }
private:
    size_t big;
};

using bb_seqentry_t = std::vector<std::pair<int64, int64>>;
struct bb_entry_t {
    bb_seqentry_t bb_seq;
    ULA bb_cols;
    size_t iv;
};

void addUniqueSegments(std::vector<bb_seqentry_t>& bb_seq_list, size_t min_length = 20);
void mergeAdjacentSegments(std::vector<bb_seqentry_t>& bb_seq_list);

class BbSeqEntrySorter {
public:
    explicit BbSeqEntrySorter(size_t seqI) : m_seq(seqI) {}
    bool operator()(const bb_seqentry_t& a, const bb_seqentry_t& b) const {
        return genome::absolut(a[m_seq].first) < genome::absolut(b[m_seq].first);
    }
    size_t m_seq;
};

// Inline utilities
inline void printBbSeq(std::ostream& os, const bb_seqentry_t& bbseq) {
    for (size_t i = 0; i < bbseq.size(); ++i) {
        if (i > 0)
            os << '\t';
        os << "(" << bbseq[i].first << ", " << bbseq[i].second << ")";
    }
}

inline void readBackboneSeqFile(std::istream& bbseq_input, std::vector<bb_seqentry_t>& backbone) {
    std::string cur_line;
    std::getline(bbseq_input, cur_line); // read off the header line
    while (std::getline(bbseq_input, cur_line)) {
        bb_seqentry_t bb;
        std::stringstream line_str(cur_line);
        int64 lpos = 0;
        while (line_str >> lpos) {
            int64 rpos = 0;
            line_str >> rpos;
            bb.emplace_back(lpos, rpos);
        }
        backbone.push_back(bb);
    }
}

inline void writeBackboneSeqFile(std::ostream& bbseq_out, std::vector<bb_seqentry_t>& backbone) {
    if (backbone.empty())
        return; // can't write if there's no backbone!
    for (size_t seqI = 0; seqI < backbone[0].size(); seqI++) {
        if (seqI > 0)
            bbseq_out << '\t';
        std::stringstream ss;
        ss << "seq" << seqI;
        bbseq_out << ss.str() << "_leftend\t" << ss.str() << "_rightend";
    }
    bbseq_out << std::endl;
    for (const auto& bb : backbone) {
        for (size_t seqI = 0; seqI < bb.size(); seqI++) {
            if (seqI > 0)
                bbseq_out << '\t';
            bbseq_out << bb[seqI].first << '\t' << bb[seqI].second;
        }
        bbseq_out << std::endl;
    }
}

inline void readBackboneColsFile(std::istream& bbcol_input, std::vector<std::pair<size_t, ULA>>& bb_list) {
    std::string cur_line;
    while (std::getline(bbcol_input, cur_line)) {
        ULA tmp_ula;
        size_t ivI;
        std::stringstream ss(cur_line);
        ss >> ivI;
        size_t left_col;
        size_t len;
        ss >> left_col;
        ss >> len;
        gnSeqI bbseq;
        while (ss >> bbseq) {
            tmp_ula.SetStart(bbseq, left_col);
        }
        tmp_ula.SetLength(len);
        bb_list.emplace_back(ivI, tmp_ula);
    }
}

void makeAllPairwiseGenomeHSS(
    IntervalList& iv_list,
    std::vector<CompactGappedAlignment<>*>& iv_ptrs,
    std::vector<CompactGappedAlignment<>*>& iv_orig_ptrs,
    pairwise_genome_hss_t& hss_cols,
    const HssDetector* detector
);

void mergePairwiseHomologyPredictions(
    std::vector<CompactGappedAlignment<>*>& iv_orig_ptrs,
    pairwise_genome_hss_t& hss_cols,
    std::vector<std::vector<ULA*>>& ula_list
);

} // namespace mems

#endif // __Backbone_h__
