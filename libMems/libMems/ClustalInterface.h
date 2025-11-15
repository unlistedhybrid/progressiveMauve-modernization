/*******************************************************************************
 * $Id: ClustalInterface.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * Copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Licensed under the GPL. See COPYING for details.
 ******************************************************************************/

#ifndef _ClustalInterface_h_
#define _ClustalInterface_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/NumericMatrix.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnSequence.h"
#include "libMems/GappedAlignment.h"
#include "libMems/GappedAligner.h"

// attempt to auto-link the ClustalW library on windows
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalW64omp.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalW64fdomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalWomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalWfdomp.lib")
#endif
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW64.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW64fd.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalWfd.lib")
#endif

namespace mems {

class ClustalInterface : public GappedAligner {
public:
	static ClustalInterface& getClustalInterface();
	bool Align(GappedAlignment& cr, Match* r_begin, Match* r_end, std::vector< genome::gnSequence* >& seq_table);
	void SetDistanceMatrix(NumericMatrix<double>& distance_matrix, std::string& tree_filename);
	void SetMinFlankSize(gnSeqI min_flank) { min_flank_size = min_flank; }

	void setGuideTree(std::string& tree_filename, NumericMatrix<double>& dist_mat, uint seq_count);

	bool guideTreeLoaded() const { return distance_matrix.cols() > 0; };

	void SetDistanceMatrix(NumericMatrix<double>& distance_matrix, std::string& tree_filename, bool reread_tree);

protected:
	bool CallClustal(std::vector< std::string >& seq_table);
	NumericMatrix<double> distance_matrix;
	gnSeqI min_flank_size;
	int clustal_score_cutoff;
	bool allocated_aln;
private:
	ClustalInterface(const ClustalInterface& ci) { *this = ci; }
	ClustalInterface& operator=(const ClustalInterface& ci);
	ClustalInterface();
};

}

#endif // _ClustalInterface_h_
