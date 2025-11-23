#include "libMUSCLE/muscle.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/msa.h"

namespace muscle {

/***
Compute weights by the CLUSTALW method.
Thompson, Higgins and Gibson (1994), CABIOS (10) 19-29;
see also CLUSTALW paper.

Weights are computed from the edge lengths of a rooted tree.

Define the strength of an edge to be its length divided by the number
of leaves under that edge. The weight of a sequence is then the sum
of edge strengths on the path from the root to the leaf.

Example.

        0.2
       -----A     0.1
	 -x         ------- B     0.7
	   --------y           ----------- C
	    0.3     ----------z
                    0.4    -------------- D
                                 0.8

Edge	Length	Leaves	Strength
----	-----	------	--------
xy		0.3		3		0.1
xA		0.2		1		0.2
yz		0.4		2		0.2
yB		0.1		1		0.1
zC		0.7		1		0.7
zD		0.8		1		0.8

Leaf	Path		Strengths			Weight
----	----		---------			------
A		xA			0.2					0.2
B		xy-yB		0.1 + 0.1			0.2
C		xy-yz-zC	0.1 + 0.2 + 0.7		1.0
D		xy-yz-zD	0.1 + 0.2 + 0.8		1.1

***/

#define TRACE 0

static unsigned CountLeaves(const Tree &tree, unsigned uNodeIndex,
  unsigned LeavesUnderNode[])
	{
	if (tree.IsLeaf(uNodeIndex))
		{
		LeavesUnderNode[uNodeIndex] = 1;
		return 1;
		}

	const unsigned uLeft = tree.GetLeft(uNodeIndex);
	const unsigned uRight = tree.GetRight(uNodeIndex);
	const unsigned uRightCount = CountLeaves(tree, uRight, LeavesUnderNode);
	const unsigned uLeftCount = CountLeaves(tree, uLeft, LeavesUnderNode);
	const unsigned uCount = uRightCount + uLeftCount;
	LeavesUnderNode[uNodeIndex] = uCount;
	return uCount;
	}

void CalcClustalWWeights(const Tree &tree, WEIGHT Weights[])
	{
	std::cerr << "DEBUG CalcClustalWWeights: Starting" << std::endl;

	const unsigned uLeafCount = tree.GetLeafCount();
	std::cerr << "DEBUG CalcClustalWWeights: uLeafCount=" << uLeafCount << std::endl;
	
	if (0 == uLeafCount)
		{
		std::cerr << "DEBUG CalcClustalWWeights: uLeafCount==0, returning" << std::endl;
		return;
		}
	else if (1 == uLeafCount)
		{
		std::cerr << "DEBUG CalcClustalWWeights: uLeafCount==1, setting weight to 1.0" << std::endl;
		Weights[0] = (WEIGHT) 1.0;
		return;
		}
	else if (2 == uLeafCount)
		{
		std::cerr << "DEBUG CalcClustalWWeights: uLeafCount==2, setting weights to 0.5 each" << std::endl;
		Weights[0] = (WEIGHT) 0.5;
		Weights[1] = (WEIGHT) 0.5;
		return;
		}

	std::cerr << "DEBUG CalcClustalWWeights: uLeafCount>2, checking if rooted" << std::endl;
	if (!tree.IsRooted())
		{
		std::cerr << "ERROR CalcClustalWWeights: Tree is not rooted!" << std::endl;
		Quit("CalcClustalWWeights requires rooted tree");
		}
	std::cerr << "DEBUG CalcClustalWWeights: Tree is rooted" << std::endl;

	const unsigned uNodeCount = tree.GetNodeCount();
	std::cerr << "DEBUG CalcClustalWWeights: uNodeCount=" << uNodeCount << std::endl;
	
	unsigned *LeavesUnderNode = new unsigned[uNodeCount];
	memset(LeavesUnderNode, 0, uNodeCount*sizeof(unsigned));

	const unsigned uRootNodeIndex = tree.GetRootNodeIndex();
	std::cerr << "DEBUG CalcClustalWWeights: uRootNodeIndex=" << uRootNodeIndex << std::endl;
	std::cerr << "DEBUG CalcClustalWWeights: Calling CountLeaves" << std::endl;
	
	unsigned uLeavesUnderRoot = CountLeaves(tree, uRootNodeIndex, LeavesUnderNode);
	std::cerr << "DEBUG CalcClustalWWeights: uLeavesUnderRoot=" << uLeavesUnderRoot << std::endl;
	
	if (uLeavesUnderRoot != uLeafCount)
		Quit("WeightsFromTreee: Internal error, root count %u %u",
		  uLeavesUnderRoot, uLeafCount);

	double *Strengths = new double[uNodeCount];
	std::cerr << "DEBUG CalcClustalWWeights: Calculating strengths" << std::endl;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (tree.IsRoot(uNodeIndex))
			{
			Strengths[uNodeIndex] = 0.0;
			continue;
			}
		const unsigned uParent = tree.GetParent(uNodeIndex);
		const double dLength = tree.GetEdgeLength(uNodeIndex, uParent);
		const unsigned uLeaves = LeavesUnderNode[uNodeIndex];
		const double dStrength = dLength / (double) uLeaves;
		Strengths[uNodeIndex] = dStrength;
		}

	std::cerr << "DEBUG CalcClustalWWeights: Calculating weights for each leaf" << std::endl;
	for (unsigned n = 0; n < uLeafCount; ++n)
		{
		const unsigned uLeafNodeIndex = tree.LeafIndexToNodeIndex(n);
		if (!tree.IsLeaf(uLeafNodeIndex))
			Quit("CalcClustalWWeights: leaf");

		double dWeight = 0;
		unsigned uNode = uLeafNodeIndex;
		while (!tree.IsRoot(uNode))
			{
			dWeight += Strengths[uNode];
			uNode = tree.GetParent(uNode);
			}
		if (dWeight < 0.0001)
			dWeight = 1.0;
		Weights[n] = (WEIGHT) dWeight;
		}

	std::cerr << "DEBUG CalcClustalWWeights: Cleaning up" << std::endl;
	delete[] Strengths;
	delete[] LeavesUnderNode;

	std::cerr << "DEBUG CalcClustalWWeights: Normalizing weights" << std::endl;
	Normalize(Weights, uLeafCount);
	std::cerr << "DEBUG CalcClustalWWeights: Completed successfully" << std::endl;
	}
#if	TRACE
		Log(" = %g\n", dWeight);
#endif
		}

	delete[] Strengths;
	delete[] LeavesUnderNode;

	Normalize(Weights, uLeafCount);
	}

void MSA::SetClustalWWeights(const Tree &tree)
	{
	const unsigned uSeqCount = GetSeqCount();
	const unsigned uLeafCount = tree.GetLeafCount();

	WEIGHT *Weights = new WEIGHT[uSeqCount];

	CalcClustalWWeights(tree, Weights);

	for (unsigned n = 0; n < uLeafCount; ++n)
		{
		const WEIGHT w = Weights[n];
		const unsigned uLeafNodeIndex = tree.LeafIndexToNodeIndex(n);
		const unsigned uId = tree.GetLeafId(uLeafNodeIndex);
		const unsigned uSeqIndex = GetSeqIndex(uId);
#if	DEBUG
		if (GetSeqName(uSeqIndex) != tree.GetLeafName(uLeafNodeIndex))
			Quit("MSA::SetClustalWWeights: names don't match");
#endif
		SetSeqWeight(uSeqIndex, w);
		}
	NormalizeWeights((WEIGHT) 1.0);

	delete[] Weights;
	}
} 
