#include "libMUSCLE/muscle.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/profile.h"
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <cstdio>
#include <cstring>

namespace muscle {

#define TRACE 0

static unsigned SubFamRecurse(const Tree &tree, unsigned uNodeIndex, unsigned uMaxLeafCount,
  unsigned SubFams[], unsigned &uSubFamCount)
{
    if (tree.IsLeaf(uNodeIndex))
        return 1;

    unsigned uLeft = tree.GetLeft(uNodeIndex);
    unsigned uRight = tree.GetRight(uNodeIndex);
    unsigned uLeftCount = SubFamRecurse(tree, uLeft, uMaxLeafCount, SubFams, uSubFamCount);
    unsigned uRightCount = SubFamRecurse(tree, uRight, uMaxLeafCount, SubFams, uSubFamCount);

    unsigned uLeafCount = uLeftCount + uRightCount;
    if (uLeftCount + uRightCount > uMaxLeafCount)
    {
        if (uLeftCount <= uMaxLeafCount)
            SubFams[uSubFamCount++] = uLeft;
        if (uRightCount <= uMaxLeafCount)
            SubFams[uSubFamCount++] = uRight;
    }
    else if (tree.IsRoot(uNodeIndex))
    {
        if (uSubFamCount != 0)
            Quit("Error in SubFamRecurse");
        SubFams[uSubFamCount++] = uNodeIndex;
    }

    return uLeafCount;
}

void SubFam(const Tree &tree, unsigned uMaxLeafCount, unsigned SubFams[], unsigned *ptruSubFamCount)
{
    *ptruSubFamCount = 0;
    SubFamRecurse(tree, tree.GetRootNodeIndex(), uMaxLeafCount, SubFams, *ptruSubFamCount);
}

static void SetInFam(const Tree &tree, unsigned uNodeIndex, bool NodeInSubFam[])
{
    if (tree.IsLeaf(uNodeIndex))
        return;
    unsigned uLeft = tree.GetLeft(uNodeIndex);
    unsigned uRight = tree.GetRight(uNodeIndex);
    NodeInSubFam[uLeft] = true;
    NodeInSubFam[uRight] = true;

    SetInFam(tree, uLeft, NodeInSubFam);
    SetInFam(tree, uRight, NodeInSubFam);
}

void AlignSubFam(SeqVect &vAll, const Tree &GuideTree, unsigned uNodeIndex,
  MSA &msaOut)
{
    const unsigned uSeqCount = vAll.GetSeqCount();

    const char *InTmp = "asf_in.tmp";
    const char *OutTmp = "asf_out.tmp";

    unsigned *Leaves = new unsigned[uSeqCount];
    unsigned uLeafCount;
    GetLeaves(GuideTree, uNodeIndex, Leaves, &uLeafCount);

    SeqVect v;
    for (unsigned i = 0; i < uLeafCount; ++i)
    {
        unsigned uLeafNodeIndex = Leaves[i];
        unsigned uId = GuideTree.GetLeafId(uLeafNodeIndex);
        Seq &s = vAll.GetSeqById(uId);
        v.AppendSeq(s);
    }

    TextFile fIn(InTmp, true);
    v.ToFASTAFile(fIn);
    fIn.Close();

    char CmdLine[4096];
    snprintf(CmdLine, sizeof(CmdLine), "probcons %s > %s 2> /dev/null", InTmp, OutTmp);
    int ignore = system(CmdLine);
    (void)ignore;

    TextFile fOut(OutTmp);
    msaOut.FromFile(fOut);

    for (unsigned uSeqIndex = 0; uSeqIndex < uLeafCount; ++uSeqIndex)
    {
        const char *Name = msaOut.GetSeqName(uSeqIndex);
        unsigned uId = vAll.GetSeqIdFromName(Name);
        msaOut.SetSeqId(uSeqIndex, uId);
    }

#ifndef _MSC_VER
    unlink(InTmp);
    unlink(OutTmp);
#else
    remove(InTmp);
    remove(OutTmp);
#endif

    delete[] Leaves;
}

}
