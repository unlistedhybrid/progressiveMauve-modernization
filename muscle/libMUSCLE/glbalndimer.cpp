#include "libMUSCLE/muscle.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "libMUSCLE/pwpath.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/gapscoredimer.h"
#include "libMUSCLE/threadstorage.h"

namespace muscle {

#define TRACE 0

static SCORE TraceBackDimer(const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  const char *TBM_, const char *TBD_, const char *TBI_,
  unsigned uLengthA, unsigned uLengthB, PWPath &Path);

static const char *LocalScoreToStr(SCORE s)
{
    static TLS<char[16]> str;
    if (MINUS_INFINITY == s)
        return "     *";
    snprintf(str.get(), 16, "%6.3g", s);
    return str.get();
}

static ProfPos getInitedPPTerm()
{
    ProfPos pp;
    pp.m_bAllGaps = false;
    pp.m_LL = 1;
    pp.m_LG = 0;
    pp.m_GL = 0;
    pp.m_GG = 0;
    pp.m_fOcc = 1;
    return pp;
}
static TLS<ProfPos> PPTerm(getInitedPPTerm());

static SCORE ScoreProfPosDimerLE(const ProfPos &PPA, const ProfPos &PPB)
{
    SCORE Score = 0;
    for (unsigned n = 0; n < 20; ++n)
    {
        const unsigned uLetter = PPA.m_uSortOrder[n];
        const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
        if (0 == fcLetter)
            break;
        Score += fcLetter * PPB.m_AAScores[uLetter];
    }
    if (0 == Score)
        return -2.5;
    SCORE logScore = logf(Score);
    return (SCORE)(logScore * (PPA.m_fOcc * PPB.m_fOcc));
}

static SCORE ScoreProfPosDimerPSP(const ProfPos &PPA, const ProfPos &PPB)
{
    SCORE Score = 0;
    for (unsigned n = 0; n < 20; ++n)
    {
        const unsigned uLetter = PPA.m_uSortOrder[n];
        const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
        if (0 == fcLetter)
            break;
        Score += fcLetter * PPB.m_AAScores[uLetter];
    }
    return Score;
}

static SCORE ScoreProfPosDimer(const ProfPos &PPA, const ProfPos &PPB)
{
    switch (g_PPScore.get())
    {
    case PPSCORE_LE:
        return ScoreProfPosDimerLE(PPA, PPB);
    case PPSCORE_SP:
    case PPSCORE_SV:
        return ScoreProfPosDimerPSP(PPA, PPB);
    case PPSCORE_Undefined:
    default:
        Quit("Invalid g_PPScore.get()");
        return 0;
    }
}

SCORE GlobalAlignDimer(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
{
    assert(uLengthB > 0 && uLengthA > 0);

    const unsigned uPrefixCountA = uLengthA + 1;
    const unsigned uPrefixCountB = uLengthB + 1;

    const size_t LM = uPrefixCountA * uPrefixCountB;
    SCORE *DPM_ = new SCORE[LM];
    SCORE *DPD_ = new SCORE[LM];
    SCORE *DPI_ = new SCORE[LM];

    char *TBM_ = new char[LM];
    char *TBD_ = new char[LM];
    char *TBI_ = new char[LM];

    DPM(0, 0) = 0;
    DPD(0, 0) = MINUS_INFINITY;
    DPI(0, 0) = MINUS_INFINITY;

    TBM(0, 0) = 'S';
    TBD(0, 0) = '?';
    TBI(0, 0) = '?';

    DPM(1, 0) = MINUS_INFINITY;
    DPD(1, 0) = GapScoreMD(PA[0], PPTerm.get());
    DPI(1, 0) = MINUS_INFINITY;

    TBM(1, 0) = '?';
    TBD(1, 0) = 'S';
    TBI(1, 0) = '?';

    DPM(0, 1) = MINUS_INFINITY;
    DPD(0, 1) = MINUS_INFINITY;
    DPI(0, 1) = GapScoreMI(PPTerm.get(), PB[0]);

    TBM(0, 1) = '?';
    TBD(0, 1) = '?';
    TBI(0, 1) = 'S';

    for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
    {
        DPM(uPrefixLengthA, 0) = MINUS_INFINITY;
        TBM(uPrefixLengthA, 0) = '?';

        DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0) +
            GapScoreDD(PA[uPrefixLengthA - 1], PPTerm.get());
        TBD(uPrefixLengthA, 0) = 'D';

        DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
        TBI(uPrefixLengthA, 0) = '?';
    }

    for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
    {
        DPM(0, uPrefixLengthB) = MINUS_INFINITY;
        TBM(0, uPrefixLengthB) = '?';

        DPD(0, uPrefixLengthB) = MINUS_INFINITY;
        TBD(0, uPrefixLengthB) = '?';

        DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1) +
            GapScoreII(PPTerm.get(), PB[uPrefixLengthB - 1]);
        TBI(0, uPrefixLengthB) = 'I';
    }

    for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
    {
        const ProfPos &PPB = PB[uPrefixLengthB - 1];
        for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
        {
            const ProfPos &PPA = PA[uPrefixLengthA - 1];
            {
                SCORE scoreLL = ScoreProfPosDimer(PPA, PPB);

                SCORE scoreMM = DPM(uPrefixLengthA - 1, uPrefixLengthB - 1) + GapScoreMM(PPA, PPB);
                SCORE scoreDM = DPD(uPrefixLengthA - 1, uPrefixLengthB - 1) + GapScoreDM(PPA, PPB);
                SCORE scoreIM = DPI(uPrefixLengthA - 1, uPrefixLengthB - 1) + GapScoreIM(PPA, PPB);

                SCORE scoreBest = scoreMM;
                char c = 'M';
                if (scoreDM > scoreBest)
                {
                    scoreBest = scoreDM;
                    c = 'D';
                }
                if (scoreIM > scoreBest)
                {
                    scoreBest = scoreIM;
                    c = 'I';
                }

                DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest + scoreLL;
                TBM(uPrefixLengthA, uPrefixLengthB) = c;
            }
            {
                SCORE scoreMD = DPM(uPrefixLengthA - 1, uPrefixLengthB) + GapScoreMD(PPA, PPB);
                SCORE scoreDD = DPD(uPrefixLengthA - 1, uPrefixLengthB) + GapScoreDD(PPA, PPB);
                SCORE scoreID = DPI(uPrefixLengthA - 1, uPrefixLengthB) + GapScoreID(PPA, PPB);

                SCORE scoreBest = scoreMD;
                char c = 'M';
                if (scoreDD > scoreBest)
                {
                    scoreBest = scoreDD;
                    c = 'D';
                }
                if (scoreID > scoreBest)
                {
                    scoreBest = scoreID;
                    c = 'I';
                }

                DPD(uPrefixLengthA, uPrefixLengthB) = scoreBest;
                TBD(uPrefixLengthA, uPrefixLengthB) = c;
            }
            {
                SCORE scoreMI = DPM(uPrefixLengthA, uPrefixLengthB - 1) + GapScoreMI(PPA, PPB);
                SCORE scoreDI = DPD(uPrefixLengthA, uPrefixLengthB - 1) + GapScoreDI(PPA, PPB);
                SCORE scoreII = DPI(uPrefixLengthA, uPrefixLengthB - 1) + GapScoreII(PPA, PPB);

                SCORE scoreBest = scoreMI;
                char c = 'M';
                if (scoreDI > scoreBest)
                {
                    scoreBest = scoreDI;
                    c = 'D';
                }
                if (scoreII > scoreBest)
                {
                    scoreBest = scoreII;
                    c = 'I';
                }

                DPI(uPrefixLengthA, uPrefixLengthB) = scoreBest;
                TBI(uPrefixLengthA, uPrefixLengthB) = c;
            }
        }
    }

    SCORE Score = TraceBackDimer(DPM_, DPD_, DPI_, TBM_, TBD_, TBI_,
        uLengthA, uLengthB, Path);

    delete[] DPM_;
    delete[] DPD_;
    delete[] DPI_;

    delete[] TBM_;
    delete[] TBD_;
    delete[] TBI_;

    return Score;
}

static SCORE TraceBackDimer(const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
    const char *TBM_, const char *TBD_, const char *TBI_,
    unsigned uLengthA, unsigned uLengthB, PWPath &Path)
{
    const unsigned uPrefixCountA = uLengthA + 1;

    unsigned uPrefixLengthA = uLengthA;
    unsigned uPrefixLengthB = uLengthB;

    char cEdge = 'M';
    SCORE scoreMax = DPM(uLengthA, uLengthB);
    if (DPD(uLengthA, uLengthB) > scoreMax)
    {
        scoreMax = DPD(uLengthA, uLengthB);
        cEdge = 'D';
    }
    if (DPI(uLengthA, uLengthB) > scoreMax)
    {
        scoreMax = DPI(uLengthA, uLengthB);
        cEdge = 'I';
    }

    for (;;)
    {
        if (0 == uPrefixLengthA && 0 == uPrefixLengthB)
            break;

        PWEdge Edge;
        Edge.cType = cEdge;
        Edge.uPrefixLengthA = uPrefixLengthA;
        Edge.uPrefixLengthB = uPrefixLengthB;
        Path.PrependEdge(Edge);

        switch (cEdge)
        {
        case 'M':
            assert(uPrefixLengthA > 0 && uPrefixLengthB > 0);
            cEdge = TBM(uPrefixLengthA, uPrefixLengthB);
            --uPrefixLengthA;
            --uPrefixLengthB;
            break;
        case 'D':
            assert(uPrefixLengthA > 0);
            cEdge = TBD(uPrefixLengthA, uPrefixLengthB);
            --uPrefixLengthA;
            break;
        case 'I':
            assert(uPrefixLengthB > 0);
            cEdge = TBI(uPrefixLengthA, uPrefixLengthB);
            --uPrefixLengthB;
            break;
        default:
            Quit("Invalid edge PLA=%u PLB=%u %c", uPrefixLengthA, uPrefixLengthB, cEdge);
        }
    }
    return scoreMax;
}

}
