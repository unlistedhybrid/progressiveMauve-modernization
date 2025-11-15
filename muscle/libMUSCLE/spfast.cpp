#include "libMUSCLE/muscle.h"
#include "libMUSCLE/profile.h"

namespace muscle {

#define TRACE	0

enum
{
	LL = 0,
	LG = 1,
	GL = 2,
	GG = 3,
};

static const char *GapTypeToStr(int GapType)
{
	switch (GapType)
	{
	case LL: return "LL";
	case LG: return "LG";
	case GL: return "GL";
	case GG: return "GG";
	}
	Quit("Invalid gap type");
	return "?";
}

static TLS<SCORE[4][4]> GapScoreMatrix;

static void InitGapScoreMatrix()
{
	const SCORE t = (SCORE) 0.2;

	GapScoreMatrix.get()[LL][LL] = 0;
	GapScoreMatrix.get()[LL][LG] = g_scoreGapOpen.get();
	GapScoreMatrix.get()[LL][GL] = 0;
	GapScoreMatrix.get()[LL][GG] = 0;

	GapScoreMatrix.get()[LG][LL] = g_scoreGapOpen.get();
	GapScoreMatrix.get()[LG][LG] = 0;
	GapScoreMatrix.get()[LG][GL] = g_scoreGapOpen.get();
	GapScoreMatrix.get()[LG][GG] = t * g_scoreGapOpen.get();

	GapScoreMatrix.get()[GL][LL] = 0;
	GapScoreMatrix.get()[GL][LG] = g_scoreGapOpen.get();
	GapScoreMatrix.get()[GL][GL] = 0;
	GapScoreMatrix.get()[GL][GG] = 0;

	GapScoreMatrix.get()[GG][LL] = 0;
	GapScoreMatrix.get()[GG][LG] = t * g_scoreGapOpen.get();
	GapScoreMatrix.get()[GG][GL] = 0;
	GapScoreMatrix.get()[GG][GG] = 0;

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < i; ++j)
			if (GapScoreMatrix.get()[i][j] != GapScoreMatrix.get()[j][i])
				Quit("GapScoreMatrix.get() not symmetrical");
}

static SCORE SPColBrute(const MSA &msa, unsigned uColIndex)
{
	SCORE Sum = 0;
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
	{
		const WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		unsigned uLetter1 = msa.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter1 >= 20)
			continue;
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
		{
			const WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			unsigned uLetter2 = msa.GetLetterEx(uSeqIndex2, uColIndex);
			if (uLetter2 >= 20)
				continue;
			SCORE t = w1 * w2 * (*g_ptrScoreMatrix.get())[uLetter1][uLetter2];
#if	TRACE
			Log("Check %c %c w1=%.3g w2=%.3g Mx=%.3g t=%.3g\n",
				LetterToCharAmino(uLetter1),
				LetterToCharAmino(uLetter2),
				w1,
				w2,
				(*g_ptrScoreMatrix.get())[uLetter1][uLetter2],
				t);
#endif
			Sum += t;
		}
	}
	return Sum;
}

static SCORE SPGapFreqs(const FCOUNT Freqs[])
{
#if TRACE
	Log("Freqs=");
	for (unsigned i = 0; i < 4; ++i)
		if (Freqs[i] != 0)
			Log(" %s=%.3g", GapTypeToStr(i), Freqs[i]);
	Log("\n");
#endif

	SCORE TotalOffDiag = 0;
	SCORE TotalDiag = 0;
	for (unsigned i = 0; i < 4; ++i)
	{
		const FCOUNT fi = Freqs[i];
		if (0 == fi)
			continue;
		const float *Row = GapScoreMatrix.get()[i];
		SCORE diagt = fi * fi * Row[i];
		TotalDiag += diagt;
#if	TRACE
		Log("SPFGaps %s %s + Mx=%.3g TotalDiag += %.3g\n",
			GapTypeToStr(i),
			GapTypeToStr(i),
			Row[i],
			diagt);
#endif
		SCORE Sum = 0;
		for (unsigned j = 0; j < i; ++j)
		{
			SCORE t = Freqs[j] * Row[j];
#if	TRACE
			if (Freqs[j] != 0)
				Log("SPFGaps %s %s + Mx=%.3g Sum += %.3g\n",
					GapTypeToStr(i),
					GapTypeToStr(j),
					Row[j],
					fi * t);
#endif
			Sum += t;
		}
		TotalOffDiag += fi * Sum;
	}
#if TRACE
	Log("SPFGap TotalOffDiag=%.3g + TotalDiag=%.3g = %.3g\n",
		TotalOffDiag, TotalDiag, TotalOffDiag + TotalDiag);
#endif
	return TotalOffDiag * 2 + TotalDiag;
}

static SCORE SPFreqs(const FCOUNT Freqs[])
{
#if TRACE
	Log("Freqs=");
	for (unsigned i = 0; i < 20; ++i)
		if (Freqs[i] != 0)
			Log(" %c=%.3g", LetterToCharAmino(i), Freqs[i]);
	Log("\n");
#endif

	SCORE TotalOffDiag = 0;
	SCORE TotalDiag = 0;
	for (unsigned i = 0; i < 20; ++i)
	{
		const FCOUNT fi = Freqs[i];
		if (0 == fi)
			continue;
		const float *Row = (*g_ptrScoreMatrix.get())[i];
		SCORE diagt = fi * fi * Row[i];
		TotalDiag += diagt;
#if	TRACE
		Log("SPF %c %c + Mx=%.3g TotalDiag += %.3g\n",
			LetterToCharAmino(i),
			LetterToCharAmino(i),
			Row[i],
			diagt);
#endif
		SCORE Sum = 0;
		for (unsigned j = 0; j < i; ++j)
		{
			SCORE t = Freqs[j] * Row[j];
#if	TRACE
			if (Freqs[j] != 0)
				Log("SPF %c %c + Mx=%.3g Sum += %.3g\n",
					LetterToCharAmino(i),
					LetterToCharAmino(j),
					Row[j],
					fi * t);
#endif
			Sum += t;
		}
		TotalOffDiag += fi * Sum;
	}
#if TRACE
	Log("SPF TotalOffDiag=%.3g + TotalDiag=%.3g = %.3g\n",
		TotalOffDiag, TotalDiag, TotalOffDiag + TotalDiag);
#endif
	return TotalOffDiag * 2 + TotalDiag;
}

} // namespace muscle
