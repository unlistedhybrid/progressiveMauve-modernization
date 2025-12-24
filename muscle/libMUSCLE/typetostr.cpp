#include "libMUSCLE/muscle.h"
#include <stdio.h>

namespace muscle {

const char *SecsToStr(unsigned long Secs)
	{
	static TLS<char[32]> Str; // Increased buffer size for safety
	long hh, mm, ss;

	hh = Secs/(60*60);
	mm = (Secs/60)%60;
	ss = Secs%60;

	snprintf(Str.get(), 32, "%02ld:%02ld:%02ld", hh, mm, ss);
	return Str.get();
	}

const char *BoolToStr(bool b)
	{
	return b ? "True" : "False";
	}

const char *ScoreToStr(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return "       *";
// Hack to use "circular" buffer so when called multiple
// times in a printf-like argument list it works OK.
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static TLS<char[iBufferCount*iBufferLength]> szStr;
	static TLS<int> iBufferIndex(0);
	iBufferIndex.get() = (iBufferIndex.get() + 1)%iBufferCount;
	char *pStr = szStr.get() + iBufferIndex.get()*iBufferLength;
	snprintf(pStr, iBufferLength, "%8g", Score);
	return pStr;
	}

// Left-justified version of ScoreToStr
const char *ScoreToStrL(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return "*";
// Hack to use "circular" buffer so when called multiple
// times in a printf-like argument list it works OK.
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static TLS<char[iBufferCount*iBufferLength]> szStr;
	static TLS<int> iBufferIndex(0);
	iBufferIndex.get() = (iBufferIndex.get() + 1)%iBufferCount;
	char *pStr = szStr.get() + iBufferIndex.get()*iBufferLength;
	snprintf(pStr, iBufferLength, "%.3g", Score);
	return pStr;
	}

const char *WeightToStr(WEIGHT w)
	{
	return ScoreToStr(w);
	}
} 
