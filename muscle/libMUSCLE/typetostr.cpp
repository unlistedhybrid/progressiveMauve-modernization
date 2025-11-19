#include "libMUSCLE/muscle.h"
#include <stdio.h>

namespace muscle {

const char *SecsToStr(unsigned long Secs)
	{
	// Fixed: Increased buffer size from 16 to 32 to prevent truncation warnings.
	// 16 bytes causes -Wformat-truncation on 64-bit systems because very large 
	// 'hh' values would overflow the formatted string.
	static TLS<char[32]> Str;
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
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static TLS<char[iBufferCount*iBufferLength]> szStr;
	static TLS<int> iBufferIndex(0);
	iBufferIndex.get() = (iBufferIndex.get() + 1)%iBufferCount;
	char *pStr = szStr.get() + iBufferIndex.get()*iBufferLength;
	snprintf(pStr, 16, "%8g", Score);
	return pStr;
	}

const char *ScoreToStrL(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return "*";
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static TLS<char[iBufferCount*iBufferLength]> szStr;
	static TLS<int> iBufferIndex(0);
	iBufferIndex.get() = (iBufferIndex.get() + 1)%iBufferCount;
	char *pStr = szStr.get() + iBufferIndex.get()*iBufferLength;
	snprintf(pStr, 16, "%.3g", Score);
	return pStr;
	}

const char *WeightToStr(WEIGHT w)
	{
	return ScoreToStr(w);
	}
}
