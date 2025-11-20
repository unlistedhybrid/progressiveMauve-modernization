#include "libMUSCLE/muscle.h"
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cerrno>
#include "libMUSCLE/threadstorage.h"

#if defined(_WIN32)
#include <windows.h>
#include <share.h>
#define strdup _strdup
#endif

namespace muscle {

#ifndef MAX_PATH
#define MAX_PATH 260
#endif

static TLS<char[MAX_PATH]> g_strListFileName;
static TLS<bool> g_bListFileAppend(false);

TLS<int> g_argc;
TLS<char **> g_argv;
TLS<MSA*> ptrBestMSA;

static TLS<SEQWEIGHT> g_SeqWeight(SEQWEIGHT_Undefined);

void SetSeqWeightMethod(SEQWEIGHT Method)
{
    g_SeqWeight.get() = Method;
}

SEQWEIGHT GetSeqWeightMethod()
{
    return g_SeqWeight.get();
}

void SetListFileName(const char *ptrListFileName, bool bAppend)
{
    size_t len = std::strlen(ptrListFileName);
    if (len >= MAX_PATH) len = MAX_PATH - 1;
    std::strncpy(g_strListFileName.get(), ptrListFileName, len);
    g_strListFileName.get()[len] = '\0';
    g_bListFileAppend.get() = bAppend;
}

void Log(const char szFormat[], ...)
{
    if (g_strListFileName.get()[0] == 0)
        return;

    static TLS<FILE *> f(nullptr);
    const char *mode = g_bListFileAppend.get() ? "a" : "w";
#if defined(_WIN32)
    if (f.get() == nullptr)
        f.get() = _fsopen(g_strListFileName.get(), mode, _SH_DENYNO);
#else
    if (f.get() == nullptr)
        f.get() = fopen(g_strListFileName.get(), mode);
#endif
    if (f.get() == nullptr) {
        perror(g_strListFileName.get());
        exit(EXIT_NotStarted);
    }

    char szStr[4096];
    va_list ArgList;
    va_start(ArgList, szFormat);
    vsnprintf(szStr, sizeof(szStr), szFormat, ArgList);
    va_end(ArgList);
    fprintf(f.get(), "%s", szStr);
    fflush(f.get());
}

const char *GetTimeAsStr()
{
    static TLS<char[32]> szStr;
    time_t t;
    time(&t);
    struct tm *ptmCurrentTime = localtime(&t);
    std::strftime(szStr.get(), 32, "%a %b %d %H:%M:%S %Y", ptmCurrentTime);
    return szStr.get();
}

void Quit(const char szFormat[], ...)
{
    va_list ArgList;
    char szStr[4096];
    va_start(ArgList, szFormat);
    vsnprintf(szStr, sizeof(szStr), szFormat, ArgList);
    va_end(ArgList);

    fprintf(stderr, "\n*** ERROR ***  %s\n", szStr);

    Log("\n*** FATAL ERROR ***  ");
    Log("%s\n", szStr);
    Log("Stopped %s\n", GetTimeAsStr());

#if defined(_WIN32)
    if (IsDebuggerPresent()) {
        int iBtn = MessageBoxA(NULL, szStr, "muscle", MB_ICONERROR | MB_OKCANCEL);
        if (IDCANCEL == iBtn)
            Break();
    }
#endif
    exit(EXIT_FatalError);
}

void Warning(const char szFormat[], ...)
{
    va_list ArgList;
    char szStr[4096];
    va_start(ArgList, szFormat);
    vsnprintf(szStr, sizeof(szStr), szFormat, ArgList);
    va_end(ArgList);

    fprintf(stderr, "\n*** WARNING *** %s\n", szStr);
    Log("\n*** WARNING ***  %s\n", szStr);
}

void TrimBlanks(char szStr[])
{
    TrimLeadingBlanks(szStr);
    TrimTrailingBlanks(szStr);
}

void TrimLeadingBlanks(char szStr[])
{
    size_t n = std::strlen(szStr);
    while (szStr[0] == ' ' && n > 0) {
        std::memmove(szStr, szStr + 1, n);
        szStr[--n] = 0;
    }
}

void TrimTrailingBlanks(char szStr[])
{
    size_t n = std::strlen(szStr);
    while (n > 0 && szStr[n - 1] == ' ') {
        szStr[--n] = 0;
    }
}

bool Verbose()
{
    return true;
}

SCORE StrToScore(const char *pszStr)
{
    return (SCORE) std::atof(pszStr);
}

void StripWhitespace(char szStr[])
{
    unsigned uOutPos = 0;
    unsigned uInPos = 0;
    char c;
    while ((c = szStr[uInPos++]))
        if (c != ' ' && c != '\t' && c != '\n' && c != '\r')
            szStr[uOutPos++] = c;
    szStr[uOutPos] = 0;
}

void StripGaps(char szStr[])
{
    unsigned uOutPos = 0;
    unsigned uInPos = 0;
    char c;
    while ((c = szStr[uInPos++]))
        if (c != '-')
            szStr[uOutPos++] = c;
    szStr[uOutPos] = 0;
}

bool IsValidSignedInteger(const char *Str)
{
    size_t len = std::strlen(Str);
    if (len == 0) return false;
    if (*Str == '+' || *Str == '-') ++Str;
    while (*Str)
        if (!isdigit(*Str++)) return false;
    return true;
}

bool IsValidInteger(const char *Str)
{
    size_t len = std::strlen(Str);
    if (len == 0) return false;
    while (*Str)
        if (!isdigit(*Str++)) return false;
    return true;
}

bool isidentf(char c)
{
    return isalpha((unsigned char)(c)) || c == '_';
}

bool isident(char c)
{
    return isalpha((unsigned char)(c)) ||
           isdigit((unsigned char)(c)) ||
           c == '_';
}

bool IsValidIdentifier(const char *Str)
{
    if (!isidentf(Str[0]))
        return false;
    const char *p = Str + 1;
    while (*p)
        if (!isident(*p++)) return false;
    return true;
}

void SetLogFile()
{
    const char *strFileName = ValueOpt("loga");
    if (strFileName != nullptr)
        g_bListFileAppend.get() = true;
    else
        strFileName = ValueOpt("log");
    if (strFileName == nullptr)
        return;
    size_t len = std::strlen(strFileName);
    if (len >= MAX_PATH) len = MAX_PATH - 1;
    std::strncpy(g_strListFileName.get(), strFileName, len);
    g_strListFileName.get()[len] = '\0';
}

void NameFromPath(const char szPath[], char szName[], unsigned uBytes)
{
    if (uBytes == 0) return;
    const char *pstrLastSlash = std::strrchr(szPath, '/');
    const char *pstrLastBackslash = std::strrchr(szPath, '\\');
    const char *pstrLastDot = std::strrchr(szPath, '.');
    const char *pstrLastSep = (pstrLastSlash > pstrLastBackslash) ? pstrLastSlash : pstrLastBackslash;
    const char *pstrBegin = pstrLastSep ? pstrLastSep + 1 : szPath;
    const char *pstrEnd = pstrLastDot ? pstrLastDot - 1 : szPath + std::strlen(szPath);
    unsigned uNameLength = (unsigned)(pstrEnd - pstrBegin + 1);
    if (uNameLength > uBytes - 1)
        uNameLength = uBytes - 1;
    std::memcpy(szName, pstrBegin, uNameLength);
    szName[uNameLength] = 0;
}

char *strsave(const char *s)
{
    char *ptrCopy = strdup(s);
    if (ptrCopy == nullptr)
        Quit("Out of memory");
    return ptrCopy;
}

bool IsValidFloatChar(char c)
{
    return isdigit((unsigned char)(c)) || c == '.' || c == 'e' || c == 'E' ||
           c == 'd' || c == 'D' || c == '+' || c == '-';
}

void Call_MY_ASSERT(const char *file, int line, bool b, const char *msg)
{
    if (b)
        return;
    Quit("%s(%d): MY_ASSERT(%s)", file, line, msg);
}

static size_t g_MemTotal = 0;

void MemPlus(size_t Bytes, char *Where)
{
    g_MemTotal += Bytes;
    Log("+%10u  %6u  %6u  %s\n",
        (unsigned)Bytes,
        (unsigned)GetMemUseMB(),
        (unsigned)(g_MemTotal / 1000000),
        Where);
}

void MemMinus(size_t Bytes, char *Where)
{
    g_MemTotal -= Bytes;
    Log("-%10u  %6u  %6u  %s\n",
        (unsigned)Bytes,
        (unsigned)GetMemUseMB(),
        (unsigned)(g_MemTotal / 1000000),
        Where);
}

}
