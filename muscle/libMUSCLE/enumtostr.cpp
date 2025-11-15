#include "libMUSCLE/muscle.h"
#include "libMUSCLE/enums.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

namespace muscle {

static TLS<char[64]> szMsg;

static inline int icmp(const char* a, const char* b)
{
#if defined(_WIN32)
    return _stricmp(a, b);
#else
    return strcasecmp(a, b);
#endif
}

static const char* writeEnumValue(const char* enumName, int x)
{
    snprintf(szMsg.get(), 64, "%s_%d", enumName, x);
    return szMsg.get();
}

#define BEGIN_ENUM_TO_STR(t) const char* t##ToStr(t x) { switch (x) {
#define ENUM_CASE(t, v) case t##_##v: return #v;
#define ENUM_DEFAULT(t) default: return writeEnumValue(#t, static_cast<int>(x)); }}

#define BEGIN_STR_TO_ENUM(t) t StrTo##t(const char* Str) {
#define STR_ENUM_CASE(t, v) if (icmp(Str, #v) == 0) return t##_##v;
#define STR_ENUM_END(t) return t##_Undefined; }

#undef s
#undef c
#undef e

#define s(t)  BEGIN_ENUM_TO_STR(t)
#define c(t, x) ENUM_CASE(t, x)
#define e(t)  ENUM_DEFAULT(t)
#include "libMUSCLE/enums.h"

#undef s
#undef c
#undef e

#define s(t) BEGIN_STR_TO_ENUM(t)
#define c(t, x) STR_ENUM_CASE(t, x)
#define e(t) STR_ENUM_END(t)
#include "libMUSCLE/enums.h"

}
