#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"

// Guard this definition so it only compiles on systems that actually need it
// (Matches the logic in gnDefs.h)
#if (defined(__GNUG__) && ( __GNUC__ <= 2 )) || defined(__INTEL_COMPILER) || (defined _MSC_VER && defined __cplusplus)

// some compilers don't have abs() for 64 bit ints
int64 abs( int64 a ){
	return a < 0 ? -a : a;
}

#endif