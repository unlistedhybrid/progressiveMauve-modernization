#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"

#ifndef HAVE_LLABS
// some compilers don't have abs() for 64 bit ints, so we define it here
int64 abs( int64 a ) noexcept {
	return a < 0 ? -a : a;
}
#endif
