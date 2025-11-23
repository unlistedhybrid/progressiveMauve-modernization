#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"

namespace genome {

// some compilers don't have abs() for 64 bit ints
int64 abs_int64(int64 a) {
	return a < 0 ? -a : a;
}

} // namespace genome
