#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDefs.h"

namespace libGenome {
    // some compilers don't have abs() for 64 bit ints
    int64 abs( int64 a ) noexcept {
        return a < 0 ? -a : a;
    }
}
