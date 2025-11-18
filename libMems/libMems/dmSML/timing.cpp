#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef WIN32
#include <sys/time.h>
#include <unistd.h>
#endif

#include <cstdlib>
#include <cstdio>


#include "libMems/dmSML/util.h"
#include "libMems/dmSML/timing.h"


struct dmtimer_s {
#ifdef WIN32    
    unsigned int last;
#else
    struct timeval tv;
#endif
};



typedef int Int;
typedef unsigned int UInt;
typedef double Float64;

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
#define NUM_FREQ_BITS   (14)
static Int ShiftAmt;
static Int TicksPerSecond;
static Int LastReadValue;
static Int BaseTime;
#endif /* WIN32 */


dmtimer_t * StartTimer() {
#ifdef WIN32
    dmtimer_t * t = static_cast<dmtimer_t*>(malloc( sizeof( *t ) ));
    t->last = timeGetTime();
    return( t );
#else
    dmtimer_t * t = static_cast<dmtimer_t*>(malloc( sizeof( *t ) ));
    gettimeofday( &(t->tv), nullptr );
    return( t );
#endif /* WIN32 */
}



unsigned int ReadTimer( dmtimer_t * t ) {
#ifdef WIN32
    unsigned int cur = timeGetTime();
    return( cur - t->last );
#else
    struct timeval current;
    struct timezone dummy;
    unsigned int begintime, endtime;
    gettimeofday( &current, &dummy );
    begintime = 1000 * t->tv.tv_sec + (t->tv.tv_usec/1000);
    endtime = 1000 * current.tv_sec + (current.tv_usec/1000);
    return( endtime - begintime );
#endif
}



void StopTimer( dmtimer_t * t ) {
    free( t );
}



#ifdef WIN32
static void InitTimeWIN32() {
    timeBeginPeriod( 1 );
}
#endif /* WIN32 */


void InitTime() {
#ifdef WIN32    
    InitTimeWIN32();
#endif
}
