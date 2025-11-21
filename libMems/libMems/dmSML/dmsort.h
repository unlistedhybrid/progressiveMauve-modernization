#ifndef __DMSORT_H__
#define __DMSORT_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "libMems/dmSML/util.h"
#include "libMems/dmSML/timing.h"
#include "libMems/dmSML/asyncio.h"
#include "libMems/dmSML/buffer.h"
#include "libMems/dmSML/sorting.h"
#include "libMems/dmSML/sml.h"

//#define ASCII_KEYBYTES
//#define NNNNN_KEYBYTES
//#define NO_SORT_PERF_TEST
//#define NO_WRITE_PERF_TEST
//#define NO_BINNING_PERF_TEST
//#define NO_BIN_WRITE_PERF_TEST
//#define NO_RESTRUCTURE_PERF_TEST

#ifndef NELEMS
#define NELEMS(x) \
    ( sizeof((x)) / sizeof((x)[0]) )
#endif

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define MINRECS     (1311)
#define MAXRECS     (1311)

typedef struct device_s {
    const char      *devname;
    const char      *path;
    iodevice_t      dev;
} device_t;

#define BIN_SPECIAL     (-10000)

typedef struct bin_s {
    aFILE               *file;
    int                 dev;
    offset_t            nrecs;
    buffer_list_t       bufs;
    char*				fname;
} bin_t;

typedef struct seqbuf_s {
	aFILE				*file;
	int					dev;
	offset_t			bufpos;
	uint64				seq_pos;
	buffer_list_t		bufs;
} seqbuf_t;

enum dm_errors {
	SUCCESS,
	TOO_FEW_BINS,
	TOO_MANY_BINS,
	INPUT_NOT_OPENED,
	INVALID_WS_SIZE,
	SEQUENCE_TOO_SHORT,
	OUTPUT_NOT_OPENED,
	INVALID_NUMRECS,
	NO_FREE_BUFFERS,
	BIN_NOT_OPENED,
};


void FinishBinning();
offset_t CalculateDataReadSize( buffer_t* b );

#define ALPHA_BITS 2

void RestructureReadSMLBins( void );
int InitdmSML( long working_mb, long buffer_size, const char* input_filename, const char* output_filename, const char* const* scratch_paths, uint64 seed );
void DisplayStatusHeader( void );
void DisplayStatus( void );
void UpdateIOState( void );
void EnsureAllOperationsComplete( void );
void BinningPhase( void );
void SortReading( void );

#ifdef USE_QSORT_ONLY

int comp_keys( record_t a, record_t b );
void QBrute( record_t a[], int lo, int hi );
void QSort( record_t a[], int lo0, int hi0 );
void RecSort( record_t a[], int nelems );
int SortBuffer( buffer_t * buf );
void SortSorting( void );

#elif defined NO_SORT_PERF_TEST

void SortSorting( void );

#else 

sort_buf_t* CurrentSortBuf;
buffer_t* SortScratchBuffer;
void SortSorting( void );

#endif

void RestructureSMLBinsForWrite( void );
int CalculateSortWriteSize( int sortI );
void SortWriting( void );
void SortHandleCompletions( void );
void SortUpdateIOState();
void SortingEnsureAllOperationsComplete();
void SortingPhase( void );
int dmsort( void );

#ifdef __cplusplus
extern "C" {
#endif

int dmSML( const char* input_file, const char* output_file, const char* const* scratch_paths, uint64 seed );

#ifdef __cplusplus
}
#endif

#endif
