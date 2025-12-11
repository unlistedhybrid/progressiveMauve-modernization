#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <time.h>
#include <stddef.h>
#include "libMems/dmSML/buffer.h"
#include <string.h>

// --- FIX: PORTABLE RANDOM NUMBER GENERATOR ---
static unsigned long pm_rng_next = 1;

static void pm_srand(unsigned int seed) {
    pm_rng_next = seed;
}

static int pm_rand(void) {
    pm_rng_next = pm_rng_next * 1103515245 + 12345;
    return((unsigned)(pm_rng_next/65536) % 32768);
}
// ---------------------------------------------

static int BigRandom() {
    static char firsttime = 1;
    int i, result;
    if( firsttime ) {
        firsttime = 0;
        pm_srand( 0 );
    }
    
    result = 0;
    for( i = 0; i < sizeof( result ); i++ ) {
        result <<= sizeof( result );
        result ^= pm_rand();
    }
    return( result < 0 ? (-result < 0 ? 0 : -result) : result );
}


// Working Set support.
// returns resulting size of the entire structure.
int MakeWorkingSet( working_set_t * ws, offset_t goalsize, offset_t minrecs, offset_t maxrecs ) {
// wrap the memory allocation loop with an outer loop
// that will attempt smaller working set sizes if large ones fail to allocate
	while( 1 ){

	    offset_t cursize = 0;
	    offset_t overhead = sizeof( ws->bufs[0] );
	    offset_t minsize = overhead + minrecs * sizeof( record_t );
	    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
	    offset_t nbufs = 0;      // number of real buffers pleged to the working set
	    offset_t maxbufs = 256;  // the max number of buffers we track (this grows if necessary)
	    offset_t *buflist = malloc( sizeof( *buflist ) * maxbufs ); // grows when necessary
	    
	    record_t *recordptr;
	    offset_t i;
	    // if we can't possibly do anything useful
	    if( goalsize < minsize || maxrecs < minrecs || !buflist ) {
	    	if( buflist )
		        free( buflist );
	        return( 0 );
	    }

	    // just start allocating buffers until we can't anymore
	    while( goalsize - cursize >= maxsize ) {
	        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
	        if( nbufs == maxbufs ) {
	            // resize the array
	            maxbufs *= 2;
	            buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
	        }
	        buflist[nbufs++] = randrecs;
	        // update the number of bytes we've currently decided to allocate.
	        cursize += overhead + randrecs * sizeof( record_t );
	    }
		// now we have nbufs buffers, and the number of records they should
		// store is in the buflist list.
		// allocate one big chunk of memory
		printf( "allocating %llu bytes for working set (%llu bufs)\n", cursize, nbufs );

		ws->bufs = malloc( cursize );
		// if it failed to allocate try a smaller size
		if( !ws->bufs ){
			goalsize /= 2;
			continue;
		}

		ws->size = cursize;
		ws->nbufs = nbufs;
		// clear it out
		memset( ws->bufs, 0, cursize );

		recordptr = (record_t *)( ((ptrdiff_t)ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
		for( i = 0; i < nbufs; i++ ) {
		    ws->bufs[i].totalrecs = buflist[i];
		    ws->bufs[i].recs = recordptr;
		    recordptr += ws->bufs[i].totalrecs;
		}

		free( buflist );
	    return( cursize );
	}
    return 0;
}


// Working Set support.
// Reorganize the working set with a different distribution of buffers.
void ReorganizeWorkingSet( working_set_t * ws, offset_t minrecs, offset_t maxrecs ) {
    offset_t goalsize = ws->size;
    offset_t cursize = 0;
    offset_t overhead = sizeof( ws->bufs[0] );
    offset_t minsize = overhead + minrecs * sizeof( record_t );
    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
    offset_t nbufs = 0;      // number of real buffers pledged to the working set
    offset_t maxbufs = 256;  // the max number of buffers we're tracking (this grows if necessary)
    offset_t *buflist = malloc( sizeof( *buflist ) * maxbufs ); // grows when necessary
    offset_t leftovers;
    record_t *recordptr;
    offset_t i;
    
    // if we can't possibly do anything useful
    if( maxrecs < minrecs ) {
        free( buflist );
        return;
    }

    if( goalsize < minsize ) {
        minsize = goalsize;
        minrecs = (minsize-overhead) / sizeof( record_t );
    }
    
    // just start allocating buffers until we can't anymore
    while( goalsize - cursize >= maxsize ) {
        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
        if( nbufs == maxbufs ) {
            // resize the array
            maxbufs *= 2;
            buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
        }
        buflist[nbufs++] = randrecs;
        // update the number of bytes we've currently decided to allocate.
        cursize += overhead + randrecs * sizeof( record_t );
    }
    
    // clean up the last bit
    if( goalsize - cursize > overhead ) {
        leftovers = (goalsize - cursize - overhead) / sizeof( record_t );
        if( leftovers ) {
            if( nbufs == maxbufs ) {
                // resize the array
                maxbufs *= 2;
                buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
            }
            buflist[nbufs++] = leftovers;
            cursize += overhead + leftovers * sizeof( record_t );
        }
    }

    ws->nbufs = nbufs;
    // clear it out
    memset( ws->bufs, 0, cursize );

    recordptr = (record_t *)( ((ptrdiff_t)ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
    for( i = 0; i < nbufs; i++ ) {
        ws->bufs[i].totalrecs = buflist[i];
        ws->bufs[i].recs = recordptr;
        recordptr += ws->bufs[i].totalrecs;
    }
    
    free( buflist );
    return;
}

void UpdateWSIOFinishedState( working_set_t * ws ) {
    buffer_t *b;
    for( b = ws->bufs; b - ws->bufs < ws->nbufs; b++ ) {
        if( b->operation > OP_NONE ) {
            if( aOperationComplete( b->file, b->operation ) ) {
                b->operation = OP_FINISHED;
            } else {
            }
        }
    }
}

buffer_list_t * InitList( buffer_list_t * list ) {
    list->head = NULL;
    list->nitems = 0;
    return( list );
}


void PushHead( buffer_list_t * list, buffer_t * item ) {
    if( list->head == NULL ) {
        list->head = item;
        list->nitems = 1;
        list->head->next = list->head;
        list->head->last = list->head;
        return;
    }
    item->last = list->head->last;
    item->next = list->head;
    list->head->last->next = item;
    list->head->last = item;
    list->head = item;
    list->nitems++;
}

void PushTail( buffer_list_t * list, buffer_t * item ) {
    PushHead( list, item );
    list->head = list->head->last;
}

buffer_t * PopHead( buffer_list_t * list ) {
    buffer_t *ret;
    if( list->head == NULL ) {
        return( NULL );
    }
    list->head->next->last = list->head->last;
    list->head->last->next = list->head->next;
    ret = list->head;
    list->head = list->head->next;
    ret->next = ret->last = NULL;
    list->nitems--;
    if( list->nitems == 0 ) {
        list->head = NULL;
    }
    return( ret );
}

buffer_t * PopTail( buffer_list_t * list ) {
    if( list->head == NULL ) {
        return( list->head );
    }
    list->head = list->head->last;
    return( PopHead( list ) );
}

// returns second argument
buffer_t * RemoveItem( buffer_list_t * list, buffer_t * item ) {
    // FIX: Handle NULL cases safely
    if ( !list || !item ) {
        return NULL;
    }

    if( item == list->head ) {
        return( PopHead( list ) );
    }
    
    // Safe removal logic
    if ( item->next ) item->next->last = item->last;
    if ( item->last ) item->last->next = item->next;
    
    item->next = item->last = NULL;
    
    if ( list->nitems > 0 ) list->nitems--;
    if( list->nitems == 0 ) {
        list->head = NULL;
    }
    return( item );
}


int CompareKeys_qsort_wrapper( const void *r1, const void *r2 ) {
    return( CompareKeys( (record_t *)r1, (record_t *)r2 ) );
}

int CompareKeys( const record_t *r1, const record_t *r2 ) {
    // FIX: TIE-BREAKER FOR STABLE SORT
    // If the primary keys are equal, qsort is not guaranteed to preserve order.
    // Platforms differ in qsort implementation (e.g. Mac vs Linux).
    // We break ties by comparing the memory content of the records.
    int ret = COMPARE_KEYS( *r1, *r2 );
    if ( ret == 0 ) {
        return memcmp(r1, r2, sizeof(record_t));
    }
    return ret;
}

void UpdateDeviceIOExecuteState( working_set_t * ws, iodevice_t * dev ) {
    if( !dev->buf || dev->state == DEV_FREE || dev->buf->operation == OP_FINISHED ) {
        buffer_t *b;
        buffer_t *found_buf = NULL;
        dev->state = DEV_FREE;
        dev->buf = NULL;

        for( b = ws->bufs; b - ws->bufs < ws->nbufs; b++ ) {
            if( b->operation == OP_PENDING && b->device == dev ) {
                if( !found_buf ) {
                    found_buf = b;
                } else if( (b->file == found_buf->file) && 
                    (b->fileop < found_buf->fileop) ) {
                    found_buf = b;
                }
            }
        }
        
        if( found_buf ) {
            dev->buf = found_buf;
            found_buf->operation = found_buf->file->mode == A_READ 
                ? aRead( found_buf->recs, 1, found_buf->io_size, found_buf->file, found_buf->io_pos )
                : aWrite( found_buf->recs, 1, found_buf->io_size, found_buf->file, found_buf->io_pos );
            dev->state = DEV_BUSY;
        }
    }
}

void ReadBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev ) {
	buffer->io_size = num_recs * sizeof( record_t );
    buffer->numrecs = num_recs;
    buffer->device = dev;
    buffer->fileop = buffer->file->op++;
    buffer->io_pos = CURRENT_POS;
    if( buffer->operation != OP_NONE ) {
        printf( "weird!\n" );
    } else {
        buffer->operation = OP_PENDING;
    }
}

void WriteBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev ) {
    ReadBuffer( buffer, num_recs, dev );
}