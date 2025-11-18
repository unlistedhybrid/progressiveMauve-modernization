#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctime>
#include <cstddef>
#include "libMems/dmSML/buffer.h"
#include <cstring>

static int BigRandom() {
    static char firsttime = 1;
    int i, result;
    if( firsttime ) {
        firsttime = 0;
        srand( 0 );
    }
    
    result = 0;
    for( i = 0; i < static_cast<int>(sizeof( result )); i++ ) {
        result <<= sizeof( result );
        result ^= rand();
    }
    return( result < 0 ? (-result < 0 ? 0 : -result) : result );
}


int MakeWorkingSet( working_set_t * ws, offset_t goalsize, offset_t minrecs, offset_t maxrecs ) {
	while( 1 ){

	    offset_t cursize = 0;
	    offset_t overhead = sizeof( ws->bufs[0] );
	    offset_t minsize = overhead + minrecs * sizeof( record_t );
	    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
	    offset_t nbufs = 0;
	    offset_t maxbufs = 256;
	    offset_t *buflist = static_cast<offset_t*>(malloc( sizeof( *buflist ) * maxbufs ));
	    
	    record_t *recordptr;
	    offset_t i;
	    if( goalsize < minsize || maxrecs < minrecs || !buflist ) {
	    	if( buflist )
		        free( buflist );
	        return( 0 );
	    }

	    while( goalsize - cursize >= maxsize ) {
	        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
	        if( nbufs == maxbufs ) {
	            maxbufs *= 2;
	            buflist = static_cast<offset_t*>(realloc( buflist, sizeof( *buflist ) * maxbufs ));
	        }
	        buflist[nbufs++] = randrecs;
	        cursize += overhead + randrecs * sizeof( record_t );
	    }
		printf( "allocating %llu bytes for working set (%llu bufs)\n", cursize, nbufs );

		ws->bufs = static_cast<buffer_t*>(malloc( cursize ));
		if( !ws->bufs ){
			goalsize /= 2;
			continue;
		}

		ws->size = cursize;
		ws->nbufs = nbufs;
		std::memset( ws->bufs, 0, cursize );

		recordptr = reinterpret_cast<record_t*>( reinterpret_cast<std::ptrdiff_t>(ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
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


void ReorganizeWorkingSet( working_set_t * ws, offset_t minrecs, offset_t maxrecs ) {
    offset_t goalsize = ws->size;
    offset_t cursize = 0;
    offset_t overhead = sizeof( ws->bufs[0] );
    offset_t minsize = overhead + minrecs * sizeof( record_t );
    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
    offset_t nbufs = 0;
    offset_t maxbufs = 256;
    offset_t *buflist = static_cast<offset_t*>(malloc( sizeof( *buflist ) * maxbufs ));
    offset_t leftovers;
    record_t *recordptr;
    offset_t i;
    
    if( maxrecs < minrecs ) {
        free( buflist );
        return;
    }

    if( goalsize < minsize ) {
        minsize = goalsize;
        minrecs = (minsize-overhead) / sizeof( record_t );
    }
    
    while( goalsize - cursize >= maxsize ) {
        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
        if( nbufs == maxbufs ) {
            maxbufs *= 2;
            buflist = static_cast<offset_t*>(realloc( buflist, sizeof( *buflist ) * maxbufs ));
        }
        buflist[nbufs++] = randrecs;
        cursize += overhead + randrecs * sizeof( record_t );
    }
    
    if( goalsize - cursize > overhead ) {
        leftovers = (goalsize - cursize - overhead) / sizeof( record_t );
        if( leftovers ) {
            if( nbufs == maxbufs ) {
                maxbufs *= 2;
                buflist = static_cast<offset_t*>(realloc( buflist, sizeof( *buflist ) * maxbufs ));
            }
            buflist[nbufs++] = leftovers;
            cursize += overhead + leftovers * sizeof( record_t );
        }
    }

    ws->nbufs = nbufs;
    std::memset( ws->bufs, 0, cursize );
    recordptr = reinterpret_cast<record_t*>( reinterpret_cast<std::ptrdiff_t>(ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
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
            }
        }
    }
}


buffer_list_t * InitList( buffer_list_t * list ) {
    list->head = nullptr;
    list->nitems = 0;
    return( list );
}


void PushHead( buffer_list_t * list, buffer_t * item ) {
    if( list->head == nullptr ) {
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
    if( list->head == nullptr ) {
        return( nullptr );
    }
    list->head->next->last = list->head->last;
    list->head->last->next = list->head->next;
    ret = list->head;
    list->head = list->head->next;
    ret->next = ret->last = nullptr;
    list->nitems--;
    if( list->nitems == 0 ) {
        list->head = nullptr;
    }
    return( ret );
}

buffer_t * PopTail( buffer_list_t * list ) {
    if( list->head == nullptr ) {
        return( list->head );
    }
    list->head = list->head->last;
    return( PopHead( list ) );
}

buffer_t * RemoveItem( buffer_list_t * list, buffer_t * item ) {
    if( item == list->head ) {
        return( PopHead( list ) );
    }
    item->next->last = item->last;
    item->last->next = item->next;
    item->next = item->last = nullptr;
    list->nitems--;
    if( list->nitems == 0 ) {
        list->head = nullptr;
    }
    return( item );
}


int CompareKeys_qsort_wrapper( const void *r1, const void *r2 ) {
    return( CompareKeys( static_cast<const record_t*>(r1), static_cast<const record_t*>(r2) ) );
}


int CompareKeys( const record_t *r1, const record_t *r2 ) {
    return( COMPARE_KEYS( *r1, *r2 ) );
}


void UpdateDeviceIOExecuteState( working_set_t * ws, iodevice_t * dev ) {
    if( !dev->buf || dev->state == DEV_FREE || dev->buf->operation == OP_FINISHED ) {
        buffer_t *b;
        buffer_t *found_buf = nullptr;
        dev->state = DEV_FREE;
        dev->buf = nullptr;
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
