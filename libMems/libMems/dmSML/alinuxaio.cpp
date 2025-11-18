#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/alinuxaio.h"
#ifdef USE_LINUX_AIO

#include <libaio.h>

#include "libMems/dmSML/asyncio.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

io_context_t ctx_id = nullptr;

#ifndef __u64
typedef unsigned long long __u64;
#endif

__u64 current_id = 0;

unsigned event_max = 10000;

typedef struct completion_id_s {
	__u64 data;
	struct completion_id_s* next;
	struct completion_id_s* last;
} completion_id_t;

typedef struct completion_id_list_s { 
    int nitems;
    completion_id_t * head;
} completion_id_list_t;

completion_id_list_t * InitListComp( completion_id_list_t * list );
void PushHeadComp( completion_id_list_t * list, completion_id_t * item );
void PushTailComp( completion_id_list_t * list, completion_id_t * item );
completion_id_t * PopHeadComp( completion_id_list_t * list );
completion_id_t * PopTailComp( completion_id_list_t * list );
completion_id_t * RemoveItemComp( completion_id_list_t * list, completion_id_t * item );


completion_id_list_t * InitListComp( completion_id_list_t * list ) {
    list->head = nullptr;
    list->nitems = 0;
    return( list );
}


void PushHeadComp( completion_id_list_t * list, completion_id_t * item ) {
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

void PushTailComp( completion_id_list_t * list, completion_id_t * item ) {
    PushHeadComp( list, item );
    list->head = list->head->last;
}

completion_id_t * PopHeadComp( completion_id_list_t * list ) {
    completion_id_t *ret;
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

completion_id_t * PopTailComp( completion_id_list_t * list ) {
    if( list->head == nullptr ) {
        return( list->head );
    }
    list->head = list->head->last;
    return( PopHeadComp( list ) );
}

completion_id_t * RemoveItemComp( completion_id_list_t * list, completion_id_t * item ) {
    if( item == list->head ) {
        return( PopHeadComp( list ) );
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


completion_id_list_t *completion_list = nullptr;

int OpenLinux( aFILE * file, const char *path, int mode ){
	long error;
	if( ctx_id == nullptr ){
		error = io_queue_init( event_max, &ctx_id );
		if( error != 0 )
			perror( "io_setup" );
	}
	if( completion_list == nullptr ){
		completion_list = static_cast<completion_id_list_t*>(malloc( sizeof( completion_id_list_t ) ));
		completion_list = InitListComp( completion_list );
	}
		
	if(mode == A_READ){
		file->file_descriptor = open(path, O_LARGEFILE | O_RDONLY, S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP );
	}else{
		file->file_descriptor = open(path, O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE,  S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP);
	}
	if(file->file_descriptor < 0){
		
		perror(path);
	}
	return file->file_descriptor >= 0;
}

int CloseLinux( aFILE * file ){	
	return close( file->file_descriptor ) == 0;
}

void CleanupLinux(){
	free( completion_list );
	completion_list = nullptr;
	ctx_id = nullptr;
}

int FillAIOStruct( aFILE * file, aIORec * rec ){
	rec->aio_cb = static_cast<iocb_t*>( malloc( sizeof(iocb_t)));
	if(rec->aio_cb == nullptr)
		return 0;

	std::memset(rec->aio_cb, 0, sizeof(iocb_t));
	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		tmppos >>= 32;
		file->filep_high = tmppos;
		tmppos = rec->pos;
		tmppos <<= 32;
		tmppos >>= 32;
		file->filep_low = tmppos;
	}

	rec->aio_cb->aio_fildes = file->file_descriptor;
	rec->aio_cb->u.c.offset = file->filep_high;
	rec->aio_cb->u.c.offset <<= 32;
	rec->aio_cb->u.c.offset |= file->filep_low;
	rec->aio_cb->u.c.buf = rec->buf;
	rec->aio_cb->u.c.nbytes = rec->size * rec->count;
	
	return 1;
}

int WriteLinux( aFILE * file, aIORec * rec ){
        int req_error;
	struct iocb *request_array[] = { rec->aio_cb };
	if( FillAIOStruct( file, rec ) ){
		rec->aio_cb->aio_lio_opcode = IO_CMD_PWRITE;
		req_error = io_submit( ctx_id, 1, &rec->aio_cb );
		if(req_error != 1){
			printf("write_submit: io_submit res=%d [%s]\n", req_error, strerror(-req_error));
            printf( "aiocb->aio_fildes = %d\n", rec->aio_cb->aio_fildes );
            printf( "aiocb->u.c.offset = %llu\n", rec->aio_cb->u.c.offset );
            printf( "aiocb->u.c.buf = %lx\n", rec->aio_cb->u.c.buf );
            printf( "aiocb->u.c.nbytes = %llu\n", rec->aio_cb->u.c.nbytes );
            printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
		}
		return req_error == 1;
	}
	return 0;
}

int ReadLinux( aFILE * file, aIORec * rec ){
	int req_error;
	struct iocb *request_array[] = { rec->aio_cb };
	if( FillAIOStruct( file, rec ) ){
		rec->aio_cb->aio_lio_opcode = IO_CMD_PREAD;
		req_error = io_submit( ctx_id, 1, &rec->aio_cb );
        if(req_error != 1){
			printf("read_submit: io_submit res=%d [%s]\n", req_error, strerror(-req_error));
                printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
        }
		return req_error == 1;
	}
	return 0;
}


int QueryLastCompleteLinux( aFILE * file ){
	int rval;
	int compI;
	completion_id_t *comp;
	struct io_event ioe;
	struct timespec zero_wait;

	zero_wait.tv_sec = 0;
	zero_wait.tv_nsec = 10000000;
	
	rval = io_getevents( ctx_id, 0, 1, &ioe, &zero_wait );
	if( rval == 1 ){
		completion_id_t *completion = static_cast<completion_id_t*>(malloc( sizeof(completion_id_t) ));
		completion->data = ioe.data;
		PushTailComp( completion_list, completion );
	}
	comp = completion_list->head;
	for( compI = 0; compI < completion_list->nitems; compI++ ){
		if( comp->data == ioe.data )
			break;
	}
	if( compI != completion_list->nitems ){
		RemoveItemComp( completion_list, comp );
		return 1;
	}
	return 0;
}

#endif
