#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/aPOSIXaio.h"
#ifdef USE_POSIX_AIO

#include "libMems/dmSML/asyncio.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

int OpenPAIO( aFILE * file, const char *path, int mode ){
	int flags = 0;
#ifdef O_LARGEFILE
	flags |= O_LARGEFILE;
#endif
	if(mode == A_READ){
		file->file_descriptor = open(path, flags | O_RDONLY, S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP );
	}else{
		file->file_descriptor = open(path, flags | O_RDWR | O_CREAT | O_TRUNC,  S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP);
	}
	if(file->file_descriptor < 0){
		
		perror(path);
	}
	return file->file_descriptor >= 0;
}

int ClosePAIO( aFILE * file ){
	return close( file->file_descriptor ) == 0;
}

int FillAIOStruct( aFILE * file, aIORec * rec ){
	rec->aio_cb = static_cast<aiocb_t*>( malloc( sizeof(aiocb_t)));
	if(rec->aio_cb == nullptr)
		return 0;

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
	rec->aio_cb->aio_offset = file->filep_high;
	rec->aio_cb->aio_offset <<= 32;
	rec->aio_cb->aio_offset |= file->filep_low;
	rec->aio_cb->aio_buf = rec->buf;
	rec->aio_cb->aio_nbytes = rec->size * rec->count;
	rec->aio_cb->aio_reqprio = 0;
	std::memset(&(rec->aio_cb->aio_sigevent), 0, sizeof(struct sigevent) );
	return 1;
}

int WritePAIO( aFILE * file, aIORec * rec ){
        int req_error;
	if( FillAIOStruct( file, rec ) ){
		rec->aio_cb->aio_lio_opcode = LIO_WRITE;
		req_error = aio_write(rec->aio_cb);
		if(req_error == -1){
			perror("write");
            printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
		}
		return req_error == 0;
	}
	return 0;
}

int ReadPAIO( aFILE * file, aIORec * rec ){
	int req_error;
	if( FillAIOStruct( file, rec ) ){
		rec->aio_cb->aio_lio_opcode = LIO_READ;
		req_error = aio_read(rec->aio_cb);
        if(req_error == -1){
                perror("write");
                printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
        }
		return req_error == 0;
	}
	return 0;
}

int QueryLastCompletePAIO( aFILE * file ){
	int rval;
	struct aiocb *request_array[] = { file->queuetail->aio_cb };
	struct timespec zero_wait;

	zero_wait.tv_sec = 0;
	zero_wait.tv_nsec = 0;
	
	rval = aio_suspend(request_array, 1, &zero_wait);
	if(rval == 0){
		return 1;
	}else if(rval == -1)
		;
	return 0;
}

#endif /* USE_POSIX_AIO */
