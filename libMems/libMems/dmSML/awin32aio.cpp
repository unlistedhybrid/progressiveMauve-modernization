#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/awin32aio.h"
#include "libMems/dmSML/util.h"

// Added for printf
#include <cstdio>

#ifdef USE_WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <cstring>
#include <cstdlib>

static VOID CALLBACK DummyCompletionRoutine( DWORD err, DWORD nbytes, LPOVERLAPPED lpo ) {
    printf( "completion routine!\n" );
}


int OpenWIN32( aFILE * file, const char *path, int mode ) {
    HANDLE result;
    DWORD access = mode == A_READ ? GENERIC_READ : GENERIC_WRITE;
    DWORD disposition = mode == A_READ ? OPEN_EXISTING : CREATE_ALWAYS;
    result = CreateFile( 
        path, 
        access, 
        FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
        nullptr,
        disposition,
        FILE_FLAG_OVERLAPPED,
        nullptr );
    if( result == INVALID_HANDLE_VALUE ) {
    	access = GetLastError();
    	printf( "Error opening %s, code %u\n", path, access );
        return( 0 );
    }
    file->w32handle = result;
    return( 1 );
}


int CloseWIN32( aFILE * file ) {
    return( CloseHandle( file->w32handle ) );
}


int WriteWIN32( aFILE * file, aIORec * rec ) {

    static offset_t total_bytes = 0;
    DWORD err;
    if( file->mode != A_WRITE ) {
        return( 0 );
    }

    rec->w32overlapped = static_cast<OVERLAPPED*>(malloc( sizeof( *(rec->w32overlapped) ) ));
    std::memset( rec->w32overlapped, 0, sizeof( *(rec->w32overlapped) ) );
	
	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		tmppos >>= 32;
		file->filep_high = (DWORD)tmppos; 
		tmppos = rec->pos;
		tmppos <<= 32;
		tmppos >>= 32;
		file->filep_low = (DWORD)tmppos;
	}

    rec->w32overlapped->OffsetHigh = file->filep_high;
    rec->w32overlapped->Offset = file->filep_low;

    total_bytes += rec->size * rec->count;
    if( WriteFileEx( 
        file->w32handle, 
        rec->buf, 
        (DWORD)(rec->size*rec->count), 
        rec->w32overlapped,
        DummyCompletionRoutine ) == 0 ) {
        err = GetLastError();
        printf( "error with WriteFileEx: %u\n", err );
        return( 0 );
    }
    return( 1 );
}


int ReadWIN32( aFILE * file, aIORec * rec ) {
    DWORD err;
    if( file->mode != A_READ ) {
        return( 0 );
    }
    rec->w32overlapped = static_cast<OVERLAPPED*>(malloc( sizeof( *(rec->w32overlapped) ) ));
    std::memset( rec->w32overlapped, 0, sizeof( *(rec->w32overlapped) ) );

	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		tmppos >>= 32;
		file->filep_high = (DWORD)tmppos;
		tmppos = rec->pos;
		tmppos <<= 32;
		tmppos >>= 32;
		file->filep_low = (DWORD)tmppos;
	}

    rec->w32overlapped->OffsetHigh = file->filep_high;
    rec->w32overlapped->Offset = file->filep_low;
    if( ReadFileEx( 
        file->w32handle, 
        rec->buf, 
        (DWORD)(rec->size*rec->count), 
        rec->w32overlapped,
        DummyCompletionRoutine ) == 0 ) {
        err = GetLastError();
        switch( err ) {
        case ERROR_HANDLE_EOF:
            printf( "readfileex says EOF -- we'll pretend it worked\n" );
            return( 1 );
        default:
            printf( "error with ReadFileEx -- Last Error: %u\n", GetLastError() );
            printf( "called:  ReadFileEx( %p, %p, %u, %p, %p )\n", 
                (void*)file->w32handle, 
                rec->buf, 
                (unsigned int)(rec->size*rec->count), 
                (void*)rec->w32overlapped,
                (void*)DummyCompletionRoutine );
            return( 0 );
        }
    }
    return( 1 );
}


int QueryLastCompleteWIN32( aFILE * file ) {
    DWORD result;
    if( file->queuetail && file->queuetail->w32overlapped ) {
        result = WaitForSingleObject( file->w32handle, 0 );
        if( result != WAIT_TIMEOUT ) {
            return( 1 );
        } else {
            return( 0 );
        }
    } else {
        return( 0 );
    }
}

#endif /* USE_WIN32 */
