#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/awin32aio.h"
#include "libMems/dmSML/util.h"
#include <cstdio>

#ifdef USE_WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <cstring>

// No completion routine needed for synchronous I/O

int OpenWIN32( aFILE * file, const char *path, int mode ) {
    HANDLE result;
    DWORD access = mode == A_READ ? GENERIC_READ : GENERIC_WRITE;
    DWORD disposition = mode == A_READ ? OPEN_EXISTING : CREATE_ALWAYS;
    
    // CHANGED: Removed FILE_FLAG_OVERLAPPED to force synchronous (blocking) mode.
    // This ensures data is fully read/written before the function returns.
    result = CreateFile( 
        path, 
        access, 
        FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
        nullptr,
        disposition,
        FILE_ATTRIBUTE_NORMAL, // Was FILE_FLAG_OVERLAPPED
        nullptr );

    if( result == INVALID_HANDLE_VALUE ) {
    	DWORD err = GetLastError();
    	printf( "Error opening %s, code %u\n", path, err );
        return( 0 );
    }
    file->w32handle = result;
    return( 1 );
}


int CloseWIN32( aFILE * file ) {
    return( CloseHandle( file->w32handle ) );
}


int WriteWIN32( aFILE * file, aIORec * rec ) {
    DWORD bytesWritten = 0;
    OVERLAPPED ov = {0}; // Use stack-allocated OVERLAPPED for offset specification

    if( file->mode != A_WRITE ) {
        return( 0 );
    }
	
    // Handle 64-bit offsets
	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		ov.OffsetHigh = (DWORD)(tmppos >> 32);
		ov.Offset = (DWORD)(tmppos & 0xFFFFFFFF);
	} else {
        // If CURRENT_POS, we generally let the file pointer be. 
        // However, when using OVERLAPPED with blocking files, it updates the file ptr 
        // if we don't specify an offset. But to be safe with libMems logic:
        ov.Offset = 0xFFFFFFFF; 
        ov.OffsetHigh = 0xFFFFFFFF;
    }

    // WriteFile in blocking mode returns non-zero on success
    if( WriteFile( 
        file->w32handle, 
        rec->buf, 
        (DWORD)(rec->size * rec->count), 
        &bytesWritten,
        (rec->pos == CURRENT_POS) ? nullptr : &ov ) == 0 ) {
        
        DWORD err = GetLastError();
        printf( "error with WriteFile: %u\n", err );
        return( 0 );
    }
    return( 1 );
}


int ReadWIN32( aFILE * file, aIORec * rec ) {
    DWORD bytesRead = 0;
    OVERLAPPED ov = {0};

    if( file->mode != A_READ ) {
        return( 0 );
    }

    // Handle 64-bit offsets
	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		ov.OffsetHigh = (DWORD)(tmppos >> 32);
		ov.Offset = (DWORD)(tmppos & 0xFFFFFFFF);
	}

    // ReadFile in blocking mode halts execution until data is in the buffer.
    // This guarantees 'rec->buf' is full when we return.
    if( ReadFile( 
        file->w32handle, 
        rec->buf, 
        (DWORD)(rec->size * rec->count), 
        &bytesRead,
        (rec->pos == CURRENT_POS) ? nullptr : &ov ) == 0 ) {
        
        DWORD err = GetLastError();
        if (err == ERROR_HANDLE_EOF) {
            // EOF is fine, just return success with 0 bytes read (handled by caller)
            return 1; 
        }
        printf( "error with ReadFile -- Last Error: %u\n", err );
        return( 0 );
    }
    return( 1 );
}


int QueryLastCompleteWIN32( aFILE * file ) {
    // Since we are now using Synchronous I/O, the operation is ALWAYS complete
    // by the time the Read/Write function returns.
    return( 1 );
}

#endif /* USE_WIN32 */
