#ifndef _alinuxaio_h_
#define _alinuxaio_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"

#ifdef __cplusplus
extern "C" {
#endif

int OpenLinux( aFILE * file, const char *path, int mode );
int CloseLinux( aFILE * file );

int WriteLinux( aFILE * file, aIORec * rec );
int ReadLinux( aFILE * file, aIORec * rec );

int QueryLastCompleteLinux( aFILE * file );

#ifdef __cplusplus
}
#endif

#endif /* _alinuxaio_h_ */
