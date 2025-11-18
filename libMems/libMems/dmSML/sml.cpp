#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/sml.h"
#include "libMems/SeedMasks.h"
#include <cstring>


SMLHeader_t InitSML( aFILE* file, uint64 file_size, uint64 seed ){
	SMLHeader_t header;
	int retcode;
	
	header.version = 5;
	header.alphabet_bits = 2;
	header.seed = seed;
	header.seed_length = getSeedLength( seed );
	header.seed_weight = getSeedWeight( seed );
	header.length = file_size;
	header.unique_mers = -1;
	header.word_size = 32;
	header.little_endian = 1;
	header.id = 0;
	header.circular = 0;
	std::memcpy(header.translation_table, CreateBasicDNATable(), UINT8_MAX);
	header.description[ 0 ] = 0;
	
	retcode = aWrite( static_cast<void*>(&header), sizeof( header ), 1, file, 0 );
	if( retcode == 0 )
		printf( "Error writing to SML\n" );
	aWaitComplete( file, retcode );
	return header;
}
