#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>     // For std::string (used in InitdmSML)
#include <algorithm>  // For std::memset (used in InitdmSML and dmSML)

// FIX: Ensure this type is declared as a global variable outside of the sml namespace
// because it is assigned using a function (sml::CreateBasicDNATable) and used globally.
unsigned char* DNA_TABLE = nullptr; 

// Include necessary project headers
#include "libMems/dmSML/util.h"
#include "libMems/dmSML/timing.h"
#include "libMems/dmSML/asyncio.h"
#include "libMems/dmSML/buffer.h"
#include "libMems/dmSML/sorting.h"
#include "libMems/dmSML/sml.h"
#include "libMems/dmSML/dmsort.h"


device_t *Devices;
int NumDevices;

int NSortBufs;
sort_buf_t *SortBufs;

offset_t BufferSizeMin;
offset_t BufferSizeMax;

bin_t   *Bins;
int     NumBins;
int     NumBinDevs;

seqbuf_t Seqbuf;
    
aFILE   *Data;
int     DataDev;

const char *OutFileName = "unset";
aFILE   *Output;
int     OutputDev;

int BinToRead, BinToWrite, BinToSort;

working_set_t   WS;

offset_t NumRecs;
offset_t RecsProcessed;
offset_t RecsRead;
offset_t RecsUnread;
offset_t RecsCommitted;
offset_t RecsWritten;

double RunningTime;
dmtimer_t *RunningTimer;
double BinningTime;
dmtimer_t *BinningTimer;
double SortingTime;
dmtimer_t *SortingTimer;

double QSortTime;
dmtimer_t *QSortTimer;

double ReadIdleTime;
dmtimer_t *ReadIdleTimer;
double SortIdleTime;
dmtimer_t *SortIdleTimer;
double WriteIdleTime;
dmtimer_t *WriteIdleTimer;

buffer_list_t   Free;
buffer_list_t   ToProcess;
buffer_list_t   Reading;
buffer_list_t   Restructure;

static buffer_t * AllocateFree( void ) {
    buffer_t * ret;
    if( Free.nitems ) {
        ret = PopHead( &Free );
    } else {
        printf( "error: called AllocateFree but free list is empty\n" );
        return( nullptr );
    }
    ret->device = nullptr;
    ret->file = nullptr;
    ret->last = ret->next = nullptr;
    ret->numrecs = 0;
    ret->operation = OP_NONE;
    return( ret );
}

static size_t divisor = 0;

static int ComputeBinNumber( const unsigned char key[10] ) {
    int i;
    size_t keyval = 0;
    
    if( divisor == 0 ) {
        divisor = 16777216 / (size_t)NumBins;
        divisor += (16777216 % (size_t)NumBins) ? 1 : 0;
        printf( "Divisor is: %zu\n", divisor );
    }
    
    for( i = 0; i < 3; i++ ) {
        keyval <<= 8;
        keyval += key[i];
    }
    
    return (int)(keyval / divisor);
}

static offset_t         consumed_recs = 0;
static buffer_t         *toprocess = nullptr;

static void DoBinning( void ) {
    while( 1 ) {
        int bin = -1;
        if( toprocess == nullptr ) {
            if( ToProcess.nitems ) {
                toprocess = PopHead( &(ToProcess) );
                consumed_recs = 0;
            } else {
                return;
            }
        }
        for( ; consumed_recs < toprocess->numrecs; consumed_recs++, RecsProcessed++ ) {
            
            buffer_t *headbuf;
            record_t *rec = &(toprocess->recs[consumed_recs]);
            
#ifdef ASCII_KEYBYTES
            bin = ComputeAsciiBinNumber( rec->key );
#else
#ifdef NNNNN_KEYBYTES
			bin = ComputeNNNNNBinNumber( rec->key );
#else
            bin = ComputeBinNumber( rec->key );
#endif
#endif
            if( (bin >= static_cast<int>(NumBins)) || (bin < 0) ) {
                printf( "error: invalid bin from ComputeBinNumber: %d\n", bin );
            }

            headbuf = Bins[bin].bufs.head;
            if( !headbuf || 
                headbuf->numrecs == headbuf->totalrecs || 
                headbuf->operation != OP_NONE ) {
                if( headbuf->operation == BIN_SPECIAL ) {
                    headbuf->numrecs = 0;
                    headbuf->operation = OP_NONE;
                } else {
                    if( Free.nitems ) {
                        PushHead( &(Bins[bin].bufs), AllocateFree() );
                        headbuf = Bins[bin].bufs.head;
                    } else {
                        return;
                    }
                }
            }
            headbuf->recs[headbuf->numrecs++] = *rec;
            Bins[bin].nrecs++;
            if( headbuf->numrecs >= headbuf->totalrecs ) {
                headbuf->file = Bins[bin].file;
                headbuf->device = &(Devices[Bins[bin].dev].dev);
                RecsCommitted += headbuf->numrecs;
#ifdef NO_BIN_WRITE_PERF_TEST
				headbuf->operation = OP_FINISHED;
#else
                WriteBuffer( headbuf, headbuf->numrecs, headbuf->device );
#endif
                headbuf = nullptr;
            }
            
        }
        
        if( consumed_recs >= toprocess->numrecs ) {
            PushTail( &Free, toprocess );
            toprocess = nullptr;
        }
    }
}

void FinishBinning() {
    int i;
    buffer_t *b;
    offset_t recs = 0;
    for( i = 0; i < NumBins; i++ ) {
        while( Bins[i].bufs.nitems ) {
            b = PopHead( &(Bins[i].bufs) );
            if( b->operation == OP_NONE && b->numrecs ) {
                recs += b->numrecs;
                b->file = Bins[i].file;
                b->device = &(Devices[Bins[i].dev].dev);
#ifdef NO_BIN_WRITE_PERF_TEST
				b->operation = OP_FINISHED;
#else
                WriteBuffer( b, b->numrecs, b->device );
#endif
            }
        }
    }
    RecsCommitted += recs;
}

offset_t CalculateDataReadSize( buffer_t* b ){
	return (b->totalrecs + sml::mask_length - 1) < (RecsUnread + sml::mask_length - 1) ? (b->totalrecs + sml::mask_length - 1) : (RecsUnread + sml::mask_length - 1);
}

static void DoReading( void ) {
    buffer_t * b;
    if( RecsUnread && Free.nitems ) {
        b = AllocateFree();
        
        b->file = Data;
        ReadBuffer( b, (b->totalrecs < RecsUnread) ? b->totalrecs : RecsUnread, &(Devices[DataDev].dev) );

        b->input_pos = NumRecs - RecsUnread;
        b->io_pos = b->input_pos;
        b->io_size = CalculateDataReadSize( b );
        RecsUnread -= ((b->totalrecs < RecsUnread) ? b->totalrecs : RecsUnread);
        
        PushTail( &Reading, b );
    }
}

static void HandleBinWriteCompletions( void ) {
    int i;
    buffer_t *b, *tmpnext;
    for( i = 0; i < NumBins; i++ ) {
        b = Bins[i].bufs.head;
        do {
            if( !b ) {
                break;
            }
            tmpnext = b->next;
            if( b->operation == OP_FINISHED ) {
                RecsWritten += b->numrecs;
                if( Bins[i].bufs.nitems > 1 ) {
                    b->operation = OP_NONE;
                    PushHead( &Free, RemoveItem( &(Bins[i].bufs), b ) );
                } else {
                    b->operation = BIN_SPECIAL;
                }
            }
            b = tmpnext;
        } while( b != Bins[i].bufs.head && Bins[i].bufs.nitems > 1 );
    }
}

static void HandleSeqbufWriteCompletions( void ) {
    buffer_t *b, *tmpnext;
    b = Seqbuf.bufs.head;
    do {
        if( !b ) {
            break;
        }
        tmpnext = b->next;
        if( b->operation == OP_FINISHED ) {
            if( Seqbuf.bufs.nitems > 1 ) {
                b->operation = OP_NONE;
                PushHead( &Free, RemoveItem( &(Seqbuf.bufs), b ) );
            } 
        }
        b = tmpnext;
    } while( b != Seqbuf.bufs.head && Seqbuf.bufs.nitems > 1 );
}

#define ALPHA_BITS 2

static void Translate32(uint32* dest, const char* src, const unsigned len){
	uint8 start_bit = 0;
	unsigned cur_word = 0;
	uint32 word_mer = 0;
	uint32 i = 0;
	if( len == 0 )
		return;
	for(i=0; i < len; i++){
		if(start_bit + ALPHA_BITS <= 32){
			word_mer <<= ALPHA_BITS;
			// FIX: Use the global DNA_TABLE
			word_mer |= DNA_TABLE[static_cast<unsigned char>(src[i])];
			dest[cur_word] = word_mer;
			start_bit += ALPHA_BITS;
			if(start_bit >= 32 && i < len - 1){
				word_mer = 0;
				start_bit %= 32;
				cur_word++;
			}
		}else{
			printf("Error, this should never happen with DNA sequence\n" );
		}
	}
	if( start_bit != 0 ){
		dest[cur_word] <<= 32 - start_bit;
	}
}

void RestructureReadSMLBins( void ) {
	char little_endian = 1;
	sml::mask_t bit;
	sml::mask_t mer, rc_mer;
	record_t forward, reverse;
	record_t begin[6];
	int i;
	offset_t seqI, extras;
	char* sequence;
	sml::sml_t *sml;

    buffer_t *b, *tmpnext;
	
    buffer_t *headbuf;
	int seq_bit;
	int seq_word;
	int word_remainder;
	offset_t translate_length;
	
    b = Restructure.head;
    do {
        if( !b ) {
            break;
        }
		if( static_cast<uint64>(b->input_pos) != Seqbuf.seq_pos ){
			b = b->next;
			continue;
		}
		
        tmpnext = b->next;
		sequence = reinterpret_cast<char*>(b->recs);
		sml = reinterpret_cast<sml::sml_t*>(b->recs);

        headbuf = Seqbuf.bufs.head;
        if( !headbuf || 
            headbuf->operation != OP_NONE ) {
            if( headbuf->operation == OP_FINISHED ) {
                headbuf->numrecs = 0;
                headbuf->operation = OP_NONE;
	            Seqbuf.bufpos = 0;
            } else {
                if( Free.nitems ) {
                    PushHead( &(Seqbuf.bufs), AllocateFree() );
                    headbuf = Seqbuf.bufs.head;
		            Seqbuf.bufpos = 0;
                } else {
                    return;
                }
            }
        }

		seq_bit = Seqbuf.bufpos * 2;
		seq_word = seq_bit / 32;
		word_remainder = seq_bit % 32;
		if( word_remainder != 0 ){
			seq_word++;
		}
		
		translate_length = b->io_size - static_cast<offset_t>(sml::mask_length) + 1 - (word_remainder / 2);
		if( b->io_size + b->input_pos >= NumRecs ){
			translate_length += static_cast<offset_t>(sml::mask_length) - 1;
		}
		
#ifndef NO_RESTRUCTURE_PERF_TEST
		Translate32( reinterpret_cast<uint32*>(headbuf->recs) + seq_word, reinterpret_cast<char*>(b->recs) + (word_remainder / 2), translate_length );
#endif
		
		if( word_remainder != 0 ){
			int begin_mer = 0;
			for( seqI = 0; seqI < static_cast<offset_t>(word_remainder / 2); seqI++ ){
				begin_mer <<= 2;
				// FIX: Use the global DNA_TABLE
				begin_mer |= DNA_TABLE[static_cast<unsigned char>(sequence[ seqI ])];
			}
			(reinterpret_cast<uint32*>(headbuf->recs))[ seq_word - 1 ] |= begin_mer;
		}
		
		Seqbuf.bufpos += translate_length + (word_remainder / 2);
		Seqbuf.seq_pos += translate_length + (word_remainder / 2);

        if( Seqbuf.bufpos == headbuf->totalrecs * sizeof( record_t ) * 4 ||
			b->io_size + b->input_pos >= NumRecs ) {
            headbuf->file = Seqbuf.file;
            headbuf->device = &(Devices[Seqbuf.dev].dev);
            WriteBuffer( headbuf, headbuf->totalrecs, headbuf->device );
            headbuf->io_size = Seqbuf.bufpos / 4;
            if( b->io_size + b->input_pos >= NumRecs ){
            	offset_t offI = 0;
            	offI = headbuf->io_size % 4;
            	if( offI != 0 )
	            	headbuf->io_size += 4 - offI;
            	for( offI = 0; offI < 8; offI++ )
            		(reinterpret_cast<char*>(headbuf->recs))[ headbuf->io_size + offI ] = 0;
            	headbuf->io_size += 8;
            }
            headbuf = nullptr;
        }else if( Seqbuf.bufpos > headbuf->totalrecs * sizeof( record_t ) * 4 ){
        	printf( "Error.  Over filled Seqbuf\n" );
        }

#ifndef NO_RESTRUCTURE_PERF_TEST
		for( seqI = b->io_size - static_cast<offset_t>(sml::mask_length) + 1; seqI > 0; seqI-- ){
			bit = 1;
			bit <<= sml::mask_length - 1;
			mer = 0;
			for( i = 0; i < sml::mask_length; i++ ){
				if( bit & sml::seed_mask ){
					mer <<= 2;
					// FIX: Use the global DNA_TABLE
					mer |= DNA_TABLE[static_cast<unsigned char>(sequence[ seqI + i - 1 ])];
				}
				bit >>= 1;
			}
			// FIX: Use sml::mask_weight
			mer <<= 64 - (2 * sml::mask_weight);
			if( little_endian ){
				// FIX: Use sml::MASK_T_BYTES
				for( i = 0; i < sml::MASK_T_BYTES; i++ )
					forward.key[i] = (reinterpret_cast<char*>(&mer))[ sizeof( mer ) - i - 1 ];

			}else{
				// FIX: Use sml::MASK_T_BYTES
				for( i = 0; i < sml::MASK_T_BYTES; i++ )
					forward.key[i] = (reinterpret_cast<char*>(&mer))[ i ];
			}

			mer = ~mer;
			for( i = 0; i < 64; i += 2 ){
				rc_mer <<= 2;
				rc_mer |= mer & 3;
				mer >>= 2;
			}
			// FIX: Use sml::mask_weight
			rc_mer <<= 64 - (2 * sml::mask_weight);
			if( little_endian ){
				// FIX: Use sml::MASK_T_BYTES
				for( i = 0; i < sml::MASK_T_BYTES; i++ )
					reverse.key[i] = (reinterpret_cast<char*>(&rc_mer))[ sizeof( mer ) - i - 1 ];
			}else{
				// FIX: Use sml::MASK_T_BYTES
				for( i = 0; i < sml::MASK_T_BYTES; i++ )
					reverse.key[i] = (reinterpret_cast<char*>(&rc_mer))[i];
			}
			if( COMPARE_KEYS( forward, reverse ) > 0)
				forward = reverse;
			
			if( seqI <= 6 ){
				begin[ seqI - 1] = forward;
			}else{
				b->recs[ seqI - 1 ] = forward;
				sml[ seqI - 1 ].pos = b->input_pos + seqI - 1;
			}
		}

		extras = b->io_size - static_cast<offset_t>(sml::mask_length) + 1 < 6 ? b->io_size - static_cast<offset_t>(sml::mask_length) + 1 : 6;
		
		for(; seqI < static_cast<offset_t>(extras); seqI++ ){
			b->recs[ seqI ] = begin[ seqI ];
			sml[ seqI ].pos = b->input_pos + seqI;
		}
#else
	if(1){
    int i;
    size_t keyval = 0;
	size_t tmpval = 0;
    if( divisor == 0 ) {
        divisor = (unsigned)16777216 / (unsigned)NumBins;
        divisor += (unsigned)16777216 % (unsigned)NumBins ? 1 : 0;
        printf( "Divisor is: %u\n", divisor );
    }
	for( seqI = 0; seqI < b->numrecs; seqI++ ){
		tmpval = keyval;
	    for( i = 3; i > 0; i-- ) {
			b->recs[ seqI ].key[ i - 1 ] = (tmpval & 0xFF);
			b->recs[ seqI ].key[ i - 1 ] = 0;
			tmpval >>= 8;
		}
		keyval += divisor;
	}
	}
#endif
		
        PushTail( &ToProcess, RemoveItem( &Restructure, b ) );
        b = tmpnext;
    } while( b != Restructure.head && Restructure.nitems );
}

static void HandleReadingCompletions( void ) {
    buffer_t *b, *tmpnext;
    b = Reading.head;
    do {
        if( !b ) {
            break;
        }
        tmpnext = b->next;
        if( b->operation == OP_FINISHED ) {
            b->operation = OP_NONE;
            PushTail( &Restructure, RemoveItem( &Reading, b ) );
            RecsRead += b->numrecs;
        }
        b = tmpnext;
    } while( b != Reading.head && Reading.nitems );
}

extern "C" {
int InitdmSML( long working_mb, long buffer_size, const char* input_filename, const char* output_filename, const char* const* scratch_paths, sml::uint64 seed ) {
    int i, j;
    offset_t desired_ws_size, actual_ws_size;
    // FIX: Use sml::SMLHeader_t
    sml::SMLHeader_t header;
    struct {
        const char * bin_dev;
        int devnum;
        int nbins;
    } bins[8];

	char *bin_name;
	int scratchI = 0;

    InitTime();

    RunningTime = 0;
    RunningTimer = StartTimer();

	if( working_mb != 0 ){
	desired_ws_size = working_mb;
	desired_ws_size *= 1024 * 1024;
	}else{
#ifdef WIN32
	{
	MEMORYSTATUS ms;
	std::memset( &ms, 0, sizeof( MEMORYSTATUS ) );
	GlobalMemoryStatus( &ms );
	desired_ws_size = ms.dwTotalPhys / 2;
	}
#else
    {
	// Use std::string for safer file handling if possible, but stick to C-style for compatibility here
	FILE *fp = fopen("/proc/meminfo", "r");
	if ( fp )
	{
		long memTotal;
		char buf[1024];
		if ( fgets(buf, sizeof(buf), fp) )
		{
			sscanf(buf, "MemTotal: %ld kB", &memTotal);
			fprintf( stderr, "%s", buf );
		}
		fclose(fp);
		desired_ws_size = memTotal * 512;
	}
	}
#endif
	if( desired_ws_size / 1024  > 2048 * 1024 ){
		desired_ws_size = 1024 * 1024;
		desired_ws_size *= 2048;
	}
	}
	
	if( buffer_size == 0 ){
		buffer_size = 1;
		while( desired_ws_size / (buffer_size*sizeof(record_t)) > 2048 ){
			buffer_size *= 2;
		}
	}

	BufferSizeMin = BufferSizeMax = buffer_size;
	OutFileName = output_filename;
	
	for( ; ; scratchI++ ){
		if( !scratch_paths || scratch_paths[ scratchI ] == nullptr )
			break;
	}
	
    
	NumBinDevs = scratchI;
	NumDevices = 2 + NumBinDevs;
	Devices = static_cast<device_t*>(malloc( NumDevices * sizeof(device_t) ));
	DataDev = 0;
	OutputDev = 1;
	Devices[DataDev].devname = "Input device";
	Devices[DataDev].path = input_filename;
	Devices[DataDev].dev.buf = nullptr;
	Devices[OutputDev].devname = "Output device";
	Devices[OutputDev].path = OutFileName;
	Devices[OutputDev].dev.buf = nullptr;
    
    if( NumBinDevs == 0 ) {
    	return TOO_FEW_BINS;
    } else if( NumBinDevs > 8 ) {
    	return TOO_MANY_BINS;
    }
	
	NumRecs = aStatFileSize( input_filename );

	NumBins = desired_ws_size / (200 * NumBinDevs);
	NumBins = NumRecs / NumBins;
	NumBins = NumBins < 5 * NumBinDevs ? 5 * NumBinDevs : NumBins;
	if( NumBins % NumBinDevs != 0 )
		NumBins = ( (NumBins / NumBinDevs) + 1 ) * NumBinDevs;
	printf( "Creating %d bin files\n", NumBins );
	for( i = 2; i < NumDevices; i++ ){
		bin_name = static_cast<char*>(malloc( 10 ));
		strcpy( bin_name, "bin dev__" );
		bin_name[8] = 0x40 + i - 2;
		Devices[i].devname = bin_name;
		Devices[i].path = scratch_paths[ i - 2 ];
		Devices[i].dev.buf = nullptr;
		bins[i - 2].bin_dev = bin_name;
		bins[i - 2].nbins = NumBins / NumBinDevs;
		bins[i - 2].devnum = i;
	}
	
    if( BufferSizeMin == 0 ) {
        BufferSizeMin = MINRECS;
        BufferSizeMax = MAXRECS;
    }

    Data = aOpen( input_filename, A_READ );
	if( Data == nullptr ) {
	        printf( "couldn't open data file\n" );
		return INPUT_NOT_OPENED;
	}
   
    if( desired_ws_size == 0 ) {
        printf( "invalid working set size (%llu) -- must be at least 0\n", desired_ws_size );
    	return INVALID_WS_SIZE;
    }
	
	// FIX: Use sml::CreateBasicDNATable() since it's now in the namespace
	DNA_TABLE = sml::CreateBasicDNATable();

    Output = aOpen( OutFileName, A_WRITE );
    if( !Output ) {
        printf( "couldn't open output file!\n" );
    	return OUTPUT_NOT_OPENED;
    }
	
	header = sml::InitSML( Output, NumRecs, seed );
	sml::seed_mask = header.seed;
	sml::mask_length = header.seed_length;
	sml::mask_weight = header.seed_weight;
	
	if( NumRecs <= static_cast<offset_t>(sml::mask_length) - 1 ){
	        printf( "Sequence must be at least %d characters in length\n", sml::mask_length );
		return SEQUENCE_TOO_SHORT;
	}

	NumRecs -= static_cast<offset_t>(sml::mask_length) - 1;
	printf( "NumRecs is: %llu \n", NumRecs );
    RecsProcessed = 0;
    RecsUnread = NumRecs;
    if( NumRecs <= 0 ) {
    	return INVALID_NUMRECS;
    }
    
    actual_ws_size = MakeWorkingSet( &WS, desired_ws_size, BufferSizeMin, BufferSizeMax );
    printf( "desired working set: %llu, actual working set: %llu\n", 
        desired_ws_size, actual_ws_size );

    for( i = 0; i < WS.nbufs; i++ ) {
        PushHead( &Free, &(WS.bufs[i]) );
    }

    printf( "working set size        : %llu\n", actual_ws_size );
    printf( "total buffers           : %d\n", WS.nbufs );
    ToProcess.nitems = Reading.nitems = 0;
    ToProcess.head = Reading.head = nullptr;
	Restructure.nitems = 0;
	Restructure.head = nullptr;
		
	Seqbuf.file = Output;
	Seqbuf.dev = OutputDev;
	Seqbuf.bufpos = 0;
	Seqbuf.seq_pos = 0;
    if( Free.nitems ) {
        PushHead( &(Seqbuf.bufs), AllocateFree() );
    } else {
        printf( "error: could not give a buffer to Seqbuf\n" );
        return NO_FREE_BUFFERS;
    }

    Bins = static_cast<bin_t*>(malloc( sizeof( *Bins ) * NumBins ));
    std::memset( Bins, 0, sizeof( *Bins ) * NumBins );

    printf( "opening %d bins\n", NumBins );
    j = -1;
    for( i = 0; i < NumBins; i++ ) {
        while( 1 ) {
            j = (j+1) % NumBinDevs;
            if( bins[j].nbins ) {
                const char *fname = Fmt("%sout%05d.binned",Devices[bins[j].devnum].path,i);
                Bins[i].dev = bins[j].devnum;
                Bins[i].fname = static_cast<char*>(malloc( strlen( fname ) + 1 ));
                strcpy( Bins[i].fname, fname );

#ifndef NO_BINNING_PERF_TEST
                Bins[i].file = aOpen( fname, A_WRITE );
                if( Bins[i].file == nullptr ) {
                    printf( "couldn't open output bin file '%s'\n", fname );
					return BIN_NOT_OPENED;
                }
#else
                Bins[i].nrecs = aStatSize( fname );
		if( Bins[i].nrecs == 0 ){
	                Bins[i].file = aOpen( fname, A_WRITE );
			aClose( Bins[i].file );
			Bins[i].file = nullptr;
		}
#endif
                bins[j].nbins--;
                break;
            }
        }
    }

    for( i = 0; i < NumBins; i++ ) {
        if( Free.nitems ) {
            PushHead( &(Bins[i].bufs), AllocateFree() );
        } else {
            printf( "error: could not give one buffer to each bin\n" );
	        return NO_FREE_BUFFERS;
        }
    }
	
	return 0;
}
}

void DisplayStatusHeader( void ) {
    printf( "time recs_read recs_processed recs_committed recs_written binning_rate free reading toprocess bins restructure\n" );
}

void DisplayStatus( void ) {

    printf( "%f %llu %llu %llu %llu %f %d %d %d %d %d\n",
        RunningTime, RecsRead, RecsProcessed, RecsCommitted, RecsWritten, 
        RecsProcessed/RunningTime, Free.nitems, Reading.nitems, ToProcess.nitems, 
        WS.nbufs - Free.nitems - Reading.nitems - ToProcess.nitems - Restructure.nitems, Restructure.nitems );
}

void UpdateIOState( void ) {
    int i;
    aUpdateOperations( Data );
    for( i = 0; i < NumBins; i++ ) {
        aUpdateOperations( Bins[i].file );
    }
    aUpdateOperations( Output );
    UpdateWSIOFinishedState( &WS );
    for( i = 0; i < NumDevices; i++ ) {
        UpdateDeviceIOExecuteState( &WS, &(Devices[i].dev) );
    }
}

void EnsureAllOperationsComplete( void ) {
    int i;
    int not_complete = 1;
    dmtimer_t *wait;
    wait = StartTimer();
    while( not_complete ) {
        UpdateIOState();
        not_complete = 0;
        for( i = 0; i < WS.nbufs; i++ ) {
            if( WS.bufs[i].device &&
                WS.bufs[i].file &&
                (WS.bufs[i].operation == OP_PENDING || WS.bufs[i].operation > OP_NONE) ) {
                not_complete = 1;
                break;
            }
        }
    }
    printf( "Ensure All Operations Complete: %d msec\n", ReadTimer( wait ) );
    StopTimer( wait );
}

static double lasttime = 0;

void BinningPhase( void ) {

    int i;
    int iter;

    printf( "----------------- Starting -----------------\n" );
    printf( "working set buffers : %d\n", WS.nbufs );
    printf( "number of bins      : %d\n", NumBins );
    iter = 0;
    DisplayStatusHeader();
    while( RecsProcessed < NumRecs ) {

        if( (RunningTime - lasttime) >= 2.0f ) {
            DisplayStatus();
            lasttime = RunningTime;
        }
        
        UpdateIOState();
        
        HandleReadingCompletions();
		HandleSeqbufWriteCompletions();
		RestructureReadSMLBins();
        HandleBinWriteCompletions();

        DoReading();
        DoBinning();
        
        iter++;

        RunningTime = static_cast<double>(ReadTimer( RunningTimer )) / 1000.0;

    }
    
    printf( "total iters: %d\n", iter );
    FinishBinning();
    EnsureAllOperationsComplete();

    aClose( Data );
    Data = nullptr;
    for( i = 0; i < NumBins; i++ ) {
        aClose( Bins[i].file );
        Bins[i].file = nullptr;
    }
    printf( "Finally, RecsCommitted: %llu\n", RecsCommitted );

    DisplayStatus();

}

void SortReading( void ) {

    int i;

    for( i = 0; i < NSortBufs; i++ ) {
        if( BinToRead >= NumBins ) {
            return;
        }
        if( SortBufs[i].state == WAIT_READ ) {
            const char *fname = Fmt("%sout%05d.binned",Devices[Bins[BinToRead].dev].path,BinToRead);
            aFILE *in = aOpen( fname, A_READ );
            if( !in ) {
                printf( "couldn't open '%s' to read!\n", fname );
            }
            if( Bins[BinToRead].nrecs > SortBufs[i].buf->totalrecs ) {
                printf( "buffer not big enough to hold bin!\n" );
            }
            SortBufs[i].bin = BinToRead;
            SortBufs[i].dev = &(Devices[Bins[BinToRead].dev].dev);
            SortBufs[i].state = BUSY_READ;
            SortBufs[i].buf->file = in;
            ReadBuffer( SortBufs[i].buf, Bins[BinToRead].nrecs, SortBufs[i].dev );
            printf( "scheduled read of bin %d\n", BinToRead );
            BinToRead++;
            return;
        }
    }

}

#ifdef USE_QSORT_ONLY

int comp_keys( record_t a, record_t b ){
	int compval;
	compval = COMPARE_KEYS( a, b );
	return compval;
}

void QBrute( record_t a[], int lo, int hi ) {
    if ((hi-lo) == 1) {
        if( comp_keys( a[hi], a[lo] ) < 0 ) {
            record_t T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }
    if ((hi-lo) == 2) {
        int pmin = comp_keys( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = comp_keys( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        if (pmin != lo) {
            record_t T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        QBrute(a, lo+1, hi);
    }
    if ((hi-lo) == 3) {
        int pmin, pmax;
        pmin = comp_keys( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = comp_keys( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        pmin = comp_keys( a[pmin], a[lo+3] ) < 0 ? pmin : lo+3;
        if (pmin != lo) {
            record_t T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        pmax = comp_keys( a[hi], a[hi-1] ) > 0 ? hi : hi-1;
        pmax = comp_keys( a[pmax], a[hi-2] ) > 0 ? pmax : hi-2;
        if (pmax != hi) {
            record_t T = a[hi];
            a[hi] = a[pmax];
            a[pmax] = T;
        }
        QBrute(a, lo+1, hi-1);
    }
}

void QSort( record_t a[], int lo0, int hi0 ) {
    
    int lo = lo0;
    int hi = hi0;
    
    record_t pivot;

    if ((hi-lo) <= 3) {
        QBrute(a, lo, hi);
        return;
    }
    
    pivot = a[(lo + hi) / 2];
    a[(lo + hi) / 2] = a[hi];
    a[hi] = pivot;
    
    while( lo < hi ) {

        while( (comp_keys( a[lo], pivot ) <= 0) && lo < hi ) {
            lo++;
        }
        
        while( (comp_keys( pivot, a[hi] ) <= 0) && lo < hi ) {
            hi--;
        }
        
        if( lo < hi ) {
            record_t T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }
    
    a[hi0] = a[hi];
    a[hi] = pivot;
    
    QSort( a, lo0, lo-1 );
    QSort( a, hi+1, hi0 );
}

void RecSort( record_t a[], int nelems ) {

    QSort( a, 0, nelems-1 );

}

int SortBuffer( buffer_t * buf ) {

    RecSort( buf->recs, buf->numrecs );
    return( 1 );

}

void SortSorting( void ) {

    int i, finished;
    int lowest = -1;
    QSortTimer = StartTimer();

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == SORTING ) {
            if( lowest == -1 || SortBufs[i].bin < SortBufs[lowest].bin ) {
                lowest = i;
            }
        }
    }

    if( lowest != -1 ) {
        printf( "sorting bin %d\n", SortBufs[lowest].bin );
        finished = SortBuffer( SortBufs[lowest].buf );
        if( finished ) {
            SortBufs[lowest].state = WRITE_RESTRUCTURE;
        }
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}

#elif defined NO_SORT_PERF_TEST

void SortSorting( void ) {
    
    int i;

    QSortTimer = StartTimer();

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == SORTING ) {
            SortBufs[i].state = WAIT_WRITE;
        }
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}

#else 

sort_buf_t* CurrentSortBuf;
buffer_t* SortScratchBuffer;

void SortSorting( void ) {

    int i;

    QSortTimer = StartTimer();

	if( CurrentSortBuf == nullptr ){
	    for( i = 0; i < NSortBufs; i++ ) {
	        if( SortBufs[i].state == SORTING && SortBufs[i].bin == BinToSort ) {
	        	CurrentSortBuf = &SortBufs[i];
	        	InitRadixSort( CurrentSortBuf, SortScratchBuffer );
	            printf( "scheduling sort of bin %d\n", BinToSort );
	            break;
	        }
	    }
	}
	
	if( CurrentSortBuf != nullptr ){
		if( CurrentSortBuf->state != WRITE_RESTRUCTURE ){

			RadixSort( CurrentSortBuf );

			if( CurrentSortBuf->state == WRITE_RESTRUCTURE ){
				CurrentSortBuf = nullptr;
	            BinToSort++;
	        }
		}
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}

#endif

void RestructureSMLBinsForWrite( void ) {
    int i;
    offset_t j;
	sml::position_t* positions;
	sml::sml_t *sml;

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == WRITE_RESTRUCTURE ) {
            printf( "restructuring bin %d\n", SortBufs[i].bin );
            positions = reinterpret_cast<sml::position_t*>(SortBufs[i].buf->recs);
            sml = reinterpret_cast<sml::sml_t*>(SortBufs[i].buf->recs);
            for( j = 0; j < Bins[SortBufs[i].bin].nrecs; j++ ){
            	positions[ j ] = sml[ j ].pos;
            }
            
            SortBufs[i].state = WAIT_WRITE;
        }
    }
}

int CalculateSortWriteSize( int sortI ){
     // FIX: Use sml::position_t size
     return Bins[SortBufs[sortI].bin].nrecs * sizeof( sml::position_t );
}

void SortWriting( void ) {

    int i;

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == WAIT_WRITE && SortBufs[i].bin == BinToWrite ) {
#ifdef NO_WRITE_PERF_TEST
			SortBufs[i].state = WAIT_READ;
#else
            printf( "scheduling write of bin %d\n", BinToWrite );
            SortBufs[i].dev = &(Devices[OutputDev].dev);
            SortBufs[i].state = BUSY_WRITE;
            SortBufs[i].buf->file = Output;
            WriteBuffer( SortBufs[i].buf, Bins[SortBufs[i].bin].nrecs, &(Devices[OutputDev].dev) );
			SortBufs[i].buf->io_size = CalculateSortWriteSize( i );
#endif
            BinToWrite++;
        }
    }

}

void SortHandleCompletions( void ) {

    int i;

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == BUSY_READ || SortBufs[i].state == BUSY_WRITE ) {
            if( SortBufs[i].buf->operation == OP_FINISHED ) {
                SortBufs[i].buf->operation = OP_NONE;
                SortBufs[i].state = SortBufs[i].state == BUSY_READ ? SORTING : WAIT_READ;
#ifdef NNNNN_KEYBYTES
				if( SortBufs[i].bin == 0 && SortBufs[i].state == SORTING )
					SortBufs[i].state = WAIT_WRITE;
#endif
            }
        }
    }

}

void SortUpdateIOState() {

    int i;
    
    aUpdateOperations( Output );
    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].buf->file ) {
            aUpdateOperations( SortBufs[i].buf->file );
        }
    }
    UpdateWSIOFinishedState( &WS );
    for( i = 0; i < NumDevices; i++ ) {
        UpdateDeviceIOExecuteState( &WS, &(Devices[i].dev) );
    }

}

void SortingEnsureAllOperationsComplete() {
    int i;
    int not_complete = 1;
    dmtimer_t *wait;
    wait = StartTimer();
    while( not_complete ) {
        SortUpdateIOState();
        not_complete = 0;
        for( i = 0; i < WS.nbufs; i++ ) {
            if( WS.bufs[i].device &&
                WS.bufs[i].file &&
                (WS.bufs[i].operation == OP_PENDING || WS.bufs[i].operation > OP_NONE) ) {
                not_complete = 1;
                break;
            }
        }
    }

    aFlush( Output );

    printf( "Sort Ensure All Operations Complete: %d msec\n", ReadTimer( wait ) );
    StopTimer( wait );
}

void SortingPhase( void ) {

    int i;
    offset_t recs_per_buffer;
    offset_t biggest_nrecs = 0;
    
    NSortBufs = NumBinDevs;

    for( i = 0; i < NumBins; i++ ) {
        if( Bins[i].nrecs > biggest_nrecs ) {
            biggest_nrecs = Bins[i].nrecs;
        }
    }
    
    recs_per_buffer = biggest_nrecs;

    if( (WS.size / sizeof( record_t )) < static_cast<unsigned>(recs_per_buffer) ) {
        printf( "working set holds %llu recs, but we need %llu\n", 
            (WS.size / sizeof( record_t )), recs_per_buffer );
    }

    NSortBufs = (WS.size / sizeof( record_t )) / recs_per_buffer;

    printf( "NSortBufs = %d\n", NSortBufs );

    BinToRead = 0;
    BinToWrite = 0;
	BinToSort = 0;
	
    printf( "reorganizing working set: %llu recs per buffer, %d sort bufs\n", recs_per_buffer, NSortBufs );
    ReorganizeWorkingSet( &WS, recs_per_buffer, recs_per_buffer );

#if !defined USE_QSORT_ONLY && !defined NO_SORT_PERF_TEST
    NSortBufs--;
	SortScratchBuffer = &(WS.bufs[NSortBufs]);
	SortScratchBuffer->operation = SORTING_SCRATCH;
#endif
    
    printf( "reorganized working set has %d buffers of %llu bytes\n", WS.nbufs, recs_per_buffer * sizeof(record_t) );
    SortBufs = static_cast<sort_buf_t*>(malloc( sizeof( *SortBufs ) * NSortBufs ));
    std::memset( SortBufs, 0, sizeof( *SortBufs ) * NSortBufs );

    for( i = 0; i < NSortBufs; i++ ) {

        SortBufs[i].state = WAIT_READ;
        SortBufs[i].buf = &(WS.bufs[i]);
        SortBufs[i].dev = nullptr;

    }
#ifdef NNNNN_KEYBYTES
    while( BinToWrite < 1 ) {
        SortReading();
        SortSorting();
        RestructureSMLBinsForWrite();
        SortWriting();
        SortUpdateIOState();
        SortHandleCompletions();
    }
    SortingEnsureAllOperationsComplete();

    for( i = 1; i < NumBins; i++ ) {
        if( Bins[i].nrecs > biggest_nrecs ) {
            biggest_nrecs = Bins[i].nrecs;
        }
    }
    recs_per_buffer = biggest_nrecs;
    if( (WS.size / sizeof( record_t )) < static_cast<unsigned>(recs_per_buffer) ) {
        printf( "working set holds %llu recs, but we need %llu\n", 
            (WS.size / sizeof( record_t )), recs_per_buffer );
    }
    NSortBufs = (WS.size / sizeof( record_t )) / recs_per_buffer;
    printf( "NSortBufs = %d\n", NSortBufs );
    BinToRead = 1;
    BinToWrite = 1;
	BinToSort = 1;
	
    printf( "reorganizing working set: %llu recs per buffer, %d sort bufs\n", recs_per_buffer, NSortBufs );
    ReorganizeWorkingSet( &WS, recs_per_buffer, recs_per_buffer );

#if !defined USE_QSORT_ONLY && !defined NO_SORT_PERF_TEST
    NSortBufs--;
	SortScratchBuffer = &(WS.bufs[NSortBufs]);
	SortScratchBuffer->operation = SORTING_SCRATCH;
#endif
    
    printf( "reorganized working set has %d buffers of %llu bytes\n", WS.nbufs, recs_per_buffer * sizeof(record_t) );
    SortBufs = static_cast<sort_buf_t*>(malloc( sizeof( *SortBufs ) * NSortBufs ));
    std::memset( SortBufs, 0, sizeof( *SortBufs ) * NSortBufs );

    for( i = 0; i < NSortBufs; i++ ) {
        SortBufs[i].state = WAIT_READ;
        SortBufs[i].buf = &(WS.bufs[i]);
        SortBufs[i].dev = nullptr;
    }
#endif    
    

    while( BinToWrite < NumBins ) {
        
        SortReading();
        
        SortSorting();
        
        RestructureSMLBinsForWrite();
        SortWriting();
        
        SortUpdateIOState();

        SortHandleCompletions();
        
    }

    SortingEnsureAllOperationsComplete();

    printf( "QSort took %f seconds\n", QSortTime );

}

int dmsort() {

    BinningTimer = StartTimer();
#ifndef NO_BINNING_PERF_TEST
    BinningPhase();
    BinningTime = ReadTimer( BinningTimer ) / 1000.0;
#endif
    StopTimer( BinningTimer );

    SortingTimer = StartTimer();
    SortingPhase();
    SortingTime = ReadTimer( SortingTimer ) / 1000.0;
    StopTimer( SortingTimer );

    RunningTime = ReadTimer( RunningTimer ) / 1000.0;
    StopTimer( RunningTimer );

    printf( "total time      : %f sec\n", RunningTime );
    printf( "binning time    : %f sec (%f%%)\n", BinningTime, BinningTime/RunningTime * sizeof(record_t) );
    printf( "sorting time    : %f sec (%f%%)\n", SortingTime, SortingTime/RunningTime * sizeof(record_t) );
    
    printf( "total rate      : %f MB/sec\n", (((double)NumRecs)/10485.760)/RunningTime );
    printf( "total bin rate  : %f MB/sec\n", (((double)NumRecs)/10485.760)/BinningTime );
    printf( "total sort rate : %f MB/sec\n", (((double)NumRecs)/10485.760)/SortingTime );

    return 0;
}

extern "C" {
int dmSML( const char* input_file, const char* output_file, const char* const* scratch_paths, sml::uint64 seed ) {
	int rval = 0;
	int i = 0;
	rval = InitdmSML( 0, 0, input_file, output_file, scratch_paths, seed );
	if( rval != 0 )
		return rval;
	rval = dmsort();
	
	// FIX: Must free DNA_TABLE memory allocated by sml::CreateBasicDNATable()
	if( DNA_TABLE )
		free( DNA_TABLE );
	DNA_TABLE = nullptr;

	for( i = 0; i < NumBins; i++ ){
		removeFile( Bins[ i ].fname, FALSE );
		free( Bins[ i ].fname );
	}
	if( Bins )
		free( Bins );
	Bins = nullptr;
	NumBins = 0;
	NumDevices = 0;
	if( Devices )
		free( Devices );
	Devices = nullptr;
	if( SortBufs )
		free( SortBufs );
	SortBufs = nullptr;

	NSortBufs = 0;

	BufferSizeMin = 0;
	BufferSizeMax = 0;
	
    std::memset( &Seqbuf, 0, sizeof( seqbuf_t ) );

	DataDev = 0;

	OutFileName = "unset";

    aClose( Output );
    Output = nullptr;
	OutputDev = 0;

	BinToRead = 0;
	BinToWrite = 0;
	BinToSort = 0;

	free( WS.bufs );
	std::memset( &WS, 0, sizeof( working_set_t ) );

	NumRecs = 0;
	RecsProcessed = 0;
	RecsRead = 0;
	RecsUnread = 0;
	RecsCommitted = 0;
	RecsWritten = 0;
	

	RunningTime = 0;
	RunningTimer= nullptr;
	BinningTime = 0;
	BinningTimer= nullptr;
	SortingTime = 0;
	SortingTimer = nullptr;

	QSortTime = 0;
	QSortTimer = nullptr;

	ReadIdleTime = 0;
	ReadIdleTimer = nullptr;
	SortIdleTime = 0;
	SortIdleTimer = nullptr;
	WriteIdleTime = 0;
	WriteIdleTimer = nullptr;
	
	
	std::memset( &Free, 0, sizeof( buffer_list_t ) );
	std::memset( &ToProcess, 0, sizeof( buffer_list_t ) );
	std::memset( &Reading, 0, sizeof( buffer_list_t ) );
	std::memset( &Restructure, 0, sizeof( buffer_list_t ) );

	divisor = 0;
	consumed_recs = 0;
	toprocess = nullptr;
	lasttime = 0;
	
	return rval;
}
}
