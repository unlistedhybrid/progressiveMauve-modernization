#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/sorting.h"
#include <cmath>
#include <cstring>

#ifndef USE_QSORT_ONLY

void RadixHistogram( sort_buf_t* sortbuf );
void RadixCopy( sort_buf_t* sortbuf );
void QSortPointers( sort_buf_t* sortbuf );
void QBrute( record_t* a[], int lo, int hi );
void QSort( record_t* a[], int lo0, int hi0 );
void CopySortedData ( sort_buf_t* sortbuf );


void InitRadixSort( sort_buf_t* sortbuf, buffer_t* scratch_buffer )
{
	unsigned int bin_divisor;
	unsigned int i, keyval = 0;
	sortbuf->histogram_size = 1;
	sortbuf->histogram_size <<= RADIX_BITS;
	sortbuf->histogram = static_cast<unsigned*>( malloc( sortbuf->histogram_size * sizeof(unsigned) ) );
	sortbuf->cur_ptr_offsets = static_cast<unsigned*>( malloc( sortbuf->histogram_size * sizeof(unsigned) ) );

	std::memset( sortbuf->histogram, 0, sortbuf->histogram_size * sizeof(unsigned) );
	

    bin_divisor = (unsigned)16777216 / (unsigned)NumBins;
    bin_divisor += (unsigned)16777216 % (unsigned)NumBins ? 1 : 0;

    for( i = 0; i < 3; i++ ) {
        keyval <<= 8;
        keyval += sortbuf->buf->recs[0].key[i];
    }
	sortbuf->base_number = (keyval / bin_divisor) * bin_divisor;
	
	sortbuf->divisor = (unsigned)bin_divisor / (unsigned)sortbuf->histogram_size;
	sortbuf->divisor += (unsigned)bin_divisor % (unsigned)sortbuf->histogram_size ? 1 : 0;
	
	sortbuf->cur_position = 0;
	sortbuf->sort_state = CalculateHistogram;
	sortbuf->radix_tmp = scratch_buffer;
	
	sortbuf->rec_ptrs = static_cast<record_t**>( malloc( sortbuf->buf->numrecs * sizeof(record_t*) ) );

}

void RadixSort( sort_buf_t* sortbuf )
{
	switch(sortbuf->sort_state){
		case CalculateHistogram:
			RadixHistogram( sortbuf );
			break;
		case CopyPointers:
			RadixCopy( sortbuf );
			break;
		case QsortPointers:
			QSortPointers( sortbuf );
			break;
		case CopyData:
			CopySortedData( sortbuf );
			break;
		default:
			printf("Error in sort_state\n");
	}
}

void RadixHistogram( sort_buf_t* sortbuf ){
	unsigned data_bucket;
	unsigned maxI;
	unsigned histI;
	unsigned cur_offset;
	unsigned tmp;
	record_t* cur_rec;
	
	maxI = sortbuf->cur_position + HISTOGRAM_CHUNK_SIZE;
	maxI = maxI < static_cast<unsigned>(sortbuf->buf->numrecs) ? maxI : static_cast<unsigned>(sortbuf->buf->numrecs);

	for(; sortbuf->cur_position < maxI; sortbuf->cur_position++){
		cur_rec = &(sortbuf->buf->recs[ sortbuf->cur_position ]);
		data_bucket = cur_rec->key[0];
		data_bucket <<= 8;
		data_bucket += cur_rec->key[1];
		data_bucket <<= 8;
		data_bucket += cur_rec->key[2];

		data_bucket -= sortbuf->base_number;
		data_bucket /= sortbuf->divisor;
		sortbuf->histogram[data_bucket]++;
	}
	
	if( sortbuf->cur_position == static_cast<unsigned>(sortbuf->buf->numrecs) ){

		cur_offset = 0;
		for( histI = 0; histI < sortbuf->histogram_size; histI++){
			tmp = sortbuf->histogram[ histI ];
			sortbuf->histogram[ histI ] = cur_offset;
			cur_offset += tmp;
		}

		sortbuf->sort_state = CopyPointers;
		sortbuf->cur_position = 0;
	}
}

void RadixCopy( sort_buf_t* sortbuf ){

	unsigned data_bucket;

	unsigned maxI;
	record_t* cur_rec;

	maxI = sortbuf->cur_position + PTR_COPY_CHUNK_SIZE;
	maxI = maxI < static_cast<unsigned>(sortbuf->buf->numrecs) ? maxI : static_cast<unsigned>(sortbuf->buf->numrecs);

	if(sortbuf->cur_position == 0 )
		std::memcpy(sortbuf->cur_ptr_offsets, sortbuf->histogram, sortbuf->histogram_size * sizeof(unsigned) );

	for(; sortbuf->cur_position < maxI; sortbuf->cur_position++){
		cur_rec = &(sortbuf->buf->recs[ sortbuf->cur_position ]);
		data_bucket = cur_rec->key[0];
		data_bucket <<= 8;
		data_bucket += cur_rec->key[1];
		data_bucket <<= 8;
		data_bucket += cur_rec->key[2];

		data_bucket -= sortbuf->base_number;
		data_bucket /= sortbuf->divisor;
		
		sortbuf->rec_ptrs[ sortbuf->cur_ptr_offsets[ data_bucket ] ] = cur_rec;
		sortbuf->cur_ptr_offsets[ data_bucket ]++;
	}
	
	if( sortbuf->cur_position == static_cast<unsigned>(sortbuf->buf->numrecs) ){
		sortbuf->sort_state = QsortPointers;
		sortbuf->cur_position = 0;
	}
	
}

void QSortPointers( sort_buf_t* sortbuf )
{
	unsigned binI = sortbuf->cur_position;
	unsigned maxI = binI + SORT_BINS_SIZE;

	maxI = maxI < sortbuf->histogram_size ? maxI : sortbuf->histogram_size - 1;

	for(; binI < maxI; binI++){
		if( sortbuf->histogram[binI + 1] - sortbuf->histogram[binI] > 1 )
			QSort( sortbuf->rec_ptrs, sortbuf->histogram[binI], sortbuf->histogram[binI + 1] - 1 );
	}
	sortbuf->cur_position = binI;

	if( binI == sortbuf->histogram_size - 1 ){
		if( (sortbuf->buf->numrecs - 1) - sortbuf->histogram[binI] > 1 )
			QSort( sortbuf->rec_ptrs, sortbuf->histogram[binI], sortbuf->buf->numrecs - 1 );
		sortbuf->sort_state = CopyData;
		sortbuf->cur_position = 0;
	}
}


void CopySortedData ( sort_buf_t* sortbuf ){
	unsigned recordI = sortbuf->cur_position;
	unsigned maxI = recordI + COPY_CHUNK_SIZE;
	record_t* tmp;
	
	maxI = maxI < static_cast<unsigned>(sortbuf->buf->numrecs) ? maxI : static_cast<unsigned>(sortbuf->buf->numrecs);

	for(; recordI < maxI; recordI++ )
		sortbuf->radix_tmp->recs[recordI] = *(sortbuf->rec_ptrs[recordI]);

	sortbuf->cur_position = recordI;

	if(recordI == static_cast<unsigned>(sortbuf->buf->numrecs)){
		tmp = sortbuf->radix_tmp->recs;
		sortbuf->radix_tmp->recs = sortbuf->buf->recs;
		sortbuf->buf->recs = tmp;

		sortbuf->state = WRITE_RESTRUCTURE;
		
		free( sortbuf->rec_ptrs );
		free( sortbuf->histogram );
		free( sortbuf->cur_ptr_offsets );
	}


}


void QBrute( record_t* a[], int lo, int hi ) {
    if ((hi-lo) == 1) {
        if( CompareKeyPtrs( a[hi], a[lo] ) < 0 ) {
            record_t* T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }else
    if ((hi-lo) == 2) {
        int pmin = CompareKeyPtrs( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = CompareKeyPtrs( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        if (pmin != lo) {
            record_t* T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        QBrute(a, lo+1, hi);
    }else
    if ((hi-lo) == 3) {
        int pmin, pmax;
        pmin = CompareKeyPtrs( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = CompareKeyPtrs( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        pmin = CompareKeyPtrs( a[pmin], a[lo+3] ) < 0 ? pmin : lo+3;
        if (pmin != lo) {
            record_t* T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        pmax = CompareKeyPtrs( a[hi], a[hi-1] ) > 0 ? hi : hi-1;
        pmax = CompareKeyPtrs( a[pmax], a[hi-2] ) > 0 ? pmax : hi-2;
        if (pmax != hi) {
            record_t* T = a[hi];
            a[hi] = a[pmax];
            a[pmax] = T;
        }
        QBrute(a, lo+1, hi-1);
    }
}



void QSort( record_t* a[], int lo0, int hi0 ) {
    
    int lo = lo0;
    int hi = hi0;
    
    record_t* pivot;

    if ((hi-lo) <= 3) {
        QBrute(a, lo, hi);
        return;
    }
    
    pivot = a[(lo + hi) / 2];
    a[(lo + hi) / 2] = a[hi];
    a[hi] = pivot;
    
    while( lo < hi ) {

        while( (CompareKeyPtrs( a[lo], pivot ) <= 0) && lo < hi ) {
            lo++;
        }
        
        while( (CompareKeyPtrs( pivot, a[hi] ) <= 0) && lo < hi ) {
            hi--;
        }
        
        if( lo < hi ) {
            record_t* T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }
    
    a[hi0] = a[hi];
    a[hi] = pivot;
    
    QSort( a, lo0, lo-1 );
    QSort( a, hi+1, hi0 );
}


#endif /* !USE_QSORT_ONLY */
