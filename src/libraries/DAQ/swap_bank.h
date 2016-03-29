

// This contains code for byte swapping EVIO banks. The code
// is mostly copied from HDEVIO, but modified to just print
// any errors rather than use the error flag/message system
// built into HDEVIO. It can therefore be used outside of
// HDEVIO.


#include <stdint.h>

#undef swap64
#undef swap32
#undef swap16
#undef swap_block

// ----- Stolen from evio.h -----------
#define swap64(x) ( (((x) >> 56) & 0x00000000000000FFL) | \
                         (((x) >> 40) & 0x000000000000FF00L) | \
                         (((x) >> 24) & 0x0000000000FF0000L) | \
                         (((x) >> 8)  & 0x00000000FF000000L) | \
                         (((x) << 8)  & 0x000000FF00000000L) | \
                         (((x) << 24) & 0x0000FF0000000000L) | \
                         (((x) << 40) & 0x00FF000000000000L) | \
                         (((x) << 56) & 0xFF00000000000000L) )

#define swap32(x) ( (((x) >> 24) & 0x000000FF) | \
                         (((x) >> 8)  & 0x0000FF00) | \
                         (((x) << 8)  & 0x00FF0000) | \
                         (((x) << 24) & 0xFF000000) )

#define swap16(x) ( (((x) >> 8) & 0x00FF) | \
                         (((x) << 8) & 0xFF00) )
//---------------------------------------


uint32_t swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);


//---------------------------------
// swap_block
//---------------------------------
inline void swap_block(uint16_t *inbuff, uint16_t len, uint16_t *outbuff)
{
	for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
		*outbuff = swap16(*inbuff);
	}
}

//---------------------------------
// swap_block
//---------------------------------
inline void swap_block(uint32_t *inbuff, uint32_t len, uint32_t *outbuff)
{
	for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
		*outbuff = swap32(*inbuff);
	}
}

//---------------------------------
// swap_block
//---------------------------------
inline void swap_block(uint64_t *inbuff, uint64_t len, uint64_t *outbuff)
{
	for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
		*outbuff = swap64(*inbuff);
	}
}

