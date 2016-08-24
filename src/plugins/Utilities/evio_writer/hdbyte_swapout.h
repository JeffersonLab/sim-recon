
#include <stdint.h>

void swap_block_out(uint16_t *inbuff, uint16_t len, uint16_t *outbuff);
void swap_block_out(uint32_t *inbuff, uint32_t len, uint32_t *outbuff);
void swap_block_out(uint64_t *inbuff, uint64_t len, uint64_t *outbuff);

uint32_t swap_bank_out(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_tagsegment_out(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_segment_out(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);

