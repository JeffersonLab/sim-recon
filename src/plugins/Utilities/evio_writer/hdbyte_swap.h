
#include <stdint.h>

void swap_block(uint16_t *inbuff, uint16_t len, uint16_t *outbuff);
void swap_block(uint32_t *inbuff, uint32_t len, uint32_t *outbuff);
void swap_block(uint64_t *inbuff, uint64_t len, uint64_t *outbuff);

uint32_t swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
uint32_t swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);

