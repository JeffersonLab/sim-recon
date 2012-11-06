// 
// align_16.h
//
//	Macros to support allocation of variables and structures
//	that are guaranteed to be aligned on a 16-byte boundary.
//	Note that ONLY the first variable in the structure is
//	guaranteed to align on a 16-byte boundary; the alignment
//	of the remaining elements follows the standard packing
//	rules for c/c++ structure members.
//
// Richard Jones, University of Connecticut
// March 5, 2011
//

#ifndef ALIGNED_16_BLOCK

#include <stdint.h>

#define ALIGNED_16_BLOCK(TYPE,NUM,PTR) \
uint32_t padded_space_for_pointer_##PTR \
	[ ( ( ( NUM ) * sizeof(TYPE)) + 18 ) / 4 ]; \
TYPE * const PTR;

#define ALIGNED_16_BLOCK_PTR(TYPE,NUM,PTR) \
reinterpret_cast< TYPE *> ((reinterpret_cast<uintptr_t>(padded_space_for_pointer_##PTR) + 15) & (~ 15))

#define ALIGNED_16_BLOCK_WITH_PTR(TYPE,NUM,PTR) \
uint32_t padded_space_for_pointer_##PTR \
	[ ( ( ( NUM ) * sizeof(TYPE)) + 18 ) / 4 ]; \
TYPE * const PTR = ALIGNED_16_BLOCK_PTR(TYPE,NUM,PTR);

#endif
