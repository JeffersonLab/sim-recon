
/// Base class for all Hall-D in-memory data banks.
/// This is essentially a container class for holding
/// lists of objects that actually hold the data. This
/// could have been done using something like the STL
/// vector class. However, we want to recycle the objects
/// without constantly deleting and allocating memory.
///
/// Memory allocation/deallocation can take significant
/// CPU time. A test on a 1.9GHz Xenon showed about 120ns 
/// per allocation/deallocation cycle. For a 5kb = 1250 word event,
/// one can assume an average of 1000 banks being allocated/
/// deallocated per event. This would result in 0.08ms per
/// event. A 200 CPU farm processing data coming in at
/// 200kHz will have about 1ms per event to make a Level-3
/// analysis of the data. This would mean about 8% of the
/// time would be spent in simply allocating/deallocating
/// memory!

#include <string.h>
#include <derror.h>

#ifndef _DBANK_H_
#define _DBANK_H_

class DBank{
	public:
		int nrows;
		int maxrows;
		int rowsize;
		void **bank_ptr;
		char *name;
		
		DBank(void** ptr,int my_rowsize, char *my_name){
			bank_ptr = ptr;
			*bank_ptr = NULL;
			maxrows = 0;
			rowsize = my_rowsize;
			name = strdup(my_name);
		}

		derror_t Grow(void);
		derror_t Grow(int my_maxrows);		
};

#endif // _DBANK_H_
