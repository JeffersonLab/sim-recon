
/// Base class for all Hall-D in-memory data containers.
/// This is essentially a container class for holding
/// lists of objects that actually hold the data. This
/// could have been done using something like the STL
/// vector class. However, we want to recycle the objects
/// without constantly deleting and allocating memory.
///
/// Memory allocation/deallocation can take significant
/// CPU time. A test on a 1.9GHz Xenon showed about 120ns 
/// per allocation/deallocation cycle. For a 5kb = 1250 word event,
/// one can assume an average of 1000 containers being allocated/
/// deallocated per event. This would result in 0.08ms per
/// event. A 200 CPU farm processing data coming in at
/// 200kHz will have about 1ms per event to make a Level-3
/// analysis of the data. This would mean about 8% of the
/// time would be spent in simply allocating/deallocating
/// memory!

#include <iostream>
#include <iomanip>
using namespace std;

#include <string.h>
#include <stdlib.h>
#include "derror.h"

#ifndef _DCONTAINER_H_
#define _DCONTAINER_H_

class DContainer{
	public:
		int nrows;
		int maxrows;
		int rowsize;
		void **container_ptr;
		char *name;
		
		enum DContainer_Flags_t{
			DCONTAINER_NULL			=0x00,
			PERSISTANT			=0x01,
			WRITE_TO_OUTPUT	=0x02
		};
		DContainer_Flags_t flags;
		
		DContainer(void** ptr,int my_rowsize, char *my_name){
			container_ptr = ptr;
			*container_ptr = NULL;
			maxrows = 0;
			rowsize = my_rowsize;
			name = strdup(my_name);
			flags = DCONTAINER_NULL;
		}

		inline derror_t Grow(void){
			if(nrows>maxrows)Grow(nrows);
		}
		
		inline derror_t Grow(int my_maxrows)
		{
			if(my_maxrows<=maxrows)return NOERROR;
			maxrows = my_maxrows+10;
			unsigned long size = maxrows*rowsize;
			if(*container_ptr == NULL){
				*container_ptr = malloc(size);
			}else{
				*container_ptr = realloc(*container_ptr, size);
			}
			return NOERROR;
		}
		
		inline void* Add(void){return Add(1);}
		
		void* Add(int N){
			nrows+=N;
			Grow();
			char *ptr = (char*)*container_ptr;
			ptr += rowsize*(nrows-N);
			return (void*)ptr;
		}
};

#endif // _DCONTAINER_H_
