
/// Base class for all Hall-D in-memory data containers.
///
/// This is essentially a container class for holding
/// lists of objects that actually hold the data. It
/// is designed so that it can hold the data itself or
/// of course, a set of pointers to the data.


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
			DCONTAINER_NULL	=0x00,
			PERSISTANT			=0x01,
			WRITE_TO_OUTPUT	=0x02
		};
		DContainer_Flags_t flags;
		
		DContainer(void** ptr,int my_rowsize, char *my_name){
			/// ptr can either be NULL, or can point to the location
			/// where the array pointer should be kept. This allows
			/// one to access the data through a mechanism external
			/// to the DContainer class. If NULL is passed, the pointer
			/// value is kept internally in a private member.
			container_ptr = ptr==NULL ? &private_ptr:ptr;
			*container_ptr = NULL;
			nrows = 0;
			maxrows = 0;
			rowsize = my_rowsize;
			name = my_name==NULL ? (char*)"":strdup(my_name);
			flags = DCONTAINER_NULL;
		}
		
		inline derror_t Set(int my_rows, void** my_ptr){
			maxrows = nrows = my_rows;
			container_ptr = my_ptr;
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
		
		derror_t Delete(int idx){
			if(idx<0 || idx>=nrows)return NOERROR;
			if(idx>=0 && idx<nrows-1){
				char *ptr = (char*)*container_ptr;
				memmove(ptr+rowsize*idx, ptr+rowsize*(idx+1), rowsize*(nrows-idx-1));
			}
			nrows--;
		}
		
		inline void ResetNrows(void){nrows=0;}
		inline void* first(void){return (void*)*container_ptr;}
		inline void* last(void){return (void*)(((char*)*container_ptr)+rowsize*(nrows-1));}
		inline void* index(int idx){return (void*)(((char*)*container_ptr)+rowsize*idx);}

	private:
		void *private_ptr;
};

#endif // _DCONTAINER_H_
