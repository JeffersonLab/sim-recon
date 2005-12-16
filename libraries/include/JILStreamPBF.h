// $Id$

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <JILStream.h>

#ifndef _JILSTREAMPBF_H_
#define _JILSTREAMPBF_H_


class JILStreamPBF: public JILStream {
   public:
	
	JILStreamPBF(string filename="", string mode="w"):JILStream(filename,mode)
	{
		pbfout = NULL;
		pbfin = NULL;
		buff = NULL;
		buff_start = NULL;
		max_buff_size = 1024*250; // default is 250kB
	
		// Open the output or input file		
		if(iotype == STREAM_OUTPUT){
			// Output to file
			if(filename != ""){
				pbfout = new ofstream(filename.c_str());
			}else{
				cerr<<"A filename MUST be supplied when instantiating a JILStreamPBF!"<<endl;
				exit(-1);
			}
			(*pbfout)<<"JILStreamPBF version=0.1"<<std::endl;
			(*pbfout)<<JILMyDictionary();
		}else{
			// Input from file
			if(filename != ""){
				pbfin = new ifstream(filename.c_str());
			}else{
				cerr<<"A filename MUST be supplied when instantiating a JILStreamPBF!"<<endl;
				exit(-1);
			}
			
			// Read in the XML dictionary
			GetDictionary(pbfin);
		}
	}

	~JILStreamPBF(){
		FreeNamed();
		if(pbfout)delete pbfout;
		if(pbfin)delete pbfin;
	}
	
	//-------------------- Serialization -------------------------
	// Return a token that can be used to refer to the current
	// stream position. For the basic JILStreamASCII, we just use the
	// file position.
	std::streamoff GetStreamPosition(void){
		return (std::streamoff)(buff - buff_start);
	}

	// Atomic types. All atomic types (plus STL types)
	JILStreamPBF& operator<<(short i){*(short*)buff=i; buff+=sizeof(short); return *this;}
   JILStreamPBF& operator<<(int i){*(int*)buff=i; buff+=sizeof(int); return *this;}
   JILStreamPBF& operator<<(long i){*(long*)buff=i; buff+=sizeof(long); return *this;}
   JILStreamPBF& operator<<(unsigned short i){*(unsigned short*)buff=i; buff+=sizeof(unsigned short); return *this;}
   JILStreamPBF& operator<<(unsigned int i){*(unsigned int*)buff=i; buff+=sizeof(unsigned int); return *this;}
   JILStreamPBF& operator<<(unsigned long i){*(unsigned long*)buff=i; buff+=sizeof(unsigned long); return *this;}
   JILStreamPBF& operator<<(float f){*(float*)buff=f; buff+=sizeof(float); return *this;}
   JILStreamPBF& operator<<(double f){*(double*)buff=f; buff+=sizeof(double); return *this;}
   JILStreamPBF& operator<<(std::string s){
		// strings are prefixed by size
		*(unsigned int*)buff=s.size();
		buff+=sizeof(unsigned int);
		memcpy((char*)buff, s.c_str(), s.size());
		buff+=s.size();
		return *this;
	}

	// Catch stream manipulators
	JILStreamPBF& operator<<(JILStreamManipulator_t m){
		unsigned int *size_ptr;
		unsigned int size;
		switch(m){
			case END_NAMED:
				named_depth--;
				if(names.size()>0){
					string s = *names.begin();
					names.pop_front();
					
					// Write the total size of this buffer to the front
					// and then write it to the file
					section_size = (unsigned int)buff - (unsigned int)buff_start;
					*(unsigned int*)buff_start = section_size - sizeof(unsigned int); // don't include size word
					pbfout->write(buff_start, section_size);
					
					// Delete the buffer
					delete buff_start;
					buff = NULL;
					buff_start = NULL;
					
				}else{
					std::cerr<<"END_NAMED manipulator reached without matching name!"<<std::endl;
				}
				break;
			case END_OBJECT:
				type_depth--;
				// Update size of object in buffer for top-level objects
				if(type_depth == 0){
					size_ptr = *object_sizes.begin();
					size = (unsigned int)buff - (unsigned int)size_ptr;
					*size_ptr = size;
					object_sizes.pop_front();
				}
				break;
			case END_VECTOR:
				vector_depth--;
				break;
			case END_LIST:
				list_depth--;
				break;
			case END_ARRAY:
				array_depth--;
				break;
			case END_POINTER:
				pointer_depth--;
				break;
		}
		return *this;
	}
	
	/// The StartObjectWrite() method is called just before the
	/// data members are streamed. Returning true tells the
	/// serializer method to go ahead and stream the members.
	/// Returning false tells it to skip that and just send the
	/// END_OBJECT manipulator right away.
	bool StartObjectWrite(const std::type_info *t, void *ptr){
		type_depth++;

		// Write out object header only for top-level objects
		if(type_depth==1){
			object_sizes.push_front((unsigned int*)buff);
			(*this)<<(unsigned int)0;		// place holder for object size
			(*this)<<JILtypeid2name(t);	// name(type) of object
			(*this)<<tag;						// tag of object
		}

		// Now write out the object as though it were streamed to us
		// as a pointer.
		return StartPointerWrite(t, ptr);
	}

	// Allow user to insert object-like tags
	void StartNamedWrite(const char *name){
		names.push_front(string(name));
		if(names.size() > 1){
			cerr<<__FILE__<<":"<<__LINE__<<" JILStreamPBF does not support nested named sections! This may not work how you want it to ..."<<std::endl;
		}else{
			buff_start = new char[max_buff_size];
			buff = buff_start;
			if(!buff){
				std::cerr<<"Unable to allocate "<<max_buff_size<<" byte buffer!"<<std::endl;
				exit(-1);
			}
			
			// First sizeof(unsigned int) bytes is reserved for buffer size
			unsigned int tmp=0;
			(*this)<<tmp;
			
			// Next item is name of the named section
			string namestr(name);
			(*this)<<namestr;
		}
	}

	/// Called before all items of a vector are streamed
	void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	/// Called before all items of a list are streamed
	void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	/// Called before all items of an array are streamed
	void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	/// Called for all pointers. If it returns "true", then the members of
	/// the object are streamed. Otherwise, they are not.
	bool StartPointerWrite(const std::type_info* &t, void* ptr){
		pointer_depth++;
	
		// CacheObjectPointerWrite() will keep track of pointers for
		// the current named section modifying pos if needed according
		// to the current pointer tracking model set in the base class.
		std::streamoff pos = GetStreamPosition();
		bool write_object = CacheObjectPointerWrite(pos, t, ptr);
		(*this)<<(uint)pos;				// position of this object if writing out
												// position of other buffer location if not
		return write_object;
	}
	
	//-------------------- Deserialization -------------------------

	/// Set the stream position to that referred to by the given reference.
	/// For the basic JILStreamASCII, we just use the file position.
	bool SetStreamPosition(std::streamoff pos){
		buff = (char*)((std::streamoff)buff_start + pos);
		return true;
	}

	JILStreamPBF& operator>>(short &i){i=*(short*)buff; buff+=sizeof(short); return *this;}
   JILStreamPBF& operator>>(int &i){i=*(int*)buff; buff+=sizeof(int); return *this;}
   JILStreamPBF& operator>>(long &i){i=*(long*)buff; buff+=sizeof(long); return *this;}
   JILStreamPBF& operator>>(unsigned short &i){i=*(unsigned short*)buff; buff+=sizeof(unsigned short); return *this;}
   JILStreamPBF& operator>>(unsigned int &i){i=*(unsigned int*)buff; buff+=sizeof(unsigned int); return *this;}
   JILStreamPBF& operator>>(unsigned long &i){i=*(unsigned long*)buff; buff+=sizeof(unsigned long); return *this;}
   JILStreamPBF& operator>>(float &f){f=*(float*)buff; buff+=sizeof(float); return *this;}
   JILStreamPBF& operator>>(double &f){f=*(double*)buff; buff+=sizeof(double); return *this;}
   JILStreamPBF& operator>>(std::string &s){
		// strings are prefixed by size
		unsigned int size = *(unsigned int*)buff;
		buff+=sizeof(unsigned int);
		char *cstr = new char[size+1];
		memcpy(cstr, (char*)buff, size);
		cstr[size] = 0;
		s=cstr;
		delete cstr;
		buff+=size;
		return *this;
	}

	/// This is called just before an object is read in from the
	/// stream. The next item in the stream is the position of the
	/// object in the stream which is used to decide whether or not
	/// to go ahead and read or just set to the pointer to one
	/// already read in. Returning false means the deserializer
	/// routine should NOT try filling from the stream.
	bool StartPointerRead(const type_info* type, void* &ptr){
		// First, read in the stream position of the object. There
		// are 3 possibilities: 1.) The position is the current
		// stream position so the object is actually stored here 
		// and should be read in. 2.) The position is NULL which
		// indicates the pointer value was NULL when the object
		// was written. 3.) The position refers to another place
		// in the event and should therefore be in the cache.
		std::streamoff curr_pos = GetStreamPosition();
		unsigned int ipos;
		(*this)>>ipos;
		std::streamoff pos = (unsigned int)ipos;
		
		// Check if pointer was NULL when object was written
		if(pos == 0){ptr=NULL; return false;}
		
		// Check if object is stored at this point in stream
		if(pos == curr_pos){ return true;}
		
		// Pointer should be in cache. Look for it there
		if(CacheFindObjectPointer(pos, ptr)){return false;}
		
		// Uh-oh, somthin' jus' ain't right!
		std::cerr<<__FILE__<<":"<<__LINE__<<" Object position in stream invalid."<<std::endl;
		std::cerr<<"curr_pos="<<curr_pos<<" pos="<<pos<<std::endl;
		
		return false;
	}

	unsigned int StartVectorRead(const type_info &t){
		unsigned int size;
		(*this)>>size;
		return size;
	}
	
	unsigned int StartListRead(const type_info &t){
		unsigned int size;
		(*this)>>size;
		return size;
	}
	
	unsigned int StartArrayRead(const type_info &t, unsigned int size){
		unsigned int my_size;
		(*this)>>my_size;
		return my_size;
	}

	bool GetNamed(const char *name){
	
		do{
			FreeNamed();
			
			// First sizeof(unsigned int) bytes is buffer size
			section_size=0;
			pbfin->read((char*)&section_size, sizeof(unsigned int));
			if(section_size==0)return false;

			// Create buffer to hold section
			buff_start = new char[section_size+sizeof(unsigned int)];
			buff = buff_start;
			if(!buff){
				std::cerr<<"Unable to allocate "<<section_size<<" byte buffer!"<<std::endl;
				exit(-1);
			}
			
			// Copy section size into front of buffer so that the buffer
			// is identical to when it was written
			*(unsigned int*)buff = section_size;
			buff += sizeof(unsigned int);
		
			// Read in named section into buffer
			pbfin->read(buff, section_size);
			if(pbfin->gcount() != (int)section_size)return false;

			// Next item is name of the named section
			string namestr;
			(*this)>>namestr;
			if(namestr==name)break;
		}while(1);
		
		// Read in the top-level objects
		while(buff < buff_start+section_size){
			// Read in object size,type, and tag
			char *object_start = buff;
			unsigned int size = 0;
			string type="", tag="";
			(*this)>>size;
			(*this)>>type;
			(*this)>>tag;
			JILObjectRecord *rec = JILMakeObject(type.c_str(), this, tag.c_str());
			if(rec){
				// If this happens to be pointing to an object already
				// owned by another JILObjectRecord, tell this one he's
				// not the owner
				if(FindObjectRecord(objects, rec->ptr))rec->am_owner = false;
				objects.push_back(rec);
			}

			// If we couldn't read in the object, then "buff" will not have
			// been updated to point to the next object. Force it to be
			// the correct value here so programs that don't know about
			// this object type can still work.
			buff = object_start + size;
				
			// Record some stat info about the object
			AddToObjectStats(type, tag, rec ? rec->type:NULL, size);
		}
			
		return true;
	}
	
	/// Delete all objects allocated on last call to GetNamed() and
	/// all strings in the "lines" list.
	void FreeNamed(void){
		FreeNamedRecords();
		object_sizes.clear();
		delete buff_start;
		buff = NULL;
		buff_start = NULL;
	}
	private:
		ostream *pbfout;
		istream *pbfin;
		list<string> names;
		char *buff;
		char *buff_start;
		unsigned int max_buff_size;
		unsigned int section_size;
		list<unsigned int*> object_sizes;
		
		JILStreamPBF(){} /// Don't allow calls to default constructor. Force a filename
};

#endif //_JILSTREAMPBF_H_

