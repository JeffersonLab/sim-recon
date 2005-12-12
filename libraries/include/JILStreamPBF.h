// $Id$

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifndef _JILSTREAMPBF_H_
#define _JILSTREAMPBF_H_

#include <JILStream.h>

class JILStreamPBF: public JILStream {
   public:
	
	JILStreamPBF(string filename="", string mode="w"){

		buff = NULL;
		buff_start = NULL;
		max_buff_size = 1024*250; // default is 250kB
	
		this->filename = filename;
		if(mode == "r")JILStreamInit(STREAM_INPUT);
		else if(mode == "w")JILStreamInit(STREAM_OUTPUT);
		else{
			cerr<<"Unknown mode for JILStreamPBF \""<<mode<<"\"  should be \"r\" or \"w\""<<endl;
			exit(-1);
		}
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
		//ClearNamed();
		if(iotype == STREAM_OUTPUT){
			delete pbfout;
		}else{
			delete pbfin;
		}
	};
	
	//-------------------- Serialization -------------------------

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
	JILStreamPBF& operator<<(const std::type_info *t){
		// This is called when an object is about to be written out
		// We write out a place holder for the object size first, then
		// the object type and tags as strings. The actual object
		// size is updated when an END_OBJECT is sent.
		object_sizes.push_front((unsigned int*)buff);
		return (*this)<<(unsigned int)0<<JILtypeid2name(t)<<tag;
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
				// Update size of object in buffer
				size_ptr = *object_sizes.begin();
				size = (unsigned int)buff - (unsigned int)size_ptr;
				*size_ptr = size;
				object_sizes.pop_front();
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
	
	JILStream& UnknownOut(const type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}
	
	JILStream& UnknownIn(const type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
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

	// Called before all items of a vector are streamed
	void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	// Called before all items of a list are streamed
	void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	// Called before all items of an array are streamed
	void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*this)<<size;
	}

	// Called for all pointers. If it returns "true", then the thing being pointed
	// to is streamed. Otherwise, it is ignored and NO corresponding
	// END_POINTER manipulator will be streamed.
	bool StartPointerWrite(const type_info &t, void* ptr){
		(*this)<<(unsigned int)ptr;
		return ptr!=NULL;
	}
	
	//-------------------- Deserialization -------------------------

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

	bool GetNamed(const char *name){
	
		do{
			ClearNamed();
			
			// First sizeof(unsigned int) bytes is buffer size
			section_size=0;
			pbfin->read((char*)&section_size, sizeof(unsigned int));
			if(section_size==0)return false;

			// Create buffer to hold section
			buff_start = new char[section_size];
			buff = buff_start;
			if(!buff){
				std::cerr<<"Unable to allocate "<<section_size<<" byte buffer!"<<std::endl;
				exit(-1);
			}
		
			// Read in named section into buffer
			pbfin->read(buff_start, section_size);
			if(pbfin->gcount() != (int)section_size)return false;

			// Next item is name of the named section
			string namestr;
			(*this)>>namestr;
			if(namestr==name)break;
		}while(1);
		
		// Read in the objects
		while(buff < buff_start+section_size){
			// Read in object size,type, and tag
			char *object_start = buff;
			unsigned int size = 0;
			string type="", tag="";
			(*this)>>size;
			(*this)>>type;
			//(*this)>>tag;
			JILObjectRecord *rec = JILMakeObject(type.c_str(), this, tag.c_str());
			if(rec)objects.push_back(rec);
			
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
	
	void ClearNamed(void){
		DeleteObjectRecords(objects);
		object_stats.clear();
		object_sizes.clear();
		delete buff_start;
		buff = NULL;
		buff_start = NULL;
	}

	bool GetPointerFromStream(const type_info* type, void* &ptr){
		// returning false means we have handled  this object and
		// the deserializer routine should NOT try filling from the
		// stream.
		unsigned int my_ptr;
		(*this)>>my_ptr;
		return my_ptr!=0;
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

