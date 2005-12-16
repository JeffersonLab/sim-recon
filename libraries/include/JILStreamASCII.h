// $Id$


#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <JILStream.h>

#ifndef _JILSTREAMASCII_H_
#define _JILSTREAMASCII_H_

class JILStreamASCII:public JILStream{

   public:
	
	/// Constructor that opens basic JILStreamASCII file. 
	JILStreamASCII(string filename="", string mode="w"):JILStream(filename,mode)
	{
		fin = NULL;
		fout = NULL;

		// Open the output or input file		
		if(iotype == STREAM_OUTPUT){
			// Output to file
			if(filename != ""){
				fout = new ofstream(filename.c_str());
				(*fout)<<"JILStreamASCII"<<std::endl;
				(*fout)<<JILMyDictionary();
			}else{
				fout = &std::cout;
			}
		}else{
			// Input from file
			if(filename != ""){
				fin = new ifstream(filename.c_str());
			}else{
				fin = &std::cin;
			}
			char line[256];
			fin->getline(line, 256);
			if(strncmp(line, "JILStreamASCII", strlen("JILStreamASCII"))){
				std::cerr<<"File is not a JILStreamASCII file!!"<<std::endl;
				exit(-1);
			}
			
			// Read in the XML dictionary
			GetDictionary(fin);
		}
	}

	/// Destructor of basic JILStreamASCII object.
	~JILStreamASCII(){
		FreeNamed();
		if(fout!=NULL && fout != &std::cout)delete fout;
		if(fin!=NULL && fin != &std::cin)delete fin;
	}
	
	
	//-------------------- Serialization -------------------------
	/// Pointer tracking not supported here. Warn user
	virtual void SetPointerTracking(JILStreamPointerTracking_t ptmode){
		std::cerr<<"WARNING: pointer tracking not supported in JILStreamASCII!!"<<std::endl;
		pointer_tracking = ptmode;
	}

	// Return a token that can be used to refer to the current
	// stream position. For the basic JILStreamASCII, we just use the
	// file position.
	std::streamoff GetStreamPosition(void){
		return iotype==STREAM_OUTPUT ? fout->tellp():fin->tellg();
	}
	
	/// Atomic types. All atomic types (plus strings) are explicitly defined.
	JILStream& operator<<(short i){(*fout)<<i<<" short"<<std::endl; return *this;}
	JILStream& operator<<(int i){(*fout)<<i<<" int"<<std::endl; return *this;}
	JILStream& operator<<(long i){(*fout)<<i<<" long"<<std::endl; return *this;}
	JILStream& operator<<(unsigned short i){(*fout)<<i<<" unsigned short"<<std::endl; return *this;}
	JILStream& operator<<(unsigned int i){(*fout)<<i<<" unsigned int"<<std::endl; return *this;}
	JILStream& operator<<(unsigned long i){(*fout)<<i<<" unsigned long"<<std::endl; return *this;}
	JILStream& operator<<(float f){(*fout)<<f<<" float"<<std::endl; return *this;}
	JILStream& operator<<(double f){(*fout)<<f<<" double"<<std::endl; return *this;}
	JILStream& operator<<(std::string s){(*fout)<<s<<std::endl; return *this;}

	/// Handle stream manipulators
	JILStream& operator<<(JILStreamManipulator_t m){
		switch(m){
			case END_NAMED:
				named_depth--;
				(*fout)<<"end_named"<<std::endl;
				break;
			case END_OBJECT:
				type_depth--;
				(*fout)<<"end_type"<<std::endl;
				break;
			case END_VECTOR:
				vector_depth--;
				(*fout)<<"end_vector"<<std::endl;
				break;
			case END_LIST:
				list_depth--;
				(*fout)<<"end_list"<<std::endl;
				break;
			case END_ARRAY:
				array_depth--;
				(*fout)<<"end_array"<<std::endl;
				break;
			case END_POINTER:
				pointer_depth--;
				(*fout)<<"end_pointer"<<std::endl;
				break;
		}
		return (*this);
	}
	
	// The StartObjectWrite() method is called just before the
	// data members are streamed. If pointer_tracking is set to
	// PTR_AUTOMATIC, then this is called only when an object that
	// has not been written to this named section is about to be
	// written. Otherwise, this is always called allowing the subclass
	// to implement a pointer tracking scheme.
	bool StartObjectWrite(const std::type_info *t, void *ptr){
		type_depth++;
		(*fout)<<"type="<<JILtypeid2name(t);
		if(tag.size()>0)(*fout)<<" tag="<<tag;
		(*fout)<<std::endl;
		return StartPointerWrite(t, ptr);
	}

	/// Allow user to insert object-like tags
	void StartNamedWrite(const char *name){
		named_depth++;
		(*fout)<<"named:"<<name<<std::endl;
	}
	
	/// Called before all items of a vector are streamed
	void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		vector_depth++;
		(*fout)<<"vector: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}

	/// Called before all items of a list are streamed
	void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		list_depth++;
		(*fout)<<"list: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}

	/// Called before all items of an array are streamed
	void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		array_depth++;
		(*fout)<<"array: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}
	
	/// Called for all pointers. If it returns "true", then the thing being pointed
	/// to is streamed. Otherwise, it is ignored and NO corresponding
	/// END_POINTER manipulator will be streamed.
	bool StartPointerWrite(const type_info* t, void* ptr){
		pointer_depth++;

		// CacheObjectPointerWrite() will keep track of pointers for
		// the current named section modifying pos if needed according
		// to the current pointer tracking model set in the base class.
		(*fout)<<"0x"<<hex<<(unsigned int)ptr<<dec<<" "<<JILtypeid2name(t)<<std::endl;

		return true;
	}
	
	//-------------------- Deserialization -------------------------

	// Set the stream position to that referred to by the given reference.
	// For the basic JILStreamASCII, we just use the file position.
	bool SetStreamPosition(std::streamoff pos){
		switch(iotype){
			case STREAM_OUTPUT: fout->seekp(pos); break;
			case STREAM_INPUT: fin->seekg(pos); break;
		}
		return true;
	}
		
	JILStream& operator>>(short &i){if(!lines.empty()){i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(int &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(long &i){if(!lines.empty()){ i=atol(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(unsigned short &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(unsigned int &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;};
	JILStream& operator>>(unsigned long &i){if(!lines.empty()){ i=atol(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(float &f){if(!lines.empty()){ f=atof(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(double &f){if(!lines.empty()){ f=atof(lines.front().c_str()); lines.pop_front();} return *this;}
	JILStream& operator>>(std::string &s){if(!lines.empty()){ s = lines.front(); lines.pop_front();} return *this;}

	/// returning false means we have handled  this object and
	/// the deserializer routine should NOT try filling from the
	/// stream.
	bool StartPointerRead(const type_info* type, void* &ptr){
		if(lines.empty())return false;
		unsigned int my_ptr = strtol(lines.front().c_str(), NULL, 0);
		lines.pop_front();

		// Check if pointer was NULL when object was written
		return my_ptr != 0;
	}
	
	/// Get vector info from stream
	unsigned int StartVectorRead(const type_info &t){
		if(lines.empty())return 0;
		// next item in stream should be vector info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("vector: size=");
		return atoi(ptr);
	}
	
	/// Get list info from stream
	unsigned int StartListRead(const type_info &t){
		if(lines.empty())return 0;
		// next item in stream should be list info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("list: size=");
		return atoi(ptr);
	}

	/// Get array info from stream
	unsigned int StartArrayRead(const type_info &t, unsigned int size){
		if(lines.empty())return 0;
		// next item in stream should be array info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("array: size=");
		return atoi(ptr);
	}

	/// When deserializing, the caller doesn't really know what objects are in the
	/// file so they call this to read a "map" into memory. In this simple case,
	/// the whole named section is read in and then parsed to create the objects
	/// right away. A more sophisticated deserializer might read in only a header
	/// for the section and save the I/O calls for when the objects are actually
	/// requested. This returns true on success and false on failure.
	bool GetNamed(const char *name){
		
		// Clear anything we currently have in memory.
		FreeNamed();
	
		// loop through lines of input, storing them when we come across the named section
		char line[1024];
		bool in_named = false;
		int depth = 0;
		char named_str[256];
		sprintf(named_str, "named:%s", name);
		do{
			line[0] = 0;
			fin->getline(line, 1024);
			if(in_named==true){
				if(!strcmp(line, "end_named"))depth--;
				if(depth==0)break;
				lines.push_back(string(line));
			}
			if(!strcmp(line, named_str))in_named=true;
			if(!strncmp(line,"named:", strlen("named:")))depth++;
		}while(strlen(line) != 0);
		if(lines.empty())return false;

		// Loop over the top-level objects in this named section. We need to remove
		// all of the "type=XXX" lines since the object deserializer
		// already has that information.
		do{
			if(lines.empty())break;
			string s = lines.front();
			lines.pop_front();
			string::size_type type_pos = s.find("type=",0);
			if(type_pos != string::npos){
				// Remove all type=XXX and end_XXX
				int depth = 1;
				list<string>::iterator iter = lines.begin();
				unsigned int stream_bytes_for_object = 0;
				for(; iter!=lines.end(); iter++){
					stream_bytes_for_object += (*iter).size();
					if((*iter) == "end_type")depth--;
					bool erase_line = false;
					if(!strncmp((*iter).c_str(), "type=", strlen("type="))){
						depth++;
						erase_line = true;
					}
					if(!strncmp((*iter).c_str(), "end_", strlen("end_")))erase_line = true;
					if(!strncmp((*iter).c_str(), "pointer:", strlen("pointer:")))erase_line = true;
					if(erase_line){
						lines.erase(iter);
						iter--;
					}
					if(depth==0)break;
				}
				
				// If a tag is present, extract it
				type_pos += strlen("type=");
				string type_str = s.substr(type_pos, s.size()-type_pos);
				string tag_str = "";
				string::size_type tag_pos = s.find("tag=",type_pos);
				if(tag_pos != string::npos){
					tag_pos += strlen("tag=");
					tag_str = s.substr(tag_pos, s.size()-tag_pos);
					tag_pos -= strlen("tag=")+1;
					type_str = s.substr(0, tag_pos); 
				}

				// At this point all type= and end_XXX should be removed from
				// the "lines" list, including the end_type for this section.
				// Ergo, "lines" should have all of the data for this object
				JILObjectRecord *rec = JILMakeObject(type_str.c_str(), this, tag_str.c_str());
				if(rec){
					// If this happens to be pointing to an object already
					// owned by another JILObjectRecord, tell this one he's
					// not the owner
					if(FindObjectRecord(objects, rec->ptr))rec->am_owner = false;
					objects.push_back(rec);
				}
				
				// Record some stat info about the object
				AddToObjectStats(type_str, tag_str, rec ? rec->type:NULL, stream_bytes_for_object);
			}
		}while(1);
		
		return true;
	}
	
	/// Delete all objects allocated on last call to GetNamed() and
	/// all strings in the "lines" list.
	void FreeNamed(void){
		FreeNamedRecords();
		lines.clear();
	}

	private:
		ostream *fout;
		istream *fin;
		list<string> lines;
};

#endif // _JILSTREAMASCII_H_

