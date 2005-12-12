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

#include <JILObjectRecord.h>

#ifndef _JILSTREAM_H_
#define _JILSTREAM_H_

// First pass compilation will not have this
// defined so we can provide a dummy function
const char* JILMyDictionary(void);
#ifndef _JIL_SERIALIZERS_H_
inline const char* JILMyDictionary(void){return "";}
#endif

/// The JILStream class is the base class used to derive classes
/// capable of reading/writing objects in specialized file formats. The
/// JILStream class itself is fully functional, but the file format
/// it uses is uncompressed ASCII with lots of redundant information
/// which leads to very large files.
///
/// When writing a subclass of JILStream, one typically needs to
/// override the default virtual methods for the following methods
/// in addition to the methods for converting atomic types. See the
/// descriptions below for more details.
///
/// For serializing to a file:
/// <ul>
/// <li>JILStream& operator<<(const std::type_info *t)
/// <li>JILStream& operator<<(JILStreamManipulator_t m)
/// <li>void StartNamedWrite(const char *name)
/// <li>void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item)
/// <li>void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item)
/// <li>void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item)
/// <li>bool StartPointerWrite(const type_info &t, void* ptr)
/// </ul>
///
/// For deserializing objects from a file:
/// <ul>
/// <li>bool GetNamed(const char *name)
/// <li>void ClearNamed(void)
/// <li>bool GetPointerFromStream(const type_info* type, void* &ptr)
/// <li>unsigned int StartVectorRead(const type_info &t)
/// <li>unsigned int StartListRead(const type_info &t)
/// <li>unsigned int StartArrayRead(const type_info &t, unsigned int size)
/// <li>void GetObjectsFromSource(const type_info *t, string tag, list<JILObjectRecord*> *objects_ptr)
/// </ul>
/// <hr>

class JILStream{

   public:
	
	enum JILStreamIOType_t{
		STREAM_INPUT,
		STREAM_OUTPUT
	};

	/// These "manipulators" are sent into the stream to flag the
	/// end of the various types of items.
	enum JILStreamManipulator_t{
		END_NAMED,
		END_OBJECT,
		END_VECTOR,
		END_LIST,
		END_ARRAY,
		END_POINTER
	};
	
	/// Constructor that opens basic JILStream file. This should NOT
	/// be called when instantiating a subclass object. See JILStreamInit
	/// below.
	JILStream(string filename="", string mode="w"){

		this->filename = filename;
		if(mode == "r")JILStreamInit(STREAM_INPUT);
		else if(mode == "w")JILStreamInit(STREAM_OUTPUT);
		else{
			cerr<<"Unknown mode for JILStream \""<<mode<<"\"  should be \"r\" or \"w\""<<endl;
			exit(-1);
		}
		
		if(iotype == STREAM_OUTPUT){
			// Output to file
			if(filename != ""){
				fout = new ofstream(filename.c_str());
				(*fout)<<"BasicJILStream"<<std::endl;
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
			if(strncmp(line, "BasicJILStream", strlen("BasicJILStream"))){
				std::cerr<<"File is not a Basic JILStream file!!"<<std::endl;
				exit(-1);
			}
			
			// Read in the XML dictionary
			GetDictionary(fin);
		}
	}

	/// This initializes the JILStream object and should
	/// be called by all subclasses.
	void JILStreamInit(JILStreamIOType_t iotype){
		this->iotype = iotype;
		pointer_depth = 0;
		type_depth = 0;
		named_depth = 0;
		vector_depth = 0;
		list_depth = 0;
	}

	/// Destructor of basic JILStream object.
	virtual ~JILStream(){
		ClearNamed();
		if(!fout && fout != &std::cout)delete fout;
		if(!fin && fin != &std::cin)delete fin;
	}
	
	/// Read in dictionary information from the stream. This will
	/// read in the XML formatted dictionary information from the
	/// stream. It is assumed that most subclasses will use this
	/// routine, but it is made virtual in case some need to implement
	/// dictionary info in something other than the default XML.
	virtual void GetDictionary(istream *instream){
		// For now, we just read in until we find the end tag
		char line[1024];
		string endtag="</JILDictionary>";
		do{
			line[0] = 0;
			instream->getline(line,1024);
			if(endtag == line)break;
		}while(instream->gcount()!=0);
	}
	
	//-------------------- Serialization -------------------------

	/// Set the current tag to be written out with subsequent objects
	void SetTag(string tag){
		this->tag = tag;
	}
	
	/// Atomic types. All atomic types (plus strings) are explicitly defined.
	virtual JILStream& operator<<(short i){(*fout)<<i<<" short"<<std::endl; return *this;}
	virtual JILStream& operator<<(int i){(*fout)<<i<<" int"<<std::endl; return *this;}
	virtual JILStream& operator<<(long i){(*fout)<<i<<" long"<<std::endl; return *this;}
	virtual JILStream& operator<<(unsigned short i){(*fout)<<i<<" unsigned short"<<std::endl; return *this;}
	virtual JILStream& operator<<(unsigned int i){(*fout)<<i<<" unsigned int"<<std::endl; return *this;}
	virtual JILStream& operator<<(unsigned long i){(*fout)<<i<<" unsigned long"<<std::endl; return *this;}
	virtual JILStream& operator<<(float f){(*fout)<<f<<" float"<<std::endl; return *this;}
	virtual JILStream& operator<<(double f){(*fout)<<f<<" double"<<std::endl; return *this;}
	virtual JILStream& operator<<(std::string s){(*fout)<<s<<std::endl; return *this;}
	virtual JILStream& operator<<(const std::type_info *t){
		(*fout)<<"type="<<JILtypeid2name(t);
		if(tag.size()>0)(*fout)<<" tag="<<tag;
		(*fout)<<std::endl;
		type_depth++;
		return *this;
	}

	/// Handle stream manipulators
	virtual JILStream& operator<<(JILStreamManipulator_t m){
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

	/// The first pass compilation will use this "catch-all" template
	/// since no specific operator<< methods are defined for the custom
	/// classes. Executables using this template shouldn't be run and are
	/// just used to provide the introspection information.
   template<typename T>
   JILStream& operator<<(T t){
		return UnknownOut(&typeid(t));
	}

	/// Pointer types can be caught by this template
	/// which can then either do something fancy, or just stream
	/// the object pointed to.
	template<typename T>
	JILStream& operator<<(const T* tptr){
		if(StartPointerWrite(typeid(T), (void*)tptr)){
			if(tptr)(*this)<<*tptr;
			(*this)<<END_POINTER;
		}
		return (*this);
	}
	
	/// Pass this over to the const version
   template<typename T>
   JILStream& operator<<(T* tptr){
		(*this)<<(const T*)tptr;
		return (*this);
   }

	/// STL container classes are handled by templates. They basically
	/// convert the call to the template method into a call to a
	/// non-templated virtual method. This means the template only
	/// has to be here and the subclass does not need to deal with
	/// templates at all. We first call the StartVectorWrite() method so we
	/// can tell it(the subclass) what type of data, how many items, and how big
	/// each item is before streaming the individual items. After all
	/// items are streamed, we send an END_VECTOR manipulator, which may
	/// be a little redundant, but is easily enough ignored if it is not
	/// needed by the subclass.
   template<typename T>
   JILStream& operator<<(std::vector<T> v){
		StartVectorWrite(typeid(T), v.size(), sizeof(T));
		for(unsigned int i=0; i<v.size(); i++)(*this)<<v[i];
		(*this)<<END_VECTOR;
		return (*this);
	}

   template<typename T>
   JILStream& operator<<(std::list<T> v){
		StartListWrite(typeid(T), v.size(), sizeof(T));
		//for(list<T>::iterator iter=v.begin(); iter!=v.end(); iter++)(*this)<<(*iter);
		(*this)<<END_LIST;
		return (*this);
	}

   template<typename T>
   JILStream& WriteArray(T* tpr, unsigned int size){
		StartArrayWrite(typeid(T), size, sizeof(T));
		for(unsigned int i=0; i<size; i++, tpr++)(*this)<<*tpr;
		(*this)<<END_ARRAY;
		return (*this);
	}

	/// Allow user to insert object-like tags
	virtual void StartNamedWrite(const char *name){
		named_depth++;
		(*fout)<<"named:"<<name<<std::endl;
	}
	
	/// Called before all items of a vector are streamed
	virtual void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		vector_depth++;
		(*fout)<<"vector: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}

	/// Called before all items of a list are streamed
	virtual void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		list_depth++;
		(*fout)<<"list: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}

	/// Called before all items of an array are streamed
	virtual void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		array_depth++;
		(*fout)<<"array: size="<<size<<" bytes_per_item="<<bytes_per_item<<std::endl;
	}
	
	/// Called for all pointers. If it returns "true", then the thing being pointed
	/// to is streamed. Otherwise, it is ignored and NO corresponding
	/// END_POINTER manipulator will be streamed.
	virtual bool StartPointerWrite(const type_info &t, void* ptr){
		pointer_depth++;
		(*fout)<<"0x"<<hex<<(unsigned long)ptr<<dec<<" "<<JILtypeid2name(&t)<<std::endl;
		return true;
	}

	/// Handle unknown data types. Generally, if this gets called
	/// there's a problem since the user is trying to serialize
	/// an object that has no serializer routine.
	virtual JILStream& UnknownOut(const type_info* type){
		(*fout)<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}
	
	//-------------------- Deserialization -------------------------
		
	virtual JILStream& operator>>(short &i){if(!lines.empty()){i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(int &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(long &i){if(!lines.empty()){ i=atol(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(unsigned short &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(unsigned int &i){if(!lines.empty()){ i=atoi(lines.front().c_str()); lines.pop_front();} return *this;};
	virtual JILStream& operator>>(unsigned long &i){if(!lines.empty()){ i=atol(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(float &f){if(!lines.empty()){ f=atof(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(double &f){if(!lines.empty()){ f=atof(lines.front().c_str()); lines.pop_front();} return *this;}
	virtual JILStream& operator>>(std::string &s){if(!lines.empty()){ s = lines.front(); lines.pop_front();} return *this;}

	/// Here we do something similar to the above, but for input streams.
	/// There is an additional trick required here. Since the input
	/// operators take a pointer, we need to catch the stream calls
	/// that use the object directly (not the pointer) and convert them
	/// into calls using the pointer. The second template routine below
	/// will catch all of the undefined objects for the first pass compilation.
   template<typename T>
   JILStream& operator>>(T &t){
		T* tt= &t;
		return (*this)>>tt;
	}
	
   template<typename T>
   JILStream& operator>>(T* &t){
		return UnknownIn(&typeid(t));
	}

	/// Read a vector in from the stream
	template<typename T>
   JILStream& operator>>(std::vector<T> &v){
		unsigned int size=StartVectorRead(typeid(T));
		for(unsigned int i=0;i<size;i++){
			T t;
			(*this)>>t;
			v.push_back(t);
		}
		
		// end_vector should already be removed.
		return *this;
	}

	/// Read a list in from the stream
	template<typename T>
   JILStream& operator>>(std::list<T> &v){
		unsigned int size=StartListRead(typeid(T));
		for(unsigned int i=0;i<size;i++){
			T t;
			(*this)>>t;
			v.push_back(t);
		}
		
		// end_list should already be removed.
		return *this;
	}
	
	/// This template just converts the typed call to an un-typed one
	/// to GetPointerFromStream (which can be overridden by the subclass).
	template <typename T>
	bool GetPointerFromStreamT(T* &t){
		return GetPointerFromStream(&typeid(T), (void*&)t);
	}
	
	/// returning false means we have handled  this object and
	/// the deserializer routine should NOT try filling from the
	/// stream.
	virtual bool GetPointerFromStream(const type_info* type, void* &ptr){
		if(lines.empty())return false;
		unsigned long myptr = strtol(lines.front().c_str(), NULL, 0);
		lines.pop_front();
		return myptr!=0x0;
	}
	
	/// Get vector info from stream
	virtual unsigned int StartVectorRead(const type_info &t){
		if(lines.empty())return 0;
		// next item in stream should be vector info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("vector: size=");
		return atoi(ptr);
	}
	
	/// Get list info from stream
	virtual unsigned int StartListRead(const type_info &t){
		if(lines.empty())return 0;
		// next item in stream should be list info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("list: size=");
		return atoi(ptr);
	}

	/// Template method to read in an array of items
   template<typename T>
   JILStream& ReadArray(T* tpr, unsigned int size){
		unsigned int mysize=StartArrayRead(typeid(T), size);
		if(mysize<size)size = mysize; // only read in what is in file
		for(unsigned int i=0; i<size; i++, tpr++)(*this)>>*tpr;
		return *this;
	}

	/// Get array info from stream
	virtual unsigned int StartArrayRead(const type_info &t, unsigned int size){
		if(lines.empty())return 0;
		// next item in stream should be array info
		const char *ptr = lines.front().c_str();
		lines.pop_front();
		ptr += strlen("array: size=");
		return atoi(ptr);
	}

	virtual JILStream& UnknownIn(const type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}

	/// When deserializing, the caller doesn't really know what objects are in the
	/// file so they call this to read a "map" into memory. In this simple case,
	/// the whole named section is read in and then parsed to create the objects
	/// right away. A more sophisticated deserializer might read in only a header
	/// for the section and save the I/O calls for when the objects are actually
	/// requested. This returns true on success and false on failure.
	virtual bool GetNamed(const char *name){
		
		// Clear anything we currently have in memory.
		ClearNamed();
	
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
			char *ptr = strstr(s.c_str(), "type=");
			if(ptr){
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
				ptr += strlen("type=");
				char *tag_ptr = strstr(ptr,"tag=");
				if(tag_ptr){
					tag_ptr[-1] = 0; // change "space" between object type and tag into null 
					tag_ptr+=strlen("tag=");
				}else{
					tag_ptr = "";
				}

				// At this point all type= and end_XXX should be removed from
				// the "lines" list, including the end_type for this section.
				// Ergo, "lines" should have all of the data for this object
				JILObjectRecord *rec = JILMakeObject(ptr, this, tag_ptr);
				if(rec)objects.push_back(rec);
				
				// Record some stat info about the object
				AddToObjectStats(ptr, tag_ptr, rec ? rec->type:NULL, stream_bytes_for_object);
			}
		}while(1);
		
		return true;
	}
	
	/// Get a list of the objects of type T in the current named section.
	/// The objects list is passed in so lists that have been "adopted" can be
	/// accessed in the same way that the internal list is.
	template <typename T>
	void GetObjects(vector<T*> &v, string tag = "", list<JILObjectRecord*> *objects_ptr=NULL){
		// Look through all object records and copy the pointers of
		// type T into the vector v.
		if(objects_ptr==NULL)objects_ptr = &objects;
		const type_info *t = &typeid(T);
		GetObjectsFromSource(t, tag, objects_ptr);
		list<JILObjectRecord*>::iterator iter;
		for(iter=objects_ptr->begin(); iter!=objects_ptr->end(); iter++){
			if((*iter)->type == t && (*iter)->tag == tag)v.push_back((T*)(*iter)->ptr);
		}
	}
	
	/// Read the objects of the specified type and tag in and add them to the 
	/// JILObjectRecord list. This is not needed for the file format used by
	/// the base class, but is here for subclasses that implement a 
	/// rich header format that can read in objects only when requested.
	virtual void GetObjectsFromSource(const type_info *t, string tag, list<JILObjectRecord*> *objects_ptr){
		return;
	}
	
	/// Delete all objects allocated on last call to GetNamed() and
	/// all strings in the "lines" list. Subclasses should include
	/// a call to DeleteObjectRecords(objects);
	virtual void ClearNamed(void){
		DeleteObjectRecords(objects);
		lines.clear();
		object_stats.clear();
	}
	
	/// Delete all objects in the given list. The list is passed in
	/// so lists that have been "adopted" can be easily freed.
	virtual void DeleteObjectRecords(list<JILObjectRecord*> &objects){
		list<JILObjectRecord*>::iterator iter;
		for(iter=objects.begin(); iter!=objects.end(); iter++)delete *iter;
		objects.clear();
	}
	
	/// Allow object records (and the objects they point to)
	/// to be adopted so that the objects are not deleted
	/// by the JILStream and the stream forgets about them
	void AdoptObjectRecords(list<JILObjectRecord*> &objects)
	{
		objects = this->objects;
		this->objects.clear();
	}

	void PrintObjectHisto(list<JILObjectRecord*> *objects_ptr=NULL){
		if(!objects_ptr)objects_ptr = &objects;
		map<string, int> h;
		list<JILObjectRecord*>::iterator iter;
		for(iter=objects_ptr->begin(); iter!=objects_ptr->end(); iter++){
			h[string((*iter)->type_name)]++;
		}
		
		cout<<endl;
		cout<<" N  Type"<<endl;
		cout<<"--  ----"<<endl;
		map<string, int>::iterator mi;
		for(mi=h.begin(); mi!=h.end(); mi++){
			cout<<mi->second<<"   "<<mi->first<<std::endl;
		}
		cout<<endl;
	}
	
	/// Print statistics about the objects found in the last call to
	/// GetNamed(). This can be useful in cases where PrintObjectHisto()
	/// is not since that requires the class definitions be known to
	/// the current program whereas this does not. If a subclass of
	/// JILStream is used, it must call AddToObjectStats() for each
	// object it encounters in the file in order for this to be useful.
	void PrintObjectStats(void){
			list<object_stat_t>::iterator iter;
			std::cout<<"Object Type               Num.    Total Bytes"<<std::endl;
			std::cout<<"------------------------- ----    -----------"<<std::endl;
			for(iter = object_stats.begin(); iter !=object_stats.end(); iter++){
				string str(79,' ');
				str.replace(0, (*iter).name.size(), (*iter).name);
				str.replace(20, (*iter).tag.size(), (*iter).tag);
				stringstream ss,st;
				ss<<(*iter).num;
				str.replace(30-ss.str().size(), ss.str().size(), ss.str());
				st<<(*iter).total_size;
				str.replace(45-st.str().size(), st.str().size(), st.str());
				
				std::cout<<str<<std::endl;
			}
	}

	protected:
		JILStreamIOType_t iotype;
		string filename;
		int pointer_depth;
		int type_depth;
		int named_depth;
		int vector_depth;
		int list_depth;
		int array_depth;
		string tag;
		
		typedef struct{
			string name;
			string tag;
			const std::type_info *type;
			unsigned int num;
			unsigned int total_size;
		}object_stat_t;
		
		list<JILObjectRecord*> objects;
		list<object_stat_t> object_stats;
		
		/// Used to keep track of objects found in file, even if the
		/// classes aren't known (i.e. compiled into) the current program
		void AddToObjectStats(string name, string tag, const std::type_info *type, unsigned int bytes){
			// First, see if an object_stat_t exists for this object
			list<object_stat_t>::iterator iter;
			for(iter = object_stats.begin(); iter !=object_stats.end(); iter++){
				if((*iter).name == name && (*iter).tag == tag && (*iter).type == type){
					(*iter).num++;
					(*iter).total_size += bytes;
					return;
				}
			}
			
			object_stat_t object_stat;
			object_stat.name=name;
			object_stat.tag=tag;
			object_stat.type=type;
			object_stat.num=1;
			object_stat.total_size=bytes;
			object_stats.push_back(object_stat);
		}
		
	private:
		ostream *fout;
		istream *fin;
		list<string> lines;
};

#endif // _JILSTREAM_H_

