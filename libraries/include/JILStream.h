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
/// <li>bool StartPointerWrite(const type_info* &t, void* ptr)
/// </ul>
///
/// For deserializing objects from a file:
/// <ul>
/// <li>bool GetNamed(const char *name)
/// <li>void FreeNamed(void)
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

	enum JILStreamPointerTracking_t{
		PTR_AUTOMATIC,
		PTR_NONE
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

		if(mode == "w")iotype = STREAM_OUTPUT;
		else if(mode == "r")iotype = STREAM_INPUT;
		else {std::cerr<<"Unknown mode \""<<mode<<"\". Aborting.."<<endl; exit(-1);}

		this->filename = filename;
		pointer_depth = 0;
		type_depth = 0;
		named_depth = 0;
		vector_depth = 0;
		list_depth = 0;
		pointer_tracking = PTR_AUTOMATIC;
	}

	/// Destructor of basic JILStream object.
	virtual ~JILStream(){
		FreeNamed();
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
	
	/// Set the pointer tracking mode
	virtual void SetPointerTracking(JILStreamPointerTracking_t ptmode){
		pointer_tracking = ptmode;
	}
	
	// Return a value that can be used to refer to the current
	// stream position. What this position is relative to is
	// determined by the subclass
	virtual std::streamoff GetStreamPosition(void)=0;
	
	/// Atomic types. All atomic types (plus strings) are explicitly defined.
	virtual JILStream& operator<<(short i)=0;
	virtual JILStream& operator<<(int i)=0;
	virtual JILStream& operator<<(long i)=0;
	virtual JILStream& operator<<(unsigned short i)=0;
	virtual JILStream& operator<<(unsigned int i)=0;
	virtual JILStream& operator<<(unsigned long i)=0;
	virtual JILStream& operator<<(float f)=0;
	virtual JILStream& operator<<(double f)=0;
	virtual JILStream& operator<<(std::string s)=0;

	/// Handle stream manipulators
	virtual JILStream& operator<<(JILStreamManipulator_t m)=0;

	/// The first pass compilation will use this "catch-all" template
	/// since no specific operator<< methods are defined for the custom
	/// classes. Executables using this template shouldn't be run and are
	/// just used to provide the introspection information.
   template<typename T>
   JILStream& operator<<(T t){
		return UnknownOut(&typeid(t));
	}

	/// This is called when an object is about to be written to the stream.
	/// If it returns "true", then the object data is sent, otherwise, it
	/// is not. This (optionally) keeps track of the pointers of the objects
	/// written so that each object is recorded only once in the output stream.
	template <typename T>
	bool WriteObject(T *ptr){
		return StartObjectWrite(&typeid(T), (void*)ptr);
	}
	
	/// Keep and check a cache of pointers for the current named section.
	/// This routine is used when serializing the object.
	bool CacheObjectPointerWrite(std::streamoff &pos, const std::type_info *t, void *ptr)
	{
		// If the pointer_tracking model is not PTR_AUTOMATIC, then
		// just tell the subclass to write out the object
		if(pointer_tracking != PTR_AUTOMATIC)return true;
		
		// Look for pointer in cache
		list<pointer_cache_t>::iterator iter = pointer_cache.begin();
		for(;iter != pointer_cache.end(); iter++){
			if((*iter).ptr == ptr && (*iter).type == t){
				pos = (*iter).pos;
				return false;
			}
		}
		
		AddObjectToCache(pos,t,ptr);
		return true;
	}
	
	/// Keep and check a cache of pointers for the current named section.
	/// This routine is used when deserializing the object.
	bool CacheFindObjectPointer(std::streamoff pos, void* &ptr)
	{
		// pointer tracking model only affects how objects are written
		// so we ignore it here.
		
		// Look for pointer in cache
		list<pointer_cache_t>::iterator iter = pointer_cache.begin();
		for(;iter != pointer_cache.end(); iter++){
			if((*iter).pos == pos){
				ptr = (*iter).ptr;
				return true; // pointer found in cache.
			}
		}
		
		return false; // pointer not in cache.
	}

	/// Add a pointer with corresponding stream position to the cache
	void AddObjectToCache(std::streamoff pos, const std::type_info *t, void *ptr){
		// Pointer wasn't found in cache. Add it and write out object
		pointer_cache_t p;
		p.ptr = ptr;
		p.type = t;
		p.pos = pos;
		pointer_cache.push_back(p);
	}
	
	// The StartObjectWrite() method is called just before the
	// data members are streamed. If pointer_tracking is set to
	// PTR_AUTOMATIC, then this is called only when an object that
	// has not been written to this named section is about to be
	// written. Otherwise, this is always called allowing the subclass
	// to implement a pointer tracking scheme.
	virtual bool StartObjectWrite(const std::type_info *t, void *ptr)=0;

	/// Pointer types can be caught by this template
	/// which can then either do something fancy, or just stream
	/// the object pointed to.
	template<typename T>
	JILStream& operator<<(const T* tptr){
		if(StartPointerWrite(&typeid(T), (void*)tptr)){
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

	/// Open a new named section
	void StartNamed(const char *name){
		pointer_cache.clear();
		StartNamedWrite(name);
	}

	/// Open a new named section in the subclass
	virtual void StartNamedWrite(const char *name)=0;
	
	/// Called before all items of a vector are streamed
	virtual void StartVectorWrite(const std::type_info &t, unsigned int size, unsigned int bytes_per_item)=0;

	/// Called before all items of a list are streamed
	virtual void StartListWrite(const std::type_info &t, unsigned int size, unsigned int bytes_per_item)=0;

	/// Called before all items of an array are streamed
	virtual void StartArrayWrite(const std::type_info &t, unsigned int size, unsigned int bytes_per_item)=0;
	
	/// Called for all pointers. If it returns "true", then the thing being pointed
	/// to is streamed. Otherwise, it is ignored and NO corresponding
	/// END_POINTER manipulator will be streamed.
	virtual bool StartPointerWrite(const std::type_info* t, void* ptr)=0;

	/// Handle unknown data types. Generally, if this gets called
	/// there's a problem since the user is trying to serialize
	/// an object that has no serializer routine.
	virtual JILStream& UnknownOut(const std::type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}
	
	//-------------------- Deserialization -------------------------

	// Set the stream position to that referred to by the given reference.
	// Like GetStreamPosition() the exact meaning of this values is
	// determined by the subclass
	virtual bool SetStreamPosition(std::streamoff pos)=0;
		
	virtual JILStream& operator>>(short &i)=0;
	virtual JILStream& operator>>(int &i)=0;
	virtual JILStream& operator>>(long &i)=0;
	virtual JILStream& operator>>(unsigned short &i)=0;
	virtual JILStream& operator>>(unsigned int &i)=0;
	virtual JILStream& operator>>(unsigned long &i)=0;
	virtual JILStream& operator>>(float &f)=0;
	virtual JILStream& operator>>(double &f)=0;
	virtual JILStream& operator>>(std::string &s)=0;

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
	/// to StartPointerRead (which can be overridden by the subclass).
	template <typename T>
	bool StartPointerReadT(T* &t){
		std::streamoff pos = GetStreamPosition();
		if(!StartPointerRead(&typeid(T), (void*&)t))return false;		

		if(t == NULL) t = new T();
		AddObjectToCache(pos, &typeid(T), (void*)t);

		return true;
	}
	
	/// returning false means we have handled  this object and
	/// the deserializer routine should NOT try filling from the
	/// stream.
	virtual bool StartPointerRead(const type_info* type, void* &ptr)=0;
	
	/// Get vector info from stream
	virtual unsigned int StartVectorRead(const type_info &t)=0;
	
	/// Get list info from stream
	virtual unsigned int StartListRead(const type_info &t)=0;

	/// Template method to read in an array of items
   template<typename T>
   JILStream& ReadArray(T* tpr, unsigned int size){
		unsigned int mysize=StartArrayRead(typeid(T), size);
		if(mysize<size)size = mysize; // only read in what is in file
		for(unsigned int i=0; i<size; i++, tpr++)(*this)>>*tpr;
		return *this;
	}

	/// Get array info from stream
	virtual unsigned int StartArrayRead(const type_info &t, unsigned int size)=0;

	virtual JILStream& UnknownIn(const type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}

	/// When deserializing, the caller doesn't really know what objects are in the
	/// file so they call this to read a "map" into memory. In the simple case,
	/// the whole named section is read in and then parsed to create the objects
	/// right away. A more sophisticated method might read in only a header
	/// for the section and save the I/O calls for when the objects are actually
	/// requested. This returns true on success and false on failure. Failure
	/// can mean there are simply no more of the requested named sections in the
	/// file.
	virtual bool GetNamed(const char *name)=0;
	
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
	/// JILObjectRecord list. This may not be needed by all subclasses and
	/// so a default empty routine is provided.
	virtual void GetObjectsFromSource(const type_info *t, string tag, list<JILObjectRecord*> *objects_ptr){}
	
	/// Delete all objects allocated on last call to GetNamed()
	void FreeNamedRecords(void){
		DeleteObjectRecords(objects);
		object_stats.clear();
		pointer_cache.clear();
	}
	
	/// This should be called by the subclass at the begining of
	/// GetNamed() to free any memory from the last named section
	/// that was read in. If the subclass implements this method,
	/// it needs to call FreeNamedRecords().
	virtual void FreeNamed(void){
		FreeNamedRecords();
	}
	
	/// Delete all objects in the given list. The list is passed in
	/// so lists that have been "adopted" can be easily freed.
	void DeleteObjectRecords(list<JILObjectRecord*> &objects){
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

	/// Search the object records list for one with the given pointer.
	/// if found, return it. Otherwise, return NULL
	JILObjectRecord* FindObjectRecord(list<JILObjectRecord*> &objects, void* ptr){
		list<JILObjectRecord*>::iterator iter;
		for(iter=objects.begin(); iter!=objects.end(); iter++){
			if((*iter)->ptr == ptr)return *iter;
		}
		return NULL;
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
		JILStreamPointerTracking_t pointer_tracking;
		
		typedef struct{
			string name;
			string tag;
			const std::type_info *type;
			unsigned int num;
			unsigned int total_size;
		}object_stat_t;
		
		typedef struct{
			std::streamoff pos;
			const std::type_info *type;
			void *ptr;
		}pointer_cache_t;

		list<JILObjectRecord*> objects;
		list<object_stat_t> object_stats;
		list<pointer_cache_t> pointer_cache;

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
};

#endif // _JILSTREAM_H_

