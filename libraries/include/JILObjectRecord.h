// $Id$

#include <typeinfo>

#ifndef _JILOBJECTRECORD_
#define _JILOBJECTRECORD_

// This is a place holder for the real function that is placed in the
// *_serializers.cc file. By using a template, it allows the initial
// executable to be linked. When re-compiled with the serializers.h
// file, that takes precendence and this is ignored. The only drawback
// is that this has to use the template type in the argument to avoid
// multiply-defined errors. This means the user could get confusing errors
// if they pass something other than a type_info& in as T.
template<typename T>
const char* JILtypeid2name(const T* t){return t->name();}


// Likewise, this is a template place holder for the actual generated
// routine.
class JILObjectRecord;
class JILStream;
template<typename T>
JILObjectRecord* JILMakeObject(const char* type_name, T *s, string tag=""){return new JILObjectRecord();};

// JILObjectRecord
//
/// This keeps info about an object that has been created
/// in memory by JIL. It is a base class and only subclasses
/// of this are uesful. See JILObjectRecordT.
class JILObjectRecord{
	public:
	
		JILObjectRecord(){
			ptr = NULL;
			am_owner = false;
			type = NULL;
			type_name = "";
			tag="";
		}
		virtual ~JILObjectRecord(){}

		void *ptr;
		bool am_owner;
		const std::type_info* type;
		const char *type_name;
		string tag;
};

/// This class actually creates and deletes the objects. When
/// an object is read in, an JILObjectRecordT<T> is created
/// and stored in the JILStream object records list as a
/// JILObjectRecord. It is necessary to do it this way because
/// we need to delete the object using a pointer of type T
/// but the JILObjectRecord objects only have a void* . In order
/// to store a list of records, they must be stored as the base
/// class.
template <typename T>
class JILObjectRecordT:public JILObjectRecord{
	public:
		JILObjectRecordT(JILStream *s, string tag=""){
			tptr = NULL; // Deserializer will allocate and set tptr
			am_owner = true;
			type= &typeid(T);
			type_name = JILtypeid2name(type);
			this->tag = tag;

			// Fill object with contents from stream
			(*s)>>tptr;

			ptr = (void*)tptr;
		}

		~JILObjectRecordT(){
			if(am_owner)delete tptr;
		}
		
		T* tptr;
};

#endif //_JILOBJECTRECORD_


