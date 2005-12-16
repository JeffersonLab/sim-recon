// $Id$

#include <typeinfo>

#ifndef _JILOBJECTRECORD_
#define _JILOBJECTRECORD_

// Declare a couple of classes
class JILObjectRecord;
class JILStream;

// There are default, conditionally compiled versions of these
// below , after the definition of class JILObjectRecord.
const char* JILtypeid2name(const std::type_info* t);
JILObjectRecord* JILMakeObject(const char* type_name, JILStream *s, string tag="");


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

#ifndef _JIL_SERIALIZERS_H_
// These are dummy routines that are compiled only when the
// serializer routines are not present to allow linking of
// an initial exectuable.
inline const char* JILtypeid2name(const std::type_info* t){return t->name();}
inline JILObjectRecord* JILMakeObject(const char* type_name, JILStream *s, string tag){return new JILObjectRecord();};
#endif


#endif //_JILOBJECTRECORD_


