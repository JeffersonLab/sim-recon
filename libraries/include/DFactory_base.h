// $Id$

///
/// Hall-D Data Factory
///
/// All data (except that read from the data source) must
/// be derived from the data that was read from the source.
/// One DFactory object should exist for each type of data
/// that is "generated". For example: Clusters in the FCAL
/// are generated from Hits in the FCAL. Asking the factory
/// for a list of FCAL clusters will make the factory ask
/// for a list of FCAL Hits from the FCAL Hits factory which
/// may, in turn, ask for a list of raw data values from yet
/// another factory).
/// A list of DFactories is kept in the DEvent object so
/// all data may be accessed through it.

#ifndef _DFACTORY_BASE_H_
#define _DFACTORY_BASE_H_

#include <vector>
#include <list>
#include <string>
using namespace std;

#include "DEventProcessor.h"

class DEventLoop;
class JILStream;
class JILObjectRecord;
#include "hddm_s.h"

//-----------------------
// class DFactory_base
//-----------------------
class DFactory_base:public DEventProcessor{
	/// This class is used as a base class for all factory classes.
	/// Typically, a factory object will be an instance of the
	/// the DFactory template class, (based on the class type of the
	/// objects it provides). In order for the DEvent object to
	/// treat all factories alike (i.e. keep an array of them) they
	/// must all be derived from a common base class. This is that
	/// base class.
	
	friend class DEvent;

	public:
	
		/// This gets typecast in the template member function
		/// also named "Get()" in the DEventLoop class. Since there
		/// is no base class for vector objects, we give it something
		/// that should at least be the same size i.e. vector<int*>
		virtual vector<void*>& Get()=0;
		
		/// Returns the number of rows.
		virtual const int GetNrows(void)=0;
		
		/// Delete the factory's data depending on the flags
		virtual derror_t Reset(void)=0;
		
		/// Delete the factory's data regardless of the flags
		virtual derror_t HardReset(void)=0;
		
		/// Return the pointer to the class's name which this
		/// factory provides.
		virtual const char* dataClassName(void)=0;
		
		/// Returns the size of the data class on which this factory is based
		virtual int dataClassSize(void)=0;
		
		/// Returns a string object with a nicely formatted ASCII table of the data
		virtual const string toString(void){return string(" <Print method undefined for ")+dataClassName()+string("> ");}

		/// Extracts data from the given hddm_s structure and places it
		/// in objects whose pointers are placed in the passed vector
		/// (Most factories won't need this)
		virtual derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v){return OBJECT_NOT_AVAILABLE;}

		/// The data tag string associated with this factory. Most factories
		/// will not overide this.
		virtual inline const char* Tag(void){return "";}

#ifdef JILIO		
		/// Access method to have the DFactory template class stream
		/// the data objects to the output
		virtual void StreamToOutput(JILStream *jilstream)=0;

		/// Access method to have the DFactory template class read
		/// the data objects from the input
		virtual void StreamFromInput(JILStream *jilstream, list<JILObjectRecord*> &objects, vector<void*> &v)=0;
#endif // JILIO

		/// Used by DEventLoop to give a pointer back to itself
		void SetDEventLoop(DEventLoop *loop){this->eventLoop=loop;}
		
		enum DFactory_Flags_t{
			DFACTORY_NULL		=0x00,
			PERSISTANT			=0x01,
			WRITE_TO_OUTPUT	=0x02,
			NOT_OBJECT_OWNER	=0x04
		};
		
		const DFactory_Flags_t& GetFactoryFlags(void){return flags;}
	
	protected:
		DEventLoop *eventLoop;
		DFactory_Flags_t flags;
		int debug_level;
		int busy;
		string _table;
		string _row;
		int _icol;
		int _columns[100];
		int header_width;

		inline void SetObjectOwner(void){flags = (DFactory_Flags_t)(flags & ~NOT_OBJECT_OWNER);}
		inline void SetNotObjectOwner(void){flags = (DFactory_Flags_t)(flags | NOT_OBJECT_OWNER);}

		// Methods useful in help produce nicely formatted ASCII
		void printheader(const char *header);
		void printnewrow(void);
		void printcol(const char *str);
		template<typename T> void printcol(const char* format, T);
		void printrow(void);
};


//-------------
// printcol
//-------------
template<typename T>
void DFactory_base::printcol(const char *format, T val)
{
	/// Print a formatted value to "str". Used by Print()
	char str[32];
	sprintf(str, format, val);
	_row.replace(_columns[_icol++]-strlen(str), strlen(str), str);
}

#endif // _DFACTORY_BASE_H_
