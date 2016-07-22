// $Id$
//
//    File: DBORptrs.h
// Created: Tue Apr 26 14:52:22 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

// This is used to keep a list of pointers to various BOR config
// objects. These are kept globally so all events get copies of
// of the pointers. The DBORptrs objects are kept in the
// JEventSource_EVIOpp object.

#ifndef _DBORptrs_
#define _DBORptrs_

#include <DAQ/Df250BORConfig.h>
#include <DAQ/Df125BORConfig.h>
#include <DAQ/DF1TDCBORConfig.h>
#include <DAQ/DCAEN1290TDCBORConfig.h>

#define MyBORTypes(X) \
	X(Df250BORConfig) \
	X(Df125BORConfig) \
	X(DF1TDCBORConfig) \
	X(DCAEN1290TDCBORConfig)


#include <JANA/jerror.h>

#include <DAQ/LinkAssociations.h>

class JEventSource_EVIOpp;

class DBORptrs{
	public:
		DBORptrs(void){}
		virtual ~DBORptrs(){ Delete(); }
		

		// For each type defined in "MyTypes" above, define a vector of
		// pointers to it with a name made by prepending a "v" to the classname
		// The following expands to things like e.g.
		//
		//       vector<Df250BORConfig*> vDf250BORConfig;
		//
		#define makevector(A) vector<A*>  v##A;
		MyBORTypes(makevector)
		#undef makevector
		
		// Method to delete all objects in all vectors. This should
		// usually only be called from the destructor
		#define deletevector(A)     for(auto p : v##A) delete p;
		#define clearvectors(A)     v##A.clear();
		void Delete(void){
			MyBORTypes(deletevector)
			MyBORTypes(clearvectors)
		}
		#undef deletevector
		#undef clearvectors


		// Sort all vectors by rocid then slot (use sort from LinkAssociations.h)
		#define sortvector(A) if( v##A.size()>1 ) sort(v##A.begin(), v##A.end(), SortByModule<A>);
		void Sort(void){ MyBORTypes(sortvector) }
		#undef sortvector
};

#endif // _DBORptrs_

