// $Id$
//
//    File: DClassDef.h
// Created: Tue Oct 16 04:22:52 EDT 2012
// Creator: davidl (on Darwin harriet.local 11.4.2 i386)
//

#ifndef _DClassDef_
#define _DClassDef_

#include <map>
#include <set>
#include <string>
using namespace std;


class DClassDef{
	public:
		DClassDef();
		virtual ~DClassDef();
		
		string name;
		map<string, string> members;
		set<string> include_types;
		unsigned int depth;
		
	protected:
	
	
	private:

};

#endif // _DClassDef_

