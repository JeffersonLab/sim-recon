// $Id$
//
//    File: DException.h
// Created: Wed Sep 14 07:42:51 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DException_
#define _DException_

#include <string>
using namespace std;


class DException{
	public:
		DException(const char* filename, int line, string message);
		virtual ~DException();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DException";}
		
		const char* filename;
		int line;
		string message;
		
	protected:
	
	
	private:

};

#endif // _DException_

