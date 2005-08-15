// $Id$
//
//    File: DParameter.h
// Created: Fri Aug 12 12:36:42 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DParameter_
#define _DParameter_

#include <string>
#include <stdlib.h>
using namespace std;

#include "derror.h"

class DParameter{
	public:
		DParameter(string my_key, string my_value);
		virtual ~DParameter();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DParameter";}
		
		inline void SetValue(string val){value = val;}
		inline const string& GetKey(void){return key;}					///< Return key as a STL string
		inline const string& GetValue(void){return value;}				///< Return value as a STL string
		inline float f(void){return (float)atof(value.c_str());}		///< Return value as a float
		inline double d(void){return (double)atof(value.c_str());}	///< Return value as a double
		inline int i(void){return (int)atoi(value.c_str());}			///< Return value as an int
		bool isdefault;
		
	protected:
		string key;
		string value;
	
	private:

};

#endif // _DParameter_

