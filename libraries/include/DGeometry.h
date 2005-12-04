// $Id$
//
//    File: DGeometry.h
// Created: Thu Aug 11 22:17:29 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DGeometry_
#define _DGeometry_

#include "derror.h"

class DGeometry{
	public:
		DGeometry(){}
		DGeometry(unsigned int run_number);
		virtual ~DGeometry();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DGeometry";}

		inline bool IsInRange(unsigned int run){return run>=min_run_number && run<=max_run_number;}
		
	protected:
		unsigned int min_run_number;	
		unsigned int max_run_number;	
	
	private:

};

#endif // _DGeometry_

