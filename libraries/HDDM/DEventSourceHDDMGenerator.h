// $Id$
//
//    File: DEventSourceHDDMGenerator.h
// Created: Sat Jul  1 19:23:54 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#ifndef _DEventSourceHDDMGenerator_
#define _DEventSourceHDDMGenerator_

#include <JANA/JEventSourceGenerator.h>
using namespace jana;

class DEventSourceHDDMGenerator:public JEventSourceGenerator{
	public:
		DEventSourceHDDMGenerator(){}
		~DEventSourceHDDMGenerator(){}
		const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSourceHDDMGenerator";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		JEventSource* MakeJEventSource(string source);
		
};

#endif // _DEventSourceHDDMGenerator_

