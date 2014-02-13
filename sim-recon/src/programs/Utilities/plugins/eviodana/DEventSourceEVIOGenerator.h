// $Id$
//
//    File: DEventSourceEVIOGenerator.h
// Created: Sat May  8 13:54:46 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DEventSourceEVIOGenerator_
#define _DEventSourceEVIOGenerator_

#include <string>
using namespace std;

#include <JANA/JApplication.h>
#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>
using namespace jana;

class DEventSourceEVIOGenerator:public JEventSourceGenerator{
	public:
		DEventSourceEVIOGenerator();
		virtual ~DEventSourceEVIOGenerator();
		
		const char* Description(void);
		double CheckOpenable(string source);
		JEventSource* MakeJEventSource(string source);
		
	protected:
	
	
	private:

};

#endif // _DEventSourceEVIOGenerator_

