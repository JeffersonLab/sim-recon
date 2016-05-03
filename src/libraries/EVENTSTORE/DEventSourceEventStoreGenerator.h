// $Id$
//
//    File: DEventSourceEventStoreGenerator.h
// Created: Sat May  8 13:54:46 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DEventSourceEventStoreGenerator_
#define _DEventSourceEventStoreGenerator_

#include <string>
using namespace std;

#include <JANA/JApplication.h>
#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>
using namespace jana;

class DEventSourceEventStoreGenerator:public JEventSourceGenerator{
	public:
		DEventSourceEventStoreGenerator();
		virtual ~DEventSourceEventStoreGenerator();
		
		const char* Description(void);
		double CheckOpenable(string source);
		JEventSource* MakeJEventSource(string source);
		
	protected:
	
	
	private:

};

#endif // _DEventSourceEventStoreGenerator_

