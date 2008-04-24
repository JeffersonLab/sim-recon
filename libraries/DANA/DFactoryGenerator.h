// $Id$
//
//    File: DFactoryGenerator.h
// Created: Mon Jul  3 21:46:40 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#ifndef _DFactoryGenerator_
#define _DFactoryGenerator_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>
using namespace jana;

class DFactoryGenerator: public JFactoryGenerator{
	public:
		DFactoryGenerator();
		virtual ~DFactoryGenerator();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator";}
		
		jerror_t GenerateFactories(JEventLoop*);

	protected:
	
	
	private:

};

#endif // _DFactoryGenerator_

