// $Id$
//
//    File: JEventSourceGenerator_ETEVIO.h
// Created: Mon Nov 26 10:56:40 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _JEventSourceGenerator_ETEVIO_
#define _JEventSourceGenerator_ETEVIO_

#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>

#include "JEventSource_ETEVIO.h"

class JEventSourceGenerator_ETEVIO: public jana::JEventSourceGenerator{
	public:
		JEventSourceGenerator_ETEVIO(){}
		virtual ~JEventSourceGenerator_ETEVIO(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSourceGenerator_ETEVIO";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		jana::JEventSource* MakeJEventSource(string source);
};

#endif // _JEventSourceGenerator_DAQ_

