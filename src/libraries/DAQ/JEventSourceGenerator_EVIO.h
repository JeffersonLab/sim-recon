// $Id$
// $HeadURL$
//
//    File: JEventSourceGenerator_EVIO.h
// Created: Tue May 21 14:05:48 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _JEventSourceGenerator_EVIO_
#define _JEventSourceGenerator_EVIO_

#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>

#include "JEventSource_EVIO.h"

class JEventSourceGenerator_EVIO: public jana::JEventSourceGenerator{
	public:
		JEventSourceGenerator_EVIO(){}
		virtual ~JEventSourceGenerator_EVIO(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSourceGenerator_EVIO";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		jana::JEventSource* MakeJEventSource(string source);
};

#endif // _JEventSourceGenerator_EVIO_

