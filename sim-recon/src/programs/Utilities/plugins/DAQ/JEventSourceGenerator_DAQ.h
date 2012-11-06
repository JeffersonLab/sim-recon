// $Id$
//
//    File: JEventSourceGenerator_DAQ.h
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _JEventSourceGenerator_DAQ_
#define _JEventSourceGenerator_DAQ_

#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>

#include "JEventSource_DAQ.h"

class JEventSourceGenerator_DAQ: public jana::JEventSourceGenerator{
	public:
		JEventSourceGenerator_DAQ(){}
		virtual ~JEventSourceGenerator_DAQ(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSourceGenerator_DAQ";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		jana::JEventSource* MakeJEventSource(string source);
};

#endif // _JEventSourceGenerator_DAQ_

