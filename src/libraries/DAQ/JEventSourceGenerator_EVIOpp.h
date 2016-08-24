// $Id$
//
//    File: JEventSourceGenerator_EVIOpp.h
// Created: Tue Mar 29 08:14:42 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _JEventSourceGenerator_EVIOpp_
#define _JEventSourceGenerator_EVIOpp_

#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>

#include <DAQ/HDEVIO.h>

#include "JEventSource_EVIOpp.h"

class JEventSourceGenerator_EVIOpp: public jana::JEventSourceGenerator{
	public:
		JEventSourceGenerator_EVIOpp(){}
		virtual ~JEventSourceGenerator_EVIOpp(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSourceGenerator_EVIOpp";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		jana::JEventSource* MakeJEventSource(string source);
};

#endif // _JEventSourceGenerator_EVIOpp_

