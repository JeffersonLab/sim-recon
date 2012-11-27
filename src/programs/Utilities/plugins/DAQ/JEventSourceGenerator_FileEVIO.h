// $Id$
//
//    File: JEventSourceGenerator_FileEVIO.h
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _JEventSourceGenerator_FileEVIO_
#define _JEventSourceGenerator_FileEVIO_

#include <JANA/jerror.h>
#include <JANA/JEventSourceGenerator.h>

#include "JEventSource_FileEVIO.h"

class JEventSourceGenerator_FileEVIO: public jana::JEventSourceGenerator{
	public:
		JEventSourceGenerator_FileEVIO(){}
		virtual ~JEventSourceGenerator_FileEVIO(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSourceGenerator_FileEVIO";}
		
		const char* Description(void);
		double CheckOpenable(string source);
		jana::JEventSource* MakeJEventSource(string source);
};

#endif // _JEventSourceGenerator_FileEVIO_

