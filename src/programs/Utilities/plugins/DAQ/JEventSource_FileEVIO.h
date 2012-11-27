// $Id$
//
//    File: JEventSource_FileEVIO.h
// Created: Mon Nov 26 10:48:42 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2  x86_64)
//

#ifndef _JEventSource_FileEVIO_
#define _JEventSource_FileEVIO_


#include "JEventSource_EVIO.h"


//-----------------------------------------------------------------------
/// The JEventSource_FileEVIO class implements a JEventSource capable of reading in
/// EVIO data from an ET system. The heavy lifting is all done in the JEventSource_DAQ
/// class so look there for details. 

class JEventSource_FileEVIO: public JEventSource_EVIO{
	public:
		                    JEventSource_FileEVIO(const char* source_name);
		            virtual ~JEventSource_FileEVIO();
		virtual const char* className(void){return static_className();}
		 static const char* static_className(void){return "JEventSource_FileEVIO";}

		           jerror_t ReadEVIOEvent(void);
};

#endif // _JEventSourceGenerator_DAQ_

