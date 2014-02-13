// $Id$
//
//    File: JFactoryGenerator_ThreadCancelHandler.h
// Created: Wed Dec 19 15:52:58 EST 2012
// Creator: davidl (on Linux ifarm1101 2.6.18-274.3.1.el5 x86_64)
//


// The entire purpose of this class is to replace the
// signal handler for HUP signals installed by JANA
// with out own. It does NOT actually generate any
// factories! See the comments at the top of MyProcessor.cc
// for details.

#ifndef _JFactoryGenerator_ThreadCancelHandler_
#define _JFactoryGenerator_ThreadCancelHandler_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

extern void mcsmear_thread_HUP_sighandler(int sig);

class JFactoryGenerator_ThreadCancelHandler: public jana::JFactoryGenerator{
	public:
		JFactoryGenerator_ThreadCancelHandler(){}
		virtual ~JFactoryGenerator_ThreadCancelHandler(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JFactoryGenerator_ThreadCancelHandler";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop){
		
			jout<<"Installing special signal handler for mcsmear..."<<endl;
			signal(SIGHUP, mcsmear_thread_HUP_sighandler);
		
			return NOERROR;
		}

};

#endif // _JFactoryGenerator_ThreadCancelHandler_

