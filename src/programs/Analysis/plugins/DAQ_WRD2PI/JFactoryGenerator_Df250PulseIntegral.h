// $Id$
//
//    File: JFactoryGenerator_Df250PulseIntegral.h
// Created: Thu Feb 13 12:49:12 EST 2014
// Creator: dalton (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JFactoryGenerator_Df250PulseIntegral_
#define _JFactoryGenerator_Df250PulseIntegral_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "Df250PulseIntegral_factory.h"

class JFactoryGenerator_Df250PulseIntegral: public jana::JFactoryGenerator{
	public:
		JFactoryGenerator_Df250PulseIntegral(){}
		virtual ~JFactoryGenerator_Df250PulseIntegral(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JFactoryGenerator_Df250PulseIntegral";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop){
			loop->AddFactory(new Df250PulseIntegral_factory());
			return NOERROR;
		}

};

#endif // _JFactoryGenerator_Df250PulseIntegral_

