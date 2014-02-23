// $Id$
//
//    File: JFactoryGenerator_DTranslationTable.h
// Created: Thu Jun 27 15:33:38 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _JFactoryGenerator_DTranslationTable_
#define _JFactoryGenerator_DTranslationTable_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DTranslationTable_factory.h"

class JFactoryGenerator_DTranslationTable: public jana::JFactoryGenerator{
	public:
		JFactoryGenerator_DTranslationTable(){}
		virtual ~JFactoryGenerator_DTranslationTable(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JFactoryGenerator_DTranslationTable";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop){
			loop->AddFactory(new DTranslationTable_factory());
			return NOERROR;
		}

};

#endif // _JFactoryGenerator_DTranslationTable_

