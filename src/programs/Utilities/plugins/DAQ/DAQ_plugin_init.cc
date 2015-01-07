// $Id$
//
//    File: DAQ_plugin_init.cc
// Created: Tue Jan  6 09:23:34 EST 2015
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <JANA/JApplication.h>
#include <DAQ/JEventSourceGenerator_EVIO.h>
#include <DAQ/JFactoryGenerator_DAQ.h>
using namespace jana;


// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddEventSourceGenerator(new JEventSourceGenerator_EVIO());
	app->AddFactoryGenerator(new JFactoryGenerator_DAQ());
}
} // "C"

