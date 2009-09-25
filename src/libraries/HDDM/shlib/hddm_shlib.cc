// $id
//
// hddm_shlib.cc
//
/// This file provides is used to generated a shared library
/// that can handle a specific version of HDDM. It can be dynamically
/// loaded by an an existing executable using the --so= or
/// --sodir= command line switches.
///
/// THIS CAN BE DANGEROUS!!! The routines here assume a certain format
/// for the C++ classes which may not be consistent with the
/// format compiled into the executable. This is intended as a quick
/// solution to the problem of reading in older simulation files
/// with newer executables where changes to HDDM were made in between.
/// It can only be expected to work if the executable is looking
/// at a portion of the data where the C++ class definitions didn't
/// change.
///
/// This has not been thoroughly tested so your mileage may vary.
///


#include <iostream>
using namespace std;

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include "HDDM/JEventSource_HDDM.h"
#include "BCAL/JFactory_DBCALHit.h"
#include "CDC/JFactory_DCDCHit.h"
#include "CHERENKOV/JFactory_DCHERENKOVHit.h"
#include "FCAL/JFactory_DFCALHit.h"
#include "FDC/JFactory_DFDCHit.h"
#include "TRACKING/JFactory_DMCTrackHit.h"
#include "TRACKING/JFactory_DMCThrown.h"
#include "TOF/JFactory_DTOFHit.h"
#include "UPV/JFactory_DUPVHit.h"

extern "C" {
extern InitPlugin_t InitPlugin;
}

void InitPlugin(JApplication *app)
{
	cout<<"Adding factories from shared object ..."<<endl;	
	eventLoop->AddFactory(new JFactory_DBCALHit());
	eventLoop->AddFactory(new JFactory_DCDCHit());
	eventLoop->AddFactory(new JFactory_DCHERENKOVHit());
	eventLoop->AddFactory(new JFactory_DFCALHit());
	eventLoop->AddFactory(new JFactory_DFDCHit());
	eventLoop->AddFactory(new JFactory_DMCTrackHit());
	eventLoop->AddFactory(new JFactory_DMCThrown());
	eventLoop->AddFactory(new JFactory_DTOFHit());
	eventLoop->AddFactory(new JFactory_DUPVHit());
}
