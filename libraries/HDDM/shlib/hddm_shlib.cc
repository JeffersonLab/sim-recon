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

#include <DApplication.h>
#include <DEventSourceHDDM.h>
#include <DEventLoop.h>

#include "DFactory_DBCALHit.h"
#include "DFactory_DCDCHit.h"
#include "DFactory_DCHERENKOVHit.h"
#include "DFactory_DFCALHit.h"
#include "DFactory_DFDCHit.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DTOFHit.h"
#include "DFactory_DUPVHit.h"

extern "C" {
extern GetDEventSourceType_t GetDEventSourceType;
extern MakeDEventSource_t MakeDEventSource;
extern InitFactories_t InitFactories;
}

const char* GetDEventSourceType(void){return "DEventSourceHDDM";}

DEventSource* MakeDEventSource(const char* name)
{
	return new DEventSourceHDDM(name);
}

void InitFactories(DEventLoop* eventLoop)
{
	cout<<"Adding factories from shared object ..."<<endl;	
	eventLoop->AddFactory(new DFactory_DBCALHit());
	eventLoop->AddFactory(new DFactory_DCDCHit());
	eventLoop->AddFactory(new DFactory_DCHERENKOVHit());
	eventLoop->AddFactory(new DFactory_DFCALHit());
	eventLoop->AddFactory(new DFactory_DFDCHit());
	eventLoop->AddFactory(new DFactory_DMCCheatHit());
	eventLoop->AddFactory(new DFactory_DMCThrown());
	eventLoop->AddFactory(new DFactory_DTOFHit());
	eventLoop->AddFactory(new DFactory_DUPVHit());
}
