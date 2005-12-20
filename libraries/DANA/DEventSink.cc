// $Id$
//
//    File: DEventSink.cc
// Created: Mon Dec 19 16:15:49 EST 2005
// Creator: davidl (on Linux phecda 2.6.9-11.ELsmp athlon)
//

#pragma implementation "DEventSink.h"

#include "DEventSink.h"
#include "DEventLoop.h"

//---------------------------------
// DEventSink    (Constructor)
//---------------------------------
DEventSink::DEventSink()
{
	initialized = false;
	pthread_mutex_init(&sink_mutex, NULL);
}

//---------------------------------
// brun
//---------------------------------
derror_t DEventSink::brun(DEventLoop *loop, int runnumber)
{
	derror_t err=NOERROR;

	// We want to make sure that the factory write list is
	// generated only once.
	if(initialized)return NOERROR;
	LockSink();
	if(!initialized){
	
		err = brun_sink(loop, runnumber);
		initialized = true;
	}
	UnlockSink();

	return err;
}

//---------------------------------
// AddToWriteList
//---------------------------------
void DEventSink::AddToWriteList(string name, string tag)
{
	// We don't want to add a factory to the list twice!
	if(IsInWriteList(name,tag))return;

	factory_name_spec_t f;
	f.name = name;
	f.tag = tag;
	factories_to_write.push_back(f);
}

//---------------------------------
// AddAllToWriteList
//---------------------------------
void DEventSink::AddAllToWriteList(DEventLoop *loop)
{
	// Get list of all factories
	vector<DFactory_base*> factories = loop->GetFactories();
	for(unsigned int i=0; i<factories.size(); i++){
		if(factories[i]->GetNrows()==0)continue;
		if(!factories[i]->TestFactoryFlag(DFactory_base::WRITE_TO_OUTPUT))continue;
		AddToWriteList(factories[i]->dataClassName(), factories[i]->Tag());
	}
}

//---------------------------------
// RemoveFromWriteList
//---------------------------------
void DEventSink::RemoveFromWriteList(string name, string tag)
{
	vector<factory_name_spec_t>::iterator iter = factories_to_write.begin();
	for(; iter != factories_to_write.end(); iter++){
		if((*iter).name == name){
			if((*iter).tag == tag){
				factories_to_write.erase(iter);
				return;
			}
		}
	}
}

//---------------------------------
// ClearWriteList
//---------------------------------
void DEventSink::ClearWriteList(void)
{
	factories_to_write.clear();
}

