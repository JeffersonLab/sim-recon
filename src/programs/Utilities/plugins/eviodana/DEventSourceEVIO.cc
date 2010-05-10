// $Id$
//
//    File: DEventSourceEVIO.cc
// Created: Sat May  8 13:54:35 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#include <iostream>
using namespace std;

#include <DANA/DApplication.h>

#include "DEventSourceEVIO.h"


//---------------------------------
// DEventSourceEVIO    (Constructor)
//---------------------------------
DEventSourceEVIO::DEventSourceEVIO(const char* source_name):JEventSource(source_name)
{
	// Open the EVIO file, catching any exceptions
	try {
		chan = new evioFileChannel(source_name,"r", 65536);
		chan->open();
		jout << "Opened EVIO file \""<<source_name<<"\" for reading"<<endl;
		
		// Fill in tagMap
		// We do this by hand for now, but will eventually get it from the XML dictionary file
		tagMap["DTrackTimeBased"]			= pair<int, int>(11500,0);
		tagMap["DTrackTimeBased.objId"]	= pair<int, int>(11500, 1);
		tagMap["DTrackTimeBased.chisq"]	= pair<int, int>(11500, 2);
		tagMap["DTrackTimeBased.Ndof"]	= pair<int, int>(11500, 3);
		tagMap["DTrackTimeBased.FOM"]		= pair<int, int>(11500, 4);
		tagMap["DTrackTimeBased.x"]		= pair<int, int>(11500, 5);
		tagMap["DTrackTimeBased.y"]		= pair<int, int>(11500, 6);
		tagMap["DTrackTimeBased.z"]		= pair<int, int>(11500, 7);
		tagMap["DTrackTimeBased.px"]		= pair<int, int>(11500, 8);
		tagMap["DTrackTimeBased.py"]		= pair<int, int>(11500, 9);
		tagMap["DTrackTimeBased.pz"]		= pair<int, int>(11500, 10);
		tagMap["DTrackTimeBased.q"]		= pair<int, int>(11500, 11);
		tagMap["DTrackTimeBased.E"]		= pair<int, int>(11500, 12);
		tagMap["DTrackTimeBased.mass"]	= pair<int, int>(11500, 13);
		tagMap["DTrackTimeBased.t0"]		= pair<int, int>(11500, 14);
		tagMap["DTrackTimeBased.assocObjectBanks"] = pair<int, int>(11500, 254);
		tagMap["DTrackTimeBased.assocObjects"] = pair<int, int>(11500, 255);

	}catch (evioException e) {
		jerr << e.toString() << endl;
		chan = NULL;
	} catch (...) {
		jerr << "?unknown exception" << endl;
		chan = NULL;
	}
	
	// Used for reference trajectories
	pthread_mutex_init(&rt_mutex, NULL);
}

//---------------------------------
// ~DEventSourceEVIO    (Destructor)
//---------------------------------
DEventSourceEVIO::~DEventSourceEVIO()
{
	// Close EVIO file, catching any exceptions
	try {
		if(chan){
			chan->close();
			delete(chan);
		}
	}catch (evioException e) {
		jerr << e.toString() << endl;
		chan = NULL;
	} catch (...) {
		jerr << "?unknown exception" << endl;
		chan = NULL;
	}
}

//---------------------------------
// GetEvent
//---------------------------------
jerror_t DEventSourceEVIO::GetEvent(JEvent &event)
{
	if(chan==NULL)return NO_MORE_EVENTS_IN_SOURCE;
	if(chan->read()){
		
		evioDOMTree *evt = new evioDOMTree(chan);

		// Copy the reference info into the JEvent object
		event.SetJEventSource(this);
		event.SetEventNumber(++Nevents_read);
		event.SetRunNumber(1);
		event.SetRef(evt);
		
	}else{
		return NO_MORE_EVENTS_IN_SOURCE;
	}

	return NOERROR;
}

//---------------------------------
// FreeEvent
//---------------------------------
void DEventSourceEVIO::FreeEvent(JEvent &event)
{
	evioDOMTree *evt = (evioDOMTree*)event.GetRef();
	if(evt)delete evt;

	// Check for DReferenceTrajectory objects we need to delete
	pthread_mutex_lock(&rt_mutex);
	map<evioDOMTree*, vector<DReferenceTrajectory*> >::iterator iter = rt_by_event.find(evt);
	if(iter != rt_by_event.end()){
		vector<DReferenceTrajectory*> &rts = iter->second;
		for(unsigned int i=0; i<rts.size(); i++)rt_pool.push_back(rts[i]);
		rt_by_event.erase(iter);
	}
	pthread_mutex_unlock(&rt_mutex);
}

//---------------------------------
// GetObjects
//---------------------------------
jerror_t DEventSourceEVIO::GetObjects(JEvent &event, JFactory_base *factory)
{
	/// This gets called through the virtual method of the
	/// JEventSource base class. It creates the objects of the type
	/// on which factory is based. It uses the evioDOMTree object
	/// kept in the ref field of the JEvent object passed.

	// We must have a factory to hold the data
	if(!factory)throw RESOURCE_UNAVAILABLE;

	evioDOMTree *evt = (evioDOMTree*)event.GetRef();
	if(!evt)throw RESOURCE_UNAVAILABLE;

	// Get pointer to the B-field object and Geometry object
	JEventLoop *loop = event.GetJEventLoop();
	if(loop){
		DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
		if(dapp){
			bfield = dapp->GetBfield();
			geom = dapp->GetDGeometry(event.GetRunNumber());
		}
	}

	// Get name of data class we're trying to extract
	string dataClassName = factory->GetDataClassName();

	if(dataClassName =="DTrackTimeBased")
	  return Extract_DTrackTimeBased(evt, dynamic_cast<JFactory<DTrackTimeBased>*>(factory));	

	return OBJECT_NOT_AVAILABLE;
}

//---------------------------------
// Extract_DTrackTimeBased
//---------------------------------
jerror_t DEventSourceEVIO::Extract_DTrackTimeBased(evioDOMTree *evt,  JFactory<DTrackTimeBased> *factory)
{
	// Note: Since this is a reconstructed factory, we want to generally return OBJECT_NOT_AVAILABLE
	// rather than NOERROR. The reason being that the caller interprets "NOERROR" to mean "yes I
	// usually can provide objects of that type, but this event has none." This will cause it to
	// skip any attempt at reconstruction. On the other hand, a value of "OBJECT_NOT_AVAILABLE" tells
	// it "I cannot provide those type of objects for this event.

  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

	vector<DTrackTimeBased*> data;
	vector<DReferenceTrajectory*> rts;
  
	bool event_had_tracktimebaseds = false;

	// Get list of DTrackTimeBased banks. At this level, there should usually be
	// only one but we get the list and loop over them since that is how it is formatted.
	// The actual DTrackTimeBased objects are contained in child banks of this (these).
	evioDOMNodeListP bList = evt->getNodeList(tagNumEquals(tagMap["DTrackTimeBased"]));

	// loop over all banks in list (should only be one)
	evioDOMNodeList::const_iterator iter1;
	for(iter1=bList->begin(); iter1!=bList->end(); iter1++) {

		evioDOMNodeList *members = (*iter1)->getChildList();

		try{
			const vector<uint64_t> &v_objId = GetVector<uint64_t>(members, "DTrackTimeBased.objId");
			const vector<float> &v_chisq = GetVector<float>(members, "DTrackTimeBased.chisq");
			const vector<int> &v_Ndof = GetVector<int>(members, "DTrackTimeBased.Ndof");
			const vector<float> &v_FOM = GetVector<float>(members, "DTrackTimeBased.FOM");
			const vector<float> &v_x = GetVector<float>(members, "DTrackTimeBased.x");
			const vector<float> &v_y = GetVector<float>(members, "DTrackTimeBased.y");
			const vector<float> &v_z = GetVector<float>(members, "DTrackTimeBased.z");
			const vector<float> &v_px = GetVector<float>(members, "DTrackTimeBased.px");
			const vector<float> &v_py = GetVector<float>(members, "DTrackTimeBased.py");
			const vector<float> &v_pz = GetVector<float>(members, "DTrackTimeBased.pz");
			const vector<float> &v_q = GetVector<float>(members, "DTrackTimeBased.q");
			//const vector<float> &v_E = GetVector<float>(members, "DTrackTimeBased.E");
			const vector<float> &v_mass = GetVector<float>(members, "DTrackTimeBased.mass");
			//const vector<float> &v_t0 = GetVector<float>(members, "DTrackTimeBased.t0");
			
			// Get enough DReferenceTrajectory objects for all of the DTrackTimeBased Objects
			// we're about to read in. This seems a little complicated, but that's because it
			// is expensive to allocate these things so we recycle as much as possible.
			list<DReferenceTrajectory*> my_rts;
			pthread_mutex_lock(&rt_mutex);
			for(unsigned int i=0; i<v_objId.size(); i++){
				while(my_rts.size() < v_objId.size()){
					if(rt_pool.size()>0){
						my_rts.push_back(rt_pool.back());
						rt_pool.pop_back();
					}else{
						my_rts.push_back(new DReferenceTrajectory(bfield));
					}
				}
			}
			pthread_mutex_unlock(&rt_mutex);
			
			// Loop over DTrackTimeBased objects
			event_had_tracktimebaseds = true;
			for(unsigned int i=0; i<v_objId.size(); i++){

				DVector3 pos(v_x[i], v_y[i], v_z[i]);
				DVector3 mom(v_px[i], v_py[i], v_pz[i]);

				DTrackTimeBased *track = new DTrackTimeBased();

				track->setMomentum(mom);
				track->setPosition(pos);
				track->setCharge(v_q[i]);
				track->setMass(v_mass[i]);
				track->chisq = v_chisq[i];
				track->Ndof = v_Ndof[i];
				track->FOM = v_FOM[i];
				track->id = v_objId[i];
				
				// We need to create a pool of DReferenceTrajectory objects that we can assign 
				// to the track as is done in DTrackTimeBased_factory.cc. For now though, we
				// skip that in order to focus on getting the rest of this working.
				DReferenceTrajectory *rt = my_rts.back();
				my_rts.pop_back();
				if(rt){
					rt->SetMass(track->mass());
					rt->SetDGeometry(geom);
					rt->Swim(pos, mom, track->charge());
					rts.push_back(rt);
				}
				track->rt = rt;

				data.push_back(track);
			}
			
		}catch(JException &e){
			cout<<e.toString()<<endl;
		}
	}

	// Copy into factory
	if(event_had_tracktimebaseds){
		factory->CopyTo(data);
		
		// Add DReferenceTrajectory objects to rt_by_event so they can be deleted later.
		// The rt_by_event maintains lists indexed by the hddm_s pointer since multiple
		// threads may be calling us. Note that we first look to see if a list already
		// exists for this event and append to it if it does. This is so we can the
		// same list for all objects that use DReferenceTrajectories.
		pthread_mutex_lock(&rt_mutex);
		map<evioDOMTree*, vector<DReferenceTrajectory*> >::iterator iter = rt_by_event.find(evt);
		if(iter != rt_by_event.end()){
			vector<DReferenceTrajectory*> &my_rts = iter->second;
			my_rts.insert(my_rts.end(), rts.begin(), rts.end());
		}else{
			rt_by_event[evt] = rts;
		}
		pthread_mutex_unlock(&rt_mutex);

		// If the event had a s_Tracktimebased_t pointer, then report back that
		// we read them in from the file. Otherwise, report OBJECT_NOT_AVAILABLE
		return NOERROR;
	}

	// If we get to here then there was not even a placeholder in the HDDM file.
	// Return OBJECT_NOT_AVAILABLE to indicate reconstruction should be tried.
	return OBJECT_NOT_AVAILABLE;
}
		
