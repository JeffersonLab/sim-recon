// $Id$
//
//    File: JEventProcessor_recon2mc.cc
// Created: Tue Nov 10 13:07:57 EST 2015
// Creator: davidl (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include <sstream>

#include "JEventProcessor_recon2mc.h"
using namespace jana;

#include <TRACKING/DTrackTimeBased.h>
#include <particleType.h>

// These parameters are exposed as JANA config.
// parameters below in init. 
string OUTFILENAME = "recon2mc.hddm";
double MIN_FOM = 1.0E-2; // minimum FOM to accept
double MIN_P =  0.0; // minimum momentum in GeV/c
double MAX_P = 20.0; // maximum momentum in GeV/c
vector<int> pids_to_keep;
double VX = -1000.0; // vertex to override reconstructed one
double VY = -1000.0; // This is only used if the VERTEX config.
double VZ = -1000.0; // parameter is set.
bool OVERRIDE_VERTEX = false;

using namespace hddm_s;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_recon2mc());
}
} // "C"


//------------------
// JEventProcessor_recon2mc (Constructor)
//------------------
JEventProcessor_recon2mc::JEventProcessor_recon2mc()
{
	ostr_s = NULL;
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~JEventProcessor_recon2mc (Destructor)
//------------------
JEventProcessor_recon2mc::~JEventProcessor_recon2mc()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_recon2mc::init(void)
{
	// Open output HDDM file
	ofs.open(OUTFILENAME.c_str());
	if (! ofs.is_open()) {
		std::cout << " Error opening output file \"" << OUTFILENAME << "\"!" << std::endl;
		exit(-1);
	}
	ostr_s = new hddm_s::ostream(ofs);
	//ostr_s->setCompression(hddm_s::k_bz2_compression);     // hdgeant can't handle compressed files
	//ostr_s->setIntegrityChecks(hddm_s::k_crc32_integrity); // hdgeant can't handle integrity checks

	string pidlist = "8,9";
	string vertex = "";
	gPARMS->SetDefaultParameter("OUTFILENAME", OUTFILENAME, "Filename for output HDDM file");
	gPARMS->SetDefaultParameter("MIN_FOM", MIN_FOM, "Minimum tracking FOM for track to be passed to output");
	gPARMS->SetDefaultParameter("MIN_P", MIN_P, "Minimum reconstructed track momentum in GeV/c for track to be passed to output");
	gPARMS->SetDefaultParameter("MAX_P", MAX_P, "Maximum reconstructed track momentum in GeV/c for track to be passed to output");
	gPARMS->SetDefaultParameter("PIDLIST", pidlist, "Comma separated list of GEANT particle numbers indicating types of particles to keep. Empty string means keep them all (probably not what you want)");
	gPARMS->SetDefaultParameter("VERTEX", vertex, "Comma separated vertex coordinates in cm. If empty (default) the reconstructed vertex is used.");
	
	// Parse PIDLIST
	if(pidlist.length()>0){
		stringstream ss(pidlist);
		int i;
		while (ss >> i){
			pids_to_keep.push_back(i);
			if (ss.peek() == ',') ss.ignore();
		}
	}

	// Parse VERTEX
	if(vertex.length()>0){
		stringstream ss(vertex);
		
		ss >> VX;
		if (ss.peek() == ',') ss.ignore();
		ss >> VY;
		if (ss.peek() == ',') ss.ignore();
		ss >> VZ;
		
		OVERRIDE_VERTEX = true;		
	}
	
	jout << "=========================================" << endl;
	jout << "recon2mc settings:" << endl;
	jout << "-------------------" << endl;
	jout << " OUTFILENAME: " << OUTFILENAME << endl;
	jout << "     MIN_FOM: " << MIN_FOM << endl;
	jout << "       MIN_P: " << MIN_P << " GeV/c" << endl;
	jout << "       MAX_P: " << MAX_P << " GeV/c" << endl;
	jout << "PIDs to keep: ";
	for(uint32_t i=0; i<pids_to_keep.size(); i++) jout << pids_to_keep[i] << ", ";
	jout << endl;
	jout << "      vertex: ";
	if(OVERRIDE_VERTEX){
		jout << VX << ", " << VY << ", " << VZ << " cm" << endl;
	}else{
		jout << "<from reconstructed track>" << endl;
	}
	
	jout << "=========================================" << endl;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_recon2mc::brun(JEventLoop *eventLoop, int runnumber)
{
	runNumber = runnumber;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_recon2mc::evnt(JEventLoop *loop, int eventnumber)
{
	// Get list of tracks
	vector<const DTrackTimeBased*> tbts;
	loop->Get(tbts);
	
	// Copy tracks we want to keep to special list
	vector<const DTrackTimeBased*> tbts_to_keep;
	for(uint32_t i=0; i<tbts.size(); i++){
		const DTrackTimeBased *tbt = tbts[i];

		// Filter poorly fit tracks out
		if(tbt->FOM < MIN_FOM) continue;
		
		// Filter tracks outside acceptable momentum range
		double p = tbt->pmag();
		if( p<MIN_P || p>MAX_P) continue;
		
		// Optionally filter out undesired types
		if( !pids_to_keep.empty() ){
			int pid = tbt->PID();
			bool keep = false;
			for(uint32_t i=0; i<pids_to_keep.size(); i++){
				if( pid == pids_to_keep[i] ){
					keep = true;
					break;
				}
			}
			if(!keep) continue;
		}

		tbts_to_keep.push_back(tbt);
	}
	
	// Don't write out events where we found no tracks of interest
	if(tbts_to_keep.empty()) return NOERROR;
	
	//==================================================
	//    PLACE FILTER CODE HERE IF YOU WANT TO
	//    CUT ON OTHER CONDITIONS. TO IGNORE AN
	//    EVENT, JUST RETURN "NOERROR" IF YOU WANT
	//    TO THROW AWAY THE EVENT. IF YOU WANT TO
	//    FILTER OUT SPECIFIC TRACKS, THEN EITHER
	//    DO IT IN THE LOOP ABOVE OR REMOVE THE
	//    ENTRY FROM tbts_to_keep.
	//==================================================

	// Create new HDDM event
	hddm_s::HDDM hddm;
	hddm.addPhysicsEvents();
	PhysicsEvent &pe = hddm.getPhysicsEvent();
	pe.setEventNo(eventnumber);
	pe.setRunNo(runNumber);
	
	pe.addReactions();
	Reaction &re = pe.getReaction();
	re.setType(0);
	re.addVertices(tbts_to_keep.size()); // add separate vertex for each track

	for(uint32_t i=0; i<tbts_to_keep.size(); i++){
		const DTrackTimeBased *tbt = tbts_to_keep[i];
		Vertex &vtx = re.getVertex(i);

		vtx.addOrigins();		
		Origin &origin = vtx.getOrigin();		
		origin.setT(  tbt->time() );
		origin.setVx( OVERRIDE_VERTEX ? VX:tbt->x()    );
		origin.setVy( OVERRIDE_VERTEX ? VY:tbt->y()    );
		origin.setVz( OVERRIDE_VERTEX ? VZ:tbt->z()    );
		
		vtx.addProducts();
		Product &prod = vtx.getProduct();
		prod.setMech(0);
		prod.setParentid(0);
		prod.setPdgtype( PDGtype(tbt->PID()) );
		prod.setType(tbt->PID());

		prod.addMomenta();
		Momentum &mom = prod.getMomentum();
		mom.setE(  tbt->energy() );
		mom.setPx( tbt->px()     );
		mom.setPy( tbt->py()     );
		mom.setPz( tbt->pz()     );

		prod.addPropertiesList();
		Properties &prop = prod.getProperties();
		prop.setCharge( tbt->charge() );
		prop.setMass(   tbt->mass()   );
		
	}
	
	// Write hddm event to output
	pthread_mutex_lock(&mutex);
	*ostr_s << hddm;
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_recon2mc::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_recon2mc::fini(void)
{

	// Close output HDDM file
	delete ostr_s;
	ofs.close();

	return NOERROR;
}

