// $Id$
//
//    File: JEventProcessor_dumpcandidates.cc
// Created: Tue Feb  4 09:29:35 EST 2014
// Creator: davidl (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_dumpcandidates.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>

#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackHitSelector.h>
#include <TRACKING/DTrackFitter.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_dumpcandidates());
}
} // "C"


//------------------
// JEventProcessor_dumpcandidates (Constructor)
//------------------
JEventProcessor_dumpcandidates::JEventProcessor_dumpcandidates()
{

}

//------------------
// ~JEventProcessor_dumpcandidates (Destructor)
//------------------
JEventProcessor_dumpcandidates::~JEventProcessor_dumpcandidates()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_dumpcandidates::init(void)
{
	
	// Open output file
	ofs = new ofstream("gluex_candidates.txt");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_dumpcandidates::brun(JEventLoop *eventLoop, int runnumber)
{
	// Get pointer to geometry object
	DApplication* dapp=dynamic_cast<DApplication*>(japp);
	dgeom = dapp->GetDGeometry(1);
	if(!dgeom){
		jerr << "Couldn't get DGeometry pointer!!" << endl;
		return UNKNOWN_ERROR;
	}

	// Get CDC wires
	vector<vector<DCDCWire *> > cdcwires;
	dgeom->GetCDCWires(cdcwires);

	// Get FDC wires
	vector<vector<DFDCWire *> > fdcwires;
	dgeom->GetFDCWires(fdcwires);
	
	// Generate map keyed by wire object's address
	// (in the form of an unsigned long) and whose
	// value is the index of the wire (i.e. position
	// as produced by dumpwires)
	int index = 0;
	for(unsigned int i=0; i<cdcwires.size(); i++){
		vector<DCDCWire *> &wires = cdcwires[i];
		for(unsigned int j=0; j<wires.size(); j++){

			DCDCWire *w = wires[j];
			unsigned long addr = (unsigned long)w;
			wireID[addr] = index++;
		}
	}

	for(unsigned int i=0; i<fdcwires.size(); i++){
		vector<DFDCWire *> &wires = fdcwires[i];
		for(unsigned int j=0; j<wires.size(); j++){

			DFDCWire *w = wires[j];
			unsigned long addr = (unsigned long)w;
			wireID[addr] = index++;
		}
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_dumpcandidates::evnt(JEventLoop *loop, int eventnumber)
{
	// Get track candidates
	vector<const DTrackCandidate*> candidates;
	loop->Get(candidates);

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"Unable to get a DTrackHitSelector object!"<<endl;
		return UNKNOWN_ERROR;
	}
	const DTrackHitSelector * hitselector = hitselectors[0];

	// Get pointer to DTrackFitter object
	vector<const DTrackFitter *> fitters;
	loop->Get(fitters);
	if(fitters.size()<1){
		_DBG_<<"Unable to get a DTrackFitter object!"<<endl;
		return UNKNOWN_ERROR;
	}
	DTrackFitter *fitter = const_cast<DTrackFitter*>(fitters[0]);
	
	for(unsigned int i=0; i<candidates.size(); i++){
		const DTrackCandidate *can = candidates[i];
		
		// Create DReferenceTrajectory for this candidate
		DReferenceTrajectory *rt = new DReferenceTrajectory(fitter->GetDMagneticFieldMap());
		rt->SetDGeometry(dgeom);
      rt->q = can->charge();
		rt->SetMass(0.1396);
		rt->Swim(can->position(),can->momentum(),can->charge());
		
		// Get hits to be used for the fit
		vector<const DCDCTrackHit*> cdctrackhits;
		vector<const DFDCPseudo*> fdcpseudos;
		loop->Get(cdctrackhits);
		loop->Get(fdcpseudos);
		hitselector->GetAllHits(DTrackHitSelector::kHelical, rt, cdctrackhits, fdcpseudos, fitter);

		// Get list of wire ids
		vector<int> wire_ids;
		for(unsigned int j=0; j<cdctrackhits.size(); j++){
			wire_ids.push_back(GetWireIndex(cdctrackhits[j]->wire));
		}
		for(unsigned int j=0; j<fdcpseudos.size(); j++){
			wire_ids.push_back(GetWireIndex(fdcpseudos[j]->wire));
		}

		// Write candidate to string
		stringstream ss;
		ss << can->charge();
		ss << " " << can->x() << " " << can->y() << " " << can->z();
		ss << " " << can->px() << " " << can->py() << " " << can->pz();
		for(unsigned int j=0; j<wire_ids.size(); j++){
			ss << " " << wire_ids[j];
		}
		
		// Write candidate string to file
		(*ofs) << ss.str() << endl;
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_dumpcandidates::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_dumpcandidates::fini(void)
{
	if(ofs){
		ofs->close();
		delete ofs;
		ofs =NULL;
	}

	return NOERROR;
}



