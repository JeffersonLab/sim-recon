// $Id$
//
//    File: DL3Trigger_factory.cc
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DL3Trigger_factory.h"
using namespace jana;

#include <START_COUNTER/DSCDigiHit.h>
#include <TOF/DTOFDigiHit.h>
#include <BCAL/DBCALPoint.h>
#include <BCAL/DBCALCluster.h>
#include <FCAL/DFCALCluster.h>
#include <TRACKING/DTrackCandidate.h>

#include <TRACKING/DTrackWireBased.h>
#include <BCAL/DBCALCluster.h>
#include <TRIGGER/DL1Trigger.h>

//------------------
// init
//------------------
jerror_t DL3Trigger_factory::init(void)
{
	FRACTION_TO_KEEP = 1.0;
	DO_WIRE_BASED_TRACKING = false;
	DO_BCAL_CLUSTER = false;
	L1_TRIG_MASK    = 0xffffffff;
	L1_FP_TRIG_MASK = 0xffffffff;
	MVA_WEIGHTS = "";
	MVA_CUT = -0.2;
	mvareader = NULL;

	gPARMS->SetDefaultParameter("L3:FRACTION_TO_KEEP", FRACTION_TO_KEEP ,"Random Fraction of event L3 should keep. (Only used for debugging).");
	gPARMS->SetDefaultParameter("L3:DO_WIRE_BASED_TRACKING", DO_WIRE_BASED_TRACKING ,"Activate wire-based tracking for every event");
	gPARMS->SetDefaultParameter("L3:DO_BCAL_CLUSTER", DO_BCAL_CLUSTER ,"Activate BCAL clusters for every event");
	gPARMS->SetDefaultParameter("L3:L1_TRIG_MASK", L1_TRIG_MASK ,"Discard events that don't have one of these bits set in DL1Trigger::trig_mask (or in L1_FP_TRIG_MASK)");
	gPARMS->SetDefaultParameter("L3:L1_FP_TRIG_MASK", L1_FP_TRIG_MASK ,"Discard events that don't have one of these bits set in DL1Trigger::fp_trig_mask (or in L1_TRIG_MASK)");
	gPARMS->SetDefaultParameter("L3:MVA_WEIGHTS", MVA_WEIGHTS ,"TMVA weights file");
	gPARMS->SetDefaultParameter("L3:MVA_CUT", MVA_CUT ,"Cut on MVA response function. Event with values less than this are discarded.");

	if(MVA_WEIGHTS != ""){
#ifdef HAVE_TMVA
		mvareader = new TMVA::Reader();
		mvareader->AddVariable("Nstart_counter",      &Nstart_counter);
		mvareader->AddVariable("Ntof",                &Ntof);
		mvareader->AddVariable("Nbcal_points",        &Nbcal_points);
		mvareader->AddVariable("Nbcal_clusters",      &Nbcal_clusters);
		mvareader->AddVariable("Ebcal_points",        &Ebcal_points);
		mvareader->AddVariable("Ebcal_clusters",      &Ebcal_clusters);
		mvareader->AddVariable("Nfcal_clusters",      &Nfcal_clusters);
		mvareader->AddVariable("Efcal_clusters",      &Efcal_clusters);
		mvareader->AddVariable("Ntrack_candidates",   &Ntrack_candidates);
		mvareader->AddVariable("Ptot_candidates",     &Ptot_candidates);
		
		mvareader->BookMVA("MVA", MVA_WEIGHTS);
#endif
	}

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DL3Trigger_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL3Trigger_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// Simple pass-through L3 trigger
	// algorithm = 0x1
	
	

	DL3Trigger *l3trig = new DL3Trigger(DL3Trigger::kKEEP_EVENT, 0x0L, 0x1);
	_data.push_back(l3trig);

	if(FRACTION_TO_KEEP!=1.0){
		double r = (double)random()/(double)RAND_MAX;
		if(r > FRACTION_TO_KEEP) l3trig->L3_decision = DL3Trigger::kDISCARD_EVENT;
	}
	
	// If L1 trigger filter is being applied, do that here
	if( (L1_TRIG_MASK!=0xffffffff) || (L1_FP_TRIG_MASK!=0xffffffff) ){

		vector<const DL1Trigger*> l1triggers;
		loop->Get(l1triggers);
		bool trig_bit_is_set = false;
		for(auto t : l1triggers){
			if( t->trig_mask&L1_TRIG_MASK       ) trig_bit_is_set = true;
			if( t->fp_trig_mask&L1_FP_TRIG_MASK ) trig_bit_is_set = true;
		}
		if(!trig_bit_is_set) l3trig->L3_decision = DL3Trigger::kDISCARD_EVENT;
	}
	
#ifdef HAVE_TMVA
	if(mvareader){
		vector<const DSCDigiHit*> scdigihits;
		vector<const DTOFDigiHit*> tofdigihits;
		vector<const DBCALPoint*> bcalpoints;
		vector<const DBCALCluster*> bcalclusters;
		vector<const DFCALCluster*> fcalclusters;
		vector<const DTrackCandidate*> trackcandidates;
		loop->Get(scdigihits);
		loop->Get(tofdigihits);
		loop->Get(bcalpoints);
		loop->Get(bcalclusters);
		loop->Get(fcalclusters);
		loop->Get(trackcandidates);

		// Calorimeter energies
		double Ebcal_points   = 0.0;
		double Ebcal_clusters = 0.0;
		double Efcal_clusters = 0.0;
		for(auto bp : bcalpoints  ) Ebcal_points   += bp->E();
		for(auto bc : bcalclusters) Ebcal_clusters += bc->E();
		for(auto fc : fcalclusters) Efcal_clusters += fc->getEnergy();

		// Ptot for candidates
		double Ptot_candidates = 0.0;
		for(auto tc : trackcandidates) Ptot_candidates += tc->momentum().Mag();

		Nstart_counter    = scdigihits.size();
		Ntof              = tofdigihits.size();
		Nbcal_points      = bcalpoints.size();
		Nbcal_clusters    = bcalclusters.size();
		Ebcal_points      = Ebcal_points;
		Ebcal_clusters    = Ebcal_clusters;
		Nfcal_clusters    = fcalclusters.size();
		Efcal_clusters    = Efcal_clusters;
		Ntrack_candidates = trackcandidates.size();
		Ptot_candidates   = Ptot_candidates;

		double v = mvareader->EvaluateMVA("MVA");
		if( v < MVA_CUT ) l3trig->L3_decision = DL3Trigger::kDISCARD_EVENT;
	}
#endif

	
	if(DO_WIRE_BASED_TRACKING){
		vector<const DTrackWireBased*> wbt;
		loop->Get(wbt);
	}

	if(DO_BCAL_CLUSTER){
		vector<const DBCALCluster*> bcalclusters;
		loop->Get(bcalclusters);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL3Trigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL3Trigger_factory::fini(void)
{
	return NOERROR;
}

