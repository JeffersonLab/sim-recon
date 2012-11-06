// $Id$
//
// File: DEventProcessor_pidstudies_tree.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_pidstudies_tree.h"
#include "particleType.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_pidstudies_tree());
}
} // "C"

bool Compare_TrackMatches(DEventProcessor_pidstudies_tree::plugin_trackmatch_t *locTrackMatch1, DEventProcessor_pidstudies_tree::plugin_trackmatch_t *locTrackMatch2){
	return (locTrackMatch1->dFOM > locTrackMatch2->dFOM);
};

//------------------
// init
//------------------
jerror_t DEventProcessor_pidstudies_tree::init(void)
{

	dMCReconstructionStatuses = new MCReconstructionStatuses();
	dPluginTree_MCReconstructionStatuses = new TTree("dPluginTree_MCReconstructionStatuses", "MC Reconstruction Statuses");
	dPluginTree_MCReconstructionStatuses->Branch("dPluginBranch_MCReconstructionStatuses", "MCReconstructionStatuses", &dMCReconstructionStatuses);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_pidstudies_tree::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_pidstudies_tree::evnt(JEventLoop *loop, int eventnumber)
{
	unsigned int loc_i, loc_j;
	DLorentzVector locThrownFourMomentum, locThrownSpacetimeVertex;
	DLorentzVector locReconstructedFourMomentum, locReconstructedSpacetimeVertex;
	const DChargedTrackHypothesis *locChargedTrackHypothesis;
	MCReconstructionStatus *locMCReconstructionStatus;
	ReconstructedHypothesis *locReconstructedHypothesis;
	map<const DMCThrown*, const DChargedTrack*> locTrackMap;
	map<const DMCThrown*, const DChargedTrack*>::iterator locTrackMapIterator;

	vector<const DMCThrown*> locMCThrownVector;
	const DMCThrown *locMCThrown;
	loop->Get(locMCThrownVector);

	vector<const DChargedTrack*> locChargedTrackVector;
	const DChargedTrack *locChargedTrack;
	loop->Get(locChargedTrackVector);

	// Find track matches
	//Assume reconstructed track p/theta not much different between hypotheses: just choose the first hypothesis
	double locFOM, locMinimumFOM = 0.1;
	plugin_trackmatch_t *locTrackMatch;
	dTrackMatches.resize(0);
	for(loc_i = 0; loc_i < locMCThrownVector.size(); loc_i++){
		locMCThrown = locMCThrownVector[loc_i];
		for(loc_j = 0; loc_j < locChargedTrackVector.size(); loc_j++){
			locChargedTrack = locChargedTrackVector[loc_j];
			locFOM = Calc_MatchFOM(locMCThrown->momentum(), locChargedTrack->dChargedTrackHypotheses[0]->dTrackTimeBased->momentum());
			if(locFOM > locMinimumFOM){
				locTrackMatch = new plugin_trackmatch_t();
				locTrackMatch->dMCThrown = locMCThrown;
				locTrackMatch->dChargedTrack = locChargedTrack;
				locTrackMatch->dFOM = locFOM;
				dTrackMatches.push_back(locTrackMatch);
			}
		}
	}
	//sort track matches by fom
	sort(dTrackMatches.begin(), dTrackMatches.end(), Compare_TrackMatches);

	//fill track map in order of best matches
	while(dTrackMatches.size() > 0){
		locMCThrown = dTrackMatches[0]->dMCThrown;
		locChargedTrack = dTrackMatches[0]->dChargedTrack;
		locTrackMap[locMCThrown] = locChargedTrack;
		//erase all matches that include these objects
		for(loc_i = dTrackMatches.size() - 1; loc_i >= 0; loc_i--){
			locTrackMatch = dTrackMatches[loc_i];
			if(locTrackMatch->dMCThrown == locMCThrown)
				dTrackMatches.erase(dTrackMatches.begin() + loc_i);
			else if(locTrackMatch->dChargedTrack == locChargedTrack)
				dTrackMatches.erase(dTrackMatches.begin() + loc_i);
			if(loc_i == 0) //unsigned int!
				break;
		}
	}

	//fill information into trees
	dMCReconstructionStatuses->dMCReconstructionStatusVector.resize(0);
	for(loc_i = 0; loc_i < locMCThrownVector.size(); loc_i++){
		locMCThrown = locMCThrownVector[loc_i];
		locMCReconstructionStatus = new MCReconstructionStatus();
		//set thrown track information
		locMCReconstructionStatus->dThrownID = Particle_t(locMCThrown->type);
		locThrownFourMomentum.SetVect(locMCThrown->momentum());
		locThrownFourMomentum.SetT(locMCThrown->energy());
		locMCReconstructionStatus->dThrownFourMomentum = locThrownFourMomentum;
		locThrownSpacetimeVertex.SetVect(locMCThrown->position());
		locThrownSpacetimeVertex.SetT(locMCThrown->t0());
		locMCReconstructionStatus->dThrownSpacetimeVertex = locThrownSpacetimeVertex;
		//see if any match to charged particles
		locChargedTrack = NULL;
		for(locTrackMapIterator = locTrackMap.begin(); locTrackMapIterator != locTrackMap.end(); locTrackMapIterator++){
			if(locMCThrown == (*locTrackMapIterator).first){
				locChargedTrack = (*locTrackMapIterator).second;
				break;
			}
		}
		if(locChargedTrack == NULL){
			//set reconstructed hypothesis information
			dMCReconstructionStatuses->dMCReconstructionStatusVector.push_back(locMCReconstructionStatus);
			continue;
		}

		//set reconstructed hypothesis information
		for(loc_i = 0; loc_i < locChargedTrack->dChargedTrackHypotheses.size(); loc_i++){
			locChargedTrackHypothesis = locChargedTrack->dChargedTrackHypotheses[loc_i];
			locReconstructedHypothesis = new ReconstructedHypothesis();

			locReconstructedHypothesis->dPID = locChargedTrackHypothesis->dPID;
			locReconstructedFourMomentum.SetVect(locChargedTrackHypothesis->dTrackTimeBased->momentum());
			locReconstructedFourMomentum.SetT(locChargedTrackHypothesis->dTrackTimeBased->energy());
			locReconstructedSpacetimeVertex.SetVect(locChargedTrackHypothesis->dTrackTimeBased->position());
			locReconstructedSpacetimeVertex.SetT(locChargedTrackHypothesis->dTrackTimeBased->t0());
			locReconstructedHypothesis->dFourMomentum = locReconstructedFourMomentum;
			locReconstructedHypothesis->dSpacetimeVertex = locReconstructedSpacetimeVertex;

			locReconstructedHypothesis->dChiSq_Overall = locChargedTrackHypothesis->dChiSq;
			locReconstructedHypothesis->dNDF_Overall = locChargedTrackHypothesis->dNDF;

			locReconstructedHypothesis->dChiSq_Tracking = locChargedTrackHypothesis->dTrackTimeBased->chisq;
			locReconstructedHypothesis->dNDF_Tracking = locChargedTrackHypothesis->dTrackTimeBased->Ndof;

			locReconstructedHypothesis->dChiSq_DCdEdx = locChargedTrackHypothesis->dTrackTimeBased->chi2_dedx;
			locReconstructedHypothesis->dNDF_DCdEdx = 1;

			locReconstructedHypothesis->dChiSq_Timing = locChargedTrackHypothesis->dChiSq_Timing;
			locReconstructedHypothesis->dNDF_Timing = locChargedTrackHypothesis->dNDF_Timing;

			locReconstructedHypothesis->dChiSq_TOFdEdx = 0.0;
			locReconstructedHypothesis->dNDF_TOFdEdx = 0;

			locReconstructedHypothesis->dChiSq_BCALdEdx = 0.0;
			locReconstructedHypothesis->dNDF_BCALdEdx = 0;

			locMCReconstructionStatus->dReconstructedHypothesisVector.push_back(locReconstructedHypothesis);
		}

		//set reconstructed hypothesis information
		dMCReconstructionStatuses->dMCReconstructionStatusVector.push_back(locMCReconstructionStatus);
	}
	dPluginTree_MCReconstructionStatuses->Fill();

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_pidstudies_tree::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_pidstudies_tree::fini(void)
{
	return NOERROR;
}

//Copied from DEventProcessor_phys_tree.cc
double DEventProcessor_pidstudies_tree::Calc_MatchFOM(const DVector3& locMomentum_Track1, const DVector3& locMomentum_Track2) const
{
	// This is a kind of brain-dead algorithm. It wants to use both the
	// momentum direction and magnitude to determine the FOM. For the
	// magnitude, we use the curvature which becomes close to zero for
	// high momentum tracks (good because a 4GeV/c and an 8GeV/c track
	// both look more or less like straight lines).
	//
	// For the direction, we just use the relative angle between the
	// two tracks.
	//
	// For both, we take a reciprocal so that the closer the match,
	// the higher the FOM. We take a product of the 2 FOM components
	// so that both components must have a high value in order for the
	// total FOM to be large.

	double epsilon = 1.0E-6; // prevent reciprocals from resulting in infinity

	double curature_a = 1.0/locMomentum_Track1.Mag();
	double curature_b = 1.0/locMomentum_Track1.Mag();
	double curvature_diff = fabs(curature_a - curature_b);
	double curvature_fom = 1.0/(curvature_diff + epsilon);

	double theta_rel = fabs(locMomentum_Track1.Angle(locMomentum_Track2));
	double theta_fom = 1.0/(theta_rel + epsilon);
	
	double fom = curvature_fom*theta_fom;

	return fom;
}

