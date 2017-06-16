// $Id$
//
//    File: DCustomAction_ee_convert_cuts.cc
// Created: Wed Jun 14 06:26:48 EDT 2017
// Creator: jrsteven (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#include "DCustomAction_ee_convert_cuts.h"

void DCustomAction_ee_convert_cuts::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		// Optional: Create a ROOT subfolder.
			//If another thread has already created the folder, it just changes to it. 
		// CreateAndChangeTo_Directory("MyDirName", "MyDirTitle");
			//make sub-directory content here
		// gDirectory->cd(".."); //return to the action directory

		//	(Optional) Example: Create a histogram. 
			// This function will return the histogram if already created by another thread. If not pre-existing, it will create and return it. 
			// Function arguments are identical to those used for the histogram constructors
		dHist_Vertex = GetOrCreate_Histogram<TH2I>("eeVertex", "; Vertex Z (cm); Vertex R (cm)", 2000, 0.0, 200.0, 500, 0., 25.);
		dHist_Vertex_DOCA = GetOrCreate_Histogram<TH2I>("eeVertex_DOCA", "; Vertex Z (cm); Vertex R (cm)", 2000, 0.0, 200.0, 500, 0., 25.);
		dHist_Vertex_LowMass = GetOrCreate_Histogram<TH2I>("eeVertex_LowMass", "; Vertex Z (cm); Vertex R (cm)", 2000, 0.0, 200.0, 500, 0., 25.);
		dHist_Vertex_Final = GetOrCreate_Histogram<TH2I>("eeVertex_Final", "; Vertex Z (cm); Vertex R (cm)", 2000, 0.0, 200.0, 500, 0., 25.);

		dHist_DOCAVertex_DOCA = GetOrCreate_Histogram<TH2I>("eeDOCAVertex_DOCA", "; Vertex Z (cm); Vertex R (cm)", 2000, 0.0, 200.0, 500, 0., 25.);
		dHist_DOCA_DOCAVertex = GetOrCreate_Histogram<TH2I>("DOCA_DOCAVertex", "; DOCA Vertex (cm); DOCA (cm)", 200, 0.0, 10.0, 200, 0., 10.);
		dHist_DOCA_DeltaPhi = GetOrCreate_Histogram<TH2I>("DOCA_DeltaPhi", "; #Delta#phi (degrees); DOCA (cm)", 360, -180.0, 180.0, 200, 0., 10.);
		dHist_DOCA_DeltaTheta = GetOrCreate_Histogram<TH2I>("DOCA_DeltaTheta", "; #Delta#theta (degrees); DOCA (cm)", 100, -30.0, 30.0, 200, 0., 10.);
		dHist_DeltaPhi_DeltaTheta = GetOrCreate_Histogram<TH2I>("DeltaPhi_DeltaTheta", "; #Delta#theta (degrees); #Delta#phi (degrees)", 100, -30.0, 30.0, 360, -180., 180.);
		dHist_DOCA_PhiV = GetOrCreate_Histogram<TH2I>("DOCA_PhiV", "; #phi_{V} (rad); DOCA (cm)", 100, 0, 3.14, 200, 0., 10.);

		dHist_VertexZ_DOCAVertexZ = GetOrCreate_Histogram<TH2I>("VertexZ_DOCAVertexZ", "; DOCA Vertex Z (cm); Vertex Z (cm)", 200, 0.0, 200.0, 200, 0., 200.);
		dHist_VertexR_DOCAVertexR = GetOrCreate_Histogram<TH2I>("VertexR_DOCAVertexR", "; DOCA Vertex R (cm); Vertex R (cm)", 100, 0.0, 25.0, 100, 0., 25.);
		dHist_Mee = GetOrCreate_Histogram<TH1I>("Mee", "; M_{ee} (GeV)", 200, 0., 1.0);

		dHist_EOverP = GetOrCreate_Histogram<TH2I>("EOverP", "; Electron E/p; Positron E/p", 100, 0., 2., 100, 0., 2.);
		dHist_EOverP_BothFCAL = GetOrCreate_Histogram<TH2I>("EOverP_BothFCAL", "; Electron E/p; Positron E/p", 100, 0., 2., 100, 0., 2.);
		dHist_EOverP_BothBCAL = GetOrCreate_Histogram<TH2I>("EOverP_BothBCAL", "; Electron E/p; Positron E/p", 100, 0., 2., 100, 0., 2.);
		dHist_EOverP_FCALBCAL = GetOrCreate_Histogram<TH2I>("EOverP_FCALBCAL", "; Electron E/p; Positron E/p", 100, 0., 2., 100, 0., 2.);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_ee_convert_cuts::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	// Optional: Useful utility functions.
	const DAnalysisUtilities* locAnalysisUtilities;
	locEventLoop->GetSingle(locAnalysisUtilities);

	//Optional: check whether the user wanted to use the kinematic fit results when performing this action
	bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();

	// Get q+q- pair from combo
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);
	deque<const DKinematicData*> locParticles;
	if(locUseKinFitResultsFlag)
		locParticleComboStep->Get_FinalParticles(locParticles);
	else 
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	deque<const DKinematicData*> locParticlesMeasured;
	locParticleComboStep->Get_FinalParticles_Measured(locParticlesMeasured);

	// Charged track properties (require good FOM)
	deque<const DTrackTimeBased*> locTracksTimeBased;
	double locEOverP[2] = {0., 0.};
	bool locFCAL[2] = {false, false};
	bool locBCAL[2] = {false, false};
	for(size_t loc_i = 0; loc_i < 2; ++loc_i) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
		if(locChargedTrack == NULL) continue; // should never happen
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		locTracksTimeBased.push_back(locTrackTimeBased);

		// track FOM cut
		//if(locTrackTimeBased->FOM < 0.0001)
		//	return false;

		// BCAL/FCAL shower matches
		const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
		const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();

		if(locBCALShowerMatchParams != NULL) {
			const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
			locEOverP[loc_i] = locBCALShower->E/locParticles[loc_i]->lorentzMomentum().P();
			locBCAL[loc_i] = true;
		}
		if(locFCALShowerMatchParams != NULL) {
			const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
			locEOverP[loc_i] = locFCALShower->getEnergy()/locParticles[loc_i]->lorentzMomentum().P();
			locFCAL[loc_i] = true;
		}	
	}
	if(locTracksTimeBased.size() < 2) 
		return false;

	// Vertex position from KinFit
	TVector3 locVertex = locParticles[0]->position();

	// DOCA between tracks and vertex from POCA
	double locDOCA = locAnalysisUtilities->Calc_DOCA(locParticlesMeasured[0], locParticlesMeasured[1]);
	TVector3 locDOCAVertex;
	double locDOCAVertexVal = locAnalysisUtilities->Calc_DOCAVertex(locParticlesMeasured[0], locParticlesMeasured[1], locDOCAVertex);

	// Angles between tracks
	double locDeltaPhi = (locParticlesMeasured[0]->lorentzMomentum().Phi() - locParticlesMeasured[1]->lorentzMomentum().Phi())*180./TMath::Pi();
	if(locDeltaPhi > 180.) locDeltaPhi -= 360.;
	if(locDeltaPhi < -180.) locDeltaPhi += 360.;
	double locDeltaTheta = (locParticlesMeasured[0]->lorentzMomentum().Theta() - locParticlesMeasured[1]->lorentzMomentum().Theta())*180./TMath::Pi();

	// P4 for mass calculation?
	set<pair<const JObject*, unsigned int> > locSourceObjects;
	DLorentzVector locPairP4 = locAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 0, locSourceObjects, locUseKinFitResultsFlag);
	locPairP4 = locParticles[0]->lorentzMomentum();
	locPairP4 += locParticles[1]->lorentzMomentum();

	// Simon's track intersection via reference trajectories (only from EVIO files)
	//DVector3 pos;
	//double doca;
	double var_doca;
	DKinematicData kd1=*locTracksTimeBased[0],kd2=*locTracksTimeBased[1];
	if (locTracksTimeBased[0]->rt!=NULL) {
		//cout<<"found rt"<<endl;
		//locTracksTimeBased[0]->rt->IntersectTracks(locTracksTimeBased[1]->rt,&kd1,&kd2,pos,doca,var_doca);
		locTracksTimeBased[0]->rt->IntersectTracks(locTracksTimeBased[1]->rt,&kd1,&kd2,locVertex,locDOCA,var_doca);
	}

	// Vertex variables for filling histograms and cuts
	double locVertexZ = locVertex.Z();
	double locVertexR = locVertex.Perp();
	double locSlope = -6./20.;
	double locMinVertexR = 28.0 + locSlope*locVertexZ; // tight 30.25
	double locMaxVertexR = 34.0 + locSlope*locVertexZ; // tight 32

	// Conversion phi variable
	TVector3 negP3 = locParticlesMeasured[0]->lorentzMomentum().Vect();
	TVector3 posP3 = locParticlesMeasured[1]->lorentzMomentum().Vect();
	TVector3 uTot = posP3 + negP3;
	TVector3 u = locPairP4.Vect().Unit();
	TVector3 v = (posP3.Unit().Cross(negP3.Unit())).Unit();
	TVector3 z(0, 0, 1.);
	TVector3 w = (u.Cross(v)).Unit();
	TVector3 wc = (u.Cross(z)).Unit();
	double phiv = acos(w.Dot(wc));

	// only keep unique q+q- pairs since inclusive analysis
	if(dPreviousSourceObjects.find(locSourceObjects) == dPreviousSourceObjects.end())
		dPreviousSourceObjects.insert(locSourceObjects);
	else 
		return false;

	//Optional: FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// Fill any histograms here
		dHist_Vertex->Fill(locVertexZ, locVertexR);
		
		if(phiv > 3.0 && locDOCA < 3.0) {
			dHist_Vertex_DOCA->Fill(locVertexZ, locVertexR);
			dHist_DOCAVertex_DOCA->Fill(locDOCAVertex.Z(), locDOCAVertex.Perp());
			dHist_VertexZ_DOCAVertexZ->Fill(locDOCAVertex.Z(), locVertexZ);
			dHist_VertexR_DOCAVertexR->Fill(locDOCAVertex.Perp(), locVertexR);
			
			if(locPairP4.M() < 0.1)
				dHist_Vertex_LowMass->Fill(locVertexZ, locVertexR);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	//Require events to come from region dominated by conversions
	if(locVertexR < 2.) // exclude target
		return false;
	if(locVertexZ < 60. || locVertexZ > 100.) // start counter limit
		return false;
	if(locVertexZ < 80. && (locVertexR < 4. || locVertexR > 10.)) // straight section
		return false;
	if(locVertexZ > 80. && (locVertexR < locMinVertexR || locVertexR > locMaxVertexR)) // nose section
		return false;

	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHist_DOCA_DeltaTheta->Fill(locDeltaTheta, locDOCA);
		dHist_DOCA_DeltaPhi->Fill(locDeltaPhi, locDOCA);
		dHist_DOCA_PhiV->Fill(phiv, locDOCA);
		dHist_DOCA_DOCAVertex->Fill(locDOCAVertexVal, locDOCA);

		// Histogram mass after vertex requirement
		if(phiv > 3.0 && locDOCA < 3.0) {
			dHist_DeltaPhi_DeltaTheta->Fill(locDeltaTheta, locDeltaPhi);
			dHist_Mee->Fill(locPairP4.M());
			dHist_Vertex_Final->Fill(locVertexZ, locVertexR);
			dHist_EOverP->Fill(locEOverP[0], locEOverP[1]);
			if(locFCAL[0] && locFCAL[1]) 
				dHist_EOverP_BothFCAL->Fill(locEOverP[0], locEOverP[1]);
			if(locBCAL[0] && locBCAL[1]) 
				dHist_EOverP_BothBCAL->Fill(locEOverP[0], locEOverP[1]);
			if(locBCAL[0] && locFCAL[1])
				dHist_EOverP_FCALBCAL->Fill(locEOverP[0], locEOverP[1]);
			if(locBCAL[1] && locFCAL[0]) 
				dHist_EOverP_FCALBCAL->Fill(locEOverP[1], locEOverP[0]);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	// Require small DOCA between track pair
	if(locDOCA > 3.0) 
		return false;

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
