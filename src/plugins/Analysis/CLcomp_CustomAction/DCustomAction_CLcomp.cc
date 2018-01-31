// $Id$
//
//    File: DCustomAction_CLcomp.cc
// Created: Tue Jan 30 11:16:36 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#include "DCustomAction_CLcomp.h"

void DCustomAction_CLcomp::Initialize(JEventLoop* locEventLoop)
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
		// dMyHist = GetOrCreate_Histogram<TH1I>("MyHistName", "MyHistTitle", 100, 0.0, 1.0);
		dHist_IM = GetOrCreate_Histogram<TH1I>("IM", "Invariant Mass;Invariant Mass (GeV/c^{2})", 1000, 0.9, 2.4); 
		dHist_IM_kept = GetOrCreate_Histogram<TH1I>("IM_kept", "Invariant Mass of kept;Invariant Mass (GeV/c^{2})", 1000, 0.9, 2.4); 
		dHist_IM_cut = GetOrCreate_Histogram<TH1I>("IM_cut", "Invariant Mass of cut;Invariant Mass (GeV/c^{2})", 1000, 0.9, 2.4); 
		dHist_CL_KK = GetOrCreate_Histogram<TH1I>("CL_KK", "log10(CL_{KK});log10(CL_{KK})", 100, -50.0, 0.0);
		dHist_CL_KK_kept = GetOrCreate_Histogram<TH1I>("CL_KK_kept", "log10(CL_{KK}) kept after ratio cut;log10(CL_{KK})", 100, -50.0, 0.0);
		dHist_CL_KK_cut = GetOrCreate_Histogram<TH1I>("CL_KK_cut", "log10(CL_{KK}) cut;log10(CL_{KK})", 100, -50.0, 0.0);
		dHist_CL_pipi = GetOrCreate_Histogram<TH1I>("CL_pipi", "log10(CL_{#pi#pi});log10(CL_{#pi#pi})", 100, -50.0, 0.0);
		dHist_CL_pipi_kept = GetOrCreate_Histogram<TH1I>("CL_pipi_kept", "log10(CL_{#pi#pi}) kept after ratio cut;log10(CL_{#pi#pi})", 100, -50.0, 0.0);
		dHist_CL_pipi_cut = GetOrCreate_Histogram<TH1I>("CL_pipi_cut", "log10(CL_{#pi#pi}) cut;log10(CL_{#pi#pi})", 100, -50.0, 0.0);
		dHist_log10_ratio_vs_E = GetOrCreate_Histogram<TH2I>("log10_ratio_vs_E", "log10(CL_{KK}/CL_{#pi#pi}) vs beam energy;Beam Energy (GeV);Log10(CL_{KK}/CL_{#pi#pi})", 90, 3.0, 12.0, 100, -50.0, 50.0);
		dHist_log10_ratio_vs_mass = GetOrCreate_Histogram<TH2I>("log10_ratio_vs_mass", "log10(CL_{KK}/CL_{#pi#pi}) vs KK mass;IM_{KK} (GeV/c^{2});Log10(CL_{KK}/CL_{#pi#pi})", 1000, 0.9, 2.4, 100, -50.0, 50.0);
		dHist_mass_vs_E = GetOrCreate_Histogram<TH2I>("mass_vs_E", "KK mass vs beam energy;Beam Energy (GeV);IM_{KK} (GeV/c^{2})", 90, 3.0, 12.0, 300, 0.9, 2.4);
		dHist_mass_vs_E_kept = GetOrCreate_Histogram<TH2I>("mass_vs_E_kept", "KK mass vs beam energy for kept combos;Beam Energy (GeV);IM_{KK} (GeV/c^{2})", 90, 3.0, 12.0, 1000, 0.9, 2.4);
		dHist_mass_vs_E_cut = GetOrCreate_Histogram<TH2I>("mass_vs_E_cut", "KK mass vs beam energy for cut combos;Beam Energy (GeV);IM_{KK} (GeV/c^{2})", 90, 3.0, 12.0, 1000, 0.9, 2.4);

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	dAnalysisUtils = new DAnalysisUtilities(locEventLoop);
        dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
        dKinFitter = new DKinFitter(dKinFitUtils);
}

bool DCustomAction_CLcomp::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

	// Optional: Useful utility functions.
	// const DAnalysisUtilities* locAnalysisUtilities;
	// locEventLoop->GetSingle(locAnalysisUtilities);

	//Optional: check whether the user wanted to use the kinematic fit results when performing this action
	bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();
	if (!locUseKinFitResultsFlag)
		return true;

	// Get kinfit results from original reaction
	const DKinFitResults* locKinFitResults_KK = locParticleCombo->Get_KinFitResults();
        double locCL_KK = 0;
        if (locKinFitResults_KK == NULL) return false;
        locCL_KK = locKinFitResults_KK->Get_ConfidenceLevel();

	dKinFitter->Reset_NewEvent();

        //CREATE DKINFITPARTICLE OBJECTS FOR EACH PARTICLE
        	// First get the original hypothesis
        const DBeamPhoton* locBeamPhoton = static_cast<const DBeamPhoton*>(locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured());
        shared_ptr<DKinFitParticle> locKinFitParticle_BeamPhoton = dKinFitUtils->Make_BeamParticle(locBeamPhoton);
        shared_ptr<DKinFitParticle> locKinFitParticle_Target = dKinFitUtils->Make_TargetParticle(Proton);
	shared_ptr<DKinFitParticle> locKinFitParticle_KPlus = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(0));
        shared_ptr<DKinFitParticle> locKinFitParticle_KMinus = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(1));
	shared_ptr<DKinFitParticle> locKinFitParticle_Proton = dKinFitUtils->Make_DetectedParticle(locParticleCombo->Get_ParticleComboStep(0)->Get_FinalParticle_Measured(2)); // WARNING!!! This index depends on the order of your DReaction!!!

		// Get the tracks to set up the new hypotheses
	shared_ptr<DKinFitParticle> locKinFitParticle_PiPlus;
        shared_ptr<DKinFitParticle> locKinFitParticle_PiMinus;
        vector<const JObject*> FinalObjs = locParticleCombo->Get_FinalParticle_SourceObjects();
        for(size_t loc_i = 0; loc_i < FinalObjs.size() - 1; ++loc_i) // size() - 1 to get steps 0 and 1 which correspond to the kaon steps. Skip proton
	{
		auto particle = dynamic_cast<const DChargedTrack*>(FinalObjs[loc_i]);
                if(particle == nullptr) continue;

                if(particle->Get_Charge() > 0)
		{
			const DKinematicData *locKinData_PiPlus = dynamic_cast<const DKinematicData*>(particle->Get_Hypothesis(PiPlus));
                        if(locKinData_PiPlus == nullptr) continue;
                        locKinFitParticle_PiPlus = dKinFitUtils->Make_DetectedParticle(locKinData_PiPlus);
		}
		else if(particle->Get_Charge() < 0)
		{
			const DKinematicData *locKinData_PiMinus = dynamic_cast<const DKinematicData*>(particle->Get_Hypothesis(PiMinus));
                        if(locKinData_PiMinus == nullptr) continue;
			locKinFitParticle_PiMinus = dKinFitUtils->Make_DetectedParticle(locKinData_PiMinus);
		}
	}
        if (locKinFitParticle_PiPlus == nullptr || locKinFitParticle_PiMinus == nullptr) return true; // keep the event

        // SETUP THE CONSTRAINTS
        dKinFitter->Reset_NewFit(); //discards everything from the previous kinematic fit

        // P4 Constraint
        set<shared_ptr<DKinFitParticle>> locInitialParticles, locFinalParticles;
        locInitialParticles.insert(locKinFitParticle_BeamPhoton);
	locInitialParticles.insert(locKinFitParticle_Target);
        locFinalParticles.insert(locKinFitParticle_PiPlus);
	locFinalParticles.insert(locKinFitParticle_PiMinus);
	locFinalParticles.insert(locKinFitParticle_Proton);

	shared_ptr<DKinFitConstraint_P4> locP4Constraint = dKinFitUtils->Make_P4Constraint(locInitialParticles, locFinalParticles);
        TVector3 locMomentum = locKinFitParticle_BeamPhoton->Get_Momentum() + locKinFitParticle_Target->Get_Momentum();
        locMomentum -= (locKinFitParticle_PiPlus->Get_Momentum() + locKinFitParticle_PiMinus->Get_Momentum() + locKinFitParticle_Proton->Get_Momentum());
        locP4Constraint->Set_InitP3Guess(locMomentum);
        dKinFitter->Add_Constraint(locP4Constraint);

        // Vertex constraint
        set<shared_ptr<DKinFitParticle>> locFullConstrainParticles, locNoConstrainParticles;
        locFullConstrainParticles.insert(locKinFitParticle_PiPlus);
        locFullConstrainParticles.insert(locKinFitParticle_PiMinus);
        locFullConstrainParticles.insert(locKinFitParticle_Proton);
        locNoConstrainParticles.insert(locKinFitParticle_BeamPhoton);
	locNoConstrainParticles.insert(locKinFitParticle_Target);

        vector<shared_ptr<DKinFitParticle>> locVectFinalParticles;
        locVectFinalParticles.push_back(locKinFitParticle_PiPlus);
        locVectFinalParticles.push_back(locKinFitParticle_PiMinus);
        locVectFinalParticles.push_back(locKinFitParticle_Proton);
        TVector3 locVertexGuess = dAnalysisUtils->Calc_CrudeVertex(locVectFinalParticles);
        shared_ptr<DKinFitConstraint_Vertex> locDecayVertexConstraint = dKinFitUtils->Make_VertexConstraint(locFullConstrainParticles, locNoConstrainParticles, locVertexGuess);
        dKinFitter->Add_Constraint(locDecayVertexConstraint);

        // PERFORM THE KINEMATIC FIT
        dKinFitter->Set_DebugLevel(0); //e.g. 500 for everything
        dKinFitter->Fit_Reaction();

        // GET THE FIT RESULTS
        double locCL_pipi = dKinFitter->Get_ConfidenceLevel();

	// Get quantities for the histograms
	double IM_KK = (locKinFitParticle_KPlus->Get_P4() + locKinFitParticle_KMinus->Get_P4()).M();
	double E = locKinFitParticle_BeamPhoton->Get_P4().E();
	double CL_ratio = locCL_KK / locCL_pipi;
	double log10_ratio = log10( CL_ratio );

	//Optional: FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// These histograms are for EVERY combo, i.e. the same KK mass can be used multiple times!!!	

		// Fill any histograms here
		dHist_IM->Fill( IM_KK );
		dHist_CL_KK->Fill( log10(locCL_KK) );
		dHist_CL_pipi->Fill( log10(locCL_pipi) );
		dHist_log10_ratio_vs_E->Fill( E, log10_ratio );
		dHist_log10_ratio_vs_mass->Fill( IM_KK, log10_ratio );
		dHist_mass_vs_E->Fill( E, IM_KK);
		if( CL_ratio > 1)
		{
			dHist_IM_kept->Fill( IM_KK );
			dHist_CL_KK_kept->Fill( log10(locCL_KK) );
			dHist_CL_pipi_kept->Fill( log10(locCL_pipi) );
			dHist_mass_vs_E_kept->Fill( E, IM_KK);
		}
		else
		{
			dHist_IM_cut->Fill( IM_KK );
			dHist_CL_KK_cut->Fill( log10(locCL_KK) );
			dHist_CL_pipi_cut->Fill( log10(locCL_pipi) );
			dHist_mass_vs_E_cut->Fill( E, IM_KK);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return ( CL_ratio > 1);
	//return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
