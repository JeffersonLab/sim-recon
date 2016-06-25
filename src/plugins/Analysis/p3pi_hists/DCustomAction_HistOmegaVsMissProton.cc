// $Id$
//
//    File: DCustomAction_HistOmegaVsMissProton.cc
// Created: Sun Jun 28 22:48:32 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#include "DCustomAction_HistOmegaVsMissProton.h"

void DCustomAction_HistOmegaVsMissProton::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Optional: Useful utility functions.
		locEventLoop->GetSingle(dAnalysisUtilities);

		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		//	(Optional) Example: Create a histogram. 
			// This function will return the histogram if already created by another thread. If not pre-existing, it will create and return it. 
			// Function arguments are identical to those used for the histogram constructors
		string locHistTitle = ";#it{#gamma}#it{p}#rightarrow#it{#pi}^{+}#it{#pi}^{-}#it{#gamma}#it{#gamma} Missing Mass (GeV/c^{2})";
		locHistTitle += string(";#it{#pi}^{+}#it{#pi}^{-}#it{#gamma}#it{#gamma} Invariant Mass (GeV/c^{2})");
		dHist_OmegaVsMissProton = GetOrCreate_Histogram<TH2I>("OmegaVsMissProton", locHistTitle, 325, 0.3, 1.6, 300, 0.5, 1.1);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_HistOmegaVsMissProton::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//no duplicate entries: missing p4 is unique for each combo
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, false);
	DLorentzVector locOmegaP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 1, false);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHist_OmegaVsMissProton->Fill(locMissingP4.M(), locOmegaP4.M());
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
