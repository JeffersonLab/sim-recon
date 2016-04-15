// $Id$
//
//    File: DEventProcessor_track_skimmer.cc
// Created: Tue Jan 13 11:08:16 EST 2015
// Creator: Paul (on Darwin Pauls-MacBook-Pro.local 14.0.0 i386)
//

#include "DEventProcessor_track_skimmer.h"

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_track_skimmer()); //register this plugin
		locApplication->AddFactoryGenerator(new DFactoryGenerator_track_skimmer()); //register the factory generator
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_track_skimmer::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... create historgrams or trees ...
	// japp->RootUnLock();
	//

	//CREATE STREAMS

	//# Tracks
	dIDXAStreamMap["q+"] = new ofstream("q+.idxa");
	dIDXAStreamMap["2q+"] = new ofstream("2q+.idxa");
	dIDXAStreamMap["3q+"] = new ofstream("3q+.idxa");
	dIDXAStreamMap["q-"] = new ofstream("q-.idxa");
	dIDXAStreamMap["2q-"] = new ofstream("2q-.idxa");

	//# Showers
	dIDXAStreamMap["q0"] = new ofstream("q0.idxa");
	dIDXAStreamMap["2q0"] = new ofstream("2q0.idxa");
	dIDXAStreamMap["3q0"] = new ofstream("3q0.idxa");
	dIDXAStreamMap["4q0"] = new ofstream("4q0.idxa");
	dIDXAStreamMap["5q0"] = new ofstream("5q0.idxa");
	dIDXAStreamMap["6q0"] = new ofstream("6q0.idxa");

	//# pi0s
	dIDXAStreamMap["pi0"] = new ofstream("pi0.idxa");
	dIDXAStreamMap["2pi0"] = new ofstream("2pi0.idxa");
	dIDXAStreamMap["3pi0"] = new ofstream("3pi0.idxa");

	//KShort
	dIDXAStreamMap["ks_2piq"] = new ofstream("ks_2piq.idxa");
	dIDXAStreamMap["ks_2pi0"] = new ofstream("ks_2pi0.idxa");

	//Eta
	dIDXAStreamMap["eta_2g"] = new ofstream("eta_2g.idxa");
	dIDXAStreamMap["eta_3pi0"] = new ofstream("eta_3pi0.idxa");
	dIDXAStreamMap["eta_3piq"] = new ofstream("eta_3piq.idxa");

	//Omega
	dIDXAStreamMap["omega_3piq"] = new ofstream("omega_3piq.idxa");
	dIDXAStreamMap["omega_pi0g"] = new ofstream("omega_pi0g.idxa");

	//EtaPrime
	dIDXAStreamMap["etaprm_2piqeta_2g"] = new ofstream("etaprm_2piqeta_2g.idxa");
	dIDXAStreamMap["etaprm_2piqeta_3pi0"] = new ofstream("etaprm_2piqeta_3pi0.idxa");
	dIDXAStreamMap["etaprm_2piqeta_3piq"] = new ofstream("etaprm_2piqeta_3piq.idxa");
	dIDXAStreamMap["etaprm_2piqg"] = new ofstream("etaprm_2piqg.idxa");
	dIDXAStreamMap["etaprm_2pi0eta_2g"] = new ofstream("etaprm_2pi0eta_2g.idxa");
	dIDXAStreamMap["etaprm_2pi0eta_3pi0"] = new ofstream("etaprm_2pi0eta_3pi0.idxa");
	dIDXAStreamMap["etaprm_2pi0eta_3piq"] = new ofstream("etaprm_2pi0eta_3piq.idxa");

	//Phi
	dIDXAStreamMap["phi_2kq"] = new ofstream("phi_2kq.idxa");
	dIDXAStreamMap["phi_3piq"] = new ofstream("phi_3piq.idxa");

	//Strange Baryons
	dIDXAStreamMap["lambda"] = new ofstream("lambda.idxa");
	dIDXAStreamMap["sigma0"] = new ofstream("sigma0.idxa");
	dIDXAStreamMap["sigma+"] = new ofstream("sigma+.idxa");
	dIDXAStreamMap["xi-"] = new ofstream("xi-.idxa");
	dIDXAStreamMap["xi0"] = new ofstream("xi0.idxa");

	//First line
	map<string, ofstream*>::iterator locStreamIterator = dIDXAStreamMap.begin();
	for(; locStreamIterator != dIDXAStreamMap.end(); ++locStreamIterator)
		*(locStreamIterator->second) << "IDXA" << endl;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_track_skimmer::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_track_skimmer::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software

	// See whether this is MC data or real data
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	unsigned int locRunNumber = locEventLoop->GetJEvent().GetRunNumber();
	unsigned int locUniqueID = locMCThrowns.empty() ? 1 : Get_FileNumber(locEventLoop);

	//Default Tag: No PreSelect!! // Definition of PreSelect may change, and these skims are intended to be universal
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);
	size_t locNumQPlus = 0, locNumQMinus = 0;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locChargedTracks[loc_i]->Contains_Charge(1))
			++locNumQPlus;
		if(locChargedTracks[loc_i]->Contains_Charge(-1)) //intentionally NOT "else if"!!! //can have a track with both charges!
			++locNumQMinus;
	}

	//Default Tag: No PreSelect!! // Definition of PreSelect may change, and these skims are intended to be universal
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	//Find which skims are satisfied by this event
	set<ofstream*> locSaveSkimStreams;

	//# Tracks
	if(locNumQPlus >= 1)
		locSaveSkimStreams.insert(dIDXAStreamMap["q+"]);
	if(locNumQPlus >= 2)
		locSaveSkimStreams.insert(dIDXAStreamMap["2q+"]);
	if(locNumQPlus >= 3)
		locSaveSkimStreams.insert(dIDXAStreamMap["3q+"]);
	if(locNumQMinus >= 1)
		locSaveSkimStreams.insert(dIDXAStreamMap["q-"]);
	if(locNumQMinus >= 2)
		locSaveSkimStreams.insert(dIDXAStreamMap["2q-"]);

	//# Showers
	if(locNeutralShowers.size() >= 1)
		locSaveSkimStreams.insert(dIDXAStreamMap["q0"]);
	if(locNeutralShowers.size() >= 2)
		locSaveSkimStreams.insert(dIDXAStreamMap["2q0"]);
	if(locNeutralShowers.size() >= 3)
		locSaveSkimStreams.insert(dIDXAStreamMap["3q0"]);
	if(locNeutralShowers.size() >= 4)
		locSaveSkimStreams.insert(dIDXAStreamMap["4q0"]);
	if(locNeutralShowers.size() >= 5)
		locSaveSkimStreams.insert(dIDXAStreamMap["5q0"]);
	if(locNeutralShowers.size() >= 6)
		locSaveSkimStreams.insert(dIDXAStreamMap["6q0"]);

	//Get the analysis results for all DReactions. 
		//Getting these objects triggers the analysis, if it wasn't performed already. 
		//These objects contain the DParticleCombo objects that survived the DAnalysisAction cuts that were added to the DReactions
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	//DECAYING PARTICLES
	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		const DAnalysisResults* locAnalysisResults = locAnalysisResultsVector[loc_i];
		if(locAnalysisResults->Get_NumPassedParticleCombos() == 0)
			continue; // no combos passed

		string locReactionName = locAnalysisResults->Get_Reaction()->Get_ReactionName();

		//# Pi0
		if(locReactionName == "skim_pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["pi0"]);
		else if(locReactionName == "skim_2pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["2pi0"]);
		else if(locReactionName == "skim_3pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["3pi0"]);

		//KShort
		else if(locReactionName == "skim_ks_2piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["ks_2piq"]);
		else if(locReactionName == "skim_ks_2pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["ks_2pi0"]);

		//Eta
		else if(locReactionName == "skim_eta_2g")
			locSaveSkimStreams.insert(dIDXAStreamMap["eta_2g"]);
		else if(locReactionName == "skim_eta_3pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["eta_3pi0"]);
		else if(locReactionName == "skim_eta_3piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["eta_3piq"]);

		//Omega
		else if(locReactionName == "skim_omega_3piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["omega_3piq"]);
		else if(locReactionName == "skim_omega_pi0g")
			locSaveSkimStreams.insert(dIDXAStreamMap["omega_pi0g"]);

		//EtaPrime
		else if(locReactionName == "skim_etaprm_2piqeta_2g")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2piqeta_2g"]);
		else if(locReactionName == "skim_etaprm_2piqeta_3pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2piqeta_3pi0"]);
		else if(locReactionName == "skim_etaprm_2piqeta_3piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2piqeta_3piq"]);
		else if(locReactionName == "skim_etaprm_2piqg")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2piqg"]);
		else if(locReactionName == "skim_etaprm_2pi0eta_2g")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2pi0eta_2g"]);
		else if(locReactionName == "skim_etaprm_2pi0eta_3pi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2pi0eta_3pi0"]);
		else if(locReactionName == "skim_etaprm_2pi0eta_3piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["etaprm_2pi0eta_3piq"]);

		//Phi
		else if(locReactionName == "skim_phi_2kq")
			locSaveSkimStreams.insert(dIDXAStreamMap["phi_2kq"]);
		else if(locReactionName == "skim_phi_3piq")
			locSaveSkimStreams.insert(dIDXAStreamMap["phi_3piq"]);

		//Strange Baryons
		else if(locReactionName == "skim_lambda")
			locSaveSkimStreams.insert(dIDXAStreamMap["lambda"]);
		else if(locReactionName == "skim_sigma0")
			locSaveSkimStreams.insert(dIDXAStreamMap["sigma0"]);
		else if(locReactionName == "skim_sigma+")
			locSaveSkimStreams.insert(dIDXAStreamMap["sigma+"]);
		else if(locReactionName == "skim_xi-")
			locSaveSkimStreams.insert(dIDXAStreamMap["xi-"]);
		else if(locReactionName == "skim_xi0")
			locSaveSkimStreams.insert(dIDXAStreamMap["xi0"]);
	}

	//BUILD TO-PRINT STRING
	ostringstream locToPrintStream;
	locToPrintStream << locRunNumber << " " << locEventNumber << " " << locUniqueID;
	string locToPrintString = locToPrintStream.str();

	//LOCK AND PRINT
	set<ofstream*>::iterator locSetIterator = locSaveSkimStreams.begin();
	LockState();
	{
		for(; locSetIterator != locSaveSkimStreams.end(); ++locSetIterator)
			(**locSetIterator) << locToPrintString << endl;
	}
	UnlockState();

	return NOERROR;
}

int DEventProcessor_track_skimmer::Get_FileNumber(JEventLoop* locEventLoop) const
{
	//Assume that the file name is in the format: *_X.ext, where:
		//X is the file number (a string of numbers of any length)
		//ext is the file extension (probably .evio or .hddm)

	//get the event source
	JEventSource* locEventSource = locEventLoop->GetJEvent().GetJEventSource();
	if(locEventSource == NULL)
		return -1;

	//get the source file name (strip the path)
	string locSourceFileName = locEventSource->GetSourceName();

	//find the last "_" & "." indices
	size_t locUnderscoreIndex = locSourceFileName.rfind("_");
	size_t locDotIndex = locSourceFileName.rfind(".");
	if((locUnderscoreIndex == string::npos) || (locDotIndex == string::npos))
		return -1;

	size_t locNumberLength = locDotIndex - locUnderscoreIndex - 1;
	string locFileNumberString = locSourceFileName.substr(locUnderscoreIndex + 1, locNumberLength);

	int locFileNumber = -1;
	istringstream locFileNumberStream(locFileNumberString);
	locFileNumberStream >> locFileNumber;

	return locFileNumber;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_track_skimmer::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_track_skimmer::fini(void)
{
	// Called before program exit after event processing is finished.

	//Close streams
	map<string, ofstream*>::iterator locStreamIterator = dIDXAStreamMap.begin();
	for(; locStreamIterator != dIDXAStreamMap.end(); ++locStreamIterator)
		locStreamIterator->second->close();

	return NOERROR;
}

