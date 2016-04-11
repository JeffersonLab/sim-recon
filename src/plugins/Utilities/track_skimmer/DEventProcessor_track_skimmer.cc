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

	dIDXAStream_2track.open("eventlist_2track.idxa");
	dIDXAStream_2track << "IDXA" << endl;

	dIDXAStream_2track1pi0.open("eventlist_2track1pi0.idxa");
	dIDXAStream_2track1pi0 << "IDXA" << endl;

	dIDXAStream_3track.open("eventlist_3track.idxa");
	dIDXAStream_3track << "IDXA" << endl;

	dIDXAStream_3track1pi0.open("eventlist_3track1pi0.idxa");
	dIDXAStream_3track1pi0 << "IDXA" << endl;

	dIDXAStream_4track.open("eventlist_4track.idxa");
	dIDXAStream_4track << "IDXA" << endl;

	dIDXAStream_4track1pi0.open("eventlist_4track1pi0.idxa");
	dIDXAStream_4track1pi0 << "IDXA" << endl;

	dIDXAStream_5track.open("eventlist_5track.idxa");
	dIDXAStream_5track << "IDXA" << endl;

	dIDXAStream_5track1pi0.open("eventlist_5track1pi0.idxa");
	dIDXAStream_5track1pi0 << "IDXA" << endl;

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

	//Get the analysis results for all DReactions. 
		//Getting these objects triggers the analysis, if it wasn't performed already. 
		//These objects contain the DParticleCombo objects that survived the DAnalysisAction cuts that were added to the DReactions
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	//MUST LOCK AROUND MODIFICATION OF MEMBER VARIABLES IN brun() or evnt().
	LockState();
	{
		if(locChargedTracks.size() >= 2)
			dIDXAStream_2track << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
		if(locChargedTracks.size() >= 3)
			dIDXAStream_3track << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
		if(locChargedTracks.size() >= 4)
			dIDXAStream_4track << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
		if(locChargedTracks.size() >= 5)
			dIDXAStream_5track << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;

		// If a particle combo passed the cuts, save the event info in the output file
		for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
		{
			const DAnalysisResults* locAnalysisResults = locAnalysisResultsVector[loc_i];
			if(locAnalysisResults->Get_Reaction()->Get_ReactionName() != "pi0")
				continue; // analysis results were for a different reaction
			if(locAnalysisResults->Get_NumPassedParticleCombos() == 0)
				continue; // no combos passed

			if(locChargedTracks.size() >= 2)
				dIDXAStream_2track1pi0 << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
			if(locChargedTracks.size() >= 3)
				dIDXAStream_3track1pi0 << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
			if(locChargedTracks.size() >= 4)
				dIDXAStream_4track1pi0 << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
			if(locChargedTracks.size() >= 5)
				dIDXAStream_5track1pi0 << locRunNumber << " " << locEventNumber << " " << locUniqueID << endl;
		}
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
	dIDXAStream_2track.close();
	dIDXAStream_2track1pi0.close();

	dIDXAStream_3track.close();
	dIDXAStream_3track1pi0.close();

	dIDXAStream_4track.close();
	dIDXAStream_4track1pi0.close();

	dIDXAStream_5track.close();
	dIDXAStream_5track1pi0.close();

	return NOERROR;
}

