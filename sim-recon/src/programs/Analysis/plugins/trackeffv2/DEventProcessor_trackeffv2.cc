#include <iostream>
#include <cmath>
using namespace std;

#include "DEventProcessor_trackeffv2.h"

bool gHistWireBasedFlag = true;

// Routine used to create our DEventProcessor
extern "C"
{
	void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new DEventProcessor_trackeffv2());
	}
} // "C"

void DEventProcessor_trackeffv2::MakeAndChangeTo_SubDirectory(string locDirName)
{
	TDirectoryFile* locSubDirectory = static_cast<TDirectoryFile*>(gDirectory->Get(locDirName.c_str()));
	if(locSubDirectory == NULL) //else folder already created earlier
		locSubDirectory = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
	locSubDirectory->cd();
}

//------------------
// init
//------------------
jerror_t DEventProcessor_trackeffv2::init(void)
{
	//MODIFY THE BELOW VARIABLES TO CHANGE THE HISTOGRAM PARTICLES & BINNING
	dFinalStatePIDs.push_back(Proton);  dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KPlus);

	dMinimumMatchFOM = -1.0;

	dKinematicsHistBinRange_MinP = 0.1;
	dKinematicsHistBinRange_MaxP = 1.0;
	dKinematicsHistBins_NumPBins = 9;
	double locPBinSize = (dKinematicsHistBinRange_MaxP - dKinematicsHistBinRange_MinP)/(double(dKinematicsHistBins_NumPBins));

	dKinematicsHistBinRange_MinTheta = 0.0;
	dKinematicsHistBinRange_MaxTheta = 140.0;
	dKinematicsHistBins_NumThetaBins = 14;
	double locThetaBinSize = (dKinematicsHistBinRange_MaxTheta - dKinematicsHistBinRange_MinTheta)/(double(dKinematicsHistBins_NumThetaBins));

	dEfficiencyHists_MinP = 0.0;
	dEfficiencyHists_MaxP = 3.0;
	dEfficiencyHists_NumPBins = 60;

	dEfficiencyHists_MinTheta = 0.0;
	dEfficiencyHists_MaxTheta = 140.0;
	dEfficiencyHists_NumThetaBins = 28;


	// Histogram Creation
	ostringstream locHistNameStream, locHistTitleStream;
	string locHistName, locHistTitle;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!! //This function will be called only one time, so it is safe to assume that none of the hists have been created yet

	locHistName = "NumCandidates";
	locHistTitle = "# Candidates";
	dHist_NumCandidates = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 21, -0.5, 20.5);

	if(gHistWireBasedFlag)
	{
		locHistName = "NumWireBased";
		locHistTitle = "# Wire-Based Particles";
		dHist_NumWireBased = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 21, -0.5, 20.5);
	}

	locHistName = "NumTimeBased";
	locHistTitle = "# Time-Based Particles";
	dHist_NumTimeBased = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 21, -0.5, 20.5);

	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
	{
		Particle_t locPID = dFinalStatePIDs[loc_i];
		string locParticleName = ParticleType(locPID);
		string locParticleROOTName = ParticleName_ROOT(locPID);

		gDirectory->cd("/");
		MakeAndChangeTo_SubDirectory(locParticleName);

		//setup counts maps
		deque<unsigned int> locTempIntDeque(dEfficiencyHists_NumPBins, 0);
		dNumThrown[locPID].assign(dEfficiencyHists_NumThetaBins, locTempIntDeque);
		dNumTimesClose_Candidates[locPID].assign(dEfficiencyHists_NumThetaBins, locTempIntDeque);
		if(gHistWireBasedFlag)
			dNumTimesClose_WireBased[locPID].assign(dEfficiencyHists_NumThetaBins, locTempIntDeque);
		dNumTimesClose_TimeBased[locPID].assign(dEfficiencyHists_NumThetaBins, locTempIntDeque);

		//setup hist maps
		deque<TH1D*> locTempHistDeque;
		locTempHistDeque.resize(dKinematicsHistBins_NumPBins);
		dHistMap_DeltaPOverP_Candidates[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaTheta_Candidates[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaPhi_Candidates[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaVertexZ_Candidates[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		if(gHistWireBasedFlag)
		{
			dHistMap_DeltaPOverP_WireBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
			dHistMap_DeltaTheta_WireBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
			dHistMap_DeltaPhi_WireBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
			dHistMap_DeltaVertexZ_WireBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		}
		dHistMap_DeltaPOverP_TimeBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaTheta_TimeBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaPhi_TimeBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);
		dHistMap_DeltaVertexZ_TimeBased[locPID].assign(dKinematicsHistBins_NumThetaBins, locTempHistDeque);

		// create counts & efficiency hists
		locHistName = "NumThrown";
		locHistTitle = string("# Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
		dHist_NumThrown[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dEfficiencyHists_NumThetaBins, dEfficiencyHists_MinTheta, dEfficiencyHists_MaxTheta, dEfficiencyHists_NumPBins, dEfficiencyHists_MinP, dEfficiencyHists_MaxP);

		locHistName = "Efficiencies_Candidates";
		locHistTitle = locParticleROOTName + string(" Track Candidate Efficiencies;#theta#circ;p (GeV/c)");
		dHist_CloseEfficiencies_Candidates[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dEfficiencyHists_NumThetaBins, dEfficiencyHists_MinTheta, dEfficiencyHists_MaxTheta, dEfficiencyHists_NumPBins, dEfficiencyHists_MinP, dEfficiencyHists_MaxP);

		if(gHistWireBasedFlag)
		{
			locHistName = "Efficiencies_WireBased";
			locHistTitle = locParticleROOTName + string(" Wire-Based Track Efficiencies;#theta#circ;p (GeV/c)");
			dHist_CloseEfficiencies_WireBased[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dEfficiencyHists_NumThetaBins, dEfficiencyHists_MinTheta, dEfficiencyHists_MaxTheta, dEfficiencyHists_NumPBins, dEfficiencyHists_MinP, dEfficiencyHists_MaxP);
		}

		locHistName = "Efficiencies_TimeBased";
		locHistTitle = locParticleROOTName + string(" Time-Based Track Efficiencies;#theta#circ;p (GeV/c)");
		dHist_CloseEfficiencies_TimeBased[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dEfficiencyHists_NumThetaBins, dEfficiencyHists_MinTheta, dEfficiencyHists_MaxTheta, dEfficiencyHists_NumPBins, dEfficiencyHists_MinP, dEfficiencyHists_MaxP);

		//create hists
		for(unsigned int locThetaBin = 0; locThetaBin < dKinematicsHistBins_NumThetaBins; ++locThetaBin)
		{
			double locMinTheta = dKinematicsHistBinRange_MinTheta + (double(locThetaBin))*locThetaBinSize;
			double locMaxTheta = locMinTheta + locThetaBinSize;
			for(unsigned int locPBin = 0; locPBin < dKinematicsHistBins_NumPBins; ++locPBin)
			{
				double locMinP = dKinematicsHistBinRange_MinP + (double(locPBin))*locPBinSize;
				double locMaxP = locMinP + locPBinSize;

				locHistTitleStream.str("");
				locHistTitleStream << locParticleROOTName << ", #theta #in [" << locMinTheta << "#circ, " << locMaxTheta << "#circ), p #in [" << locMinP << ", " << locMaxP << ") GeV/c";

				MakeAndChangeTo_SubDirectory("Candidates");
				MakeAndChangeTo_SubDirectory("DeltaPOverP");

				// Candidates: DeltaP/P
				locHistNameStream.str("");
				locHistNameStream << "DeltaPOverP_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Candidate ") + locHistTitleStream.str() + string(";#Deltap/p (Reconstructed - Thrown)");
				(dHistMap_DeltaPOverP_Candidates[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 500, -1.0, 1.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaTheta");

				// Candidates: DeltaTheta
				locHistNameStream.str("");
				locHistNameStream << "DeltaTheta_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Candidate ") + locHistTitleStream.str() + string(";#Delta#theta#circ (Reconstructed - Thrown)");
				(dHistMap_DeltaTheta_Candidates[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaPhi");

				// Candidates: DeltaPhi
				locHistNameStream.str("");
				locHistNameStream << "DeltaPhi_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Candidate ") + locHistTitleStream.str() + string(";#Delta#phi#circ (Reconstructed - Thrown)");
				(dHistMap_DeltaPhi_Candidates[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaVertexZ");

				// Candidates: DeltaVertexZ
				locHistNameStream.str("");
				locHistNameStream << "DeltaVertexZ_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Candidate ") + locHistTitleStream.str() + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				(dHistMap_DeltaVertexZ_Candidates[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 500, -100.0, 100.0);

				if(gHistWireBasedFlag)
				{
					gDirectory->cd("..");
					gDirectory->cd("..");
					MakeAndChangeTo_SubDirectory("WireBased");
					MakeAndChangeTo_SubDirectory("DeltaPOverP");

					// Wire-Based: DeltaP/P
					locHistNameStream.str("");
					locHistNameStream << "DeltaPOverP_" << locThetaBin << "_" << locPBin;
					locHistTitle = string("Wire-Based ") + locHistTitleStream.str() + string(";#Deltap/p (Reconstructed - Thrown)");
					(dHistMap_DeltaPOverP_WireBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 1000, -1.0, 1.0);

					gDirectory->cd("..");
					MakeAndChangeTo_SubDirectory("DeltaTheta");

					// Wire-Based: DeltaTheta
					locHistNameStream.str("");
					locHistNameStream << "DeltaTheta_" << locThetaBin << "_" << locPBin;
					locHistTitle = string("Wire-Based ") + locHistTitleStream.str() + string(";#Delta#theta#circ (Reconstructed - Thrown)");
					(dHistMap_DeltaTheta_WireBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

					gDirectory->cd("..");
					MakeAndChangeTo_SubDirectory("DeltaPhi");

					// Wire-Based: DeltaPhi
					locHistNameStream.str("");
					locHistNameStream << "DeltaPhi_" << locThetaBin << "_" << locPBin;
					locHistTitle = string("Wire-Based ") + locHistTitleStream.str() + string(";#Delta#phi#circ (Reconstructed - Thrown)");
					(dHistMap_DeltaPhi_WireBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

					gDirectory->cd("..");
					MakeAndChangeTo_SubDirectory("DeltaVertexZ");

					// Wire-Based: DeltaVertexZ
					locHistNameStream.str("");
					locHistNameStream << "DeltaVertexZ_" << locThetaBin << "_" << locPBin;
					locHistTitle = string("Wire-Based ") + locHistTitleStream.str() + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
					(dHistMap_DeltaVertexZ_WireBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 500, -100.0, 100.0);
				}

				gDirectory->cd("..");
				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("TimeBased");
				MakeAndChangeTo_SubDirectory("DeltaPOverP");

				// Time-Based: DeltaP/P
				locHistNameStream.str("");
				locHistNameStream << "DeltaPOverP_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Time-Based ") + locHistTitleStream.str() + string(";#Deltap/p (Reconstructed - Thrown)");
				(dHistMap_DeltaPOverP_TimeBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 500, -1.0, 1.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaTheta");

				// Time-Based: DeltaTheta
				locHistNameStream.str("");
				locHistNameStream << "DeltaTheta_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Time-Based ") + locHistTitleStream.str() + string(";#Delta#theta#circ (Reconstructed - Thrown)");
				(dHistMap_DeltaTheta_TimeBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaPhi");

				// Time-Based: DeltaPhi
				locHistNameStream.str("");
				locHistNameStream << "DeltaPhi_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Time-Based ") + locHistTitleStream.str() + string(";#Delta#phi#circ (Reconstructed - Thrown)");
				(dHistMap_DeltaPhi_TimeBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 2000, -180.0, 180.0);

				gDirectory->cd("..");
				MakeAndChangeTo_SubDirectory("DeltaVertexZ");

				// Time-Based: DeltaVertexZ
				locHistNameStream.str("");
				locHistNameStream << "DeltaVertexZ_" << locThetaBin << "_" << locPBin;
				locHistTitle = string("Time-Based ") + locHistTitleStream.str() + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				(dHistMap_DeltaVertexZ_TimeBased[locPID])[locThetaBin][locPBin] = new TH1D(locHistNameStream.str().c_str(), locHistTitle.c_str(), 500, -100.0, 100.0);

				gDirectory->cd("..");
				gDirectory->cd("..");
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trackeffv2::brun(JEventLoop* locEventLoop, int runnumber)
{
	dMCThrownMatchingFactory = static_cast<DMCThrownMatching_factory*>(locEventLoop->GetFactory("DMCThrownMatching"));
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trackeffv2::evnt(JEventLoop* locEventLoop, int eventnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	vector<const DTrackCandidate*> locTrackCandidates;
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	
	locEventLoop->Get(locMCThrowns, "FinalState");
	locEventLoop->Get(locTrackCandidates);
	if(gHistWireBasedFlag)
		locEventLoop->Get(locTrackWireBasedVector);
	locEventLoop->Get(locTrackTimeBasedVector);

	// Thrown Matching: Candidates
	vector<const DKinematicData*> locInputKinematicDataVector;
	locInputKinematicDataVector.clear();
	map<const DKinematicData*, const DMCThrown*> locDataToThrownMap_Candidates;
	map<const DMCThrown*, const DKinematicData*> locThrownToDataMap_Candidates;
	for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
		locInputKinematicDataVector.push_back(locTrackCandidates[loc_i]);
	Find_GenReconMatches(locMCThrowns, locInputKinematicDataVector, locDataToThrownMap_Candidates, locThrownToDataMap_Candidates, false);

	// Thrown Matching: Wire-Based
	locInputKinematicDataVector.clear();
	map<const DKinematicData*, const DMCThrown*> locDataToThrownMap_WireBased;
	map<const DMCThrown*, const DKinematicData*> locThrownToDataMap_WireBased;
	if(gHistWireBasedFlag)
	{
		for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
			locInputKinematicDataVector.push_back(locTrackWireBasedVector[loc_i]);
		Find_GenReconMatches(locMCThrowns, locInputKinematicDataVector, locDataToThrownMap_WireBased, locThrownToDataMap_WireBased, true);
	}

	// Thrown Matching: Time-Based
	locInputKinematicDataVector.clear();
	map<const DKinematicData*, const DMCThrown*> locDataToThrownMap_TimeBased;
	map<const DMCThrown*, const DKinematicData*> locThrownToDataMap_TimeBased;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		locInputKinematicDataVector.push_back(locTrackTimeBasedVector[loc_i]);
	Find_GenReconMatches(locMCThrowns, locInputKinematicDataVector, locDataToThrownMap_TimeBased, locThrownToDataMap_TimeBased, true);

	// # tracks at each stage
	set<unsigned long> locUniqueCandidateIDs_WireBased;
	if(gHistWireBasedFlag)
	{
		for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
			locUniqueCandidateIDs_WireBased.insert(locTrackWireBasedVector[loc_i]->candidateid);
	}
	set<unsigned long> locUniqueCandidateIDs_TimeBased;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		locUniqueCandidateIDs_TimeBased.insert(locTrackTimeBasedVector[loc_i]->candidateid);

	//Fill histograms
	japp->RootWriteLock();
	{
		//counts
		dHist_NumCandidates->Fill(locTrackCandidates.size());
		if(gHistWireBasedFlag)
			dHist_NumWireBased->Fill(locUniqueCandidateIDs_WireBased.size());
		dHist_NumTimeBased->Fill(locUniqueCandidateIDs_TimeBased.size());

		//Fill Delta Hists
		double locPBinSize_Kinematics = (dKinematicsHistBinRange_MaxP - dKinematicsHistBinRange_MinP)/(double(dKinematicsHistBins_NumPBins));
		double locThetaBinSize_Kinematics = (dKinematicsHistBinRange_MaxTheta - dKinematicsHistBinRange_MinTheta)/(double(dKinematicsHistBins_NumThetaBins));
		double locPBinSize_Efficiencies = (dEfficiencyHists_MaxP - dEfficiencyHists_MinP)/(double(dEfficiencyHists_NumPBins));
		double locThetaBinSize_Efficiencies = (dEfficiencyHists_MaxTheta - dEfficiencyHists_MinTheta)/(double(dEfficiencyHists_NumThetaBins));
		map<Particle_t, size_t> locNumCloseCandidates;
		map<Particle_t, size_t> locNumCloseWireBased;
		map<Particle_t, size_t> locNumCloseTimeBased;
		for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		{
			const DMCThrown* locMCThrown = locMCThrowns[loc_i];
			Particle_t locPID = locMCThrown->PID();
			if(dHistMap_DeltaPOverP_Candidates.find(locPID) == dHistMap_DeltaPOverP_Candidates.end())
				continue; //not histogramming info for this PID

			DVector3 locThrownP3 = locMCThrown->momentum();
			double locThrownP = locThrownP3.Mag();
			double locThrownTheta = locThrownP3.Theta()*180.0/TMath::Pi();
			double locThrownPhi = locThrownP3.Phi()*180.0/TMath::Pi();

			// Determine Kinematics P & Theta Bins
			bool locKinematicsBinsOKFlag = true;
			if((locThrownP < dKinematicsHistBinRange_MinP) || (locThrownP > dKinematicsHistBinRange_MaxP))
				locKinematicsBinsOKFlag = false;
			int locPBin = int((locThrownP - dKinematicsHistBinRange_MinP)/locPBinSize_Kinematics);
			if((locThrownTheta < dKinematicsHistBinRange_MinTheta) || (locThrownTheta > dKinematicsHistBinRange_MaxTheta))
				locKinematicsBinsOKFlag = false;
			int locThetaBin = int((locThrownTheta - dKinematicsHistBinRange_MinTheta)/locThetaBinSize_Kinematics);

			// Determine Efficiency P & Theta Bins
			bool locEfficiencyBinsOKFlag = true;
			if((locThrownP < dEfficiencyHists_MinP) || (locThrownP > dEfficiencyHists_MaxP))
				locEfficiencyBinsOKFlag = false;
			int locPBin_Efficiencies = int((locThrownP - dEfficiencyHists_MinP)/locPBinSize_Efficiencies);
			if((locThrownTheta < dEfficiencyHists_MinTheta) || (locThrownTheta > dEfficiencyHists_MaxTheta))
				locEfficiencyBinsOKFlag = false;
			int locThetaBin_Efficiencies = int((locThrownTheta - dEfficiencyHists_MinTheta)/locThetaBinSize_Efficiencies);

			// Add to num thrown
			if(locEfficiencyBinsOKFlag)
				++(dNumThrown[locPID])[locThetaBin_Efficiencies][locPBin_Efficiencies];

			// Hist Candidates
			if(locThrownToDataMap_Candidates.find(locMCThrown) != locThrownToDataMap_Candidates.end())
			{
				const DKinematicData* locKinematicData = locThrownToDataMap_Candidates[locMCThrown];
				if(Check_IsCloseMatch(locKinematicData, locMCThrown, true) && locEfficiencyBinsOKFlag)
					++(dNumTimesClose_Candidates[locPID])[locThetaBin_Efficiencies][locPBin_Efficiencies];

				if(locKinematicsBinsOKFlag)
				{
					double locDeltaPFraction = (locKinematicData->momentum().Mag() - locThrownP)/locThrownP;
					double locDeltaTheta = locKinematicData->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
					double locDeltaPhi = locKinematicData->momentum().Phi()*180.0/TMath::Pi() - locThrownPhi;
					if(locDeltaPhi > 180.0)
						locDeltaPhi -= 360.0;
					if(locDeltaPhi < -180.0)
						locDeltaPhi += 360.0;
					double locDeltaZ = locKinematicData->position().Z() - locMCThrown->position().Z();
					(dHistMap_DeltaPOverP_Candidates[locPID])[locThetaBin][locPBin]->Fill(locDeltaPFraction);
					(dHistMap_DeltaTheta_Candidates[locPID])[locThetaBin][locPBin]->Fill(locDeltaTheta);
					(dHistMap_DeltaPhi_Candidates[locPID])[locThetaBin][locPBin]->Fill(locDeltaPhi);
					(dHistMap_DeltaVertexZ_Candidates[locPID])[locThetaBin][locPBin]->Fill(locDeltaZ);
				}
			}

			// Hist WireBased
			if(gHistWireBasedFlag && (locThrownToDataMap_WireBased.find(locMCThrown) != locThrownToDataMap_WireBased.end()))
			{
				const DKinematicData* locKinematicData = locThrownToDataMap_WireBased[locMCThrown];
				if(Check_IsCloseMatch(locKinematicData, locMCThrown, false) && locEfficiencyBinsOKFlag)
					++(dNumTimesClose_WireBased[locPID])[locThetaBin_Efficiencies][locPBin_Efficiencies];

				if(locKinematicsBinsOKFlag)
				{
					double locDeltaPFraction = (locKinematicData->momentum().Mag() - locThrownP)/locThrownP;
					double locDeltaTheta = locKinematicData->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
					double locDeltaPhi = locKinematicData->momentum().Phi()*180.0/TMath::Pi() - locThrownPhi;
					if(locDeltaPhi > 180.0)
						locDeltaPhi -= 360.0;
					if(locDeltaPhi < -180.0)
						locDeltaPhi += 360.0;
					double locDeltaZ = locKinematicData->position().Z() - locMCThrown->position().Z();
					(dHistMap_DeltaPOverP_WireBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaPFraction);
					(dHistMap_DeltaTheta_WireBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaTheta);
					(dHistMap_DeltaPhi_WireBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaPhi);
					(dHistMap_DeltaVertexZ_WireBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaZ);
				}
			}

			// Hist TimeBased
			if(locThrownToDataMap_TimeBased.find(locMCThrown) != locThrownToDataMap_TimeBased.end())
			{
				const DKinematicData* locKinematicData = locThrownToDataMap_TimeBased[locMCThrown];
				if(Check_IsCloseMatch(locKinematicData, locMCThrown, false) && locEfficiencyBinsOKFlag)
					++(dNumTimesClose_TimeBased[locPID])[locThetaBin_Efficiencies][locPBin_Efficiencies];

				if(locKinematicsBinsOKFlag)
				{
					double locDeltaPFraction = (locKinematicData->momentum().Mag() - locThrownP)/locThrownP;
					double locDeltaTheta = locKinematicData->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
					double locDeltaPhi = locKinematicData->momentum().Phi()*180.0/TMath::Pi() - locThrownPhi;
					if(locDeltaPhi > 180.0)
						locDeltaPhi -= 360.0;
					if(locDeltaPhi < -180.0)
						locDeltaPhi += 360.0;
					double locDeltaZ = locKinematicData->position().Z() - locMCThrown->position().Z();
					(dHistMap_DeltaPOverP_TimeBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaPFraction);
					(dHistMap_DeltaTheta_TimeBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaTheta);
					(dHistMap_DeltaPhi_TimeBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaPhi);
					(dHistMap_DeltaVertexZ_TimeBased[locPID])[locThetaBin][locPBin]->Fill(locDeltaZ);
				}
			}
		}
	}
	japp->RootUnLock();

	return NOERROR;
}

bool DEventProcessor_trackeffv2::Check_IsCloseMatch(const DKinematicData* locKinematicData, const DMCThrown* locMCThrown, bool locIsCandidateFlag)
{
	//ignores PID and charge
	DVector3 locThrownP3 = locMCThrown->momentum();
	double locThrownP = locThrownP3.Mag();

	double locDeltaPFraction = (locKinematicData->momentum().Mag() - locThrownP)/locThrownP;
	double locMaxDeltaPFraction = locIsCandidateFlag ? 0.2 : 0.1;
	if(fabs(locDeltaPFraction) > locMaxDeltaPFraction)
		return false;

	double locThrownTheta = locThrownP3.Theta()*180.0/TMath::Pi();
	double locDeltaTheta = locKinematicData->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
	double locMaxDeltaTheta = locIsCandidateFlag ? 40.0 : 15.0;
	if(fabs(locDeltaTheta) > locMaxDeltaTheta)
		return false;

	double locDeltaPhi = fabs(locKinematicData->momentum().Phi() - locThrownP3.Phi());
	if(locDeltaPhi > TMath::Pi())
		locDeltaPhi -= 2.0*TMath::Pi();
	locDeltaPhi *= 180.0/TMath::Pi();
	double locMaxDeltaPhi = locIsCandidateFlag ? 40.0 : 15.0;
	if(locIsCandidateFlag && (locThrownTheta < 5.0))
		locMaxDeltaPhi = 180.0; //phase-space is so tight, ignore phi here
	if(fabs(locDeltaPhi) > locMaxDeltaPhi)
		return false;

	double locDeltaZ = locKinematicData->position().Z() - locMCThrown->position().Z();
	double locMaxDeltaZ = locIsCandidateFlag ? 100000.0 : 10.0;
	if(fabs(locDeltaZ) > locMaxDeltaZ)
		return false;

	return true;
}

void DEventProcessor_trackeffv2::Find_GenReconMatches(const vector<const DMCThrown*>& locInputMCThrowns, const vector<const DKinematicData*>& locInputKinematicDataVector, map<const DKinematicData*, const DMCThrown*>& locDataToThrownMap, map<const DMCThrown*, const DKinematicData*>& locThrownToDataMap, bool locRequirePIDMatchFlag)
{
	// Loop over all pairs of thrown and reconstructed tracks, finding the best match and setting it aside (and saving it in the maps).
		// Do this until either all throwns or all reconstructed particles have been matched
	const DKinematicData* locKinematicData;
	const DMCThrown* locMCThrown;
	size_t locBestDataIndex = 0, locBestMCThrownIndex = 0;
	double locMatchFOM, locBestMatchFOM;
	vector<const DMCThrown*> locMCThrownVector = locInputMCThrowns;
	vector<const DKinematicData*> locKinematicDataVector = locInputKinematicDataVector;
	while((!locMCThrownVector.empty()) && (!locKinematicDataVector.empty()))
	{
		bool locMatchFoundFlag = false;
		locBestMatchFOM = dMinimumMatchFOM;
		for(size_t loc_i = 0; loc_i < locMCThrownVector.size(); ++loc_i)
		{
			locMCThrown = locMCThrownVector[loc_i];
			for(size_t loc_j = 0; loc_j < locKinematicDataVector.size(); ++loc_j)
			{
				locKinematicData = locKinematicDataVector[loc_j];
				if(locRequirePIDMatchFlag)
				{
					if(locKinematicData->PID() != locMCThrown->PID())
						continue; //wrong pid
				}
				locMatchFOM = dMCThrownMatchingFactory->Calc_MatchFOM(locMCThrown->momentum(), locKinematicData->momentum());

				if(locMatchFOM >= locBestMatchFOM)
				{
					locMatchFoundFlag = true;
					locBestMatchFOM = locMatchFOM;
					locBestMCThrownIndex = loc_i;
					locBestDataIndex = loc_j;
				}
			}
		}

		if(!locMatchFoundFlag) //no more good matches!
			break;

		// save match to maps
		locMCThrown = locMCThrownVector[locBestMCThrownIndex];
		locKinematicData = locKinematicDataVector[locBestDataIndex];

		locDataToThrownMap[locKinematicData] = locMCThrown;
		locThrownToDataMap[locMCThrown] = locKinematicData;
		locKinematicDataVector.erase(locKinematicDataVector.begin() + locBestDataIndex);
		locMCThrownVector.erase(locMCThrownVector.begin() + locBestMCThrownIndex);
	}
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trackeffv2::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trackeffv2::fini(void)
{
	// Calculate efficiencies for each particle in each theta/p bin
	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
	{
		Particle_t locPID = dFinalStatePIDs[loc_i];
		unsigned int locEntries_Thrown = 0;
		unsigned int locEntries_Candidates = 0;
		unsigned int locEntries_WireBased = 0;
		unsigned int locEntries_TimeBased = 0;
		for(unsigned int locThetaBin = 0; locThetaBin < dEfficiencyHists_NumThetaBins; ++locThetaBin)
		{
			for(unsigned int locPBin = 0; locPBin < dEfficiencyHists_NumPBins; ++locPBin)
			{
				unsigned int locNumThrown = (dNumThrown[locPID])[locThetaBin][locPBin];
				if(locNumThrown == 0)
					continue;
				locEntries_Thrown += locNumThrown;
				dHist_NumThrown[locPID]->SetBinContent(locThetaBin + 1, locPBin + 1, locNumThrown);

				locEntries_Candidates += (dNumTimesClose_Candidates[locPID])[locThetaBin][locPBin];
				unsigned int locNumClose_Candidates = (dNumTimesClose_Candidates[locPID])[locThetaBin][locPBin];
				double locPercentage_Candidates = 100.0*double(locNumClose_Candidates)/double(locNumThrown);
				dHist_CloseEfficiencies_Candidates[locPID]->SetBinContent(locThetaBin + 1, locPBin + 1, locPercentage_Candidates/100.0);

				if(gHistWireBasedFlag)
				{
					locEntries_WireBased += (dNumTimesClose_WireBased[locPID])[locThetaBin][locPBin];
					unsigned int locNumClose_WireBased = (dNumTimesClose_WireBased[locPID])[locThetaBin][locPBin];
					double locPercentage_WireBased = 100.0*double(locNumClose_WireBased)/double(locNumThrown);
					dHist_CloseEfficiencies_WireBased[locPID]->SetBinContent(locThetaBin + 1, locPBin + 1, locPercentage_WireBased/100.0);
				}

				locEntries_TimeBased += (dNumTimesClose_TimeBased[locPID])[locThetaBin][locPBin];
				unsigned int locNumClose_TimeBased = (dNumTimesClose_TimeBased[locPID])[locThetaBin][locPBin];
				double locPercentage_TimeBased = 100.0*double(locNumClose_TimeBased)/double(locNumThrown);
				dHist_CloseEfficiencies_TimeBased[locPID]->SetBinContent(locThetaBin + 1, locPBin + 1, locPercentage_TimeBased/100.0);
			}
		}
		dHist_NumThrown[locPID]->SetEntries(locEntries_Thrown);
		dHist_CloseEfficiencies_Candidates[locPID]->SetEntries(locEntries_Candidates);
		if(gHistWireBasedFlag)
			dHist_CloseEfficiencies_WireBased[locPID]->SetEntries(locEntries_WireBased);
		dHist_CloseEfficiencies_TimeBased[locPID]->SetEntries(locEntries_TimeBased);
	}
	return NOERROR;
}

