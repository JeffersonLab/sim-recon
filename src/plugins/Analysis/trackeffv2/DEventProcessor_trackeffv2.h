#ifndef _DEventProcessor_trackeffv2_
#define _DEventProcessor_trackeffv2_

#include <map>
#include <set>
#include <deque>

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>

#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DMCThrown.h>

#include "ANALYSIS/DMCThrownMatching_factory.h"

using namespace std;
using namespace jana;

class DEventProcessor_trackeffv2 : public JEventProcessor
{
	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop* locEventLoop, int32_t runnumber);
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void MakeAndChangeTo_SubDirectory(string locDirName);
		void Find_GenReconMatches(const vector<const DMCThrown*>& locInputMCThrowns, const vector<const DKinematicData*>& locInputKinematicDataVector, map<const DKinematicData*, const DMCThrown*>& locDataToThrownMap, map<const DMCThrown*, const DKinematicData*>& locThrownToDataMap, bool locRequirePIDMatchFlag);
		bool Check_IsCloseMatch(const DKinematicData* locKinematicData, const DMCThrown* locMCThrown, bool locIsCandidateFlag);

		DMCThrownMatching_factory* dMCThrownMatchingFactory;

		double dMinimumMatchFOM;
		deque<Particle_t> dFinalStatePIDs;

		double dKinematicsHistBinRange_MinP;
		double dKinematicsHistBinRange_MaxP;
		size_t dKinematicsHistBins_NumPBins;

		double dKinematicsHistBinRange_MinTheta;
		double dKinematicsHistBinRange_MaxTheta;
		size_t dKinematicsHistBins_NumThetaBins;

		double dEfficiencyHists_MinP;
		double dEfficiencyHists_MaxP;
		size_t dEfficiencyHists_NumPBins;

		double dEfficiencyHists_MinTheta;
		double dEfficiencyHists_MaxTheta;
		size_t dEfficiencyHists_NumThetaBins;

		TH1D* dHist_NumCandidates;
		TH1D* dHist_NumWireBased;
		TH1D* dHist_NumTimeBased;

		map<Particle_t, TH2D*> dHist_NumThrown;
		map<Particle_t, TH2D*> dHist_CloseEfficiencies_Candidates;
		map<Particle_t, TH2D*> dHist_CloseEfficiencies_WireBased;
		map<Particle_t, TH2D*> dHist_CloseEfficiencies_TimeBased;

		//1st deque index is thrown p-bin, 2nd is thrown theta-bin
		map<Particle_t, deque<deque<unsigned int> > > dNumThrown;
		map<Particle_t, deque<deque<unsigned int> > > dNumTimesClose_Candidates;
		map<Particle_t, deque<deque<unsigned int> > > dNumTimesClose_WireBased;
		map<Particle_t, deque<deque<unsigned int> > > dNumTimesClose_TimeBased;

		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPOverP_Candidates;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaTheta_Candidates;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPhi_Candidates;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaVertexZ_Candidates;

		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPOverP_WireBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaTheta_WireBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPhi_WireBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaVertexZ_WireBased;

		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPOverP_TimeBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaTheta_TimeBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaPhi_TimeBased;
		map<Particle_t, deque<deque<TH1D*> > > dHistMap_DeltaVertexZ_TimeBased;
};

#endif // _DEventProcessor_trackeffv2_

