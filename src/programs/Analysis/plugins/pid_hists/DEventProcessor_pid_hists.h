#ifndef _DEventProcessor_pid_hists_
#define _DEventProcessor_pid_hists_

#include <vector>
using namespace std;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DMCThrown.h>

class DEventProcessor_pid_hists:public JEventProcessor{
	public:
		DEventProcessor_pid_hists();
		~DEventProcessor_pid_hists();
		const char* className(void){return "DEventProcessor_pid_hists";}
		
		int NgoodPid;
		int NbadRecon;
		int NpipTruth_protonPid;
		int NpipTruth_KpPid;
		int NprotonTruth_pipPid;
		int NprotonTruth_KpPid;
		int NKpTruth_pipPid;
		int NKpTruth_protonPid;
		int NpimTruth_KmPid;
		int NKmTruth_pimPid;
		TH1F* pid;
		TH2F* kine;
		TH2F* kineBadRec;
		TH1F* tof;
		TH1F* deltaBeta;
		TH1F* dEdx;
		TH1F* tdiff_Kpi;
		TH2F* tdiff_Kpi_mom;
		TH1F* tdiff_ppi;
		TH2F* tdiff_ppi_mom;
		TH1F* tdiff_pK;
		TH2F* tdiff_pK_mom;

		pthread_mutex_t mutex;
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double GetTruthMatchingFOM(int trackIndex,const DTrackTimeBased *track,vector<const DMCThrown*>mcthrowns);
		void GetThrownIndex(const DKinematicData *kd, int &MAX_TRACKS, double &f, int &track);
		void GetThrownIndex2(const DKinematicData *kd,vector<const DMCThrown*>mcthrowns,int &track);
		double GetMomChisq(const DKinematicData *kd,const DKinematicData *kd_thrown);
		void printInfo(const DKinematicData *kd,int event,int track,int N,TString string);
		int GetType(const DTrackTimeBased *track);
			
}; 

#endif // _DEventProcessor_pid_hists_

