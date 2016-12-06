// $Id$
//
// File: DEventProcessor_monitoring_hists.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_monitoring_hists.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"
{
	void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new DEventProcessor_monitoring_hists());
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_monitoring_hists::init(void)
{
	dNumMemoryMonitorEvents = 0;
	gPARMS->SetDefaultParameter("MONITOR:MEMORY_EVENTS", dNumMemoryMonitorEvents);

	string locOutputFileName = "hd_root.root";
	if(gPARMS->Exists("OUTPUT_FILENAME"))
		gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

	//go to file
	TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
	if(locFile != NULL)
		locFile->cd("");
	else
		gDirectory->Cd("/");

	//go to directory
	TDirectoryFile* locSubDirectory = static_cast<TDirectoryFile*>(gDirectory->Get("Independent"));
	if(locSubDirectory == NULL) //else folder already created
		locSubDirectory = new TDirectoryFile("Independent", "Independent");
	locSubDirectory->cd();

	dHist_IsEvent = new TH1D("IsEvent", "Is the event an event?", 2, -0.5, 1.5);
	dHist_IsEvent->GetXaxis()->SetBinLabel(1, "False");
	dHist_IsEvent->GetXaxis()->SetBinLabel(2, "True");

	if(dNumMemoryMonitorEvents > 0)
	{
	        dVirtualMemoryVsEventNumber = new TH1D("VirtualMemoryVsEventNumber", ";Event # ;Virtual Memory (MB)", dNumMemoryMonitorEvents, 0.5, (double)dNumMemoryMonitorEvents + 0.5);
        	dResidentMemoryVsEventNumber = new TH1D("ResidentMemoryVsEventNumber", ";Event # ;Resident Memory (MB)", dNumMemoryMonitorEvents, 0.5, (double)dNumMemoryMonitorEvents + 0.5);
	}
	gDirectory->cd("..");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_monitoring_hists::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//Initialize Actions
	dHistogramAction_NumReconstructedObjects.Initialize(locEventLoop);
	dHistogramAction_Reconstruction.Initialize(locEventLoop);
	dHistogramAction_EventVertex.Initialize(locEventLoop);

	dHistogramAction_DetectorMatching.Initialize(locEventLoop);
	dHistogramAction_DetectorMatchParams.Initialize(locEventLoop);
	dHistogramAction_Neutrals.Initialize(locEventLoop);
	dHistogramAction_DetectorPID.Initialize(locEventLoop);

	dHistogramAction_TrackMultiplicity.Initialize(locEventLoop);
	dHistogramAction_DetectedParticleKinematics.Initialize(locEventLoop);
	dHistogramAction_TrackShowerErrors.Initialize(locEventLoop);

	if(dNumMemoryMonitorEvents > 0)
	{
		dHistogramAction_ObjectMemory.dMaxNumEvents = dNumMemoryMonitorEvents;
		dHistogramAction_ObjectMemory.Initialize(locEventLoop);
	}

	if(!locMCThrowns.empty())
	{
		dHistogramAction_ThrownParticleKinematics.Initialize(locEventLoop);
		dHistogramAction_ReconnedThrownKinematics.Initialize(locEventLoop);
		dHistogramAction_GenReconTrackComparison.Initialize(locEventLoop);
	}

	return NOERROR;
}

void DEventProcessor_monitoring_hists::Read_MemoryUsage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage     = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_monitoring_hists::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock

	//CHECK TRIGGER TYPE
	const DTrigger* locTrigger = NULL;
	locEventLoop->GetSingle(locTrigger);
	if(!locTrigger->Get_IsPhysicsEvent())
		return NOERROR;

	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	{
		dHist_IsEvent->Fill(1);
		if(dNumMemoryMonitorEvents > 0)
		{
			double vm, rss;
			Read_MemoryUsage(vm, rss);
			dVirtualMemoryVsEventNumber->SetBinContent(eventnumber, vm / 1024.0);
			dResidentMemoryVsEventNumber->SetBinContent(eventnumber, rss / 1024.0);
		}
	}
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//Fill reaction-independent histograms.
	dHistogramAction_NumReconstructedObjects(locEventLoop);
	dHistogramAction_Reconstruction(locEventLoop);
	dHistogramAction_EventVertex(locEventLoop);

	dHistogramAction_DetectorMatching(locEventLoop);
	dHistogramAction_DetectorMatchParams(locEventLoop);
	dHistogramAction_Neutrals(locEventLoop);
	dHistogramAction_DetectorPID(locEventLoop);

	dHistogramAction_TrackMultiplicity(locEventLoop);
	dHistogramAction_DetectedParticleKinematics(locEventLoop);
	dHistogramAction_TrackShowerErrors(locEventLoop);
	if(dNumMemoryMonitorEvents > 0)
		dHistogramAction_ObjectMemory(locEventLoop);

	if(!locMCThrowns.empty())
	{
		dHistogramAction_ThrownParticleKinematics(locEventLoop);
		dHistogramAction_ReconnedThrownKinematics(locEventLoop);
		dHistogramAction_GenReconTrackComparison(locEventLoop);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_monitoring_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_monitoring_hists::fini(void)
{
	return NOERROR;
}

