// $Id$
//
//    File: JEventProcessor_PS_flux.cc
//
//
#include <iostream>
#include <sstream>
#include "JEventProcessor_PS_flux.h"

using namespace std;
using namespace jana;

#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <PAIR_SPECTROMETER/DPSPair.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGHGeometry.h>
#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGMGeometry.h>
#include <PID/DBeamPhoton.h>
#include <DAQ/DBeamCurrent.h>

#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>

const int NC_PSC = DPSGeometry::NUM_COARSE_COLUMNS; // 8: number of PSC modules (counters) per arm
const int NC_PS = DPSGeometry::NUM_FINE_COLUMNS; // 145: number of PS columns (tiles) per arm
const int NC_TAGH = DTAGHGeometry::kCounterCount; // 274: number of TAGH counters
const int NC_TAGM = DTAGMGeometry::kColumnCount; // 102: number of TAGM columns
// PS,PSC coincidences
static TH1I *psflux_num_events;
static TH1F *hFiducialTime;
static TH1I *hBeamCurrentTime;
static TH1I *hBeamCurrentTimeFiducial;
static TH2I *hPS_tdiffVsE;
static TH2I *hPSC_tdiffVsE;
static TH2I *hPSPSC_tdiffVsE;
static TH1I *hPS_E;
// PSC,PS,TAGH coincidences
static TH2I *hPSTAGH_tdiffVsEdiff;
static TH2I *hPSTAGH_tdiffVsEtagh;
static TH2I *hPSTAGH_tdiffVsCounter;
static TH2I *hPSTAGH_EdiffVsTAGHCounterID;
static TH2I *hPSTAGH_EdiffVsEtagh;
// PSC,PS,TAGM coincidences
static TH2I *hPSTAGM_tdiffVsEdiff;
static TH2I *hPSTAGM_tdiffVsEtagm;
static TH2I *hPSTAGM_tdiffVsColumn;
static TH2I *hPSTAGM_EdiffVsTAGMColumn;
static TH2I *hPSTAGM_EdiffVsEtagm;
// PSC,PS,TAG coincidences
static TH2I *hPSTAG_tdiffVsEtag;

//-------------------------
// Routine used to create our JEventProcessor
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_PS_flux());

    }
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_PS_flux::dTreeFillData;

//------------------
// JEventProcessor_PS_flux (Constructor)
//------------------
JEventProcessor_PS_flux::JEventProcessor_PS_flux()
{

}

//------------------
// ~JEventProcessor_PS_flux (Destructor)
//------------------
JEventProcessor_PS_flux::~JEventProcessor_PS_flux()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_PS_flux::init(void)
{
    // energy binning
    const double Ebw_PS = 0.05;
    const double Ebl_PS = 5.0; 
    const double Ebh_PS = 13.0;
    const int NEb_PS = int((Ebh_PS-Ebl_PS)/Ebw_PS);
    const double Ebw_TAG = 0.1;
    const double Ebl_TAG = 2.5;
    const double Ebh_TAG = 10.5;
    const int NEb_TAG = int((Ebh_TAG-Ebl_TAG)/Ebw_TAG);
    // time binning
    const double Tbl_TAG = -100.0;
    const double Tbh_TAG = 100.0;
    const int NTb_TAG = 200;
    // energy difference binning
    const int NEdb = 200; const double Edbl = -1.0; const double Edbh = 1.0;

    t_start = -999.;
    t_end = 0;

    // create root folder for pspair and cd to it, store main dir
    TDirectory *psFluxDir = gDirectory->mkdir("PS_flux");
    psFluxDir->cd();
    // book hists
    psflux_num_events = new TH1I("psflux_num_events","PS flux number of events",2,0.5,2.5);
    hFiducialTime = new TH1F("fiducialTime", "Fiducial time", 1, 0, 1);
    hBeamCurrentTime = new TH1I("beamCurrentTime", "; Event time from DBeamCurrent for all events (seconds)", 1000, 0, 500);
    hBeamCurrentTimeFiducial = new TH1I("beamCurrentTimeFiducial", "; Event time from DBeamCurrent for fiducial events (seconds)", 1000, 0, 500);

    //
    gDirectory->mkdir("PSC_PS")->cd();
    hPS_E = new TH1I("PS_E","PS pair energy; PS pair energy [GeV]",NEb_PS,Ebl_PS,Ebh_PS);
    hPS_tdiffVsE = new TH2I("PS_tdiffVsE","PS pair time difference vs. PS pair energy;PS pair energy [GeV];PS pair time difference [ns]",NEb_PS,Ebl_PS,Ebh_PS,200,-10.0,10.0);
    hPSC_tdiffVsE = new TH2I("PSC_tdiffVsE","PSC pair time difference vs. PS pair energy;PS pair energy [GeV];PSC pair time difference [ns]",NEb_PS,Ebl_PS,Ebh_PS,200,-10.0,10.0);
    hPSPSC_tdiffVsE = new TH2I("PSPSC_tdiffVsE","PSC/PS pair time difference of differences vs. PS pair energy;PS pair energy [GeV];PSC/PS pair time difference of differences [ns]",NEb_PS,Ebl_PS,Ebh_PS,200,-10.0,10.0);
    psFluxDir->cd();
    //
    gDirectory->mkdir("PSC_PS_TAGH")->cd();
    hPSTAGH_tdiffVsEdiff = new TH2I("PSTAGH_tdiffVsEdiff","PS pair - TAGH: PS-TAGH time difference vs. PS-TAGH energy difference;E(PS) - E(TAGH) [GeV];PSC/TAGH time difference [ns]",NEdb,Edbl,Edbh,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGH_tdiffVsEtagh = new TH2I("PSTAGH_tdiffVsEtagh","PSC/TAGH time difference vs. TAGH energy; TAGH energy [GeV];PSC/TAGH time difference [ns]",NEb_PS,Ebl_PS,Ebh_PS,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGH_tdiffVsCounter = new TH2I("PSTAGH_tdiffVsCounter","PSC/TAGH time difference vs. TAGH Counter; TAGH Counter;PSC/TAGH time difference [ns]",NC_TAGH,0.5,0.5+NC_TAGH,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGH_EdiffVsTAGHCounterID = new TH2I("PSTAGH_EdiffVsTAGHCounterID","PS pair - TAGH: PS-TAGH energy difference vs. TAGH counter ID;TAGH counter ID;E(PS) - E(TAGH) [GeV]",NC_TAGH,0.5,0.5+NC_TAGH,NEdb,Edbl,Edbh);
    hPSTAGH_EdiffVsEtagh = new TH2I("PSTAGH_EdiffVsEtagh","PS pair - TAGH: PS-TAGH energy difference vs. TAGH energy;TAGH energy [GeV];E(PS) - E(TAGH) [GeV]",NEb_TAG,Ebl_TAG,Ebh_TAG,NEdb,Edbl,Edbh);
    psFluxDir->cd();
    //
    gDirectory->mkdir("PSC_PS_TAGM")->cd();
    hPSTAGM_tdiffVsEdiff = new TH2I("PSTAGM_tdiffVsEdiff","PS pair - TAGM: PS-TAGM time difference vs. PS-TAGM energy difference;E(PS) - E(TAGM) [GeV];PSC/TAGM time difference [ns]",NEdb,Edbl,Edbh,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGM_tdiffVsEtagm = new TH2I("PSTAGM_tdiffVsEtagm","PSC/TAGM time difference vs. TAGM energy; TAGM energy [GeV];PSC/TAGM time difference [ns]",NEb_PS,Ebl_PS,Ebh_PS,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGM_tdiffVsColumn = new TH2I("PSTAGM_tdiffVsColumn","PSC/TAGM time difference vs. TAGM Column; TAGM Column;PSC/TAGM time difference [ns]",NC_TAGM,0.5,0.5+NC_TAGM,NTb_TAG,Tbl_TAG,Tbh_TAG);
    hPSTAGM_EdiffVsTAGMColumn = new TH2I("PSTAGM_EdiffVsTAGMColumn","PS pair - TAGM: PS-TAGM energy difference vs. TAGM column;TAGM column;E(PS) - E(TAGM) [GeV]",NC_TAGM,0.5,0.5+NC_TAGM,NEdb,Edbl,Edbh);
    hPSTAGM_EdiffVsEtagm = new TH2I("PSTAGM_EdiffVsEtagm","PS pair - TAGM: PS-TAGM energy difference vs. TAGM energy;TAGM energy [GeV];E(PS) - E(TAGM) [GeV]",NEb_PS,Ebl_PS,Ebh_PS,NEdb,Edbl,Edbh);
    psFluxDir->cd();
    gDirectory->mkdir("PSC_PS_TAG")->cd();
    hPSTAG_tdiffVsEtag = new TH2I("PSTAG_tdiffVsEtag","PSC/TAG time difference vs. TAG energy; TAG energy [GeV];PSC/TAG time difference [ns]",NEb_PS,Ebl_PS,Ebh_PS,NTb_TAG,Tbl_TAG,Tbh_TAG);
    psFluxDir->cd();
    gDirectory->cd("../");

    double locNumTAGHhits = 500;
    double locNumTAGMhits = 200;

    //TTREE INTERFACE
    //MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
    dTreeInterface = DTreeInterface::Create_DTreeInterface("PSFlux_Tree", "tree_PSFlux.root");

    //TTREE BRANCHES
    DTreeBranchRegister locTreeBranchRegister;
    
    locTreeBranchRegister.Register_Single<ULong64_t>("EventNumber");
    locTreeBranchRegister.Register_Single<Bool_t>("IsFiducial");
    locTreeBranchRegister.Register_Single<Double_t>("PSCtimeL");
    locTreeBranchRegister.Register_Single<Double_t>("PSCtimeR");
    locTreeBranchRegister.Register_Single<Double_t>("PStimeL");
    locTreeBranchRegister.Register_Single<Double_t>("PStimeR");
    locTreeBranchRegister.Register_Single<Double_t>("PSenergyL");
    locTreeBranchRegister.Register_Single<Double_t>("PSenergyR");
    locTreeBranchRegister.Register_Single<Double_t>("PSPairEnergy");

    locTreeBranchRegister.Register_Single<Int_t>("NumTAGMhits");
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMenergy", "NumTAGMhits", locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMtime", "NumTAGMhits", locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMintegral", "NumTAGMhits", locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMpulse_peak", "NumTAGMhits", locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<bool>("TAGMhasTDC", "NumTAGMhits", locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Int_t>("TAGMcolumn", "NumTAGMhits", locNumTAGMhits);

    locTreeBranchRegister.Register_Single<Int_t>("NumTAGHhits");
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHenergy", "NumTAGHhits", locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHtime", "NumTAGHhits", locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHintegral", "NumTAGHhits", locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHpulse_peak", "NumTAGHhits", locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<bool>("TAGHhasTDC", "NumTAGHhits", locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Int_t>("TAGHcounter", "NumTAGHhits", locNumTAGHhits);

    //REGISTER BRANCHES
    dTreeInterface->Create_Branches(locTreeBranchRegister);

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_PS_flux::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    // extract the PS geometry
    vector<const DPSGeometry*> psGeomVect;
    eventLoop->Get(psGeomVect);
    if (psGeomVect.size() < 1)
    return OBJECT_NOT_AVAILABLE;
    const DPSGeometry& psGeom = *(psGeomVect[0]);

    dBeamCurrentFactory = new DBeamCurrent_factory();
    dBeamCurrentFactory->init();
    dBeamCurrentFactory->brun(eventLoop, runnumber);
   
    // PS pair energy binning
    double wl_min=0.05,wr_min = 0.05;
    double Ebw_PS = wl_min + wr_min;
    const double Ebl_PS = psGeom.getElow(0,1) + psGeom.getElow(1,1);
    const double Ebh_PS = psGeom.getEhigh(0,NC_PS) + psGeom.getEhigh(1,NC_PS);
    double range = fabs(Ebh_PS-Ebl_PS);
    int NEb_PS = range/Ebw_PS-int(range/Ebw_PS) < 0.5 ? int(range/Ebw_PS) : int(range/Ebw_PS) + 1;
    double Elows_PS[NEb_PS+1];
    for (int i=0;i<NEb_PS+1;i++) {
        Elows_PS[i] = Ebl_PS + i*Ebw_PS;
    }
    // extract the TAGH geometry
    vector<const DTAGHGeometry*> taghGeomVect;
    eventLoop->Get(taghGeomVect);
    if (taghGeomVect.size() < 1)
    return OBJECT_NOT_AVAILABLE;
    const DTAGHGeometry& taghGeom = *(taghGeomVect[0]);
    // get photon energy bin low of each counter for energy histogram binning
    double Elows_TAGH[NC_TAGH+1];
    for (int i=0;i<NC_TAGH;i++) {
        Elows_TAGH[i] = taghGeom.getElow(NC_TAGH-i);
    }
    // add the upper limit
    Elows_TAGH[NC_TAGH] = taghGeom.getEhigh(1);
    // extract the TAGM geometry
    vector<const DTAGMGeometry*> tagmGeomVect;
    eventLoop->Get(tagmGeomVect);
    if (tagmGeomVect.size() < 1)
    return OBJECT_NOT_AVAILABLE;
    const DTAGMGeometry& tagmGeom = *(tagmGeomVect[0]);
    // get photon energy bin low of each counter for energy histogram binning
    double Elows_TAGM[NC_TAGM+1];
    for (int i=0;i<NC_TAGM;i++) {
        Elows_TAGM[i] = tagmGeom.getElow(NC_TAGM-i);
    }
    // add the upper limit
    Elows_TAGM[NC_TAGM] = tagmGeom.getEhigh(1);
    //
    const int NEdiff = 200;
    double EdiffLows[NEdiff+1];
    for (int i=0;i<NEdiff+1;i++) {
        EdiffLows[i] = -1.0 + i*0.01;
    }
    const int NTb_PS = 200;
    double Tlows_PS[NTb_PS+1];
    for (int i=0;i<=NTb_PS;i++) {
        Tlows_PS[i] = -10.0+i*0.1;
    }

    vector<double> locBeamPeriodVector;
    eventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
    dBeamBunchPeriod = locBeamPeriodVector[0];

    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

    // set variable-width energy bins if histogram is empty
    // PS
    if (hPS_tdiffVsE->GetEntries()==0) hPS_tdiffVsE->SetBins(NEb_PS,Elows_PS,NTb_PS,Tlows_PS);
    if (hPSC_tdiffVsE->GetEntries()==0) hPSC_tdiffVsE->SetBins(NEb_PS,Elows_PS,NTb_PS,Tlows_PS);
    if (hPSPSC_tdiffVsE->GetEntries()==0) hPSPSC_tdiffVsE->SetBins(NEb_PS,Elows_PS,NTb_PS,Tlows_PS);
    // TAGH
    if (hPSTAGH_EdiffVsEtagh->GetEntries()==0) hPSTAGH_EdiffVsEtagh->SetBins(NC_TAGH,Elows_TAGH,NEdiff,EdiffLows);
    // TAGM
    if (hPSTAGM_EdiffVsEtagm->GetEntries()==0) hPSTAGM_EdiffVsEtagm->SetBins(NC_TAGM,Elows_TAGM,NEdiff,EdiffLows);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    //
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_PS_flux::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // This is called for every event. Use of common resources like writing
    // to a file or filling a histogram should be mutex protected. Using
    // loop->Get(...) to get reconstructed objects (and thereby activating the
    // reconstruction algorithm) should be done outside of any mutex lock
    // since multiple threads may call this method at the same time.
    // coarse PS pairs
    vector<const DPSCPair*> cpairs;
    loop->Get(cpairs);
    // fine PS pairs
    vector<const DPSPair*> fpairs;
    loop->Get(fpairs);

    // tagger hits
    vector<const DBeamPhoton*> beamPhotons;
    loop->Get(beamPhotons);

    // beam current and fiducial definition
    vector<const DBeamCurrent*> beamCurrent;
    loop->Get(beamCurrent);
    if( beamCurrent.empty() ) 
	return NOERROR;

    // get start and end time of events relative to start of run from DBeamCurrent
    if(t_start < 0.) {
	t_start = beamCurrent[0]->t;
    }
    t_end = beamCurrent[0]->t;

    // FILL HISTOGRAMS
    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    psflux_num_events->Fill(1);
    hBeamCurrentTime->Fill(beamCurrent[0]->t);
    if(beamCurrent[0]->is_fiducial) {
       int timeBin = hBeamCurrentTimeFiducial->FindBin(beamCurrent[0]->t);
       hBeamCurrentTimeFiducial->SetBinContent(timeBin, 1);
    }

    // PSC coincidences
    if (cpairs.size()>=1) {
        // take pair with smallest time difference from sorted vector
        const DPSCHit* clhit = cpairs[0]->ee.first; // left hit in coarse PS
        const DPSCHit* crhit = cpairs[0]->ee.second;// right hit in coarse PS
        double PSC_tdiff = clhit->t-crhit->t;
	if (fabs(PSC_tdiff) > 6.) {japp->RootFillUnLock(this); return NOERROR;}

        // PSC,PS coincidences
        if (fpairs.size()>=1) {
            psflux_num_events->Fill(2);
            // take pair with smallest time difference from sorted vector
            const DPSHit* flhit = fpairs[0]->ee.first;  // left hit in fine PS
            const DPSHit* frhit = fpairs[0]->ee.second; // right hit in fine PS
            double E_pair = flhit->E+frhit->E;
            double PS_tdiff = flhit->t-frhit->t;
            
            hPS_tdiffVsE->Fill(E_pair,PS_tdiff);
            hPSC_tdiffVsE->Fill(E_pair,PSC_tdiff);
            hPSPSC_tdiffVsE->Fill(E_pair,PSC_tdiff-PS_tdiff);
	    if(fabs(PS_tdiff) < 0.5*dBeamBunchPeriod && fabs(PSC_tdiff) < 0.5*dBeamBunchPeriod)
		    hPS_E->Fill(E_pair);

	    // Fill PS variables
	    dTreeFillData.Fill_Single<ULong64_t>("EventNumber", eventnumber);
	    dTreeFillData.Fill_Single<Bool_t>("IsFiducial", beamCurrent[0]->is_fiducial);
	    dTreeFillData.Fill_Single<Double_t>("PSCtimeL", clhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PSCtimeR", crhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PStimeL", flhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PStimeR", frhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PSenergyL", flhit->E);
	    dTreeFillData.Fill_Single<Double_t>("PSenergyR", frhit->E);
	    dTreeFillData.Fill_Single<Double_t>("PSPairEnergy", E_pair);

	    int itagm = 0;
	    int itagh = 0;
	    int locTAGMhits = 0;
	    int locTAGHhits = 0;

	    // PSC,PS,TAGH coincidences
	    for(unsigned int i=0; i < beamPhotons.size(); i++) {
	        const DTAGHHit* tag;
		beamPhotons[i]->GetSingle(tag);
		if(!tag) continue;
		
		// loose cuts on matching
		double Ediff = E_pair-tag->E;
		double tdiff_left = clhit->t-tag->t;
		hPSTAGH_tdiffVsEdiff->Fill(Ediff,tdiff_left);

		if(fabs(tdiff_left) > 10.*dBeamBunchPeriod || fabs(Ediff) > 4.0) 
			continue;

		dTreeFillData.Fill_Array<Double_t>("TAGHenergy", tag->E, itagh);
		dTreeFillData.Fill_Array<Double_t>("TAGHtime", tag->t, itagh);
		dTreeFillData.Fill_Array<Double_t>("TAGHintegral", tag->integral, itagh);
		dTreeFillData.Fill_Array<Double_t>("TAGHpulse_peak", tag->pulse_peak, itagh);
		dTreeFillData.Fill_Array<bool>("TAGHhasTDC", tag->has_TDC, itagh);
		dTreeFillData.Fill_Array<Int_t>("TAGHcounter", tag->counter_id, itagh);
		itagh++;
		locTAGHhits++;

                hPSTAGH_EdiffVsEtagh->Fill(tag->E,Ediff);
                hPSTAGH_EdiffVsTAGHCounterID->Fill(tag->counter_id,Ediff);
		if (Ediff < -0.45 || Ediff > 0.15) // cut on energy difference
			continue; 
		hPSTAGH_tdiffVsEtagh->Fill(tag->E,tdiff_left);
		hPSTAGH_tdiffVsCounter->Fill(tag->counter_id, tdiff_left);
		hPSTAG_tdiffVsEtag->Fill(tag->E,tdiff_left);
            }
	    // PSC,PS,TAGM coincidences
	    for(unsigned int i=0; i < beamPhotons.size(); i++) {
	        const DTAGMHit* tag;
		beamPhotons[i]->GetSingle(tag);
		if(!tag) continue;

		// loose cuts on matching
		double Ediff = E_pair-tag->E;
		double tdiff_left = clhit->t-tag->t;
		hPSTAGM_tdiffVsEdiff->Fill(Ediff,tdiff_left);

		if(fabs(tdiff_left) > 10.*dBeamBunchPeriod || fabs(Ediff) > 4.0) 
			continue;

		dTreeFillData.Fill_Array<Double_t>("TAGMenergy", tag->E, itagm);
		dTreeFillData.Fill_Array<Double_t>("TAGMtime", tag->t, itagm);
		dTreeFillData.Fill_Array<Double_t>("TAGMintegral", tag->integral, itagm);
		dTreeFillData.Fill_Array<Double_t>("TAGMpulse_peak", tag->pulse_peak, itagm);
		dTreeFillData.Fill_Array<bool>("TAGMhasTDC", tag->has_TDC, itagm);
		dTreeFillData.Fill_Array<Int_t>("TAGMcolumn", tag->column, itagm);
		itagm++;
		locTAGMhits++;

                hPSTAGM_EdiffVsEtagm->Fill(tag->E,Ediff);
                hPSTAGM_EdiffVsTAGMColumn->Fill(tag->column,Ediff);
		if (Ediff < -0.25 || Ediff > 0.15) // loose cut on energy difference
			continue; 
		hPSTAGM_tdiffVsEtagm->Fill(tag->E,tdiff_left);
		hPSTAGM_tdiffVsColumn->Fill(tag->column, tdiff_left);
		hPSTAG_tdiffVsEtag->Fill(tag->E,tdiff_left);
            }

	    //FILL TTREE
	    dTreeFillData.Fill_Single<Int_t>("NumTAGMhits", locTAGMhits);
	    dTreeFillData.Fill_Single<Int_t>("NumTAGHhits", locTAGHhits);

	    dTreeInterface->Fill(dTreeFillData);
        }
    }
    //
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_PS_flux::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.

    //cout<<"t_start = "<<t_start<<" t_end = "<<t_end<<endl;
    double finalFiducialTime = dBeamCurrentFactory->IntegratedFiducialTime(t_start, t_end);
    //cout<<endl<<"Fiducial time: "<<finalFiducialTime<<endl;
    hFiducialTime->SetBinContent(1., finalFiducialTime);

    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_PS_flux::fini(void)
{
    // Called before program exit after event processing is finished.
	
    delete dTreeInterface; //saves trees to file, closes file

    return NOERROR;
}

