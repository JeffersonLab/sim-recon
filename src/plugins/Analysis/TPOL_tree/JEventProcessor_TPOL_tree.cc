// $Id$
//
//    File: JEventProcessor_TPOL_tree.cc
// Created: Thu Feb  4 16:11:54 EST 2016
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.4.4.el7.x86_64 x86_64)
//
#include <iostream>
#include <cmath>
#include <stdint.h>
#include "JEventProcessor_TPOL_tree.h"
using namespace jana;
using namespace std;

#include <TRIGGER/DL1Trigger.h>
#include <TPOL/DTPOLHit_factory.h>
#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <PAIR_SPECTROMETER/DPSPair.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGHGeometry.h>
#include <TAGGER/DTAGMGeometry.h>
#include <PID/DBeamPhoton.h>
#include <DAQ/DBeamCurrent.h>

const int NSECTORS = DTPOLHit_factory::NSECTORS;
const double SECTOR_DIVISION = DTPOLHit_factory::SECTOR_DIVISION;
const int NC_PSC = DPSGeometry::NUM_COARSE_COLUMNS;
const int NC_PS = DPSGeometry::NUM_FINE_COLUMNS;
const int NC_TAGH = DTAGHGeometry::kCounterCount;
const int NC_TAGM = DTAGMGeometry::kColumnCount;
const bool VERBOSE = false;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TVector3.h>

extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_TPOL_tree());
    }
} // "C"

thread_local DTreeFillData JEventProcessor_TPOL_tree::dTreeFillData;

//------------------
// JEventProcessor_TPOL_tree (Constructor)
//------------------
JEventProcessor_TPOL_tree::JEventProcessor_TPOL_tree()
{

}

//------------------
// ~JEventProcessor_TPOL_tree (Destructor)
//------------------
JEventProcessor_TPOL_tree::~JEventProcessor_TPOL_tree()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TPOL_tree::init(void)
{
    // This is called once at program startup. If you are creating
    // and filling historgrams in this plugin, you should lock the
    // ROOT mutex like this:
    //


    // Construct DTreeInterface and register branches for TPOL tree
    double locNumTAGHhits = 500;
    double locNumTAGMhits = 200;

    dTreeInterface = DTreeInterface::Create_DTreeInterface("TPOL_tree","tree_TPOL.root");
    DTreeBranchRegister locTreeBranchRegister;

    locTreeBranchRegister.Register_Single<UShort_t>("nadc");
    locTreeBranchRegister.Register_Single<ULong64_t>("eventnum");
    locTreeBranchRegister.Register_Single<Bool_t>("isFiducial");
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("rocid","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("slot","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("channel","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("itrigger","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<ULong64_t>("w_integral","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<ULong64_t>("w_max","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<ULong64_t>("w_min","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<ULong64_t>("w_samp1","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("sector","nadc",NSECTORS);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("phi","nadc",NSECTORS);
    locTreeBranchRegister.Register_Single<ULong64_t>("ntpol");
    locTreeBranchRegister.Register_FundamentalArray<UShort_t>("waveform","ntpol",NSECTORS*150);
    locTreeBranchRegister.Register_Single<Double_t>("PSenergy_lhit");
    locTreeBranchRegister.Register_Single<Double_t>("PSenergy_rhit");
    locTreeBranchRegister.Register_Single<Double_t>("PSCtime_lhit");
    locTreeBranchRegister.Register_Single<Double_t>("PSCtime_rhit");
    locTreeBranchRegister.Register_Single<Double_t>("PStime_lhit");
    locTreeBranchRegister.Register_Single<Double_t>("PStime_rhit");
    locTreeBranchRegister.Register_Single<UShort_t>("ntagh");
    locTreeBranchRegister.Register_FundamentalArray<Bool_t>("TAGH_DBeam","ntagh",locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHenergy","ntagh",locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGHtime","ntagh",locNumTAGHhits);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("TAGHcounter","ntagh",locNumTAGHhits);
    locTreeBranchRegister.Register_Single<UShort_t>("ntagm");
    locTreeBranchRegister.Register_FundamentalArray<Bool_t>("TAGM_DBeam","ntagm",locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMenergy","ntagm",locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<Double_t>("TAGMtime","ntagm",locNumTAGMhits);
    locTreeBranchRegister.Register_FundamentalArray<UInt_t>("TAGMcolumn","ntagm",locNumTAGMhits);
    dTreeInterface->Create_Branches(locTreeBranchRegister);

    //
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TPOL_tree::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    // Set up beam current factory for is_Fiducial
    dBeamCurrentFactory = new DBeamCurrent_factory();
    dBeamCurrentFactory->init();
    dBeamCurrentFactory->brun(eventLoop, runnumber);

    // Set up beam period for beam bunches
    vector<double> locBeamPeriodVector;
    eventLoop->GetCalib("PHOTON_BEAM/RF/beam_period",locBeamPeriodVector);
    dBeamBunchPeriod = locBeamPeriodVector[0];
    
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TPOL_tree::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // This is called for every event. Use of common resources like writing
    // to a file or filling a histogram should be mutex protected. Using
    // loop->Get(...) to get reconstructed objects (and thereby activating the
    // reconstruction algorithm) should be done outside of any mutex lock
    // since multiple threads may call this method at the same time.
    
    // Construct trigger information
    const DL1Trigger *trig_words = NULL;
    uint32_t trig_mask, fp_trig_mask;
    try {
        loop->GetSingle(trig_words);
    } catch(...) {};
    if (trig_words) {
        trig_mask = trig_words->trig_mask;
        fp_trig_mask = trig_words->fp_trig_mask;
    }
    else {
        trig_mask = 0;
        fp_trig_mask = 0;
    }
    int trig_bits = fp_trig_mask > 0 ? 10 + fp_trig_mask:trig_mask;
    // skim PS triggers
    if (trig_bits!=8) {
        return NOERROR;
    }
    
    // Get fADC 250 windowraws 
    vector<const Df250WindowRawData*> windowraws;
    loop->Get(windowraws);

    // Get coarse PS pairs
    vector<const DPSCPair*> cpairs;
    loop->Get(cpairs);

    // Get fine PS pairs
    vector<const DPSPair*> fpairs;
    loop->Get(fpairs);

    // Get TAGH hits
    vector<const DTAGHHit*> taghhits;
    loop->Get(taghhits);

    // Get TAGM hits
    vector<const DTAGMHit*> tagmhits;
    loop->Get(tagmhits);

    // Get beam photons
    vector<const DBeamPhoton*> beamPhotons;
    loop->Get(beamPhotons);

    // Get beam current
    vector<const DBeamCurrent*> beamCurrent;
    loop->Get(beamCurrent);

    japp->RootFillLock(this);
    if (!beamCurrent.empty())
    {
	// Check that photons are is_Fiducial
    	Bool_t isFiducial = beamCurrent[0]->is_fiducial;
    	dTreeFillData.Fill_Single<Bool_t>("isFiducial",isFiducial);
    }

    // PSC coincidences
    if (cpairs.size()>=1) {
        // take pair with smallest time difference from sorted vector
        const DPSCHit* clhit = cpairs[0]->ee.first; // left hit in coarse PS
        const DPSCHit* crhit = cpairs[0]->ee.second;// right hit in coarse PS
 	double PSC_tdiff = fabs(clhit->t-crhit->t);
	if (PSC_tdiff > 6.0)
	{
		japp->RootFillUnLock(this);
		return NOERROR;
	}
        // PSC,PS coincidences
        if (fpairs.size()>=1) {
            // take pair with smallest time difference from sorted vector
            const DPSPair::PSClust* flhit = fpairs[0]->ee.first;  // left hit in fine PS
            const DPSPair::PSClust* frhit = fpairs[0]->ee.second; // right hit in fine PS
            if(flhit->column < geomModuleColumn[clhit->module-1][0] || flhit->column > geomModuleColumn[clhit->module-1][1])
	    {
    	        japp->RootFillUnLock(this);
	        return NOERROR;
	    }
	    if(frhit->column < geomModuleColumn[crhit->module-1][0] || frhit->column > geomModuleColumn[crhit->module-1][1])
	    {	
		japp->RootFillUnLock(this);
		return NOERROR;
	    }
            double t_lhit = clhit->t; 
            double E_pair = flhit->E+frhit->E;
     
	    dTreeFillData.Fill_Single<ULong64_t>("eventnum",eventnumber);
	    dTreeFillData.Fill_Single<Double_t>("PSenergy_lhit",flhit->E);
            dTreeFillData.Fill_Single<Double_t>("PSenergy_rhit",frhit->E);
	    dTreeFillData.Fill_Single<Double_t>("PSCtime_lhit",clhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PSCtime_rhit",crhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PStime_lhit",flhit->t);
	    dTreeFillData.Fill_Single<Double_t>("PStime_rhit",frhit->t);

	    // PSC, PS, and TAGX coincidences
	    // Loop over TAGH hits and match to DBeamPhotons
	    unsigned int htag = 0;
	    unsigned int htag_DBeam = 0;
            double EdiffMax = 0.5; double tdiffMax = 10.0*dBeamBunchPeriod;
            for (unsigned int i=0; i < taghhits.size(); i++) {
                const DTAGHHit* tag = taghhits[i];
                if (!tag->has_TDC||!tag->has_fADC) continue;
		if (std::isnan(tag->t) || std::isnan(tag->E)) continue;
                if (fabs(E_pair-tag->E) > EdiffMax || fabs(t_lhit-tag->t) > tdiffMax) continue;
		dTreeFillData.Fill_Array<Double_t>("TAGHenergy",tag->E,htag);
                dTreeFillData.Fill_Array<Double_t>("TAGHtime",tag->t,htag);
                dTreeFillData.Fill_Array<UInt_t>("TAGHcounter",tag->counter_id,htag);
		unsigned int same = 0;
                for (unsigned int j=0; j < beamPhotons.size(); j++)
		{
			const DTAGHHit* tagh;
			beamPhotons[j]->GetSingle(tagh);
			if (!tagh) continue;
			if (!tagh->has_TDC || !tagh->has_fADC) continue;
			if (std::isnan(tagh->t) || std::isnan(tagh->E)) 
			{
				if (VERBOSE) jerr<<"Found TAGH with NAN."<<tagh->counter_id<<endl;
				continue;
			}
			if (fabs(E_pair-tagh->E) > EdiffMax || fabs(t_lhit-tagh->t) > tdiffMax) continue;
			if (tagh->t != tag->t || tagh->E != tag->E || tagh->counter_id != tag->counter_id) continue;
			same++;
			htag_DBeam++;
		}
		if (same > 1 && VERBOSE) jerr<<"Found more than one match for TAGH."<<endl;
		if (same == 0) dTreeFillData.Fill_Array<Bool_t>("TAGH_DBeam",false,htag);
		else dTreeFillData.Fill_Array<Bool_t>("TAGH_DBeam",true,htag);
		htag++;
            }

	    // Loop over TAGM hits and match to DBeamPhotons
	    unsigned int mtag = 0;
	    unsigned int mtag_DBeam = 0;
            for (unsigned int i=0; i < tagmhits.size(); i++) {
                const DTAGMHit* tag = tagmhits[i];
                if (!tag->has_TDC||!tag->has_fADC) continue;
                if (tag->row!=0) continue;
		if (std::isnan(tag->t) || std::isnan(tag->E)) continue;
                if (fabs(E_pair-tag->E) > EdiffMax || fabs(t_lhit-tag->t) > tdiffMax) continue;
		dTreeFillData.Fill_Array<Double_t>("TAGMenergy",tag->E,mtag);
                dTreeFillData.Fill_Array<Double_t>("TAGMtime",tag->t,mtag);
                dTreeFillData.Fill_Array<UInt_t>("TAGMcolumn",tag->column,mtag);
		unsigned int same = 0;
                for (unsigned int j=0; j < beamPhotons.size(); j++)
                {
                        const DTAGMHit* tagm;
                        beamPhotons[j]->GetSingle(tagm);
                        if (!tagm) continue;
                        if (!tagm->has_TDC || !tagm->has_fADC) continue;
			if (tagm->row != 0) continue;
			if (std::isnan(tagm->t) || std::isnan(tagm->E)) 
			{
				if (VERBOSE) jerr<<"Found TAGM with NAN."<<tagm->column<<endl;
				continue;
			}
                        if (fabs(E_pair-tagm->E) > EdiffMax || fabs(t_lhit-tagm->t) > tdiffMax) continue;
			if (tagm->t != tag->t || tagm->E != tag->E || tagm->column != tag->column) continue;
                        same++;
                        mtag_DBeam++;
                }
		if (same > 1 && VERBOSE) jerr<<"Found more than one match for TAGM."<<endl;
                if (same == 0) dTreeFillData.Fill_Array<Bool_t>("TAGM_DBeam",false,mtag);
                else dTreeFillData.Fill_Array<Bool_t>("TAGM_DBeam",true,mtag);
		mtag++;
	    }

            // Ensure TAGH hits match DBeamPhotons
            unsigned int htag_Check = 0;
            for (unsigned int i=0; i < beamPhotons.size(); i++) {
                const DTAGHHit* tag;
                beamPhotons[i]->GetSingle(tag);
		if (!tag) continue;
		if (!tag->has_TDC||!tag->has_fADC) continue;
		if (std::isnan(tag->t) || std::isnan(tag->E)) continue;
		if (fabs(E_pair-tag->E) > EdiffMax || fabs(t_lhit-tag->t) > tdiffMax) continue;	
                htag_Check++;
            }
	    if (htag_Check != htag_DBeam && VERBOSE) jerr<<"TAGH: "<<htag_DBeam<<" , "<<htag_Check<<endl;

	    // Ensure TAGM hits match DBeamPhotons
	    unsigned int mtag_Check = 0;
            for (unsigned int i=0; i < beamPhotons.size(); i++) {
                const DTAGMHit* tag;
		beamPhotons[i]->GetSingle(tag);
		if(!tag) continue;
                if (!tag->has_TDC||!tag->has_fADC) continue;
                if (tag->row!=0) continue;
		if (std::isnan(tag->t) || std::isnan(tag->E)) continue;
                if (fabs(E_pair-tag->E) > EdiffMax || fabs(t_lhit-tag->t) > tdiffMax) continue;
                mtag_Check++;
            }
	    if (mtag_Check != mtag_DBeam && VERBOSE) jerr<<"TAGM: "<<mtag_DBeam<<" , "<<mtag_Check<<endl;

	    dTreeFillData.Fill_Single<UShort_t>("ntagh",htag);
	    dTreeFillData.Fill_Single<UShort_t>("ntagm",mtag);

	    // Loop over windowraws to collect TPOL hits
	    // No cuts are applied to the TPOL hits
            unsigned int hit = 0;
	    unsigned int ntpol = 0;
            for(unsigned int i=0; i< windowraws.size(); i++) {
                const Df250WindowRawData *windowraw = windowraws[i];
                if (windowraw->rocid!=84) continue;
                if (!(windowraw->slot==13||windowraw->slot==14)) continue; // azimuthal sectors, rings: 15,16
                unsigned int rocid = windowraw->rocid;
                unsigned int slot = windowraw->slot;
                unsigned int channel = windowraw->channel;
                unsigned int itrigger = windowraw->itrigger;
                // Get a vector of the samples for this channel
                const vector<uint16_t> &samplesvector = windowraw->samples;
                unsigned int nsamples = samplesvector.size();

                // loop over the samples to calculate integral, min, max
                if (nsamples==0 && VERBOSE) jerr << "Raw samples vector is empty." << endl;

                unsigned int w_integral = 0;
		unsigned int w_max = 0;
		unsigned int w_min = 0;
		unsigned int w_samp1 = 0;
                for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
		    dTreeFillData.Fill_Array<UShort_t>("waveform",samplesvector[c_samp],ntpol);
		    ntpol++;
                    if (c_samp==0) {  // use first sample for initialization
                        w_integral = samplesvector[0];
                        w_min = samplesvector[0];
                        w_max = samplesvector[0];
                        w_samp1 = samplesvector[0];
                    } else {
                        w_integral += samplesvector[c_samp];
                        if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
                        if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
                    }
                }

                unsigned int sector = GetSector(slot,channel);
                double phi = GetPhi(sector);
	
		dTreeFillData.Fill_Array<UInt_t>("rocid",rocid,hit);
		dTreeFillData.Fill_Array<UInt_t>("slot",slot,hit);
		dTreeFillData.Fill_Array<UInt_t>("channel",channel,hit);
		dTreeFillData.Fill_Array<UInt_t>("itrigger",itrigger,hit);
		dTreeFillData.Fill_Array<ULong64_t>("w_integral",w_integral,hit);
		dTreeFillData.Fill_Array<ULong64_t>("w_max",w_max,hit);
		dTreeFillData.Fill_Array<ULong64_t>("w_min",w_min,hit);
		dTreeFillData.Fill_Array<ULong64_t>("w_samp1",w_samp1,hit);
		dTreeFillData.Fill_Array<UInt_t>("sector",sector,hit);
		dTreeFillData.Fill_Array<ULong64_t>("phi",phi,hit);
		hit++;
            }
            unsigned int nadc = hit;
            if (nadc>NSECTORS && VERBOSE) jerr << "TPOL_tree plugin error: nadc exceeds nmax(" << NSECTORS << ")." << endl;
            dTreeFillData.Fill_Single<UShort_t>("nadc",nadc);
	    dTreeFillData.Fill_Single<ULong64_t>("ntpol",ntpol);
	    dTreeInterface->Fill(dTreeFillData);
        }
    }
    japp->RootFillUnLock(this);
    //
    return NOERROR;
}

int JEventProcessor_TPOL_tree::GetSector(int slot,int channel)
{
    int sector = 0;
    if (slot == 13) sector = 25 - channel;
    if (slot == 14) {
        if (channel <= 8) sector = 9 - channel;
        else sector = NSECTORS + 9 - channel;
    }
    // fix cable swap
    if (sector == 9) sector = 6;
    else if (sector == 6) sector = 9;
    if (sector == 0 && VERBOSE) jerr << "sector did not change from initial value (0)." << endl;
    return sector;
}
double JEventProcessor_TPOL_tree::GetPhi(int sector)
{
    double phi = -10.0;
    if(sector <= 8) phi = (sector + 23)*SECTOR_DIVISION + 0.5*SECTOR_DIVISION;
    if(sector >= 9) phi = (sector - 9)*SECTOR_DIVISION + 0.5*SECTOR_DIVISION;
    return phi;
}
double JEventProcessor_TPOL_tree::GetPulseTime(const vector<uint16_t> waveform,double w_min,double w_max,double minpeakheight)
{
    // find the time to cross half peak height
    int lastbelowsamp=0; double peakheight = w_max-w_min;
    double threshold = w_min + peakheight/2.0;
    double  firstaboveheight=0, lastbelowheight=0;
    double w_time=0;
    if (peakheight > minpeakheight) {
        for (uint16_t c_samp=0; c_samp<waveform.size(); c_samp++) {
            if (waveform[c_samp]>threshold) {
                firstaboveheight = waveform[c_samp];
                lastbelowsamp = c_samp-1;
                lastbelowheight = waveform[c_samp-1];
                break;
            }
        }
        w_time =  lastbelowsamp + (threshold-lastbelowheight)/(firstaboveheight-lastbelowheight);
    }
    return 64.0*w_time;
}
//------------------
// erun
//------------------
jerror_t JEventProcessor_TPOL_tree::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TPOL_tree::fini(void)
{
    delete dTreeInterface;
    // Called before program exit after event processing is finished.
    return NOERROR;
}
