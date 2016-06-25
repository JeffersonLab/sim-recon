// $Id$
//
//    File: JEventProcessor_TPOL_tree.h
// Created: Thu Feb  4 16:11:54 EST 2016
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.4.4.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_TPOL_tree_
#define _JEventProcessor_TPOL_tree_

#include <JANA/JEventProcessor.h>

#include <DAQ/Df250WindowRawData.h>
#include <TTree.h>

class JEventProcessor_TPOL_tree:public jana::JEventProcessor{
public:
    JEventProcessor_TPOL_tree();
    ~JEventProcessor_TPOL_tree();
    const char* className(void){return "JEventProcessor_TPOL_tree";}

    TTree *TPOL;
    static const UInt_t nmax = 32;
    UInt_t nadc;
    UInt_t eventnum;                 ///< Event number
    UInt_t rocid[nmax];              ///< (from DDAQAddress) Crate number
    UInt_t slot[nmax];               ///< (from DDAQAddress) Slot number in crate
    UInt_t channel[nmax];            ///< (from DDAQAddress) Channel number in slot
    UInt_t itrigger[nmax];           ///< (from DDAQAddress) Trigger number for cases when this hit was read in a multi-event block (from DDAQAddress)
    UInt_t nsamples;                 ///< Number of samples in the waveform
    UInt_t waveform[nmax][150];      ///< array of samples in the waveform for the event
    UInt_t w_integral[nmax];         ///< Sum of all samples in the waveform
    UInt_t w_min[nmax];              ///< Minimum sample in the waveform
    UInt_t w_max[nmax];              ///< Maximum sample in the waveform
    UInt_t w_samp1[nmax];            ///< First sample in the waveform
    Double_t w_time[nmax];               ///< Half-pulse-height time in ns
    UInt_t sector[nmax];
    Double_t phi[nmax];
    Double_t E_lhit,E_rhit,t_lhit,t_rhit;
    UInt_t ntag;
    static const UInt_t ntag_max = 12;
    Double_t E_tag[ntag_max],t_tag[ntag_max];
    Bool_t is_tagm[ntag_max];

    int GetSector(int slot,int channel);
    double GetPhi(int sector);
    double GetPulseTime(const vector<uint16_t> waveform,double w_min,double w_max,double minpeakheight);

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TPOL_tree_
