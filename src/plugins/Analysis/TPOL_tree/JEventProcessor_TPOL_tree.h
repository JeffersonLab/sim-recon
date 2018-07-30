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
//#include <TTree.h>

#include "PAIR_SPECTROMETER/DPSGeometry.h"
#include "ANALYSIS/DTreeInterface.h"
#include "DAQ/DBeamCurrent.h"
#include "DAQ/DBeamCurrent_factory.h"

class JEventProcessor_TPOL_tree:public jana::JEventProcessor{
public:
    JEventProcessor_TPOL_tree();
    ~JEventProcessor_TPOL_tree();
    const char* className(void){return "JEventProcessor_TPOL_tree";}

    int GetSector(int slot,int channel);
    double GetPhi(int sector);
    double GetPulseTime(const vector<uint16_t> waveform,double w_min,double w_max,double minpeakheight);

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.

    DBeamCurrent_factory *dBeamCurrentFactory;
    double dBeamBunchPeriod;

    DTreeInterface* dTreeInterface;
    static thread_local DTreeFillData dTreeFillData;

    int geomModuleColumn[8][2] = {{110, 145}, {90, 115}, {73, 93}, {56, 76}, {40, 60}, {24, 45}, {8, 28}, {0, 12}};

};

#endif // _JEventProcessor_TPOL_tree_
