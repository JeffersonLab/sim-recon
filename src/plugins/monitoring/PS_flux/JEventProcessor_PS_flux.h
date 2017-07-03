// $Id$
//
//    File: JEventProcessor_PS_flux.h
//
//

#ifndef _JEventProcessor_PS_flux_
#define _JEventProcessor_PS_flux_

#include "TRandom3.h"

#include <JANA/JEventProcessor.h>

#include "PAIR_SPECTROMETER/DPSGeometry.h"
#include "ANALYSIS/DTreeInterface.h"
#include "DAQ/DBeamCurrent.h"
#include "DAQ/DBeamCurrent_factory.h"

class JEventProcessor_PS_flux:public jana::JEventProcessor{
public:
    JEventProcessor_PS_flux();
    ~JEventProcessor_PS_flux();
    const char* className(void){return "JEventProcessor_PS_flux";}

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.

    DBeamCurrent_factory *dBeamCurrentFactory;
    double dBeamBunchPeriod;
    double t_start;
    double t_end;
    double t_fiducial;

    //TREE
    DTreeInterface* dTreeInterface;
    //thread_local: Each thread has its own object: no lock needed
    //important: manages it's own data internally: don't want to call new/delete every event!
    static thread_local DTreeFillData dTreeFillData;

    const DPSGeometry *dPSGeom;
    int geomModuleColumn[8][2] = {{110, 145}, {90, 115}, {73, 93}, {56, 76}, {40, 60}, {24, 45}, {8, 28}, {0, 12}}; // {min, max}
    TRandom3 *dRandom;
};

#endif // _JEventProcessor_PS_flux_

