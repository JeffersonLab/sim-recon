
#ifndef _JEventProcessor_FCAL_TimingOffsets_
#define _JEventProcessor_FCAL_TimingOffsets_

#include <JANA/JEventProcessor.h>
#include <FCAL/DFCALGeometry.h>
#include "ANALYSIS/DTreeInterface.h"
#include <thread>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALShower.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250PulsePedestal.h"
#include "DAQ/DEPICSvalue.h"
#include "TRIGGER/DL1Trigger.h"
#include <DANA/DStatusBits.h>
#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "HDGEOMETRY/DGeometry.h"
#include "DANA/DApplication.h"
#include <TTree.h>


class DFCALHit;

class JEventProcessor_FCAL_TimingOffsets:public jana::JEventProcessor{
 public:

  JEventProcessor_FCAL_TimingOffsets();
  ~JEventProcessor_FCAL_TimingOffsets();
  const char* className(void){return "JEventProcessor_FCAL_TimingOffsets";}

 private:

  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  double m_targetZ;

  int m_column;
  int m_row;
  int m_chan;


  const DFCALGeometry* m_fcalGeom;
  double m_FCALfront;
 DVector3 dTargetCenter;
 DVector3 dFCALblockcenter;


  
  
 
};

#endif // _JEventProcessor_FCAL_TimingOffsets_

