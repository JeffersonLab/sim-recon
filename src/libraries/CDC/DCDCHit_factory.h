// $Id$
//
//    File: DCDCHit_factory.h
// Created: Tue Aug  6 11:29:56 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DCDCHit_factory_
#define _DCDCHit_factory_

#include <JANA/JFactory.h>
#include <DAQ/Df125CDCPulse.h>
#include <TTAB/DTranslationTable.h>

#include "DCDCHit.h"

using namespace std;
using namespace jana;

class DCDCHit_factory: public jana::JFactory<DCDCHit>{
 public:
  DCDCHit_factory(){};
  ~DCDCHit_factory(){};
  
  // we need to store information on the hits with respect to their readout channels in order to look for correlated hits
  struct cdchit_info_t{
    uint32_t rocid;
    uint32_t slot;
    uint32_t connector;
    
    double time;
    double max;
    
    inline bool operator==(const struct cdchit_info_t &rhs) const {
      return (rocid==rhs.rocid) && (slot==rhs.slot) && (connector==rhs.connector);
    }
  };
  
  
  int RemoveCorrelationHits;
  double CorrelationHitsCut;
  double CorrelatedHitPeak;
  int Disable_CDC_TimingCuts;

  // timing cut limits
  double LowTCut;
  double HighTCut;
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  vector<const DTranslationTable *> ttab;
};

#endif // _DCDCHit_factory_

