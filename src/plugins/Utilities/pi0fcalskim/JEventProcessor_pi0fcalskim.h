// $Id$
//
//    File: JEventProcessor_pi0fcalskim.h
// Created: Mon Dec  1 14:57:11 EST 2014
// Creator: shepherd (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_pi0fcalskim_
#define _JEventProcessor_pi0fcalskim_

#include <JANA/JEventProcessor.h>
#include "evio_writer/DEventWriterEVIO.h"

class TTree;
class DFCALCluster;

#include <vector>

class JEventProcessor_pi0fcalskim:public jana::JEventProcessor{
 public:

  enum { kMaxHits = 500, kMaxClus = 20 };

  JEventProcessor_pi0fcalskim();
  ~JEventProcessor_pi0fcalskim();
  const char* className(void){return "JEventProcessor_pi0fcalskim";}

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  void writeClustersToRoot( const vector< const DFCALCluster* > clusVec );

  const DEventWriterEVIO* dEventWriterEVIO;		

  double MIN_MASS;
  double MAX_MASS;
  double MIN_E;
  double MIN_R;
  double MAX_DT;
  double MAX_ETOT;
  int    MIN_BLOCKS;

  int WRITE_ROOT;
  int WRITE_EVIO;

  TTree* m_tree;
  int m_nClus;
  int m_hit0[kMaxClus];
  float m_px[kMaxClus];
  float m_py[kMaxClus];
  float m_pz[kMaxClus];

  int m_nHit;
  int m_chan[kMaxHits];
  float m_e[kMaxHits];
};

#endif // _JEventProcessor_pi0fcalskim_

