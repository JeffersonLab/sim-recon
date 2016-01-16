// $Id$
//
//    File: JEventProcessor_FCAL_online.h
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FCAL_online_
#define _JEventProcessor_FCAL_online_

#include <JANA/JEventProcessor.h>

class TH1D;
class TH1I;
class TH2I;
class TH2F;
class TH1F;
class TProfile;

class JEventProcessor_FCAL_online:public jana::JEventProcessor{
 public:
  JEventProcessor_FCAL_online();
  ~JEventProcessor_FCAL_online();
  const char* className(void){return "JEventProcessor_FCAL_online";}

 private:

  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  double m_targetZ;

  TH1D *fcal_num_events;

  TH1I* m_digInt;
  TH1I* m_digCoarseT;
  TProfile* m_digCoarseTChan;
  TH1I* m_digPreciseT;
  TProfile* m_digPreciseTChan;
  TH1I* m_digT;
  TH1I* m_digT0;
  TH1I* m_digTmT0;
  TH2F* m_digTmT02D;
  TH1I* m_digPed;
  TProfile* m_digPedChan;
  TH2F* m_digPed2D;
  TH2F* m_digPedSq2D;
  TH1I* m_digQual;
  TH1I* m_digNUnder;
  TH1I* m_digNOver;
  TH1I* m_digN;
  TH1I* m_digPeakV;
  TH2F* m_digPeakV2D;
  TH2F* m_digOcc2D;
  TH2I* m_digIntVsPeak;
  TH1I* m_digIntToPeak;

  TH1I* m_hitN;
  TH1I* m_hitE;
  TH1I* m_hitETot;
  TH1I* m_hitT;
  TH1I* m_hitT0;
  TH1I* m_hitTmT0;
  TH2F* m_hitE2D;
  TH2F* m_hitTmT02D;
  TH2F* m_hitTmT0Sq2D;
  TH2F* m_hitOcc2D;
  
  TH1I* m_clusN;
  TH1I* m_clusE;
  TH1I* m_clusETot;
  TH1I* m_clusT;
  TH1I* m_clusT0;
  TH1I* m_clusTmT0;
  TH2I* m_clusXYHigh;
  TH2I* m_clusXYLow;
  TH1I* m_clusPhi;
  TH1I* m_clus2GMass;

  TH1I* m_show2GMass;
  TH2I* m_showZvsE;
  TH2I* m_showECorVsE;
  TH2I* m_showTsMTcVsZ;

};

#endif // _JEventProcessor_FCAL_online_

