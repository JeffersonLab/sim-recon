// -----------------------------------------
// DEventProcessor_truth_dirc.h
// created on: 07.04.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DEVENTPROCESSOR_TRUTH_DIRC_H_
#define DEVENTPROCESSOR_TRUTH_DIRC_H_

#include <iostream>
#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>
#include <DANA/DApplication.h>
#include <HDGEOMETRY/DGeometry.h>

using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <PID/DKinematicData.h>
#include <PID/DBeamPhoton.h>
#include <DIRC/DDIRCTruthBarHit.h>
#include <DIRC/DDIRCTruthPmtHit.h>

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectoryFile.h>
#include <TThread.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TClonesArray.h>


class DEventProcessor_truth_dirc: public JEventProcessor {

public:
  DEventProcessor_truth_dirc();
  ~DEventProcessor_truth_dirc();

  pthread_mutex_t mutex;

private:
  jerror_t init(void);
  jerror_t brun(jana::JEventLoop *loop, int32_t runnumber);
  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);
  jerror_t erun(void);
  jerror_t fini(void); // called after last event

  TH1F *hTruthBarHitBar;
  TH2F *hTruthBarHitXY;
  TH2F *hTruthPmtHitZY_North, *hTruthPmtHitZY_South;
  TH2F *hTruthPmtHit_North, *hTruthPmtHit_South;
  TH2F *hTruthPixelHit_North, *hTruthPixelHit_South;

};

#endif /* DEVENTPROCESSOR_TRUTH_DIRC_H_ */
