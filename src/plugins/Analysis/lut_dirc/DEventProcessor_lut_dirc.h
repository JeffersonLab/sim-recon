// -----------------------------------------
// DEventProcessor_lut_dirc.h
// created on: 29.11.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DEVENTPROCESSOR_LUT_DIRC_H_
#define DEVENTPROCESSOR_LUT_DIRC_H_

#include <iostream>
#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>
using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <PID/DKinematicData.h>
#include <PID/DBeamPhoton.h>
#include <DIRC/DDIRCTruthBarHit.h>
#include <DIRC/DDIRCTruthPmtHit.h>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TThread.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TClonesArray.h>

#include "DrcLutNode.h"

class DEventProcessor_lut_dirc: public JEventProcessor {

public:
  DEventProcessor_lut_dirc();
  ~DEventProcessor_lut_dirc();

  TClonesArray *fLut[48];
  DrcLutNode *fLutNode;
  TTree *fTree;

  pthread_mutex_t mutex;

private:
  jerror_t init(void);
  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);
  jerror_t erun(void);
  jerror_t fini(void); // called after last event

};

#endif /* DEVENTPROCESSOR_LUT_DIRC_H_ */
