// MyProcessor.cc
//

#include <iostream>
#include <vector>
#include <string>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

#include "MyProcessor.h"
#include "FDC/DFDCHit.h"
#include "FDC/DFDCPseudo_factory.h"
#include "FDC/DFDCSegment_factory.h"
#include "CDC/DCDCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DMCTrackHit.h"
#include "DANA/DApplication.h"

#include "MyTrajectory.h"
#include "MyTrajectoryBfield.h"
#include "chisqMin.h"
#include "combinedResidFunc.h"

#define TARGET_POSITION 65.0
#define THETA_FORWARD_CUT 1.178097
#define THETA_BACKWARD_CUT 1.963495
#define THREEPIOVER4 2.356194

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
}

//------------------------------------------------------------------
// ~MyProcessor 
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{
  cout << "MyProcessor destructor called\n";
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
  rootfile = new TFile("fitter.root", "RECREATE");
  nFitter = new TNtuple("nFitter", "fitter results", "event:status:xp0:z0:theta:phi:ptinv:chisq:nFDC:nCDC:xp0True:z0True:thetaTrue:phiTrue:ptinvTrue");
  lsfitter.init();
  return NOERROR;
}

//------------------------------------------------------------------
// fini
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
  cout << "fini called\n";
  rootfile->Write();
  return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{
  jerror_t status = lsfitter.brun(eventLoop, runnumber);
  return status;
}

//------------------------------------------------------------------
// erun 
//------------------------------------------------------------------
jerror_t MyProcessor::erun()
{
  jerror_t status = lsfitter.erun();
  return status;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------

jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{

  vector<const DMCThrown*>thrown;
  eventLoop->Get(thrown);
  double thetaTrue = thrown[0]->momentum().Theta();
  double phiTrue = thrown[0]->momentum().Phi();
  double z0True = TARGET_POSITION;
  double xp0True = 0.0;
  double ptinvTrue = 1.0/(thrown[0]->momentum().Mag()*sin(thetaTrue));

  //if (eventnumber != 1) return NOERROR;

  jerror_t status = lsfitter.evnt(eventLoop, eventnumber);
  int lsfitterStatus = lsfitter.getStatus();

  // fill the ntuple
  HepVector ppEnd(5);
  ppEnd = lsfitter.getParams();
  double chisq;
  chisq = lsfitter.getChiSquared(); // this will be garbage for now
  int size_fdc = lsfitter.getSizeFDC();
  int size_cdc = lsfitter.getSizeCDC();
  nFitter->Fill(eventnumber, (float)lsfitterStatus, ppEnd(1), ppEnd(2), ppEnd(3), ppEnd(4), ppEnd(5),
		chisq, size_fdc, size_cdc, xp0True, z0True,
		thetaTrue, phiTrue, ptinvTrue);
  return status;
 
}

