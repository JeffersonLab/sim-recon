// Author: David Lawrence  June 25, 2004
//
//
// DTrackLSFitter.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#ifndef _DTRACKLSFITTER_H_
#define _DTRACKLSFITTER_H_

#define DTRACKLSFITTER_UNDEFINED -2
#define DTRACKLSFITTER_FIT_NOT_ATTEMPTED -1
#define DTRACKLSFITTER_NOMINAL 0
#define DTRACKLSFITTER_EXCEPTION_THROWN 1

#include <JANA/JEventProcessor.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <FDC/DFDCSegment_factory.h>
#include <CDC/DCDCTrackHit.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackFitter.h>
#include <TH1F.h>
#include <TFile.h>
#include <TNtuple.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>
#include "chisqMin.h"
#include "hitDetails.h"
#include "combinedResidFunc.h"
#include "MyTrajectoryGrkuta.h"

class DReferenceTrajectory;

class DTrackLSFitter:public DTrackFitter
{
 public:
  DTrackLSFitter(JEventLoop *loop);
  ~DTrackLSFitter();
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);			///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);		///< Called every event.
  int eventNo;
  jerror_t erun(void);
  jerror_t fini(void);
  int debug_level;
  HepVector getParams();
  double getChiSquared();
  int getSizeFDC();
  int getSizeCDC();
  int getStatus();
  //void writeFDCHitsHddm(fitter_Event_t &eventHddm);

		// Virtual methods from TrackFitter base class
		string Name(void) const {return string("MMI");}
		fit_status_t FitTrack(void);
		double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);
 private:
  
  //const DMagneticFieldMap *bfield; // supplied by DTrackFitter base class
  const DLorentzDeflections *lorentz_def;
  
  DFDCSegment_factory *segment_factory;
  //ofstream *signatureFile; 
  //ifstream *configFile; 
  HepVector ppEnd;
  chisqMin *fitPtr;
  int size_fdc, size_cdc;
  //vector<const DFDCPseudo*>pseudopoints;
  //vector<const DCDCTrackHit*>trackhits;
  //vector<const DMCThrown*>thrown;
  double xpInitial, zInitial, thetaInitial, phiInitial, ptinvInitial;
  void setFitterStartParams();
  //fitter_iostream_t* ios;
  //void writeTrajectoryHddm(MyTrajectoryBfield &traj, int tag, fitter_Trajectorys_t *trajsHddm);
  //void writeCDCDetailsHddm(vector<CDCHitDetails*> *CDCDetailsPtr, fitter_Trajectorys_t *trajsHddm);
  //void writeCDCHitsHddm(fitter_Event_t &event);
  int status; // status code
  //void writeResidsHddm(const HepVector &params,
	//		vector<CDCHitDetails*>* &CDCDetailsPtr,
	//		fitter_Trajectorys_t *trajsHddm,
	//		combinedResidFunc &prf,
	///		MyTrajectoryBfield &trajectory);
  //void writeFitHddm(fitter_Trajectorys_t* trajsHddm);
};

#endif // _DTRACKLSFITTER_H_
