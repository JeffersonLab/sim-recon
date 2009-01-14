// DTrackLSFitter.cc
//
/*******************************************************************************
To Do:

Allow for multiple tracks, either by declaring this as a one track
class, or by handling multiple hit lists.

Distinquish actions between Monte Carlo and real data.

*******************************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

#include "FDC/DFDCHit.h"
#include "FDC/DFDCPseudo_factory.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DMCTrackHit.h"
#include "DANA/DApplication.h"

#include "MyTrajectory.h"
#include "MyTrajectoryBfield.h"
#include "chisqMin.h"
#include "combinedResidFunc.h"
#include "globalGslFuncs.h"
#include "DTrackLSFitter.h"

#define TARGET_POSITION 65.0
#define THETA_FORWARD_CUT 1.178097
#define THETA_BACKWARD_CUT 1.963495
#define THREEPIOVER4 2.356194
#define MAX_TRAJS 1000
#define MAX_POINTS 10000
#define MAX_TRACKS 200
#define PI 3.141592653589793238

double radialDist2(const DCDCTrackHit *trkhit) {
  double x, y, r2;
  const DCDCWire *wire;
  wire = trkhit->wire;
  x = wire->origin.X();
  y = wire->origin.Y();
  r2 = x*x + y*y;
  return r2;
}

class compareTrackHits {
public:
  int operator() (const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) const
  {return radialDist2(hit1) < radialDist2(hit2);}
};

//------------------------------------------------------------------
// DTrackLSFitter 
//------------------------------------------------------------------
DTrackLSFitter::DTrackLSFitter():debug_level(1), ppEnd(5)
{
  ifstream configFile("config.txt");
  if (!configFile) {
    cout << "error: config file not found, terminating" << endl;
    exit(1);
  }
  char hddmFileName[128] = "fitter.hddm";
  ios = init_fitter_HDDM(hddmFileName);
  int debug_level_in;
  configFile >> debug_level_in;
  cout << "debug_level_in = " << debug_level_in << endl;
  if (debug_level_in >= -10 && debug_level_in <= 10) {
    debug_level = debug_level_in;
  }
}

//------------------------------------------------------------------
// ~DTrackLSFitter 
//------------------------------------------------------------------
DTrackLSFitter::~DTrackLSFitter()
{
  cout << "DTrackLSFitter destructor called\n";
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t DTrackLSFitter::init(void)
{
  signatureFile = new ofstream("fitter_signature.txt");
  return NOERROR;
}

//------------------------------------------------------------------
// fini
//------------------------------------------------------------------
jerror_t DTrackLSFitter::fini(void)
{
  signatureFile->close();
  delete signatureFile;
  close_fitter_HDDM(ios);
  return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t DTrackLSFitter::brun(JEventLoop *eventLoop, int runnumber)
{
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  lorentz_def = dapp->GetLorentzDeflections();
  JFactory_base *base = eventLoop->GetFactory("DFDCSegment");
  segment_factory = dynamic_cast<DFDCSegment_factory*>(base);
  return NOERROR;
}

//------------------------------------------------------------------
// erun 
//------------------------------------------------------------------
jerror_t DTrackLSFitter::erun()
{
 return NOERROR;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------

jerror_t DTrackLSFitter::evnt(JEventLoop *eventLoop, int eventnumber)
{

  status = DTRACKLSFITTER_UNDEFINED; // mark status as undefined

  if (debug_level > 1) cout << "----------- New Event "
			     << eventnumber << " -------------" << endl;

  fitPtr = NULL;
  //
  // set up hddm stuff
  //
  fitter_HDDM_t* hddm;
  hddm = make_fitter_HDDM();
  fitter_Events_t* eventsHddm;
  eventsHddm = make_fitter_Events(1);
  hddm->events = eventsHddm;
  eventsHddm->mult = 1;
  fitter_Event_t eventHddm;
  eventHddm.run = -1;
  eventHddm.number = eventnumber;
  fitter_Tracks_t* tracksHddm = make_fitter_Tracks(MAX_TRACKS);
  eventHddm.tracks = tracksHddm;
  tracksHddm->mult = 1;
  fitter_Track_t* trackHddm = tracksHddm->in;
  trackHddm->chisq = -1.0;
  fitter_Trajectorys_t* trajsHddm = make_fitter_Trajectorys(MAX_TRAJS);
  trackHddm->trajectorys = trajsHddm;
  //
  // get vectors from factories
  //
  vector<const DTrack*>tracks;
  eventLoop->Get(tracks, "THROWN");
  if (debug_level > 1) cout << "number of tracks = " << tracks.size() << endl;

  vector<const DMCTrackHit*> mctrackhits;
  eventLoop->Get(mctrackhits);
  if (debug_level > 1) cout << "number of mctrackhits = " << mctrackhits.size() << endl;

  pseudopoints.clear();
  //eventLoop->Get(pseudopoints);  
  tracks[0]->Get(pseudopoints);  

  trackhits.clear();
  //eventLoop->Get(trackhits);  
  tracks[0]->Get(trackhits);  
  
  thrown.clear();
  eventLoop->Get(thrown);

  size_fdc =  pseudopoints.size();
  if (debug_level > 1) cout << "number of pseudopoints = " << size_fdc << endl;
  trackHddm->nFDC = size_fdc;

  size_cdc = trackhits.size();
  if (debug_level > 1) cout << "number of trackhits = " << size_cdc << endl;
  trackHddm->nCDC = size_cdc;
  stable_sort(trackhits.begin(), trackhits.end(), compareTrackHits());
  //
  // write raw hit info to hddm file
  //
  writeFDCHitsHddm(eventHddm);
  writeCDCHitsHddm(eventHddm);
  //
  // check for minimum hits to attempt fit
  //
  if (size_fdc + size_cdc < 10) {
    trackHddm->fitStatus = DTRACKLSFITTER_FIT_NOT_ATTEMPTED;
    flush_fitter_HDDM(hddm, ios);
    return NOERROR;
  }
  //
  // find the starting parameters
  //
  setFitterStartParams();

#ifdef GRKUTA
  MyTrajectoryGrkuta trajectory(bfield, debug_level);
#else
  MyTrajectoryBfield trajectory(bfield, debug_level);
#endif

  // instantiate the residual function class
  combinedResidFunc prf(&pseudopoints, &trackhits, &trajectory, lorentz_def,
			debug_level);

  residFuncPtr = &prf; // set the global residual point to refer to this class
  fitPtr = new chisqMin(&prf, debug_level);
  HepVector ppStart(5), paramTrue(5);
  //
  // swim using MC truth points
  //
  trajectory.swimMC(mctrackhits);
  int tag_mc = 3;
  writeTrajectoryHddm(trajectory, tag_mc, trajsHddm);
  trajectory.clear();

  try {
    //
    // swim with MC parameters
    //
    paramTrue(1) = 0.0;
    paramTrue(2) = TARGET_POSITION;
    paramTrue(3) = thrown[0]->momentum().Theta();
    paramTrue(4) = thrown[0]->momentum().Phi();
    paramTrue(5) = 1.0/thrown[0]->momentum().Pt();
    trajectory.swim(paramTrue);
    int tag_true = 2;
    writeTrajectoryHddm(trajectory, tag_true, trajsHddm);
    trajectory.clear();
    vector<CDCHitDetails*> *CDCDetailsPtr;
    writeResidsHddm(paramTrue, CDCDetailsPtr, trajsHddm, prf, trajectory);
    prf.clearDetails();
    //
    // swim with initial parameters
    //
    ppStart(1) = xpInitial;
    ppStart(2) = zInitial;
    ppStart(3) = thetaInitial;
    ppStart(4) = phiInitial;
    ppStart(5) = ptinvInitial;
    trajectory.swim(ppStart); // swim with initial parameters
    int tag_initial = 0;
    writeTrajectoryHddm(trajectory, tag_initial, trajsHddm);
    trajectory.clear();
    //
    // do the fit
    //
    for (double df = 0.0; df < 1.0001; df += 1.0) {
      prf.setInnerResidFrac(df);
      fitPtr->setStartParams(ppStart);
      fitPtr->fit();
      fitPtr->getParams(ppEnd);
      // cout << df << ' ' << fitPtr->getChi2() << endl;
      ppStart = ppEnd;
    }
    //
    // store the results
    //
    trackHddm->chisq = fitPtr->getChi2();
    trajectory.swim(ppEnd); // swim with final parameters
    int tag_final = 1;
    writeTrajectoryHddm(trajectory, tag_final, trajsHddm);
    writeFitHddm(trajsHddm);
    *signatureFile << eventnumber << ' '
		   << size_fdc << ' '
		   << size_cdc << ' '
		   << thetaInitial << ' '
		   << phiInitial << ' '
		   << ppEnd(1) << ' ' 
		   << ppEnd(2) << ' '
		   << ppEnd(3) << ' '
		   << ppEnd(4) << ' '
		   << ppEnd(5) << ' '
		   << fitPtr->getChi2() << endl;
    // get details of hits on track
    writeResidsHddm(ppEnd, CDCDetailsPtr, trajsHddm, prf, trajectory);

    status = DTRACKLSFITTER_NOMINAL;
  } catch (int code) {

    cout << "==========fit error = " << code << "===========" << endl;
    status = DTRACKLSFITTER_EXCEPTION_THROWN;
    cout << "= at event " << eventnumber << endl;
    if (debug_level >= 3) {
      cout << "= trajectory" << endl;
      trajectory.print();
    }

  }

  // clean up
  prf.clearDetails();
  prf.setStoreDetails(false);
  trajectory.clear();

  if (debug_level >= 3) {
    cout << "FDCHitDetails: " << FDCHitDetails::getInstanceCount() << endl;
    cout << "CDCHitDetails: " << CDCHitDetails::getInstanceCount() << endl;
  }

  trackHddm->fitStatus = status;
  //
  // load hddm event into buffer
  //
  eventsHddm->in[0] = eventHddm;
  flush_fitter_HDDM(hddm, ios);

  if (fitPtr != NULL) {
    delete fitPtr;
    fitPtr = NULL;
  }

  return NOERROR;
 
}

HepVector DTrackLSFitter::getParams() {
  return ppEnd;
}

double DTrackLSFitter::getChiSquared() {
  double chisq = -1.0;
  if (fitPtr) chisq = fitPtr->getChi2();
  return chisq;
}

int DTrackLSFitter::getSizeFDC() {
  return size_fdc;
}

int DTrackLSFitter::getSizeCDC() {
  return size_cdc;
}

void DTrackLSFitter::setFitterStartParams() {
  zInitial = TARGET_POSITION;
  ptinvInitial = 0.0;
  if (size_fdc > 1) {
    double xfirst, yfirst, zfirst, xlast, ylast, zlast;
    xfirst = (pseudopoints[0])->x;
    yfirst = (pseudopoints[0])->y;
    zfirst = (pseudopoints[0])->wire->origin(2);
    xlast = (pseudopoints[size_fdc - 1])->x;
    ylast = (pseudopoints[size_fdc - 1])->y;
    zlast = (pseudopoints[size_fdc - 1])->wire->origin(2);
    double rfirst = sqrt(xfirst*xfirst + yfirst*yfirst);
    double delta_z = zlast - zfirst;
    double rPerp = sqrt((xlast - xfirst)*(xlast - xfirst) + (ylast - yfirst)*(ylast - yfirst));
    thetaInitial = atan2(rPerp, delta_z);
    phiInitial = atan2(ylast - yfirst, xlast - xfirst);
    xpInitial = (xfirst - yfirst/tan(phiInitial))*sin(phiInitial);
    zInitial = zfirst - sqrt(rfirst*rfirst - xpInitial*xpInitial)/tan(thetaInitial);
  } else if (size_cdc > 0) {

    float thetaInitialTrue = thrown[0]->momentum().Theta();
    if (thetaInitialTrue < PIOVER2) {
      thetaInitial = thetaInitialTrue - 5*PI/180;
    } else {
      thetaInitial = thetaInitialTrue + 5*PI/180;
    }

    const DCDCTrackHit *trkhit_0 = trackhits[0];
    const DCDCWire *wire_0 = trkhit_0->wire;
    const DCDCTrackHit *trkhit_n = trackhits[size_cdc - 1];
    const DCDCWire *wire_n = trkhit_n->wire;
    double deltaX = wire_n->origin.X() - wire_0->origin.X();
    double deltaY = wire_n->origin.Y() - wire_0->origin.Y();
    phiInitial = atan2(deltaY, deltaX);
    double alpha = phiInitial - PIOVER2;
    xpInitial = wire_0->origin.X()*cos(alpha) + wire_0->origin.Y()*sin(alpha); 
  } else {
    int code = 32;
    throw code;
  }
}

void DTrackLSFitter::writeTrajectoryHddm
(MyTrajectoryBfield &trajectory, int tag, fitter_Trajectorys_t *trajsHddm) {
  // get the next trajectory pointer
  fitter_Trajectory_t* trajHddm = &(trajsHddm->in[trajsHddm->mult++]);
  // set the tag
  trajHddm->index = tag;
  // create the parameters structure
  unsigned int nparams = trajectory.getNumberOfParams();
  fitter_Parameters_t* parametersHddm = make_fitter_Parameters(nparams);
  parametersHddm->mult = nparams; // set the number of parameters for output
  // put parameters structure in the trajectory
  trajHddm->parameters = parametersHddm;
  // get the parameters and put them in the parameters structure
  HepVector params = trajectory.getParams();
  fitter_Parameter_t* thisParamHddmPtr;
  for (int i = 0; i < (int)nparams; i++) {
    thisParamHddmPtr = parametersHddm->in + i;
    thisParamHddmPtr->label = i;
    thisParamHddmPtr->value = params(i+1);
  }
  // do the points
  fitter_Points_t* pointsHddm = make_fitter_Points(MAX_POINTS);
  trajHddm->points = pointsHddm;
  vector<HepLorentzVector*> *traj;
  traj = trajectory.getTrajectory();
  HepLorentzVector trajPoint;
  fitter_Point_t* pointHddm;
  for (unsigned int i = 0; i < traj->size(); i++) {
    trajPoint = *((*traj)[i]);
    pointHddm = &(pointsHddm->in[pointsHddm->mult++]);
    pointHddm->x = trajPoint.getX();;
    pointHddm->y = trajPoint.getY();
    pointHddm->z = trajPoint.getZ();
    pointHddm->t = trajPoint.getT();
  }
}

void DTrackLSFitter::writeCDCDetailsHddm(vector<CDCHitDetails*> *CDCDetailsPtr, fitter_Trajectorys_t *trajsHddm) {
  int nCDC = CDCDetailsPtr->size();
  fitter_Trajectory_t* trajHddm = &(trajsHddm->in[trajsHddm->mult - 1]);
  fitter_Residuals_t* residsHddm = make_fitter_Residuals(nCDC);
  trajHddm->residuals = residsHddm;
  double dist, doca, resid;
  fitter_Residual_t* thisResidPtr;
  fitter_ResidInfoCdc_t* infoCdcPtr;
  for (int i = 0; i < nCDC; i++) {
    dist = ((*CDCDetailsPtr)[i])->dist;
    doca = ((*CDCDetailsPtr)[i])->doca;
    resid = doca - dist;
    thisResidPtr = residsHddm->in + i;
    thisResidPtr->detector_index = i;
    thisResidPtr->value = resid;
    infoCdcPtr = make_fitter_ResidInfoCdc();
    thisResidPtr->residInfoCdc = infoCdcPtr;
    infoCdcPtr->id = i;
    infoCdcPtr->dist = dist;
    infoCdcPtr->doca = doca;
    infoCdcPtr->xWire = ((*CDCDetailsPtr)[i])->posWire(1);
    infoCdcPtr->yWire = ((*CDCDetailsPtr)[i])->posWire(2);
    infoCdcPtr->zWire = ((*CDCDetailsPtr)[i])->posWire(3);
    infoCdcPtr->xTraj = ((*CDCDetailsPtr)[i])->poca.getX();
    infoCdcPtr->yTraj = ((*CDCDetailsPtr)[i])->poca.getY();
    infoCdcPtr->zTraj = ((*CDCDetailsPtr)[i])->poca.getZ();
    infoCdcPtr->tTraj = ((*CDCDetailsPtr)[i])->poca.getT();
  }
  residsHddm->mult = nCDC;
}

int DTrackLSFitter::getStatus(void) {
  return status;
}

void DTrackLSFitter::writeFDCHitsHddm(fitter_Event_t &event) {
  fitter_Fdcdata_t* fdcdataPtr = make_fitter_Fdcdata();
  event.fdcdata = fdcdataPtr;
  fitter_Fdchits_t* fdchitsPtr = make_fitter_Fdchits(size_fdc);
  fdcdataPtr->fdchits = fdchitsPtr;
  fdchitsPtr->mult = size_fdc;
  for (int i = 0; i < size_fdc; i++) {
    cout << i << ' ' << (pseudopoints[i])->x << ' ' << (pseudopoints[i])->y << ' ' << (pseudopoints[i])->wire->origin(2) << endl;
    fdchitsPtr->in[i].id = i;
    fdchitsPtr->in[i].x = (pseudopoints[i])->x;
    fdchitsPtr->in[i].y = (pseudopoints[i])->y;
    fdchitsPtr->in[i].z = (pseudopoints[i])->wire->origin(2);
  }
}

void DTrackLSFitter::writeCDCHitsHddm(fitter_Event_t &event) {
  fitter_Cdcdata_t* cdcdataPtr = make_fitter_Cdcdata();
  event.cdcdata = cdcdataPtr;
  fitter_Cdchits_t* cdchitsPtr = make_fitter_Cdchits(size_cdc);
  cdcdataPtr->cdchits = cdchitsPtr;
  cdchitsPtr->mult = size_cdc;
  //double x, y, phi_wire;
  double phi, theta;
  for (int i = 0; i < size_cdc; i++) {
    const DCDCTrackHit *trkhit = trackhits[i];
    const DCDCWire *wire = trkhit->wire;
    //x = wire->origin.X();
    //y = wire->origin.Y();
    //phi_wire = atan2(y, x);
    //phi = phi_wire + PIOVER2;
    phi = atan2(wire->udir.y(), wire->udir.x());
    //theta = wire->stereo;
    theta = acos(wire->udir.z());
    cdchitsPtr->in[i].id = i;
    cdchitsPtr->in[i].x = wire->origin.x();
    cdchitsPtr->in[i].y = wire->origin.y();
    cdchitsPtr->in[i].z= wire->origin.z();
    cdchitsPtr->in[i].L = wire->L;
    cdchitsPtr->in[i].theta = theta;
    cdchitsPtr->in[i].phi = phi;
    cdchitsPtr->in[i].tdrift = trkhit->tdrift;
    cdchitsPtr->in[i].ring = wire->ring;
    cdchitsPtr->in[i].straw = wire->straw;
  }
}

void DTrackLSFitter::writeResidsHddm(const HepVector &params,
				      vector<CDCHitDetails*>* &CDCDetailsPtr,
				      fitter_Trajectorys_t *trajsHddm,
				      combinedResidFunc &prf,
				      MyTrajectoryBfield &trajectory) {
  trajectory.clear();
  prf.setStoreDetails(true);
  int dummy_data;
  HepVector dummy_vector(size_fdc + size_cdc);
  prf.resid(&params, &dummy_data, &dummy_vector);
  trajsHddm->in[trajsHddm->mult - 1].chisq = prf.getChiSquared();
  CDCDetailsPtr = prf.getCDCDetails();
  writeCDCDetailsHddm(CDCDetailsPtr, trajsHddm);
  prf.setStoreDetails(false);
  trajectory.clear();
}

void DTrackLSFitter::writeFitHddm(fitter_Trajectorys_t* trajsHddm) {
  int nparams = fitPtr->getP();
  HepMatrix covar(nparams, nparams);
  covar = fitPtr->getCovar();
  // get the pointer to the current trajectory
  fitter_Trajectory_t* trajHddmPtr = trajsHddm->in + trajsHddm->mult - 1;
  // create the fit structure
  fitter_Fit_t* fitHddmPtr = make_fitter_Fit();
  // store it in the current trajectory
  trajHddmPtr->fit = fitHddmPtr;
  // store attributes for fit
  fitHddmPtr->type = 1;
  fitHddmPtr->iterations = fitPtr->getIter();
  // create the matrix structure
  fitter_CovarianceMatrix_t* matrixHddmPtr = make_fitter_CovarianceMatrix();
  // store it in the fit
  fitHddmPtr->covarianceMatrix = matrixHddmPtr;
  // store the number of parameters
  matrixHddmPtr->params = nparams;
  // create the elements array
  fitter_CovarianceElements_t* elementsHddmPtr = make_fitter_CovarianceElements(100);
  // store it in the matrix
  matrixHddmPtr->covarianceElements = elementsHddmPtr;
  // loop over all matrix elements
  int ielement = 0;
  for (int i = 1; i <= nparams; i++) {
    for (int j = i; j <= nparams; j++) {
      // get the point to this element
      fitter_CovarianceElement_t* elementHddmPtr = elementsHddmPtr->in + ielement;
      // fill the element with info
      elementHddmPtr->row = i;
      elementHddmPtr->column = j;
      elementHddmPtr->value = covar(i, j);
      // store the complete element in the matrix
      elementsHddmPtr->in[ielement] = *elementHddmPtr;
      // increment the element counter
      ielement++;
    }
  }
  // store the element count in the matrix
  elementsHddmPtr->mult = ielement;
}
