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

#include <FDC/DFDCHit.h>
#include <FDC/DFDCPseudo_factory.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <DANA/DApplication.h>

#include "MyTrajectory.h"
#include "MyTrajectoryBfield.h"
#include "chisqMin.h"
#include "combinedResidFunc.h"
#include "globalGslFuncs.h"
#include "DTrackLSFitter.h"
#include "DFactoryGeneratorLSLM.h"

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

// Routine used If we're a plugin
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddFactoryGenerator(new DFactoryGeneratorLSLM());
}
} // "C"


//------------------------------------------------------------------
// DTrackLSFitter 
//------------------------------------------------------------------
DTrackLSFitter::DTrackLSFitter(JEventLoop *loop):DTrackFitter(loop),debug_level(1), ppEnd(5)
{
	// Get the DApplication pointer so we can get pointers to the Lorentz deflections
	// table for the FDC and the DFDCSegment factory pointer.
	// NOTE: The bfield member is supplied by the DTrackFitter base class and
	// set in the DTrackFitter(loop) constructor.
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
  lorentz_def = dapp->GetLorentzDeflections();
  JFactory_base *base = loop->GetFactory("DFDCSegment");
  segment_factory = dynamic_cast<DFDCSegment_factory*>(base);
}

//------------------------------------------------------------------
// ~DTrackLSFitter 
//------------------------------------------------------------------
DTrackLSFitter::~DTrackLSFitter()
{
  cout << "DTrackLSFitter destructor called\n";
}

//------------------------------------------------------------------
// ChiSq 
//------------------------------------------------------------------
double DTrackLSFitter::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr, vector<pull_t> *pulls_ptr)
{
	// This will need to be filled in by Mark
	if(chisq_ptr)*chisq_ptr=0.0;
	if(dof_ptr)*dof_ptr=0.0;
	return 0.0; // return the chisq/Ndof
}

//------------------------------------------------------------------
// FitTrack 
//------------------------------------------------------------------
DTrackFitter::fit_status_t DTrackLSFitter::FitTrack(void)
{

	// Initialize result to indicate fit in progress
	this->fit_status = kFitNotDone;

  status = DTRACKLSFITTER_UNDEFINED; // mark status as undefined

  fitPtr = NULL;

	// The CDC and FDC hits are already filled out by the DTrackFitter base class
	// This code used different names though so we just set a couple of aliases.
	vector<const DCDCTrackHit*> &trackhits = cdchits;
	vector<const DFDCPseudo*> &pseudopoints = fdchits;

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
  try {

    fitPtr = new chisqMin(&prf, debug_level);
    HepVector ppStart(5), paramTrue(5);
  //
  // swim using MC truth points
  //

//  trajectory.swimMC(mctrackhits);
    trajectory.clear();


    //
    // swim with MC parameters
    //
//    paramTrue(1) = 0.0;
//    paramTrue(2) = TARGET_POSITION;
//    paramTrue(3) = thrown[0]->momentum().Theta();
//    paramTrue(4) = thrown[0]->momentum().Phi();
//    paramTrue(5) = 1.0/thrown[0]->momentum().Pt();
//    trajectory.swim(paramTrue);

//    trajectory.clear();
 
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
    trajectory.swim(ppEnd); // swim with final parameters
 
	 /*
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
		*/
		
    // get details of hits on track

    status = DTRACKLSFITTER_NOMINAL;
  } catch (int code) {

    cout << "==========fit error = " << code << "===========" << endl;
    status = DTRACKLSFITTER_EXCEPTION_THROWN;
    cout << "= at event " << loop->GetJEvent().GetEventNumber() << endl;
    if (debug_level >= 3) {
      cout << "= trajectory" << endl;
      trajectory.print();
    }
	 this->fit_status = kFitFailed;

  }

  // clean up
  prf.clearDetails();
  prf.setStoreDetails(false);
  trajectory.clear();

  if (debug_level >= 3) {
    cout << "FDCHitDetails: " << FDCHitDetails::getInstanceCount() << endl;
    cout << "CDCHitDetails: " << CDCHitDetails::getInstanceCount() << endl;
  }

  if (status == DTRACKLSFITTER_NOMINAL) {

	// Copy results into DTrackFitter data members
	double r0 = ppEnd(1);
	double z0 = ppEnd(2);
	double theta = ppEnd(3);
	double phi = ppEnd(4);
	double ptot = fabs(1.0/ppEnd(5)/sin(theta));
	
	double x0 = r0*cos(phi);
	double y0 = r0*sin(phi);
	DVector3 pos(x0, y0, z0);
	DVector3 mom;
	mom.SetMagThetaPhi(ptot, theta, phi);
	fit_params.setMomentum(mom);
	fit_params.setPosition(pos);
	fit_params.setCharge(ppEnd(5)>=0.0 ? 1.0:-1.0);
	
	this->chisq = fitPtr->getChi2();
	this->Ndof = fitPtr->getN() - fitPtr->getP();
	if(this->fit_status==kFitNotDone)this->fit_status = kFitSuccess; // honor it if fit_status has already been set (e.g. to kFitFailed)
	this->cdchits_used_in_fit = cdchits;
	this->fdchits_used_in_fit = fdchits;
  }

  if (fitPtr != NULL) {
    delete fitPtr;
    fitPtr = NULL;
  }

  return this->fit_status;
 
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
	xpInitial = input_params.position().Perp(); // distance from beam line at point of closest approach
	zInitial = input_params.position().Z();
	thetaInitial = input_params.momentum().Theta();
	phiInitial = input_params.momentum().Phi();
	ptinvInitial = 1.0/input_params.momentum().Perp();

#if 0
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
#endif
}


int DTrackLSFitter::getStatus(void) {
  return status;
}

