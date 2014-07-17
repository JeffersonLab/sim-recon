// $Id$
//
//    File: DTrackCandidate_factory_FDCCathodes.h
// Created: Tue Nov  6 13:37:08 EST 2007
// Creator: staylor (on Linux ifarml1.jlab.org 2.4.21-47.0.1.ELsmp i686)
//

#ifndef _DTrackCandidate_factory_FDCCathodes_
#define _DTrackCandidate_factory_FDCCathodes_

#include <JANA/JFactory.h>
#include "DTrackCandidate.h"
#include <DMatrix.h>
#include "FDC/DFDCSegment_factory.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DMagneticFieldStepper.h"
#include <TH1F.h>
#include <TH2F.h>
#include "DHelicalFit.h"

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="CORE.png" width="100">
///	</A>
/// \endhtmlonly

/// Find track candidates in the Forward Drift Chambers (FDC).
/// This will try to form candidates from hits found in the FDC.
/// Roughly, clusters of hits in single packages (6 wire planes)
/// are found and fit. The fit is propagated to adjacent packages
/// to link them up with hit clusters there. Eventually, all hits
/// are used to fit using a Reimann method to obtain the parameters
/// of the candidate.
///
/// These candidates will be merged with those from the FDC in
/// the DTrackCandidate_factory class.

class DTrackCandidate_factory_FDCCathodes:public JFactory<DTrackCandidate>{
 public:
  DTrackCandidate_factory_FDCCathodes(){
    DEBUG_HISTS = false;
    //DEBUG_HISTS = true;
  };
  ~DTrackCandidate_factory_FDCCathodes(){};
  const char* Tag(void){return "FDCCathodes";}
  
 private:  
  const DMagneticFieldMap *bfield;
  DMagneticFieldStepper *stepper;
		
  //jerror_t init(void);						///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  jerror_t GetPositionAndMomentum(const DFDCSegment *segment);
  jerror_t GetPositionAndMomentum(DFDCSegment *segment,
				  DVector3 &pos, DVector3 &mom);
  jerror_t GetPositionAndMomentum(DVector3 &pos,DVector3 &mom);
  jerror_t GetPositionAndMomentum(vector<const DFDCSegment *>segments,
				  DVector3 &pos,DVector3 &mom);
  
  double GetCharge(const DVector3 &pos,const DFDCSegment *segment);
  double GetCharge(const DVector3 &pos,vector<const DFDCSegment *>segments);

  double DocaSqToHelix(const DFDCPseudo *hit);
  DFDCSegment *GetTrackMatch(DFDCSegment *segment,vector<DFDCSegment*>package,
			     unsigned int &match_id);
  void LinkSegments(unsigned int pack1,vector<DFDCSegment *>packages[4]); 
  double Match(double p);

  bool GetTrackMatch(double q,DVector3 &pos,DVector3 &mom,
		     const DFDCSegment *segment);
  bool LinkStraySegment(const DFDCSegment *segment);

  bool DEBUG_HISTS,USE_FDC,APPLY_MOMENTUM_CORRECTION;
  double p_factor1,p_factor2;
  TH2F *match_dist_fdc,*match_center_dist2;
  TH1F *Hcircle_fit_prob;
  vector<double>z_wires;
  double endplate_z;
  double TARGET_Z;
  double MAX_R_VERTEX_LIMIT;
  double zpack[4];
  double FactorForSenseOfRotation;
  
  // Fit parameters
  double xc,yc,rc,z_vertex,q,phi0,tanl;
 
  // vector of fit results
  vector<DHelicalFit>fit_results;
 
  // Parameters at the end of the segment
  double xs,ys,zs;
  double p,cosphi,sinphi,twokappa,one_over_twokappa,cotl;
};

inline double DTrackCandidate_factory_FDCCathodes::Match(double p){
  double cut=-1.0+13.5/p;
  if (cut>16.0) cut=16.0;
  if (cut<1.) cut=1.0;
  return cut;
}


#endif // _DTrackCandidate_factory_FDCCathodes_

