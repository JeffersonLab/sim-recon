// $Id$
//
//    File: DParticleID.h
// Created: Mon Feb 28 13:47:49 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_
#define _DParticleID_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <DVector3.h>
#include "HDGEOMETRY/DRootGeom.h"
#include <TRACKING/DTrackTimeBased_factory.h>
#include <BCAL/DBCALShower.h>
//#include <FCAL/DFCALCluster.h>
#include <FCAL/DFCALShower.h>
#include <TOF/DTOFPoint.h>
#include <START_COUNTER/DSCHit.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <particleType.h>
#include <PID/DChargedTrackHypothesis.h>
#include <TRACKING/DMagneticFieldStepper.h>

class DTrackTimeBased;
class DCDCTrackHit;

class DParticleID:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DParticleID);
 
  // Constructor and destructor
  DParticleID(JEventLoop *loop); // require JEventLoop in constructor
  virtual ~DParticleID();

  class dedx_t{
  public:
    dedx_t(double dE,double dx, double p):dE(dE),dx(dx),p(p){dEdx = dE/dx;}
      double dE; // energy loss in layer
      double dx; // path length in layer
	  double dEdx; // ratio dE/dx
      double p;  // momentum at this dE/dx measurement
      
  };

  virtual jerror_t CalcDCdEdxChiSq(DChargedTrackHypothesis *locChargedTrackHypothesis) const = 0;

  jerror_t GetDCdEdxHits(const DTrackTimeBased *track, vector<dedx_t>& dEdxHits_CDC, vector<dedx_t>& dEdxHits_FDC) const;
  jerror_t CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const;
  jerror_t CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, const vector<dedx_t>& locdEdxHits_CDC, const vector<dedx_t>& locdEdxHits_FDC, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const;

  jerror_t CalcdEdxHit(const DVector3 &mom, const DVector3 &pos, const DCDCTrackHit *hit, pair <double,double> &dedx) const;
  jerror_t GroupTracks(vector<const DTrackTimeBased *> &tracks, vector<vector<const DTrackTimeBased*> >&grouped_tracks) const;

  jerror_t MatchToTOF(const DReferenceTrajectory *rt, DTrackFitter::fit_type_t fit_type, vector<const DTOFPoint*>&tof_points, double &tproj, unsigned int &tof_match_id, double &locPathLength, double &locFlightTime) const;
  jerror_t MatchToBCAL(const DReferenceTrajectory *rt, const vector<const DBCALShower*>& locInputBCALShowers, vector<const DBCALShower*>& locMatchedBCALShowers, double& locProjectedTime, double& locPathLength, double& locFlightTime) const;
  jerror_t MatchToFCAL(const DReferenceTrajectory *rt, const vector<const DFCALShower*>& locInputFCALShowers, vector<const DFCALShower*>& locMatchedFCALShowers, double& locProjectedTime, double& locPathLength, double& locFlightTime) const;
  jerror_t MatchToSC(const DReferenceTrajectory *rt, DTrackFitter::fit_type_t fit_type, vector<const DSCHit*>&sc_hits, double &tproj,unsigned int &sc_match_id, double &locPathLength, double &locFlightTime) const;
  jerror_t MatchToSC(const DKinematicData &parms, vector<const DSCHit*>&sc_hits, double &tproj,unsigned int &sc_match_id) const;

  virtual Particle_t IDTrack(float locCharge, float locMass) const;
  void Calc_TimingChiSq(DChargedTrackHypothesis* locChargedTrackHypothesis, double locRFTime, double locRFBunchFrequency) const;
  void Calc_ChargedPIDFOM(DChargedTrackHypothesis* locChargedTrackHypothesis, double locRFTime, double locRFBunchFrequency) const;

  protected:
		// gas material properties
		double dKRhoZoverA_FDC, dRhoZoverA_FDC, dLnI_FDC;
		double dKRhoZoverA_CDC, dRhoZoverA_CDC, dLnI_CDC;
		double dDensity_FDC;
		double dDensity_CDC;
		double dA_CDC;
		double dA_FDC;

 private: 
 
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID();
  
  const DGeometry *geom;
  const DMagneticFieldMap *bfield;
  DMagneticFieldStepper *stepper;
  double dTargetZCenter;


  // start counter geometry parameters
  double sc_leg_tcor;
  double sc_angle_cor;
  vector<DVector3>sc_pos;
  vector<DVector3>sc_norm;  

  double DELTA_R_BCAL;
  double DELTA_R_FCAL;
};

#endif // _DParticleID_

