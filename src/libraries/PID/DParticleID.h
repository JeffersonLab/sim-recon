// $Id$
//
//    File: DParticleID.h
// Created: Mon Feb 28 13:47:49 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_
#define _DParticleID_

#include <deque>
#include <limits>
#include <cmath>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <DVector3.h>
#include "HDGEOMETRY/DRootGeom.h"
#include <TRACKING/DTrackTimeBased_factory.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TOF/DTOFPoint.h>
#include <START_COUNTER/DSCHit.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackFinder.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <particleType.h>
#include <PID/DChargedTrackHypothesis.h>
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DEventRFBunch.h>
#include <PID/DDetectorMatches.h>
#include <TRACKING/DMagneticFieldStepper.h>

#include <TMath.h>

class DTrackTimeBased;
class DCDCTrackHit;

class DParticleID:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DParticleID);
 
  // Constructor and destructor
  DParticleID(JEventLoop *loop); // require JEventLoop in constructor
  virtual ~DParticleID(void){}

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

  void GetScintMPdEandSigma(double p,double M,double x,double &most_prob_dE,
			    double &sigma_dE) const;

	//called by track reconstruction
	bool MatchToTOF(const DReferenceTrajectory* rt, const vector<const DTOFPoint*>& locTOFPoints, double& locStartTime, double& locTimeVariance) const;
	bool MatchToBCAL(const DReferenceTrajectory* rt, const vector<const DBCALShower*>& locBCALShowers, double& locStartTime, double& locTimeVariance) const;
	bool MatchToFCAL(const DReferenceTrajectory* rt, const vector<const DFCALShower*>& locFCALShowers, double& locStartTime, double& locTimeVariance) const;
	bool MatchToSC(const DReferenceTrajectory* rt, const vector<const DSCHit*>& locSCHits, double& locStartTime, double& locTimeVariance) const;

	//matching tracks to hits/showers routines (can be called by DDetectorMatches factory)
	bool MatchToBCAL(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DBCALShower* locBCALShower, double locInputStartTime, DShowerMatchParams& locShowerMatchParams) const;
	bool MatchToTOF(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DTOFPoint* locTOFPoint, double locInputStartTime, DTOFHitMatchParams& locTOFHitMatchParams) const;
	bool MatchToFCAL(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DFCALShower* locFCALShower, double locInputStartTime, DShowerMatchParams& locShowerMatchParams) const;
	bool MatchToSC(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DSCHit* locSCHit, double locInputStartTime, DSCHitMatchParams& locSCHitMatchParams) const;

	// Alternate SC matching algorithm for straight line tracks
	bool MatchToSC(const DKinematicData *kd,
		       const vector<const DSCHit*>& locSCHits,
		       vector<DSCHitMatchParams>& locSCHitMatchParams) const;

	//select "best" matches //called by several factories
	bool Get_BestSCMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DSCHitMatchParams& locBestMatchParams) const;
	bool Get_BestBCALMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DShowerMatchParams& locBestMatchParams) const;
	bool Get_BestTOFMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DTOFHitMatchParams& locBestMatchParams) const;
	bool Get_BestFCALMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DShowerMatchParams& locBestMatchParams) const;

	//matching showers to tracks routines (called by DDetectorMatches factory)
	bool Distance_ToTrack(const DBCALShower* locBCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance, double& locDeltaPhi, double& locDeltaZ) const;
	bool Distance_ToTrack(const DFCALShower* locFCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const;
	bool Distance_ToTrack(const DTOFPoint* locTOFPoint, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const;
	bool Distance_ToTrack(const DSCHit* locSCHit, const DReferenceTrajectory* rt, double locInputStartTime, double& locDeltaPhi) const;

	double Calc_BCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
	double Calc_FCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
	double Calc_TOFFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
	double Calc_SCFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches) const;

  virtual Particle_t IDTrack(float locCharge, float locMass) const;

  void Calc_PropagatedRFTime(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, double& locPropagatedRFTime, bool locRFTimeFixedFlag) const;
  bool Calc_TrackStartTime(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, double& locStartTime, double& locStartTimeVariance, bool locRFTimeFixedFlag) const;
  double Calc_TimingChiSq(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag, unsigned int &locNDF, double& locTimingPull) const;
  double Calc_TimingChiSq(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DEventRFBunch* locEventRFBunch, unsigned int &locNDF, double& locPull) const;
  void Calc_ChargedPIDFOM(DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag) const;

  protected:
		// gas material properties
		double dKRhoZoverA_FDC, dRhoZoverA_FDC, dLnI_FDC;	
		double dKRhoZoverA_Scint, dRhoZoverA_Scint, dLnI_Scint;
		double dKRhoZoverA_CDC, dRhoZoverA_CDC, dLnI_CDC;
		double dDensity_FDC;
		double dDensity_CDC;
		double dA_CDC;
		double dA_FDC;

		double BCAL_Z_CUT,BCAL_PHI_CUT_PAR1,BCAL_PHI_CUT_PAR2;
		double FCAL_CUT_PAR1,FCAL_CUT_PAR2;

		double DELTA_R_FCAL;
		double C_EFFECTIVE; // start counter light propagation speed
		double ATTEN_LENGTH; // Start counter attenuation length
		double OUT_OF_TIME_CUT; //for all matches

 private: 
 
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID();
  
	// start counter geometry parameters
	double sc_leg_tcor;
	double sc_angle_cor;
	vector<DVector3> sc_pos;
	vector<DVector3> sc_norm;
	double dSCdphi;
	double dSCphi0;

  double dTargetZCenter;
  double dRFBunchFrequency;

  const DTrackFinder *finder;
};

#endif // _DParticleID_

