// $Id$
//
//    File: DParticleID.h
// Created: Mon Feb 28 13:47:49 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_
#define _DParticleID_

#include <deque>
#include <map>
#include <limits>
#include <cmath>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <DVector3.h>
#include "HDGEOMETRY/DRootGeom.h"
#include <TRACKING/DTrackTimeBased_factory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <BCAL/DBCALShower.h>
#include <BCAL/DBCALCluster.h>
#include <FCAL/DFCALShower.h>
#include <FCAL/DFCALGeometry_factory.h>
#include <TOF/DTOFPoint.h>
#include <TOF/DTOFPaddleHit.h>
#include <TOF/DTOFGeometry_factory.h>
#include <TOF/DTOFPoint_factory.h>
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
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackCandidate.h>

#include <TMath.h>

class DTrackTimeBased;

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

  void GetScintMPdEandSigma(double p,double M,double x,double &most_prob_dE, double &sigma_dE) const;
  double GetMostProbabledEdx_DC(double p,double mass,double dx, bool locIsCDCFlag) const; //bool is false for FDC
  double GetdEdxSigma_DC(double num_hits,double p,double mass, double mean_path_length, bool locIsCDCFlag) const; //bool is false for FDC

	//called by track reconstruction
	bool MatchToTOF(const DReferenceTrajectory* rt, const vector<const DTOFPoint*>& locTOFPoints, double& locStartTime, double& locTimeVariance) const;
	bool MatchToBCAL(const DReferenceTrajectory* rt, const vector<const DBCALShower*>& locBCALShowers, double& locStartTime, double& locTimeVariance) const;
	bool MatchToFCAL(const DReferenceTrajectory* rt, const vector<const DFCALShower*>& locFCALShowers, double& locStartTime, double& locTimeVariance) const;
	bool MatchToSC(const DReferenceTrajectory* rt, const vector<const DSCHit*>& locSCHits, double& locStartTime, double& locTimeVariance, bool locIsTimeBased=false) const;

	// Routines to predict which detector elements should fire given a track
		//SC locProjBarrelRegion = true/false for barrel / (bend / nose) regions
	unsigned int PredictSCSector(const DReferenceTrajectory* rt, const double dphi_cut, DVector3* locProjPos = NULL, bool* locProjBarrelRegion = NULL, double* locMinDPhi = NULL) const;
	bool PredictTOFPaddles(const DReferenceTrajectory *rt,
			       unsigned int &hbar,unsigned int &vbar,
			       DVector3 *intersection=NULL) const;
	bool PredictFCALHit(const DReferenceTrajectory *rt,
			    unsigned int &row, unsigned int &col,
			    DVector3 *intersection=NULL) const;
	bool PredictBCALWedge(const DReferenceTrajectory *rt,
			      unsigned int &module,unsigned int &sector,
			      DVector3 *intersection=NULL) const;

	//matching tracks to hits/showers routines (can be called by DDetectorMatches factory)
	bool MatchToBCAL(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DBCALShower* locBCALShower, double locInputStartTime, DBCALShowerMatchParams& locShowerMatchParams) const;
	bool MatchToTOF(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DTOFPoint* locTOFPoint, double locInputStartTime, DTOFHitMatchParams& locTOFHitMatchParams) const;
	bool MatchToFCAL(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DFCALShower* locFCALShower, double locInputStartTime, DFCALShowerMatchParams& locShowerMatchParams) const;
	bool MatchToSC(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DSCHit* locSCHit, double locInputStartTime, DSCHitMatchParams& locSCHitMatchParams, bool locIsTimeBased=false, const float* locInputDeltaPhiCut=NULL, DVector3 *IntersectionPoint=NULL, DVector3 *IntersectionDir=NULL) const;

	//select "best" matches //called by several factories
	bool Get_BestSCMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DSCHitMatchParams& locBestMatchParams) const;
	void Get_BestSCMatchParams(vector<DSCHitMatchParams>& locSCHitMatchParams, DSCHitMatchParams& locBestMatchParams) const;
	bool Get_BestBCALMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DBCALShowerMatchParams& locBestMatchParams) const;
	void Get_BestBCALMatchParams(DVector3 locMomentum, vector<DBCALShowerMatchParams>& locShowerMatchParams, DBCALShowerMatchParams& locBestMatchParams) const;
	bool Get_BestTOFMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DTOFHitMatchParams& locBestMatchParams) const;
	void Get_BestTOFMatchParams(vector<DTOFHitMatchParams>& locTOFHitMatchParams, DTOFHitMatchParams& locBestMatchParams) const;
	bool Get_BestFCALMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DFCALShowerMatchParams& locBestMatchParams) const;
	void Get_BestFCALMatchParams(vector<DFCALShowerMatchParams>& locShowerMatchParams, DFCALShowerMatchParams& locBestMatchParams) const;

	//matching showers to tracks routines (mostly called by DDetectorMatches factory)
	bool Distance_ToTrack(const DBCALShower* locBCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance, double& locDeltaPhi, double& locDeltaZ) const;
	bool Distance_ToTrack(const DFCALShower* locFCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const;
	bool Distance_ToTrack(const DTOFPoint* locTOFPoint, const DReferenceTrajectory* rt, double locInputStartTime, double& locDeltaX, double& locDeltaY) const;
	bool Distance_ToTrack(const DSCHit* locSCHit, const DReferenceTrajectory* rt, double locInputStartTime, double& locDeltaPhi) const;

	//select closest shower/hit to track //track MUST be either DTrackTimeBased or DTrackWireBased //NULL if no possibly-valid matches
	const DBCALShower* Get_ClosestToTrack_BCAL(const DKinematicData* locTrack, vector<const DBCALShower*>& locBCALShowers, double& locBestMatchDeltaPhi, double& locBestMatchDeltaZ) const;
	const DFCALShower* Get_ClosestToTrack_FCAL(const DKinematicData* locTrack, vector<const DFCALShower*>& locFCALShowers, double& locBestDistance) const;
	const DSCHit* Get_ClosestToTrack_SC(const DKinematicData* locTrack, vector<const DSCHit*>& locSCHits, double& locBestDeltaPhi) const;
	const DTOFPoint* Get_ClosestToTrack_TOFPoint(const DKinematicData* locTrack, vector<const DTOFPoint*>& locTOFPoints, double& locBestDeltaX, double& locBestDeltaY) const;

	//first in pair is vertical, second is horizontal // NULL If none / doesn't hit TOF
	pair<const DTOFPaddleHit*, const DTOFPaddleHit*> Get_ClosestToTrack_TOFPaddles(const DKinematicData* locTrack, vector<const DTOFPaddleHit*>& locTOFPaddleHits, double& locBestDeltaX, double& locBestDeltaY, double locMaxDeltaT = -1.0) const;
	const DTOFPaddleHit* Get_ClosestTOFPaddleHit_Horizontal(const DReferenceTrajectory* locReferenceTrajectory, const vector<const DTOFPaddleHit*>& locTOFPaddleHits, double locInputStartTime, double& locBestDeltaY) const;
	const DTOFPaddleHit* Get_ClosestTOFPaddleHit_Vertical(const DReferenceTrajectory* locReferenceTrajectory, const vector<const DTOFPaddleHit*>& locTOFPaddleHits, double locInputStartTime, double& locBestDeltaX) const;

	double Calc_BCALFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const;
	double Calc_FCALFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const;
	double Calc_TOFFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const;
	double Calc_SCFlightTimePCorrelation(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches) const;

	virtual Particle_t IDTrack(float locCharge, float locMass) const;

	double Calc_PropagatedRFTime(const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch) const;
	double Calc_TimingChiSq(const DKinematicData* locKinematicData, unsigned int &locNDF, double& locTimingPull) const;
	void Calc_ChargedPIDFOM(DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch) const;

	unsigned int Get_CDCRingBitPattern(vector<const DCDCTrackHit*>& locCDCTrackHits) const;
	unsigned int Get_FDCPlaneBitPattern(vector<const DFDCPseudo*>& locFDCPseudos) const;
	void Get_CDCRings(int locBitPattern, set<int>& locCDCRings) const;
	void Get_FDCPlanes(int locBitPattern, set<int>& locFDCPlanes) const;

	void Get_CDCNumHitRingsPerSuperlayer(int locBitPattern, map<int, int>& locNumHitRingsPerSuperlayer) const;
	void Get_CDCNumHitRingsPerSuperlayer(const set<int>& locCDCRings, map<int, int>& locNumHitRingsPerSuperlayer) const;
	void Get_FDCNumHitPlanesPerPackage(int locBitPattern, map<int, int>& locNumHitPlanesPerPackage) const;
	void Get_FDCNumHitPlanesPerPackage(const set<int>& locFDCPlanes, map<int, int>& locNumHitPlanesPerPackage) const;

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
		double BCAL_PHI_CUT_PAR3;
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
	vector<vector<DVector3> >sc_dir; // direction vector in plane of plastic
	vector<vector<DVector3> >sc_pos;
	vector<vector<DVector3> >sc_norm;
	double dSCdphi;
	double dSCphi0;
	// start counter calibration parameters
	// Propagation time (pt) parameters
	enum sc_region_t{
	  SC_STRAIGHT,
	  SC_BEND,
	  SC_NOSE,
	};
	vector<double>sc_veff[3];
	vector<double>sc_pt_yint[3];
	vector<double>sc_pt_slope[3];
	// Attenuation (attn) calibration paramerters
	enum sc_region_attn{
	  SC_STRAIGHT_ATTN,
	  SC_BENDNOSE_ATTN,
	};
	vector<double> sc_attn_A[2];
	vector<double> sc_attn_B[2];
	vector<double> sc_attn_C[2];

    vector<double> sc_paddle_resols;

	// FCAL geometry
	double dFCALz;
	const DFCALGeometry *dFCALGeometry;

	// TOF calibration constants
		// used to update hit energy & time when matching to un-matched, position-ill-defined bars
	const DTOFGeometry* dTOFGeometry;
	vector<double> propagation_speed;
	double TOF_HALFPADDLE;
	double dHalfPaddle_OneSided;
	double TOF_ATTEN_LENGTH;
	double TOF_E_THRESHOLD;
	double ONESIDED_PADDLE_MIDPOINT_MAG; //+/- this number for North/South

  double dTargetZCenter;
  double SC_DPHI_CUT,SC_DPHI_CUT_WB,SC_DPHI_CUT_SLOPE;

  const DTrackFinder *finder;
  DTOFPoint_factory* dTOFPointFactory;
};

#endif // _DParticleID_

