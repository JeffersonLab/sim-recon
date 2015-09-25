// $Id$
//
//    File: DNeutralParticleHypothesis_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TMath.h>

#include "DNeutralParticleHypothesis_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	dTargetLength = 30.0;
	double locTargetCenterZ = 65.0;
	dTargetRadius = 1.5; //FIX: grab from database!!!

	locEventLoop->GetSingle(dParticleID);

	if(locGeometry)
	{
		locGeometry->GetTargetZ(locTargetCenterZ);
		locGeometry->GetTargetLength(dTargetLength);
	}
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<Particle_t> locPIDHypotheses;
	locPIDHypotheses.push_back(Gamma);
	locPIDHypotheses.push_back(Neutron);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	// Loop over DNeutralShowers
	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
	{
		const DNeutralShower *locNeutralShower = locNeutralShowers[loc_i];
		// Loop over vertices and PID hypotheses & create DNeutralParticleHypotheses for each combination
		for(size_t loc_j = 0; loc_j < locPIDHypotheses.size(); ++loc_j)
		{
			DNeutralParticleHypothesis* locNeutralParticleHypothesis = Create_DNeutralParticleHypothesis(locNeutralShower, locPIDHypotheses[loc_j], locEventRFBunch, locVertex);
			if(locNeutralParticleHypothesis != NULL)
				_data.push_back(locNeutralParticleHypothesis);	
		}
	}

	return NOERROR;
}

DNeutralParticleHypothesis* DNeutralParticleHypothesis_factory::Create_DNeutralParticleHypothesis(const DNeutralShower* locNeutralShower, Particle_t locPID, const DEventRFBunch* locEventRFBunch, const DVertex* locVertex) const
{
	DVector3 locVertexGuess = locVertex->dSpacetimeVertex.Vect();
	double locStartTime = locVertex->dSpacetimeVertex.T();

	double locHitTime = locNeutralShower->dSpacetimeVertex.T();
	DVector3 locHitPoint = locNeutralShower->dSpacetimeVertex.Vect();

	// Calculate DNeutralParticleHypothesis Quantities (projected time at vertex for given id, etc.)
	double locMass = ParticleMass(locPID);

	DVector3 locPath = locHitPoint - locVertexGuess;
	double locPathLength = locPath.Mag();
	if(!(locPathLength > 0.0))
		return NULL; //invalid, will divide by zero when creating error matrix, so skip!

	DVector3 locMomentum(locPath);
	DMatrixDSym locParticleCovariance;

	double locProjectedTime = 0.0, locPMag = 0.0;
	if(locPID != Gamma)
	{
		double locDeltaT = locHitTime - locStartTime;
		double locBeta = locPathLength/(locDeltaT*29.9792458);
		if(locBeta >= 1.0)
			locBeta = 0.9999;
		if(locBeta < 0.0)
			locBeta = 0.0;
		double locGamma = 1.0/sqrt(1.0 - locBeta*locBeta);
		locPMag = locGamma*locBeta*locMass;
		locMomentum.SetMag(locPMag);
		locProjectedTime = locStartTime + (locVertexGuess.Z() - dTargetCenter.Z())/29.9792458;
		Calc_ParticleCovariance_Massive(locNeutralShower, locMass, locDeltaT, locEventRFBunch->dTimeVariance, locMomentum, locPath, locParticleCovariance);
	}
	else
	{
		locPMag = locNeutralShower->dEnergy;
		double locFlightTime = locPathLength/29.9792458;
		locProjectedTime = locHitTime - locFlightTime;
		locMomentum.SetMag(locPMag);
		Calc_ParticleCovariance_Photon(locNeutralShower, locMomentum, locPath, locParticleCovariance);
	}

	// Build DNeutralParticleHypothesis // dEdx not set
	DNeutralParticleHypothesis* locNeutralParticleHypothesis = new DNeutralParticleHypothesis;
	locNeutralParticleHypothesis->AddAssociatedObject(locNeutralShower);

	locNeutralParticleHypothesis->dNeutralShowerID = locNeutralShower->dShowerID;
	locNeutralParticleHypothesis->setPID(locPID);
	locNeutralParticleHypothesis->setMass(locMass);
	locNeutralParticleHypothesis->setCharge(0.0);
	locNeutralParticleHypothesis->setMomentum(locMomentum);
	locNeutralParticleHypothesis->setPosition(locVertexGuess);
	locNeutralParticleHypothesis->setT0(locStartTime, sqrt(locEventRFBunch->dTimeVariance), locEventRFBunch->dTimeSource);
	locNeutralParticleHypothesis->setTime(locProjectedTime);
	locNeutralParticleHypothesis->setT1(locNeutralShower->dSpacetimeVertex.T(), sqrt(locNeutralShower->dCovarianceMatrix(4, 4)), locNeutralShower->dDetectorSystem);
	locNeutralParticleHypothesis->setPathLength(locPathLength, 0.0); //zero uncertainty (for now)
	locNeutralParticleHypothesis->setErrorMatrix(locParticleCovariance);

	// Calculate DNeutralParticleHypothesis FOM
	unsigned int locNDF = 0;
	double locChiSq = 0.0;
	double locFOM = -1.0; //undefined for non-photons
	if(locPID == Gamma)
	{
		double locTimePull = 0.0;
		locChiSq = dParticleID->Calc_TimingChiSq(locNeutralParticleHypothesis, locNDF, locTimePull);
		locFOM = TMath::Prob(locChiSq, locNDF);
	}
	locNeutralParticleHypothesis->dChiSq = locChiSq;
	locNeutralParticleHypothesis->dNDF = locNDF;
	locNeutralParticleHypothesis->dFOM = locFOM;
	return locNeutralParticleHypothesis;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory::fini(void)
{
	return NOERROR;
}

void DNeutralParticleHypothesis_factory::Calc_ParticleCovariance_Photon(const DNeutralShower* locNeutralShower, const DVector3& locMomentum, const DVector3& locPathVector, DMatrixDSym& locParticleCovariance) const
{
	//build 8x8 matrix: 5x5 shower, 3x3 vertex position
	DMatrixDSym locShowerPlusVertCovariance(8);
	for(unsigned int loc_l = 0; loc_l < 5; ++loc_l) //shower: e, x, y, z, t
	{
		for(unsigned int loc_m = 0; loc_m < 5; ++loc_m)
			locShowerPlusVertCovariance(loc_l, loc_m) = locNeutralShower->dCovarianceMatrix(loc_l, loc_m);
	}
	locShowerPlusVertCovariance(5, 5) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position x
	locShowerPlusVertCovariance(6, 6) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position y
	locShowerPlusVertCovariance(7, 7) = dTargetLength*dTargetLength/12.0; //vertex position z

	DVector3 locDeltaX = -1.0*locPathVector; //defined oppositely in document! //delta_x defined here as common_vertex - hit_point
	DVector3 locDeltaXOverDeltaXSq = (1.0/locDeltaX.Mag2())*locDeltaX;
	DVector3 locUnitP = (1.0/locNeutralShower->dEnergy)*locMomentum;
	DVector3 locUnitDeltaXOverC = (1.0/(29.9792458*locDeltaX.Mag()))*locDeltaX;

	//build transform matrix
	DMatrix locTransformMatrix(7, 8);

	locTransformMatrix(0, 0) = locUnitP.X(); //partial deriv of px wrst shower-e
	locTransformMatrix(0, 1) = locMomentum.Px()*(locDeltaXOverDeltaXSq.X() - 1.0/locDeltaX.X()); //partial deriv of px wrst shower-x
	locTransformMatrix(0, 2) = locMomentum.Px()*locDeltaX.Y()/locDeltaX.Mag2(); //partial deriv of px wrst shower-y
	locTransformMatrix(0, 3) = locMomentum.Px()*locDeltaX.Z()/locDeltaX.Mag2(); //partial deriv of px wrst shower-z
	locTransformMatrix(0, 5) = -1.0*locTransformMatrix(0, 1); //partial deriv of px wrst vert-x
	locTransformMatrix(0, 6) = -1.0*locTransformMatrix(0, 2); //partial deriv of px wrst vert-y
	locTransformMatrix(0, 7) = -1.0*locTransformMatrix(0, 3); //partial deriv of px wrst vert-z

	locTransformMatrix(1, 0) = locUnitP.Y(); //partial deriv of py wrst shower-e
	locTransformMatrix(1, 1) = locMomentum.Py()*locDeltaX.X()/locDeltaX.Mag2(); //partial deriv of py wrst shower-x
	locTransformMatrix(1, 2) = locMomentum.Py()*(locDeltaXOverDeltaXSq.Y() - 1.0/locDeltaX.Y()); //partial deriv of py wrst shower-y
	locTransformMatrix(1, 3) = locMomentum.Py()*locDeltaX.Z()/locDeltaX.Mag2(); //partial deriv of py wrst shower-z
	locTransformMatrix(1, 5) = -1.0*locTransformMatrix(1, 1); //partial deriv of py wrst vert-x
	locTransformMatrix(1, 6) = -1.0*locTransformMatrix(1, 2); //partial deriv of py wrst vert-y
	locTransformMatrix(1, 7) = -1.0*locTransformMatrix(1, 3); //partial deriv of py wrst vert-z

	locTransformMatrix(2, 0) = locUnitP.Z(); //partial deriv of pz wrst shower-e
	locTransformMatrix(2, 1) = locMomentum.Pz()*locDeltaX.X()/locDeltaX.Mag2(); //partial deriv of pz wrst shower-x
	locTransformMatrix(2, 2) = locMomentum.Pz()*locDeltaX.Y()/locDeltaX.Mag2(); //partial deriv of pz wrst shower-y
	locTransformMatrix(2, 3) = locMomentum.Pz()*(locDeltaXOverDeltaXSq.Z() - 1.0/locDeltaX.Z()); //partial deriv of pz wrst shower-z
	locTransformMatrix(2, 5) = -1.0*locTransformMatrix(2, 1); //partial deriv of pz wrst vert-x
	locTransformMatrix(2, 6) = -1.0*locTransformMatrix(2, 2); //partial deriv of pz wrst vert-y
	locTransformMatrix(2, 7) = -1.0*locTransformMatrix(2, 3); //partial deriv of pz wrst vert-z

	locTransformMatrix(3, 5) = 1.0; //partial deriv of x wrst vertex-x
	locTransformMatrix(4, 6) = 1.0; //partial deriv of y wrst vertex-y
	locTransformMatrix(5, 7) = 1.0; //partial deriv of z wrst vertex-z

	locTransformMatrix(6, 0) = 0.0; //partial deriv of t wrst shower-e //beta defined = 1
	locTransformMatrix(6, 1) = locUnitDeltaXOverC.X(); //partial deriv of t wrst shower-x
	locTransformMatrix(6, 2) = locUnitDeltaXOverC.Y(); //partial deriv of t wrst shower-y
	locTransformMatrix(6, 3) = locUnitDeltaXOverC.Z(); //partial deriv of t wrst shower-z
	locTransformMatrix(6, 4) = 1.0; //partial deriv of t wrst shower-t
	locTransformMatrix(6, 5) = -1.0*locTransformMatrix(6, 1); //partial deriv of t wrst vert-x
	locTransformMatrix(6, 6) = -1.0*locTransformMatrix(6, 2); //partial deriv of t wrst vert-y
	locTransformMatrix(6, 7) = -1.0*locTransformMatrix(6, 3); //partial deriv of t wrst vert-z

	//convert
	locParticleCovariance.ResizeTo(7, 7);
	locParticleCovariance = locShowerPlusVertCovariance.Similarity(locTransformMatrix);
}

void DNeutralParticleHypothesis_factory::Calc_ParticleCovariance_Massive(const DNeutralShower* locNeutralShower, double locMass, double locDeltaT, double locStartTimeVariance, const DVector3& locMomentum, const DVector3& locPathVector, DMatrixDSym& locParticleCovariance) const
{
	//build 9x9 matrix: 5x5 shower, 4x4 vertex position & time
	DMatrixDSym locShowerPlusVertCovariance(9);
	for(unsigned int loc_l = 0; loc_l < 5; ++loc_l) //shower: e, x, y, z, t
	{
		for(unsigned int loc_m = 0; loc_m < 5; ++loc_m)
			locShowerPlusVertCovariance(loc_l, loc_m) = locNeutralShower->dCovarianceMatrix(loc_l, loc_m);
	}
	locShowerPlusVertCovariance(5, 5) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position x
	locShowerPlusVertCovariance(6, 6) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position y
	locShowerPlusVertCovariance(7, 7) = dTargetLength*dTargetLength/12.0; //vertex position z
	locShowerPlusVertCovariance(8, 8) = locStartTimeVariance; //vertex time

	DVector3 locDeltaX = -1.0*locPathVector; //defined oppositely in document! //delta_x defined here as common_vertex - hit_point
	locDeltaT *= -1.0; //defined oppositely in document! //delta_t defined here as common_t - hit_t
	double locCSq = 29.9792458*29.9792458;
	double locCDeltaTSq = locCSq*locDeltaT*locDeltaT;
	double locDeltaX4Sq = locCDeltaTSq - locDeltaX.Mag2();
	DVector3 locDeltaXOverDeltaX4Sq = (1.0/locDeltaX4Sq)*locDeltaX;
	DVector3 locEPVecOverPSq = (locNeutralShower->dEnergy/locMomentum.Mag2())*locMomentum;
	DVector3 locEPVecOverCPMagDeltaXMag = (locNeutralShower->dEnergy/(29.9792458*locDeltaX.Mag()*locMomentum.Mag()))*locDeltaX;

	//build transform matrix
	DMatrix locTransformMatrix(7, 9);

	locTransformMatrix(0, 0) = 0.0; //partial deriv of px wrst shower-e
	locTransformMatrix(0, 1) = -1.0*locMomentum.Px()*(locDeltaXOverDeltaX4Sq.X() + 1.0/locDeltaX.X()); //partial deriv of px wrst shower-x
	locTransformMatrix(0, 2) = -1.0*locMomentum.Px()*locDeltaX.Y()/locDeltaX4Sq; //partial deriv of px wrst shower-y
	locTransformMatrix(0, 3) = -1.0*locMomentum.Px()*locDeltaX.Z()/locDeltaX4Sq; //partial deriv of px wrst shower-z
	locTransformMatrix(0, 4) = locDeltaT*locCSq*locMomentum.Px()/locDeltaX4Sq; //partial deriv of px wrst shower-t
	locTransformMatrix(0, 5) = -1.0*locTransformMatrix(0, 1); //partial deriv of px wrst vert-x
	locTransformMatrix(0, 6) = -1.0*locTransformMatrix(0, 2); //partial deriv of px wrst vert-y
	locTransformMatrix(0, 7) = -1.0*locTransformMatrix(0, 3); //partial deriv of px wrst vert-z
	locTransformMatrix(0, 8) = -1.0*locTransformMatrix(0, 4); //partial deriv of px wrst vert-t

	locTransformMatrix(1, 0) = 0.0; //partial deriv of py wrst shower-e
	locTransformMatrix(1, 1) = -1.0*locMomentum.Py()*locDeltaX.X()/locDeltaX4Sq; //partial deriv of py wrst shower-x
	locTransformMatrix(1, 2) = -1.0*locMomentum.Py()*(locDeltaXOverDeltaX4Sq.Y() + 1.0/locDeltaX.Y()); //partial deriv of py wrst shower-y
	locTransformMatrix(1, 3) = -1.0*locMomentum.Py()*locDeltaX.Z()/locDeltaX4Sq; //partial deriv of py wrst shower-z
	locTransformMatrix(1, 4) = locDeltaT*locCSq*locMomentum.Py()/locDeltaX4Sq; //partial deriv of py wrst shower-t
	locTransformMatrix(1, 5) = -1.0*locTransformMatrix(1, 1); //partial deriv of py wrst vert-x
	locTransformMatrix(1, 6) = -1.0*locTransformMatrix(1, 2); //partial deriv of py wrst vert-y
	locTransformMatrix(1, 7) = -1.0*locTransformMatrix(1, 3); //partial deriv of py wrst vert-z
	locTransformMatrix(1, 8) = -1.0*locTransformMatrix(1, 4); //partial deriv of py wrst vert-t

	locTransformMatrix(2, 0) = 0.0; //partial deriv of pz wrst shower-e
	locTransformMatrix(2, 1) = -1.0*locMomentum.Pz()*locDeltaX.X()/locDeltaX4Sq; //partial deriv of pz wrst shower-x
	locTransformMatrix(2, 2) = -1.0*locMomentum.Pz()*locDeltaX.Y()/locDeltaX4Sq; //partial deriv of pz wrst shower-y
	locTransformMatrix(2, 3) = -1.0*locMomentum.Pz()*(locDeltaXOverDeltaX4Sq.Z() + 1.0/locDeltaX.Z()); //partial deriv of pz wrst shower-z
	locTransformMatrix(2, 4) = locDeltaT*locCSq*locMomentum.Pz()/locDeltaX4Sq; //partial deriv of pz wrst shower-t
	locTransformMatrix(2, 5) = -1.0*locTransformMatrix(2, 1); //partial deriv of pz wrst vert-x
	locTransformMatrix(2, 6) = -1.0*locTransformMatrix(2, 2); //partial deriv of pz wrst vert-y
	locTransformMatrix(2, 7) = -1.0*locTransformMatrix(2, 3); //partial deriv of pz wrst vert-z
	locTransformMatrix(2, 8) = -1.0*locTransformMatrix(2, 4); //partial deriv of pz wrst vert-t

	locTransformMatrix(3, 5) = 1.0; //partial deriv of x wrst vertex-x
	locTransformMatrix(4, 6) = 1.0; //partial deriv of y wrst vertex-y
	locTransformMatrix(5, 7) = 1.0; //partial deriv of z wrst vertex-z
	locTransformMatrix(6, 8) = 1.0; //partial deriv of t wrst vertex-t

	//convert
	locParticleCovariance.ResizeTo(7, 7);
	locParticleCovariance = locShowerPlusVertCovariance.Similarity(locTransformMatrix);
}

