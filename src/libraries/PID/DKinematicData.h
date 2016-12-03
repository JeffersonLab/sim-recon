#ifndef _DKINEMATICDATA_
#define _DKINEMATICDATA_

#include <JANA/JObject.h>
using namespace jana;

#include "GlueX.h"    

#include "DVector3.h"    
#include "DLorentzVector.h" 
#include "particleType.h" 
#include "TMatrixFSym.h"

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.9792458
#endif

class DKinematicData : public JObject
{
	public:

		// constructors and destructor
		DKinematicData(void);
		DKinematicData(const DVector3& aMomentum, const DVector3& aPosition, double aMass, double aCharge, const TMatrixFSym* aErrorMatrix = nullptr);
		virtual ~DKinematicData(void) {};

		//Reset
		virtual void Reset(void);

		//GETTERS
		Particle_t PID(void) const{return m_pid;}
		double mass(void) const{return m_mass;}
		double charge(void) const{return m_charge;}
		const DVector3& momentum(void) const{return m_momentum;}
		const DVector3& position(void) const{return m_position;}
		double time(void) const{return m_time;}

		//components
		double px(void) const{return m_momentum.Px();}
		double py(void) const{return m_momentum.Py();}
		double pz(void) const{return m_momentum.Pz();}
		double x(void) const{return m_position.X();}
		double y(void) const{return m_position.Y();}
		double z(void) const{return m_position.Z();}

		//derived quantities
		double energy(void) const{return sqrt(m_mass*m_mass + pmag2());}
		double pperp(void) const{return sqrt(px()*px() + py()*py());}
		double pperp2(void) const{return px()*px() + py()*py();}
		double pmag(void) const{return m_momentum.Mag();}
		double pmag2(void) const{return m_momentum.Mag2();}
		DLorentzVector lorentzMomentum(void) const{return DLorentzVector(m_momentum, energy());}

		//matrices, tracking
		const TMatrixFSym* errorMatrix(void) const{return m_errorMatrix;}
		const TMatrixFSym* TrackingErrorMatrix(void) const{return m_TrackingErrorMatrix;}
		bool forwardParmFlag(void) const{return m_use_forward_parameters;}
		void TrackingStateVector(double aVec[5]) const;

		//detector timing
		double t0(void) const{return m_t0;}
		double t0_err(void) const{return m_t0_err;}
		DetectorSystem_t t0_detector(void) const{return m_t0_detector;}
		double t1(void) const{return m_t1;}
		double t1_err(void) const{return m_t1_err;}
		DetectorSystem_t t1_detector(void) const{return m_t1_detector;}

		//pathlength, dE/dx
		double pathLength(void) const{return m_pathLength;}
		double pathLength_err(void) const{return m_pathLength_err;}
		double dEdx(void) const{return m_dedx;}

		//derived timing, PID info
		double TOF(void) const{return ((m_t1 - m_t0));}
		double TOF_err(void) const{return sqrt(m_t1_err*m_t1_err + m_t0_err*m_t0_err);}
		double deltaInvBeta(void) const{return (1.0/lorentzMomentum().Beta() - 1.0/measuredBeta());}
		double measuredInvBeta_err(void) const;
		double deltaBeta(void) const{return (lorentzMomentum().Beta() - measuredBeta());}
		double measuredBeta(void) const{return ((m_pathLength/(m_t1 - m_t0)))/29.9792458;}
		double measuredBeta_err(void) const;

		//SETTERS
		void setPID(Particle_t locPID){m_pid = locPID;}
		void setMass(double aMass){m_mass = aMass;}
		void setCharge(double aCharge){m_charge = aCharge;}
		void setMomentum(const DVector3& aMomentum){m_momentum = aMomentum;}
		void setPosition(const DVector3& aPosition){m_position = aPosition;}
		void setTime(double locTime){m_time = locTime;}

		// Time of flight information
		void setT0(double at0, double at0_err, DetectorSystem_t at0_detector);
		void setT1(double at1, double at1_err, DetectorSystem_t at1_detector);

		//pathlength, dE/dx
		void setPathLength(double apathLength, double apathLength_err);
		void setdEdx(double adedx){m_dedx = adedx;}

		//Tracking
		void setForwardParmFlag(bool aFlag){m_use_forward_parameters = aFlag;}
		void setErrorMatrix(const TMatrixFSym* aMatrix){m_errorMatrix = aMatrix;}
		void setTrackingErrorMatrix(const TMatrixFSym* aMatrix){m_TrackingErrorMatrix = aMatrix;}
		void setTrackingStateVector(double a1, double a2, double a3, double a4, double a5);

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "PID", "%i", (int)PID());
			AddString(items, "Name", "%s", ParticleType(PID()));
			AddString(items, "q", "%+1.0f", charge());
			AddString(items, "x(cm)", "%3.1f", x());
			AddString(items, "y(cm)", "%3.1f", y());
			AddString(items, "z(cm)", "%3.1f", z());
			AddString(items, "E(GeV)", "%2.4f", energy());
			AddString(items, "t(ns)", "%2.3f", time());
			AddString(items, "p(GeV/c)", "%2.3f", momentum().Mag());
			AddString(items, "theta(deg)", "%2.3f", momentum().Theta()*180.0/M_PI);
			AddString(items, "phi(deg)", "%2.3f", momentum().Phi()*180.0/M_PI);
		}

	private:

		// data members
		Particle_t m_pid;
		double m_mass;
		double m_charge;
		DVector3 m_momentum;
		DVector3 m_position;
		double m_time; // Time of the track propagated at m_position

		// detector timing
		double m_t0; /// Measured Start time (ns)
		double m_t0_err; /// Measured Start time error
		DetectorSystem_t m_t0_detector; /// Detector used to measure the start time
		double m_t1; /// Measured End of flight time (ns)
		double m_t1_err; /// Measured End of flight time error
		DetectorSystem_t m_t1_detector; /// Detector used to measure the end of flight time

		// pathlength, dE/dx
		double m_pathLength; /// Flight path length (cm)
		double m_pathLength_err; /// Flight path length err
		double m_dedx;

		// Tracking information //The setter's are responsible for managing the matrix memory!  These are NEVER the owners.
		const TMatrixFSym* m_errorMatrix;   // Order is (px, py, pz, x, y, z, t)
		const TMatrixFSym *m_TrackingErrorMatrix;  // order is q/pt,phi,tanl,D,z
		bool m_use_forward_parameters; // Flag indicating the use of the forward parameterization (x,y,tx,ty,q/p)
		double m_TrackingStateVector[5]; // order is q/pt,phi,tanl,D,z
};

/******************************************************************** CONSTRUCTORS *********************************************************************/

inline DKinematicData::DKinematicData(void) :
m_pid(Unknown), m_mass(0.0), m_charge(0.0), m_momentum(DVector3()), m_position(DVector3()), m_time(0.0),
m_t0(0.0), m_t0_err(0.0), m_t0_detector(SYS_NULL), m_t1(0.0), m_t1_err(0.0), m_t1_detector(SYS_NULL),
m_pathLength(0.0), m_pathLength_err(0.0), m_dedx(0.0), m_errorMatrix(nullptr), m_TrackingErrorMatrix(nullptr), m_use_forward_parameters(false)
{
	for(unsigned int i = 0; i < 5; ++i)
		m_TrackingStateVector[i] = 0.0;
}

inline DKinematicData::DKinematicData(const DVector3& aMomentum, const DVector3& aPosition, double aMass, double aCharge, const TMatrixFSym* aErrorMatrix):
m_pid(Unknown), m_mass(aMass), m_charge(aCharge), m_momentum(aMomentum), m_position(aPosition), m_time(0.0),
m_t0(0.0), m_t0_err(0.0), m_t0_detector(SYS_NULL), m_t1(0.0), m_t1_err(0.0), m_t1_detector(SYS_NULL),
m_pathLength(0.0), m_pathLength_err(0.0), m_dedx(0.0), m_errorMatrix(aErrorMatrix), m_TrackingErrorMatrix(nullptr), m_use_forward_parameters(false)
{
	for(unsigned int i = 0; i < 5; ++i)
		m_TrackingStateVector[i] = 0.0;
}

/********************************************************************** GETTERS ************************************************************************/

inline void DKinematicData::TrackingStateVector(double aVec[5]) const
{
	for (unsigned int i = 0; i < 5; ++i)
		aVec[i] = m_TrackingStateVector[i];
}

inline double DKinematicData::measuredBeta_err(void) const
{
	double a = pow(m_pathLength_err/m_pathLength, 2);
	double b = (pow(m_t1_err,2) + pow(m_t0_err,2))/pow(m_t1 - m_t0,2);

	double err = measuredBeta()*sqrt(a + b);
	return err;
}

inline double DKinematicData::measuredInvBeta_err(void) const
{
	double a = (m_t1_err*m_t1_err + m_t0_err*m_t0_err)/pow(m_t1 - m_t0,2);
	double b = pow(m_pathLength_err/m_pathLength,2);

	double err = (1.0/measuredBeta())*sqrt(a+b);
	return err;
}

/********************************************************************** SETTERS ************************************************************************/

inline void DKinematicData::setTrackingStateVector(double a1, double a2, double a3, double a4, double a5)
{
	m_TrackingStateVector[0]=a1;
	m_TrackingStateVector[1]=a2;
	m_TrackingStateVector[2]=a3;
	m_TrackingStateVector[3]=a4;
	m_TrackingStateVector[4]=a5;
}

inline void DKinematicData::setT0(double at0, double at0_err, const DetectorSystem_t at0_detector)
{
	m_t0 = at0;
	m_t0_err = at0_err;
	m_t0_detector = at0_detector;
}

inline void DKinematicData::setT1(double at1, double at1_err, const DetectorSystem_t at1_detector)
{
	m_t1 = at1;
	m_t1_err = at1_err;
	m_t1_detector = at1_detector;
}

inline void DKinematicData::setPathLength(double apathLength, double apathLength_err)
{
	m_pathLength = apathLength;
	m_pathLength_err = apathLength_err;
}

#endif /* _DKINEMATICDATA_ */
