#ifndef _DKINEMATICDATA_
#define _DKINEMATICDATA_

#include <memory>

#include <JANA/JObject.h>

#include "GlueX.h"    

#include "DVector3.h"    
#include "DLorentzVector.h" 
#include "particleType.h" 
#include "TMatrixFSym.h"

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.9792458
#endif

using namespace std;
using namespace jana;

class DKinematicData : public JObject
{
	public:

		// constructors and destructor
		DKinematicData(void);
		DKinematicData(const DKinematicData& locSourceData, bool locShareNonKinematicsFlag = false, bool locShareKinematicsFlag = false);
		virtual ~DKinematicData(void) {};

		//Assignment operator
		DKinematicData& operator=(const DKinematicData& locSourceData);

		//Reset
		virtual void Reset(void);

		//GETTERS
		Particle_t PID(void) const{return dNonKinematics->m_pid;}
		double mass(void) const{return ParticleMass(dNonKinematics->m_pid);}
		double charge(void) const{return ParticleCharge(dNonKinematics->m_pid);}
		const DVector3& momentum(void) const{return dKinematics->m_momentum;}
		const DVector3& position(void) const{return dKinematics->m_position;}
		double time(void) const{return dKinematics->m_time;}

		//components
		double px(void) const{return dKinematics->m_momentum.Px();}
		double py(void) const{return dKinematics->m_momentum.Py();}
		double pz(void) const{return dKinematics->m_momentum.Pz();}
		double x(void) const{return dKinematics->m_position.X();}
		double y(void) const{return dKinematics->m_position.Y();}
		double z(void) const{return dKinematics->m_position.Z();}

		//derived quantities
		double energy(void) const{return sqrt(mass()*mass() + pmag2());}
		double pperp(void) const{return sqrt(px()*px() + py()*py());}
		double pperp2(void) const{return px()*px() + py()*py();}
		double pmag(void) const{return dKinematics->m_momentum.Mag();}
		double pmag2(void) const{return dKinematics->m_momentum.Mag2();}
		DLorentzVector lorentzMomentum(void) const{return DLorentzVector(momentum(), energy());}

		//matrices, tracking
		const TMatrixFSym* errorMatrix(void) const{return dKinematics->m_errorMatrix;}
		const TMatrixFSym* TrackingErrorMatrix(void) const{return dNonKinematics->m_TrackingErrorMatrix;}
		bool forwardParmFlag(void) const{return dNonKinematics->m_use_forward_parameters;}
		void TrackingStateVector(double aVec[5]) const;

		//SETTERS
		void setPID(Particle_t locPID){dNonKinematics->m_pid = locPID;}
		void setMomentum(const DVector3& aMomentum){dKinematics->m_momentum = aMomentum;}
		void setPosition(const DVector3& aPosition){dKinematics->m_position = aPosition;}
		void setTime(double locTime){dKinematics->m_time = locTime;}

		//Tracking
		void setForwardParmFlag(bool aFlag){dNonKinematics->m_use_forward_parameters = aFlag;}
		void setErrorMatrix(const TMatrixFSym* aMatrix){dKinematics->m_errorMatrix = aMatrix;}
		void setTrackingErrorMatrix(const TMatrixFSym* aMatrix){dNonKinematics->m_TrackingErrorMatrix = aMatrix;}
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

		struct DKinematics
		{
			DKinematics(void);

			DVector3 m_momentum;
			DVector3 m_position;
			double m_time; // Time of the track propagated at m_position

			// Tracking information //The setter's are responsible for managing the matrix memory!  These are NEVER the owners.
			const TMatrixFSym* m_errorMatrix;   // Order is (px, py, pz, x, y, z, t)
		};

		struct DNonKinematics
		{
			// NONE OF THIS DEPENDS ON THE KINEMATIC FIT
			// so, this can be stored as a pointer & shared between multiple objects instead of stored separately for each
			DNonKinematics(void);

			Particle_t m_pid;

			const TMatrixFSym *m_TrackingErrorMatrix;  // order is q/pt,phi,tanl,D,z
			bool m_use_forward_parameters; // Flag indicating the use of the forward parameterization (x,y,tx,ty,q/p)
			double m_TrackingStateVector[5]; // order is q/pt,phi,tanl,D,z
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		shared_ptr<DKinematics> dKinematics;
		shared_ptr<DNonKinematics> dNonKinematics;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DKinematicData::DKinematicData(void) :
dKinematics(make_shared<DKinematics>()), dNonKinematics(make_shared<DNonKinematics>())
{}

inline DKinematicData::DKinematicData(const DKinematicData& locSourceData, bool locShareNonKinematicsFlag, bool locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dKinematics = locShareKinematicsFlag ? locSourceData->dKinematics : make_shared<DKinematics>(*(locSourceData->dKinematics));
	dNonKinematics = locShareNonKinematicsFlag ? locSourceData->dNonKinematics : make_shared<DNonKinematics>(*(locSourceData->dNonKinematics));
}

inline DKinematicData& DKinematicData::operator=(const DKinematicData& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dKinematics = make_shared<DKinematics>(*(locSourceData->dKinematics));
	dNonKinematics = make_shared<DNonKinematics>(*(locSourceData->dNonKinematics));
}

inline DKinematicData::DKinematics::DKinematics(void) :
m_momentum(DVector3()), m_position(DVector3()), m_time(0.0), m_errorMatrix(nullptr)
{}

inline DKinematicData::DNonKinematics::DNonKinematics(void) :
m_pid(Unknown), m_TrackingErrorMatrix(nullptr), m_use_forward_parameters(false), m_TrackingStateVector{0.0, 0.0, 0.0, 0.0, 0.0}
{}

/********************************************************************** GETTERS ************************************************************************/

inline void DKinematicData::TrackingStateVector(double aVec[5]) const
{
	for (unsigned int i = 0; i < 5; ++i)
		aVec[i] = dNonKinematics->m_TrackingStateVector[i];
}

/********************************************************************** SETTERS ************************************************************************/

inline void DKinematicData::setTrackingStateVector(double a1, double a2, double a3, double a4, double a5)
{
	dNonKinematics->m_TrackingStateVector[0]=a1;
	dNonKinematics->m_TrackingStateVector[1]=a2;
	dNonKinematics->m_TrackingStateVector[2]=a3;
	dNonKinematics->m_TrackingStateVector[3]=a4;
	dNonKinematics->m_TrackingStateVector[4]=a5;
}

#endif /* _DKINEMATICDATA_ */
