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
		virtual ~DKinematicData(void) {};

		//Reset
		virtual void Reset(void);

		//GETTERS
		Particle_t PID(void) const{return dNonKinFitInfo->m_pid;}
		double mass(void) const{return ParticleMass(dNonKinFitInfo->m_pid);}
		double charge(void) const{return ParticleCharge(dNonKinFitInfo->m_pid);}
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
		double energy(void) const{return sqrt(mass()*mass() + pmag2());}
		double pperp(void) const{return sqrt(px()*px() + py()*py());}
		double pperp2(void) const{return px()*px() + py()*py();}
		double pmag(void) const{return m_momentum.Mag();}
		double pmag2(void) const{return m_momentum.Mag2();}
		DLorentzVector lorentzMomentum(void) const{return DLorentzVector(m_momentum, energy());}

		//matrices, tracking
		const TMatrixFSym* errorMatrix(void) const{return m_errorMatrix;}
		const TMatrixFSym* TrackingErrorMatrix(void) const{return dNonKinFitInfo->m_TrackingErrorMatrix;}
		bool forwardParmFlag(void) const{return dNonKinFitInfo->m_use_forward_parameters;}
		void TrackingStateVector(double aVec[5]) const;

		//pathlength
		double pathLength(void) const{return dNonKinFitInfo->m_pathLength;}

		//SETTERS
		void setPID(Particle_t locPID){dNonKinFitInfo->m_pid = locPID;}
		void setMomentum(const DVector3& aMomentum){m_momentum = aMomentum;}
		void setPosition(const DVector3& aPosition){m_position = aPosition;}
		void setTime(double locTime){m_time = locTime;}

		//pathlength
		void setPathLength(double apathLength){dNonKinFitInfo->m_pathLength = apathLength;}

		//Tracking
		void setForwardParmFlag(bool aFlag){dNonKinFitInfo->m_use_forward_parameters = aFlag;}
		void setErrorMatrix(const TMatrixFSym* aMatrix){m_errorMatrix = aMatrix;}
		void setTrackingErrorMatrix(const TMatrixFSym* aMatrix){dNonKinFitInfo->m_TrackingErrorMatrix = aMatrix;}
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

		struct DNonKinFitInfo
		{
			// NONE OF THIS DEPENDS ON THE KINEMATIC FIT
			// so, this can be stored as a pointer & shared between multiple objects instead of stored separately for each

			Particle_t m_pid;
			double m_pathLength; /// Flight path length (cm)

			const TMatrixFSym *m_TrackingErrorMatrix;  // order is q/pt,phi,tanl,D,z
			bool m_use_forward_parameters; // Flag indicating the use of the forward parameterization (x,y,tx,ty,q/p)
			double m_TrackingStateVector[5]; // order is q/pt,phi,tanl,D,z
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		shared_ptr<DNonKinFitInfo> dNonKinFitInfo;

		DVector3 m_momentum;
		DVector3 m_position;
		double m_time; // Time of the track propagated at m_position

		// Tracking information //The setter's are responsible for managing the matrix memory!  These are NEVER the owners.
		const TMatrixFSym* m_errorMatrix;   // Order is (px, py, pz, x, y, z, t)
};

/******************************************************************** CONSTRUCTORS *********************************************************************/

inline DKinematicData::DKinematicData(void) :
dNonKinFitInfo(make_shared<DNonKinFitInfo>()), m_momentum(DVector3()), m_position(DVector3()), m_time(0.0), m_errorMatrix(nullptr)
{}

inline DKinematicData::DNonKinFitInfo::DNonKinFitInfo(void) :
m_pid(Unknown), m_pathLength(0.0), m_TrackingErrorMatrix(nullptr), m_use_forward_parameters(false), m_TrackingStateVector{0.0, 0.0, 0.0, 0.0, 0.0}
{}

/********************************************************************** GETTERS ************************************************************************/

inline void DKinematicData::TrackingStateVector(double aVec[5]) const
{
	for (unsigned int i = 0; i < 5; ++i)
		aVec[i] = dNonKinFitInfo->m_TrackingStateVector[i];
}

/********************************************************************** SETTERS ************************************************************************/

inline void DKinematicData::setTrackingStateVector(double a1, double a2, double a3, double a4, double a5)
{
	dNonKinFitInfo->m_TrackingStateVector[0]=a1;
	dNonKinFitInfo->m_TrackingStateVector[1]=a2;
	dNonKinFitInfo->m_TrackingStateVector[2]=a3;
	dNonKinFitInfo->m_TrackingStateVector[3]=a4;
	dNonKinFitInfo->m_TrackingStateVector[4]=a5;
}

#endif /* _DKINEMATICDATA_ */
