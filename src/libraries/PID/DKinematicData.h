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
		DKinematicData(const DKinematicData& locSourceData, bool locShareKinematicsFlag = false);
		virtual ~DKinematicData(void) {};

		//Assignment operator
		DKinematicData& operator=(const DKinematicData& locSourceData);

		//Reset
		virtual void Reset(void);

		//GETTERS
		Particle_t PID(void) const{return dKinematicInfo->m_pid;}
		double mass(void) const{return ParticleMass(dKinematicInfo->m_pid);}
		double charge(void) const{return ParticleCharge(dKinematicInfo->m_pid);}
		const DVector3& momentum(void) const{return dKinematicInfo->m_momentum;}
		const DVector3& position(void) const{return dKinematicInfo->m_position;}
		double time(void) const{return dKinematicInfo->m_time;}
		const TMatrixFSym* errorMatrix(void) const{return dKinematicInfo->m_errorMatrix;}

		//components
		double px(void) const{return dKinematicInfo->m_momentum.Px();}
		double py(void) const{return dKinematicInfo->m_momentum.Py();}
		double pz(void) const{return dKinematicInfo->m_momentum.Pz();}
		double x(void) const{return dKinematicInfo->m_position.X();}
		double y(void) const{return dKinematicInfo->m_position.Y();}
		double z(void) const{return dKinematicInfo->m_position.Z();}

		//derived quantities
		double energy(void) const{return sqrt(mass()*mass() + pmag2());}
		double pperp(void) const{return sqrt(px()*px() + py()*py());}
		double pperp2(void) const{return px()*px() + py()*py();}
		double pmag(void) const{return dKinematicInfo->m_momentum.Mag();}
		double pmag2(void) const{return dKinematicInfo->m_momentum.Mag2();}
		DLorentzVector lorentzMomentum(void) const{return DLorentzVector(momentum(), energy());}

		//SETTERS
		void setPID(Particle_t locPID){dKinematicInfo->m_pid = locPID;}
		void setMomentum(const DVector3& aMomentum){dKinematicInfo->m_momentum = aMomentum;}
		void setPosition(const DVector3& aPosition){dKinematicInfo->m_position = aPosition;}
		void setTime(double locTime){dKinematicInfo->m_time = locTime;}
		void setErrorMatrix(const TMatrixFSym* aMatrix){dKinematicInfo->m_errorMatrix = aMatrix;}

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

		struct DKinematicInfo
		{
			DKinematicInfo(void);

			Particle_t m_pid;
			DVector3 m_momentum;
			DVector3 m_position;
			double m_time; // Time of the track propagated at m_position

			// Tracking information //The setter's are responsible for managing the matrix memory!  These are NEVER the owners.
			const TMatrixFSym* m_errorMatrix;   // Order is (px, py, pz, x, y, z, t)
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		//This is done because sometimes a new object is needed (e.g. DChargedTrackHypothesis) for which this info hasn't changed (from DTrackTimeBased)
		//Thus, just share this between the two objects, instead of doubling the memory usage
		shared_ptr<DKinematicInfo> dKinematicInfo;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DKinematicData::DKinematicData(void) : dKinematicInfo(make_shared<DKinematicInfo>()) {}

inline DKinematicData::DKinematicData(const DKinematicData& locSourceData, bool locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dKinematicInfo = locShareKinematicsFlag ? locSourceData.dKinematicInfo : make_shared<DKinematicInfo>(*(locSourceData.dKinematicInfo));
}

inline DKinematicData& DKinematicData::operator=(const DKinematicData& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dKinematicInfo = make_shared<DKinematicInfo>(*(locSourceData.dKinematicInfo));
	return *this;
}

inline DKinematicData::DKinematicInfo::DKinematicInfo(void) :
m_pid(Unknown), m_momentum(DVector3()), m_position(DVector3()), m_time(0.0), m_errorMatrix(nullptr) {}

/*********************************************************************** RESET *************************************************************************/

inline void DKinematicData::Reset(void)
{
	dKinematicInfo = make_shared<DKinematicInfo>(); //not safe to reset individually, since you don't know what it's shared with
	ClearAssociatedObjects();
}

#endif /* _DKINEMATICDATA_ */
