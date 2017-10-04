#ifndef _DKINEMATICDATA_
#define _DKINEMATICDATA_

#include <memory>

#include <JANA/JObject.h>

#include "GlueX.h"    

#include "DVector3.h"    
#include "DLorentzVector.h" 
#include "particleType.h" 
#include "TMatrixFSym.h"
#include "DResettable.h"
#include "DResourcePool.h"

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.9792458
#endif

using namespace std;
using namespace jana;

class DKinematicData : public JObject, public DResettable
{
	public:

		// constructors and destructor
		DKinematicData(void);
		DKinematicData(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition = DVector3(), double locTime = 0.0, const shared_ptr<const TMatrixFSym>& locErrorMatrix = nullptr);
		DKinematicData(const DKinematicData& locSourceData, bool locShareKinematicsFlag = false);
		virtual ~DKinematicData(void) {};

		//Assignment operator
		DKinematicData& operator=(const DKinematicData& locSourceData);
		void Share_FromInput_Kinematics(const DKinematicData* locSourceData);

		//Reset & Release
		virtual void Reset(void);
		virtual void Release(void);

		//GETTERS
		Particle_t PID(void) const{return dKinematicInfo->dPID;}
		double mass(void) const{return ParticleMass(dKinematicInfo->dPID);}
		double charge(void) const{return ParticleCharge(dKinematicInfo->dPID);}
		const DVector3& momentum(void) const{return dKinematicInfo->dMomentum;}
		const DVector3& position(void) const{return dKinematicInfo->dPosition;}
		double time(void) const{return dKinematicInfo->dTime;}
		shared_ptr<const TMatrixFSym> errorMatrix(void) const{return dErrorMatrix;}

		//components
		double px(void) const{return dKinematicInfo->dMomentum.Px();}
		double py(void) const{return dKinematicInfo->dMomentum.Py();}
		double pz(void) const{return dKinematicInfo->dMomentum.Pz();}
		double x(void) const{return dKinematicInfo->dPosition.X();}
		double y(void) const{return dKinematicInfo->dPosition.Y();}
		double z(void) const{return dKinematicInfo->dPosition.Z();}

		//derived quantities
		double energy(void) const{return sqrt(mass()*mass() + pmag2());}
		double pperp(void) const{return sqrt(px()*px() + py()*py());}
		double pperp2(void) const{return px()*px() + py()*py();}
		double pmag(void) const{return dKinematicInfo->dMomentum.Mag();}
		double pmag2(void) const{return dKinematicInfo->dMomentum.Mag2();}
		DLorentzVector lorentzMomentum(void) const{return DLorentzVector(momentum(), energy());}
		DLorentzVector x4(void) const{return DLorentzVector(position(), time());}

		//SETTERS
		void Set_Members(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition = DVector3(), double locTime = 0.0, const shared_ptr<const TMatrixFSym>& locErrorMatrix = nullptr);
		void setPID(Particle_t locPID){dKinematicInfo->dPID = locPID;}
		void setMomentum(const DVector3& aMomentum){dKinematicInfo->dMomentum = aMomentum;}
		void setPosition(const DVector3& aPosition){dKinematicInfo->dPosition = aPosition;}
		void setTime(double locTime){dKinematicInfo->dTime = locTime;}
		void setErrorMatrix(const shared_ptr<const TMatrixFSym>& aMatrix){dErrorMatrix = aMatrix;}
		void setErrorMatrix(const shared_ptr<TMatrixFSym>& aMatrix){dErrorMatrix = std::const_pointer_cast<const TMatrixFSym>(aMatrix);}

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

		class DKinematicInfo : public DResettable
		{
			public:
				//CONSTRUCTORS
				DKinematicInfo(void) = default;
				DKinematicInfo(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition = DVector3(), double locTime = 0.0);

				void Set_Members(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition = DVector3(), double locTime = 0.0);
				void Reset(void);
				void Release(void){};

				//MEMBERS
				Particle_t dPID = Unknown;
				DVector3 dMomentum;
				DVector3 dPosition;
				double dTime = 0.0; // Time of the track propagated at dPosition
		};

	private:

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		//This is done because sometimes a new object is needed (e.g. DChargedTrackHypothesis) for which this info hasn't changed (from DTrackTimeBased)
		//Thus, just share this between the two objects, instead of doubling the memory usage
		//By inheriting this class, you also get to share the same interface
		shared_ptr<DKinematicInfo> dKinematicInfo = nullptr;

		//This member must be separate from DKinematicInfo or else it may be recycled on a different thread from which it was created, causing a race condition!
		shared_ptr<const TMatrixFSym> dErrorMatrix = nullptr; // Order is (px, py, pz, x, y, z, t)

		static thread_local shared_ptr<DResourcePool<DKinematicInfo>> dResourcePool_KinematicInfo;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DKinematicData::DKinematicData(void) : dKinematicInfo(dResourcePool_KinematicInfo->Get_SharedResource()), dErrorMatrix(nullptr){}

inline DKinematicData::DKinematicData(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition, double locTime, const shared_ptr<const TMatrixFSym>& locErrorMatrix) :
		dKinematicInfo(dResourcePool_KinematicInfo->Get_SharedResource()), dErrorMatrix(locErrorMatrix)
{
	dKinematicInfo->Set_Members(locPID, locMomentum, locPosition, locTime);
}

inline void DKinematicData::Share_FromInput_Kinematics(const DKinematicData* locSourceData)
{
	dKinematicInfo = const_cast<DKinematicData*>(locSourceData)->dKinematicInfo;
	dErrorMatrix = locSourceData->dErrorMatrix;
}

inline DKinematicData::DKinematicData(const DKinematicData& locSourceData, bool locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	if(locShareKinematicsFlag)
		dKinematicInfo = locSourceData.dKinematicInfo;
	else
	{
		dKinematicInfo = dResourcePool_KinematicInfo->Get_SharedResource();
		*dKinematicInfo = *(locSourceData.dKinematicInfo);
	}
	dErrorMatrix = locSourceData.dErrorMatrix; //it's const so it can't be modified: share and then replace if desired
}

inline DKinematicData& DKinematicData::operator=(const DKinematicData& locSourceData)
{
	//Guard against self-assignment
	if(dKinematicInfo == locSourceData.dKinematicInfo)
		return *this;

	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dKinematicInfo = dResourcePool_KinematicInfo->Get_SharedResource();
	*dKinematicInfo = *(locSourceData.dKinematicInfo);

	dErrorMatrix = locSourceData.dErrorMatrix; //it's const so it can't be modified: share and then replace if desired
	return *this;
}

inline DKinematicData::DKinematicInfo::DKinematicInfo(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition, double locTime) :
dPID(locPID), dMomentum(locMomentum), dPosition(locPosition), dTime(locTime) {}

/*********************************************************************** RESET & RELEASE *************************************************************************/

inline void DKinematicData::Set_Members(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition, double locTime, const shared_ptr<const TMatrixFSym>& locErrorMatrix)
{
	dKinematicInfo->Set_Members(locPID, locMomentum, locPosition, locTime);
	dErrorMatrix = locErrorMatrix;
}

inline void DKinematicData::DKinematicInfo::Set_Members(Particle_t locPID, const DVector3& locMomentum, DVector3 locPosition, double locTime)
{
	dPID = locPID;
	dMomentum = locMomentum;
	dPosition = locPosition;
	dTime = locTime;
}

inline void DKinematicData::Reset(void)
{
	dKinematicInfo = dResourcePool_KinematicInfo->Get_SharedResource(); //not safe to reset individually, since you don't know what it's shared with
	dErrorMatrix = nullptr;
	ClearAssociatedObjects();
}

inline void DKinematicData::Release(void)
{
	dKinematicInfo = nullptr;
	dErrorMatrix = nullptr;
	ClearAssociatedObjects();
}

inline void DKinematicData::DKinematicInfo::Reset(void)
{
	dPID = Unknown;
	dMomentum = DVector3();
	dPosition = DVector3();
	dTime = 0.0;
}

#endif /* _DKINEMATICDATA_ */
