#ifndef _DKinFitParticle_
#define _DKinFitParticle_

#include <set>
#include <math.h>
#include <memory>

#include "DResettable.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixFSym.h"

using namespace std;

enum DKinFitParticleType
{
	d_DetectedParticle = 0, 
	d_BeamParticle,
	d_TargetParticle,
	d_DecayingParticle, 
	d_MissingParticle
};

enum DKinFitPullType
{
	d_EPull = 0,
	d_PxPull,
	d_PyPull,
	d_PzPull,
	d_XxPull,
	d_XyPull,
	d_XzPull,
	d_TPull
};

class DKinFitter;
class DKinFitUtils;
class DKinFitConstraint_Spacetime;
class DKinFitConstraint_Vertex;
class DKinFitConstraint_P4;
class DKinFitConstraint_Mass;
class DKinFitUtils_GlueX;

class DKinFitParticle : public DResettable
{
	friend class DKinFitter;
	friend class DKinFitUtils;
	friend class DKinFitConstraint_Spacetime;
	friend class DKinFitConstraint_Vertex;
	friend class DKinFitConstraint_P4;
	friend class DKinFitConstraint_Mass;
	friend class DKinFitUtils_GlueX;

	public:

		//STRUCTORS
		DKinFitParticle(void);
		~DKinFitParticle(void){}

		//RESET AND PRINT
		void Reset(void);
		void Release(void);
		void Print_ParticleParams(void) const;

		//SETTERS
		void Set_KinFitParticleType(DKinFitParticleType locKinFitParticleType){dKinFitParticleType = locKinFitParticleType;}
		void Set_PID(int locPID){dPID = locPID;}
		void Set_Charge(char locCharge){dCharge = locCharge;}
		void Set_Mass(double locMass){dMass = locMass;}

		void Set_Position(TVector3 locPosition){dSpacetimeVertex.SetVect(locPosition);}
		void Set_Time(double locTime){dSpacetimeVertex.SetT(locTime);}
		void Set_SpacetimeVertex(TLorentzVector locSpacetimeVertex){dSpacetimeVertex = locSpacetimeVertex;}

		void Set_Momentum(TVector3 locMomentum){dMomentum = locMomentum;}
		void Set_CovarianceMatrix(const shared_ptr<TMatrixFSym>& locCovarianceMatrix){dCovarianceMatrix = std::const_pointer_cast<const TMatrixFSym>(locCovarianceMatrix);}
		void Set_CovarianceMatrix(const shared_ptr<const TMatrixFSym>& locCovarianceMatrix){dCovarianceMatrix = locCovarianceMatrix;}
		void Set_ShowerEnergy(double locShowerEnergy){dShowerEnergy = locShowerEnergy;}
		void Set_PathLength(double locPathLength){dPathLength = locPathLength;}
		void Set_PathLengthUncertainty(double locPathLengthUncertainty){dPathLengthUncertainty = locPathLengthUncertainty;}
		void Set_RestFrameLifetimeUncertainty(double locRestFrameLifetimeUncertainty){dRestFrameLifetimeUncertainty = locRestFrameLifetimeUncertainty;}
		void Set_RestFrameLifetime(double locRestFrameLifetime){dRestFrameLifetime = locRestFrameLifetime;}

		void Set_CommonVertex(TVector3 locCommonVertex){dCommonSpacetimeVertex.SetVect(locCommonVertex);}
		void Set_CommonTime(double locCommonTime){dCommonSpacetimeVertex.SetT(locCommonTime);}
		void Set_CommonSpacetimeVertex(TLorentzVector locCommonSpacetimeVertex){dCommonSpacetimeVertex = locCommonSpacetimeVertex;}

		void Set_PxParamIndex(char locPxParamIndex){dPxParamIndex = locPxParamIndex;}
		void Set_VxParamIndex(char locVxParamIndex){dVxParamIndex = locVxParamIndex;}
		void Set_TParamIndex(char locTParamIndex){dTParamIndex = locTParamIndex;}
		void Set_CommonVxParamIndex(char locCommonVxParamIndex){dCommonVxParamIndex = locCommonVxParamIndex;}
		void Set_CommonTParamIndex(char locCommonTParamIndex){dCommonTParamIndex = locCommonTParamIndex;}
		void Set_EParamIndex(char locEParamIndex){dEParamIndex = locEParamIndex;}

		void Set_IsNeutralShowerFlag(bool locIsNeutralShowerFlag){dIsNeutralShowerFlag = locIsNeutralShowerFlag;}
		void Set_VertexConstraintFlag(unsigned char locVertexConstraintFlag){dVertexConstraintFlag = locVertexConstraintFlag;}

		void Set_VertexP4AtProductionVertex(bool locVertexP4AtProductionVertex){dVertexP4AtProductionVertex = locVertexP4AtProductionVertex;}
		void Set_FromInitialState(const set<shared_ptr<DKinFitParticle>>& locFromInitialState){dFromInitialState = locFromInitialState;}
		void Set_FromFinalState(const set<shared_ptr<DKinFitParticle>>& locFromFinalState){dFromFinalState = locFromFinalState;}

		//GETTERS
		int Get_PID(void) const{return dPID;}
		double Get_Energy(void) const{return sqrt(dMass*dMass + dMomentum.Mag2());}
		char Get_Charge(void) const{return dCharge;}
		TLorentzVector Get_P4(void) const{return TLorentzVector(dMomentum, Get_Energy());}
		TVector3 Get_Momentum(void) const{return dMomentum;}
		TVector3 Get_Position(void) const{return dSpacetimeVertex.Vect();}
		double Get_Mass(void) const{return dMass;}
		double Get_Beta(void) const{return dMomentum.Mag()/(Get_Energy());}
		double Get_Time(void) const{return dSpacetimeVertex.T();}
		double Get_ShowerEnergy(void) const{return dShowerEnergy;}
		double Get_PathLength(void) const{return dPathLength;}
		double Get_PathLengthUncertainty(void) const{return dPathLengthUncertainty;}
		double Get_RestFrameLifetimeUncertainty(void) const{return dRestFrameLifetimeUncertainty;}
		double Get_RestFrameLifetime(void) const{return dRestFrameLifetime;}
		shared_ptr<const TMatrixFSym> Get_CovarianceMatrix(void) const{return dCovarianceMatrix;}
		TLorentzVector Get_SpacetimeVertex(void) const{return dSpacetimeVertex;}
		TVector3 Get_CommonVertex(void) const{return dCommonSpacetimeVertex.Vect();}
		double Get_CommonTime(void) const{return dCommonSpacetimeVertex.T();}
		TLorentzVector Get_CommonSpacetimeVertex(void) const{return dCommonSpacetimeVertex;}

		unsigned char Get_VertexConstraintFlag(void) const{return dVertexConstraintFlag;}
		bool Get_FitCommonVertexFlag(void) const{return (dCommonVxParamIndex >= 0);}
		bool Get_FitCommonTimeFlag(void) const{return (dCommonTParamIndex >= 0);}
		DKinFitParticleType Get_KinFitParticleType(void) const{return dKinFitParticleType;}

		set<shared_ptr<DKinFitParticle>> Get_FromInitialState(void) const{return dFromInitialState;}
		set<shared_ptr<DKinFitParticle>> Get_FromFinalState(void) const{return dFromFinalState;}
		set<shared_ptr<DKinFitParticle>> Get_FromAllParticles(void) const;

		bool Get_VertexP4AtProductionVertex(void) const{return dVertexP4AtProductionVertex;}

		char Get_PxParamIndex(void) const{return dPxParamIndex;}
		char Get_VxParamIndex(void) const{return dVxParamIndex;}
		char Get_TParamIndex(void) const{return dTParamIndex;}
		char Get_CommonVxParamIndex(void) const{return dCommonVxParamIndex;}
		char Get_CommonTParamIndex(void) const{return dCommonTParamIndex;}
		char Get_EParamIndex(void) const{return dEParamIndex;}

		int Get_CovMatrixEParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() >= 7) ? -1 : 0));}
		int Get_CovMatrixPxParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() >= 7) ? 0 : -1));}
		int Get_CovMatrixVxParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() >= 7) ? 3 : 1));}
		int Get_CovMatrixTParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() >= 7) ? 6 : 4));}

		bool Get_IsNeutralShowerFlag(void) const{return dIsNeutralShowerFlag;}

	private:

		DKinFitParticleType dKinFitParticleType;

		int dPID; //PDG PID
		char dCharge;
		double dMass;

		//p, x, & t must all coincide: t & p at point x (for charged tracks p is a function of x in a b-field!)
		TLorentzVector dSpacetimeVertex;
		TLorentzVector dCommonSpacetimeVertex; //if not in vertex/time fit, will be same as dSpacetimeVertex

		double dShowerEnergy;
		TVector3 dMomentum; //must be the value of the momentum at dSpacetimeVertex

		//is 7x7 for detected charged particles, either 7x7 (particles) or 5x5 (showers) for neutrals
		//for decaying particles, is 11x11 if involved in 2 vertex fits (otherwise 7x7): includes common vertex
		shared_ptr<const TMatrixFSym> dCovarianceMatrix; 
		//5x5 format: E, x, y, z, t
		//7x7 format: px, py, pz, x, y, z, t
		//11x11 format: px, py, pz, x, y, z, t, common x, common y, common z, common t

		double dPathLength;
		double dPathLengthUncertainty;
		double dRestFrameLifetime; //is 0 unless decaying particle in 2 vertex fits
		double dRestFrameLifetimeUncertainty; //is 0 unless decaying particle in 2 vertex fits

		unsigned char dVertexConstraintFlag; //unused unless in vertex fit //can choose between equations //only for non-accelerating particles not constrained in time

		char dEParamIndex;
		char dPxParamIndex;
		char dVxParamIndex;
		char dTParamIndex;
		char dCommonVxParamIndex;
		char dCommonTParamIndex;

		//Decaying particles are reconstructed from these particles //ignored if not decaying
		set<shared_ptr<DKinFitParticle>> dFromInitialState;
		set<shared_ptr<DKinFitParticle>> dFromFinalState;

		//true if the object's p3, v3, & t are defined at its production vertex (& common v & t are at decay vertex). else at it's decay vertex (& common v & t are at production vertex)
			//note: if a decaying particle is not in a vertex fit, then this quantity doesn't matter (default true)
		bool dVertexP4AtProductionVertex;

		bool dIsNeutralShowerFlag;
};

inline DKinFitParticle::DKinFitParticle(void)
{
	Reset();
}

inline void DKinFitParticle::Reset(void)
{
	dPID = 0;
	dCharge = 0;
	dMass = 0.0;
	dSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dCommonSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dShowerEnergy = 0.0;
	dMomentum.SetXYZ(0.0, 0.0, 0.0);
	dCovarianceMatrix = NULL;
	dPathLength = 0.0;
	dPathLengthUncertainty = 0.0;
	dRestFrameLifetime = 0.0;
	dRestFrameLifetimeUncertainty = 0.0;

	dVertexConstraintFlag = 0;

	dKinFitParticleType = d_DetectedParticle;

	dEParamIndex = -1;
	dPxParamIndex = -1;
	dVxParamIndex = -1;
	dTParamIndex = -1;
	dCommonVxParamIndex = -1;
	dCommonTParamIndex = -1;

	dFromInitialState.clear();
	dFromFinalState.clear();

	dVertexP4AtProductionVertex = false;
	dIsNeutralShowerFlag = false;
}

inline void DKinFitParticle::Release(void)
{
	dCovarianceMatrix = nullptr;
	dFromInitialState.clear();
	dFromFinalState.clear();
}

inline set<shared_ptr<DKinFitParticle>> DKinFitParticle::Get_FromAllParticles(void) const
{
	//get all of the particles this particle is derived from, excluding decaying particles
	set<shared_ptr<DKinFitParticle>> locFromAllParticles;

	//from initial state
	for(auto& locKinFitParticle : dFromInitialState)
	{
		locFromAllParticles.insert(locKinFitParticle);
		if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			auto locNewParticles = locKinFitParticle->Get_FromAllParticles();
			locFromAllParticles.insert(locNewParticles.begin(), locNewParticles.end());
		}
	}

	//from initial state
	for(auto& locKinFitParticle : dFromFinalState)
	{
		locFromAllParticles.insert(locKinFitParticle);
		if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			auto locNewParticles = locKinFitParticle->Get_FromAllParticles();
			locFromAllParticles.insert(locNewParticles.begin(), locNewParticles.end());
		}
	}

	return locFromAllParticles;
}

inline void DKinFitParticle::Print_ParticleParams(void) const
{
	cout << "DKinFitParticle: Particle Type Enum, pointer: " << dKinFitParticleType << ", " << this << endl;

	cout << "DKinFitParticle: Particle PID, Q, Mass = " << dPID << ", " << int(dCharge) << ", " << dMass << endl;
	cout << "DKinFitParticle: Particle P3, V3, T, path length = " << dMomentum.Px() << ", " << dMomentum.Py() << ", " << dMomentum.Pz() << ", ";
	cout << dSpacetimeVertex.X() << ", " << dSpacetimeVertex.Y() << ", " << dSpacetimeVertex.Z() << ", " << dSpacetimeVertex.T() << ", " << dPathLength << endl;
	cout << "DKinFitParticle: Particle Common V3, Common T, ShowerE = " << dCommonSpacetimeVertex.X() << ", " << dCommonSpacetimeVertex.Y();
	cout << ", " << dCommonSpacetimeVertex.Z() << ", " << dCommonSpacetimeVertex.T() << ", " << dShowerEnergy << endl;

	cout << "DKinFitParticle: FitCommonVertexFlag, FitCommonTimeFlag, dVertexP4AtProductionVertex, dIsNeutralShowerFlag = ";
	cout << Get_FitCommonVertexFlag() << ", " << Get_FitCommonTimeFlag() << ", " << dVertexP4AtProductionVertex << ", " << dIsNeutralShowerFlag << endl;
	if(dCovarianceMatrix != NULL)
	{
		cout << "DKinFitParticle: CovMatrix Diagonal Terms: ";
		for(int loc_i = 0; loc_i < dCovarianceMatrix->GetNcols(); ++loc_i)
			cout << (*dCovarianceMatrix)(loc_i, loc_i) << ", ";
		cout << endl;
	}

	cout << "DKinFitParticle: Particle E, Px, Vx, Common Vx, T, Common T indices = " << int(dEParamIndex) << ", " << int(dPxParamIndex) << ", ";
	cout << int(dVxParamIndex) << ", " << int(dCommonVxParamIndex) << ", " << int(dTParamIndex) << ", " << int(dCommonTParamIndex) << endl;

	cout << "dFromInitialState size, PIDs: " << dFromInitialState.size();
	for(auto& locKinFitParticle : dFromInitialState)
		cout << ", " << locKinFitParticle->Get_PID();
	cout << endl;

	cout << "dFromFinalState size, PIDs: " << dFromFinalState.size();
	for(auto& locKinFitParticle : dFromFinalState)
		cout << ", " << locKinFitParticle->Get_PID();
	cout << endl;

	auto locFromAllParticles = Get_FromAllParticles();
	cout << "FromAllButDecaying size, PIDs: " << locFromAllParticles.size();
	for(auto& locKinFitParticle : locFromAllParticles)
		cout << ", " << locKinFitParticle->Get_PID();
	cout << endl;
}

#endif // _DKinFitParticle_

