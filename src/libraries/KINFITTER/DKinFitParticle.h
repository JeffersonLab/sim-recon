#ifndef _DKinFitParticle_
#define _DKinFitParticle_

#include <set>
#include <math.h>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"

//#include "DKinFitConstraints.h"

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

// dCovarianceMatrix is owned by DKinFitUtils (DKinFitUtils is responsible for new/delete)
class DKinFitter;
class DKinFitUtils;
class DKinFitConstraint_Spacetime;
class DKinFitConstraint_Vertex;
class DKinFitConstraint_P4;
class DKinFitConstraint_Mass;

class DKinFitParticle
{
	friend class DKinFitter;
	friend class DKinFitUtils;
	friend class DKinFitConstraint_Spacetime;
	friend class DKinFitConstraint_Vertex;
	friend class DKinFitConstraint_P4;
	friend class DKinFitConstraint_Mass;

	public:

		int Get_PID(void) const{return dPID;}
		double Get_Energy(void) const{return sqrt(dMass*dMass + dMomentum.Mag2());}
		int Get_Charge(void) const{return dCharge;}
		TLorentzVector Get_P4(void) const{return TLorentzVector(dMomentum, Get_Energy());}
		TVector3 Get_Momentum(void) const{return dMomentum;}
		TVector3 Get_Position(void) const{return dSpacetimeVertex.Vect();}
		double Get_Mass(void) const{return dMass;}
		double Get_Beta(void) const{return dMomentum.Mag()/(Get_Energy());}
		double Get_Time(void) const{return dSpacetimeVertex.T();}
		double Get_ShowerEnergy(void) const{return dShowerEnergy;}
		double Get_PathLength(void) const{return dPathLength;}
		double Get_PathLengthUncertainty(void) const{return dPathLengthUncertainty;}
		const TMatrixDSym* Get_CovarianceMatrix(void) const{return dCovarianceMatrix;}
		TLorentzVector Get_SpacetimeVertex(void) const{return dSpacetimeVertex;}
		TVector3 Get_CommonVertex(void) const{return dCommonSpacetimeVertex.Vect();}
		double Get_CommonTime(void) const{return dCommonSpacetimeVertex.T();}
		TLorentzVector Get_CommonSpacetimeVertex(void) const{return dCommonSpacetimeVertex;}

		unsigned short int Get_VertexConstraintFlag(void) const{return dVertexConstraintFlag;}
		bool Get_FitCommonVertexFlag(void) const{return (dCommonVxParamIndex >= 0);}
		bool Get_FitCommonTimeFlag(void) const{return (dCommonTParamIndex >= 0);}
		DKinFitParticleType Get_KinFitParticleType(void) const{return dKinFitParticleType;}

		set<DKinFitParticle*> Get_FromInitialState(void) const{return dFromInitialState;}
		set<DKinFitParticle*> Get_FromFinalState(void) const{return dFromFinalState;}
		set<DKinFitParticle*> Get_FromAllButDecaying(void) const;

		bool Get_VertexP4AtProductionVertex(void) const{return dVertexP4AtProductionVertex;}

		int Get_PxParamIndex(void) const{return dPxParamIndex;}
		int Get_VxParamIndex(void) const{return dVxParamIndex;}
		int Get_TParamIndex(void) const{return dTParamIndex;}
		int Get_CommonVxParamIndex(void) const{return dCommonVxParamIndex;}
		int Get_CommonTParamIndex(void) const{return dCommonTParamIndex;}
		int Get_EParamIndex(void) const{return dEParamIndex;}

		int Get_CovMatrixEParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() == 7) ? -1 : 0));}
		int Get_CovMatrixPxParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() == 7) ? 0 : -1));}
		int Get_CovMatrixVxParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() == 7) ? 3 : 1));}
		int Get_CovMatrixTParamIndex(void) const{return ((dCovarianceMatrix == NULL) ? -1 : ((dCovarianceMatrix->GetNcols() == 7) ? 6 : 4));}

		bool Get_IsNeutralShowerFlag(void) const{return dIsNeutralShowerFlag;}

		void Print_ParticleParams(void) const;

	private:

		DKinFitParticle(void);
		~DKinFitParticle(void){}

		void Set_KinFitParticleType(DKinFitParticleType locKinFitParticleType){dKinFitParticleType = locKinFitParticleType;}
		void Set_PID(int locPID){dPID = locPID;}
		void Set_Charge(int locCharge){dCharge = locCharge;}
		void Set_Mass(double locMass){dMass = locMass;}

		void Set_Position(TVector3 locPosition){dSpacetimeVertex.SetVect(locPosition);}
		void Set_Time(double locTime){dSpacetimeVertex.SetT(locTime);}
		void Set_SpacetimeVertex(TLorentzVector locSpacetimeVertex){dSpacetimeVertex = locSpacetimeVertex;}

		void Set_Momentum(TVector3 locMomentum){dMomentum = locMomentum;}
		void Set_CovarianceMatrix(TMatrixDSym* locCovarianceMatrix){dCovarianceMatrix = locCovarianceMatrix;}
		void Set_ShowerEnergy(double locShowerEnergy){dShowerEnergy = locShowerEnergy;}
		void Set_PathLength(double locPathLength){dPathLength = locPathLength;}
		void Set_PathLengthUncertainty(double locPathLengthUncertainty){dPathLengthUncertainty = locPathLengthUncertainty;}

		void Set_CommonVertex(TVector3 locCommonVertex){dCommonSpacetimeVertex.SetVect(locCommonVertex);}
		void Set_CommonTime(double locCommonTime){dCommonSpacetimeVertex.SetT(locCommonTime);}
		void Set_CommonSpacetimeVertex(TLorentzVector locCommonSpacetimeVertex){dCommonSpacetimeVertex = locCommonSpacetimeVertex;}

		void Set_PxParamIndex(int locPxParamIndex){dPxParamIndex = locPxParamIndex;}
		void Set_VxParamIndex(int locVxParamIndex){dVxParamIndex = locVxParamIndex;}
		void Set_TParamIndex(int locTParamIndex){dTParamIndex = locTParamIndex;}
		void Set_CommonVxParamIndex(int locCommonVxParamIndex){dCommonVxParamIndex = locCommonVxParamIndex;}
		void Set_CommonTParamIndex(int locCommonTParamIndex){dCommonTParamIndex = locCommonTParamIndex;}
		void Set_EParamIndex(int locEParamIndex){dEParamIndex = locEParamIndex;}

		void Set_IsNeutralShowerFlag(bool locIsNeutralShowerFlag){dIsNeutralShowerFlag = locIsNeutralShowerFlag;}
		void Set_VertexConstraintFlag(unsigned short int locVertexConstraintFlag){dVertexConstraintFlag = locVertexConstraintFlag;}

		void Set_VertexP4AtProductionVertex(bool locVertexP4AtProductionVertex){dVertexP4AtProductionVertex = locVertexP4AtProductionVertex;}
		void Set_FromInitialState(const set<DKinFitParticle*>& locFromInitialState){dFromInitialState = locFromInitialState;}
		void Set_FromFinalState(const set<DKinFitParticle*>& locFromFinalState){dFromFinalState = locFromFinalState;}

		void Reset(void);

		DKinFitParticleType dKinFitParticleType;

		int dPID; //PDG PID
		int dCharge;
		double dMass;

		//p, x, & t must all coincide: t & p at point x (for charged tracks p is a function of x in a b-field!)
		TLorentzVector dSpacetimeVertex;
		TLorentzVector dCommonSpacetimeVertex; //if not in vertex/time fit, will be same as dSpacetimeVertex

		double dShowerEnergy;
		TVector3 dMomentum; //must be the value of the momentum at dSpacetimeVertex
		TMatrixDSym* dCovarianceMatrix; //is 7x7 for charged particles, and is either 7x7 or 5x5 for neutrals
		//7x7 format: px, py, pz, x, y, z, t
		//5x5 format: E, x, y, z, t

		//not sure if still needed
		double dPathLength;
		double dPathLengthUncertainty;

		unsigned short int dVertexConstraintFlag; //unused unless in vertex fit //can choose between equations //only for non-accelerating particles not constrained in time

		int dEParamIndex;
		int dPxParamIndex;
		int dVxParamIndex;
		int dTParamIndex;
		int dCommonVxParamIndex;
		int dCommonTParamIndex;

		//Decaying particles are reconstructed from these particles //ignored if not decaying
		set<DKinFitParticle*> dFromInitialState;
		set<DKinFitParticle*> dFromFinalState;

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

inline set<DKinFitParticle*> DKinFitParticle::Get_FromAllButDecaying(void) const
{
	//get all of the particles this particle is derived from, excluding decaying particles
	set<DKinFitParticle*> locFromAllButDecaying;

	//from initial state
	set<DKinFitParticle*>::const_iterator locIterator = dFromInitialState.begin();
	for(; locIterator != dFromInitialState.end(); ++locIterator)
	{
		DKinFitParticle* locKinFitParticle = *locIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			locFromAllButDecaying.insert(locKinFitParticle);
		else
		{
			set<DKinFitParticle*> locNewParticles = locKinFitParticle->Get_FromAllButDecaying();
			locFromAllButDecaying.insert(locNewParticles.begin(), locNewParticles.end());
		}
	}

	//from initial state
	for(locIterator = dFromFinalState.begin(); locIterator != dFromFinalState.end(); ++locIterator)
	{
		DKinFitParticle* locKinFitParticle = *locIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			locFromAllButDecaying.insert(locKinFitParticle);
		else
		{
			set<DKinFitParticle*> locNewParticles = locKinFitParticle->Get_FromAllButDecaying();
			locFromAllButDecaying.insert(locNewParticles.begin(), locNewParticles.end());
		}
	}

	return locFromAllButDecaying;
}

inline void DKinFitParticle::Print_ParticleParams(void) const
{
	cout << "DKinFitParticle: Particle Type Enum: " << dKinFitParticleType << endl;

	cout << "DKinFitParticle: Particle PID, Q, Mass = " << dPID << ", " << dCharge << ", " << dMass << endl;
	cout << "DKinFitParticle: Particle P3, V3, T = " << dMomentum.Px() << ", " << dMomentum.Py() << ", " << dMomentum.Pz() << ", ";
	cout << dSpacetimeVertex.X() << ", " << dSpacetimeVertex.Y() << ", " << dSpacetimeVertex.Z() << ", " << dSpacetimeVertex.T() << endl;
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

	cout << "DKinFitParticle: Particle E, Px, Vx, Common Vx, T, Common T indices = " << dEParamIndex << ", " << dPxParamIndex << ", " << dVxParamIndex << ", " << dCommonVxParamIndex << ", " << dTParamIndex << ", " << dCommonTParamIndex << endl;

	cout << "dFromInitialState size, PIDs: " << dFromInitialState.size();
	set<DKinFitParticle*>::const_iterator locIterator = dFromInitialState.begin();
	for(; locIterator != dFromInitialState.end(); ++locIterator)
		cout << ", " << (*locIterator)->Get_PID();
	cout << endl;

	cout << "dFromFinalState size, PIDs: " << dFromFinalState.size();
	for(locIterator = dFromFinalState.begin(); locIterator != dFromFinalState.end(); ++locIterator)
		cout << ", " << (*locIterator)->Get_PID();
	cout << endl;
}

#endif // _DKinFitParticle_

