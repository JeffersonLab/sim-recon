#ifndef _DKinFitConstraints_Spacetime_
#define _DKinFitConstraints_Spacetime_

#include <set>
#include <algorithm>

#include "TVector3.h"
#include "TLorentzVector.h"

#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"

using namespace std;

class DKinFitter;
class DKinFitUtils;

class DKinFitConstraint_Spacetime : public DKinFitConstraint_Vertex
{
	friend class DKinFitter;
	friend class DKinFitUtils;

	public:
		double Get_InitTimeGuess(void) const{return dInitTimeGuess;};
		void Set_InitTimeGuess(double locInitTimeGuess){dInitTimeGuess = locInitTimeGuess;};

		TLorentzVector Get_CommonSpacetime(void) const;
		double Get_CommonTime(void) const;

		int Get_CommonTParamIndex(void) const;

		set<DKinFitParticle*> Get_OnlyConstrainTimeParticles(void) const{return dOnlyConstrainTimeParticles;}
		set<DKinFitParticle*> Get_AllConstrainedParticles(void) const;
		set<DKinFitParticle*> Get_AllParticles(void) const;

		void Print_ConstraintInfo(void) const;

	private:
		DKinFitConstraint_Spacetime(void);
		~DKinFitConstraint_Spacetime(void){}

		void Reset(void);

		void Set_CommonTime(double locTime);
		void Set_CommonVertex(TVector3& locVertex);
		void Set_CommonSpacetime(TLorentzVector& locSpacetime);

		void Set_CommonTParamIndex(int locCommonTParamIndex);
		void Set_OnlyConstrainTimeParticles(const set<DKinFitParticle*>& locOnlyConstrainTimeParticles){dOnlyConstrainTimeParticles = locOnlyConstrainTimeParticles;}

		set<DKinFitParticle*> dOnlyConstrainTimeParticles; //neutral showers //not used to constrain vertex, but fit vertex is used for time constraint

		double dInitTimeGuess;
};

inline DKinFitConstraint_Spacetime::DKinFitConstraint_Spacetime(void)
{
	Reset();
}

inline void DKinFitConstraint_Spacetime::Reset(void)
{
	DKinFitConstraint_Vertex::Reset();
	dInitTimeGuess = 0.0;
	dOnlyConstrainTimeParticles.clear();
}

inline set<DKinFitParticle*> DKinFitConstraint_Spacetime::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	set<DKinFitParticle*> locBaseParticles = DKinFitConstraint_Vertex::Get_AllParticles();
	set_union(locBaseParticles.begin(), locBaseParticles.end(), dOnlyConstrainTimeParticles.begin(), 
		dOnlyConstrainTimeParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

inline set<DKinFitParticle*> DKinFitConstraint_Spacetime::Get_AllConstrainedParticles(void) const
{
	set<DKinFitParticle*> locAllConstrainedParticles;
	set_union(dFullConstrainParticles.begin(), dFullConstrainParticles.end(), dOnlyConstrainTimeParticles.begin(), 
		dOnlyConstrainTimeParticles.end(), inserter(locAllConstrainedParticles, locAllConstrainedParticles.begin()));
	return locAllConstrainedParticles;
}

inline TLorentzVector DKinFitConstraint_Spacetime::Get_CommonSpacetime(void) const
{
	return TLorentzVector(Get_CommonVertex(), Get_CommonTime());
}

inline double DKinFitConstraint_Spacetime::Get_CommonTime(void) const
{
	if(dFullConstrainParticles.empty())
		return -1;
	return (*dFullConstrainParticles.begin())->Get_CommonVxParamIndex();
}

inline void DKinFitConstraint_Spacetime::Set_CommonTime(double locTime)
{
	set<DKinFitParticle*>::iterator locIterator = dFullConstrainParticles.begin();
	for(; locIterator != dFullConstrainParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonTime(locTime);
	for(locIterator = dNoConstrainParticles.begin(); locIterator != dNoConstrainParticles.end(); ++locIterator)
	{
		DKinFitParticleType locKinFitParticleType = (*locIterator)->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			(*locIterator)->Set_Time(locTime);
		else
			(*locIterator)->Set_CommonTime(locTime);
	}
	for(locIterator = dOnlyConstrainTimeParticles.begin(); locIterator != dOnlyConstrainTimeParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonTime(locTime);
}

inline void DKinFitConstraint_Spacetime::Set_CommonVertex(TVector3& locVertex)
{
	DKinFitConstraint_Vertex::Set_CommonVertex(locVertex);
	set<DKinFitParticle*>::iterator locIterator = dOnlyConstrainTimeParticles.begin();
	for(; locIterator != dOnlyConstrainTimeParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonVertex(locVertex);
}

inline void DKinFitConstraint_Spacetime::Set_CommonSpacetime(TLorentzVector& locSpacetime)
{
	DKinFitConstraint_Vertex::Set_CommonVertex(locSpacetime.Vect());
	Set_CommonTime(locSpacetime.T());
}

inline void DKinFitConstraint_Spacetime::Set_CommonTParamIndex(int locCommonTParamIndex)
{
	set<DKinFitParticle*>::iterator locIterator = dFullConstrainParticles.begin();
	for(; locIterator != dFullConstrainParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonTParamIndex(locCommonTParamIndex);
	for(locIterator = dOnlyConstrainTimeParticles.begin(); locIterator != dOnlyConstrainTimeParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonTParamIndex(locCommonTParamIndex);
	for(locIterator = dNoConstrainParticles.begin(); locIterator != dNoConstrainParticles.end(); ++locIterator)
	{
		DKinFitParticleType locKinFitParticleType = (*locIterator)->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			(*locIterator)->Set_TParamIndex(locCommonTParamIndex); //not included in fit, but particle vertex is defined by the fit result
		else
			(*locIterator)->Set_CommonTParamIndex(locCommonTParamIndex);
	}
}

inline int DKinFitConstraint_Spacetime::Get_CommonTParamIndex(void) const
{
	if(dFullConstrainParticles.empty())
		return -1;
	return (*dFullConstrainParticles.begin())->Get_CommonTParamIndex();
}

inline void DKinFitConstraint_Spacetime::Print_ConstraintInfo(void) const
{
	DKinFitConstraint_Vertex::Print_ConstraintInfo();

	cout << "DKinFitConstraint_Spacetime: Only-time-constrained particle PID's, pointers: " << endl;
	set<DKinFitParticle*>::const_iterator locIterator = dOnlyConstrainTimeParticles.begin();
	for(; locIterator != dOnlyConstrainTimeParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << (*locIterator) << endl;
}

#endif // _DKinFitConstraint_Spacetime_

