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
		DKinFitConstraint_Spacetime(void);
		~DKinFitConstraint_Spacetime(void){}

		double Get_InitTimeGuess(void) const{return dInitTimeGuess;};
		void Set_InitTimeGuess(double locInitTimeGuess){dInitTimeGuess = locInitTimeGuess;};

		TLorentzVector Get_CommonSpacetime(void) const;
		double Get_CommonTime(void) const;

		char Get_CommonTParamIndex(void) const;

		set<shared_ptr<DKinFitParticle>> Get_OnlyConstrainTimeParticles(void) const{return dOnlyConstrainTimeParticles;}
		set<shared_ptr<DKinFitParticle>> Get_AllConstrainingParticles(void) const;
		set<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const;

		void Print_ConstraintInfo(void) const;

		void Reset(void);
		void Release(void);

	private:

		void Set_CommonTime(double locTime);
		void Set_CommonVertex(const TVector3& locVertex);
		void Set_CommonSpacetime(TLorentzVector& locSpacetime);

		void Set_CommonTParamIndex(char locCommonTParamIndex);
		void Set_OnlyConstrainTimeParticles(const set<shared_ptr<DKinFitParticle>>& locOnlyConstrainTimeParticles){dOnlyConstrainTimeParticles = locOnlyConstrainTimeParticles;}

		set<shared_ptr<DKinFitParticle>> dOnlyConstrainTimeParticles; //neutral showers //not used to constrain vertex, but fit vertex is used for time constraint

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

inline void DKinFitConstraint_Spacetime::Release(void)
{
	DKinFitConstraint_Vertex::Release();
	dOnlyConstrainTimeParticles.clear();
}

inline set<shared_ptr<DKinFitParticle>> DKinFitConstraint_Spacetime::Get_AllParticles(void) const
{
	set<shared_ptr<DKinFitParticle>> locAllParticles;
	auto locBaseParticles = DKinFitConstraint_Vertex::Get_AllParticles();
	set_union(locBaseParticles.begin(), locBaseParticles.end(), dOnlyConstrainTimeParticles.begin(), dOnlyConstrainTimeParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

inline set<shared_ptr<DKinFitParticle>> DKinFitConstraint_Spacetime::Get_AllConstrainingParticles(void) const
{
	set<shared_ptr<DKinFitParticle>> locAllConstrainedParticles;
	set_union(dFullConstrainParticles.begin(), dFullConstrainParticles.end(), dOnlyConstrainTimeParticles.begin(), dOnlyConstrainTimeParticles.end(), inserter(locAllConstrainedParticles, locAllConstrainedParticles.begin()));
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
	for(auto& locParticle : dFullConstrainParticles)
		locParticle->Set_CommonTime(locTime);
	for(auto& locParticle : dNoConstrainParticles)
	{
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locParticle->Set_Time(locTime);
		else
			locParticle->Set_CommonTime(locTime);
	}
	for(auto& locParticle : dOnlyConstrainTimeParticles)
		locParticle->Set_CommonTime(locTime);
}

inline void DKinFitConstraint_Spacetime::Set_CommonVertex(const TVector3& locVertex)
{
	DKinFitConstraint_Vertex::Set_CommonVertex(locVertex);
	for(auto& locParticle : dOnlyConstrainTimeParticles)
		locParticle->Set_CommonVertex(locVertex);
}

inline void DKinFitConstraint_Spacetime::Set_CommonSpacetime(TLorentzVector& locSpacetime)
{
	DKinFitConstraint_Vertex::Set_CommonVertex(locSpacetime.Vect());
	Set_CommonTime(locSpacetime.T());
}

inline void DKinFitConstraint_Spacetime::Set_CommonTParamIndex(char locCommonTParamIndex)
{
	for(auto& locParticle : dFullConstrainParticles)
		locParticle->Set_CommonTParamIndex(locCommonTParamIndex);
	for(auto& locParticle : dOnlyConstrainTimeParticles)
		locParticle->Set_CommonTParamIndex(locCommonTParamIndex);
	for(auto& locParticle : dNoConstrainParticles)
	{
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locParticle->Set_TParamIndex(locCommonTParamIndex); //not included in fit, but particle vertex is defined by the fit result
		else
			locParticle->Set_CommonTParamIndex(locCommonTParamIndex);
	}
}

inline char DKinFitConstraint_Spacetime::Get_CommonTParamIndex(void) const
{
	if(dFullConstrainParticles.empty())
		return -1;
	return (*dFullConstrainParticles.begin())->Get_CommonTParamIndex();
}

inline void DKinFitConstraint_Spacetime::Print_ConstraintInfo(void) const
{
	DKinFitConstraint_Vertex::Print_ConstraintInfo();

	cout << "DKinFitConstraint_Spacetime: Only-time-constrained particle PID's, pointers: " << endl;
	for(auto& locParticle : dOnlyConstrainTimeParticles)
		cout << locParticle->Get_PID() << ", " << locParticle << endl;
}

#endif // _DKinFitConstraint_Spacetime_
