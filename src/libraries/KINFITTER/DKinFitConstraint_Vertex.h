#ifndef _DKinFitConstraints_Vertex_
#define _DKinFitConstraints_Vertex_

#include <set>
#include <algorithm>

#include "TVector3.h"

#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"

using namespace std;

class DKinFitter;
class DKinFitUtils;

class DKinFitConstraint_Vertex : public DKinFitConstraint
{
	friend class DKinFitter;
	friend class DKinFitUtils;

	public:

		TVector3 Get_InitVertexGuess(void) const{return dInitVertexGuess;};
		void Set_InitVertexGuess(const TVector3& locInitVertexGuess){dInitVertexGuess = locInitVertexGuess;};

		TVector3 Get_CommonVertex(void) const;
		char Get_CommonVxParamIndex(void) const;
		char Get_FIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle) const;

		virtual set<shared_ptr<DKinFitParticle>> Get_AllConstrainingParticles(void) const{return dFullConstrainParticles;}
		set<shared_ptr<DKinFitParticle>> Get_FullConstrainParticles(void) const{return dFullConstrainParticles;}
		set<shared_ptr<DKinFitParticle>> Get_NoConstrainParticles(void) const{return dNoConstrainParticles;}
		virtual set<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const;

		void Print_ConstraintInfo(void) const;

	protected:

		DKinFitConstraint_Vertex(void);
		virtual ~DKinFitConstraint_Vertex(void){}

		virtual void Reset(void);

		void Set_FIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle, char locFIndex){dConstraintEquationParticleMap[locKinFitParticle] = locFIndex;}
		void Set_CommonVxParamIndex(char locCommonVxParamIndex);
		virtual void Set_CommonVertex(const TVector3& locVertex);

		void Set_FullConstrainParticles(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles){dFullConstrainParticles = locFullConstrainParticles;}
		void Set_NoConstrainParticles(const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles){dNoConstrainParticles = locNoConstrainParticles;}

		set<shared_ptr<DKinFitParticle>> dFullConstrainParticles; //charged particles, decaying particles, beam particles (not neutral showers!)
		set<shared_ptr<DKinFitParticle>> dNoConstrainParticles; //missing particles, decaying particles, & neutral showers //fit vertex is set for these

		//key is particle, value is the constraint equation index
		map<shared_ptr<DKinFitParticle>, char> dConstraintEquationParticleMap;

		TVector3 dInitVertexGuess;
};

inline DKinFitConstraint_Vertex::DKinFitConstraint_Vertex(void)
{
	Reset();
}

inline void DKinFitConstraint_Vertex::Reset(void)
{
	dInitVertexGuess = TVector3(0.0, 0.0, 0.0);
	dFullConstrainParticles.clear();
	dNoConstrainParticles.clear();
	dConstraintEquationParticleMap.clear();
}

inline char DKinFitConstraint_Vertex::Get_FIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle) const
{
	auto locIterator = dConstraintEquationParticleMap.find(locKinFitParticle);
	if(locIterator == dConstraintEquationParticleMap.end())
		return -1;
	return locIterator->second;
}

inline set<shared_ptr<DKinFitParticle>> DKinFitConstraint_Vertex::Get_AllParticles(void) const
{
	set<shared_ptr<DKinFitParticle>> locAllParticles;
	set_union(dFullConstrainParticles.begin(), dFullConstrainParticles.end(), dNoConstrainParticles.begin(), dNoConstrainParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

inline TVector3 DKinFitConstraint_Vertex::Get_CommonVertex(void) const
{
	if(dFullConstrainParticles.empty())
		return TVector3();
	return (*dFullConstrainParticles.begin())->Get_CommonVertex();
}

inline void DKinFitConstraint_Vertex::Set_CommonVertex(const TVector3& locVertex)
{
	for(auto& locParticle : dFullConstrainParticles)
		locParticle->Set_CommonVertex(locVertex);
	for(auto& locParticle : dNoConstrainParticles)
	{
		auto locKinFitParticleType = locParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locParticle->Set_Position(locVertex);
		else
			locParticle->Set_CommonVertex(locVertex);
	}
}

inline char DKinFitConstraint_Vertex::Get_CommonVxParamIndex(void) const
{
	if(dFullConstrainParticles.empty())
		return -1;
	return (*dFullConstrainParticles.begin())->Get_CommonVxParamIndex();
}

inline void DKinFitConstraint_Vertex::Set_CommonVxParamIndex(char locCommonVxParamIndex)
{
	for(auto& locParticle : dFullConstrainParticles)
		locParticle->Set_CommonVxParamIndex(locCommonVxParamIndex);
	for(auto& locParticle : dNoConstrainParticles)
	{
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locParticle->Set_VxParamIndex(locCommonVxParamIndex);
		else
			locParticle->Set_CommonVxParamIndex(locCommonVxParamIndex);
	}
}

inline void DKinFitConstraint_Vertex::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_Vertex: Full-constrained particle PID's, pointers: " << endl;
	for(auto& locParticle : dFullConstrainParticles)
		cout << locParticle->Get_PID() << ", " << locParticle << endl;

	cout << "DKinFitConstraint_Vertex: No-constrain particle PID's, pointers: " << endl;
	for(auto& locParticle : dNoConstrainParticles)
		cout << locParticle->Get_PID() << ", " << locParticle << endl;
}

#endif // _DKinFitConstraint_Vertex_

