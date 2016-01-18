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
		int Get_CommonVxParamIndex(void) const;
		int Get_FIndex(DKinFitParticle* locKinFitParticle) const;

		virtual set<DKinFitParticle*> Get_AllConstrainedParticles(void) const{return dFullConstrainParticles;}
		set<DKinFitParticle*> Get_FullConstrainParticles(void) const{return dFullConstrainParticles;}
		set<DKinFitParticle*> Get_NoConstrainParticles(void) const{return dNoConstrainParticles;}
		virtual set<DKinFitParticle*> Get_AllParticles(void) const;

		void Print_ConstraintInfo(void) const;

	protected:

		DKinFitConstraint_Vertex(void);
		virtual ~DKinFitConstraint_Vertex(void){}

		virtual void Reset(void);

		void Set_FIndex(DKinFitParticle* locKinFitParticle, unsigned int locFIndex){dConstraintEquationParticleMap[locKinFitParticle] = locFIndex;}
		void Set_CommonVxParamIndex(int locCommonVxParamIndex);
		virtual void Set_CommonVertex(const TVector3& locVertex);

		void Set_FullConstrainParticles(const set<DKinFitParticle*>& locFullConstrainParticles){dFullConstrainParticles = locFullConstrainParticles;}
		void Set_NoConstrainParticles(const set<DKinFitParticle*>& locNoConstrainParticles){dNoConstrainParticles = locNoConstrainParticles;}

		set<DKinFitParticle*> dFullConstrainParticles; //charged particles, decaying particles, beam particles (not neutral showers!)
		set<DKinFitParticle*> dNoConstrainParticles; //missing particles, decaying particles, & neutral showers //fit vertex is set for these

		//key is particle, value is the constraint equation index
		map<DKinFitParticle*, unsigned int> dConstraintEquationParticleMap;

		TVector3 dInitVertexGuess;
};

DKinFitConstraint_Vertex::DKinFitConstraint_Vertex(void)
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

inline int DKinFitConstraint_Vertex::Get_FIndex(DKinFitParticle* locKinFitParticle) const
{
	map<DKinFitParticle*, unsigned int>::const_iterator locIterator = dConstraintEquationParticleMap.find(locKinFitParticle);
	if(locIterator == dConstraintEquationParticleMap.end())
		return -1;
	return locIterator->second;
}

inline set<DKinFitParticle*> DKinFitConstraint_Vertex::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	set_union(dFullConstrainParticles.begin(), dFullConstrainParticles.end(), dNoConstrainParticles.begin(), 
		dNoConstrainParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
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
	set<DKinFitParticle*>::iterator locIterator = dFullConstrainParticles.begin();
	for(; locIterator != dFullConstrainParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonVertex(locVertex);
	for(locIterator = dNoConstrainParticles.begin(); locIterator != dNoConstrainParticles.end(); ++locIterator)
	{
		DKinFitParticleType locKinFitParticleType = (*locIterator)->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			(*locIterator)->Set_Position(locVertex);
		else
			(*locIterator)->Set_CommonVertex(locVertex);
	}
}

inline int DKinFitConstraint_Vertex::Get_CommonVxParamIndex(void) const
{
	if(dFullConstrainParticles.empty())
		return -1;
	return (*dFullConstrainParticles.begin())->Get_CommonVxParamIndex();
}

inline void DKinFitConstraint_Vertex::Set_CommonVxParamIndex(int locCommonVxParamIndex)
{
	set<DKinFitParticle*>::iterator locIterator = dFullConstrainParticles.begin();
	for(; locIterator != dFullConstrainParticles.end(); ++locIterator)
		(*locIterator)->Set_CommonVxParamIndex(locCommonVxParamIndex);
	for(locIterator = dNoConstrainParticles.begin(); locIterator != dNoConstrainParticles.end(); ++locIterator)
	{
		DKinFitParticleType locKinFitParticleType = (*locIterator)->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			(*locIterator)->Set_VxParamIndex(locCommonVxParamIndex);
		else
			(*locIterator)->Set_CommonVxParamIndex(locCommonVxParamIndex);
	}
}

inline void DKinFitConstraint_Vertex::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_Vertex: Full-constrained particle PID's, q's, masses: " << endl;
	set<DKinFitParticle*>::const_iterator locIterator = dFullConstrainParticles.begin();
	for(; locIterator != dFullConstrainParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << (*locIterator)->Get_Charge() << ", " << (*locIterator)->Get_Mass() << endl;

	cout << "DKinFitConstraint_Vertex: No-constrain particle PID's, q's, masses: " << endl;
	for(locIterator = dNoConstrainParticles.begin(); locIterator != dNoConstrainParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << (*locIterator)->Get_Charge() << ", " << (*locIterator)->Get_Mass() << endl;
}

#endif // _DKinFitConstraint_Vertex_

