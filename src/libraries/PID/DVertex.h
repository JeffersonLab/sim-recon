// $Id$
//
//    File: DVertex.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_
#define _DVertex_

#include <map>

#include <JANA/JObject.h>

#include <DLorentzVector.h>

#include "PID/DKinematicData.h"
#include "KINFITTER/DKinFitParticle.h"

using namespace std;
using namespace jana;

class DVertex: public jana::JObject
{
	public:    
		JOBJECT_PUBLIC(DVertex);

		DLorentzVector dSpacetimeVertex; // vertex position in cm + vertex time in ns
		TMatrixFSym dCovarianceMatrix; //xyzt order

		unsigned int dKinFitNDF;
		double dKinFitChiSq;
		map<const JObject*, map<DKinFitPullType, double> > dKinFitPulls;

		// Objects used to calculate this added as Associated Objects
		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "x", "%3.2f", dSpacetimeVertex.X());
			AddString(items, "y", "%3.2f", dSpacetimeVertex.Y());
			AddString(items, "z", "%3.2f", dSpacetimeVertex.Z());
			AddString(items, "t", "%3.2f", dSpacetimeVertex.T());
		}
};

#endif // _DVertex_

