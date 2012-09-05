// $Id$
//
//    File: DNeutralShower.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralShower_
#define _DNeutralShower_

#include <vector>
#include <utility>
#include <string>

#include "JANA/JObject.h"

#include "GlueX.h"
#include "DLorentzVector.h"
#include "DMatrixDSym.h"

#include "BCAL/DBCALShower.h"
#include "FCAL/DFCALShower.h"

using namespace std;

class DNeutralShower : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DNeutralShower);

		DLorentzVector dSpacetimeVertex;
		float dEnergy;
		DetectorSystem_t dDetectorSystem;
		DMatrixDSym dCovarianceMatrix; //E, x, y, z, t

		DNeutralShower(const DBCALShower *locBCALShower)
		{
			dDetectorSystem = SYS_BCAL;
			dEnergy = locBCALShower->E;
			double locEnergyUncertainty = (dEnergy >= 0.0) ? dEnergy*sqrt( 0.0598*0.0598/dEnergy + 0.0094*0.0094 ) : 1e-3; //last updated at svn revision 9242
			dSpacetimeVertex.SetXYZT(locBCALShower->x, locBCALShower->y, locBCALShower->z, locBCALShower->t);
			dCovarianceMatrix.ResizeTo(5, 5);
			dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
			dCovarianceMatrix(1, 1) = locBCALShower->xErr*locBCALShower->xErr;
			dCovarianceMatrix(2, 2) = locBCALShower->xErr*locBCALShower->yErr;
			dCovarianceMatrix(3, 3) = locBCALShower->xErr*locBCALShower->zErr;
			dCovarianceMatrix(4, 4) = locBCALShower->xErr*locBCALShower->tErr;
			//NEED CORRELATIONS!
		}

		DNeutralShower(const DFCALShower *locFCALShower)
		{
			dDetectorSystem = SYS_FCAL;
			dEnergy = locFCALShower->getEnergy();
			double locEnergyUncertainty = (dEnergy >= 0.0) ? 0.042*sqrt(dEnergy) + 0.0001 : 1e-3; //from old DPhoton_factory::makeFCalPhoton() function
			dSpacetimeVertex.SetVect(locFCALShower->getPosition());
			dSpacetimeVertex.SetT(locFCALShower->getTime());
			dCovarianceMatrix.ResizeTo(5, 5);
			dCovarianceMatrix(0, 0) = locEnergyUncertainty*locEnergyUncertainty;
			dCovarianceMatrix(1, 1) = locFCALShower->getPositionError().X()*locFCALShower->getPositionError().X();
			dCovarianceMatrix(2, 2) = locFCALShower->getPositionError().Y()*locFCALShower->getPositionError().Y();
			dCovarianceMatrix(3, 3) = locFCALShower->getPositionError().Z()*locFCALShower->getPositionError().Z();
			dCovarianceMatrix(4, 4) = 0.0; //not stored in DFCALShower
			//NEED CORRELATIONS!
		}

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "E", "%3.5f", dEnergy);
			AddString(items, "x", "%3.2f", dSpacetimeVertex.X());
			AddString(items, "y", "%3.2f", dSpacetimeVertex.Y());
			AddString(items, "z", "%3.2f", dSpacetimeVertex.Z());
			AddString(items, "t", "%3.2f", dSpacetimeVertex.T());
		}

};

#endif // _DNeutralShower_

