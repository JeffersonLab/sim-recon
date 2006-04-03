//********************************************************
// DFactory_DFDCPseudo.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//********************************************************

#include "DFactory_DFDCPseudo.h"

bool unsorted(vector<const DFDCHit*>& v) {
	for (unsigned int i=0; i < v.size(); i++)
		if (v[i+1]->globalPlane < v[i]->globalPlane)	
			return true;
	return false; 
}

bool DFDCHit_dE_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->dE < b->dE;
}

bool DFDCHit_layer_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->globalLayer < b->globalLayer;
}

bool DFDCHit_plane_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->globalPlane < b->globalPlane;
}


DFactory_DFDCPseudo::DFactory_DFDCPseudo() {}

DFactory_DFDCPseudo::~DFactory_DFDCPseudo() {}

derror_t DFactory_DFDCPseudo::evnt(DEventLoop* eventLoop, int eventNo) {
	vector<const DFDCHit*> fdcHits;
	vector<const DFDCHit*> uHits;
	vector<const DFDCHit*> vHits;
	map<const int, const DFDCHit*> anodeHits;
	float angle = 0.0;
	float pi	= 3.1415926;
	cout << "here" << endl;
	try {
		if (unsorted(fdcHits))
			std::sort(fdcHits.begin(), fdcHits.end(), DFDCHit_layer_cmp);
		
		vector<const DFDCHit*>::iterator iHit = fdcHits.begin(); 
		for (int iLayer = 1; iLayer <= 24; ++iLayer) {
			while ((iHit != fdcHits.end()) && (iLayer == (*iHit)->layer)) {
				const DFDCHit* aHit = *iHit;
				if (aHit->type) {
					if (aHit->plane == 1)
						uHits.push_back(aHit);	
					else
						vHits.push_back(aHit);
				}
				else
					anodeHits.insert(pair<const int, const DFDCHit*>(aHit->element, aHit));
			}
			
			std::stable_sort(uHits.begin(), uHits.end(), DFDCHit_dE_cmp);
			std::stable_sort(vHits.begin(), vHits.end(), DFDCHit_dE_cmp);
			
			conjure(uHits, vHits, anodeHits, angle);
			uHits.clear();
			vHits.clear();
			anodeHits.clear();
			angle += (pi/3);
		}	
	}
	catch (exception e) {
		cerr << e.what() << endl;
	}
	
	return NOERROR;
}

void DFactory_DFDCPseudo::conjure(	const vector<const DFDCHit*>& u, 
									const vector<const DFDCHit*>& v,
									const map<const int, const DFDCHit*>& x,
									float angle) {	
	float uYdist(0.0), vYdist(0.0);
	float uYint(0.0), vYint(0.0); 	
	float xCoordUV(0.0);
	int wireCandidateLow(0), wireCandidateHigh(0); 

	for (unsigned int i=0; i < u.size(); i++) {
		for (unsigned int j=0; j < v.size(); j++) {
			uYdist = 2*u[i]->r*u[i]->r;
			vYdist = 2*v[j]->r*v[j]->r;
			if (u[i]->r > 0)
				uYint = uYdist;
			else if (u[i]->r < 0)
				uYint = -uYdist;
			else
				uYint = 0.0;
			
			if (v[i]->r > 0)
				vYint = -vYdist;
			else if (v[i]->r < 0)
				vYint = vYdist;
			else
				vYint = 0.0;
				
			xCoordUV = (uYint - vYint) / 2;
			wireCandidateLow = static_cast<int>(floor(xCoordUV) + 60);
			wireCandidateHigh = static_cast<int>(ceil(xCoordUV) + 60);
			
			if ((x.count(wireCandidateLow)) || (x.count(wireCandidateHigh))) {
				float yCoordUV(0.0);
				float xCoordReal(0.0), yCoordReal(0.0), zCoordReal(0.0);
				yCoordUV = (vYint - uYint) / 2;
				xCoordReal = xCoordUV*cos(angle) - yCoordUV*sin(angle);
				yCoordReal = xCoordUV*sin(angle) + yCoordUV*sin(angle);
				zCoordReal = _geo.getLayerZ(u[i]);
				DFDCPseudo* newPseudo = new DFDCPseudo(xCoordReal, yCoordReal,
													 zCoordReal);
				// Add hits to DFDCPseudo::members?
				_data.push_back(newPseudo);
			}
		}
	}
}
			
