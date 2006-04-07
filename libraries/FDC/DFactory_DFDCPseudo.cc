//********************************************************
// DFactory_DFDCPseudo.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//********************************************************

#include "DFactory_DFDCPseudo.h"

void DFactory_DFDCPseudo::crap(	vector<const DFDCHit*>& u, 
		  						vector<const DFDCHit*>& v,
								map<const int, const DFDCHit*>& x,
								float angle,
								int evNo,
								int layerNo) {	
	float uYdist(0.0), vYdist(0.0);
	float uYint(0.0), vYint(0.0); 	
	float xCoordUV(0.0);
	int wireCandidateLow(0), wireCandidateHigh(0);
	ofstream UVofs("uv.dat", ios::app);
	ofstream Xofs("x.dat", ios::app);
	*_log << "In crap, layer: " << layerNo << " u: " << u.size() << " v: " 
		  << v.size() << " x: " << x.size() << endMsg;	
	for (vector<const DFDCHit*>::iterator i = u.begin(); i != u.end(); ++i) {
		for (vector<const DFDCHit*>::iterator j = v.begin(); j != v.end(); ++j) {
			uYdist = sqrt(2)*(*i)->r;
			vYdist = sqrt(2)*(*j)->r;
			
			if ((*i)->r > 0)
				uYint = uYdist;
			else if ((*i)->r < 0)
				uYint = -uYdist;
			else
				uYint = 0.0;
			
			if ((*j)->r > 0)
				vYint = -vYdist;
			else if ((*j)->r < 0)
				vYint = vYdist;
			else
				vYint = 0.0;
			xCoordUV = (uYint - vYint) / 2;
			wireCandidateLow = static_cast<int>(floor(xCoordUV) + 60);
			wireCandidateHigh = static_cast<int>(ceil(xCoordUV) + 60);
			
			UVofs << "(" << wireCandidateLow << "," << wireCandidateHigh << ")" << " ";
		}
	}
	for (map<int, const DFDCHit*>::iterator i = x.begin(); i != x.end(); ++i) {
		Xofs << (*i).second->r << "\t";
	}
	
	UVofs << endl << endl;
	Xofs << endl << endl;
	UVofs.close();
	Xofs.close();
}


bool unsorted(vector<const DFDCHit*>& v) {
	for (unsigned int i=0; i < v.size(); i++)
		if (v[i+1]->gPlane < v[i]->gPlane)	
			return true;
	return false; 
}

bool DFDCHit_dE_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->dE < b->dE;
}

bool DFDCHit_layer_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->gLayer < b->gLayer;
}


DFactory_DFDCPseudo::DFactory_DFDCPseudo() {
	logFile = new ofstream("DFactory_DFDCPseudo.log");
	_log = new DStreamLog(*logFile, "PSEUDO");
	*_log << "File initialized." << endMsg;
}

DFactory_DFDCPseudo::~DFactory_DFDCPseudo() {
	logFile->close();
	delete logFile;
	delete _log;
}

derror_t DFactory_DFDCPseudo::evnt(DEventLoop* eventLoop, int eventNo) {
	vector<const DFDCHit*> fdcHits;
	vector<const DFDCHit*> uHits;
	vector<const DFDCHit*> vHits;
	map<const int, const DFDCHit*> anodeHits;
	float angle = 0.0;
	float pi	= 3.1415926;
	unsigned int uHitsTot(0), vHitsTot(0), xHitsTot(0);
	
	eventLoop->Get(fdcHits);
	*_log << "fdcHits.size() = " << fdcHits.size() << endMsg;
	try {
		if (unsorted(fdcHits))
			std::sort(fdcHits.begin(), fdcHits.end(), DFDCHit_layer_cmp);
		
		vector<const DFDCHit*>::iterator iHit = fdcHits.begin(); 
		for (int iLayer = 1; iLayer <= 24; ++iLayer) {
			while ((iHit != fdcHits.end()) && (iLayer == (*iHit)->gLayer)) {
				const DFDCHit* aHit = *iHit;
				if (aHit->type) {
					if (aHit->plane == 1)
						uHits.push_back(aHit);	
					else
						vHits.push_back(aHit);
				}
				else
					anodeHits.insert(pair<const int, const DFDCHit*>(aHit->element, aHit));
				iHit++;
				
			}
			uHitsTot += uHits.size();
			vHitsTot += vHits.size();
			xHitsTot += anodeHits.size();
				  
			std::stable_sort(uHits.begin(), uHits.end(), DFDCHit_dE_cmp);
			std::stable_sort(vHits.begin(), vHits.end(), DFDCHit_dE_cmp);
//			conjure(uHits, vHits, anodeHits, angle);
			crap(uHits,vHits,anodeHits,angle,eventNo,iLayer);
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

void DFactory_DFDCPseudo::conjure(	vector<const DFDCHit*>& u, 
									vector<const DFDCHit*>& v,
									map<const int, const DFDCHit*>& x,
									float angle) {	
									
	// TODO: Use the pulse thingy to associate hits. For real.
	
	float uYdist(0.0), vYdist(0.0);
	float uYint(0.0), vYint(0.0); 	
	float xCoordUV(0.0);
	int wireCandidateLow(0), wireCandidateHigh(0);
	for (vector<const DFDCHit*>::iterator i = u.begin(); i != u.end(); ++i) {
		for (vector<const DFDCHit*>::iterator j = v.begin(); j != v.end(); ++j) {
			uYdist = sqrt(2)*(*i)->r;
			vYdist = sqrt(2)*(*j)->r;
			
			if ((*i)->r > 0)
				uYint = uYdist;
			else if ((*i)->r < 0)
				uYint = -uYdist;
			else
				uYint = 0.0;
			
			if ((*j)->r > 0)
				vYint = -vYdist;
			else if ((*j)->r < 0)
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
				DFDCPseudo* newPseudo = new DFDCPseudo(xCoordReal, yCoordReal,
													 (*i)->gLayer);
				// Add hits to DFDCPseudo::members?
				_data.push_back(newPseudo);
				
			}
		}
	}
}
			
