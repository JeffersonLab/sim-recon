//********************************************************
// DFDCPseudo_factory.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//********************************************************

#include "DFDCPseudo_factory.h"

///
/// DFDCCathodeCluster_gLayer_cmp(): 
/// non-member function passed to std::sort() to sort DFDCCathodeCluster pointers 
/// by their gLayer attributes.
///
bool DFDCCathodeCluster_gLayer_cmp(const DFDCCathodeCluster* a, 
								   const DFDCCathodeCluster* b) {
	return a->gLayer < b->gLayer;
}

///
/// DFDCPseudo_factory::DFDCPseudo_factory():
/// default constructor -- initializes log file
///
DFDCPseudo_factory::DFDCPseudo_factory() {
	logFile = new ofstream("DFDCPseudo_factory.log");
	_log = new JStreamLog(*logFile, "PSEUDO");
	*_log << "File initialized." << endMsg;
}

///
/// DFDCPseudo_factory::~DFDCPseudo_factory():
/// default destructor -- closes log file
///
DFDCPseudo_factory::~DFDCPseudo_factory() {
	logFile->close();
	delete logFile;
	delete _log;
}

///
/// DFDCPseudo_factory::evnt():
/// this is the place that anode hits and DFDCCathodeClusters are organized into pseudopoints.
/// For now, this is done purely by geometry, with no drift or peak-finding. See also
/// DFDCPseudo_factory::makePseudo().
///
jerror_t DFDCPseudo_factory::evnt(JEventLoop* eventLoop, int eventNo) {
	vector<const DFDCHit*> fdcHits;
	vector<const DFDCHit*> xHits;
	vector<const DFDCCathodeCluster*> cathClus;
	vector<const DFDCCathodeCluster*> uClus;
	vector<const DFDCCathodeCluster*> oneLayerU;
	vector<const DFDCCathodeCluster*> vClus;
	vector<const DFDCCathodeCluster*> oneLayerV;
	std::map<int, const DFDCHit*> oneLayerX;
	DFDCGeometry geo;
	float angle = 0.0;
	float pi	= 3.1415926;
	
	eventLoop->Get(fdcHits);
	eventLoop->Get(cathClus);
	
	// Ensure clusters are in order of ascending Z position.
	std::sort(cathClus.begin(), cathClus.end(), DFDCCathodeCluster_gLayer_cmp);
	
	// Sift through hits and select out anode hits.
	for (unsigned int i=0; i < fdcHits.size(); i++)
		if (fdcHits[i]->type == 0)
			xHits.push_back(fdcHits[i]);
			
	// Sift through clusters and put U and V clusters into respective vectors.
	for (unsigned int i=0; i < cathClus.size(); i++) {
		if (cathClus[i]->plane == 1)
			vClus.push_back(cathClus[i]);
		else
			uClus.push_back(cathClus[i]);
	}

	vector<const DFDCCathodeCluster*>::iterator uIt = uClus.begin();
	vector<const DFDCCathodeCluster*>::iterator vIt = vClus.begin();
	vector<const DFDCHit*>::iterator xIt = xHits.begin();
	
	// For each layer, get its sets of V, X, and U hits, and then pass them to the geometrical
	// organization routine, DFDCPseudo_factory::makePseudo()
	for (int iLayer=1; iLayer <= 24; iLayer++) {
		for (; ((uIt != uClus.end() && (*uIt)->gLayer == iLayer)); uIt++)
			oneLayerU.push_back(*uIt);
		for (; ((vIt != vClus.end() && (*vIt)->gLayer == iLayer)); vIt++)
			oneLayerV.push_back(*vIt);
		for (; ((xIt != xHits.end() && (*xIt)->gLayer == iLayer)); xIt++)
			oneLayerX.insert(std::pair<const int, const DFDCHit*>((*xIt)->element,
																  *xIt));
		makePseudo(oneLayerX, oneLayerU, oneLayerV, geo.getLayerRotation(iLayer), 	
					iLayer);
		angle += pi/3.0;
		oneLayerU.clear();
		oneLayerV.clear();
		oneLayerX.clear();
	}
	
	return NOERROR;
}

/// 
/// DFDCPseudo_factory::makePseudo():
/// performs UV+X matching to create pseudopoints
///
void DFDCPseudo_factory::makePseudo(	std::map<int, const DFDCHit*>& x,
										vector<const DFDCCathodeCluster*>& u,
										vector<const DFDCCathodeCluster*>& v,
										float angle,
										int layer)
{
	if ((u.size() == 0) || (v.size() == 0) || (x.size() == 0))
		return;
	std::map<int, const DFDCHit*>::iterator xIt;
	vector<const DFDCCathodeCluster*>::iterator uIt;
	vector<const DFDCCathodeCluster*>::iterator vIt;
	
	// Loop over all U and V clusters
	for (uIt = u.begin(); uIt != u.end(); uIt++) {
		for (vIt = v.begin(); vIt != v.end(); vIt++) {
			// Find the left and right edges of the cluster
			float leftX 	= intersectX((*uIt)->beginStrip, (*vIt)->beginStrip);
			float rightX 	= intersectX((*uIt)->endStrip, (*vIt)->endStrip);

			// Find the Y coordinate of the U and V strips with maxmimum dE
			float y			= intersectY((*uIt)->maxStrip, (*vIt)->maxStrip);
			float res 		= ((*uIt)->width + (*vIt)->width / 2) / 2;

			// Since wire numbers are integers, cast left and right bounds to integers.
			int left 		= static_cast<int>(floor(leftX + 60.0));
			int right 		= static_cast<int>(ceil(rightX + 60.0));
			
			// Ask the map if it has a hit for wire numbers between left and right;
			// if it does, create a new pseudopoint.			
			for (int i = left; i <= right; i++) 
			{
				xIt = x.find(i);
				if (xIt == x.end())
					continue;
				float xc = (*xIt).first - 60;
				xc = xc*cos(angle) - y*sin(angle);
				y = xc*sin(angle) + y*sin(angle);
				DFDCPseudo* newPseu = new DFDCPseudo(xc, y, layer, res);
				_data.push_back(newPseu);
			}
		}
	}
}			

///
/// DFDCPseudo_factory::intersectX():
/// finds the X coordinate of a U-V intersection
///
float DFDCPseudo_factory::intersectX(int u, int v) {	
	float uYint = (119 - u)*sqrt(2.0) + (1/sqrt(2.0));
	float vYint = (v - 119)*sqrt(2.0) - (1/sqrt(2.0)); 	

	return (uYint - vYint) / 2;
}			

///
/// DFDCPseudo_factory::intersectY():
/// finds the Y coordinate of a U-V intersection
///
float DFDCPseudo_factory::intersectY(int u, int v) {
	float uYint = (119 - u)*sqrt(2.0) + (1/sqrt(2.0));
	float vYint = (v - 119)*sqrt(2.0) - (1/sqrt(2.0));
	
	return (uYint + vYint) / 2;
}

