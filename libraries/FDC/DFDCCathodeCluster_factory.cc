//******************************************************************************
// DFDCCathodeCluster_factory.cc
//
// Author:	Craig Bookwalter (craigb at jlab.org)
// Date:	April 2006
//
//******************************************************************************

#include "DFDCCathodeCluster_factory.h"

///
///	DFDCHit_gLayer_cmp(): 
/// a non-member function passed to std::sort() for sorting DFDCHit pointers by
/// their gLayer attributes.
///
bool DFDCHit_gLayer_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->gLayer < b->gLayer;
}

///
/// DFDCHit_element_cmp():
///	a non-member function passed to std::sort() for sorting DFDCHit pointers by
/// their element (wire or strip) numbers. Typically only used for a single layer
/// of hits.
///
bool DFDCHit_element_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->element < b->element;
}

///
/// DFDCCathodeCluster_gPlane_cmp():
/// a non-member function passed to std::sort() for sorting DFDCCathodeCluster pointers
/// by their gPlane (plane number over all modules, 1-74) attributes.
///
bool DFDCCathodeCluster_gPlane_cmp(	const DFDCCathodeCluster* a, 
									const DFDCCathodeCluster* b) {
	return a->gPlane < b->gPlane;
}

///
/// DFDCCathodeCluster_factory::DFDCCathodeCluster_factory():
///	default constructor--initializes log file
///
DFDCCathodeCluster_factory::DFDCCathodeCluster_factory() {
	_logFile = new ofstream("DFDCCathodeCluster_factory.log");
	_log = new JStreamLog(*_logFile, "CLUST");
}

///
/// DFDCCathodeCluster_factory::~DFDCCathodeCluster_factory():
/// default destructor--closes log file.
///
DFDCCathodeCluster_factory::~DFDCCathodeCluster_factory() {
	_logFile->close();
	delete _logFile;
	delete _log;
}

///
/// DFDCCathodeCluster_factory::evnt():
/// This (along with DFDCCathodeCluster_factory::pique()) 
/// is the place cathode hits are associated into cathode clusters. This function 
/// should eventually be modified to do more sophisticated peak finding. 
///
jerror_t DFDCCathodeCluster_factory::evnt(JEventLoop *eventLoop, int eventNo) {
	vector<const DFDCHit*> allHits;
	vector<const DFDCHit*> uHits;
	vector<const DFDCHit*> vHits;
	vector<const DFDCHit*> thisLayer;
	//vector<const DFDCHit*>::iterator planeBegin;
	//vector<const DFDCHit*>::iterator planeEnd;
	int iLayer = 1;
	
	try {
		eventLoop->Get(allHits);

		// Sift through all hits and select out U and V hits.
		for (vector<const DFDCHit*>::iterator i = allHits.begin(); 
			 i != allHits.end(); ++i)
		{
			if ((*i)->type) {
				if ((*i)->plane == 1)
					vHits.push_back(*i);
				else
					uHits.push_back(*i);
			}
		}
		
		// Ensure all cathode hits are in order of increasing Z position.
		std::sort(uHits.begin(), uHits.end(), DFDCHit_gLayer_cmp);
		std::sort(vHits.begin(), vHits.end(), DFDCHit_gLayer_cmp);
		thisLayer.clear();

		// Layer by layer, create clusters of U hits. 		
		for (vector<const DFDCHit*>::iterator i = uHits.begin(); 
			 i != uHits.end(); ++i)
		{
			if ((*i)->gLayer == iLayer)		
				thisLayer.push_back(*i);
			else {
				++iLayer;
				pique(thisLayer);
				thisLayer.clear();
			}
		}
		
		// Set up for V layers.				
		iLayer = 1;
		thisLayer.clear();

		// Layer by layer, create clusters of V hits.		
		for (vector<const DFDCHit*>::iterator i = vHits.begin(); 
			 i != vHits.end(); ++i)
		{
			if ((*i)->gLayer == iLayer)		
				thisLayer.push_back(*i);
			else {
				++iLayer;
				pique(thisLayer);
				thisLayer.clear();
			}
		}
		
		// Ensure that the data are still in order of Z position.
		std::sort(_data.begin(), _data.end(), DFDCCathodeCluster_gPlane_cmp);
	}
	catch (JException d) {
		cout << d << endl;
	}	
	catch (...) {
		cerr << "exception caught in DFDCCathodeCluster_factory" << endl;
	}
	
	return NOERROR;	
}			

///
/// DFDCCathodeCluster_factory::pique():
/// takes a single layer's worth of cathode hits and attempts to create DFDCCathodeClusters
/// by grouping together hits with consecutive strip numbers.
///
void DFDCCathodeCluster_factory::pique(vector<const DFDCHit*>& H) {
	if (H.size() == 0)
		return;
		
	int width(1);
	int beginStrip(0);
	int maxStrip(0);
	float dEtot(0.0);
	float dEmax(0.0);
	
	// Ensure the hits are in ascending strip number order
	std::sort(H.begin(), H.end(), DFDCHit_element_cmp);
	
	beginStrip = (*(H.begin()))->element;

	// For all hits in this layer, associate consecutively-numbered strips into a 
	// DFDCCathodeCluster object. 
	for (vector<const DFDCHit*>::iterator i = H.begin(); i != H.end(); i++) {
		// If we're not at the end of the array, and the strip number of the next hit is
		// equal to the strip number + 1 of this hit, then we continue our cluster.
		if ((i+1 != H.end()) && ((*i)->element + 1 == (*(i+1))->element)) {
			width++;
			dEtot += (*i)->dE;
			if ((*i)->dE > dEmax) {
				dEmax = (*i)->dE;
				maxStrip = (*i)->element;
			}
		}
		// If not, our cluster must have ended, so we record the information into a new
		// DFDCCathodeCluster object and reset for the next cluster.
		else {
			DFDCCathodeCluster* newCluster = new DFDCCathodeCluster();
			if (width > 1) {
				newCluster->beginStrip  = beginStrip;
				newCluster->maxStrip 	= maxStrip;
				newCluster->dEtot 		= dEtot;
			}
			else {
				newCluster->beginStrip  = (*i)->element;
				newCluster->maxStrip	= (*i)->element;
				newCluster->dEtot		= (*i)->dE;
			}
			newCluster->width 		= width;
			newCluster->endStrip	= (*i)->element;
			newCluster->gLayer		= (*i)->gLayer;
			newCluster->gPlane		= (*i)->gPlane;
			newCluster->plane		= (*i)->plane;
			_data.push_back(newCluster);
			width 		= 1;
			maxStrip 	= 0;
			dEtot 		= 0.0;
			dEmax		= 0.0;
			if (i+1 != H.end())
				beginStrip  = (*(i+1))->element;
		}
	}
}

///
/// DFDCCathodeCluster_factory::toString():
/// returns a sensible std::string representation of the data contained in this 
/// factory.
///
const string DFDCCathodeCluster_factory::toString() {
	if (_data.size() <= 0)
		return "";
		
	stringstream s;
	s << (*_data.begin())->header() << endl;

	// Simply call each cluster's toString() method and stream it into s.
	for (unsigned int i=0; i < _data.size(); ++i) 
		s << _data[i]->toString() << endl;
	
	return s.str();
}
