//******************************************************************************
// DFDCCathodeCluster_factory.cc
//
// Author:	Craig Bookwalter (craigb at jlab.org)
// Date:	April 2006
//
//******************************************************************************

#include "DFDCCathodeCluster_factory.h"

bool DFDCHit_cmp(const DFDCHit* a, const DFDCHit* b) {
  if (a->gLayer==b->gLayer){
    return a->t < b->t;
  }
  return a->gLayer < b->gLayer;
}


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
	if(a->element != b->element) return a->element < b->element;
	if(a->t       != b->t      ) return a->t < b->t;
	return a->q < b->q;
}

///
/// DFDCHit_time_cmp()
///    a non-member function passed to std::stable_sort() for sorting DFDCHit 
/// pointers in order of increasing time, provided that the time difference is
/// significant.
///

bool DFDCHit_time_cmp(const DFDCHit* a, const DFDCHit* b) {
  if (fabs(a->t-b->t)>HIT_TIME_DIFF_MIN && (a->t < b->t))
    return true;
  return false;
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
	_log = new JStreamLog(std::cout, "DFDCCathodeCluster >>");
}

///
/// DFDCCathodeCluster_factory::~DFDCCathodeCluster_factory():
/// default destructor--closes log file.
///
DFDCCathodeCluster_factory::~DFDCCathodeCluster_factory() {
	delete _log;
}

///
/// Initialization
///
jerror_t DFDCCathodeCluster_factory::init(void){
  TIME_SLICE=100.0; //ns,  Changed from 10->100 4/7/16 SJT
  gPARMS->SetDefaultParameter("FDC:CLUSTER_TIME_SLICE",TIME_SLICE);
  return NOERROR;	
}

///
/// DFDCCathodeCluster_factory::evnt():
/// This (along with DFDCCathodeCluster_factory::pique()) 
/// is the place cathode hits are associated into cathode clusters.  
///
jerror_t DFDCCathodeCluster_factory::evnt(JEventLoop *eventLoop, uint64_t eventNo) {
  vector<const DFDCHit*> allHits;
  vector<const DFDCHit*> uHits;
  vector<const DFDCHit*> vHits;
  vector<vector<const DFDCHit*> >thisLayer;
  
  try {
    eventLoop->Get(allHits);

    if (allHits.size()>0) {
      // Sort hits by layer number and by time
      sort(allHits.begin(),allHits.end(),DFDCHit_cmp);
      
      // Sift through all hits and select out U and V hits.
      for (vector<const DFDCHit*>::iterator i = allHits.begin(); 
	   i != allHits.end(); ++i){
	if ((*i)->type) {
	  if ((*i)->plane == 1)
	    vHits.push_back(*i);
	  else
	    uHits.push_back(*i);
	}
      }  
    	
      // Layer by layer, create clusters of U hits.
      if (uHits.size()>0){
	thisLayer.clear();
	vector<const DFDCHit*>::iterator i = uHits.begin();
	for (int iLayer=1;iLayer<25;iLayer++){
	  if (i==uHits.end()) break;
	  
	  vector<const DFDCHit*> hits;	
	  float old_time=(*i)->t;
	  while((i!=uHits.end()) && ((*i)->gLayer == iLayer)){ 
	    // Look for hits falling within a time slice
	    if (fabs((*i)->t-old_time)>TIME_SLICE){
	      // Sort hits by element number
	      sort(hits.begin(),hits.end(),DFDCHit_element_cmp);
	      // put into the vector
	      thisLayer.push_back(hits);
	      hits.clear();
	      old_time=(*i)->t;
	    }
	    hits.push_back(*i);
	    
	    i++;
	  }
	  // Sort hits by element number
	  sort(hits.begin(),hits.end(),DFDCHit_element_cmp);
	  // add the last vector of hits
	  thisLayer.push_back(hits);
	  
	  // Create clusters from these lists of hits
	  for (unsigned int k=0;k<thisLayer.size();k++) pique(thisLayer[k]);
	  
	  // Clear the hits and layer vectors for the next ones
	  thisLayer.clear();	
	  hits.clear();
	}
      }

      // Layer by layer, create clusters of V hits.
      if (vHits.size()>0){
	thisLayer.clear();
	vector<const DFDCHit*>::iterator i = vHits.begin();
	for (int iLayer=1;iLayer<25;iLayer++){
	  if (i==vHits.end()) break;
	  
	  vector<const DFDCHit*> hits;
	  float old_time=(*i)->t;
	  while((i!=vHits.end()) && ((*i)->gLayer == iLayer)){
	    // Look for hits falling within a time slice
	    if (fabs((*i)->t-old_time)>TIME_SLICE){
	      // Sort hits by element number
	      sort(hits.begin(),hits.end(),DFDCHit_element_cmp);
	      // put into the vector
	      thisLayer.push_back(hits);
	      hits.clear();
	      old_time=(*i)->t;
	    }
	    hits.push_back(*i);
	    
	    i++;
	  }
	  // Sort hits by element number
	  sort(hits.begin(),hits.end(),DFDCHit_element_cmp);
	  // add the last vector of hits
	  thisLayer.push_back(hits);
	  
	  // Create clusters from these lists of hits	  		
	  for (unsigned int k=0;k<thisLayer.size();k++) pique(thisLayer[k]);
	  
	  // Clear the hits and layer vectors for the next ones
	  thisLayer.clear();
	  hits.clear();
	}
      }
      
      // Ensure that the data are still in order of Z position.
      std::sort(_data.begin(), _data.end(), DFDCCathodeCluster_gPlane_cmp);
    }
  }
  catch (JException d) {
    cout << d << endl;
  }	
  catch (...) {
		cerr << "exception caught in DFDCCathodeCluster_factory" << endl;
  }
  
  return NOERROR;	
}			

//-----------------------------
// pique
//-----------------------------
void DFDCCathodeCluster_factory::pique(vector<const DFDCHit*>& H)
{
	/// Find clusters within cathode plane.
	///
	/// Upon entry, the vector "H" should already be sorted
	/// by strip number and should only contains hits from
	/// the same plane that are in time with each other.
	/// This will form clusters from all contiguous strips.

	// Loop over hits
	for(uint32_t istart=0; istart<H.size(); istart++){
		const DFDCHit *first_hit = H[istart];
		
		// Find end of contiguous section
		uint32_t iend=istart+1;
		for(; iend<H.size(); iend++){
			if(iend>=H.size()) break;
			if( (H[iend]->element - H[iend-1]->element) > 1 ) break;
		}
		if( (iend-istart)<2 ) continue; // don't allow single strip clusters
		
		// istart should now point to beginning of cluster 
		// and iend to one past end of cluster
		DFDCCathodeCluster* newCluster = new DFDCCathodeCluster();
		newCluster->q_tot   = 0.0;
		newCluster->gLayer  = first_hit->gLayer;
		newCluster->gPlane  = first_hit->gPlane;
		newCluster->plane   = first_hit->plane;
		for(uint32_t i=istart; i<iend; i++){
			newCluster->q_tot += H[i]->q;
			newCluster->members.push_back(H[i]);
		}
		_data.push_back(newCluster);
		
		istart = iend-1;
	}
}

