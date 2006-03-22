
#include "DFactory_DHDDMFDCTruth.h"
#include "DHDDMFDCTruth.h"



DFactory_DHDDMFDCTruth::DFactory_DHDDMFDCTruth() {}

DFactory_DHDDMFDCTruth::~DFactory_DHDDMFDCTruth() {}

derror_t DFactory_DHDDMFDCTruth::evnt(DEventLoop *eventLoop, int eventnumber)
{
	/// Place holder for now. 

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DHDDMFDCTruth::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.		
	
	v.clear();
	
	// Loop over Physics Events
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents) {
		throw DException("Attempt to get physics events from HDDM source failed.");
		return NOERROR;
	}
	
	for (unsigned int i=0; i < allEvents->mult; i++) {
		s_HitView_t *hits = allEvents->in[i].hitView;
		
		if (hits == HDDM_NULL) {
			throw DException("HDDM source has no hits.");
			return NOERROR;
		}
		
		if (hits->forwardDC == HDDM_NULL) {
			throw DException("HDDM source has no forwardDC information.");
			return NOERROR;
		}
		
		if (hits->forwardDC->fdcChambers == HDDM_NULL) {
			throw DException("HDDM source has no hits in the FDC.");		
			return NOERROR;
		}
		
		s_FdcChambers_t* fdcChamberSet = hits->forwardDC->fdcChambers;
		
		for (unsigned int i=0; i < fdcChamberSet->mult; i++) {
			s_FdcChamber_t fdcChamber = fdcChamberSet->in[i];
			s_FdcTruthPoints_t* fdcTruthPoints = fdcChamber.fdcTruthPoints;
			if (fdcTruthPoints == HDDM_NULL)
				continue;
			for (unsigned int j=0; j < fdcTruthPoints->mult; j++) {
				s_FdcTruthPoint_t truthPoint = fdcTruthPoints->in[j];
				DHDDMFDCTruth* newPoint = new DHDDMFDCTruth();
				newPoint->x				= truthPoint.x;
				newPoint->y				= truthPoint.y;
				newPoint->z				= truthPoint.z;
				newPoint->track			= truthPoint.track;
				newPoint->primary		= truthPoint.primary;
				newPoint->dradius		= truthPoint.dradius;
				newPoint->dEdx			= truthPoint.dEdx;
				v.push_back(newPoint);
			}
		}
	}		
	std::cout << "OKAY! " << v.size() << std::endl;
	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DHDDMFDCTruth::toString(void)
{
/*	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: layer: module: tau(rad):    z(cm):  u(cm):  dE(MeV):   t(ns):   type:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFDCHit *fdchit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d", fdchit->layer);
		printcol("%d", fdchit->module);
		printcol("%3.1f", fdchit->tau);
		printcol("%3.1f", fdchit->z);
		printcol("%2.3f", fdchit->u);
//		if(!fdchit->type){
			printcol("%1.3f", fdchit->dE*1000.0);
			printcol("%4.0f", fdchit->t);
		}else{
			printcol("");
			printcol("");
		}
//		printcol("%s", fdchit->type ? "cathode":"anode");
		printrow();
	}

	return _table;
*/
	return "";
}
