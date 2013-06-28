// $Id$
//
//    File: DTranslationTable.cc
// Created: Thu Jun 27 16:07:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#include "DTranslationTable.h"



using namespace jana;
using namespace std;

// Less than operator for csc_t data types. This is used by
// the map<csc_t, XX> to order the entires by key
bool operator<(const DTranslationTable::csc_t &a, const DTranslationTable::csc_t &b){
	if(a.rocid < b.rocid) return true;
	if(a.rocid > b.rocid) return false;
	if(a.slot < b.slot) return true;
	if(a.slot > b.slot) return false;
	if(a.channel < b.channel) return true;
	return false;
}

//---------------------------------
// DTranslationTable    (Constructor)
//---------------------------------
DTranslationTable::DTranslationTable(JEventLoop *loop)
{
	// These are used to help the event source report which
	// types of data it is capable of providing. For practical
	// purposes, these types are provided by the source
	// because they are generated and placed into their
	// respective JANA factories during a call to GetObjects().
	supplied_data_types.insert("DBCALHit");
	supplied_data_types.insert("DCDCHit");
	supplied_data_types.insert("DFCALHit");
	supplied_data_types.insert("DFDCHit");
	supplied_data_types.insert("DSCHit");
	supplied_data_types.insert("DTOFRawHit");

}

//---------------------------------
// ~DTranslationTable    (Destructor)
//---------------------------------
DTranslationTable::~DTranslationTable()
{

}

//---------------------------------
// IsSuppliedType
//---------------------------------
bool DTranslationTable::IsSuppliedType(string dataClassName) const
{
	return (supplied_data_types.find(dataClassName) != supplied_data_types.end());
}

//---------------------------------
// ApplyTranslationTable
//---------------------------------
void DTranslationTable::ApplyTranslationTable(JEventLoop *loop) const
{
	/// This will get all of the low level objects and
	/// generate detector hit objects from them, placing
	/// them in the appropriate DANA factories.
	
	vector<DBCALHit*> vbcal;
	vector<DCDCHit*> vcdc;
	vector<DFCALHit*> vfcal;
	vector<DFDCHit*> vfdc;
	vector<DSCHit*> vsc;
	vector<DTOFRawHit*> vtof;
	
	vector<const Df250PulseIntegral*> pulseintegrals;
	loop->Get(pulseintegrals);
	for(uint32_t i=0; i<pulseintegrals.size(); i++){
		const Df250PulseIntegral *pi = pulseintegrals[i];

		csc_t csc = {pi->rocid, pi->slot, pi->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()) continue;
		
		const DChannelInfo &chaninfo = iter->second;
		switch(chaninfo.det_sys){
			case BCAL        : vbcal.push_back( MakeBCALHit(pi, chaninfo.bcal) ); break;
			case CDC         : vcdc.push_back ( MakeCDCHit(pi, chaninfo.cdc) ); break;
			case FCAL        : vfcal.push_back( MakeFCALHit(pi, chaninfo.fcal) ); break;
			case FDC_CATHODES: vfdc.push_back ( MakeFDCCathodeHit(pi, chaninfo.fdc_cathodes) ); break;
			case FDC_WIRES   : vfdc.push_back ( MakeFDCCathodeHit(pi, chaninfo.fdc_wires) ); break;
			case ST          : vst.push_back  ( MakeSTHit(pi, chaninfo.st) ); break;
			case TOF         : vtof.push_back ( MakeTOFHit(pi, chaninfo.tof) ); break;

			default:  break;
		}
	}
	
	JFactory<DBCALHit>   *fac_bcal = dynamic_cast<JFactory<DBCALHit>*>(loop->GetFactory("DBCALHit"));
	JFactory<DCDCHit>    *fac_cdc  = dynamic_cast<JFactory<DCDCHit>*>(loop->GetFactory("DCDCHit"));
	JFactory<DFCALHit>   *fac_fcal = dynamic_cast<JFactory<DFCALHit>*>(loop->GetFactory("DFCALHit"));
	JFactory<DFDCHit>    *fac_fdc  = dynamic_cast<JFactory<DFDCHit>*>(loop->GetFactory("DFDCHit"));
	JFactory<DSCHit>     *fac_sc   = dynamic_cast<JFactory<DSCHit>*>(loop->GetFactory("DSCHit"));
	JFactory<DTOFRawHit> *fac_tof  = dynamic_cast<JFactory<DTOFRawHit>*>(loop->GetFactory("DTOFRawHit"));
	
	if(fac_bcal) fac_bcal->CopyTo(vbcal);
	if(fac_cdc ) fac_cdc->CopyTo(vcdc);
	if(fac_fcal) fac_fcal->CopyTo(vfcal);
	if(fac_fdc ) fac_fdc->CopyTo(vfdc);
	if(fac_sc  ) fac_sc->CopyTo(vsc);
	if(fac_tof ) fac_tof->CopyTo(vtof);
}

//---------------------------------
// MakeBCALHit
//---------------------------------
DBCALHit* DTranslationTable::MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const
{
	DBCALHit *h = new DBCALHit();
	h->module = idx.module;
	h->layer = idx.layer;
	h->sector = idx.sector;
	h->end = idx.end==0 ? DBCALGeometry::kUpstream:DBCALGeometry::kDownstream;
	
	h->E = hit->integral;
	h->t = 0;
	
	h->AddAssociatedObject(hit);
	return h;
}




