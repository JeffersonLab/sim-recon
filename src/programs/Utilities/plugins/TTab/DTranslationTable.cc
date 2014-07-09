// $Id$
//
//    File: DTranslationTable.cc
// Created: Thu Jun 27 16:07:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#include "DTranslationTable.h"

#include <expat.h>
#include <sstream>

#include <DAQ/DModuleType.h>

using namespace jana;
using namespace std;

// Use one translation table for all threads
static pthread_mutex_t tt_mutex = PTHREAD_MUTEX_INITIALIZER;
static bool tt_initialized = false;
static map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo> TT;
string ROCID_MAP_FILENAME;
static map<uint32_t, uint32_t> rocid_map;     // (see ReadOptionalROCidTranslation() for details)
static map<uint32_t, uint32_t> rocid_inv_map; // (see ReadOptionalROCidTranslation() for details)


//...................................
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

//...................................
// sort functions
bool SortBCALDigiHit(const DBCALDigiHit *a, const DBCALDigiHit *b){
	if(a->module == b->module){
		if(a->layer == b->layer){
			if(a->sector == b->sector){
				if(a->end == b->end){
					return a->pulse_time < b->pulse_time;
				}else{ return a->end < b->end; }
			}else{ return a->sector < b->sector; }
		}else{ return a->layer< b->layer; }
	}else { return a->module < b->module; }
}

//---------------------------------
// DTranslationTable    (Constructor)
//---------------------------------
DTranslationTable::DTranslationTable(JEventLoop *loop)
{
	// Default is to just read translation table from CCDB. If this fails,
	// then an attempt will be made to read from a file on the local disk.
	// The filename can be specified to be anything, but if the user specifies
	// this, then we assume that they want to use it and skip using the CCDB.
	// They may also specify that they want to skip checking the CCDB via
	// the "TT:NO_CCDB" parameter. This would only be useful if they want to
	// force the use of a local file named "tt.xml".
	NO_CCDB = false;
	XML_FILENAME = "tt.xml";
	VERBOSE = 0;
	gPARMS->SetDefaultParameter("TT:NO_CCDB", NO_CCDB, "Don't try getting translation table from CCDB and just look for file. Only useful if you want to force reading tt.xml. This is automatically set if you specify a different filename via the TT:XML_FILENAME parameter.");
	JParameter *p = gPARMS->SetDefaultParameter("TT:XML_FILENAME", XML_FILENAME, "Fallback filename of translation table XML file. If set to non-default, CCDB will not be checked.");
	if(p->GetDefault() != p->GetValue()) NO_CCDB = true;
	gPARMS->SetDefaultParameter("TT:VERBOSE", VERBOSE, "Verbosity level for Applying Translation Table. 0=no messages, 10=all messages.");
	
	ROCID_MAP_FILENAME = "rocid.map";
	gPARMS->SetDefaultParameter("TT:ROCID_MAP_FILENAME", ROCID_MAP_FILENAME, "Optional rocid to rocid conversion map for use with files generated with the non-standard rocid's");

	// Initialize dedicated JStreamLog used for debugging messages
	ttout.SetTag("--- TT ---: ");
	ttout.SetTimestampFlag();
	ttout.SetThreadstampFlag();

	// Look for and read in an optional rocid <-> rocid translation table
	ReadOptionalROCidTranslation();

	// Read in Translation table. This will create DChannelInfo objects
	// and store them in the "TT" map, indexed by csc_t objects
	ReadTranslationTable(loop->GetJCalibration());

	// These are used to help the event source report which
	// types of data it is capable of providing. For practical
	// purposes, these types are "provided" by the source
	// because they are generated and placed into their
	// respective JANA factories during a call to GetObjects().
	// The source is responsible for reporting the types it is
	// directly responsible for (e.g. Df250PulseIntegral)
	supplied_data_types.insert("DBCALDigiHit");
	supplied_data_types.insert("DBCALTDCDigiHit");
	supplied_data_types.insert("DCDCDigiHit");
	supplied_data_types.insert("DFCALDigiHit");
	supplied_data_types.insert("DFDCCathodeDigiHit");
	supplied_data_types.insert("DFDCWireDigiHit");
	supplied_data_types.insert("DSCDigiHit");
	supplied_data_types.insert("DSCTDCDigiHit");
	supplied_data_types.insert("DTOFDigiHit");
	supplied_data_types.insert("DTOFTDCDigiHit");
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
// ReadOptionalROCidTranslation
//---------------------------------
void DTranslationTable::ReadOptionalROCidTranslation(void)
{
	// Some data may be taken with the ROC ID value set
	// incorrectly in CODA. For CODA 3.0 data, there is
	// actually no way to set it so it can be different
	// for every CODA configuration. A simple work-around
	// for this is to use a map file to list the translation
	// from the crate numbers used in the evio file to those
	// stored in the TT. Check here if a local file exists
	// with the name specified by the TT:ROCID_MAP_FILENAME
	// config parameter (default is "rocid.map"). If so,
	// read it in. The format is just 2 values per line.
	// The first is the rocid in the evio file, and the
	// second, what the rocid is in the TT. Note that the
	// value of the crate copied into the data objects 
	// will be what is in the EVIO file.
	ifstream ifs(ROCID_MAP_FILENAME.c_str());
	if(!ifs.is_open()) return;
	
	cout << "Opened ROC id translation map: " << ROCID_MAP_FILENAME << endl;
	while(ifs.good()){
		char line[256];
		ifs.getline(line, 256);
		if(ifs.gcount() < 1) break;
		if(line[0] == '#') continue;

		stringstream ss(line);
		uint32_t from=10000, to=10000;
		ss >> from >> to;  // from=evio  to=TT
		if( to==10000 ){
			if( from!=10000){
				cout << "unable to convert line:" << endl;
				cout << "  " << line;
			}
		}else{
			rocid_map[from] = to;
			rocid_inv_map[to] = from;
		}
	}
	ifs.close();
	
	if(rocid_map.size() == rocid_inv_map.size()){
		cout << "   Read " << rocid_map.size() << " entries" << endl;
		map<uint32_t,uint32_t>::iterator iter;
		for(iter=rocid_map.begin(); iter != rocid_map.end(); iter++){
			cout << "   rocid " << iter->first << " -> rocid " << iter->second << endl;
		}
	}else{
		cout << "Entries not unique! This can happen if there are" <<endl;
		cout << "more than one entry with the same value (either" << endl;
		cout << "two keys or two vals the same.)" << endl;
		cout << "Please fix the file \"" << ROCID_MAP_FILENAME << "\"." <<endl;
		exit(-1);
	}
}

//---------------------------------
// ApplyTranslationTable
//---------------------------------
void DTranslationTable::ApplyTranslationTable(JEventLoop *loop) const
{
	/// This will get all of the low level objects and
	/// generate detector hit objects from them, placing
	/// them in the appropriate DANA factories.

	if(VERBOSE>0) ttout << "Entering ApplyTranslationTable:" << endl;
	
	// Containers to hold all of the detector-specific "Digi"
	// objects. Once filled, these will be copied to the
	// respective factories at the end of this method.
	vector<DBCALDigiHit*> vbcal;
	vector<DBCALTDCDigiHit*> vbcaltdc;
	vector<DCDCDigiHit*> vcdc;
	vector<DFCALDigiHit*> vfcal;
	vector<DFDCCathodeDigiHit*> vfdccathode;
	vector<DFDCWireDigiHit*> vfdcwire;
	vector<DSCDigiHit*> vsc;
	vector<DSCTDCDigiHit*> vsctdc;
	vector<DTOFDigiHit*> vtof;
	vector<DTOFTDCDigiHit*> vtoftdc;
	
	// Df250PulseIntegral (will apply Df250PulseTime via associated objects)
	vector<const Df250PulseIntegral*> pulseintegrals250;
	loop->Get(pulseintegrals250);
	if(VERBOSE>2) ttout << "  Number Df250PulseIntegral objects: " << pulseintegrals250.size() << endl;
	for(uint32_t i=0; i<pulseintegrals250.size(); i++){
		const Df250PulseIntegral *pi = pulseintegrals250[i];
		
		// Apply optional rocid translation
		uint32_t rocid = pi->rocid;
		map<uint32_t, uint32_t>::iterator rocid_iter = rocid_map.find(rocid);
		if(rocid_iter != rocid_map.end()) rocid = rocid_iter->second;
		
		if(VERBOSE>4) ttout << "    Looking for rocid:" << rocid <<" slot:" << pi->slot << " chan:" << pi->channel << endl;
		
		// Create crate,slot,channel index and find entry in Translation table.
		// If none is found, then just quietly skip this hit.
		csc_t csc = {rocid, pi->slot, pi->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()){
			if(VERBOSE>6) ttout << "     - Didn't find it" << endl;
			continue;
		}
		const DChannelInfo &chaninfo = iter->second;
		if(VERBOSE>6) ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) << endl;
		
		// Check for a pulse time (this should have been added in JEventSource_EVIO.cc)
		const Df250PulseTime *pt = NULL;
		try{
			pi->GetSingle(pt);
		}catch(...){}
		
		// Create the appropriate hit type based on detector type
		switch(chaninfo.det_sys){
			case BCAL        : vbcal.push_back( MakeBCALDigiHit(chaninfo.bcal, pi, pt) ); break;
			case FCAL        : vfcal.push_back( MakeFCALDigiHit(chaninfo.fcal, pi, pt) ); break;
			case SC          : vsc.push_back  ( MakeSCDigiHit(  chaninfo.sc,   pi, pt) ); break;
			case TOF         : vtof.push_back ( MakeTOFDigiHit( chaninfo.tof,  pi, pt) ); break;

			default:
				if(VERBOSE>4) ttout << "       - Don't know how to make DigiHit objects for this detector type!" << endl;
				break;
		}
	}

	// Df125PulseIntegral (will apply Df125PulseTime via associated objects)
	vector<const Df125PulseIntegral*> pulseintegrals125;
	loop->Get(pulseintegrals125);
	if(VERBOSE>2) ttout << "  Number Df125PulseIntegral objects: " << pulseintegrals125.size() << endl;
	for(uint32_t i=0; i<pulseintegrals125.size(); i++){
		const Df125PulseIntegral *pi = pulseintegrals125[i];

		// Apply optional rocid translation
		uint32_t rocid = pi->rocid;
		map<uint32_t, uint32_t>::iterator rocid_iter = rocid_map.find(rocid);
		if(rocid_iter != rocid_map.end()) rocid = rocid_iter->second;
		
		if(VERBOSE>4) ttout << "    Looking for rocid:" << rocid <<" slot:" << pi->slot << " chan:" << pi->channel << endl;
	
		// Create crate,slot,channel index and find entry in Translation table.
		// If none is found, then just quietly skip this hit.
		csc_t csc = {pi->rocid, pi->slot, pi->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()) {
		    if(VERBOSE>6) ttout << "     - Didn't find it" << endl;
		    continue;
		}
		const DChannelInfo &chaninfo = iter->second;
		if(VERBOSE>6) ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) << endl;

		// Check for a pulse time (this should have been added in JEventSource_EVIO.cc
		const Df125PulseTime *pt = NULL;
		try{
			pi->GetSingle(pt);
		}catch(...){}

		// Create the appropriate hit type based on detector type
		switch(chaninfo.det_sys){
			case CDC         : vcdc.push_back( MakeCDCDigiHit(chaninfo.cdc, pi, pt) ); break;
			case FDC_CATHODES: vfdccathode.push_back( MakeFDCCathodeDigiHit(chaninfo.fdc_cathodes, pi, pt) ); break;

			default: 
			    if(VERBOSE>4) ttout << "       - Don't know how to make DigiHit objects for this detector type!" << endl; 
			    break;
		}
	}

	// DF1TDCHit
	vector<const DF1TDCHit*> f1tdchits;
	loop->Get(f1tdchits);
	if(VERBOSE>2) ttout << "  Number DF1TDCHit objects: " << f1tdchits.size() << endl;
	for(uint32_t i=0; i<f1tdchits.size(); i++){
		const DF1TDCHit *hit = f1tdchits[i];

		// Apply optional rocid translation
		uint32_t rocid = hit->rocid;
		map<uint32_t, uint32_t>::iterator rocid_iter = rocid_map.find(rocid);
		if(rocid_iter != rocid_map.end()) rocid = rocid_iter->second;

		if(VERBOSE>4) ttout << "    Looking for rocid:" << rocid <<" slot:" << hit->slot << " chan:" << hit->channel << endl;
		
		// Create crate,slot,channel index and find entry in Translation table.
		// If none is found, then just quietly skip this hit.
		csc_t csc = {hit->rocid, hit->slot, hit->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()) {
		    if(VERBOSE>6) ttout << "     - Didn't find it" << endl;
		    continue;
		}
		const DChannelInfo &chaninfo = iter->second;
		if(VERBOSE>6) ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) << endl;
		
		// Create the appropriate hit type based on detector type
		switch(chaninfo.det_sys){
			case BCAL     : vbcaltdc.push_back( MakeBCALTDCDigiHit(chaninfo.bcal,      hit) ); break;
			case FDC_WIRES: vfdcwire.push_back( MakeFDCWireDigiHit(chaninfo.fdc_wires, hit) ); break;
			case SC       : vsctdc.push_back( MakeSCTDCDigiHit(chaninfo.sc,      hit) ); break;

			default:  	
			    if(VERBOSE>4) ttout << "       - Don't know how to make DigiHit objects for this detector type!" << endl;
			    break;
		}
	}

	// DCAEN1290TDCHit
	vector<const DCAEN1290TDCHit*> caen1290tdchits;
	loop->Get(caen1290tdchits);
	if(VERBOSE>2) ttout << "  Number DCAEN1290TDCHit objects: " << caen1290tdchits.size() << endl;
	for(uint32_t i=0; i<caen1290tdchits.size(); i++){
		const DCAEN1290TDCHit *hit = caen1290tdchits[i];

		// Apply optional rocid translation
		uint32_t rocid = hit->rocid;
		map<uint32_t, uint32_t>::iterator rocid_iter = rocid_map.find(rocid);
		if(rocid_iter != rocid_map.end()) rocid = rocid_iter->second;

		if(VERBOSE>4) ttout << "    Looking for rocid:" << rocid <<" slot:" << hit->slot << " chan:" << hit->channel << endl;
		
		// Create crate,slot,channel index and find entry in Translation table.
		// If none is found, then just quietly skip this hit.
		csc_t csc = {hit->rocid, hit->slot, hit->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()) {
		    if(VERBOSE>6) ttout << "     - Didn't find it" << endl;
		    continue;
		}
		const DChannelInfo &chaninfo = iter->second;
		if(VERBOSE>6) ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) << endl;
		
		// Create the appropriate hit type based on detector type
		switch(chaninfo.det_sys){
			case TOF      : vtoftdc.push_back( MakeTOFTDCDigiHit(chaninfo.tof,      hit) ); break;

			default:  	
			    if(VERBOSE>4) ttout << "       - Don't know how to make DigiHit objects for this detector type!" << endl;
			    break;
		}
	}

	// Sort object order (this makes it easier to browse with hd_dump)
	sort(vbcal.begin(), vbcal.end(), SortBCALDigiHit);
	
	if(VERBOSE>3){
		ttout << "        vbcal.size() = " << vbcal.size() << endl;
		ttout << "     vbcaltdc.size() = " << vbcaltdc.size() << endl;
		ttout << "         vcdc.size() = " << vcdc.size() << endl;
		ttout << "        vfcal.size() = " << vfcal.size() << endl;
		ttout << "  vfdccathode.size() = " << vfdccathode.size() << endl;
		ttout << "     vfdcwire.size() = " << vfdcwire.size() << endl;
		ttout << "          vsc.size() = " << vsc.size() << endl;
		ttout << "       vsctdc.size() = " << vsctdc.size() << endl;
		ttout << "         vtof.size() = " << vtof.size() << endl;
		ttout << "      vtoftdc.size() = " << vtoftdc.size() << endl;
	}
	
	// Find factory for each container and copy the object pointers into it
	// (n.b. do this even if container is empty since it sets the evnt_called flag)
	CopyToFactory(loop, vbcal);
	CopyToFactory(loop, vbcaltdc);
	CopyToFactory(loop, vcdc);
	CopyToFactory(loop, vfcal);
	CopyToFactory(loop, vfdccathode);
	CopyToFactory(loop, vfdcwire);
	CopyToFactory(loop, vsc);
	CopyToFactory(loop, vsctdc);
	CopyToFactory(loop, vtof);
	CopyToFactory(loop, vtoftdc);

}

//---------------------------------
// MakeBCALDigiHit
//---------------------------------
DBCALDigiHit* DTranslationTable::MakeBCALDigiHit(const BCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const
{
	if(VERBOSE>4) ttout << "       - Making DBCALDigiHit for (mod,lay,sec,end)=(" << idx.module << "," << idx.layer << "," << idx.sector << "," << (DBCALGeometry::End)idx.end << endl;

	DBCALDigiHit *h = new DBCALDigiHit();
	CopyDf250Info(h, pi, pt);

	h->module = idx.module;
	h->layer  = idx.layer;
	h->sector = idx.sector;
	h->end    = (DBCALGeometry::End)idx.end;

	return h;
}

//---------------------------------
// MakeFCALDigiHit
//---------------------------------
DFCALDigiHit* DTranslationTable::MakeFCALDigiHit(const FCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const
{
	DFCALDigiHit *h = new DFCALDigiHit();
	CopyDf250Info(h, pi, pt);

	h->row    = idx.row;
	h->column = idx.col;
	
	return h;
}

//---------------------------------
// MakeTOFDigiHit
//---------------------------------
DTOFDigiHit* DTranslationTable::MakeTOFDigiHit(const TOFIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const
{
	DTOFDigiHit *h = new DTOFDigiHit();
	CopyDf250Info(h, pi, pt);

	h->plane = idx.plane;
	h->bar   = idx.bar;
	h->end   = idx.end;

	return h;
}

//---------------------------------
// MakeSCDigiHit
//---------------------------------
DSCDigiHit* DTranslationTable::MakeSCDigiHit(const SCIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const
{
	DSCDigiHit *h = new DSCDigiHit();
	CopyDf250Info(h, pi, pt);

	h->sector = idx.sector;

	return h;
}

//---------------------------------
// MakeCDCDigiHit
//---------------------------------
DCDCDigiHit* DTranslationTable::MakeCDCDigiHit(const CDCIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const
{
	DCDCDigiHit *h = new DCDCDigiHit();
	CopyDf125Info(h, pi, pt);

	h->ring = idx.ring;
	h->straw = idx.straw;
	
	return h;
}

//---------------------------------
// MakeFDCCathodeDigiHit
//---------------------------------
DFDCCathodeDigiHit* DTranslationTable::MakeFDCCathodeDigiHit(const FDC_CathodesIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const
{
	DFDCCathodeDigiHit *h = new DFDCCathodeDigiHit();
	CopyDf125Info(h, pi, pt);

	h->package    = idx.package;
	h->chamber    = idx.chamber;
	h->view       = idx.view;
	h->strip      = idx.strip;
	h->strip_type = idx.strip_type;

	return h;
}

//---------------------------------
// MakeBCALTDCDigiHit
//---------------------------------
DBCALTDCDigiHit* DTranslationTable::MakeBCALTDCDigiHit(const BCALIndex_t &idx, const DF1TDCHit *hit) const
{
	DBCALTDCDigiHit *h = new DBCALTDCDigiHit();
	CopyDF1TDCInfo(h, hit);

	h->module = idx.module;
	h->layer  = idx.layer;
	h->sector = idx.sector;
	h->end    = (DBCALGeometry::End)idx.end;

	return h;
}

//---------------------------------
// MakeFDCWireDigiHit
//---------------------------------
DFDCWireDigiHit* DTranslationTable::MakeFDCWireDigiHit(const FDC_WiresIndex_t &idx, const DF1TDCHit *hit) const
{
	DFDCWireDigiHit *h = new DFDCWireDigiHit();
	CopyDF1TDCInfo(h, hit);

	h->package = idx.package;
	h->chamber = idx.chamber;
	h->wire    = idx.wire;

	return h;
}

//---------------------------------
// MakeSCTDCDigiHit
//---------------------------------
DSCTDCDigiHit*  DTranslationTable::MakeSCTDCDigiHit(const SCIndex_t &idx, const DF1TDCHit *hit) const
{
	DSCTDCDigiHit *h = new DSCTDCDigiHit();
	CopyDF1TDCInfo(h, hit);

	h->sector = idx.sector;

	return h;
}

//---------------------------------
// MakeTOFTDCDigiHit
//---------------------------------
DTOFTDCDigiHit*  DTranslationTable::MakeTOFTDCDigiHit(const TOFIndex_t &idx, const DCAEN1290TDCHit *hit) const
{
	DTOFTDCDigiHit *h = new DTOFTDCDigiHit();
	CopyDCAEN1290TDCInfo(h, hit);

	h->plane = idx.plane;
	h->bar   = idx.bar;
	h->end   = idx.end;

	return h;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//  The following routines access the translation table
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

static int ModuleStr2ModID(string &type);
static DTranslationTable::Detector_t DetectorStr2DetID(string &type);
static void StartElement(void *userData, const char *xmlname, const char **atts);
static void EndElement(void *userData, const char *xmlname);


//---------------------------------
// ReadTranslationTable
//---------------------------------
void DTranslationTable::ReadTranslationTable(JCalibration *jcalib)
{
	// It seems expat is not thread safe so we lock a mutex here and
	// read in the translation table just once
	pthread_mutex_lock(&tt_mutex);
	if(tt_initialized){
		pthread_mutex_unlock(&tt_mutex);
		return;
	}

	// String to hold entire XML translation table
	string tt_xml; 

	// Try getting it from CCDB first
	if(jcalib && !NO_CCDB){
		map<string,string> tt;
		string namepath = "Translation/DAQ2detector";
		jout << "Reading translation table from calib DB: " << namepath << " ..." << endl;
		jcalib->GetCalib(namepath, tt);
		if(tt.size() != 1){
			jerr << " Error: Unexpected translation table format!" <<endl;
			jerr << "        tt.size()=" << tt.size() << " (expected 1)" <<endl;
		}else{
			// Copy table into tt string
			map<string,string>::iterator iter = tt.begin();
			tt_xml = iter->second;
		}
	}
	
	// If getting from CCDB fails, try just reading in local file
	if(tt_xml.size() == 0){
		if(!NO_CCDB) jout << "Unable to get translation table from CCDB." << endl;
		jout << "Will try reading TT from local file: " << XML_FILENAME << endl;

		// Open file
		ifstream ifs(XML_FILENAME.c_str());
		if(! ifs.is_open()){
			jerr << " Error: Cannot open file! Translation table unavailable." << endl;
			pthread_mutex_unlock(&tt_mutex);
			return;
		}
		
		// read lines into stringstream object 
		stringstream ss;
		while(ifs.good()){
			char line[4096];
			ifs.getline(line, 4096);
			ss << line;
		}

		// Close file
		ifs.close();
		
		// Copy from stringstream to tt
		tt_xml = ss.str();
	}
	
	// create parser and specify element handlers
	XML_Parser xmlParser = XML_ParserCreate(NULL);
	if(xmlParser==NULL) {
		jerr << "readTranslationTable...unable to create parser" << endl;
		exit(EXIT_FAILURE);
	}
	XML_SetElementHandler(xmlParser,StartElement,EndElement);
	XML_SetUserData(xmlParser, &TT);

	// Parse XML string
	int status=XML_Parse(xmlParser, tt_xml.c_str(), tt_xml.size(), 1); // "1" indicates this is the final piece of XML
	if(status==0) {
		jerr << "  ?readTranslationTable...parseXMLFile parse error for " << XML_FILENAME << endl;
		jerr << XML_ErrorString(XML_GetErrorCode(xmlParser)) << endl;
	}
	
	jout << TT.size() << " channels defined in translation table" << endl;
	XML_ParserFree(xmlParser);

	pthread_mutex_unlock(&tt_mutex);
	tt_initialized = true;
}

//---------------------------------
// ModuleStr2ModID
//---------------------------------
int ModuleStr2ModID(string &type)
{
	if(type=="vmecpu") {
		return(DModuleType::VMECPU);
	} else if (type=="tid") {
		return(DModuleType::TID);
	} else if (type=="fadc250") {
		return(DModuleType::FADC250);
	} else if (type=="fadc125") {
		return(DModuleType::FADC125);
	} else if (type=="f1tdcv2") {
		return(DModuleType::F1TDC32);
	} else if (type=="f1tdcv3") {
		return(DModuleType::F1TDC48);
	} else if (type=="jldisc") {
		return(DModuleType::JLAB_DISC);
	} else if (type=="vx1290a") {
		return(DModuleType::CAEN1290);
	} else {
		return(DModuleType::UNKNOWN);
	}
}

//---------------------------------
// DetectorStr2DetID
//---------------------------------
DTranslationTable::Detector_t DetectorStr2DetID(string &type)
{
	if( type=="fdc_cathodes" ) {
		return DTranslationTable::FDC_CATHODES;
	} else if( type=="fdc_wires" ) {
		return DTranslationTable::FDC_WIRES;	
	} else if( type=="bcal" ) {
		return DTranslationTable::BCAL;
	} else if( type=="cdc" ) {
		return DTranslationTable::CDC;	
	} else if( type=="fcal" ) {
		return DTranslationTable::FCAL;	
	} else if( type=="ps" ) {
		return DTranslationTable::PS;
	} else if( type=="psc" ) {
		return DTranslationTable::PSC;
	} else if( type=="st" ) {
	        // The start counter is labelled by "ST" in the translation table
	        // but we stick with the "SC" label in this plugin for consistency
	        // with the rest of the reconstruction software
		return DTranslationTable::SC;
	} else if( type=="tagh" ) {
		return DTranslationTable::TAGH;
	} else if( type=="tagm" ) {
		return DTranslationTable::TAGM;
	} else if( type=="tof" ) {
		return DTranslationTable::TOF;
	} else {
		return DTranslationTable::UNKNOWN_DETECTOR;
	}
}

//---------------------------------
// StartElement
//---------------------------------
void StartElement(void *userData, const char *xmlname, const char **atts) {
	
	static int crate=0, slot=0;
	
	static string type,Type;
	int mc2codaType= 0;
	int channel = 0;
	string Detector;
	int end=0;
	int row=0,column=0,module=0,sector=0,layer=0,chan=0;
	int ring=0,straw=0,plane=0,bar=0,gPlane=0,element=0;
	int package=0,chamber=0,view=0,strip=0,wire=0;
	int id=0, strip_type=0;

	// This complicated line just recasts the userData pointer into
	// a reference to the "TT" member of the DTranslationTable object
	// that called us.
	map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo> &TT = *((map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>*)userData);
	
	// store crate summary info, fill both maps
	if(strcasecmp(xmlname,"halld_online_translation_table")==0) {
		// do nothing
		
	} else if(strcasecmp(xmlname,"crate")==0) {
		for (int i=0; atts[i]; i+=2) {
			if(strcasecmp(atts[i],"number")==0) {
				crate = atoi(atts[i+1]);
				break;
			}
		}
		
	} else if(strcasecmp(xmlname,"slot")==0) {
		for (int i=0; atts[i]; i+=2) {
			if(strcasecmp(atts[i],"number")==0) {
				slot = atoi(atts[i+1]);
			} else if(strcasecmp(atts[i],"type")==0) {
				Type = string(atts[i+1]);
				type = string(atts[i+1]);
				std::transform(type.begin(), type.end(), type.begin(), (int(*)(int)) tolower);
			}
		}
		
		// The detID value set here shows up in the header of the Data Block Bank
		// of the output file. It should be set to one if this crate has JLab
		// made modules that output in the standard format (see document:
		// "VME Data Format Standards for JLAB Modules"). These would include
		// f250ADC, f125ADC, F1TDC, .... Slots containing other types of modules
		// (e.g. CAEN1290) should have their own unique detID. We use detID of
		// zero for non-digitizing modules like CPUs nd TIDs even though potentially,
		// one could read data from these.
		mc2codaType = ModuleStr2ModID(type);		
		
	} else if(strcasecmp(xmlname,"channel")==0) {
		
		for (int i=0; atts[i]; i+=2) {
			string tag(atts[i+0]);
			string sval(atts[i+1]);
			int ival = atoi(atts[i+1]);

			if(     tag == "number"   ) channel   = ival;
			else if(tag == "detector" ) Detector  = sval;
			else if(tag == "row"      ) row       = ival;
			else if(tag == "column"   ) column    = ival;
			else if(tag == "col"      ) column    = ival;
			else if(tag == "module"   ) module    = ival;
			else if(tag == "sector"   ) sector    = ival;
			else if(tag == "layer"    ) layer     = ival;
			else if(tag == "chan"     ) chan      = ival;
			else if(tag == "ring"     ) ring      = ival;
			else if(tag == "straw"    ) straw     = ival;
			else if(tag == "gPlane"   ) gPlane    = ival;
			else if(tag == "element"  ) element   = ival;
			else if(tag == "plane"    ) plane     = ival;
			else if(tag == "bar"      ) bar       = ival;
			else if(tag == "package"  ) package   = ival;
			else if(tag == "chamber"  ) chamber   = ival;
			else if(tag == "view"     ) {
			        if(     sval == "U"  ) view=1;
				else if(sval == "D"  ) view=3;
			}
			else if(tag == "strip"    ) strip     = ival;
			else if(tag == "wire"     ) wire      = ival;
			else if(tag == "id"       ) id        = ival;
			else if(tag == "end"      ){
				if(     sval == "U"  ) {end = DBCALGeometry::kUpstream;   view=1;}
				else if(sval == "D"  ) {end = DBCALGeometry::kDownstream; view=3;}
				else if(sval == "N"  ) end = 0; // TOF north
				else if(sval == "S"  ) end = 1; // TOF south
				else if(sval == "UP" ) end = 0; // TOF up
				else if(sval == "DW" ) end = 1; // TOF down
			}
			else if(tag == "strip_type"){
				if(     sval == "full" ) strip_type = 1;
				else if(sval == "A"    ) strip_type = 2;
				else if(sval == "B"    ) strip_type = 3;
			}


		}
		
		// ignore certain module types
		if(type == "disc") return;
		if(type == "ctp") return;
		if(type == "sd") return;
		if(type == "a1535sn") return;

		
//		// Data integrity check
//		if(crate<0 || crate>=MAXDCRATE){
//			jerr << " Crate value of "<<crate<<" is not in range 0 <= crate < " << MAXDCRATE << endl;
//			exit(-1);
//		}
//		
//		if(slot<0 || slot>=MAXDSLOT){
//			jerr << " Slot value of "<<slot<<" is not in range 0 <= slot < " << MAXDSLOT << endl;
//			exit(-1);
//		}
//		
//		if(channel<0 || channel>=MAXDCHANNEL){
//			jerr << " Crate value of "<<channel<<" is not in range 0 <= channel < " << MAXDCHANNEL << endl;
//			exit(-1);
//		}
		
		// fill maps
		
		DTranslationTable::csc_t csc = {crate,slot,channel};
		string detector = Detector;
		std::transform(detector.begin(), detector.end(), detector.begin(), (int(*)(int)) tolower);
		
		//string s="unknown::";

		// Common indexes
		DTranslationTable::DChannelInfo &ci = TT[csc];
		ci.CSC = csc;
		ci.module_type = (DModuleType::type_id_t)mc2codaType;
		ci.det_sys = DetectorStr2DetID(detector);

		// detector-specific indexes
		switch(ci.det_sys){
			case DTranslationTable::BCAL:
				ci.bcal.module = module;
				ci.bcal.layer = layer;
				ci.bcal.sector = sector;
				ci.bcal.end = end;
				break;
			case DTranslationTable::CDC:
				ci.cdc.ring = ring;
				ci.cdc.straw = straw;
				break;
			case DTranslationTable::FCAL:
				ci.fcal.row = row;
				ci.fcal.col = column;
				break;
			case DTranslationTable::FDC_CATHODES:
				ci.fdc_cathodes.package = package;
				ci.fdc_cathodes.chamber = chamber;
				ci.fdc_cathodes.view = view;
				ci.fdc_cathodes.strip = strip;
				ci.fdc_cathodes.strip_type = strip_type;
				break;
			case DTranslationTable::FDC_WIRES:
				ci.fdc_wires.package = package;
				ci.fdc_wires.chamber = chamber;
				ci.fdc_wires.wire = wire;
				break;
			case DTranslationTable::SC:
				ci.sc.sector = sector;
				break;
			case DTranslationTable::TAGH:
				ci.tagh.id = id;
				break;
			case DTranslationTable::TAGM:
				ci.tagm.col = column;
				ci.tagm.row = row;
				break;
			case DTranslationTable::TOF:
				ci.tof.plane = plane;
				ci.tof.bar = bar;
				ci.tof.end = end;
				break;
			case DTranslationTable::PS:
			case DTranslationTable::PSC:
			case DTranslationTable::UNKNOWN_DETECTOR:
				break;
		}

	} else {
		jerr << endl << endl << "?startElement...unknown xml tag " << xmlname << endl << endl;
	}
	
}


//--------------------------------------------------------------------------


void EndElement(void *userData, const char *xmlname) {
	// nothing to do yet...
}


//--------------------------------------------------------------------------




