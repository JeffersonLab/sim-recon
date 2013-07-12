// $Id$
//
//    File: DTranslationTable.cc
// Created: Thu Jun 27 16:07:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#include "DTranslationTable.h"

#include <expat.h>

#include <DAQ/DModuleType.h>

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
	// Read in Translation table. This will create DChannelInfo objects
	// and store them in the "TT" map, indexed by csc_t objects
	ReadTranslationTable();

	// These are used to help the event source report which
	// types of data it is capable of providing. For practical
	// purposes, these types are provided by the source
	// because they are generated and placed into their
	// respective JANA factories during a call to GetObjects().
	// The source is responsible for reporting the types it is
	// directly responsible for (e.g. Df250PulseIntegral)
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
	
	// These will need to get filled in stages. In particular, things like
	// DBCALHit objects contain both the amplitude from Df250PulseIntegral
	// and time from D250PulseTime objects. DPulseIntegral objects are
	// looped over first and the DPulseTime next, adding their information
	// to any existing objects if found or creating new ones if needed.
	vector<DBCALHit*> vbcal;
	vector<DCDCHit*> vcdc;
	vector<DFCALHit*> vfcal;
	vector<DFDCHit*> vfdc;
	vector<DSCHit*> vsc;
	vector<DTOFRawHit*> vtof;

	vector<DBCALTDCHit*> vbcaltdc;
	
	// Df250PulseIntegral (will apply Df250PulseTime via associated objects)
	vector<const Df250PulseIntegral*> pulseintegrals;
	loop->Get(pulseintegrals);
	for(uint32_t i=0; i<pulseintegrals.size(); i++){
		const Df250PulseIntegral *pi = pulseintegrals[i];
		
		// Create crate,slot,channel index and find entry in Translation table.
		// If none is found, then just quietly skip this hit.
		csc_t csc = {pi->rocid, pi->slot, pi->channel};
		map<csc_t, DChannelInfo>::const_iterator iter = TT.find(csc);
		if(iter == TT.end()) continue;
		const DChannelInfo &chaninfo = iter->second;
		
		// Check for a pulse time (this should have been added in JEventSource_EVIO.cc
		// If a trigger time for the module is also associated, subtract that to give
		// a time relative to the trigger.
		uint32_t t = 0;
		try{
			const Df250PulseTime *pt;
			pi->GetSingle(pt);
			t = pt->time;
			
//			const Df250TriggerTime *trigt;
//			pi->GetSingle(trigt);
//			t -= trigt->time;
		}catch(...){}
		
		// Create the appropriate hit type based on detector type
		switch(chaninfo.det_sys){
			case BCAL        : vbcal.push_back( MakeBCALHit(chaninfo.bcal, pi, t) ); break;
//			case CDC         : vcdc.push_back ( MakeCDCHit(pi, chaninfo.cdc) ); break;
//			case FCAL        : vfcal.push_back( MakeFCALHit(pi, chaninfo.fcal) ); break;
//			case FDC_CATHODES: vfdc.push_back ( MakeFDCCathodeHit(pi, chaninfo.fdc_cathodes) ); break;
//			case FDC_WIRES   : vfdc.push_back ( MakeFDCWireHit(pi, chaninfo.fdc_wires) ); break;
//			case SC          : vsc.push_back  ( MakeSCHit(pi, chaninfo.sc) ); break;
//			case TOF         : vtof.push_back ( MakeTOFHit(pi, chaninfo.tof) ); break;

			default:  break;
		}
	}
	
	JFactory<DBCALHit>   *fac_bcal = dynamic_cast<JFactory<DBCALHit>*>(loop->GetFactory("DBCALHit"));
	JFactory<DCDCHit>    *fac_cdc  = dynamic_cast<JFactory<DCDCHit>*>(loop->GetFactory("DCDCHit"));
	JFactory<DFCALHit>   *fac_fcal = dynamic_cast<JFactory<DFCALHit>*>(loop->GetFactory("DFCALHit"));
	JFactory<DFDCHit>    *fac_fdc  = dynamic_cast<JFactory<DFDCHit>*>(loop->GetFactory("DFDCHit"));
	JFactory<DSCHit>     *fac_sc   = dynamic_cast<JFactory<DSCHit>*>(loop->GetFactory("DSCHit"));
	JFactory<DTOFRawHit> *fac_tof  = dynamic_cast<JFactory<DTOFRawHit>*>(loop->GetFactory("DTOFRawHit"));

	JFactory<DBCALTDCHit>   *fac_bcaltdc = dynamic_cast<JFactory<DBCALTDCHit>*>(loop->GetFactory("DBCALTDCHit"));
	
	if(fac_bcal) fac_bcal->CopyTo(vbcal);
	if(fac_cdc ) fac_cdc->CopyTo(vcdc);
	if(fac_fcal) fac_fcal->CopyTo(vfcal);
	if(fac_fdc ) fac_fdc->CopyTo(vfdc);
	if(fac_sc  ) fac_sc->CopyTo(vsc);
	if(fac_tof ) fac_tof->CopyTo(vtof);

	if(fac_bcaltdc) fac_bcaltdc->CopyTo(vbcaltdc);
}

//---------------------------------
// MakeBCALHit
//---------------------------------
DBCALHit* DTranslationTable::MakeBCALHit(const BCALIndex_t &idx, const Df250PulseIntegral *hit, uint32_t t) const
{
	DBCALHit *h = new DBCALHit();
	h->module = idx.module;
	h->layer = idx.layer;
	h->sector = idx.sector;
	h->end = (DBCALGeometry::End)idx.end;
	
	h->E = (float)hit->integral;
	h->t = (float)t;
	
	h->AddAssociatedObject(hit);
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
void DTranslationTable::ReadTranslationTable(void)
{
	
	string translationTableName = "tt.xml";
	jout << "Reading translation table " << translationTableName << endl;
	
	// create parser and specify element handlers
	XML_Parser xmlParser = XML_ParserCreate(NULL);
	if(xmlParser==NULL) {
		jerr << "readTranslationTable...unable to create parser" << endl;
		exit(EXIT_FAILURE);
	}
	XML_SetElementHandler(xmlParser,StartElement,EndElement);
	XML_SetUserData(xmlParser, &TT);
	
	// open and parse the file
	FILE *f = fopen(translationTableName.c_str(),"r");
	if(f!=NULL) {
		int status,len;
		bool done;
		const int bufSize = 50000;
		char *buf = new char[bufSize];
		do {
			len  = fread(buf,1,bufSize,f);
			done = len!=bufSize;
			status=XML_Parse(xmlParser,buf,len,done);
			
			if((!done)&&(status==0)) {
				jerr << "  ?readTranslationTable...parseXMLFile parse error for " << translationTableName << endl;
				jerr << XML_ErrorString(XML_GetErrorCode(xmlParser)) << endl;
				fclose(f);
				delete [] buf;
				XML_ParserFree(xmlParser);
				return;
			}
		} while (!done);
		fclose(f);
		delete [] buf;
		jout << TT.size() << " channels defined in translation table" << endl;
		
	} else {
		jerr << "  ?readTranslationTable...unable to open " << translationTableName << endl;
	}
	XML_ParserFree(xmlParser);
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
	} else if( type=="sc" ) {
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
	int mc2codaType;
	int channel = 0;
	string Detector;
	int end;
	int row,column,module,sector,layer,chan;
	int ring,straw,plane,bar,gPlane,element;
	int package,chamber,view,strip,wire;
	int id, strip_type;

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
//		if(mc2codaType!=DModuleType::UNKNOWN) {
//			nModules[nCrate-1]++;
//			modules[nCrate-1][slot-1] = mc2codaType;
//			switch(mc2codaType){
//				case FADC250:
//				case FADC125:
//				case F1TDC32:
//				case F1TDC48:
//					detID[nCrate-1][slot-1]   = 1;
//					break;
//				case CAEN1190:
//				case CAEN1290:
//					detID[nCrate-1][slot-1]   = 2;
//					break;
//				default:
//					detID[nCrate-1][slot-1]   = 0;
//			}
//		}
		
		
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
			else if(tag == "view"     ) view      = ival;
			else if(tag == "strip"    ) strip     = ival;
			else if(tag == "wire"     ) wire      = ival;
			else if(tag == "id"       ) id        = ival;
			else if(tag == "end"      ){
				if(     sval == "U"  ) end = DBCALGeometry::kUpstream;
				else if(sval == "D"  ) end = DBCALGeometry::kDownstream;
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


//		
//			if(strcasecmp(atts[i],"number")==0) {
//				channel = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"detector")==0) {
//				Detector = string(atts[i+1]);
//			} else if(strcasecmp(atts[i],"row")==0) {
//				row = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"column")==0) {
//				column = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"col")==0) {
//				column = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"module")==0) {
//				module = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"sector")==0) {
//				sector = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"layer")==0) {
//				layer = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"end")==0) {
//				end = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"chan")==0) {
//				chan = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"ring")==0) {
//				ring = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"straw")==0) {
//				straw = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"gPlane")==0) {
//				gPlane = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"element")==0) {
//				element = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"plane")==0) {
//				plane = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"bar")==0) {
//				bar = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"package")==0) {
//				package = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"chamber")==0) {
//				chamber = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"view")==0) {
//				view = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"strip")==0) {
//				strip = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"wire")==0) {
//				wire = atoi(atts[i+1]);
//			} else if(strcasecmp(atts[i],"id")==0) {
//				id = atoi(atts[i+1]);
//			}
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
		
//		if(detector=="fcal") {
//			if(type!="fadc250") _DBG_ << "?startElement...illegal type for " << detector << ": " << Type << endl;
//
//			ci.module_type = mc2codaType;
//			ci.det_sys = DTranslationTable.FCAL;
//			ci.fcal.row = row;
//			ci.fcal.col = column;
//			
//		} else if(detector=="bcal") {
//			if(type=="f1tdcv2") {
//				s = "bcaltdc::";
//			} else if (type=="fadc250") {
//				s = "bcaladc::";
//			} else {
//				s = "unknownBCAL::";
//				jerr << endl << endl << "?startElement...illegal type for BCAL: " << Type << " ("<<type<<")" << endl << endl;
//			}
//			s += module + ":" + sector + ":" + layer + ":" + end;
//			cscMap[s] = csc;
//			if(type!="fadc250") _DBG_ << "?startElement...illegal type for " << detector << ": " << Type << endl;
//			
//			DChannelInfo &ci = TT[csc];
//			ci.CSC = csc;
//			ci.module_type = mc2codaType;
//			ci.det_sys = DTranslationTable.FCAL;
//			ci.fcal.row = row;
//			ci.fcal.col = column;
//			
//			
//		} else if(detector=="cdc") {
//			if(type=="fadc125") {
//				s = "cdcadc::";
//			} else {
//				s = "unknownCDC::";
//				jerr << endl << endl << "?startElement...illegal type for CDC: " << Type << endl << endl;
//			}
//			s += ring + ":" + straw;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="st") {
//			if(type=="f1tdcv2") {
//				s = "sttdc::";
//			} else if (type=="fadc250") {
//				s = "stadc::";
//			} else {
//				s = "unknownST::";
//				jerr << endl << endl << "?startElement...illegal type for ST: " << Type << endl << endl;
//			}
//			s += sector;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="fdc_cathodes") {
//			int ipackage = atoi(package.c_str());
//			int ichamber = atoi(chamber.c_str());
//			int igPlane = 1 + (ipackage-1)*6*3 + (ichamber-1)*3 + (view=="U" ? 0:2);
//			stringstream ss;
//			ss<<igPlane;
//			gPlane = ss.str();
//			if(type=="fadc125") {
//				s = "fdccathode::";
//			} else {
//				s = "unknownFDCCathode::";
//				jerr << endl << endl << "?startElement...illegal type for FDC Cathode: " << Type << endl << endl;
//			}
//			s += gPlane + ":" + strip;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="fdc_wires") {
//			int ipackage = atoi(package.c_str());
//			int ichamber = atoi(chamber.c_str());
//			int igPlane = 2 + (ipackage-1)*6*3 + (ichamber-1)*3;
//			stringstream ss;
//			ss<<igPlane;
//			gPlane = ss.str();
//			if(type=="f1tdcv3") {
//				s = "fdcanode::";
//			} else {
//				s = "unknownFDCAnode::";
//				jerr << endl << endl << "?startElement...illegal type for FDC Anode: " << Type << endl << endl;
//			}
//			s += gPlane + ":" + wire;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="tof") {
//			if(type=="vx1290a") {
//				s = "toftdc::";
//			} else if (type=="fadc250") {
//				s = "tofadc::";
//			} else {
//				s = "unknownTOF::";
//				jerr << endl << endl << "?startElement...illegal type for TOF: " << Type << endl << endl;
//			}
//			s += plane + ":" + bar + ":" + end;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="tagh") {
//			if(type=="f1tdcv2") {
//				s = "taghtdc::";
//			}else if (type=="fadc250") {
//				s = "taghadc::";
//			} else {
//				s = "unknownTagger::";
//				jerr << endl << endl << "?startElement...illegal type for TAGH: " << Type << endl << endl;
//			}
//			s += row + ":" + column;
//			cscMap[s] = csc;
//			
//		} else if(detector=="tagm") {
//			if(type=="f1tdcv2") {
//				s = "tagmtdc::";
//			} else if (type=="fadc250") {
//				s = "tagmadc::";
//			} else {
//				s = "unknownTagger::";
//				jerr << endl << endl << "?startElement...illegal type for TAGM: " << Type << endl << endl;
//			}
//			s += row + ":" + column;
//			cscMap[s] = csc;
//			
//			
//		} else if(detector=="psc") {
//			if(type=="f1tdcv2") {
//				s = "psctdc::";
//			} else if (type=="fadc250") {
//				s = "pscadc::";
//			} else {
//				s = "unknownPSC::";
//				jerr << endl << endl << "?startElement...illegal type for PSC: " << Type << endl << endl;
//			}
//			s += row + ":" + column;
//			cscMap[s] = csc;
//			
//		} else if(detector=="ps") {
//			if (type=="fadc250") {
//				s = "psadc::";
//			} else {
//				s = "unknownPS::";
//				jerr << endl << endl << "?startElement...illegal type for PS: " << Type << endl << endl;
//			}
//			s += row + ":" + column;
//			cscMap[s] = csc;
//			
//			
//		} else {
//			jerr << endl << endl << "?startElement...unknown detector " << Detector << endl << endl;
//		}
//		
//		
//		
//		// fill detector map, index is crate,slot,channel
//		detectorMap[crate][slot][channel] = s;
//		
//		
//		
	} else {
		jerr << endl << endl << "?startElement...unknown xml tag " << xmlname << endl << endl;
	}
	
}


//--------------------------------------------------------------------------


void EndElement(void *userData, const char *xmlname) {
	// nothing to do yet...
}


//--------------------------------------------------------------------------




