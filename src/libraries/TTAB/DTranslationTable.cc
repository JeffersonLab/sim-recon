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
#include <DAQ/JEventSource_EVIO.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>

using namespace jana;
using namespace std;

// Use one translation table for all threads
pthread_mutex_t& DTranslationTable::Get_TT_Mutex(void) const
{
	static pthread_mutex_t tt_mutex = PTHREAD_MUTEX_INITIALIZER;
	return tt_mutex;
}

bool& DTranslationTable::Get_TT_Initialized(void) const
{
	static bool tt_initialized = false;
	return tt_initialized;
}

map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>& DTranslationTable::Get_TT(void) const
{
	static map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo> TT;
	return TT;
}

map<uint32_t, uint32_t>& DTranslationTable::Get_ROCID_Map(void) const
{
	static map<uint32_t, uint32_t> rocid_map;     // (see ReadOptionalROCidTranslation() for details)
	return rocid_map;
}

map<uint32_t, uint32_t>& DTranslationTable::Get_ROCID_Inv_Map(void) const
{
	static map<uint32_t, uint32_t> rocid_inv_map; // (see ReadOptionalROCidTranslation() for details)
	return rocid_inv_map;
}

map<DTranslationTable::Detector_t, set<uint32_t> >& DTranslationTable::Get_ROCID_By_System(void)
{
	static map<DTranslationTable::Detector_t, set<uint32_t> > rocid_by_system;
	return rocid_by_system;
}

//...................................
// Less than operator for csc_t data types. This is used by
// the map<csc_t, XX> to order the entires by key
bool operator<(const DTranslationTable::csc_t &a, const DTranslationTable::csc_t &b) {
   if (a.rocid < b.rocid) return true;
   if (a.rocid > b.rocid) return false;
   if (a.slot < b.slot) return true;
   if (a.slot > b.slot) return false;
   if (a.channel < b.channel) return true;
   return false;
}

//...................................
// sort functions
bool SortBCALDigiHit(const DBCALDigiHit *a, const DBCALDigiHit *b) {
   if (a->module == b->module) {
      if (a->layer == b->layer) {
         if (a->sector == b->sector) {
            if (a->end == b->end) {
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
   SYSTEMS_TO_PARSE = "";
   gPARMS->SetDefaultParameter("TT:NO_CCDB", NO_CCDB, 
           "Don't try getting translation table from CCDB and just look"
           " for file. Only useful if you want to force reading tt.xml."
           " This is automatically set if you specify a different"
           " filename via the TT:XML_FILENAME parameter.");
   JParameter *p = gPARMS->SetDefaultParameter("TT:XML_FILENAME", XML_FILENAME,
           "Fallback filename of translation table XML file."
           " If set to non-default, CCDB will not be checked.");
   if (p->GetDefault() != p->GetValue())
     NO_CCDB = true;
   gPARMS->SetDefaultParameter("TT:VERBOSE", VERBOSE, 
           "Verbosity level for Applying Translation Table."
           " 0=no messages, 10=all messages.");
   
   ROCID_MAP_FILENAME = "rocid.map";
   gPARMS->SetDefaultParameter("TT:ROCID_MAP_FILENAME", ROCID_MAP_FILENAME,
           "Optional rocid to rocid conversion map for use with files"
           " generated with the non-standard rocid's");

	gPARMS->SetDefaultParameter("TT:SYSTEMS_TO_PARSE", SYSTEMS_TO_PARSE,
			"Comma separated list of systems to parse EVIO data for. "
			"Default is empty string which means to parse all. System "
			"names should be what is returned by DTranslationTable::DetectorName() .");

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
   supplied_data_types.insert("DRFDigiTime");
   supplied_data_types.insert("DSCDigiHit");
   supplied_data_types.insert("DSCTDCDigiHit");
   supplied_data_types.insert("DTOFDigiHit");
   supplied_data_types.insert("DTOFTDCDigiHit");
   supplied_data_types.insert("DTAGMDigiHit");
   supplied_data_types.insert("DTAGMTDCDigiHit");
   supplied_data_types.insert("DTAGHDigiHit");
   supplied_data_types.insert("DTAGHTDCDigiHit");
   supplied_data_types.insert("DPSDigiHit");
   supplied_data_types.insert("DPSCDigiHit");
   supplied_data_types.insert("DPSCTDCDigiHit");
   supplied_data_types.insert("DTPOLSectorDigiHit");
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
   if (!ifs.is_open()) return;
   
   std::cout << "Opened ROC id translation map: " << ROCID_MAP_FILENAME << std::endl;
   while (ifs.good()) {
      char line[256];
      ifs.getline(line, 256);
      if (ifs.gcount() < 1) break;
      if (line[0] == '#') continue;

      stringstream ss(line);
      uint32_t from=10000, to=10000;
      ss >> from >> to;  // from=evio  to=TT
      if ( to == 10000 ) {
         if ( from != 10000) {
            std::cout << "unable to convert line:" << std::endl;
            std::cout << "  " << line;
         }
      }else{
         Get_ROCID_Map()[from] = to;
         Get_ROCID_Inv_Map()[to] = from;
      }
   }
   ifs.close();
   
   if (Get_ROCID_Map().size() == Get_ROCID_Inv_Map().size()) {
      std::cout << "   Read " << Get_ROCID_Map().size() << " entries" << std::endl;
      map<uint32_t,uint32_t>::iterator iter;
      for (iter=Get_ROCID_Map().begin(); iter != Get_ROCID_Map().end(); iter++) {
         std::cout << "   rocid " << iter->first << " -> rocid "
                    << iter->second << std::endl;
      }
   }else{
      std::cout << "Entries not unique! This can happen if there are"
                << std::endl;
      std::cout << "more than one entry with the same value (either" 
                << std::endl;
      std::cout << "two keys or two vals the same.)" 
                << std::endl;
      std::cout << "Please fix the file \"" << ROCID_MAP_FILENAME << "\"." 
                << std::endl;
      exit(-1);
   }
}

//---------------------------------
// SetSystemsToParse
//---------------------------------
void DTranslationTable::SetSystemsToParse(string systems, JEventSource *eventsource)
{
	/// This takes a string of comma separated system names and
	/// identifies a list of Detector_t values from this (using
	/// strings returned by DetectorName() ). It then tries to 
	/// copy the value into the DAQ plugin so they can be used
	/// to restrict which banks to parse.
	
	if(systems == "") return; // nothing to do for empty strings

	// Make sure this is a JEventSource_EVIO object pointer
	JEventSource_EVIO *eviosource = dynamic_cast<JEventSource_EVIO*>(eventsource);
	if(!eviosource) {
		jerr << "eventsource not a JEventSource_EVIO object! Cannot restrict parsing list!" << endl;
		return;
	}

#ifndef HAVE_EVIO

    // If we don't have EVIO, then the JEventSource_EVIO objects is crippled
    // to the point where it can't be used. If this is the case, don't try
    // doing anything else.

#else // HAVE_EVIO

    // Make map of system type id by name
	map<string, Detector_t> name_to_id;
	for(uint32_t dettype=UNKNOWN_DETECTOR; dettype<NUM_DETECTOR_TYPES; dettype++){
		name_to_id[DetectorName((Detector_t)dettype)] = (Detector_t)dettype;
	}
	
	// Parse string of system names
	std::istringstream ss(systems);
	std::string token;
	while(std::getline(ss, token, ',')) {
		
		// Get initial list of rocids based on token
		set<uint32_t> rocids = Get_ROCID_By_System()[name_to_id[token]];
		
		// Let "FDC" be an alias for both cathode strips and wires
		if(token == "FDC"){
			set<uint32_t> rocids1 = Get_ROCID_By_System()[name_to_id["FDC_CATHODES"]];
			set<uint32_t> rocids2 = Get_ROCID_By_System()[name_to_id["FDC_WIRES"]];
			rocids.insert(rocids1.begin(), rocids1.end());
			rocids.insert(rocids2.begin(), rocids2.end());
		}

		// More likely than not, someone specifying "PS" will also want "PSC" 
		if(token == "PS"){
			set<uint32_t> rocids1 = Get_ROCID_By_System()[name_to_id["PSC"]];
			rocids.insert(rocids1.begin(), rocids1.end());
		}

		set<uint32_t>::iterator it;
		for(it=rocids.begin(); it!=rocids.end(); it++){

			// Add this rocid to the DAQ parsing list
			uint32_t rocid = *it;
			eviosource->AddROCIDtoParseList(rocid);	
			if(VERBOSE>0) ttout << "Added rocid " << rocid << " for system " << token << " to parse list" << endl;
		}
	}

#endif //HAVE_EVIO
}

//---------------------------------
// ApplyTranslationTable
//---------------------------------
void DTranslationTable::ApplyTranslationTable(JEventLoop *loop) const
{
   /// This will get all of the low level objects and
   /// generate detector hit objects from them, placing
   /// them in the appropriate DANA factories.

   if (VERBOSE > 0)
      ttout << "Entering ApplyTranslationTable:" << std::endl;
   
   // If the JANA call stack is being recorded, then temporarily disable it
   // so calls we make to loop->Get() here are ignored. The reason we do this
   // is because this routine is called while already in a loop->Get() call
   // so JANA will treat all other loop->Get() calls we make as being dependencies
   // of the loop->Get() call that we are already in. (Confusing eh?) 
   bool record_call_stack = loop->GetCallStackRecordingStatus();
   if (record_call_stack) loop->DisableCallStackRecording();
   
   // Containers to hold all of the detector-specific "Digi"
   // objects. Once filled, these will be copied to the
   // respective factories at the end of this method.
   vector<DBCALDigiHit*> vbcal;
   vector<DBCALTDCDigiHit*> vbcaltdc;
   vector<DCDCDigiHit*> vcdc;
   vector<DFCALDigiHit*> vfcal;
   vector<DFDCCathodeDigiHit*> vfdccathode;
   vector<DFDCWireDigiHit*> vfdcwire;
   vector<DRFDigiTime*> vrf;
   vector<DRFTDCDigiTime*> vrftdc;
   vector<DSCDigiHit*> vsc;
   vector<DSCTDCDigiHit*> vsctdc;
   vector<DTOFDigiHit*> vtof;
   vector<DTOFTDCDigiHit*> vtoftdc;
   vector<DTAGMDigiHit*> vtagm;
   vector<DTAGMTDCDigiHit*> vtagmtdc;
   vector<DTAGHDigiHit*> vtagh;
   vector<DTAGHTDCDigiHit*> vtaghtdc;
   vector<DPSDigiHit*> vps;
   vector<DPSCDigiHit*> vpsc;
   vector<DPSCTDCDigiHit*> vpsctdc;
   vector<DTPOLSectorDigiHit*> vtpolsector;

   // Df250PulseIntegral (will apply Df250PulseTime via associated objects)
   vector<const Df250PulseIntegral*> pulseintegrals250;
   loop->Get(pulseintegrals250);
   if (VERBOSE > 2)
      ttout << "  Number Df250PulseIntegral objects: " 
            << pulseintegrals250.size() << std::endl;
   for (uint32_t i=0; i<pulseintegrals250.size(); i++) {
      const Df250PulseIntegral *pi = pulseintegrals250[i];
      
      // Apply optional rocid translation
      uint32_t rocid = pi->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;
      
      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << pi->slot
               << " chan:" << pi->channel << std::endl;
      
      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, pi->slot, pi->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
         if (VERBOSE > 6)
            ttout << "     - Didn't find it" << std::endl;
         continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6)
         ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys)
               << std::endl;
      
      // Check for a pulse time (this should have been added in JEventSource_EVIO.cc)
      const Df250PulseTime *pt = NULL;
      const Df250PulsePedestal *pp = NULL;
	  pi->GetSingle(pt);
	  pi->GetSingle(pp);
      
      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
         case BCAL:
            vbcal.push_back( MakeBCALDigiHit(chaninfo.bcal, pi, pt, pp) );
            break;
         case FCAL:
            vfcal.push_back( MakeFCALDigiHit(chaninfo.fcal, pi, pt, pp) );
            break;
         case SC:
            vsc.push_back  ( MakeSCDigiHit(  chaninfo.sc, pi, pt, pp) );
            break;
         case TOF:
            vtof.push_back ( MakeTOFDigiHit( chaninfo.tof, pi, pt, pp) );
            break;
         case TAGM:
            vtagm.push_back( MakeTAGMDigiHit(chaninfo.tagm, pi, pt, pp) );
            break;
         case TAGH:
            vtagh.push_back( MakeTAGHDigiHit(chaninfo.tagh, pi, pt, pp) );
            break;
         case PS:
            vps.push_back  ( MakePSDigiHit(chaninfo.ps, pi, pt, pp) );
            break;
         case PSC:
            vpsc.push_back ( MakePSCDigiHit(chaninfo.psc, pi, pt, pp) );
            break;
         case RF:
        	vrf.push_back( MakeRFDigiTime(chaninfo.rf, pt) );
            break;
         case TPOLSECTOR:
            vtpolsector.push_back( MakeTPOLSectorDigiHit(chaninfo.tpolsector, pi, pt, pp) );
            break;
         default:
            if (VERBOSE > 4)
               ttout << "       - Don't know how to make DigiHit objects"
                     << " for this detector type!" << std::endl;
            break;
      }
   }

   // Df125PulseIntegral (will apply Df125PulseTime via associated objects)
   vector<const Df125PulseIntegral*> pulseintegrals125;
   loop->Get(pulseintegrals125);
   if (VERBOSE > 2)
      ttout << "  Number Df125PulseIntegral objects: "
            << pulseintegrals125.size() << std::endl;
   for (uint32_t i=0; i<pulseintegrals125.size(); i++) {
      const Df125PulseIntegral *pi = pulseintegrals125[i];

      // Apply optional rocid translation
      uint32_t rocid = pi->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;
      
      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << pi->slot
               << " chan:" << pi->channel << std::endl;
   
      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, pi->slot, pi->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
          if (VERBOSE > 6)
             ttout << "     - Didn't find it" << std::endl;
          continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6)
         ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) 
               << std::endl;

      // Check for a pulse time (this should have been added in JEventSource_EVIO.cc
      const Df125PulseTime *pt = NULL;
      const Df125PulsePedestal *pp = NULL;
	  pi->GetSingle(pt);
	  pi->GetSingle(pp);

      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
         case CDC:
            vcdc.push_back( MakeCDCDigiHit(chaninfo.cdc, pi, pt, pp) );
            break;
         case FDC_CATHODES:
            vfdccathode.push_back( MakeFDCCathodeDigiHit(chaninfo.fdc_cathodes, pi, pt, pp) );
            break;
         default: 
             if (VERBOSE > 4)
                ttout << "       - Don't know how to make DigiHit"
                      << " objects for this detector type!" << std::endl; 
             break;
      }
   }

   // Df125CDCPulse
   vector<const Df125CDCPulse*> cdcpulses;
   loop->Get(cdcpulses);
   if (VERBOSE > 2)
      ttout << "  Number Df125CDCPulse objects: "
            << cdcpulses.size() << std::endl;
   for (uint32_t i=0; i<cdcpulses.size(); i++) {
      const Df125CDCPulse *p = cdcpulses[i];

      // Apply optional rocid translation
      uint32_t rocid = p->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;
      
      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << p->slot
               << " chan:" << p->channel << std::endl;
   
      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, p->slot, p->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
          if (VERBOSE > 6)
             ttout << "     - Didn't find it" << std::endl;
          continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6)
         ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) 
               << std::endl;

      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
         case CDC:
            vcdc.push_back( MakeCDCDigiHit(chaninfo.cdc, p) );
            break;
         default: 
             if (VERBOSE > 4)
                ttout << "       - Don't know how to make DigiHit"
                      << " objects for this detector type!" << std::endl; 
             break;
      }
   }

   // Df125FDCPulse
   vector<const Df125FDCPulse*> fdcpulses;
   loop->Get(fdcpulses);
   if (VERBOSE > 2)
      ttout << "  Number Df125FDCPulse objects: "
            << fdcpulses.size() << std::endl;
   for (uint32_t i=0; i<fdcpulses.size(); i++) {
      const Df125FDCPulse *p = fdcpulses[i];

      // Apply optional rocid translation
      uint32_t rocid = p->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;
      
      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << p->slot
               << " chan:" << p->channel << std::endl;
   
      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, p->slot, p->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
          if (VERBOSE > 6)
             ttout << "     - Didn't find it" << std::endl;
          continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6)
         ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys) 
               << std::endl;

      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
		case FDC_CATHODES:
            vfdccathode.push_back( MakeFDCCathodeDigiHit(chaninfo.fdc_cathodes, p) );
            break;
         default: 
             if (VERBOSE > 4)
                ttout << "       - Don't know how to make DigiHit"
                      << " objects for this detector type!" << std::endl; 
             break;
      }
   }

   // DF1TDCHit
   vector<const DF1TDCHit*> f1tdchits;
   loop->Get(f1tdchits);
   if (VERBOSE > 2)
      ttout << "  Number DF1TDCHit objects: " << f1tdchits.size() << std::endl;
   for (uint32_t i=0; i<f1tdchits.size(); i++) {
      const DF1TDCHit *hit = f1tdchits[i];

      // Apply optional rocid translation
      uint32_t rocid = hit->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;

      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << hit->slot
               << " chan:" << hit->channel << std::endl;

      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, hit->slot, hit->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
          if (VERBOSE > 6)
             ttout << "     - Didn't find it" << std::endl;
          continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6) 
         ttout << "     - Found entry for: " 
               << DetectorName(chaninfo.det_sys) << std::endl;
      
      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
         case BCAL:
            vbcaltdc.push_back( MakeBCALTDCDigiHit(chaninfo.bcal, hit) );
            break;
         case FDC_WIRES:
            vfdcwire.push_back( MakeFDCWireDigiHit(chaninfo.fdc_wires, hit) );
            break;
         case RF:
        	vrftdc.push_back( MakeRFTDCDigiTime(chaninfo.rf, hit) );
            break;
         case SC:
            vsctdc.push_back( MakeSCTDCDigiHit(chaninfo.sc, hit) );
            break;
         case TAGM:
            vtagmtdc.push_back( MakeTAGMTDCDigiHit(chaninfo.tagm, hit) );
            break;
         case TAGH:
            vtaghtdc.push_back( MakeTAGHTDCDigiHit(chaninfo.tagh, hit) );
            break;
         case PSC:
            vpsctdc.push_back( MakePSCTDCDigiHit(chaninfo.psc, hit) );
            break;
         default:     
             if (VERBOSE > 4)
                ttout << "       - Don't know how to make DigiHit objects"
                      << " for this detector type!" << std::endl;
             break;
      }
   }

   // DCAEN1290TDCHit
   vector<const DCAEN1290TDCHit*> caen1290tdchits;
   loop->Get(caen1290tdchits);
   if (VERBOSE > 2)
      ttout << "  Number DCAEN1290TDCHit objects: " 
            << caen1290tdchits.size() << std::endl;
   for (uint32_t i=0; i<caen1290tdchits.size(); i++) {
      const DCAEN1290TDCHit *hit = caen1290tdchits[i];

      // Apply optional rocid translation
      uint32_t rocid = hit->rocid;
      map<uint32_t, uint32_t>::iterator rocid_iter = Get_ROCID_Map().find(rocid);
      if (rocid_iter != Get_ROCID_Map().end()) rocid = rocid_iter->second;

      if (VERBOSE > 4)
         ttout << "    Looking for rocid:" << rocid << " slot:" << hit->slot
               << " chan:" << hit->channel << std::endl;
      
      // Create crate,slot,channel index and find entry in Translation table.
      // If none is found, then just quietly skip this hit.
      csc_t csc = {rocid, hit->slot, hit->channel};
      map<csc_t, DChannelInfo>::const_iterator iter = Get_TT().find(csc);
      if (iter == Get_TT().end()) {
          if (VERBOSE > 6)
             ttout << "     - Didn't find it" << std::endl;
          continue;
      }
      const DChannelInfo &chaninfo = iter->second;
      if (VERBOSE > 6)
         ttout << "     - Found entry for: " << DetectorName(chaninfo.det_sys)
               << std::endl;
      
      // Create the appropriate hit type based on detector type
      switch (chaninfo.det_sys) {
         case TOF:
            vtoftdc.push_back( MakeTOFTDCDigiHit(chaninfo.tof, hit) );
            break;
         case RF:
        	vrftdc.push_back( MakeRFTDCDigiTime(chaninfo.rf, hit) );
            break;
         default:     
             if (VERBOSE > 4)
                ttout << "       - Don't know how to make DigiHit objects"
                      << " for this detector type!" << std::endl;
             break;
      }
   }

   // Sort object order (this makes it easier to browse with hd_dump)
   sort(vbcal.begin(), vbcal.end(), SortBCALDigiHit);
   
   if (VERBOSE > 3) {
      ttout << "        vbcal.size() = " << vbcal.size() << std::endl;
      ttout << "     vbcaltdc.size() = " << vbcaltdc.size() << std::endl;
      ttout << "         vcdc.size() = " << vcdc.size() << std::endl;
      ttout << "        vfcal.size() = " << vfcal.size() << std::endl;
      ttout << "  vfdccathode.size() = " << vfdccathode.size() << std::endl;
      ttout << "     vfdcwire.size() = " << vfdcwire.size() << std::endl;
      ttout << "          vrf.size() = " << vrf.size() << std::endl;
      ttout << "       vrftdc.size() = " << vrftdc.size() << std::endl;
      ttout << "          vsc.size() = " << vsc.size() << std::endl;
      ttout << "       vsctdc.size() = " << vsctdc.size() << std::endl;
      ttout << "         vtof.size() = " << vtof.size() << std::endl;
      ttout << "      vtoftdc.size() = " << vtoftdc.size() << std::endl;
      ttout << "        vtagm.size() = " << vtagm.size() << std::endl;
      ttout << "     vtagmtdc.size() = " << vtagmtdc.size() << std::endl;
      ttout << "        vtagh.size() = " << vtagh.size() << std::endl;
      ttout << "     vtaghtdc.size() = " << vtaghtdc.size() << std::endl;
      ttout << "          vps.size() = " << vps.size() << std::endl;
      ttout << "         vpsc.size() = " << vpsc.size() << std::endl;
      ttout << "      vpsctdc.size() = " << vpsctdc.size() << std::endl;
      ttout << "  vtpolsector.size() = " << vtpolsector.size() << std::endl;
   }
   
   // Find factory for each container and copy the object pointers into it
   // (n.b. do this even if container is empty since it sets the evnt_called flag)
   CopyToFactory(loop, vbcal);
   CopyToFactory(loop, vbcaltdc);
   CopyToFactory(loop, vcdc);
   CopyToFactory(loop, vfcal);
   CopyToFactory(loop, vfdccathode);
   CopyToFactory(loop, vfdcwire);
   CopyToFactory(loop, vrf);
   CopyToFactory(loop, vrftdc);
   CopyToFactory(loop, vsc);
   CopyToFactory(loop, vsctdc);
   CopyToFactory(loop, vtof);
   CopyToFactory(loop, vtoftdc);
   CopyToFactory(loop, vtagm);
   CopyToFactory(loop, vtagmtdc);
   CopyToFactory(loop, vtagh);
   CopyToFactory(loop, vtaghtdc);
   CopyToFactory(loop, vps);
   CopyToFactory(loop, vpsc);
   CopyToFactory(loop, vpsctdc);
   CopyToFactory(loop, vtpolsector);

   // Add to JANA's call stack some entries to make janadot draw something reasonable
   // Unfortunately, this is just us telling JANA the relationship as defined here.
   // It is not derived from the above code which would guarantee the declared relationsips
   // are correct. That would just be too complicated given how that code works.
   if (record_call_stack) {
      // re-enable call stack recording
      loop->EnableCallStackRecording();

      Addf250ObjectsToCallStack(loop, "DBCALDigiHit");
      Addf250ObjectsToCallStack(loop, "DFCALDigiHit");
      Addf250ObjectsToCallStack(loop, "DSCDigiHit");
      Addf250ObjectsToCallStack(loop, "DTOFDigiHit");
      Addf125ObjectsToCallStack(loop, "DCDCDigiHit");
      Addf125ObjectsToCallStack(loop, "DFDCCathodeDigiHit");
      AddF1TDCObjectsToCallStack(loop, "DBCALTDCDigiHit");
      AddF1TDCObjectsToCallStack(loop, "DFDCWireDigiHit");
      AddF1TDCObjectsToCallStack(loop, "DRFDigiTime");
      AddF1TDCObjectsToCallStack(loop, "DRFTDCDigiTime");
      AddF1TDCObjectsToCallStack(loop, "DSCTDCDigiHit");
      AddCAEN1290TDCObjectsToCallStack(loop, "DTOFTDCDigiHit");
   }
}

//---------------------------------
// MakeBCALDigiHit
//---------------------------------
DBCALDigiHit* DTranslationTable::MakeBCALDigiHit(const BCALIndex_t &idx,
                                                 const Df250PulseIntegral *pi,
                                                 const Df250PulseTime *pt,
                                                 const Df250PulsePedestal *pp) const
{
   if (VERBOSE > 4)
      ttout << "       - Making DBCALDigiHit for (mod,lay,sec,end)=("
            << idx.module << "," << idx.layer << "," << idx.sector 
            << "," << (DBCALGeometry::End)idx.end << std::endl;

   DBCALDigiHit *h = new DBCALDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->pulse_peak = pp==NULL ? 0 : pp->pulse_peak; // Include pulse peak information in the digihit for BCAL

   h->module = idx.module;
   h->layer  = idx.layer;
   h->sector = idx.sector;
   h->end    = (DBCALGeometry::End)idx.end;

   return h;
}

//---------------------------------
// MakeFCALDigiHit
//---------------------------------
DFCALDigiHit* DTranslationTable::MakeFCALDigiHit(const FCALIndex_t &idx,
                                                 const Df250PulseIntegral *pi,
                                                 const Df250PulseTime *pt,
                                                 const Df250PulsePedestal *pp) const
{
   DFCALDigiHit *h = new DFCALDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->row    = idx.row;
   h->column = idx.col;
   
   return h;
}

//---------------------------------
// MakeTOFDigiHit
//---------------------------------
DTOFDigiHit* DTranslationTable::MakeTOFDigiHit(const TOFIndex_t &idx,
                                               const Df250PulseIntegral *pi,
                                               const Df250PulseTime *pt,
                                               const Df250PulsePedestal *pp) const
{
   DTOFDigiHit *h = new DTOFDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->plane = idx.plane;
   h->bar   = idx.bar;
   h->end   = idx.end;

   return h;
}

//---------------------------------
// MakeSCDigiHit
//---------------------------------
DSCDigiHit* DTranslationTable::MakeSCDigiHit(const SCIndex_t &idx, 
                                             const Df250PulseIntegral *pi,
                                             const Df250PulseTime *pt,
                                             const Df250PulsePedestal *pp) const
{
   DSCDigiHit *h = new DSCDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->sector = idx.sector;

   return h;
}

//---------------------------------
// MakeTAGMDigiHit
//---------------------------------
DTAGMDigiHit* DTranslationTable::MakeTAGMDigiHit(const TAGMIndex_t &idx,
                                                 const Df250PulseIntegral *pi,
                                                 const Df250PulseTime *pt,
                                                 const Df250PulsePedestal *pp) const
{
   DTAGMDigiHit *h = new DTAGMDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->row = idx.row;
   h->column = idx.col;

   return h;
}

//---------------------------------
// MakeTAGHDigiHit
//---------------------------------
DTAGHDigiHit* DTranslationTable::MakeTAGHDigiHit(const TAGHIndex_t &idx,
                                                 const Df250PulseIntegral *pi,
                                                 const Df250PulseTime *pt,
                                                 const Df250PulsePedestal *pp) const
{
   DTAGHDigiHit *h = new DTAGHDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->counter_id = idx.id;

   return h;
}

//---------------------------------
// MakePSCDigiHit
//---------------------------------
DPSCDigiHit* DTranslationTable::MakePSCDigiHit(const PSCIndex_t &idx,
					       const Df250PulseIntegral *pi,
					       const Df250PulseTime *pt,
					       const Df250PulsePedestal *pp) const
{
   DPSCDigiHit *h = new DPSCDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->counter_id = idx.id;

   return h;
}

//---------------------------------
// MakePSDigiHit
//---------------------------------
DPSDigiHit* DTranslationTable::MakePSDigiHit(const PSIndex_t &idx,
					     const Df250PulseIntegral *pi,
					     const Df250PulseTime *pt,
					     const Df250PulsePedestal *pp) const
{
   DPSDigiHit *h = new DPSDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->arm = (DPSGeometry::Arm)idx.side;
   h->column = idx.id;

   return h;
}

//---------------------------------
// MakeCDCDigiHit
//---------------------------------
DCDCDigiHit* DTranslationTable::MakeCDCDigiHit(const CDCIndex_t &idx,
                                               const Df125PulseIntegral *pi,
                                               const Df125PulseTime *pt,
                                               const Df125PulsePedestal *pp) const
{
   DCDCDigiHit *h = new DCDCDigiHit();
   CopyDf125Info(h, pi, pt, pp);

   h->ring = idx.ring;
   h->straw = idx.straw;
   
   return h;
}

//---------------------------------
// MakeCDCDigiHit
//---------------------------------
DCDCDigiHit* DTranslationTable::MakeCDCDigiHit(const CDCIndex_t &idx,
                                               const Df125CDCPulse *p) const
{
	DCDCDigiHit *h = new DCDCDigiHit();
	h->ring              = idx.ring;
	h->straw             = idx.straw;
	h->pulse_integral    = p->integral;
	h->pulse_time        = p->le_time;
	h->pedestal          = p->pedestal;
	h->QF                = p->time_quality_bit + (p->overflow_count<<1);
	h->nsamples_integral = p->nsamples_integral;
	h->nsamples_pedestal = p->nsamples_pedestal;

	h->AddAssociatedObject(p);

	return h;
}

//---------------------------------
// MakeFDCCathodeDigiHit
//---------------------------------
DFDCCathodeDigiHit* DTranslationTable::MakeFDCCathodeDigiHit(
                                       const FDC_CathodesIndex_t &idx,
                                       const Df125PulseIntegral *pi,
                                       const Df125PulseTime *pt,
                                       const Df125PulsePedestal *pp) const
{
   DFDCCathodeDigiHit *h = new DFDCCathodeDigiHit();
   CopyDf125Info(h, pi, pt, pp);

   h->package    = idx.package;
   h->chamber    = idx.chamber;
   h->view       = idx.view;
   h->strip      = idx.strip;
   h->strip_type = idx.strip_type;

   return h;
}


//---------------------------------
// MakeFDCCathodeDigiHit
//---------------------------------
DFDCCathodeDigiHit* DTranslationTable::MakeFDCCathodeDigiHit(
                                       const FDC_CathodesIndex_t &idx,
                                       const Df125FDCPulse *p) const
{
	DFDCCathodeDigiHit *h = new DFDCCathodeDigiHit();
	h->package           = idx.package;
	h->chamber           = idx.chamber;
	h->view              = idx.view;
	h->strip             = idx.strip;
	h->strip_type        = idx.strip_type;
	h->pulse_integral    = p->integral;
	h->pulse_time        = p->le_time;
	h->pedestal          = p->pedestal;
	h->QF                = p->time_quality_bit + (p->overflow_count<<1);
	h->nsamples_integral = p->nsamples_integral;
	h->nsamples_pedestal = p->nsamples_pedestal;

	h->AddAssociatedObject(p);

	return h;
}

//---------------------------------
// MakeBCALTDCDigiHit
//---------------------------------
DBCALTDCDigiHit* DTranslationTable::MakeBCALTDCDigiHit(
                                    const BCALIndex_t &idx,
                                    const DF1TDCHit *hit) const
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
DFDCWireDigiHit* DTranslationTable::MakeFDCWireDigiHit(
                                    const FDC_WiresIndex_t &idx,
                                    const DF1TDCHit *hit) const
{
   DFDCWireDigiHit *h = new DFDCWireDigiHit();
   CopyDF1TDCInfo(h, hit);

   h->package = idx.package;
   h->chamber = idx.chamber;
   h->wire    = idx.wire;

   return h;
}

//---------------------------------
// MakeRFDigiTime
//---------------------------------
DRFTDCDigiTime*  DTranslationTable::MakeRFTDCDigiTime(
                                   const RFIndex_t &idx,
                                   const DF1TDCHit *hit) const
{
   DRFTDCDigiTime *h = new DRFTDCDigiTime();
   CopyDF1TDCInfo(h, hit);

   h->dSystem = idx.dSystem;
   h->dIsCAENTDCFlag = false;

   return h;
}

//---------------------------------
// MakeRFDigiTime
//---------------------------------
DRFTDCDigiTime*  DTranslationTable::MakeRFTDCDigiTime(
                                   const RFIndex_t &idx,
                                   const DCAEN1290TDCHit *hit) const
{
   DRFTDCDigiTime *h = new DRFTDCDigiTime();
   CopyDCAEN1290TDCInfo(h, hit);

   h->dSystem = idx.dSystem;
   h->dIsCAENTDCFlag = true;

   return h;
}

//---------------------------------
// MakeRFDigiTime
//---------------------------------
DRFDigiTime*  DTranslationTable::MakeRFDigiTime(
                                   const RFIndex_t &idx,
                                   const Df250PulseTime *hit) const
{
   DRFDigiTime *h = new DRFDigiTime();
   h->time = hit->time;

   h->dSystem = idx.dSystem;

   return h;
}

//---------------------------------
// MakeSCTDCDigiHit
//---------------------------------
DSCTDCDigiHit*  DTranslationTable::MakeSCTDCDigiHit(
                                   const SCIndex_t &idx,
                                   const DF1TDCHit *hit) const
{
   DSCTDCDigiHit *h = new DSCTDCDigiHit();
   CopyDF1TDCInfo(h, hit);

   h->sector = idx.sector;

   return h;
}

//---------------------------------
// MakeTAGMTDCDigiHit
//---------------------------------
DTAGMTDCDigiHit*  DTranslationTable::MakeTAGMTDCDigiHit(
                                     const TAGMIndex_t &idx,
                                     const DF1TDCHit *hit) const
{
   DTAGMTDCDigiHit *h = new DTAGMTDCDigiHit();
   CopyDF1TDCInfo(h, hit);

   h->row = idx.row;
   h->column = idx.col;

   return h;
}

//---------------------------------
// MakeTAGHTDCDigiHit
//---------------------------------
DTAGHTDCDigiHit*  DTranslationTable::MakeTAGHTDCDigiHit(
                                     const TAGHIndex_t &idx,
                                     const DF1TDCHit *hit) const
{
   DTAGHTDCDigiHit *h = new DTAGHTDCDigiHit();
   CopyDF1TDCInfo(h, hit);

   h->counter_id = idx.id;

   return h;
}

//---------------------------------
// MakePSCTDCDigiHit
//---------------------------------
DPSCTDCDigiHit*  DTranslationTable::MakePSCTDCDigiHit(
                                     const PSCIndex_t &idx,
                                     const DF1TDCHit *hit) const
{
   DPSCTDCDigiHit *h = new DPSCTDCDigiHit();
   CopyDF1TDCInfo(h, hit);

   h->counter_id = idx.id;

   return h;
}

//---------------------------------
// MakeTOFTDCDigiHit
//---------------------------------
DTOFTDCDigiHit*  DTranslationTable::MakeTOFTDCDigiHit(
                                    const TOFIndex_t &idx,
                                    const DCAEN1290TDCHit *hit) const
{
   DTOFTDCDigiHit *h = new DTOFTDCDigiHit();
   CopyDCAEN1290TDCInfo(h, hit);

   h->plane = idx.plane;
   h->bar   = idx.bar;
   h->end   = idx.end;

   return h;
}

//---------------------------------
// MakeTPOLSectorDigiHit
//---------------------------------
DTPOLSectorDigiHit* DTranslationTable::MakeTPOLSectorDigiHit(const TPOLSECTORIndex_t &idx,
							     const Df250PulseIntegral *pi,
							     const Df250PulseTime *pt,
							     const Df250PulsePedestal *pp) const
{
   DTPOLSectorDigiHit *h = new DTPOLSectorDigiHit();
   CopyDf250Info(h, pi, pt, pp);

   h->sector = idx.sector;
   return h;
}

//---------------------------------
// GetDetectorIndex
//---------------------------------
const DTranslationTable::DChannelInfo 
     &DTranslationTable::GetDetectorIndex(const csc_t &in_daq_index) const
{
    map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>::const_iterator detector_index_itr = Get_TT().find(in_daq_index);
    if (detector_index_itr == Get_TT().end()) {
       stringstream ss_err;
       ss_err << "Could not find detector channel in Translaton Table: "
              << "rocid = " << in_daq_index.rocid
              << "slot = " << in_daq_index.slot
              << "channel = " << in_daq_index.channel;
       throw JException(ss_err.str());
    } 

    return detector_index_itr->second;
}

//---------------------------------
// GetDAQIndex
//---------------------------------
const DTranslationTable::csc_t 
     &DTranslationTable::GetDAQIndex(const DChannelInfo &in_channel) const
{
    map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>::const_iterator tt_itr = Get_TT().begin();

    // search through the whole Table to find the key that corresponds to our detector channel
    // this is not terribly efficient - linear in the size of the table
    bool found = false;
    for (; tt_itr != Get_TT().end(); tt_itr++) {
       const DTranslationTable::DChannelInfo &det_channel = tt_itr->second;
       if ( det_channel.det_sys == in_channel.det_sys ) {
          switch ( in_channel.det_sys ) {
          case DTranslationTable::BCAL:
             if ( det_channel.bcal == in_channel.bcal ) 
                found = true;
             break;
          case DTranslationTable::CDC:
             if ( det_channel.cdc == in_channel.cdc ) 
                found = true;
             break;
          case DTranslationTable::FCAL:
             if ( det_channel.fcal == in_channel.fcal ) 
                found = true;
             break;
          case DTranslationTable::FDC_CATHODES:
             if ( det_channel.fdc_cathodes == in_channel.fdc_cathodes ) 
                found = true;
             break;
          case DTranslationTable::FDC_WIRES:
             if ( det_channel.fdc_wires == in_channel.fdc_wires ) 
                found = true;
             break;
          case DTranslationTable::PS:
             if ( det_channel.ps == in_channel.ps ) 
                found = true;
            break;
          case DTranslationTable::PSC:
             if ( det_channel.psc == in_channel.psc )
                found = true;
             break;
          case DTranslationTable::RF:
             if ( det_channel.rf == in_channel.rf )
                found = true;
             break;
          case DTranslationTable::SC:
             if ( det_channel.sc == in_channel.sc )
                found = true;
             break;
          case DTranslationTable::TAGH:
             if ( det_channel.tagh == in_channel.tagh )
                found = true;
             break;
          case DTranslationTable::TAGM:
             if ( det_channel.tagm == in_channel.tagm )
                found = true;
             break;
          case DTranslationTable::TOF:
             if ( det_channel.tof == in_channel.tof )
                found = true;
             break;
          case DTranslationTable::TPOLSECTOR:
             if ( det_channel.tpolsector == in_channel.tpolsector )
                found = true;
             break;
          default:
             jerr << "DTranslationTable::GetDAQIndex(): "
                  << "Invalid detector type = " << in_channel.det_sys 
                  << std::endl;
       }
   }

   if (found)
       break;
    }
    
    if (tt_itr == Get_TT().end()) {
       stringstream ss_err;
       ss_err << "Could not find DAQ channel in Translaton Table:  "
              << Channel2Str(in_channel) << std::endl;
       throw JException(ss_err.str());
    }

    return tt_itr->first;
}

//----------------
// Channel2Str
//----------------
string DTranslationTable::Channel2Str(const DChannelInfo &in_channel) const
{
    stringstream ss;
    
    switch ( in_channel.det_sys ) {
    case DTranslationTable::BCAL:
       ss << "module = " << in_channel.bcal.module << " layer = " << in_channel.bcal.layer 
          << " sector = " << in_channel.bcal.sector << " end = " << in_channel.bcal.end;
       break;
    case DTranslationTable::CDC:
       ss << "ring = " << in_channel.cdc.ring << " straw = " << in_channel.cdc.straw;
       break;
    case DTranslationTable::FCAL:
       ss << "row = " << in_channel.fcal.row << " column = " << in_channel.fcal.col;
       break;
    case DTranslationTable::FDC_CATHODES:
       ss << "package = " << in_channel.fdc_cathodes.package
          << " chamber = " << in_channel.fdc_cathodes.chamber
          << " view = " << in_channel.fdc_cathodes.view
          << " strip = " << in_channel.fdc_cathodes.strip 
          << " strip type = " << in_channel.fdc_cathodes.strip_type;
       break;
    case DTranslationTable::FDC_WIRES:
       ss << "package = " << in_channel.fdc_wires.package
          << " chamber = " << in_channel.fdc_wires.chamber
          << " wire = " << in_channel.fdc_wires.wire;
       break;
    case DTranslationTable::PS:
       ss << "side = " << in_channel.ps.side << " id = " << in_channel.ps.id;
       break;
    case DTranslationTable::PSC:
       ss << "id = " << in_channel.psc.id;
       break;
    case DTranslationTable::RF:
       ss << "system = " << SystemName(in_channel.rf.dSystem);
       break;
    case DTranslationTable::SC:
       ss << "sector = " << in_channel.sc.sector;
       break;
    case DTranslationTable::TAGH:
       ss << "id = " << in_channel.tagh.id;
       break;
    case DTranslationTable::TAGM:
       ss << "row = " << in_channel.tagm.row << " column = " << in_channel.tagm.col;
       break;
    case DTranslationTable::TOF:
       ss << "plane = " << in_channel.tof.plane << " bar = " << in_channel.tof.bar
          << " end = " << in_channel.tof.end;
       break;
    case DTranslationTable::TPOLSECTOR:
       ss << "sector = " << in_channel.tpolsector.sector;
       break;
    default:
       ss << "Unknown detector type" << std::endl;
    }   

    return ss.str();
}

//----------------
// Addf250ObjectsToCallStack
//----------------
void DTranslationTable::Addf250ObjectsToCallStack(JEventLoop *loop, string caller) const
{
	AddToCallStack(loop, caller, "Df250Config");
	AddToCallStack(loop, caller, "Df250PulseIntegral");
	AddToCallStack(loop, caller, "Df250PulsePedestal");
	AddToCallStack(loop, caller, "Df250PulseTime");
}

//----------------
// Addf125ObjectsToCallStack
//----------------
void DTranslationTable::Addf125ObjectsToCallStack(JEventLoop *loop, string caller) const
{
	AddToCallStack(loop, caller, "Df125Config");
	AddToCallStack(loop, caller, "Df125PulseIntegral");
	AddToCallStack(loop, caller, "Df125PulsePedestal");
	AddToCallStack(loop, caller, "Df125PulseTime");
}

//----------------
// AddF1TDCObjectsToCallStack
//----------------
void DTranslationTable::AddF1TDCObjectsToCallStack(JEventLoop *loop, string caller) const
{
	AddToCallStack(loop, caller, "DF1TDCConfig");
	AddToCallStack(loop, caller, "DF1TDCHit");
}

//----------------
// AddCAEN1290TDCObjectsToCallStack
//----------------
void DTranslationTable::AddCAEN1290TDCObjectsToCallStack(JEventLoop *loop, string caller) const
{
	AddToCallStack(loop, caller, "DCAEN1290TDCConfig");
	AddToCallStack(loop, caller, "DCAEN1290TDCHit");
}

//----------------
// AddToCallStack
//----------------
void DTranslationTable::AddToCallStack(JEventLoop *loop, 
                                       string caller, string callee) const
{
   /// This is used to give information to JANA regarding the relationship and
   /// origin of some of these data objects. This is really just needed so that
   /// the janadot program can be used to produce the correct callgraph. Because
   /// of how this plugin works, JANA can't record the correct call stack (at
   /// least not easily!) Therefore, we have to give it a little help here.

   JEventLoop::call_stack_t cs;
   cs.start_time = cs.end_time = 0.0;
   cs.caller_name = caller;
   cs.callee_name = callee;
   cs.data_source = JEventLoop::DATA_FROM_CACHE;
   loop->AddToCallStack(cs);
   cs.callee_name = cs.caller_name;
   cs.caller_name = "<ignore>";
   cs.data_source = JEventLoop::DATA_FROM_FACTORY;
   loop->AddToCallStack(cs);
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
   pthread_mutex_lock(&Get_TT_Mutex());
   if (Get_TT_Initialized()) {
      pthread_mutex_unlock(&Get_TT_Mutex());
      return;
   }

   // String to hold entire XML translation table
   string tt_xml; 

   // Try getting it from CCDB first
   if (jcalib && !NO_CCDB) {
      map<string,string> tt;
      string namepath = "Translation/DAQ2detector";
      jout << "Reading translation table from calib DB: " << namepath << " ..." << std::endl;
      jcalib->GetCalib(namepath, tt);
      if (tt.size() != 1) {
         jerr << " Error: Unexpected translation table format!" << std::endl;
         jerr << "        tt.size()=" << tt.size() << " (expected 1)" << std::endl;
      }else{
         // Copy table into tt string
         map<string,string>::iterator iter = tt.begin();
         tt_xml = iter->second;
      }
   }
   
   // If getting from CCDB fails, try just reading in local file
   if (tt_xml.size() == 0) {
      if (!NO_CCDB) jout << "Unable to get translation table from CCDB." << std::endl;
      jout << "Will try reading TT from local file: " << XML_FILENAME << std::endl;

      // Open file
      ifstream ifs(XML_FILENAME.c_str());
      if (! ifs.is_open()) {
         jerr << " Error: Cannot open file! Translation table unavailable." << std::endl;
         pthread_mutex_unlock(&Get_TT_Mutex());
         return;
      }
      
      // read lines into stringstream object 
      stringstream ss;
      while (ifs.good()) {
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
   if (xmlParser == NULL) {
      jerr << "readTranslationTable...unable to create parser" << std::endl;
      exit(EXIT_FAILURE);
   }
   XML_SetElementHandler(xmlParser,StartElement,EndElement);
   XML_SetUserData(xmlParser, &Get_TT());

   // Parse XML string
   int status=XML_Parse(xmlParser, tt_xml.c_str(), tt_xml.size(), 1); // "1" indicates this is the final piece of XML
   if (status == 0) {
      jerr << "  ?readTranslationTable...parseXMLFile parse error for " << XML_FILENAME << std::endl;
      jerr << XML_ErrorString(XML_GetErrorCode(xmlParser)) << std::endl;
   }
   
   jout << Get_TT().size() << " channels defined in translation table" << std::endl;
   XML_ParserFree(xmlParser);

   pthread_mutex_unlock(&Get_TT_Mutex());
   Get_TT_Initialized() = true;
}

//---------------------------------
// ModuleStr2ModID
//---------------------------------
int ModuleStr2ModID(string &type)
{
   if (type == "vmecpu") {
      return(DModuleType::VMECPU);
   } else if (type == "tid") {
      return(DModuleType::TID);
   } else if (type == "fadc250") {
      return(DModuleType::FADC250);
   } else if (type == "fadc125") {
      return(DModuleType::FADC125);
   } else if (type == "f1tdcv2") {
      return(DModuleType::F1TDC32);
   } else if (type == "f1tdcv3") {
      return(DModuleType::F1TDC48);
   } else if (type == "jldisc") {
      return(DModuleType::JLAB_DISC);
   } else if (type == "vx1290a") {
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
   if ( type == "fdc_cathodes" ) {
      return DTranslationTable::FDC_CATHODES;
   } else if ( type == "fdc_wires" ) {
      return DTranslationTable::FDC_WIRES;   
   } else if ( type == "bcal" ) {
      return DTranslationTable::BCAL;
   } else if ( type == "cdc" ) {
      return DTranslationTable::CDC;   
   } else if ( type == "fcal" ) {
      return DTranslationTable::FCAL;   
   } else if ( type == "ps" ) {
      return DTranslationTable::PS;
   } else if ( type == "psc" ) {
      return DTranslationTable::PSC;
   } else if ( type == "rf" ) {
      return DTranslationTable::RF;
   } else if ( type == "st" ) {
           // The start counter is labelled by "ST" in the translation table
           // but we stick with the "SC" label in this plugin for consistency
           // with the rest of the reconstruction software
      return DTranslationTable::SC;
   } else if ( type == "tagh" ) {
      return DTranslationTable::TAGH;
   } else if ( type == "tagm" ) {
      return DTranslationTable::TAGM;
   } else if ( type == "tof" ) {
      return DTranslationTable::TOF;
   } else if ( type == "tpol" ) {
      return DTranslationTable::TPOLSECTOR;
   } else {
      return DTranslationTable::UNKNOWN_DETECTOR;
   }
}

//---------------------------------
// StartElement
//---------------------------------
void StartElement(void *userData, const char *xmlname, const char **atts)
{
   static int crate=0, slot=0;
   
   static string type,Type;
   int mc2codaType= 0;
   int channel = 0;
   string Detector, locSystem;
   int end=0;
   int row=0,column=0,module=0,sector=0,layer=0;
   int ring=0,straw=0,plane=0,bar=0;
   int package=0,chamber=0,view=0,strip=0,wire=0;
   int id=0, strip_type=0;
   int side=0;

   // This complicated line just recasts the userData pointer into
   // a reference to the "TT" member of the DTranslationTable object
   // that called us.
   map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo> &TT = *((map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>*)userData);
   
   // store crate summary info, fill both maps
   if (strcasecmp(xmlname,"halld_online_translation_table") == 0) {
      // do nothing
      
   } else if (strcasecmp(xmlname,"crate") == 0) {
      for (int i=0; atts[i]; i+=2) {
         if (strcasecmp(atts[i],"number") == 0) {
            crate = atoi(atts[i+1]);
            break;
         }
      }
      
   } else if (strcasecmp(xmlname,"slot") == 0) {
      for (int i=0; atts[i]; i+=2) {
         if (strcasecmp(atts[i],"number") == 0) {
            slot = atoi(atts[i+1]);
         } else if (strcasecmp(atts[i],"type") == 0) {
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
      
   } else if (strcasecmp(xmlname,"channel") == 0) {
      
      for (int i=0; atts[i]; i+=2) {
         string tag(atts[i+0]);
         string sval(atts[i+1]);
         int ival = atoi(atts[i+1]);

         if (tag == "number")
            channel = ival;
         else if (tag == "detector")
            Detector = sval;
         else if (tag == "row")
            row = ival;
         else if (tag == "column")
            column = ival;
         else if (tag == "col")
            column = ival;
         else if (tag == "module")
            module = ival;
         else if (tag == "sector")
            sector = ival;
         else if (tag == "layer")
            layer = ival;
         else if (tag == "chan");
//            chan = ival;
         else if (tag == "ring")
            ring = ival;
         else if (tag == "straw")
            straw = ival;
         else if (tag == "gPlane");
//            gPlane = ival;
         else if (tag == "element");
//            element = ival;
         else if (tag == "plane")
            plane = ival;
         else if (tag == "bar")
            bar = ival;
         else if (tag == "package")
            package = ival;
         else if (tag == "chamber")
            chamber = ival;
         else if (tag == "view") {
            if (sval == "U")
               view=1;
            else if (sval == "D")
               view=3;
         }
         else if (tag == "strip")
            strip = ival;
         else if (tag == "wire")
            wire = ival;
         else if (tag == "side") { 
            if (sval == "A") {
               side = DPSGeometry::kNorth;
            }
            else if (sval == "B") {
               side = DPSGeometry::kSouth;
            }
	 }
         else if (tag == "id")
            id = ival;
         else if (tag == "end") {
            if (sval == "U") {
               end = DBCALGeometry::kUpstream;
               view=1;
            }
            else if (sval == "D") {
               end = DBCALGeometry::kDownstream;
               view=3;
            }
            else if (sval == "N")
               end = 0; // TOF north
            else if (sval == "S")
               end = 1; // TOF south
            else if (sval == "UP")
               end = 0; // TOF up
            else if (sval == "DW") 
               end = 1; // TOF down
         }
         else if (tag == "strip_type") {
            if (sval == "full")
               strip_type = 1;
            else if (sval == "A")
               strip_type = 2;
            else if (sval == "B")
               strip_type = 3;
         }
         else if (tag == "system")
        	 locSystem = sval;
      }
      
      // ignore certain module types
      if (type == "disc")
         return;
      if (type == "ctp")
         return;
      if (type == "sd")
         return;
      if (type == "a1535sn")
         return;

      
//      // Data integrity check
//      if (crate < 0 || crate >= MAXDCRATE) {
//         jerr << " Crate value of " << crate 
//              << " is not in range 0 <= crate < " << MAXDCRATE << std::endl;
//         exit(-1);
//      }
//      
//      if (slot < 0 || slot >= MAXDSLOT) {
//         jerr << " Slot value of " << slot 
//              << " is not in range 0 <= slot < " << MAXDSLOT << std::endl;
//         exit(-1);
//      }
//      
//      if (channel < 0 || channel >= MAXDCHANNEL) {
//         jerr << " Crate value of " << channel 
//              << " is not in range 0 <= channel < " << MAXDCHANNEL << std::endl;
//         exit(-1);
//      }
      
      // fill maps
      
      DTranslationTable::csc_t csc = {(uint32_t)crate,(uint32_t)slot,(uint32_t)channel};
      string detector = Detector;
      std::transform(detector.begin(), detector.end(), detector.begin(), (int(*)(int)) tolower);
      
      //string s="unknown::";

      // Common indexes
      DTranslationTable::DChannelInfo &ci = TT[csc];
      ci.CSC = csc;
      ci.module_type = (DModuleType::type_id_t)mc2codaType;
      ci.det_sys = DetectorStr2DetID(detector);
      DTranslationTable::Get_ROCID_By_System()[ci.det_sys].insert(crate);

      // detector-specific indexes
      switch (ci.det_sys) {
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
         case DTranslationTable::RF:
            ci.rf.dSystem = NameToSystem(locSystem.c_str());
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
	        ci.ps.side = side;
	        ci.ps.id = id;
	        break;
         case DTranslationTable::PSC:
	        ci.psc.id = id;
            break;
         case DTranslationTable::TPOLSECTOR:
            ci.tpolsector.sector = sector;
            break;
         case DTranslationTable::UNKNOWN_DETECTOR:
		 default:
            break;
      }

   } else {
      jerr << std::endl << std::endl 
           << "?startElement...unknown xml tag " << xmlname 
           << std::endl << std::endl;
   }
   
}


//--------------------------------------------------------------------------


void EndElement(void *userData, const char *xmlname) {
   // nothing to do yet...
}


//--------------------------------------------------------------------------
