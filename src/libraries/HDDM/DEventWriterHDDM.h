#ifndef _DEventWriterHDDM_
#define _DEventWriterHDDM_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
//#include <JANA/JEventProcessor.h>
#include <HDDM/hddm_s.hpp>

#include <CDC/DCDCHit.h>
#include <TOF/DTOFHit.h>
#include <FCAL/DFCALHit.h>
#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <START_COUNTER/DSCHit.h>
#include <FDC/DFDCHit.h>
#include <PAIR_SPECTROMETER/DPSHit.h>
#include <PAIR_SPECTROMETER/DPSCHit.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>
#include <TPOL/DTPOLHit.h>
using namespace std;
using namespace jana;

class DEventWriterHDDM : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterHDDM);

		DEventWriterHDDM(JEventLoop* locEventLoop, string locOutputFileBaseName);
		~DEventWriterHDDM(void);

		bool Write_HDDMEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;
		string Get_OutputFileName(string locOutputFileNameSubString) const;

	private:
		bool Write_HDDMEvent(string locOutputFileName, hddm_s::HDDM& locRecord) const;

		//contains static variables shared amongst threads
		int& Get_NumEventWriterThreads(void) const; //acquire HDDMWriter lock before modifying
		map<string, pair<ofstream*, hddm_s::ostream*> >& Get_HDDMOutputFilePointers(void) const;

		int32_t Convert_UnsignedIntToSigned(uint32_t locUnsignedInt) const;

		string dOutputFileBaseName;
		bool HDDM_USE_COMPRESSION;
		bool HDDM_USE_INTEGRITY_CHECKS;

        // metadata to save in the HDDM file
        // these should be consistent during program execution
        string HDDM_DATA_VERSION_STRING;
        string CCDB_CONTEXT_STRING;
};

#endif //_DEventWriterHDDM_
