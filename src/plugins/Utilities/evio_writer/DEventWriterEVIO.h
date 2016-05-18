#ifndef _DEventWriterEVIO_
#define _DEventWriterEVIO_

#include <map>
#include <vector>
#include <string>

#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>

#if HAVE_EVIO
#include <DAQ/JEventSource_EVIO.h>
#include <evioUtil.hxx>
#include <evioFileChannel.hxx>
#endif // HAVE_EVIO

#include <pthread.h>
#include <stdint.h>
#include <fstream>

#include <JANA/JEventLoop.h>

#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/Df125WindowRawData.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/Df125Config.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCConfig.h>
#include <DAQ/DCAEN1290TDCConfig.h>
#include <DAQ/DCAEN1290TDCHit.h>
#include <DAQ/DEPICSvalue.h>
#include <DAQ/DEventTag.h>

#include <DANA/DStatusBits.h>
#include <TTAB/DTranslationTable.h>

#include "HDEVIOWriter.h"

using namespace std;
using namespace jana;

class DEventWriterEVIO : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterEVIO);

		DEventWriterEVIO(JEventLoop* locEventLoop);
		~DEventWriterEVIO(void);

		bool Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

		string Get_OutputFileName(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

        void SetDetectorsToWriteOut(string detector_list);

		bool COMPACT;
		bool PREFER_EMULATED;
		bool DEBUG_FILES;

	protected:
		void WriteEventToBuffer(JEventLoop *locEventLoop, vector<uint32_t> &buff) const;
		bool Open_OutputFile(JEventLoop* locEventLoop, string locOutputFileName) const;
		
		void WriteCAEN1290Data(vector<uint32_t> &buff,
                               vector<const DCAEN1290TDCHit*>    &caen1290hits,
                               vector<const DCAEN1290TDCConfig*> &caen1290configs, 
                               unsigned int Nevents) const;

		void WriteF1Data(vector<uint32_t> &buff,
                         vector<const DF1TDCHit*>          &F1hits,
                         vector<const DF1TDCTriggerTime*>  &F1tts,
                         vector<const DF1TDCConfig*>       &F1configs, 
                         unsigned int Nevents) const;

		void Writef250Data(vector<uint32_t> &buff,
                           vector<const Df250PulseIntegral*> &f250pis,
                           vector<const Df250TriggerTime*>   &f250tts,
                           vector<const Df250WindowRawData*> &f250wrds,
                           unsigned int Nevents) const;
        
		void Writef125Data(vector<uint32_t> &buff,
                           vector<const Df125PulseIntegral*> &f125pis,
                           vector<const Df125CDCPulse*>      &f125cdcpulses,
                           vector<const Df125FDCPulse*>      &f125fdcpulses,
                           vector<const Df125TriggerTime*>   &f125tts,
                           vector<const Df125WindowRawData*> &f125wrds,
                           vector<const Df125Config*>        &f125configs,
                           unsigned int Nevents) const;
        
		void WriteEPICSData(vector<uint32_t> &buff,
                            vector<const DEPICSvalue*> epicsValues) const;
        
		void WriteEventTagData(vector<uint32_t> &buff,
                               uint64_t event_status,
                               const DL3Trigger* l3trigger) const;

        void WriteBORData(JEventLoop *loop, 
                          vector<uint32_t> &buff) const;
        
		std::ofstream *ofs_debug_input;
		std::ofstream *ofs_debug_output;

        const DTranslationTable *ttab;
        bool write_out_all_rocs;
        set<uint32_t> rocs_to_write_out;

	private:

		//contain static variables shared amongst threads: acquire "EVIOWriter" write lock before calling
		size_t& Get_NumEVIOOutputThreads(void) const;
		map<string, HDEVIOWriter*>& Get_EVIOOutputters(void) const;
		map<string, pthread_t>& Get_EVIOOutputThreads(void) const;

};

#endif //_DEventWriterEVIO_

