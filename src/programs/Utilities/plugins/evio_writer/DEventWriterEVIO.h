#ifndef _DEventWriterEVIO_
#define _DEventWriterEVIO_

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

using namespace std;
using namespace jana;

class DEventWriterEVIO : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterEVIO);

		DEventWriterEVIO(JEventLoop* locEventLoop);
		~DEventWriterEVIO(void);

		//It's not strictly necessary to call this prior to "Write_EVIOEvent," but it's preferred:
			//If no events are written to the output file, it will still exist
		//If the file is already open, this doesn't do anything
		bool Open_OutputFile(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

		//This method will also open the output file if it hasn't been opened yet. 
		bool Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

		string Get_OutputFileName(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

	private:
		bool Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString, uint32_t* locEVIOBuffer) const;
};

#endif //_DEventWriterEVIO_

