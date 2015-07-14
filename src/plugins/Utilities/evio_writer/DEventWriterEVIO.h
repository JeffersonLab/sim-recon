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

		bool Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

		string Get_OutputFileName(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;

	private:
		bool Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString, uint32_t* locEVIOBuffer) const;
		bool Open_OutputFile(JEventLoop* locEventLoop, string locOutputFileName) const;
};

#endif //_DEventWriterEVIO_

