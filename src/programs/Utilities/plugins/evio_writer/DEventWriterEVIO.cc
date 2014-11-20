#include "DEventWriterEVIO.h"

//file-scope so shared amongst threads; only accessed below via locks
size_t gEVIONumOutputThreads = 0;
map<string, int>* gEVIOOutputFilePointers = NULL;

DEventWriterEVIO::DEventWriterEVIO(JEventLoop* locEventLoop)
{
	japp->WriteLock("EVIOWriter");
	{
		++gEVIONumOutputThreads;
		if(gEVIOOutputFilePointers == NULL)
			gEVIOOutputFilePointers = new map<string, int>();
	}
	japp->Unlock("EVIOWriter");
}

bool DEventWriterEVIO::Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
	JEvent& locEvent = locEventLoop->GetJEvent();

	// Get pointer to JEventSource and make sure it is the right type
	JEventSource* locEventSource = locEvent.GetJEventSource();
	if(locEventSource == NULL)
		return false;

#if HAVE_EVIO
	JEventSource_EVIO* locEvioSource = dynamic_cast<JEventSource_EVIO*>(locEventSource);
	if(locEvioSource == NULL)
	{
		jerr << "WARNING!!! You MUST use this only with EVIO formatted data!!!" << endl;
		return false;
	}

	//Get the EVIO buffer
	uint32_t* locEVIOBuffer;
	uint32_t locBufferSize;
	locEvioSource->GetEVIOBuffer(locEvent, locEVIOBuffer, locBufferSize);
	if(locEVIOBuffer == NULL)
		return false;

	// write the resulting record to the output stream
	return Write_EVIOEvent(locEventLoop, locOutputFileNameSubString, locEVIOBuffer);
#else
	return false;
#endif // HAVE_EVIO
}

string DEventWriterEVIO::Get_OutputFileName(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
	//get the event source
	JEvent& locEvent = locEventLoop->GetJEvent();
	JEventSource* locEventSource = locEvent.GetJEventSource();
	if(locEventSource == NULL)
		return "no_name.evio";

	//get the source file name (strip the path)
	string locSourceFileName = locEventSource->GetSourceName();
	size_t locSlashIndex = locSourceFileName.find_last_of("/");
	string locSourceFileName_Pathless = (locSlashIndex != string::npos) ? locSourceFileName.substr(locSlashIndex + 1) : locSourceFileName;

	//strip the file extension (if present and if is a known format: .root, .evio, or .hddm)
	size_t locDotIndex = locSourceFileName_Pathless.find_last_of(".");
	if(locDotIndex != string::npos)
	{
		string locSuffix = locSourceFileName_Pathless.substr(locDotIndex + 1);
		if((locSuffix == "root") || (locSuffix == "evio") || (locSuffix == "hddm"))
			locSourceFileName_Pathless = locSourceFileName_Pathless.substr(0, locDotIndex);
	}

	return (locSourceFileName_Pathless + string(".") + locOutputFileNameSubString + string(".evio"));
}

bool DEventWriterEVIO::Open_OutputFile(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
#if HAVE_EVIO
	string locOutputFileName = Get_OutputFileName(locEventLoop, locOutputFileNameSubString);

	int locStatus = 0;
	japp->WriteLock("EVIOWriter");
	{
		if(gEVIOOutputFilePointers->find(locOutputFileName) != gEVIOOutputFilePointers->end())
		{
			//already open, don't re-open
			japp->Unlock("EVIOWriter");
			return true;
		}

		//not open: open it
		int locEVIOHandle = 0;
		const char *locFileCharName = locOutputFileName.c_str();
		locStatus = evOpen((char*)locFileCharName, (char*)"w", &locEVIOHandle);
		if(locStatus != S_SUCCESS)
			jerr << "Unable to open EVIO file." << endl;
		else //store the handle
			(*gEVIOOutputFilePointers)[locOutputFileName] = locEVIOHandle;
	}
	japp->Unlock("EVIOWriter");

	return (locStatus == S_SUCCESS);
#else
	return false;
#endif // HAVE_EVIO
}

bool DEventWriterEVIO::Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString, uint32_t* locEVIOBuffer) const
{
#if HAVE_EVIO
	string locOutputFileName = Get_OutputFileName(locEventLoop, locOutputFileNameSubString);
	japp->WriteLock("EVIOWriter");
	{
		//check to see if the EVIO file is open
		if(gEVIOOutputFilePointers->find(locOutputFileName) == gEVIOOutputFilePointers->end())
		{
			//not open, open it
			if(!Open_OutputFile(locEventLoop, locOutputFileNameSubString))
				return false; //failed to open
		}

		//open: get handle, write event
		int locEVIOHandle = (*gEVIOOutputFilePointers)[locOutputFileName];
		evWrite(locEVIOHandle, locEVIOBuffer);
	}
	japp->Unlock("EVIOWriter");

	return true;
#else
	return false;
#endif // HAVE_EVIO
}

DEventWriterEVIO::~DEventWriterEVIO(void)
{
#if HAVE_EVIO
	japp->WriteLock("EVIOWriter");
	{
		--gEVIONumOutputThreads;
		if(gEVIONumOutputThreads > 0)
		{
			japp->Unlock("EVIOWriter");
			return; //not the last thread writing to EVIO files
		}

		if(gEVIOOutputFilePointers == NULL)
		{
			japp->Unlock("EVIOWriter");
			return; //not the last thread writing to EVIO files
		}
		
		//last thread writing to EVIO files: close all files and free all memory
		map<string, int>::iterator locIterator = gEVIOOutputFilePointers->begin();
		for(; locIterator != gEVIOOutputFilePointers->end(); ++locIterator)
		{
			string locOutputFileName = locIterator->first;
			int locEVIOHandle = locIterator->second;
			evClose(locEVIOHandle);
			std::cout << "Closed EVIO file " << locOutputFileName << std::endl;
		}
		delete gEVIOOutputFilePointers;
		gEVIOOutputFilePointers = NULL;
	}
	japp->Unlock("EVIOWriter");
#endif // HAVE_EVIO
}

