#ifndef _DFDCPseudo_factory_CORRECTED_ 
#define _DFDCPseudo_factory_CORRECTED_ 

#include <JANA/JFactory.h>
using namespace jana;

#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCSegment.h>

class DFDCPseudo_factory_CORRECTED:public JFactory<DFDCPseudo>{
	public:
		DFDCPseudo_factory_CORRECTED(){};
		~DFDCPseudo_factory_CORRECTED(){}; 
		const char* Tag(void){return "CORRECTED";}

		// The values of this factory are filled from the
		// DFDCSegment factory. If we happen to be called before
		// they are (which is necessarily the case when evnt
		// is called here) then activate the DFDCSegment factory
		// to fill our _data vector.
		jerror_t evnt(JEventLoop *loop, int eventNo)
		{
			vector<const DFDCSegment*> segments;
			loop->Get(segments);
			return NOERROR;
		}

	private:

};

#endif // _DFDCPseudo_factory_CORRECTED_

