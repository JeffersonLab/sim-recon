// Smearing class for forward DIRC (FDIRC)
// needs to be filled in

#ifndef _FDIRCSMEARER_H_
#define _FDIRCSMEARER_H_

#include "Smearer.h"

//class fdirc_config_t; // forward definition for readability


class FDIRCSmearer : public Smearer
{
  public:
	FDIRCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		//fdirc_config = new fdirc_config_t(loop);
	}
	~FDIRCSmearer() {
		//delete fdirc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	//fdirc_config_t  *fdirc_config;
};


#endif // _FDIRCSMEARER_H_