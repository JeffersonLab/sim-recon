
#include "hddm.h"

//----------------------
// init_hddm_banks_t
//----------------------
derror_t init_hddm_banks_t(hddm_banks_t *hddm)
{
	/// Call constructors for all DBank derived classes
	hddm->FCALclusters	= new FCALclusters_t();
	hddm->CDChits			= new CDChits_t();
	hddm->CDCtracks		= new CDCtracks_t();
}

//----------------------
// delete_hddm_banks_t
//----------------------
derror_t delete_hddm_banks_t(hddm_banks_t *hddm)
{
	/// Call destructors for all DBank derived classes
	delete hddm->FCALclusters;
	delete hddm->CDChits;
	delete hddm->CDCtracks;
}

