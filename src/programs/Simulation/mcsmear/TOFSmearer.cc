#include "TOFSmearer.h"

//-----------
// tof_config_t  (constructor)
//-----------
tof_config_t::tof_config_t(JEventLoop *loop) 
{
	// default values
 	TOF_SIGMA = 100.*k_psec;
 	TOF_PHOTONS_PERMEV = 400.;
 	TOF_BAR_THRESHOLD    = 0.0;

	// Load data from CCDB
    cout<<"Get TOF/tof_parms parameters from CCDB..."<<endl;
    map<string, double> tofparms;
    if(loop->GetCalib("TOF/tof_parms", tofparms)) {
     	jerr << "Problem loading TOF/tof_parms from CCDB!" << endl;
     	return;
    }
     	
    TOF_SIGMA =  tofparms["TOF_SIGMA"];
    TOF_PHOTONS_PERMEV =  tofparms["TOF_PHOTONS_PERMEV"];
	
}


//-----------
// SmearEvent
//-----------
void TOFSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::FtofCounterList tofs = record->getFtofCounters();
   hddm_s::FtofCounterList::iterator iter;
   for (iter = tofs.begin(); iter != tofs.end(); ++iter) {
      // take care of hits
      iter->deleteFtofHits();
      hddm_s::FtofTruthHitList thits = iter->getFtofTruthHits();
      hddm_s::FtofTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // Smear the time
         double t = titer->getT() + SampleGaussian(tof_config->TOF_SIGMA);
         // Smear the energy
         double npe = titer->getDE() * 1000. * tof_config->TOF_PHOTONS_PERMEV;
         npe = npe +  SampleGaussian(sqrt(npe));
         float NewE = npe/tof_config->TOF_PHOTONS_PERMEV/1000.;
         if (NewE > tof_config->TOF_BAR_THRESHOLD) {
            hddm_s::FtofHitList hits = iter->addFtofHits();
            hits().setEnd(titer->getEnd());
            hits().setT(t);
            hits().setDE(NewE);
         }
      }
    
      if (config->DROP_TRUTH_HITS) {
         iter->deleteFtofTruthHits();
      }
   }
}