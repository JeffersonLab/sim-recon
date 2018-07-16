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

	// geometry
	const int TOF_NUM_PLANES = 2;
    const int TOF_NUM_BARS = 44;

	// Load data from CCDB
    cout<<"Get TOF/tof_parms parameters from CCDB..."<<endl;
    map<string, double> tofparms;
    if(loop->GetCalib("TOF/tof_parms", tofparms)) {
     	jerr << "Problem loading TOF/tof_parms from CCDB!" << endl;
     	return;
    }
     	
    TOF_SIGMA =  tofparms["TOF_SIGMA"];
    TOF_PHOTONS_PERMEV =  tofparms["TOF_PHOTONS_PERMEV"];
	
	cout<<"get TOF/paddle_resolutions from calibDB"<<endl;
    vector <double> TOF_PADDLE_TIME_RESOLUTIONS_TEMP;
    if(loop->GetCalib("TOF/paddle_resolutions", TOF_PADDLE_TIME_RESOLUTIONS_TEMP)) {
    	jerr << "Problem loading TOF/paddle_resolutions from CCDB!" << endl;
    } else {
    	for (unsigned int i = 0; i < TOF_PADDLE_TIME_RESOLUTIONS_TEMP.size(); i++) {
       		TOF_PADDLE_TIME_RESOLUTIONS.push_back(TOF_PADDLE_TIME_RESOLUTIONS_TEMP.at(i));
    	}
    }

	// load per-channel efficiencies
	vector<double> raw_table;
	if(loop->GetCalib("TOF/channel_mc_efficiency", raw_table)) {
    	jerr << "Problem loading TOF/channel_mc_efficiency from CCDB!" << endl;
    } else {
    	int channel = 0;

    	for(int plane=0; plane<TOF_NUM_PLANES; plane++) {
        	int plane_index=2*TOF_NUM_BARS*plane;
        	channel_efficiencies.push_back( vector< pair<double,double> >(TOF_NUM_BARS) );
        	for(int bar=0; bar<TOF_NUM_BARS; bar++) {
            	channel_efficiencies[plane][bar] 
                	= pair<double,double>(raw_table[plane_index+bar],
                    	    raw_table[plane_index+TOF_NUM_BARS+bar]);
            	channel+=2;
            }	      
        }
    }
    
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
         // correct simulation efficiencies 
		 if (config->APPLY_EFFICIENCY_CORRECTIONS
		 		&& !gDRandom.DecideToAcceptHit(tof_config->GetEfficiencyCorrectionFactor(titer)))
		 			continue;
		 			
         // Smear the time
         //double t = titer->getT() + gDRandom.SampleGaussian(tof_config->TOF_SIGMA);
         double t = titer->getT();
         // Smear the energy
         float NewE = titer->getDE();
         if(config->SMEAR_HITS) {
			 t += gDRandom.SampleGaussian(tof_config->GetHitTimeResolution(iter->getPlane(),iter->getBar()));
         	 double npe = titer->getDE() * 1000. * tof_config->TOF_PHOTONS_PERMEV;
         	 npe += gDRandom.SampleGaussian(sqrt(npe));
          	 NewE = npe/tof_config->TOF_PHOTONS_PERMEV/1000.;
		 }
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