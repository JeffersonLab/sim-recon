#include "FCALSmearer.h"

//-----------
// fcal_config_t  (constructor)
//-----------
fcal_config_t::fcal_config_t(JEventLoop *loop) 
{
	// default values
	FCAL_PHOT_STAT_COEF   = 0.0; //0.035;
	FCAL_BLOCK_THRESHOLD  = 0.0; //20.0*k_MeV;
	// FCAL_TSIGMA           = 0.0; // 200 ps
	FCAL_TSIGMA           = 0.2; // 200 ps - FIX THIS!!

	// Get values from CCDB
	cout << "Get FCAL/fcal_parms parameters from CCDB..." << endl;
    map<string, double> fcalparms;
    if(loop->GetCalib("FCAL/fcal_parms", fcalparms)) { 
     	jerr << "Problem loading FCAL/fcal_parms from CCDB!" << endl;
    } else {
       	FCAL_PHOT_STAT_COEF   = fcalparms["FCAL_PHOT_STAT_COEF"]; 
       	FCAL_BLOCK_THRESHOLD  = fcalparms["FCAL_BLOCK_THRESHOLD"];
	}
		
	cout<<"get FCAL/gains from calibDB"<<endl;
    vector <double> FCAL_GAINS_TEMP;

    loop->GetCalib("FCAL/gains", FCAL_GAINS_TEMP);
    for (unsigned int i = 0; i < FCAL_GAINS_TEMP.size(); i++) {
       	FCAL_GAINS.push_back(FCAL_GAINS_TEMP.at(i));
    }
     
    cout<<"get FCAL/digi_scales parameters from calibDB"<<endl;
    map<string, double> fcaldigiscales;
    loop->GetCalib("FCAL/digi_scales", fcaldigiscales);
    FCAL_MC_ESCALE = fcaldigiscales["FCAL_ADC_ASCALE"];

}
	
//-----------
// SmearEvent
//-----------
void FCALSmearer::SmearEvent(hddm_s::HDDM *record)
{
   /// Smear the FCAL hits using the nominal resolution of the individual blocks.
   /// The way this works is a little funny and warrants a little explanation.
   /// The information coming from hdgeant is truth information indexed by 
   /// row and column, but containing energy deposited and time. The mcsmear
   /// program will copy the truth information from the fcalTruthHit element
   /// to a new fcalHit element, smearing the values with the appropriate detector
   /// resolution.
   ///
   /// To access the "truth" values in DANA, get the DFCALHit objects using the
   /// "TRUTH" tag.

   //if (!fcalGeom)
   //   fcalGeom = new DFCALGeometry();

   hddm_s::FcalBlockList blocks = record->getFcalBlocks();
   hddm_s::FcalBlockList::iterator iter;
   for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
      iter->deleteFcalHits();
      hddm_s::FcalTruthHitList thits = iter->getFcalTruthHits();
      hddm_s::FcalTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // Simulation simulates a grid of blocks for simplicity. 
         // Do not bother smearing inactive blocks. They will be
         // discarded in DEventSourceHDDM.cc while being read in
         // anyway.
         if (!fcalGeom->isBlockActive(iter->getRow(), iter->getColumn()))
            continue;
         
         // Get gain constant per block
         int channelnum = fcalGeom->channel(iter->getRow(), iter->getColumn()); 
         double FCAL_gain = fcal_config->FCAL_GAINS.at(channelnum);

         // Smear the energy and timing of the hit
         double sigma = fcal_config->FCAL_PHOT_STAT_COEF/sqrt(titer->getE());
              
         // Apply constant scale factor to MC eneregy. 06/22/2016 A. Subedi
         double E = fcal_config->FCAL_MC_ESCALE * titer->getE() * (1.0 + SampleGaussian(sigma)); 
         
         
         // Smear the time by 200 ps (fixed for now) 7/2/2009 DL
         double t = titer->getT() + SampleGaussian(fcal_config->FCAL_TSIGMA); 
         // Apply a single block threshold. 
         
         // Scale threshold by gains         
         if (E >= fcal_config->FCAL_BLOCK_THRESHOLD * FCAL_gain ){
               hddm_s::FcalHitList hits = iter->addFcalHits();
               hits().setE(E);
               hits().setT(t);
         }
        
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteFcalTruthHits();
   }
}
