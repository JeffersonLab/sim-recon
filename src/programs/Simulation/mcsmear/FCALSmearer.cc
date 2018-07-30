#include "FCALSmearer.h"

//-----------
// fcal_config_t  (constructor)
//-----------
fcal_config_t::fcal_config_t(JEventLoop *loop, DFCALGeometry *fcalGeom) 
{
	// default values
	FCAL_PHOT_STAT_COEF   = 0.0; //0.035;
	FCAL_BLOCK_THRESHOLD  = 0.0; //20.0*k_MeV;
	// FCAL_TSIGMA           = 0.0; // 200 ps
	FCAL_TSIGMA           = 0.; // 400 ps 

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
    if(loop->GetCalib("FCAL/gains", FCAL_GAINS_TEMP)) {
    	jerr << "Problem loading FCAL/gains from CCDB!" << endl;
    } else {
    	for (unsigned int i = 0; i < FCAL_GAINS_TEMP.size(); i++) {
       		FCAL_GAINS.push_back(FCAL_GAINS_TEMP.at(i));
    	}
    }
     
    cout<<"get FCAL/digi_scales parameters from calibDB"<<endl;
    map<string, double> fcaldigiscales;
    if(loop->GetCalib("FCAL/digi_scales", fcaldigiscales)) {
    	jerr << "Problem loading FCAL/digi_scales from CCDB!" << endl;
    } else {
        FCAL_MC_ESCALE = fcaldigiscales["FCAL_ADC_ASCALE"];
    }

    cout<<"get FCAL/mc_timing_smear parameters from calibDB"<<endl;
    map<string, double> fcalmctimingsmear;
    if(loop->GetCalib("FCAL/mc_timing_smear", fcalmctimingsmear)) {
    	jerr << "Problem loading FCAL/mc_timing_smear from CCDB!" << endl;
    } else {
        FCAL_TSIGMA = fcalmctimingsmear["FCAL_TSIGMA"];
    }

	// initialize 2D matrix of efficiencies, indexed by (row,column)
	vector< vector<double > > new_block_efficiencies(DFCALGeometry::kBlocksTall, 
            vector<double>(DFCALGeometry::kBlocksWide));
	block_efficiencies = new_block_efficiencies;

	// load efficiencies from CCDB and fill 
	vector<double> raw_table;
	if(loop->GetCalib("FCAL/block_mc_efficiency", raw_table)) {
    	jerr << "Problem loading FCAL/block_mc_efficiency from CCDB!" << endl;
    } else {
		for (int channel=0; channel < static_cast<int>(raw_table.size()); channel++) {
    
        	// make sure that we don't try to load info for channels that don't exist
        	if (channel == fcalGeom->numActiveBlocks())
            	break;

        	int row = fcalGeom->row(channel);
        	int col = fcalGeom->column(channel);

        	// results from DFCALGeometry should be self consistent, but add in some
        	// sanity checking just to be sure
        	if (fcalGeom->isBlockActive(row,col) == false) {
        		char str[200];
            	sprintf(str, "Loading FCAL constant for inactive channel!  "
                	    "row=%d, col=%d", row, col);
            	throw JException(str);
        	}

	        block_efficiencies[row][col] = raw_table[channel];
    	}
    }

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
            
         // correct simulation efficiencies 
		 if (config->APPLY_EFFICIENCY_CORRECTIONS
             && !gDRandom.DecideToAcceptHit(fcal_config->GetEfficiencyCorrectionFactor(iter->getRow(), iter->getColumn()))) {
             continue;
         } 

         // Get gain constant per block
         int channelnum = fcalGeom->channel(iter->getRow(), iter->getColumn()); 
         double FCAL_gain = fcal_config->FCAL_GAINS.at(channelnum);
              
         double E = titer->getE();
         if(fcal_config->FCAL_ADD_LIGHTGUIDE_HITS) {
             hddm_s::FcalTruthLightGuideList lghits = titer->getFcalTruthLightGuides();
             hddm_s::FcalTruthLightGuideList::iterator lgiter;
             for (lgiter = lghits.begin(); lgiter != lghits.end(); lgiter++) {
                 E += lgiter->getE();
             }
         }
         // Apply constant scale factor to MC energy. 06/22/2016 A. Subedi
         E *= fcal_config->FCAL_MC_ESCALE; 
         
         double t = titer->getT(); 

         if(config->SMEAR_HITS) {
         	// Smear the energy and timing of the hit
         	double sigma = fcal_config->FCAL_PHOT_STAT_COEF/sqrt(titer->getE());

			t += gDRandom.SampleGaussian(fcal_config->FCAL_TSIGMA);
			E *= (1.0 + gDRandom.SampleGaussian(sigma));
		 }
		 
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
