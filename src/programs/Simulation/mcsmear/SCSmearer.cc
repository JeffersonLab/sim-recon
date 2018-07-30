#include "SCSmearer.h"

//-----------
// sc_config_t  (constructor)
//-----------
sc_config_t::sc_config_t(JEventLoop *loop) 
{
	// default values
	START_SIGMA           = 0.0; // 300ps
	START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200
	START_PADDLE_THRESHOLD  = 0.0;

    // Get the geometry
    DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
    if(!dapp){
        jerr << "Cannot get DApplication from JEventLoop!" << endl;
        return;
    }
    DGeometry* locGeometry = dapp->GetDGeometry(loop->GetJEvent().GetRunNumber());

    // Get start counter geometry
    vector<vector<DVector3> >sc_norm; 
    vector<vector<DVector3> >sc_pos;
    unsigned int MAX_SECTORS=0;
    if (locGeometry->GetStartCounterGeom(sc_pos, sc_norm))  {
		MAX_SECTORS = sc_pos.size();
        for(int sc_index=0; sc_index<sc_pos.size(); sc_index++)
            SC_START_Z.push_back( sc_pos[sc_index][0].z() );

    }

	// Load data from CCDB
    cout << "Get START_COUNTER/start_parms parameters from CCDB..." << endl;
    map<string, double> startparms;
    if(loop->GetCalib("START_COUNTER/start_parms", startparms)) {
		jerr << "Problem loading START_COUNTER/start_parms from CCDB!" << endl;
	} else {
     	START_SIGMA = startparms["START_SIGMA"] ;
     	START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];
	}
	
	cout<<"get START_COUNTER/paddle_mc_efficiency from calibDB"<<endl;
    if(loop->GetCalib("START_COUNTER/paddle_mc_efficiency", paddle_efficiencies)) {
    	jerr << "Problem loading START_COUNTER/paddle_mc_efficiency from CCDB!" << endl;
    }

    // Start counter individual paddle resolutions
    vector< vector<double> > sc_paddle_resolution_params;
    if(loop->GetCalib("START_COUNTER/time_resol_paddle_v2", sc_paddle_resolution_params))
        jout << "Error in loading START_COUNTER/time_resol_paddle_v2 !" << endl;
    else {
        if(sc_paddle_resolution_params.size() != MAX_SECTORS)
            jerr << "Start counter paddle resolutions table has wrong number of entries:" << endl
                 << "  loaded = " << sc_paddle_resolution_params.size()
                 << "  expected = " << MAX_SECTORS << endl;

        for(int i=0; i<MAX_SECTORS; i++) {
            SC_MAX_RESOLUTION.push_back( sc_paddle_resolution_params[i][0] );
            SC_BOUNDARY1.push_back( sc_paddle_resolution_params[i][1] );
            SC_BOUNDARY2.push_back( sc_paddle_resolution_params[i][2] );
            SC_SECTION1_P0.push_back( sc_paddle_resolution_params[i][3] ); 
            SC_SECTION1_P1.push_back( sc_paddle_resolution_params[i][4] );
            SC_SECTION2_P0.push_back( sc_paddle_resolution_params[i][5] ); 
            SC_SECTION2_P1.push_back( sc_paddle_resolution_params[i][6] );
            SC_SECTION3_P0.push_back( sc_paddle_resolution_params[i][7] ); 
            SC_SECTION3_P1.push_back( sc_paddle_resolution_params[i][8] );
        }
    }

    map<string,double> sc_mc_correction_factors;
    if(loop->GetCalib("START_COUNTER/mc_time_resol_corr", sc_mc_correction_factors)) {
        jout << "Error in loading START_COUNTER/mc_time_resol_corr !" << endl;
    } else {
        SC_MC_CORRECTION_P0 = sc_mc_correction_factors["P0"];
        SC_MC_CORRECTION_P1 = sc_mc_correction_factors["P1"];
    }

}


//-----------
// SmearEvent
//-----------
void SCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::StcTruthPointList truthPoints = record->getStcTruthPoints();
        
   hddm_s::StcPaddleList pads = record->getStcPaddles();
   hddm_s::StcPaddleList::iterator iter;
   for (iter = pads.begin(); iter != pads.end(); ++iter) {
      iter->deleteStcHits();
      hddm_s::StcTruthHitList thits = iter->getStcTruthHits();
      hddm_s::StcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
      	 // correct simulation efficiencies 
		 if(config->APPLY_EFFICIENCY_CORRECTIONS
		 	&& !gDRandom.DecideToAcceptHit(sc_config->GetEfficiencyCorrectionFactor(iter->getSector())))
		 	continue;

         // smear the time
         hddm_s::StcTruthPointList::iterator piter = FindMatchingTruthPoint(titer, truthPoints);
         // calculate a z-depending timing resolution
         // z is measured from the readout end of the paddles
         double z_pos = 30.;    // default value in the middle, in case we can't find a good point.  this shouldn't happen, but you never know...
         if( piter != truthPoints.end() )
             z_pos = piter->getZ() - sc_config->SC_START_Z[iter->getSector()-1];
    
         double t = titer->getT();
         double NewE = titer->getDE();
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(sc_config->GetPaddleTimeResolution(iter->getSector()-1, z_pos));
         	// smear the energy
         	double npe = titer->getDE() * 1000. *  sc_config->START_PHOTONS_PERMEV;
         	npe = npe +  gDRandom.SampleGaussian(sqrt(npe));
         	NewE = npe/sc_config->START_PHOTONS_PERMEV/1000.;
         }
         if (NewE > sc_config->START_PADDLE_THRESHOLD) {
            hddm_s::StcHitList hits = iter->addStcHits();
            hits().setT(t);
            hits().setDE(NewE);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteStcTruthHits();
   }
}

// ----------------------
// FindMatchingTruthPoint
// ----------------------
hddm_s::StcTruthPointList::iterator SCSmearer::FindMatchingTruthPoint(hddm_s::StcTruthHitList::iterator hiter, hddm_s::StcTruthPointList &truthPoints) 
{
    // Match the StcTruthHit with the most likely corresponding StcTruthPoin
    // This is needed since StcTruthHits correspond to detector hits, and so only have time and
    // energy values.   If we want to do something with a z-dependence, e.g. time resolutions,
    // we need the StcTruthPoint, which has a location in detector coordinates.
    // The only thing they have in common in the energy deposited in the scintillator paddles
    // since the StcTruthHit has a propagation time correction applied, so we use that
    // to disambiguate multiple hits in the same paddle
    hddm_s::StcTruthPointList::iterator piter;
    hddm_s::StcTruthPointList::iterator best_piter = truthPoints.end();
    double best_match_deltaE = 100.;
    for( piter = truthPoints.begin(); piter != truthPoints.end(); piter++) {
        if( hiter->getSector() == piter->getSector() ) {
            double deltaE = fabs(hiter->getDE() - piter->getDEdx());
            if(deltaE < best_match_deltaE) {
                best_piter = piter;
                best_match_deltaE = deltaE;
            }
        }
    }

    return best_piter;
}
