#include "FDCSmearer.h"

#include <FDC/DFDCGeometry.h>
#include <JANA/JException.h>
using namespace jana;

#include <sstream>

//-----------
// fdc_config_t  (constructor)
//-----------
fdc_config_t::fdc_config_t(JEventLoop *loop) 
{
	// default values
	FDC_TDRIFT_SIGMA      = 0.0;
 	FDC_CATHODE_SIGMA     = 0.0;
 	FDC_PED_NOISE         = 0.0;
 	FDC_THRESHOLD_FACTOR  = 0.0;
 	FDC_TIME_WINDOW       = 0.0;
 	FDC_THRESH_KEV        = 0.0;

	// load data from CCDB
	cout << "Get FDC/fdc_parms parameters from CCDB..." << endl;
    map<string, double> fdcparms;
     if(loop->GetCalib("FDC/fdc_parms", fdcparms)) {
     	jerr << "Problem loading FDC/fdc_parms from CCDB!" << endl;
     } else {
       	FDC_TDRIFT_SIGMA      = fdcparms["FDC_TDRIFT_SIGMA"];
       	FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];
       	FDC_THRESHOLD_FACTOR  = fdcparms["FDC_THRESHOLD_FACTOR"];
     	//FDC_PED_NOISE         = fdcparms["FDC_PED_NOISE"];  // ???
       	FDC_TIME_WINDOW       = fdcparms["FDC_TIME_WINDOW"];
       	//FDC_HIT_DROP_FRACTION = fdcparms["FDC_HIT_DROP_FRACTION"];   // ???
       	FDC_THRESH_KEV 		  = fdcparms["FDC_THRESH_KEV"]; 
	}

   	// Calculate ped noise level based on position resolution
   	//   FDC_PED_NOISE = -0.004594 + 0.008711*FDC_CATHODE_SIGMA +
   	//                    0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
   	FDC_PED_NOISE = -0.0938 + 0.0485*FDC_CATHODE_SIGMA;


	// load efficiency correction factors
	for(int package=1; package<=4; package++) {
		vector< vector<double> > new_strip_efficiencies;
		vector< vector<double> > new_wire_efficiencies;

		char ccdb_str[100];
		sprintf(ccdb_str, "/FDC/package%d/strip_mc_efficiency", package);
    	if(loop->GetCalib(ccdb_str, new_strip_efficiencies)) {
        	stringstream err_ss;
        	err_ss << "Error loading " << ccdb_str << " !";
        	throw JException(err_ss.str());
        }
		sprintf(ccdb_str,"/FDC/package%d/wire_mc_efficiency", package);
    	if(loop->GetCalib(ccdb_str, new_wire_efficiencies)) {
        	stringstream err_ss;
        	err_ss << "Error loading " << ccdb_str << " !";
        	throw JException(err_ss.str());
        }

		for(int chamber=0; chamber<6; chamber++) {
        	channel_efficiencies.push_back( new_strip_efficiencies[2*chamber+1] );
        	channel_efficiencies.push_back( new_wire_efficiencies[chamber] );
        	channel_efficiencies.push_back( new_strip_efficiencies[2*chamber] );
		}
	}
	
}


//-----------
// SmearEvent
//-----------
void FDCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   double t_max = config->TRIGGER_LOOKBACK_TIME + fdc_config->FDC_TIME_WINDOW;
   double threshold = fdc_config->FDC_THRESHOLD_FACTOR * fdc_config->FDC_PED_NOISE; // for sparsification

   hddm_s::FdcChamberList chambers = record->getFdcChambers();
   hddm_s::FdcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {

      // Add pedestal noise to strip charge data
      hddm_s::FdcCathodeStripList strips = iter->getFdcCathodeStrips();
      hddm_s::FdcCathodeStripList::iterator siter;
      for (siter = strips.begin(); siter != strips.end(); ++siter) {
          // If a fdcCathodeHit already exists delete it
          siter->deleteFdcCathodeHits();
          hddm_s::FdcCathodeTruthHitList thits = 
                                         siter->getFdcCathodeTruthHits();
          hddm_s::FdcCathodeTruthHitList::iterator titer;
          for (titer = thits.begin(); titer != thits.end(); ++titer) {
            // correct simulation efficiencies 
            if (config->APPLY_EFFICIENCY_CORRECTIONS
                  && !gDRandom.DecideToAcceptHit(fdc_config->GetEfficiencyCorrectionFactor(siter)))
              	continue;
          
            double q = titer->getQ();
            double t = titer->getT();
          	if(config->SMEAR_HITS) {
             	q += gDRandom.SampleGaussian(fdc_config->FDC_PED_NOISE);
             	t += gDRandom.SampleGaussian(fdc_config->FDC_TDRIFT_SIGMA)*1.0e9;
			}
            if (q > threshold && t > config->TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcCathodeHitList hits = siter->addFdcCathodeHits();
               hits().setQ(q);
               hits().setT(t);
            }
         }

         if (config->DROP_TRUTH_HITS)
            siter->deleteFdcCathodeTruthHits();
      }

      // Add drift time varation to the anode data 
      hddm_s::FdcAnodeWireList wires = iter->getFdcAnodeWires();
      hddm_s::FdcAnodeWireList::iterator witer;
      for (witer = wires.begin(); witer != wires.end(); ++witer) {
         // If a fdcAnodeHit exists already delete it
         witer->deleteFdcAnodeHits();
         hddm_s::FdcAnodeTruthHitList thits = witer->getFdcAnodeTruthHits();
         hddm_s::FdcAnodeTruthHitList::iterator titer;
         for (titer = thits.begin(); titer != thits.end(); ++titer) {
             // correct simulation efficiencies 
		     if (config->APPLY_EFFICIENCY_CORRECTIONS
             		&& !gDRandom.DecideToAcceptHit(fdc_config->GetEfficiencyCorrectionFactor(witer)))
             		continue;

            double t = titer->getT() + gDRandom.SampleGaussian(fdc_config->FDC_TDRIFT_SIGMA)*1.0e9;
            if (t > config->TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcAnodeHitList hits = witer->addFdcAnodeHits();
               hits().setT(t);
               hits().setDE(titer->getDE());
            }
         }

         if (config->DROP_TRUTH_HITS)
            witer->deleteFdcAnodeTruthHits();
      }
   }
}

