#ifndef _BCALSMEARER_H_
#define _BCALSMEARER_H_

#include "mcsmear_config.h"
#include "HDDM/hddm_s.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <queue>
#include <cmath>
using namespace std;

#include <DHistogram.h>
#include <BCAL/DBCALGeometry.h>

#include "units.h"
#include <TMath.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TDirectory.h>

#include "Smearer.h"

class bcal_config_t 
{
  public:
	bcal_config_t(JEventLoop *loop);
	
	void inline GetAttenuationParameters(int id, double &attenuation_length, double &attenuation_L1, double &attenuation_L2) {
   		vector<double> &parms = attenuation_parameters.at(id);
		attenuation_length = parms[0];
		attenuation_L1 = parms[1];
		attenuation_L2 = parms[2];
	}

	double inline GetEffectiveVelocity(int id) {
		return BCAL_C_EFFECTIVE;
   		//return effective_velocities.at(id);
	}
	
	// member variables
	double BCAL_SAMPLINGCOEFA;
	double BCAL_SAMPLINGCOEFB;
	double BCAL_TIMEDIFFCOEFA;
	double BCAL_TIMEDIFFCOEFB;
	double BCAL_TWO_HIT_RESO;
	double BCAL_mevPerPE;
	double BCAL_C_EFFECTIVE;

	double BCAL_LAYER1_SIGMA_SCALE;
	double BCAL_LAYER2_SIGMA_SCALE;
	double BCAL_LAYER3_SIGMA_SCALE;
	double BCAL_LAYER4_SIGMA_SCALE;

	int BCAL_NUM_MODULES;
	int BCAL_NUM_LAYERS;
	int BCAL_NUM_SECTORS;

	double BCAL_BASE_TIME_OFFSET;
	double BCAL_TDC_BASE_TIME_OFFSET;

	double BCAL_ADC_THRESHOLD_MEV;
	double BCAL_FADC_TIME_RESOLUTION;
	double BCAL_TDC_TIME_RESOLUTION;
	double BCAL_MEV_PER_ADC_COUNT;
	double BCAL_NS_PER_ADC_COUNT;
	double BCAL_NS_PER_TDC_COUNT;

	// BCAL flags
	bool NO_T_SMEAR;
	bool NO_DARK_PULSES;
	bool NO_SAMPLING_FLUCTUATIONS;
	bool NO_SAMPLING_FLOOR_TERM;
	bool NO_POISSON_STATISTICS;
	bool NO_FADC_SATURATION;

	vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
	// Assume constant effective velocity instead of channel-dependent one
	//vector<double> effective_velocities; // 16.75 (from calibDB BCAL/effective_velocities)

	vector< pair<double,double> > channel_efficiencies;
	
	double GetEfficiencyCorrectionFactor(int index, DBCALGeometry::End the_end) {
		if(the_end == DBCALGeometry::End::kUpstream)
			return channel_efficiencies.at(index).first;
		else 
			return channel_efficiencies.at(index).second;
	}

	double fADC_MinIntegral_Saturation, fADC_Saturation_Linear, fADC_Saturation_Quadratic;
};



// utility classes

//..........................
// bcal_index is a utility class that encapsulates the
// module, layer, sector, and end in a single object that
// can be used as a key to index an STL map. 
//..........................
class bcal_index{
   public:
      enum EndType{
         kUp,
         kDown
      };
      
      bcal_index(unsigned int module, unsigned int layer,
                 unsigned int sector, unsigned int incident_id,
                 EndType end)
       : module(module),
         layer(layer),
         sector(sector),
         incident_id(incident_id),
         end(end)
      {}
   
      unsigned int module;
      unsigned int layer;
      unsigned int sector;
      unsigned int incident_id;
      EndType end;
      
      bool operator<(const bcal_index &idx) const {
         if (module < idx.module)
            return true;
         if (module > idx.module)
            return false;
         if (layer < idx.layer)
            return true;
         if (layer > idx.layer)
            return false;
         if (sector < idx.sector)
            return true;
         if (sector > idx.sector)
            return false;
         if (incident_id < idx.incident_id)
            return true;
         if (incident_id > idx.incident_id)
            return false;
         if ((end==kUp) && (idx.end==kDown))
            return true;
         return false;
      }
};

//..........................
// CellHits is a utility class that holds information
// regarding the energy and time of depostions in a cell
//..........................
class CellHits{
   public:
      enum EndType{
         kUp,
         kDown
      };

      CellHits() : E(0.0), t(0.0)
      {}
      
      double E;
      double t;
      double Etruth;
      EndType end;
};

//..........................
// SumHits is a utility class that is used to hold info
// from the SiPMs contributing to that readout channel.
// This includes a list of CellHits objects, but also
// the total number of SiPMs that should be in the sum
// and the total up/downstream energies and times.
//..........................
class SumHits{
   public:
      SumHits()
      {}
      
      vector<CellHits *> cellhits;
      vector<double> EUP;
      vector<double> tUP;
      vector<double> EDN;
      vector<double> tDN;
};

//..........................
// fADCHit is a utility class that is used to hold info
// for a single fADC hit. 
//..........................
class fADCHit{
   public:
      fADCHit(double E, double t) : E(E), t(t)
      {}
      
      double E;
      double t;
};

//..........................
// fADCHitList is a utility class that is used to hold info
// for a set of fADCHit objects. 
//..........................
class fADCHitList{
   public:
      fADCHitList()
      {}
      
      int module;
      int sumlayer;
      int sumsector;
      
      vector<fADCHit> uphits;
      vector<fADCHit> dnhits;
};

//..........................
// TDCHitList is a utility class that is used to hold info
// for a single F1TDC hit
//..........................
class TDCHitList{
   public:
      TDCHitList()
      {}
      
      int module;
      int sumlayer;
      int sumsector;
      
      vector<double> uphits;
      vector<double> dnhits;
};

//..........................
// IncidentParticle_t is a utility class for holding the
// parameters of particles recorded as incident on the 
// BCAL (shower causing)
//..........................
class IncidentParticle_t{
   public:
      IncidentParticle_t(hddm_s::BcalTruthIncidentParticle &ipart) {
         x = ipart.getX();
         y = ipart.getY();
         z = ipart.getZ();
         px = ipart.getPx();
         py = ipart.getPy();
         pz = ipart.getPz();
         ptype = ipart.getPtype();
      }

      float x,y,z;
      float px, py, pz;
      int ptype, track;
};

// MAIN CLASS
class BCALSmearer : public Smearer
{
    public:
		BCALSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
			bcal_config = new bcal_config_t(loop);
			
			// pass configuration parameters
			bcal_config->NO_T_SMEAR = in_config->BCAL_NO_T_SMEAR;
			bcal_config->NO_DARK_PULSES = in_config->BCAL_NO_DARK_PULSES;
			bcal_config->NO_SAMPLING_FLUCTUATIONS = in_config->BCAL_NO_SAMPLING_FLUCTUATIONS;
			bcal_config->NO_SAMPLING_FLOOR_TERM = in_config->BCAL_NO_SAMPLING_FLOOR_TERM;
			bcal_config->NO_POISSON_STATISTICS = in_config->BCAL_NO_POISSON_STATISTICS;
			bcal_config->NO_FADC_SATURATION = in_config->BCAL_NO_FADC_SATURATION;
		}
		~BCALSmearer() {
			delete bcal_config;
		}

		void SmearEvent(hddm_s::HDDM *record);  // main smearing function

	protected:
		bcal_config_t *bcal_config;
		
		int inline GetCalibIndex(int module, int layer, int sector);

		void GetSiPMHits(hddm_s::HDDM *record,
        	             map<bcal_index, CellHits> &SiPMHits,
              	         vector<IncidentParticle_t> &incident_particles);
		void ApplySamplingFluctuations(map<bcal_index,
           		                       CellHits> &SiPMHits,
                   		               vector<IncidentParticle_t> &incident_particles);
		void MergeHits(map<bcal_index, CellHits> &SiPMHits, double Resolution);
		void ApplyPoissonStatistics(map<bcal_index, CellHits> &SiPMHits);
		void SortSiPMHits(map<bcal_index,
        	              CellHits> &SiPMHits,
             	          map<int, SumHits> &bcalfADC, double Resolution);
		void SimpleDarkHitsSmear(map<int, SumHits> &bcalfADC);
		void ApplyTimeSmearing(double sigma_ns, double sigma_ns_TDC, map<int, fADCHitList> &fADCHits, 
							   map<int, TDCHitList> &TDCHits);
		void FindHits(double thresh_MeV,
              		  map<int, SumHits> &bcalfADC,
              		  map<int, fADCHitList> &fADCHits,
              		  map<int, TDCHitList> &TDCHits);
		void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits,
                        		map<int, TDCHitList> &TDCHits,
                        		hddm_s::HDDM *record);
		
};



#endif  // _SMEAR_BCAL_H_
