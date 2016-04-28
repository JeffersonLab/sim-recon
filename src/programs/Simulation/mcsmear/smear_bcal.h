#ifndef _SMEAR_BCAL_H_
#define _SMEAR_BCAL_H_

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

namespace bcal_smearing {

class BCALSmearer
{
    public:
		BCALSmearer(mcsmear_config_t *in_mcsmear_config, bcal_config_t *in_bcal_config) {   // constructor
			mcsmear_config = in_mcsmear_config;
			bcal_config = in_bcal_config;
		}
		~BCALSmearer();  // destructor 

		Smear(hddm_s::HDDM *record);  // main smearing function

	protected:
		mcsmear_config_t *mcsmear_config;
		bcal_config_t *bcal_config;
		
		int inline GetCalibIndex(int module, int layer, int sector);
		void inline GetAttenuationParameters(int id, double &attenuation_length, double &attenuation_L1, double &attenuation_L2);
		double inline GetEffectiveVelocity(int id);

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

};

#endif  // _SMEAR_BCAL_H_