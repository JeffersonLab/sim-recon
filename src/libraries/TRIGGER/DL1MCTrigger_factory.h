#ifndef _DL1MCTrigger_factory_
#define _DL1MCTrigger_factory_

#include <JANA/JFactory.h>

#include "DL1MCTrigger.h"


#include "TTAB/DTranslationTable.h"

#include "FCAL/DFCALGeometry.h"
#include <FCAL/DFCALHit.h>

#include <BCAL/DBCALHit.h>

#include <DRandom2.h>

#include <TH1.h>
#include <TH2.h>

typedef  vector< vector<double> >  fcal_constants_t;

class DL1MCTrigger_factory:public jana::JFactory<DL1MCTrigger>{
	public:
		DL1MCTrigger_factory(){};
		~DL1MCTrigger_factory(){};

		static const int  sample        =  125;
		static const int  time_stamp    =  4.;
		static const int  max_adc_bins  =  4096;        /* number of FADC bins */

		static const int FCAL_MAX_CHANNELS = 2800;

		static const int TRIG_BASELINE  = 100;

		fcal_constants_t fcal_gains;
                fcal_constants_t fcal_pedestals;
 
		int fcal_ssp[sample];
		int fcal_gtp[sample];

		int bcal_ssp[sample];
		int bcal_gtp[sample];

		typedef struct {
		  int row;
		  int column;
		  
		  double energy;
		  double time;
		  
		  double adc_en[sample];
		  int adc_amp[sample];
		  
		  int merged;
		} fcal_signal;
		

		typedef struct {
		  int module;
		  int layer;
		  int sector;
		  int end;

		  double time;      /* Pulse time in ns */
		  double energy;    /* Pulse energy in MeV */

		  double adc_en[sample];
		  int adc_amp[sample];
		  
		  int merged;
		}  bcal_signal;
				
		typedef struct {
		  
		  int  fcal;
		  int  bcal;
		  int  en_thr;  
    
		  int  fcal_min;
		  int  fcal_max;
		  int  bcal_min;
		  int  bcal_max;
		  
		  int  st_nhit;
		  int  tagh_nhit;
		  int  tof_nhit;
		  int  ps_nhit;
		  
		  unsigned int st_pattern;
		  unsigned int tagh_pattern;
		  unsigned int tof_pattern;
		  unsigned int ps_pattern;		  
		} gtp_par;
		
		
		typedef struct {
		  
		  int  status;
		  
		  int  bit;
		  int  type;   // Bit pattern: bit 1:  FCAL
		  //              bit 2:  BCAL
		  //              bit 3:  ST
		  //              bit 4:  PS
		  //              bit 5:  TAGH
		  //              bit 6:  TOF
		  		  
		  gtp_par gtp;
		  
		  float  fcal_en;
		  float  bcal_en;
		  float  en_thr;
		  
		  float  fcal_min_en;
		  float  fcal_max_en; 
		  float  bcal_min_en;
		  float  bcal_max_en;		  
		} trigger_conf;
				
		typedef struct{
		  int roc;
		  int slot;
		  int ch;

		  int col;
		  int row;
		} fcal_mod;

		typedef struct{
		  int roc;
		  int slot;
		  int ch;

		  int module;
		  int layer; 
		  int sector;
		  int end;
		} bcal_mod;

		vector<trigger_conf> triggers_enabled;

		vector<fcal_mod> fcal_trig_mask;
		vector<bcal_mod> bcal_trig_mask;

		vector<fcal_signal> fcal_signal_hits;
		vector<bcal_signal> bcal_signal_hits;

		vector<fcal_signal> fcal_merged_hits;
		vector<bcal_signal> bcal_merged_hits;

		int    BYPASS;
		float  FCAL_ADC_PER_MEV;
		int    FCAL_CELL_THR;
		int    FCAL_EN_SC;
		int    FCAL_NSA;
		int    FCAL_NSB;
		int    FCAL_WINDOW;

		float  BCAL_ADC_PER_MEV;
		int    BCAL_CELL_THR;
		int    BCAL_EN_SC;
		int    BCAL_NSA;
		int    BCAL_NSB;
		int    BCAL_WINDOW;

		int    FCAL_BCAL_EN;
		
		float  ST_ADC_PER_MEV;
		int    ST_CELL_THR;
		int    ST_NSA;
		int    ST_NSB;
		int    ST_WINDOW;
		int    ST_NHIT;
		
		int BCAL_OFFSET;
		
		int SIMU_BASELINE;
		int SIMU_GAIN;


		double time_shift;
		double time_min;
		double time_max;
		
		int simu_baseline_fcal;
		int simu_baseline_bcal;
		double pedestal_sigma;

		int simu_gain_fcal;
		int simu_gain_bcal;


		int Read_RCDB(int32_t runnumber);		
		int SignalPulse(double en, double time, double amp_array[sample], int type);

		void AddBaseline(double adc_amp[sample], double pedestal, DRandom2 &gDRandom);

		void Digitize(double adc_amp[sample], int adc_count[sample]);

		template <typename T>  int FADC_SSP(vector<T> merged_hits, 
						    int detector);
		int GTP(int detector);
		int FindTriggers(DL1MCTrigger *trigger);
		void PrintTriggers();				

		float  BCAL_ADC_PER_MEV_CORRECT;
		
 private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		
		void LoadFCALConst( fcal_constants_t &table, 
				    const vector<double> &fcal_const_ch, 
				    const DFCALGeometry  &fcalGeom);	


		void GetSeeds(JEventLoop *loop,  uint64_t eventnumber, UInt_t &seed1, UInt_t &seed2, UInt_t &seed3);

		TH1F *hfcal_gains;
		TH2F *hfcal_gains2;
		TH1F *hfcal_ped;

		int debug;
		
};

#endif // _DL1MCTrigger_factory_

