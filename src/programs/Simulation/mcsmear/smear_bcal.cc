// $Id: smear.cc 7650 2011-03-29 22:52:30Z shepherd $
//
// Created June 22, 2005  David Lawrence
//
// Major revision March 6, 2012 David Lawrence

// Set the following to 1 to use the new timing spectrum scheme and 0 to use old scheme

#if 0

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
#include "HDDM/hddm_s.h"
#include <TMath.h>
#include <TSpline.h>
#include <TDirectory.h>

#include "DRandom2.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

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
		
		bcal_index(unsigned int module, unsigned int layer, unsigned int sector, EndType end):
			module(module),layer(layer),sector(sector),end(end){}
	
		unsigned int module;
		unsigned int layer;
		unsigned int sector;
		EndType end;
		
		// This is needed in order to use this class as the key in an STL map
		bool operator<(const bcal_index &idx) const{
			if(module<idx.module)return true;
			if(module>idx.module)return false;
			if(layer<idx.layer)return true;
			if(layer>idx.layer)return false;
			if(sector<idx.sector)return true;
			if(sector>idx.sector)return false;
			if((end==kUp) && (idx.end==kDown))return true;
			return false;
		}
};

//..........................
// CellSpectra is a utility class that holds pointers to
// histograms for each end of the cell and the total
// unsmeared, energy deposited in the cell
//..........................
class CellSpectra{
	public:
		CellSpectra():hup(NULL),hdn(NULL),Etruth(0.0){}
	
		DHistogram *hup;
		DHistogram *hdn;
		
		double Etruth;
};

//..........................
// SumSpectra is a utility class that is used to hold info
// from the SiPMs contributing to that readout channel.
// This includes a list of CellSpectra objects, but also
// the total number of SiPMs that should be in the sum
// so that dark hits for non-hit SiPMs can be added.
//..........................
class SumSpectra{
	public:
		SumSpectra():NSiPMs(0),hup(NULL),hdn(NULL){}
		
		vector<CellSpectra *> cellspectra;
		int NSiPMs;
	
		DHistogram *hup;
		DHistogram *hdn;
		
};

//..........................
// fADCHit is a utility class that is used to hold info
// for a single fADC hit. 
//..........................
class fADCHit{
	public:
		fADCHit(double E, double t):E(E),t(t){}
		
		double E;
		double t;
};

//..........................
// fADCHitList is a utility class that is used to hold info
// for a set of fADCHit objects. 
//..........................
class fADCHitList{
	public:
		fADCHitList(){}
		
		int module;
		int fADClayer;
		int fADCsector;
		
		vector<fADCHit> uphits;
		vector<fADCHit> dnhits;
		
};


// Defined in this file
void bcalInit(void);
void GetSiPMSpectra(s_HDDM_t *hddm_s, map<bcal_index, CellSpectra> &SiPMspectra);
void ApplySamplingFluctuations(map<bcal_index, CellSpectra> &SiPMspectra);
void ApplyPoissonStatistics(map<bcal_index, CellSpectra> &SiPMspectra);
void ApplySiPMTimeJitter(map<bcal_index, CellSpectra> &SiPMspectra);
void AddDarkHits(map<bcal_index, CellSpectra> &SiPMspectra);
double AddDarkHitsToOne(DHistogram *h, int NSiPMs=1);
void SortSiPMSpectra(map<bcal_index, CellSpectra> &SiPMspectra, map<int, SumSpectra> &bcalfADC);
void AddDarkHitsForNonHitSiPMs(map<int, SumSpectra> &bcalfADC);
void ApplyElectronicPulseShape(map<int, SumSpectra> &bcalfADC);
void ApplyElectronicPulseShapeOneHisto(DHistogram *h);
void FindHits(double thresh_mV, map<int, SumSpectra> &bcalfADC, map<int, fADCHitList> &fADCHits);
void FindHitsOneHisto(double thresh_mV, DHistogram *h, vector<fADCHit> &hits);
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits, s_HDDM_t *hddm_s);
DHistogram* GetHistoFromPool(void);
void ReturnHistoToPool(DHistogram*);
TSpline* MakeTSpline(void);
void SaveDebugHistos(map<bcal_index, CellSpectra> &SiPMspectra, const char *suffix="", const char *yunits="Amplitude");
void SaveDebugHistos(map<int, SumSpectra> &bcalfADC, const char *suffix="");

// The following are set in bcalInit below
bool BCAL_INITIALIZED = false;
double BCAL_mevPerPE=0.0;
double *BCAL_PULSE_SHAPE_MATRIX=NULL;

// Mutexes and histogram pool
pthread_mutex_t bcal_init_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t histo_pool_mutex = PTHREAD_MUTEX_INITIALIZER;
queue<DHistogram*> histo_pool;

// Mutex used to control accessing the ROOT global memory
extern pthread_mutex_t root_mutex;

// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn OFF the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_THRESHOLD_CUT;
extern bool BCAL_DEBUG_HISTS;


// setup response parameters
extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
extern double BCAL_XTALK_FRACT;                 // 0.157  (from calibDB BCAL/bcal_parms)
extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms) 
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)

// The following are not currently in use
//extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
//extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
//extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
//extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
//extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
//extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)


double BCAL_MIN_ENERGY_MEV = 30.0;
double BCAL_FADC_INTEGRATION_WINDOW_PRE  =  20.0; // time to integrate signal before threshold crossing in ns
double BCAL_FADC_INTEGRATION_WINDOW_POST = 180.0; // time to integrate signal after threshold crossing in ns

//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
	/// The data from HDGEANT contains attenuated, energy-weighted, timing spectra.
	/// The job of mcsmear is to use that data to create the bcalfADCUpHit and
	/// bcalfADCDownHit structures. These are made by summing signals from multiple
	/// SiPMs and are what are used as the entry points for the reconstruction.
	///
	/// The time spectra must be smeared due to sampling fluctuations, which are
	/// parameterized based on the total energy of the shower. The total
	/// energy is kept in the Etruth field of the bcalSiPMSpectrum structures.
	///
	/// In addition to the sampling fluctuations, Poisson statistics, timing jitter,
	/// dark pulses, and timewalk effects are all applied. The timewalk effect is obtained
	/// by convoluting the energy spectra with an electronic pulse shape and identifying
	/// the threshold crossing points.
	///
	
	
	// n.b. This code is slightly more complex than it might otherwise be because
	// it uses sparsified lists as opposed to full data structures for every SiPM
	// and every summed cell. It is done this way for two reasons:
	//
	// 1.) Sparsified lists are quicker to loop through and make the code faster.
	//
	// 2.) Sparsified lists are easier to keep on the stack rather than the heap
	//     avoiding expensive, large memory allocations every event. Alternatively
	//     one could keep the large structures in global variables but that would 
	//     not allow for efficient multi-threading.


	// Initialize BCAL globals on first call
	if(!BCAL_INITIALIZED)bcalInit();
	
	// First, we extract the time spectra for hit cells
	map<bcal_index, CellSpectra> SiPMspectra;
	GetSiPMSpectra(hddm_s, SiPMspectra);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(SiPMspectra, "_Raw", "Amplitude (attenuated-MeV)");

	// Sampling fluctuations
	ApplySamplingFluctuations(SiPMspectra);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(SiPMspectra, "_Sampling", "Amplitude (attenuated-MeV)");
	
	// Poisson Statistics
	ApplyPoissonStatistics(SiPMspectra);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(SiPMspectra, "_Poisson", "Amplitude (attenuated-MeV)");
	
	// Time Jitter
	ApplySiPMTimeJitter(SiPMspectra);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(SiPMspectra, "_TimeJitter", "Amplitude (attenuated-MeV)");
	
	// Add Dark Hits (for hit cells only at this point)
	AddDarkHits(SiPMspectra);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(SiPMspectra, "_DarkHits", "Amplitude (attenuated-MeV)");
	
	// Place all hit cell spectra into list indexed by fADC ID
	map<int, SumSpectra> bcalfADC;
	SortSiPMSpectra(SiPMspectra, bcalfADC);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(bcalfADC, "_SummedCell");

	// Add Dark Hits (for cells without energy deposited)
	// (n.b. this may add elements to bcalfADC for summed cells with no signal)
	AddDarkHitsForNonHitSiPMs(bcalfADC);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(bcalfADC, "_DarkHits_SummedCell");
	
	// Electronic Pulse shape
	ApplyElectronicPulseShape(bcalfADC);
	if(BCAL_DEBUG_HISTS)SaveDebugHistos(bcalfADC, "_Electronic");
	
	// Apply threshold to find hits in summed spectra
	double thresh_mV = 44.7;
	map<int, fADCHitList> fADCHits;
	FindHits(thresh_mV, bcalfADC, fADCHits);
	
	// Copy hits into HDDM tree
	CopyBCALHitsToHDDM(fADCHits, hddm_s);
	
	// Return histogram objects to the pool for recycling
	map<bcal_index, CellSpectra>::iterator iter = SiPMspectra.begin();
	for(; iter != SiPMspectra.end(); iter++){
		CellSpectra &cellspectra = iter->second;
		if(cellspectra.hup)ReturnHistoToPool(cellspectra.hup);
		if(cellspectra.hdn)ReturnHistoToPool(cellspectra.hdn);
	}
	SiPMspectra.clear();

	map<int, SumSpectra>::iterator biter = bcalfADC.begin();
	for(; biter != bcalfADC.end(); biter++){
		SumSpectra &sumspectra = biter->second;
		if(sumspectra.hup)ReturnHistoToPool(sumspectra.hup);
		if(sumspectra.hdn)ReturnHistoToPool(sumspectra.hdn);
	}
	bcalfADC.clear();
}

//-----------
// bcalInit
//-----------
void bcalInit(void)
{
	/// Pre-calculate some values that will be used every event

	// Enclose this entire routine in a mutex lock and double check 
	// the intialization flag to guarantee only one thread actually
	// runs this and that all other threads wait for it to finish
	// before continuing.
	pthread_mutex_lock(&bcal_init_mutex);

	// If we did not get here first, then initialization has already
	// been done. Return immediately.
	if(BCAL_INITIALIZED){
		pthread_mutex_unlock(&bcal_init_mutex);
		return;
	}
	
	//================ Conversion factor: MeV per PE =======================
	BCAL_mevPerPE = 1.0/( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );
	//======================================================================


	//================== Electronic Pulse Shape ============================
	// The following will initialize a table that will be used
	// for convoluting the electronic pulse_shape function. The
	// conversion is of the form of a matrix multiplication with
	// the input and output histograms being vectors and this table
	// being the transformation matrix. It is therefore a NBINSxNBINS
	// matrix with NBINS being the number of bins in the histogram
	// being convoluted.
	//
	// For simplicity, we use the same histogram definition for the 
	// input energy-weighted time histogram and the output electronic
	// pulse shape histogram. It's possible some speed optimization
	// could be done by making these different.
	
	cout<<"Pre-evaluating BCAL pulse_shape function ..."<<endl;

	// Get a histogram from the pool just so we can get it's dimensions.
	DHistogram *h = GetHistoFromPool();

	int NBINS = h->GetNbins();
	BCAL_PULSE_SHAPE_MATRIX = new double[NBINS*NBINS];
	
	// Get electronic pulse shape
	TSpline *pulse_shape = MakeTSpline();

	// Conversion factor between attenuated MeV and mV
	// (see slides from BCAL segmentation meeting 7/22/2011
	// https://halldweb1.jlab.org/wiki/images/3/3f/20110722_bcal.pdf )
	double QCD_counts_per_PE = 61.67;
	double LSB_pC_CAEN_V792 = 0.100;
	double gain_factor_for_test = 20.0;
	double SiPM_pulse_integral_pC = 1200.0;
	double SiPM_pulse_peak_mV = 2293.0;
	double TDC_gain_factor = 10.0;
	
	double mV_per_MeV = 1.0; // just a factor to start with so we can multiply/divide easily below
	mV_per_MeV *= QCD_counts_per_PE;		// QCD counts/PE
	mV_per_MeV *= LSB_pC_CAEN_V792;			// pC/PE
	mV_per_MeV /= BCAL_mevPerPE;			// pC/MeV
	mV_per_MeV /= gain_factor_for_test;		// convert to actual gain of pre-amp
	mV_per_MeV *= SiPM_pulse_peak_mV/SiPM_pulse_integral_pC; // mV/MeV
	mV_per_MeV *= TDC_gain_factor; 			// account for factor of 10 gain in TDC leg of pre-amp

	// Normalize pulse shape between t=15ns and t=100ns, excluding the
	// pre and post bumps. We scale by the amplitude to give the pulse
	// shape a height of 1 so that when multiplied by mV_per_MeV, the
	// amplitude will again be in mV.
	//double norm = pulse_shape->GetMinimum(15.0, 100.0);
	double norm = -1020.80; // mV not trivial to find minimum with TSpline3

	// The pulse shape used here represents an amplification factor
	// of 10 for the TDC leg. Elton claims this is way too high
	// and will be lowered to probably no more than 5. 
	double preamp_gain_tdc = 5.0;

	// Scale "norm" by ratio of gain used to make pulse and actual gain we expect.
	// (The data used to get the TSpline was obtained using a pre-amp with a x10 gain)
	norm *= 10.0/preamp_gain_tdc;

	// time shift pulse shape just to bring it in frame better for zoomed in picture
	double toffset = 6.0;
	
	// Loop over bins of input histo (attenuated energy)
	for(int ibin=1; ibin<=NBINS; ibin++){

		double t0 = h->GetBinCenter(ibin);

		// Loop over bins of output histo (electronic pulse)
		for(int jbin=1; jbin<=NBINS; jbin++){
			double t1 = h->GetBinCenter(jbin);
			double t = t1-t0+toffset;
			int idx = (jbin-1) + (ibin-1)*NBINS;

			// The spline data starts at -7ns and ends at 180ns. Beyond those 
			// limits it tends to diverge. Impose a cut here to keep the signal
			// at zero when evaluating the spline beyond those limits.
			if(t<-7.0 || t>180.0){
				BCAL_PULSE_SHAPE_MATRIX[idx] = 0.0;
			}else{
				BCAL_PULSE_SHAPE_MATRIX[idx] = (mV_per_MeV/norm)*pulse_shape->Eval(t-t0+toffset);
			}
		}
	}
	
	// Clean up
	delete pulse_shape; // Get rid of TSpline since we don't need it anymore
	ReturnHistoToPool(h);
	
	cout<<"Done."<<endl;
	//======================================================================
	
	BCAL_INITIALIZED = true;
	pthread_mutex_unlock(&bcal_init_mutex);

}

//-----------
// GetSiPMSpectra
//-----------
void GetSiPMSpectra(s_HDDM_t *hddm_s, map<bcal_index, CellSpectra> &SiPMspectra)
{
	/// Loop through input HDDM data and extract the timing spectrum info into
	/// CellSpectra objects.

	// Loop over PhysicsEvents
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	for(unsigned int iphysics_event=0; iphysics_event<PE->mult; iphysics_event++){
	
		// Get pointer to HDDM hits structure
		s_HitView_t *hits = PE->in[iphysics_event].hitView;
		
		// Make sure HDDM stuctures exist.
		// In the case of no real BCAL hits, we may still want to emit
		// dark hit only events. In this case, we must create the BCAL 
		// tree here.
		if(hits==HDDM_NULL)hits = PE->in[iphysics_event].hitView = make_s_HitView();
		if(hits->barrelEMcal == HDDM_NULL)hits->barrelEMcal = make_s_BarrelEMcal();
		if(hits->barrelEMcal->bcalCells == HDDM_NULL)hits->barrelEMcal->bcalCells = make_s_BcalCells(0);

		// Loop over GEANT hits in BCAL
		s_BcalSiPMSpectrums_t *bcalSiPMSpectrums = hits->barrelEMcal->bcalSiPMSpectrums;
		for(unsigned int j=0; j<bcalSiPMSpectrums->mult; j++){
			s_BcalSiPMSpectrum_t *bcalSiPMSpectrum = &bcalSiPMSpectrums->in[j];
			bcal_index idx(bcalSiPMSpectrum->module, bcalSiPMSpectrum->layer, bcalSiPMSpectrum->sector, bcalSiPMSpectrum->end==0 ? bcal_index::kUp:bcal_index::kDown);
	
			// Get reference to existing CellSpectra, or create one if it doesn't exist
			CellSpectra &cellspectra = SiPMspectra[idx];
			
			DHistogram *h=NULL;
			if(idx.end==bcal_index::kUp){
				if(!cellspectra.hup)cellspectra.hup = GetHistoFromPool();
				h = cellspectra.hup;
			}else{
				if(!cellspectra.hdn)cellspectra.hdn = GetHistoFromPool();
				h = cellspectra.hdn;
			}
			
			// The energy truth info for the cell is stored with the spectra for both
			// ends. (This is for when hdgeant cuts the spectrum for one end away
			// due to threshold). They should be the same if both spectra are here so
			// we overwrite the one data member.
			cellspectra.Etruth = bcalSiPMSpectrum->Etruth;
			
			double t = bcalSiPMSpectrum->tstart;
			double bin_width = bcalSiPMSpectrum->bin_width;
			stringstream ss(bcalSiPMSpectrum->vals);
			
			// Extract values and use them to fill histo
			int i = 0;
			double dE;
			while(ss>>dE){

				h->Fill(t, dE);
				t += bin_width;
				
				if(++i > 100000){
					_DBG_<<"More than 100k iterations parsing BCAL time spectrum string! Something is wrong!"<<endl;
					exit(-1);
				}
			}
		}
	}
}

//-----------
// ApplySamplingFluctuations
//-----------
void ApplySamplingFluctuations(map<bcal_index, CellSpectra> &SiPMspectra)
{
	/// Loop over the CellSpectra objects and apply sampling fluctuations.
	///
	/// Here we apply a statistical error due to the sampling
	/// fluctuations. The total energy (Etruth) is integrated by hdgeant.
	/// We calculate a sigma based on the deposited energy only. In
	/// reality, the sampling fluctuations are also a function of the
	/// angle of the shower particles w.r.t. the fibers. We do not include
	/// any angular dependence at this time. To do so will require more
	/// information be passed from hdgeant.
	///
	/// The error is applied by finding the ratio of the smeared
	/// cell energy to unsmeared cell energy and scaling both Eup and Edn
	/// by it.

	map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		CellSpectra &cellspectra = iter->second;
		
		// Find fractional sampling sigma based on deposited energy (whole colorimeter, not just fibers)
		double Etruth = cellspectra.Etruth;
		double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( Etruth ) + BCAL_SAMPLINGCOEFB;

		// Convert sigma into GeV
		sigmaSamp *= Etruth;

		// Randomly sample the fluctuation
		double Esmeared = gDRandom.Gaus(Etruth,sigmaSamp);

		// Calculate ratio of smeared to unsmeared
		double ratio = Esmeared/Etruth;

		// Scale attenuated energy histos for each end
		if(cellspectra.hup)cellspectra.hup->Scale(ratio);
		if(cellspectra.hdn)cellspectra.hdn->Scale(ratio);
	}
}

//-----------
// ApplyPoissonStatistics
//-----------
void ApplyPoissonStatistics(map<bcal_index, CellSpectra> &SiPMspectra)
{
	/// Loop over the CellSpectra objects and apply Poisson Statistics.
	///
	/// Because the response of the SiPM is quantized in units of photo-electrons
	/// Poisson counting statistics should be applied. This will affect the
	/// smaller energy depositions more than the larger ones.
	///
	/// We do this by converting the integral of the attenuated energy into
	/// photo-electrons and then sampling from a Poisson distribution with that
	/// mean. The ratio of the quantized, sampled value to the unquantized
	// integral (in PE) is used to scale the spectrum.

	map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		CellSpectra &cellspectra = iter->second;
		
		// Integrate hup and hdn
		double Iup = (cellspectra.hup ? cellspectra.hup->Integral():0.0); // in attenuated MeV
		double Idn = (cellspectra.hdn ? cellspectra.hdn->Integral():0.0); // in attenuated MeV

		// Upstream end
		if(Iup>0.0){
			// Convert to number of PE
			double mean_pe_up = Iup/BCAL_mevPerPE;
			
			int Npe_up = gDRandom.Poisson(mean_pe_up);
			double ratio_up = (double)Npe_up/mean_pe_up;

			cellspectra.hup->Scale(ratio_up); // only get here if cellspectra.hup is non-NULL
		}

		// Downstream end
		if(Idn>0.0){
			// Convert to number of PE
			double mean_pe_dn = Idn/BCAL_mevPerPE;
			
			int Npe_dn = gDRandom.Poisson(mean_pe_dn);
			double ratio_dn = (double)Npe_dn/mean_pe_dn;
			
			cellspectra.hdn->Scale(ratio_dn); // only get here if cellspectra.hdn is non-NULL
		}
	}
}

//-----------
// ApplySiPMTimeJitter
//-----------
void ApplySiPMTimeJitter(map<bcal_index, CellSpectra> &SiPMspectra)
{
	/// Loop over the CellSpectra objects and apply SiPM time jitter.
	///
	/// The photo detection device has an intrinsic timing resolution
	/// caused by jitter in the time between the input signal (light)
	/// and the generated electronic pulse. This is quoted in the 
	/// spec sheet as "Time Resolution (FWHM)" for a single photon.
	/// (see page 2 of the following link:
	/// http://sales.hamamatsu.com/assets/pdf/parts_S/s10362-33series_kapd1023e05.pdf)
	///
	/// n.b. XP2020 PMT's have a very similar time jitter as the SiPMs
	///
	/// To include this effect, each bin is converted into a number
	/// of photo-electrons and each of them is randomly displaced in time
	/// with the specified width. Because the photo-statistics is
	/// applied on a global scale, we do not re-quantize the bin contents
	/// into an integer number of photo-electrons. Rather, we split it
	/// into a number of pieces that is equal to the approximate number
	/// of photo-electrons coming from this bin and shift those pieces.
	/// Bins with trace amounts of energy are still kept and shifted
	/// as well to keep the energy integral constant.
	
	double fwhm_jitter = 0.600; // ns
	double sigma_jitter = fwhm_jitter/2.35; // ns
	
	// Borrow a histogram from the pool
	DHistogram *tmp = GetHistoFromPool();
	int Nbins = tmp->GetNbins();

	map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		CellSpectra &cellspectra = iter->second;
		
		// Upstream
		if(cellspectra.hup){
			tmp->Reset();
			for(int ibin=1; ibin<=Nbins; ibin++){
				double E = cellspectra.hup->GetBinContent(ibin);
				if(E==0.0)continue;
				double tbin = cellspectra.hup->GetBinCenter(ibin);

				int Npe = 1 + (int)floor(E/BCAL_mevPerPE);
				double dE = E/(double)Npe;
				
				for(int i=0; i<Npe; i++){
					double t = gDRandom.Gaus(tbin, sigma_jitter);
					tmp->Fill(t, dE);
				}
			}
			*cellspectra.hup = *tmp; // overwrite histo with tmp
		}

		// Downstream
		if(cellspectra.hdn){
			tmp->Reset();
			for(int ibin=1; ibin<=Nbins; ibin++){
				double E = cellspectra.hdn->GetBinContent(ibin);
				if(E==0.0)continue;
				double tbin = cellspectra.hdn->GetBinCenter(ibin);

				int Npe = 1 + (int)floor(E/BCAL_mevPerPE);
				double dE = E/(double)Npe;
				
				for(int i=0; i<Npe; i++){
					double t = gDRandom.Gaus(tbin, sigma_jitter);
					tmp->Fill(t, dE);
				}
			}
			*cellspectra.hdn = *tmp; // overwrite histo with tmp
		}
	}
   
   // Return histogram to pool
   ReturnHistoToPool(tmp);
}

//-----------
// AddDarkHits
//-----------
void AddDarkHits(map<bcal_index, CellSpectra> &SiPMspectra)
{
	/// Loop over the CellSpectra objects and add Dark Hits.
	///
	/// The SiPMs will randomly fire pixels at some rate, even when no
	/// light is present. These "dark hits" will add to the signal
	/// spectrum. This can affect the timing if the hits show up close
	/// to the leading edge of the real signal. They can also add up
	/// enough to occasionally cross threshold, even if no actual signal is 
	/// present. In this routine, dark hits are added to the SiPMs
	/// that already have signal. SiPM spectra with no signal will have
	/// dark hits added later in AddDarkHitsForNonHitSiPMs().
	///
	
	if(NO_DARK_PULSES)return;

	map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		CellSpectra &cellspectra = iter->second;
		
		if(cellspectra.hup)AddDarkHitsToOne(cellspectra.hup);
		if(cellspectra.hdn)AddDarkHitsToOne(cellspectra.hdn);
	}
}

//-----------
// AddDarkHitsToOne
//-----------
double AddDarkHitsToOne(DHistogram *h, int NSiPMs)
{
	/// Dark hits are added by using 1/rate as the mean time between
	/// hits. It is assumed that the timing distribution bewteen hits
	/// is exponential. The time differences are randomly sampled until
	/// a hit is found past the end of the histogram limit.
	///
	/// The second parameter is used to scale the rate up for the case 
	/// when multiple SiPMs are summed.
	///
	/// The return value is the total energy equivalent in MeV added
	/// due to dark hits. This includes primaries plus cross-talk.
	///
	/// The parameters for the dark noise should come from the CCDB. For 
	/// a reference, one can look at GlueX-doc-1754. Slide 16 has the total
	/// dark rate (17.6MHz), and slide 17 has the cross-talk (15.7%)

	int Nbins = h->GetNbins();
	double low_edge = h->GetBinLowEdge(1);
	double high_edge = h->GetBinLowEdge(Nbins+1);
	double t = low_edge;

	double tau_ns = 1.0/((double)NSiPMs*BCAL_DARKRATE_GHZ);
	
	double dE_tot = 0;
	while(t < high_edge){ // actually break inside while loop so this is superfluous

		// Loop until we get a finite time difference
		// (both s=0 and s=1 are poorly behaved).
		double delta_t = 0.0;
		do{
			double s = gDRandom.Rndm();
			delta_t = tau_ns*log(1.0/(1.0-s)); // derived by setting s equal to integral fraction of exp(t/tau)
		}while(!finite(delta_t));
		
		t += delta_t;
		if(t > high_edge) break;
		
		double dE = BCAL_mevPerPE;  // single photo electron
		if(gDRandom.Rndm() <= BCAL_XTALK_FRACT){
			dE += BCAL_mevPerPE; // single cross-talk pixel
			if(gDRandom.Rndm() <= BCAL_XTALK_FRACT){
				dE += BCAL_mevPerPE; // second cross-talk pixel
			}
		}
		h->Fill(t, dE);
		dE_tot += dE;
	}

	return dE_tot;
}

//-----------
// SortSiPMSpectra
//-----------
void SortSiPMSpectra(map<bcal_index, CellSpectra> &SiPMspectra, map<int, SumSpectra> &bcalfADC)
{
	/// Loop over the CellSpectra objects and copy pointers to them into SumSpectra objects.
	///
	/// For the BCAL, multiple SiPMs are summed together. This routine gathers individual
	/// SiPM spectra into single SumSpectra objects. Each SumSpectra represents a summed
	/// cell that is readout by an fADC channel. Since not every cell has signal in it, each
	/// SumSpectra object may not have as many input spectra as SiPMs that will actually be
	/// contributing. This is only an issue for dark hits in SiPMs with no signal. Those 
	/// non-signal SiPM dark hits are handled later in AddDarkHitsForNonHitSiPMs.
	///
	/// To save space and, more importantly, time, the first histogram is promoted to be the
	/// histo for the SumSpectra and the pointer in the CellSpectra object set to NULL. (This 
	/// prevents it from being returned to the pool twice at the end of the event.)
	
	// Loop over SiPMspectra and copy a pointer to it to the correct SumSpectra
	// element in the bcalfADC container. The bcalfADC container is an STL map which is
	// convenient since it creates a new SumSpectra object for us if it doesn't exist,
	// but otherwise, returns a reference to the existing object.
	
	map<bcal_index, CellSpectra>::iterator iter = SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		
		// Get reference to SumSpectra object
		const bcal_index &idx = iter->first;
		int fADCId = DBCALGeometry::fADCId( idx.module, idx.layer, idx.sector);
		SumSpectra &sumspectra = bcalfADC[fADCId];
		
		// Add CellSpectra object to list in SumSpectra
		CellSpectra &cellspectra = iter->second;
		sumspectra.cellspectra.push_back(&cellspectra);
		
		// If this is the first cell added to the SumSpectra, re-assign its
		// histo to the SumSpectra. Otherwise, just add to it.

		// Upstream
		if(sumspectra.hup==NULL){
			sumspectra.hup = cellspectra.hup;
			cellspectra.hup = NULL;
		}else{
			if(cellspectra.hup)sumspectra.hup->Add(cellspectra.hup);
		}
		
		// Downstream
		if(sumspectra.hdn==NULL){
			sumspectra.hdn = cellspectra.hdn;
			cellspectra.hdn = NULL;
		}else{
			if(cellspectra.hdn)sumspectra.hdn->Add(cellspectra.hdn);
		}
	}
}

//-----------
// AddDarkHitsForNonHitSiPMs
//-----------
void AddDarkHitsForNonHitSiPMs(map<int, SumSpectra> &bcalfADC)
{
	/// Loop over the SumSpectra objects and add Dark Hits for channels
	/// that do not have any signal.
	///
	/// Dark hits have already been added for cells that have signal in them
	/// (see AddDarkHits(...) ). Here, we add them in for the rest of the 
	/// channels. There are two categories of channels we need to add Dark Hits
	/// for:
	///
	/// 1. Summed cells that already have at least one SiPM with signal, but
	///    less than the number of summed SiPMs have signal
	///
	/// 2. Summed cells with no SiPMs having signal. 
	///
	/// Case 1 will add to the existing bcalfADC elements. Case 2 will
	/// need to create new elements. For dark-hit-only channels, a threshold
	/// is applied so only summed channels with at least BCAL_MIN_ENERGY_MEV
	/// equivalent deposited will have non-NULL histograms. All summed channels
	/// will have an entry in bcalfADC upon exit though, even if they have NULL
	/// histos.
	
	// In order to avoid potentially expensive mutex locks, we grab a couple
	// of histograms from the pool now and reuse them for the dark-hits-only
	// channels until we actually find one that has enough dark hits for us 
	// to keep it. Then we'll grab another from the pool. Since most channels
	// do NOT have signal and also do NOT have enough dark hits to make a
	// signal above threshold, it is worthwhile to reuse these temporary
	// histograms for multiple summed cells.

	if(NO_DARK_PULSES)return;

	DHistogram *hup_tmp = GetHistoFromPool();
	DHistogram *hdn_tmp = GetHistoFromPool();
	
	// Loop over all fADC readout cells
	for(int imodule=1; imodule<=DBCALGeometry::NBCALMODS; imodule++){

		// inner
		for(int fADC_lay=1; fADC_lay<=DBCALGeometry::NBCALLAYSIN; fADC_lay++){
			for(int fADC_sec=1; fADC_sec<=DBCALGeometry::NSUMSECSIN; fADC_sec++){
			
				// Use cellId(...) to convert fADC layer and sector into fADCId
				// (see DBCALGeometry::fADCId)
				int fADCId = DBCALGeometry::cellId(imodule, fADC_lay, fADC_sec);
				
				// Get SumSpectra object if it already exists or create new one 
				// if it doesn't.
				SumSpectra &sumspectra = bcalfADC[fADCId];
				
				// If we just created this, then grab some histos from pool.
				// n.b. we try first to use the pre-grabbed histos
            if(sumspectra.hup == NULL){
               if(hup_tmp == NULL){
                  sumspectra.hup = GetHistoFromPool();
               }else{
                  sumspectra.hup = hup_tmp;
                  hup_tmp = NULL;
               }
            }
            if(sumspectra.hdn == NULL){
               if(hdn_tmp == NULL){
                  sumspectra.hdn = GetHistoFromPool();
               }else{
                  sumspectra.hdn = hdn_tmp;
                  hdn_tmp = NULL;
               }
            }

				// Get the number of SiPMs in this readout channel.
				int NSiPMs = DBCALGeometry::NSiPMs(fADCId);
				
				// Add in dark hits for any additional SiPMs.
				int NSiPMs_no_signal =  NSiPMs - (int)sumspectra.cellspectra.size();
				if(NSiPMs_no_signal>0){
					double dEup = AddDarkHitsToOne(sumspectra.hup, NSiPMs_no_signal);
					double dEdn = AddDarkHitsToOne(sumspectra.hdn, NSiPMs_no_signal);
					
					// If this is a dark-hit-only channel, check if enough
					// equivalent energy was added to make it worthwhile to keep
					// these. Otherwise, recycle the histos for the next iteration
					// to avoid the expensive convolution and checking for thresholds
					// stages later.
					if(NSiPMs_no_signal == NSiPMs){
						if(dEup < BCAL_MIN_ENERGY_MEV){
							hup_tmp = sumspectra.hup;
							hup_tmp->Reset();
							sumspectra.hup = NULL;
						}
						if(dEdn < BCAL_MIN_ENERGY_MEV){
							hdn_tmp = sumspectra.hdn;
							hdn_tmp->Reset();
							sumspectra.hdn = NULL;
						}
					}
				}
			}
		}
	}
	
	// Return temporary histos to pool
	if(hup_tmp)ReturnHistoToPool(hup_tmp);
	if(hdn_tmp)ReturnHistoToPool(hdn_tmp);
	
}

//-----------
// ApplyElectronicPulseShape
//-----------
void ApplyElectronicPulseShape(map<int, SumSpectra> &bcalfADC)
{
	/// Loop over the SumSpectra objects and convolute any histograms
	/// with the electronic pulse shape.

	map<int, SumSpectra>::iterator iter = bcalfADC.begin();
	for(; iter!=bcalfADC.end(); iter++){
		
		SumSpectra &sumspectra = iter->second;
		
		if(sumspectra.hup)ApplyElectronicPulseShapeOneHisto(sumspectra.hup);
		if(sumspectra.hdn)ApplyElectronicPulseShapeOneHisto(sumspectra.hdn);
	}
}

//-----------
// ApplyElectronicPulseShapeOneHisto
//-----------
void ApplyElectronicPulseShapeOneHisto(DHistogram *h)
{
	/// Convolute the given histogram with the electronic pulse
	/// shape, overwriting the input histograms contents.
	
	// Make copy of input histogram and then reset it
	DHistogram tmp(*h);
	h->Reset();
	
	int Nbins = h->GetNbins();
   double bin_width = h->GetBinWidth();
   
   // Get pointers to the bin content arrays directly.
   // This will allow the two nested loops below to
   // work much more efficiently.
   // n.b. DHistogram content arrays allocate an extra
   // element at the beginning so it can use bin=1 as the
   // first bin. Thus, we need to immediately increment these
   // pointers to get to bin 1 (that first bin may eventually
   // be used for underflow).
   double* tmp_content = tmp.GetContentPointer();
   double* h_content = h->GetContentPointer();
   tmp_content++;
   h_content++;

	// Loop over bins of input histo
	for(int ibin=1; ibin<=Nbins; ibin++, tmp_content++){
		double A = *tmp_content;
      if(A==0.0)continue; // no need to continue for empty bins
		
		// Get pointer to first element to row in the 
		// pre-calculated matrix that converts the histo
		// from its input form to its convoluted form
		int idx1 = (ibin-1)*Nbins;
		
		// Loop over bins of output histo
      double *h_content_tmp = h_content;
      double *pulse_shape = &BCAL_PULSE_SHAPE_MATRIX[idx1];
		for(int jbin=1; jbin<=Nbins; jbin++, h_content_tmp++, pulse_shape++){

			double weight = A*(*pulse_shape);
         *h_content_tmp += weight;
		}
	}
}

//-----------
// FindHits
//-----------
void FindHits(double thresh_mV, map<int, SumSpectra> &bcalfADC, map<int, fADCHitList> &fADCHits)
{
	/// Loop over SumSpectra objects and find places where it crosses threshold
	map<int, SumSpectra>::iterator iter = bcalfADC.begin();
	for(; iter!=bcalfADC.end(); iter++){
		
		int fADCId = iter->first;
		SumSpectra &sumspectra = iter->second;
		
		// Find threshold crossings for both upstream and downstream readout channels
		vector<fADCHit> uphits;
		vector<fADCHit> dnhits;
		if(sumspectra.hup)FindHitsOneHisto(thresh_mV, sumspectra.hup, uphits);
		if(sumspectra.hdn)FindHitsOneHisto(thresh_mV, sumspectra.hdn, dnhits);
		
		// If at least one readout channel has a hit, add the readout cell to fADCHits
		if(uphits.size()>0 || dnhits.size()>0){
			fADCHitList &hitlist = fADCHits[fADCId];
			
			// The module, fADC layer, and fADC sector are encoded in fADCId
			// (n.b. yes, these are the same methods used for extracting 
			// similar quantities from the cellId.)
			hitlist.module = DBCALGeometry::module(fADCId);
			hitlist.fADClayer = DBCALGeometry::layer(fADCId);
			hitlist.fADCsector = DBCALGeometry::sector(fADCId);
			
			hitlist.uphits = uphits;
			hitlist.dnhits = dnhits;
		}
	}
}

//-----------
// FindHitsOneHisto
//-----------
void FindHitsOneHisto(double thresh_mV, DHistogram *h, vector<fADCHit> &hits)
{
	/// Look through the given histogram and find places where the signal
	/// size crosses the given threshold. The signal is integrated before
	/// and after the crossing time by 20ns and 180ns respectively. The 
	/// times are approximate as they are measured in number of bins by
	/// dividing by the bin width of the histo.
	///
	/// The amplitude of the signal is scaled down to account for the pre_amp
	/// gain difference between TDC and fADC signal splits. It is also converted
	/// into units of fADC counts assuming 4096 counts per 2V.
	///
	/// n.b. no quantization of the ADC counts is done at this stage since
	/// it is unclear if the conversions are correct. This can be done
	/// easily enough at a later stage after any additional scaling factors
	/// are applied.
	
	double bin_width = h->GetBinWidth();
	int Nbins_before = (int)(BCAL_FADC_INTEGRATION_WINDOW_PRE/bin_width);
	int Nbins_after = (int)(BCAL_FADC_INTEGRATION_WINDOW_POST/bin_width);

	// Loop 
	int start_bin = 1;
	int Nbins = h->GetNbins();
	while(start_bin<=Nbins){
		int ibin = h->FindFirstBinAbove(thresh_mV, start_bin);
		if(ibin<1 || ibin>Nbins)break;
		
		// Signal time. We simply use the center of the bin as the time.
		// In reality, the FPGA will implement some algorithm to extract
		// it from the surrounding samples, thus having higher accuracy
		// than the bin width. Here, the histograms should have smaller bins
		// than the 4ns the fADC uses.
		double t = h->GetBinCenter(ibin);
		
		// Calculate integration limits for signal amplitude
		int istart = ibin - Nbins_before;
		int iend = ibin + Nbins_after;
		if(istart<1)istart=1;
		if(iend>Nbins)iend = Nbins;

		// Integrate signal
		// n.b. we use the start_bin varible so it is left
		// pointing to the end of the integration window which
		// is where we start looking for the next hit on the next iteration 
		double integral = 0.0;
		for(start_bin=istart; start_bin<=iend; start_bin++){
			integral += h->GetBinContent(start_bin);
		}
		
		// Scale the integral by the ratio of bin widths to get it in
		// units of mV * fADC bins
		integral *= bin_width/4.0; // the fADC250 has 4ns wide samples
		
		// The histogram should have the signal size for the TDC, but the
		// fADC leg will actually have a smaller size since the pre-amp gain
		// will be set differently. Scale the integral down here.
		double preamp_gain_tdc = 5.0;
		integral /= preamp_gain_tdc;
		
		// Expect 4906 counts/2V
		double fADC = integral*4096.0/2000.0; // 2V = 2000mV

		// Store hit in container
		hits.push_back(fADCHit(fADC,t));
	}
}

//-----------
// CopyBCALHitsToHDDM
//-----------
void CopyBCALHitsToHDDM(map<int, fADCHitList> &fADCHits, s_HDDM_t *hddm_s)
{
	/// Loop over fADCHitList objects and copy the fADC hits into the HDDM tree.
	///
	/// This will copy all of the hits found into the first physicsEvent found
	/// in the HDDM file. Note that the hits were formed from data that may
	/// have been combined from several physicsEvent structures in the HDDM
	/// event. No attempt is made to keep track of this so all hits are thrown
	/// into only a single physicsEvent.

	// Get pointer to first physicsEvent->hitView structure
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	if(PE->mult<1)return;
	
	s_HitView_t *hits = PE->in[0].hitView;

	// Delete any existing bcalfADCCell structures (along with their bcalfADCUpHit
	// and bcalfADCDownHit structures).
	if(hits->barrelEMcal->bcalfADCCells!=HDDM_NULL){
		for(unsigned int i=0; i<hits->barrelEMcal->bcalfADCCells->mult; i++){
			s_BcalfADCUpHits_t *uphits = hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCUpHits;
			s_BcalfADCDownHits_t *downhits = hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCDownHits;
			if(uphits!=NULL && uphits!=HDDM_NULL){
				free(uphits);
			}
			if(downhits!=NULL && downhits!=HDDM_NULL){
				free(downhits);
			}
		}
		free(hits->barrelEMcal->bcalfADCCells);
	}
	
	// Make sure bcalfADCCells pointer is empty in case we return early below
	hits->barrelEMcal->bcalfADCCells = (s_BcalfADCCells_t*)HDDM_NULL;

	// If we have no cells over threshold, then bail now.
	if(fADCHits.size()==0) return;
	
	// Create bcalfADCCell structures to hold all of our hits
	hits->barrelEMcal->bcalfADCCells = make_s_BcalfADCCells(fADCHits.size());
	unsigned int &mult = hits->barrelEMcal->bcalfADCCells->mult;
	map<int, fADCHitList>::iterator it = fADCHits.begin();
	for(mult=0; mult<fADCHits.size(); mult++, it++){
	
		// Get pointer to this fADC cell in the HDDM tree
		s_BcalfADCCell_t *fADCcell = &hits->barrelEMcal->bcalfADCCells->in[mult];
		
		// Get pointer to our fADC cell information that needs to be copied to HDDM
		fADCHitList &hitlist = it->second;
		
		// Copy cell information to HDDM
		fADCcell->module = hitlist.module;
		fADCcell->layer  = hitlist.fADClayer;
		fADCcell->sector = hitlist.fADCsector;
		
		// If we have any upstream hits, copy them into HDDM
		if(hitlist.uphits.size()>0){
			fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(hitlist.uphits.size());
			fADCcell->bcalfADCUpHits->mult = 0;
			for(unsigned int i=0; i<hitlist.uphits.size(); i++){
				s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[i];
				
				fadcuphit->E = hitlist.uphits[i].E;
				fadcuphit->t = hitlist.uphits[i].t;
				fADCcell->bcalfADCUpHits->mult++;
			}
		}else{
			fADCcell->bcalfADCUpHits = (s_BcalfADCUpHits_t*)HDDM_NULL;
		}

		// If we have any downstream hits, copy them into HDDM
		if(hitlist.dnhits.size()>0){
			fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(hitlist.dnhits.size());
			fADCcell->bcalfADCDownHits->mult = 0;
			for(unsigned int i=0; i<hitlist.dnhits.size(); i++){
				s_BcalfADCDownHit_t *fadcdnhit = &fADCcell->bcalfADCDownHits->in[i];
				
				fadcdnhit->E = hitlist.dnhits[i].E;
				fadcdnhit->t = hitlist.dnhits[i].t;
				fADCcell->bcalfADCDownHits->mult++;
			}
		}else{
			fADCcell->bcalfADCDownHits = (s_BcalfADCDownHits_t*)HDDM_NULL;
		}

	} // fADCHits

}

//-----------
// GetHistoFromPool
//-----------
DHistogram* GetHistoFromPool(void)
{
	DHistogram *h = NULL;

	pthread_mutex_lock(&histo_pool_mutex);
	if(histo_pool.size()>0){
		// Take histo out of pool
		h = histo_pool.front();
		histo_pool.pop();
	}
	pthread_mutex_unlock(&histo_pool_mutex);
	
	// Do this outside of mutex lock to speed things up
	if(h==NULL){
		// There wasn't a histo in the pool. Create a new one
		h = new DHistogram(400, -100.0, 300.0);
	}else{
		// Histo came from pool. Reset it
		h->Reset();
	}

	return h;
}

//-----------
// ReturnHistoToPool
//-----------
void ReturnHistoToPool(DHistogram *h)
{
	pthread_mutex_lock(&histo_pool_mutex);
	if(h!=NULL){
		// Limit size of pool
		if(histo_pool.size() >= 200){
			delete h;
		}else{
			histo_pool.push(h);
		}
	}
	pthread_mutex_unlock(&histo_pool_mutex);
}

//---------------
// MakeTSpline
//---------------
TSpline* MakeTSpline(void)
{
	// These points are taken by mapping out figure 6
	// of GlueX-doc-1795 (the "5ns" rise time pulse).

	double t[] ={
		-7.00,
		-4.92,
		-0.75,
		4.24,
		9.37,		// 5
		
		9.71,
		10.05,
		10.39,
		10.73,
		
		11.07,
		12.03,
		12.92,
		14.08,
		16.05,	// 10
		18.56,
		21.22,
		24.31,
		28.44,
		32.92,	// 15
		38.64,
		42.07,
		//47.10,
		//52.59,
		//58.22,	// 20
		67.95,
		71.74,
		79.85,
		120.0,
		125.0,	// 25
		130.0,
		160.0,
		180.0,
		-1000.0 // flag end of data
	};
	
	double mV[]={
		0,
		0,
		0,
		0,
		0-1,					// 5

		-1.2,
		-1.4,
		-1.6,
		-1.8,

		0-2,
		-49.16039931,
		-159.0483507,
		-375.9324653,
		-711.3798959,	// 10
		-902.2379167,
		-988.9915625,
		-1020.801233,
		-960.0736806,
		-850.1857292,	// 15
		-694.0291667,
		-601.4919445,
		//-539.29,
		//-419.88,
		//-300.46,			// 20
		-142.53,
		-100.15,
		-61.63,
		-8.0,
		-4.0,				// 25
		-2.0,
		-1.0,
		0.0,
		-1000.0 // flag end of data
	};

	// Find number points
	int Npoints=0;
	while(true){
		if(t[Npoints]==-1000.0 && mV[Npoints]==-1000.0)break;
		Npoints++;
	}
	
	TSpline *s = new TSpline3("spline", t, mV, Npoints);
	
	return s;
}

//-----------
// SaveDebugHistos
//-----------
void SaveDebugHistos(map<bcal_index, CellSpectra> &SiPMspectra, const char *suffix, const char *yunits)
{
	/// Create ROOT histograms out of the CellSpectra objects. These would
	/// correspond to the individual SiPMs. These will be saved in the smear.root 
	/// file inside of a directory structure to allow debugging.
	///
	/// WARNING: This will save several histograms for every event!
	/// Usually, you will only want to run mcsmear for one or two
	/// events when this is being called.
	
	// This should probably only be called when running single-threaded
	// but better safe than sorry!
	pthread_mutex_lock(&root_mutex);

	// Save the current ROOT directory so we can restore it before returning
	TDirectory *savedir = gDirectory;
	
	map<bcal_index, CellSpectra>::iterator iter=SiPMspectra.begin();
	for(; iter!=SiPMspectra.end(); iter++){
		const bcal_index &idx = iter->first;
		CellSpectra &cellspectra = iter->second;
		
		// Skip empty spectra objects
		if(cellspectra.hup==NULL && cellspectra.hdn==NULL)continue;
		bool allempty=true;
		if(cellspectra.hup && cellspectra.hup->Integral()!=0.0)allempty = false;
		if(cellspectra.hdn && cellspectra.hdn->Integral()!=0.0)allempty = false;
		if(allempty)continue;
		
		// Extract location info
		int module = idx.module;
		int layer = idx.layer;
		int sector = idx.sector;
		
		// Create a directory to hold this readout cell's histos
		char dirname[256];
		sprintf(dirname, "SiPM_m%02dl%ds%d%s", module, layer, sector, suffix);
		
		// Directory may already exist since separate entries may be kept
		// for upstream and downstream. Check if it exists first and only
		// create it if necessary.
		savedir->cd();
		TDirectory *rcdir = (TDirectory*)savedir->FindObject(dirname);
		if(!rcdir)rcdir = savedir->mkdir(dirname);
		rcdir->cd();
		
		// Make ROOT histograms out of upstream and downstream histos
		char hname[256];
		if(cellspectra.hup && cellspectra.hup->Integral()!=0.0){
			sprintf(hname, "Upstream summed cell mod=%d lay=%d sec=%d", module, layer, sector);
			TH1D *rh = cellspectra.hup->MakeTH1D("hup", hname);
			rh->SetXTitle("Time (ns)");
			rh->SetYTitle(yunits);
		}
		if(cellspectra.hdn && cellspectra.hdn->Integral()!=0.0){
			sprintf(hname, "Downstream summed cell mod=%d lay=%d sec=%d", module, layer, sector);
			TH1D *rh = cellspectra.hdn->MakeTH1D("hdn", hname);
			rh->SetXTitle("Time (ns)");
			rh->SetYTitle(yunits);
		}
		
	}
	
	// Restore ROOT working directory
	savedir->cd();

	pthread_mutex_unlock(&root_mutex);
	
}

//-----------
// SaveDebugHistos
//-----------
void SaveDebugHistos(map<int, SumSpectra> &bcalfADC, const char *suffix)
{
	/// Create ROOT histograms out of the histograms in the SumSpectra
	/// objects. These can be either electronics pulse shapes or attenuated
	/// energy spectra (depending upon what stage of the smearing has been
	/// done to the objects given us). These will be saved in the smear.root 
	/// file inside of a directory structure to allow debugging.
	///
	/// WARNING: This will save several histograms for every event!
	/// Usually, you will only want to run mcsmear for one or two
	/// events when this is being called.
	
	// This should probably only be called when running single-threaded
	// but better safe than sorry!
	pthread_mutex_lock(&root_mutex);

	// Save the current ROOT directory so we can restore it before returning
	TDirectory *savedir = gDirectory;
	
	map<int, SumSpectra>::iterator iter=bcalfADC.begin();
	for(; iter!=bcalfADC.end(); iter++){
		int fADCId = iter->first;
		SumSpectra &sumspectra = iter->second;
		
		// Skip empty spectra objects
		if(sumspectra.hup==NULL && sumspectra.hdn==NULL)continue;
		bool allempty=true;
		if(sumspectra.hup && sumspectra.hup->Integral()!=0.0)allempty = false;
		if(sumspectra.hdn && sumspectra.hdn->Integral()!=0.0)allempty = false;
		if(allempty)continue;
		
		// Extract location info
		int module = DBCALGeometry::module(fADCId);
		int fADC_layer = DBCALGeometry::layer(fADCId);
		int fADC_sector = DBCALGeometry::sector(fADCId);
		
		// Create a directory to hold this readout cell's histos
		char dirname[256];
		sprintf(dirname, "Sum_m%02dl%ds%d%s", module, fADC_layer, fADC_sector, suffix);
		TDirectory *rcdir = savedir->mkdir(dirname);
		rcdir->cd();
		
		// Make ROOT histograms out of upstream and downstream histos
		char hname[256];
		if(sumspectra.hup && sumspectra.hup->Integral()!=0.0){
			sprintf(hname, "Upstream summed cell mod=%d lay=%d sec=%d", module, fADC_layer, fADC_sector);
			TH1D *rh = sumspectra.hup->MakeTH1D("hup", hname);
			rh->SetXTitle("Time (ns)");
			rh->SetYTitle("Amplitude (mV)");
		}
		if(sumspectra.hdn && sumspectra.hdn->Integral()!=0.0){
			sprintf(hname, "Downstream summed cell mod=%d lay=%d sec=%d", module, fADC_layer, fADC_sector);
			TH1D *rh = sumspectra.hdn->MakeTH1D("hdn", hname);
			rh->SetXTitle("Time (ns)");
			rh->SetYTitle("Amplitude (mV)");
		}
		
	}
	
	// Restore ROOT working directory
	savedir->cd();

	pthread_mutex_unlock(&root_mutex);
	
}

#else
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//                                          BELOW HERE IS DEPRECATED!!

// $Id: smear.cc 7650 2011-03-29 22:52:30Z shepherd $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
using namespace std;

#include <BCAL/DBCALGeometry.h>

#include <math.h>
#include "units.h"
#include "HDDM/hddm_s.h"
#include <TF1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TMatrix.h>

#include "DBCALReadoutChannel.h"

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

// From SampleGaussian.cc
double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);

// Defined in this file
double bcalSamplingSmear( double E );
double bcalTimeSmear(double t, double e);
int getDarkHits();
void bcalInit(void);
void CalculateThresholds(void);
TH1D* CombineSiPM_PDFs(const char *hname, TH1D *single, int NSiPM);


// The following are set in bcalInit below
bool BCAL_INITIALIZED = false;
double BCAL_mevPerPE=0.0;
double BCAL_UNATTENUATE_TO_CENTER = 0.0;
double BCAL_sigma_singlePE_sq = 0.0;
double BCAL_sigma_ped_sq = 0.0;

// since the size of a readout unit can vary even among the inner layers, we need to store the following values for each layer individually
double BCAL_inner_thresh[DBCALGeometry::NBCALLAYSIN];
double BCAL_outer_thresh[DBCALGeometry::NBCALLAYSOUT];

//number of SiPM's summed in one readout unit for each layer
int Nsum_inner[DBCALGeometry::NBCALLAYSIN];
int Nsum_outer[DBCALGeometry::NBCALLAYSOUT];


map<int, DBCALReadoutChannel> bcal_fADCs; // key is DBCALGeometry::fADCId()



// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn OFF the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_THRESHOLD_CUT;


// setup response parameters
extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
extern double BCAL_XTALK_FRACT;                 // 0.157  (from calibDB BCAL/bcal_parms)
extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms) 
extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)
extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)



//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
	DBCALGeometry bcalGeom;

	// Initialize BCAL globals on first call
	if(!BCAL_INITIALIZED)bcalInit();
	
	// Container to keep track of the elements in bcal_fADCs that
	// have hits so we can sparsely clear them at the end of this
	set<int> fADC_ids_with_hits;
	
	// Loop over PhysicsEvents
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	for(unsigned int iphysics_event=0; iphysics_event<PE->mult; iphysics_event++){
	
		// Get pointer to HDDM hits structure
		s_HitView_t *hits = PE->in[iphysics_event].hitView;
		
		// Make sure HDDM stuctures exist.
		// In the case of no real BCAL hits, we may still want to emit
		// dark hit only events. In this case, we must create the BCAL 
		// tree here.
		if(hits==HDDM_NULL)hits = PE->in[iphysics_event].hitView = make_s_HitView();
		if(hits->barrelEMcal == HDDM_NULL)hits->barrelEMcal = make_s_BarrelEMcal();
		if(hits->barrelEMcal->bcalCells == HDDM_NULL)hits->barrelEMcal->bcalCells = make_s_BcalCells(0);

		// Loop over GEANT hits in BCAL
		s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
		for(unsigned int j=0; j<cells->mult; j++){
			s_BcalCell_t *cell = &cells->in[j];

			// If data exists in HDDM tree for SiPM hits, then delete it
			// so we can write in new hits below.
			if(cell->bcalSiPMUpHits!=HDDM_NULL)free(cell->bcalSiPMUpHits);
			if(cell->bcalSiPMDownHits!=HDDM_NULL)free(cell->bcalSiPMDownHits);

			// We make one SiPM hit for every GEANT hit in both the up and
			// down stream ends
			unsigned int Nbcalhits = cell->bcalHits->mult;
			cell->bcalSiPMUpHits = make_s_BcalSiPMUpHits(Nbcalhits);
			cell->bcalSiPMDownHits = make_s_BcalSiPMDownHits(Nbcalhits);

			cell->bcalSiPMUpHits->mult = Nbcalhits;
			cell->bcalSiPMDownHits->mult = Nbcalhits;			

			// Loop over GEANT hits in BCAL for this cell
			for(unsigned int k=0; k<Nbcalhits; k++){
				s_BcalHit_t *bcalhit = &cell->bcalHits->in[k];
				s_BcalSiPMUpHit_t *bcaluphit = &cell->bcalSiPMUpHits->in[k];
				s_BcalSiPMDownHit_t *bcaldownhit = &cell->bcalSiPMDownHits->in[k];
	
				// Distance from location of energy deposition to each end of BCAL
				double upDist   = ( bcalGeom.BCALFIBERLENGTH / 2 ) + bcalhit->zLocal;
				double downDist = ( bcalGeom.BCALFIBERLENGTH / 2 ) - bcalhit->zLocal;
	
				// sampling fluctuations are correlated between ends
				double smearedE = bcalSamplingSmear( bcalhit->E );
				double upEnergy   = smearedE * exp( -upDist / bcalGeom.ATTEN_LENGTH );
				double downEnergy = smearedE * exp( -downDist / bcalGeom.ATTEN_LENGTH );

				// independently smear time for both ends -- time smearing 
				// parameters come from data taken with beam at the center of 
				// the module so there is an implicit exp( ( -L / 2 ) / lambda ) 
				// that needs to be canceled out since we are working
				// at this stage with attenuated energies
				double smearedtUp   = bcalTimeSmear( bcalhit->t,   upEnergy*BCAL_UNATTENUATE_TO_CENTER );
				double smearedtDown = bcalTimeSmear( bcalhit->t, downEnergy*BCAL_UNATTENUATE_TO_CENTER );

				// To get a realistic distribution we have to convert the smeared energy into 
				// an integer number of photoelectrons. We'll add the dark hits to that and
				// then calculate a sigma in energy based on the total number of pixels that
				// fired. The number of photoelectrons will be converted back into energy
				// and then smeared using this sigma.
				double upPE = floor(0.5 + upEnergy/(BCAL_mevPerPE * k_MeV)) + getDarkHits();
				double downPE = floor(0.5 + downEnergy/(BCAL_mevPerPE * k_MeV)) + getDarkHits();

				double sigma_up = sqrt(BCAL_sigma_ped_sq + upPE*BCAL_sigma_singlePE_sq);
				double sigma_down = sqrt(BCAL_sigma_ped_sq + downPE*BCAL_sigma_singlePE_sq);
				
				// Convert integer # of PE back into energy
				upEnergy = upPE * BCAL_mevPerPE * k_MeV;
				downEnergy = downPE * BCAL_mevPerPE * k_MeV;

				// Add in pedestal widths due to SiPMs
				upEnergy += SampleGaussian(sigma_up);
				downEnergy += SampleGaussian(sigma_down);
	
				// now offset times for propagation distance
				double upTime = smearedtUp + upDist / bcalGeom.C_EFFECTIVE;
				double downTime = smearedtDown + downDist / bcalGeom.C_EFFECTIVE;
	
				// If energy is smeared to negative or time is nan, set to 0.
				if(upEnergy <= 0 || !finite(upTime)){
					upEnergy = 0;
					upTime = 0;
				}
				if(downEnergy <= 0 || !finite(downTime)){
					downEnergy = 0;
					downTime = 0;
				}

				// Record upstream SiPM values in HDDM
				bcaluphit->E = upEnergy;
				bcaluphit->t = upTime;

				// Record downstream SiPM values in HDDM
				bcaldownhit->E = downEnergy;
				bcaldownhit->t = downTime;
			
				// Remember these hits in appropriate readout channel
				int fADC_Id = bcalGeom.fADCId( cell->module, cell->layer, cell->sector);
				bcal_fADCs[fADC_Id].uphits.push_back(bcaluphit);
				bcal_fADCs[fADC_Id].downhits.push_back(bcaldownhit);
				fADC_ids_with_hits.insert(fADC_Id);
			} //k (bcal hits)
		} //j (cells)

		// Above, we looped over GEANT hits to create SiPM hits. Below we loop
		// over fADC (readout) channels which may contain multiple SiPMs.
		// For each channel, we determine how many hits occurred in
		// BCAL_TIME_WINDOW, including cumulative dark hits. The number of 
		// cells over threshold are counted and the values recorded so that
		// we can create the correct number of bcalfADCCell stuctures in HDDM
		// later.

		// Loop over readout channels
		set<int> fADC_ids_over_thresh;
		map<int, DBCALReadoutChannel>::iterator iter;
		for(iter = bcal_fADCs.begin(); iter!=bcal_fADCs.end(); iter++){

			int fADC_id = iter->first;
			DBCALReadoutChannel &bcalfADC = iter->second;
			vector<s_BcalSiPMUpHit_t*> &uphits = bcalfADC.uphits;
			vector<s_BcalSiPMDownHit_t*> &downhits = bcalfADC.downhits;

			//--------- Upstream -------------
			// Initialize
			bcalfADC.Eup = 0.0;
			bcalfADC.tup = 0.0;

			// Sum upstream hits
			for(unsigned int j=0; j<uphits.size(); j++){
				bcalfADC.Eup += uphits[j]->E;
				bcalfADC.tup += uphits[j]->t * uphits[j]->E; // energy weighted time for sum
			}
	
			// Add in dark pulses for upstream SiPMs not hit, but that are included in sum
			unsigned int Ndark_channels_up = bcalfADC.NSiPM - uphits.size();
			if(uphits.size() > bcalfADC.NSiPM)Ndark_channels_up = 0;

			for(unsigned int j=0; j<Ndark_channels_up; j++){
				double darkPE = getDarkHits();
				double sigma_up = sqrt(BCAL_sigma_ped_sq + darkPE*BCAL_sigma_singlePE_sq);				
				double Edark = (darkPE * BCAL_mevPerPE * k_MeV) + SampleGaussian(sigma_up);
				double tdark = SampleRange( -0.25 * BCAL_INTWINDOW_NS, 0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
				bcalfADC.Eup += Edark;
				bcalfADC.tup += tdark * Edark;
			}
			
			// normalize time
			bcalfADC.tup /= bcalfADC.Eup;

			//--------- Downstream -------------
			// Initialize
			bcalfADC.Edown = 0.0;
			bcalfADC.tdown = 0.0;

			// Sum downstream hits
			for(unsigned int j=0; j<downhits.size(); j++){
				bcalfADC.Edown += downhits[j]->E;
				bcalfADC.tdown += downhits[j]->t * downhits[j]->E; // energy weighted time for sum
			}
			
			// Add in dark pulses for downstream SiPMs not hit, but that are included in sum
			unsigned int Ndark_channels_down = bcalfADC.NSiPM - downhits.size();
			if(downhits.size() > bcalfADC.NSiPM)Ndark_channels_down = 0;
			for(unsigned int j=0; j<Ndark_channels_down; j++){
				double darkPE = getDarkHits();
				double sigma_down = sqrt(BCAL_sigma_ped_sq + darkPE*BCAL_sigma_singlePE_sq);				
				double Edark = (darkPE * BCAL_mevPerPE * k_MeV) + SampleGaussian(sigma_down);
				double tdark = SampleRange( -0.25 * BCAL_INTWINDOW_NS, 0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
				bcalfADC.Edown += Edark;
				bcalfADC.tdown += tdark * Edark;
			}
			
			// normalize time
			bcalfADC.tdown /= bcalfADC.Edown;


			// If either the upstream or downstream hit is over threshold, then
			// a bcalfADCCell will need to be created so remember this id if needed.
			if((bcalfADC.Eup > bcalfADC.threshold) || (bcalfADC.Edown > bcalfADC.threshold)){
				fADC_ids_over_thresh.insert(fADC_id);
			}
		} // readout channel
		
		
		// At this point we have summed all SiPMs with dark hits included. Any cell with
		// either the upstream or downstream fADC over threshold has its id stored in
		// the fADC_ids_over_thresh container. We now need to create bcalfADCCell,
		// bcalfADCUpHit, and bcalfADCDownHit structures within the HDDM tree to hold
		// the fADC values.
		
		
		// Delete any existing bcalfADCCell structures (along with their bcalfADCUpHit
		// and bcalfADCDownHit structures).
		if(hits->barrelEMcal->bcalfADCCells!=HDDM_NULL){
			for(unsigned int i=0; i<hits->barrelEMcal->bcalfADCCells->mult; i++){
				s_BcalfADCUpHits_t *uphits = hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCUpHits;
				s_BcalfADCDownHits_t *downhits = hits->barrelEMcal->bcalfADCCells->in[i].bcalfADCDownHits;
				if(uphits!=NULL && uphits!=HDDM_NULL){
					free(uphits);
				}
				if(downhits!=NULL && downhits!=HDDM_NULL){
					free(downhits);
				}
			}
			free(hits->barrelEMcal->bcalfADCCells);
		}
		
		// Make sure bcalfADCCells pointer is empty in case we return early below
		hits->barrelEMcal->bcalfADCCells = (s_BcalfADCCells_t*)HDDM_NULL;

		// If we have no cells over threshold, then bail now.
		if(fADC_ids_over_thresh.size()==0) continue; // next iphysics_event
		
		// Create bcalfADCCell structures to hold all of our hits
		hits->barrelEMcal->bcalfADCCells = make_s_BcalfADCCells(fADC_ids_over_thresh.size());
		unsigned int &mult = hits->barrelEMcal->bcalfADCCells->mult;
		set<int>::iterator it = fADC_ids_over_thresh.begin();
		for(mult=0; mult<fADC_ids_over_thresh.size(); mult++, it++){
		
			// Get pointer to this fADC cell in the HDDM tree
			s_BcalfADCCell_t *fADCcell = &hits->barrelEMcal->bcalfADCCells->in[mult];
			
			// Get pointer to our fADC cell information that needs to be copied to HDDM
			DBCALReadoutChannel &bcalfADC = bcal_fADCs[*it];
			
			fADCcell->module = bcalfADC.module;
			fADCcell->layer  = bcalfADC.layer;
			fADCcell->sector = bcalfADC.sector;

			// Upstream hit
			if(bcalfADC.Eup > bcalfADC.threshold){
				fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
				fADCcell->bcalfADCUpHits->mult = 1;
				s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
				
				fadcuphit->E = bcalfADC.Eup;
				fadcuphit->t = bcalfADC.tup;
			}else{
				fADCcell->bcalfADCUpHits = (s_BcalfADCUpHits_t*)HDDM_NULL;
			}

			// Downstream hit
			if(bcalfADC.Edown > bcalfADC.threshold){
				fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
				fADCcell->bcalfADCDownHits->mult = 1;
				s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];
				
				fadcdownhit->E = bcalfADC.Edown;
				fadcdownhit->t = bcalfADC.tdown;
			}else{
				fADCcell->bcalfADCDownHits = (s_BcalfADCDownHits_t*)HDDM_NULL;
			}
		} // fADC_ids_over_thresh
	} // iphysics_event

	
	// Clear the bcal_fADCs elements that have hits.
	// We keep track above so we can do a sparse clear here to
	// save time clearing every fADC array for every event.
	set<int>::iterator iter=fADC_ids_with_hits.begin();
	for(; iter!=fADC_ids_with_hits.end(); iter++){
		bcal_fADCs[*iter].Clear();
	}
}

//-----------
// bcalSamplingSmear
//-----------
double bcalSamplingSmear( double E )
{
	if(NO_E_SMEAR)return E;

    double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( E ) + BCAL_SAMPLINGCOEFB;
    
    return( E * (1.0 + SampleGaussian(sigmaSamp)) );
}

//-----------
// bcalTimeSmear
//-----------
double bcalTimeSmear( double t, double E )
{
	if(NO_T_SMEAR)return t;

  double sigmaT = BCAL_TIMEDIFFCOEFA / sqrt( E ) + BCAL_TIMEDIFFCOEFB;

  return( t + SampleGaussian(sigmaT) );
}

//-----------
// getDarkHits
//-----------
int getDarkHits()
{
	if(NO_DARK_PULSES)return 0;

  int darkPulse = (int)floor(0.5+SamplePoisson( BCAL_DARKRATE_GHZ* BCAL_INTWINDOW_NS ));

  int xTalk = (int)floor(0.5+SamplePoisson( (double)darkPulse * BCAL_XTALK_FRACT ));

  return( xTalk + darkPulse );
}

//-----------
// bcalInit
//-----------
void bcalInit(void)
{
	// Calculate the inner and outer thresholds based on
	// the parameters from the calibDB. The threshold values
	// are left in the global variables BCAL_inner_thresh
	// and BCAL_outer_thresh.
	CalculateThresholds();

	for (int i=0;i<DBCALGeometry::NBCALLAYSIN;i++) {
	  if(NO_THRESHOLD_CUT) BCAL_inner_thresh[i]=0.0;
	  cout<<"BCAL threshold for layer " << i+1 << ": "<<BCAL_inner_thresh[i]*1000.0<<" MeV"<<" ("<<Nsum_inner[i]<<" SiPMs summed)"<<endl;
	}

	for (int i=0;i<DBCALGeometry::NBCALLAYSOUT;i++) {
	  if(NO_THRESHOLD_CUT) BCAL_outer_thresh[i]=0.0;
	  cout<<"BCAL threshold for layer " << i+1+DBCALGeometry::NBCALLAYSIN << ": "<<BCAL_outer_thresh[i]*1000.0<<" MeV"<<" ("<<Nsum_outer[i]<<" SiPMs summed)"<<endl;
	}
	
	// Factor to unattenuate the energy at the end to what it *would* be if it came
	// from the center of the module. This is so the proper time smearing can be
	// applied based on the 2006 beam test results which shot photons into the center
	// of the module.
	BCAL_UNATTENUATE_TO_CENTER = exp((DBCALGeometry::BCALFIBERLENGTH/2.0)/DBCALGeometry::ATTEN_LENGTH);
	
	// NOTE!!! This scheme needs to be changed if we want to get any
	// advantage from multi-threading! Otherwise, a mutex lock needs
	// to be added that encloses the entire BCAL smearing code.
	//
	// Create DBCALReadoutChannel objects for all readout channels.
	// It's easiest to do this for every cell regardless if summing is
	// turned on. Since the values are stored in a map, only one
	// DBCALReadoutChannel object will exist for every readout channel
	// in the end. That leaves the map such that it is fast to access
	// during event time since it is indexed by DBCALGeometry::fADCId(). 
	for(int i = 1; i<=48;i++){

		// inner
		for(int j = 1; j<=6;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j, k);
				int fADCId = DBCALGeometry::fADCId( i, j, k);
				int fADC_layer = DBCALGeometry::fADC_layer(cellId);
				double thresh = BCAL_inner_thresh[fADC_layer-1];
				bcal_fADCs[fADCId] = DBCALReadoutChannel(Nsum_inner[fADC_layer-1], thresh, i, fADC_layer, DBCALGeometry::fADC_sector(cellId));
			}
		}

		// outer
		for(int j = 1; j<=4;j++){
			for(int k = 1; k<=4; k++){
				int cellId = DBCALGeometry::cellId( i, j+6, k);
				int fADCId = DBCALGeometry::fADCId( i, j+6, k);
				int fADC_layer = DBCALGeometry::fADC_layer(cellId);
				double thresh = BCAL_outer_thresh[fADC_layer-1-DBCALGeometry::NBCALLAYSIN];
				bcal_fADCs[fADCId] = DBCALReadoutChannel(Nsum_outer[fADC_layer-1-DBCALGeometry::NBCALLAYSIN], thresh, i, fADC_layer, DBCALGeometry::fADC_sector(cellId));
			}
		}
	}

	// Flag that we have been called
	BCAL_INITIALIZED = true;
}

//-----------
// CalculateThreshold
//-----------
void CalculateThresholds(void)
{
	/// Calculate the BCAL thresholds based on the parameters in 
	/// the calibDB in BCAL/bcal_parms.
	///
	/// To do this, the electronic signal distribution (in MeV) is
	/// calculated for a single SiPM. This is used to determine a 
	/// cumulative distribution for multiple SiPMs in the case of
	/// summing.
	///
	/// For the course segmentation scheme, 3 inner SiPMs are summed
	/// while 4 outer SiPMs are summed. Because of this, different
	/// thresholds are calculated for the inner and outer regions.
	///
	/// The average number of "hits" per event due to dark pulses
	/// is given by the calibration parameter BCAL_AVG_DARK_DIGI_VALS_PER_EVENT.
	/// This sets a limit on the BCAL contribution to the average
	/// event size due to dark pulses and the thresholds are set 
	/// based on that.
	///
	/// A choice must be made on how to distribute the allowed,
	/// average number of hits between the inner and outer regions
	/// because they sum different numbers of SiPMs. This means
	/// they have different contributions to the signal spectrum
	/// due to the dark hits. I.e. for a given energy, there will
	/// be more dark pulses when summing 4 SiPMs together than
	/// when summing 3.
	/// 
	/// The method implemented here is to set the thresholds such that
	/// the fraction of SimPMs that are over threshold in a given event
	/// is the same for the inner SiPMs as for the outer.
	
	// Number of readout channels total
	DBCALGeometry g;
	
	int Nchan_inner=2*g.NBCALMODS*(4/g.NSUMSECSIN)*g.NBCALLAYSIN;
	int Nchan_outer=2*g.NBCALMODS*(4/g.NSUMSECSOUT)*g.NBCALLAYSOUT;
	
	// The inner BCAL region will have both an fADC and a TDC with
	// the fADC supplying both an amplitude and a time. The outer
	// region will have only an fADC so only 2 values will come
	// from those hits.
	//
	// To find the fraction of channels above threshold for a given
	// BCAL_AVG_DARK_DIGI_VALS_PER_EVENT, we use the fact that the
	// fraction of inner channels firing is equal to the fraction
	// of outer channels firing and that the sum of digitized values
	// is given by BCAL_AVG_DARK_DIGI_VALS_PER_EVENT:
	//
	// Ninner/Nchan_inner = Nouter/Nchan_outer
	//
	// BCAL_AVG_DARK_DIGI_VALS_PER_EVENT = 3*Ninner + 2*Nouter
	//
	// Solving for Nouter...
	//
	// Nouter = BCAL_AVG_DARK_DIGI_VALS_PER_EVENT/(3*Nchan_inner/Nchan_outer + 2)
	//
	// The desired fraction of readout channels firing is then:
	//
	// fraction_above_threshold = Nouter/Nchan_outer

	double Nouter = BCAL_AVG_DARK_DIGI_VALS_PER_EVENT/(3.0*(double)Nchan_inner/(double)Nchan_outer + 2.0);
	double fraction_above_threshold = Nouter/(double)Nchan_outer; // also equals Ninner/Nchan_inner

	// Some quantities are for single SiPM tiles, and some for 4x4 arrays
	double Ntiles = 16.0;

	// MeV per PhotoElectron in SiPM
	BCAL_mevPerPE = 1.0/( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE * BCAL_SAMPLING_FRACT );

	// Single PE signal and pedestal widths (values in calibDB are for single tile)
	double sigma_singlePE = BCAL_SIGMA_SIG_RELATIVE * BCAL_mevPerPE * k_MeV;
	double sigma_ped = sqrt(Ntiles) * BCAL_SIGMA_PED_RELATIVE * BCAL_mevPerPE * k_MeV;
	double sigma_gain_variation = BCAL_SIPM_GAIN_VARIATION * BCAL_mevPerPE * k_MeV;
	
	BCAL_sigma_singlePE_sq = sigma_singlePE*sigma_singlePE + sigma_gain_variation*sigma_gain_variation; // gain variation contributes in same way as single PE width
	BCAL_sigma_ped_sq = sigma_ped*sigma_ped;

	// Mean number of dark pulses from one SiPM inside integration window
	double mean_dark_pulses_array = BCAL_DARKRATE_GHZ * BCAL_INTWINDOW_NS; // BCAL_DARKRATE_GHZ is for whole 4x4 array


	// The following is code originally written by Dan Bennet at IU (I think).
	// The gist of it is that it uses the recursion relation for Poisson
	// distributions to build up a probability distribution for the number of
	// dark pulse photoelectrons per event (i.e. BCAL_INTWINDOW_NS time window).
	// This includes the cross-talk contribution, but not the cross-talk due
	// to cross talk, etc...
	// The probability of getting N photoelectrons from a SiPM array in a 
	// random BCAL_INTWINDOW_NS time window will be given by pdf[N].

	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><>
	// now we need to find n such that 
	// sum_i=1^n P(n) > ( 1 - maxOccupancy )
	// P(n) is the probability to have n pulses
	//
	// P(n) is really a convolution since we have:
	// P(n) = P( n_d + n_x ) = 
	//    P( n_d; nAvg ) * P( n_x; n_d * x )
	// 
	// n_d number of dark pulses
	// x is the cross talk rate

	// numerically build int_0^n P(n)

	// in practice the cutoff is going to be < 100 PE
	// we can build an accurate pdf up to that
	double pdf[100];
	for( int i = 0; i < 100; ++i ){ pdf[i] = 0; }

	// probability for zero is easy:
	double darkTerm = exp( -mean_dark_pulses_array );
	pdf[0] = darkTerm;

	for( int n_d = 1; n_d < 100; ++n_d ){

		darkTerm *= ( mean_dark_pulses_array / n_d );

		double xTalkAvg = n_d * BCAL_XTALK_FRACT;

		// probability for zero x-talk pulses
		double xTerm = exp( -xTalkAvg );
		pdf[n_d] += ( xTerm * darkTerm );

		// now include probability for additional
		// cross talk pulses
		for( int n_x = 1; n_x + n_d < 100; ++n_x ){

		  xTerm *= ( xTalkAvg / n_x );

		  pdf[n_d+n_x] += ( xTerm * darkTerm );
		}
	}
	//<><><><><><><><><><><><><><><><><><><><><><><><><><><><>

	// At this point we add in electronic noise using the sigma_singlePE
	// and sigma_ped calculated above. This is done by creating a histogram
	// and filling it with contributions from Npe=0 to Npe=99

	// we don't need to delete this because because ROOT owns it, I think?
	TH1D *bcal_dark_signal = new TH1D("bcal_dark_signal", "fADC signal", 1001, -0.05, 100.05);
	bcal_dark_signal->SetXTitle("fADC (MeV)");
	for(int bin=1; bin<=bcal_dark_signal->GetNbinsX(); bin++){

		double x = bcal_dark_signal->GetBinCenter(bin);
		
		
		// Loop over Npe
		double prob = 0.0;
		for(int Npe=0; Npe<100; Npe++){
			double sigma_tot = sqrt(BCAL_sigma_ped_sq + (double)Npe*BCAL_sigma_singlePE_sq)*1000.0; // in MeV
			double E = (double)Npe * BCAL_mevPerPE; // in MeV
			
			prob += TMath::Gaus(x, E , sigma_tot)*pdf[Npe];
		}

		bcal_dark_signal->SetBinContent(bin, prob);
	}

	// The histogram bcal_dark_signal now contains the probability 
	// distribution for the electronic response of a single SiPM array
	// in units of MeV. We wish to convolute that with itself
	// multiple times to get the PDF for multiple SiPMs added together
	// (in the case of summing).

	for (int i=0;i<DBCALGeometry::NBCALLAYSIN;i++) {
	  Nsum_inner[i] = g.NSUMLAYSIN[i]*g.NSUMSECSIN;

	  TH1D *bcal_dark_signal_sum_inner = CombineSiPM_PDFs("bcal_dark_signal_sum_inner", bcal_dark_signal, Nsum_inner[i]);

	  // Integrate the combined PDFs
	  TH1D *ibcal_dark_signal_sum_inner = (TH1D*)bcal_dark_signal_sum_inner->Clone("ibcal_dark_signal_sum_inner");

	  for(int bin=2; bin<=ibcal_dark_signal_sum_inner->GetNbinsX(); bin++){
	    double sum = ibcal_dark_signal_sum_inner->GetBinContent(bin-1) + ibcal_dark_signal_sum_inner->GetBinContent(bin);
	    ibcal_dark_signal_sum_inner->SetBinContent(bin, sum);
	  }

	  // Normalize
	  int last_bin = ibcal_dark_signal_sum_inner->GetNbinsX();
	  ibcal_dark_signal_sum_inner->Scale(1.0/ibcal_dark_signal_sum_inner->GetBinContent(last_bin));

	  // Finally, search the integrated, normalized histograms for the bin
	  // where the integral fraction above it is less than fraction_above_threshold
	  for(int bin=1 ; bin<=last_bin; bin++){
	    if(ibcal_dark_signal_sum_inner->GetBinContent(bin) < (1.0 - fraction_above_threshold)){
	      BCAL_inner_thresh[i] = ibcal_dark_signal_sum_inner->GetBinCenter(bin) * k_MeV;
	    }
	  }
    }

	for (int i=0;i<DBCALGeometry::NBCALLAYSOUT;i++) {
	  Nsum_outer[i] = g.NSUMLAYSOUT[i]*g.NSUMSECSOUT;

	  TH1D *bcal_dark_signal_sum_outer = CombineSiPM_PDFs("bcal_dark_signal_sum_outer", bcal_dark_signal, Nsum_outer[i]);

	  // Integrate the combined PDFs
	  TH1D *ibcal_dark_signal_sum_outer = (TH1D*)bcal_dark_signal_sum_outer->Clone("ibcal_dark_signal_sum_outer");

	  for(int bin=2; bin<=ibcal_dark_signal_sum_outer->GetNbinsX(); bin++){
	    double sum = ibcal_dark_signal_sum_outer->GetBinContent(bin-1) + ibcal_dark_signal_sum_outer->GetBinContent(bin);
	    ibcal_dark_signal_sum_outer->SetBinContent(bin, sum);
	  }

	  // Normalize
	  int last_bin = ibcal_dark_signal_sum_outer->GetNbinsX();
	  ibcal_dark_signal_sum_outer->Scale(1.0/ibcal_dark_signal_sum_outer->GetBinContent(last_bin));

	  // Finally, search the integrated, normalized histograms for the bin
	  // where the integral fraction above it is less than fraction_above_threshold
	  for(int bin=1 ; bin<=last_bin; bin++){
	    if(ibcal_dark_signal_sum_outer->GetBinContent(bin) < (1.0 - fraction_above_threshold)){
	      BCAL_outer_thresh[i] = ibcal_dark_signal_sum_outer->GetBinCenter(bin) * k_MeV;
	    }
	  }
    }
}


//-----------
// CombineSiPM_PDFs
//-----------
TH1D* CombineSiPM_PDFs(const char *hname, TH1D *single, int NSiPM)
{
	// It is quicker (and easier) to do this in stages by calculating the
	// probability of signals adding to a given energy and then 
	// subsequently calculating the probability of adding one more.
	
	TH1D *summed = (TH1D*)single->Clone(hname);
	
	for(int i=2; i<=NSiPM; i++){
	
		// Make copy of summed histogram since we'll be overwritting it
		TH1D *tmp = (TH1D*)summed->Clone("tmp");
		
		// Reset summed histo (for safety)
		summed->Reset();
		
		// Loop over all combinations of bins in the previous iteration's
		// summed probability histo and the single SiPM probability
		// histo. Fill in the new summed probability histo.
		for(int ibin=1; ibin<=tmp->GetNbinsX(); ibin++){

			double NPE1 = tmp->GetBinCenter(ibin);
			double prob1 = tmp->GetBinContent(ibin);

			for(int jbin=1; jbin<=single->GetNbinsX(); jbin++){

				double NPE2 = single->GetBinCenter(jbin);
				double prob2 = single->GetBinContent(jbin);
				
				// Find the bin in the summed histogram corresponding to
				// the sum of the energies from the bin in the previous
				// iteration's summed prob. histo and the single SiPM
				// prob. histo. If the energy sum is in range, then set
				// the bin content of the current iteration's summed prob.
				// histo.
				int kbin = summed->FindBin(NPE1 + NPE2);
				if(kbin>=1 && kbin<=summed->GetNbinsX()){
					summed->Fill(NPE1 + NPE2, prob1*prob2);
				}
			}
		}
        
        delete tmp;
	}
	
	return summed;
}


#endif
