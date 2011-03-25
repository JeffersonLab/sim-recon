// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <FCAL/DFCALGeometry.h>
#include <CCAL/DCCALGeometry.h>
#include <BCAL/DBCALGeometry.h>

#include <math.h>
#include "units.h"
#include "HDDM/hddm_s.h"
#include <TF1.h>
#include <TH2.h>
#include <TRandom3.h>

float RANDOM_MAX = (float)(0x7FFFFFFF);
#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


extern vector<vector<float> >fdc_smear_parms; 
extern TF1 *fdc_smear_function;
extern TH2F *fdc_drift_time_smear_hist;
extern TH2F *fdc_drift_dist_smear_hist;
extern TH2F *fdc_drift_time;

void SmearCDC(s_HDDM_t *hddm_s);
void AddNoiseHitsCDC(s_HDDM_t *hddm_s);
void SmearFDC(s_HDDM_t *hddm_s);
void AddNoiseHitsFDC(s_HDDM_t *hddm_s);
void SmearFCAL(s_HDDM_t *hddm_s);
void SmearCCAL(s_HDDM_t *hddm_s);
void SmearBCAL(s_HDDM_t *hddm_s);
void SmearTOF(s_HDDM_t *hddm_s);
void SmearSTC(s_HDDM_t *hddm_s);
void SmearCherenkov(s_HDDM_t *hddm_s);
void InitCDCGeometry(void);
void InitFDCGeometry(void);

void bcalInit();
float bcalSamplingSmear( float E );
float bcalTimeSmear(float t, float e);
int getDarkHits();

pthread_mutex_t mutex_fdc_smear_function = PTHREAD_MUTEX_INITIALIZER;

bool CDC_GEOMETRY_INITIALIZED = false;
int CDC_MAX_RINGS=0;
vector<unsigned int> NCDC_STRAWS;
vector<double> CDC_RING_RADIUS;

DFCALGeometry *fcalGeom = NULL;
DCCALGeometry *ccalGeom = NULL;
bool FDC_GEOMETRY_INITIALIZED = false;
unsigned int NFDC_WIRES_PER_PLANE;
vector<double> FDC_LAYER_Z;
double FDC_RATE_COEFFICIENT;

double SampleGaussian(double sigma);
int SamplePoisson(float lambda);
double SampleRange(double x1, double x2);

// Do we or do we not add noise hits
extern bool ADD_NOISE;

// Do we or do we not smear real hits
extern bool SMEAR_HITS;

// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// If the following flag is true, then include the drift-distance
// dependency on the error in the CDC position. Otherwise, use a
// flat distribution given by the CDC_TDRIFT_SIGMA below.
extern bool CDC_USE_PARAMETERIZED_SIGMA;

// The error on the drift time in the CDC. The drift times
// for the actual CDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double CDC_TDRIFT_SIGMA;

// The time window for which CDC hits are accumulated.
// This is used to determine the number of background
// hits in the CDC for a given event.
extern double CDC_TIME_WINDOW;

// The error in the energy deposition measurement in the CDC due to pedestal noise
extern double CDC_PEDESTAL_SIGMA;
 
// If the following flag is true, then include the drift-distance
// dependency on the error in the FDC position. Otherwise, use a
// flat distribution given by the FDC_TDRIFT_SIGMA below.
extern bool FDC_USE_PARAMETERIZED_SIGMA;

// The error on the drift time in the FDC. The drift times
// for the actual FDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double FDC_TDRIFT_SIGMA;

// The error in the distance along the wire as measured by
// the cathodes. This should NOT include the Lorentz
// effect which is already included in hdgeant. It
// should include any fluctuations due to ion trail density
// etc.
extern double FDC_CATHODE_SIGMA;

// The FDC pedestal noise is used to smear the cathode ADC
// values such that the position along the wire has the resolution
// specified by FDC_CATHODE_SIGMA.
extern double FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// If energy loss was turned off in the FDC then the pedestal
// noise will not be scaled properly to give the nominal 200 micron
// resolution along the wires. This flag is used to indicated
// the magic scale factor should be applied to FDC_PED_NOISE
// when it is calculated below to get the correct resolution.
extern bool FDC_ELOSS_OFF;

// Time window for acceptance of FDC hits
extern double FDC_TIME_WINDOW;

// Fraction of FDC hits to randomly drop (0=drop nothing 1=drop everything)
extern double FDC_HIT_DROP_FRACTION;

// Photon-statistics factor for smearing hit energy (from Criss's MC)
// (copied from DFCALMCResponse_factory.cc 7/2/2009 DL)
extern double FCAL_PHOT_STAT_COEF;

// Single block energy threshold (applied after smearing)
extern double FCAL_BLOCK_THRESHOLD;

// Photon-statistics factor for smearing hit energy for CompCal
// (This is just a rough estimate 11/30/2010 DL)
double CCAL_PHOT_STAT_COEF = 0.035/2.0;

// Single block energy threshold (applied after smearing)
// (This is just a rough estimate 11/30/2010 DL)
double CCAL_BLOCK_THRESHOLD = 20.0*k_MeV;

// setup response parameters
extern float BCAL_DARKRATE_GHZ        ;// 0.041;
extern float BCAL_XTALK_FRACT         ;//0.03;
extern float BCAL_INTWINDOW_NS        ;//100;
extern float BCAL_DEVICEPDE           ;//0.12;
extern float BCAL_SAMPLING_FRACT      ;//0.15;
extern float BCAL_MAXOCCUPANCY_FRACT  ;//0.05;

// GX-doc 1069, Table 1 -- try to extract back to
// photons per side per MeV in fiber
// 4.6 / PDE / attentuation  (meaurements performed in center)
// 75 = 4.6  / 0.12 / exp( -200 / 300 )
extern float BCAL_PHOTONSPERSIDEPERMEV_INFIBER; // 75

// set the sampling smearing coefficients:
// (from GlueX-doc 827 v3 Figure 13 )
extern float BCAL_SAMPLINGCOEFA ; //= 0.042;
extern float BCAL_SAMPLINGCOEFB ; //= 0.013;

// time smearing comes from beam test time difference resolution with
// beam incident on the center of the module
extern float BCAL_TIMEDIFFCOEFA ; //= 0.07 * sqrt( 2 );

// no floor term, but leave the option here:
extern float BCAL_TIMEDIFFCOEFB ; //= 0.0 * sqrt( 2 );
    
// set this low for now -- needs more thought later
extern float BCAL_CELLOUTERTHRESHOLD ; // = 1 * k_MeV;

// calculated later in based on number of photons and ph threshold
extern float Bcal_CellInnerThreshold ;

// Forward TOF resolution
extern double TOF_SIGMA;
extern double TOF_PHOTONS_PERMEV;

// Start counter resolution
extern double START_SIGMA ;
extern double START_PHOTONS_PERMEV;

// Polynomial interpolation on a grid.
// Adapted from Numerical Recipes in C (2nd Edition), pp. 121-122.
void polint(float *xa, float *ya,int n,float x, float *y,float *dy){
  int i,m,ns=0;
  float den,dif,dift,ho,hp,w;

  float *c=(float *)calloc(n,sizeof(float));
  float *d=(float *)calloc(n,sizeof(float));

  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++){
    if ((dift=fabs(x-xa[i]))<dif){
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];

  for (m=1;m<n;m++){
    for (i=1;i<=n-m;i++){
      ho=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i+1-1]-d[i-1];
      if ((den=ho-hp)==0.0) return;
      
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
      
    }
    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
}

//-----------
// Smear
//-----------
void Smear(s_HDDM_t *hddm_s)
{
	if(SMEAR_HITS)SmearCDC(hddm_s);
	if(ADD_NOISE)AddNoiseHitsCDC(hddm_s);
	if(SMEAR_HITS)SmearFDC(hddm_s);
	if(ADD_NOISE)AddNoiseHitsFDC(hddm_s);
	if(SMEAR_HITS)SmearFCAL(hddm_s);
	if(SMEAR_HITS)SmearCCAL(hddm_s);
	if(SMEAR_BCAL)SmearBCAL(hddm_s);
	if(SMEAR_HITS)SmearTOF(hddm_s);
	if(SMEAR_HITS)SmearSTC(hddm_s);
	if(SMEAR_HITS)SmearCherenkov(hddm_s);
}

//-----------
// SmearCDC
//-----------
void SmearCDC(s_HDDM_t *hddm_s)
{
	/// Smear the drift times of all CDC hits.
	/// This will add cdcStrawHit objects generated by smearing values in the
	/// cdcStrawTruthHit objects that hdgeant outputs. Any existing cdcStrawHit
	/// objects will be replaced.

	// Acquire the pointer to the physics events
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents)return;
       
	for (unsigned int m=0; m < allEvents->mult; m++) {
		// Acquire the pointer to the overall hits section of the data
		s_HitView_t *hits = allEvents->in[m].hitView;
		
		if (hits == HDDM_NULL)return;
		if (hits->centralDC == HDDM_NULL)return;
		if (hits->centralDC->cdcStraws == HDDM_NULL)return;
		for(unsigned int k=0; k<hits->centralDC->cdcStraws->mult; k++){
			s_CdcStraw_t *cdcstraw = &hits->centralDC->cdcStraws->in[k];
			
			// If the event already contains a s_CdcStrawHits_t object then free it.
			if(cdcstraw->cdcStrawHits!=HDDM_NULL){
				FREE(cdcstraw->cdcStrawHits);
				cerr<<"Values exist!"<<endl;
			}
			
			// Create new s_CdcStrawHits_t
			cdcstraw->cdcStrawHits = make_s_CdcStrawHits(cdcstraw->cdcStrawTruthHits->mult);
			cdcstraw->cdcStrawHits->mult = cdcstraw->cdcStrawTruthHits->mult;
			
			for(unsigned int j=0; j<cdcstraw->cdcStrawTruthHits->mult; j++){
				s_CdcStrawHit_t *strawhit = &cdcstraw->cdcStrawHits->in[j];
				s_CdcStrawTruthHit_t *strawtruthhit = &cdcstraw->cdcStrawTruthHits->in[j];

				double sigma_t = CDC_TDRIFT_SIGMA;
				if(CDC_USE_PARAMETERIZED_SIGMA){
					// Convert drift time back to drift distance assuming standard 55 um/ns
					//double drift_d = strawtruthhit->t*55.0E-3; // use mm since that's what the error function was paramaterized
				  double drift_d=strawtruthhit->d*10.;
					
					// The following is from a fit to Yves' numbers circa Aug 2009. The values fit were
					// resolution (microns) vs. drift distance (mm).
					// par[8] = {699.875, -559.056, 149.391, 25.6929, -22.0238, 4.75091, -0.452373, 0.0163858};
					double x = drift_d;
					//double sigma_d = (699.875) + x*((-559.056) + x*((149.391) + x*((25.6929) + x*((-22.0238) + x*((4.75091) + x*((-0.452373) + x*((0.0163858))))))));
					double sigma_d = 108.55 + 7.62391*x + 556.176*exp(-(1.12566)*pow(x,1.29645));
					sigma_t = sigma_d/55.0; // remember that sigma_d is already in microns here!
					sigma_t *= 1.0E-9; // convert sigma_t to seconds
				}

				// Smear out the CDC drift time using the specified sigma.
				// This should include both timing resolution and ion trail
				// density effects.
				double delta_t = SampleGaussian(sigma_t)*1.0E9; // delta_t is in ns
				strawhit->t = strawtruthhit->t + delta_t;
				
				// Smear energy deposition
				strawhit->dE = strawtruthhit->dE+SampleGaussian(CDC_PEDESTAL_SIGMA);
		
				// to be consistent initialize the value d=DOCA to zero
				strawhit->d = 0.;
				strawhit->itrack = strawtruthhit->itrack;
				strawhit->ptype = strawtruthhit->ptype;
				// If the time is negative, reject this smear and try again
				//if(strawhit->t<0)j--;
			}
		}
	}
}

//-----------
// AddNoiseHitsCDC
//-----------
void AddNoiseHitsCDC(s_HDDM_t *hddm_s)
{
	if(!CDC_GEOMETRY_INITIALIZED)InitCDCGeometry();
	
	// Calculate the number of noise hits for each straw and store
	// them in a sparse map. We must do it this way since we have to know
	// the total number of CdcStraw_t structures to allocate in our
	// call to make_s_CdcStraws.
	//
	// The straw rates are obtained using a parameterization done
	// to calculate the event size for the August 29, 2007 online
	// meeting. This parameterization is almost already obsolete.
	// 10/12/2007 D. L.
	vector<int> Nstraw_hits;
	vector<int> straw_number;
	vector<int> ring_number;
	int Nnoise_straws = 0;
	int Nnoise_hits = 0;
	for(unsigned int ring=1; ring<=NCDC_STRAWS.size(); ring++){
		double p[2] = {10.4705, -0.103046};
		double r_prime = (double)(ring+3);
		double N = exp(p[0] + r_prime*p[1]);
		N *= CDC_TIME_WINDOW;
		for(unsigned int straw=1; straw<=NCDC_STRAWS[ring-1]; straw++){
			// Indivdual straw rates should be way less than 1/event so
			// we just use the rate as a probablity.
			double Nhits = SampleRange(0.0, 1.0)<N ? 1.0:0.0;
			if(Nhits<1.0)continue;
			int iNhits = (int)floor(Nhits);
			Nstraw_hits.push_back(iNhits);
			straw_number.push_back(straw);
			ring_number.push_back(ring);
			Nnoise_straws++;
			Nnoise_hits+=iNhits;
		}
	}

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL)continue;
		
		// If no CDC hits were produced by HDGeant, then we need to add
		// the branches needed to the HDDM tree
		if(hits->centralDC == HDDM_NULL){
			hits->centralDC = make_s_CentralDC();
			hits->centralDC->cdcStraws = (s_CdcStraws_t*)HDDM_NULL;
			hits->centralDC->cdcTruthPoints = (s_CdcTruthPoints_t*)HDDM_NULL;
		}

		if(hits->centralDC->cdcStraws == HDDM_NULL){
			hits->centralDC->cdcStraws = make_s_CdcStraws(0);
			hits->centralDC->cdcStraws->mult=0;
		}
		
		// Get existing hits
		s_CdcStraws_t *old_cdcstraws = hits->centralDC->cdcStraws;
		unsigned int Nold = old_cdcstraws->mult;

		// Create CdcStraws structure that has enough slots for
		// both the real and noise hits.
		s_CdcStraws_t* cdcstraws = make_s_CdcStraws((unsigned int)Nnoise_straws + Nold);

		// Add real hits back in first
		cdcstraws->mult = 0;
		for(unsigned int j=0; j<Nold; j++){
			cdcstraws->in[cdcstraws->mult++] = old_cdcstraws->in[j];
			
			// We need to transfer ownership of the hits to the new cdcstraws
			// branch so they don't get deleted when old_cdcstraws is freed.
			s_CdcStraw_t *cdcstraw = &old_cdcstraws->in[j];
			cdcstraw->cdcStrawHits = (s_CdcStrawHits_t *)HDDM_NULL;
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_cdcstraws);
		hits->centralDC->cdcStraws = cdcstraws;
		
		// Loop over straws with noise hits
		for(unsigned int j=0; j<Nstraw_hits.size(); j++){
			s_CdcStraw_t *cdcstraw = &cdcstraws->in[cdcstraws->mult++];
			s_CdcStrawHits_t *strawhits = make_s_CdcStrawHits(Nstraw_hits[j]);
			cdcstraw->cdcStrawHits = strawhits;
			cdcstraw->ring = ring_number[j];
			cdcstraw->straw = straw_number[j];

			strawhits->mult = 0;
			for(int k=0; k<Nstraw_hits[j]; k++){
				s_CdcStrawHit_t *strawhit = &strawhits->in[strawhits->mult++];
				strawhit->dE = 1.0;
				strawhit->t = SampleRange(-CDC_TIME_WINDOW/2.0, +CDC_TIME_WINDOW/2.0)*1.e9;
				strawhit->d = 0.; // be consistent to initialize d=DOCA to zero
			}
		}
	}
}

//-----------
// SmearFDC
//-----------
void SmearFDC(s_HDDM_t *hddm_s)
{
	// Calculate ped noise level based on position resolution
	FDC_PED_NOISE=-0.004594+0.008711*FDC_CATHODE_SIGMA+0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
	if(FDC_ELOSS_OFF)FDC_PED_NOISE*=7.0; // empirical  4/29/2009 DL

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardDC == HDDM_NULL ||
			hits->forwardDC->fdcChambers == HDDM_NULL)continue;

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		s_FdcChamber_t *fdcChamber = fdcChambers->in;
		for(unsigned int j=0; j<fdcChambers->mult; j++, fdcChamber++){
			
			// Add pedestal noise to strip charge data
			s_FdcCathodeStrips_t *strips= fdcChamber->fdcCathodeStrips;
			if (strips!=HDDM_NULL){
			  s_FdcCathodeStrip_t *strip=strips->in;
			  for (unsigned int k=0;k<strips->mult;k++,strip++){
				
				// If a s_FdcCathodeHits_t object already exists delete it
				if(strip->fdcCathodeHits!=HDDM_NULL){
					FREE(strip->fdcCathodeHits);
					strip->fdcCathodeHits = (s_FdcCathodeHits_t*)HDDM_NULL;
				}
			  
			    s_FdcCathodeTruthHits_t *truthhits=strip->fdcCathodeTruthHits;
			    if (truthhits==HDDM_NULL)continue;
				 
				 // Allocate s_FdcCathodeHits_t objects corresponding to this s_FdcCathodeTruthHits_t
				 s_FdcCathodeHits_t *hits = strip->fdcCathodeHits = make_s_FdcCathodeHits(truthhits->mult);
				 hits->mult = truthhits->mult;
				 
			    s_FdcCathodeHit_t *hit=hits->in;
			    s_FdcCathodeTruthHit_t *truthhit=truthhits->in;
			    for (unsigned int s=0;s<hits->mult;s++,hit++,truthhit++){
					if(SampleRange(0.0, 1.0)<=FDC_HIT_DROP_FRACTION){
						hit->q = 0.0;
						hit->t = 1.0E6;
					}else{
						hit->q = truthhit->q + SampleGaussian(FDC_PED_NOISE);
						hit->t = truthhit->t;
					}
			    }
			  }
			}

			// Add drift time varation to the anode data 
			s_FdcAnodeWires_t *wires=fdcChamber->fdcAnodeWires;
			
			if (wires!=HDDM_NULL){
			  s_FdcAnodeWire_t *wire=wires->in;
			  for (unsigned int k=0;k<wires->mult;k++,wire++){

				// If a s_FdcAnodeHits_t object exists already delete it
				if(wire->fdcAnodeHits!=HDDM_NULL){
					FREE(wire->fdcAnodeHits);
					wire->fdcAnodeHits = (s_FdcAnodeHits_t*)HDDM_NULL;
				}

			    s_FdcAnodeTruthHits_t *truthhits=wire->fdcAnodeTruthHits;
			    if (truthhits==HDDM_NULL)continue;
				 
				 // Allocate s_FdcAnodeHits_t object corresponding to this s_FdcAnodeTruthHits_t
				 s_FdcAnodeHits_t *hits = wire->fdcAnodeHits = make_s_FdcAnodeHits(truthhits->mult);
				 hits->mult = truthhits->mult;
				 
			    s_FdcAnodeHit_t *hit=hits->in;
			    s_FdcAnodeTruthHit_t *truthhit=truthhits->in;
			    for (unsigned int s=0;s<hits->mult;s++, hit++, truthhit++){
			      if (FDC_USE_PARAMETERIZED_SIGMA==false){
				hit->t = truthhit->t + SampleGaussian(FDC_TDRIFT_SIGMA)*1.0E9;
			      }
			      else{		
				// The scale of the smearing in time depends 
				// on the drift distance.  We use the root TF1
				// random generator using functions representing
				// the degree of smearing in various bins in x
	  
				int ind=(truthhit->d<0.5?int(floor(truthhit->d/0.02+0.5)):25);
				double dt=0.;
				/*
				dt=SampleGaussian((5.96642*exp(-54.2745*truthhit->d)
							  +1.64413e-11*exp(56.1646*truthhit->d)
							  +2.48224));
				*/
				float xarray[26]={0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,
						  0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,
						  0.48,0.50};
	     
				float my_parms[9];
				//printf("imin %d imax %d\n",imin,imin+3);
					     
				for (int m=0;m<9;m++){
				  float dpar;
				  float parlist[26];
				  int  num=3;
				  if (ind>=23){
				    num=2;
				    if (ind==25) ind=24;
				  }
				  for (int p=0;p<num;p++){
				    parlist[p]=fdc_smear_parms[ind+p][m];
				    // printf("%f\n",parlist[p]);
				  }
				  polint(&xarray[ind],parlist,num,truthhit->d,&my_parms[m],&dpar);
				  //printf("d %f par %f\n",truthhit->d,my_parms[m]);
			
				}

				pthread_mutex_lock(&mutex_fdc_smear_function);
				for (unsigned int p=0;p<9;p++){
				  // printf("rel %f\n",truthhit->d-0.02*ind);
				  if ((p+2)%3 && my_parms[p]<0.) my_parms[p]*=-1.;
				  fdc_smear_function->SetParameter(p,my_parms[p]);
				}
				// The calib database contains the smearing in cm, so one needs to 
				// convert into nanoseconds...
				double dx=fdc_smear_function->GetRandom();
				pthread_mutex_unlock(&mutex_fdc_smear_function);
				dt=dx/0.0055;
				fdc_drift_time_smear_hist->Fill(truthhit->d,dt);
				double t_temp=90.99*tanh((truthhit->d-0.2915)/0.1626);

				fdc_drift_dist_smear_hist->Fill(truthhit->d,dx);
				fdc_drift_time->Fill(t_temp+86.67+dx*(1.-t_temp*t_temp*0.01097*0.01097)
						     /(0.01097*0.08211),truthhit->d);
				
				hit->t=truthhit->t+dt;
			      }
			      hit->dE = truthhit->dE;
			      hit->d  = 0.; // initialize d=DOCA to zero for consistency.
			    }
			  }
			}
		}
	}
}

//-----------
// AddNoiseHitsFDC
//-----------
void AddNoiseHitsFDC(s_HDDM_t *hddm_s)
{
	if(!FDC_GEOMETRY_INITIALIZED)InitFDCGeometry();
	
	// Calculate the number of noise hits for each FDC wire and store
	// them in a sparse map. We must do it this way since we have to know
	// the total number of s_FdcAnodeWire_t structures to allocate in our
	// call to make_s_FdcAnodeWires.
	//
	// We do this using the individual wire rates to calculate the probability
	// of the wire firing for a single event. For the FDC, we calculate the
	// wire rates as a function of both wire number (distance from beam line)
	// and layer (position in z). We want a roughly 1/r distribution in the
	// radial direction and a roughly exponential rise in rate in the
	// +z direction.
	//
	// The wire rates are obtained using a parameterization done
	// to calculate the event size for the August 29, 2007 online
	// meeting. This parameterization is almost already obsolete.
	// In rough terms, the layer rate (integrated over all wires)
	// is about 1 MHz. For a 24 layer chamber with a 1us time window,
	// we should have approximately 24 background hits per event.
	// 11/9/2007 D. L.
	vector<int> Nwire_hits;
	vector<int> wire_number;
	vector<int> layer_number;
	int Nnoise_wires = 0;
	int Nnoise_hits = 0;
	for(unsigned int layer=1; layer<=FDC_LAYER_Z.size(); layer++){
		double No = FDC_RATE_COEFFICIENT*exp((double)layer*log(4.0)/24.0);
		for(unsigned int wire=1; wire<=96; wire++){
			double rwire = fabs(96.0/2.0 - (double)wire);
			double N = No*log((rwire+0.5)/(rwire-0.5));

			// Indivdual wire rates should be way less than 1/event so
			// we just use the rate as a probablity.
			double Nhits = SampleRange(0.0, 1.0)<N ? 1.0:0.0;
			if(Nhits<1.0)continue;
			int iNhits = (int)floor(Nhits);
			Nwire_hits.push_back(iNhits);
			wire_number.push_back(wire);
			layer_number.push_back(layer);
			Nnoise_wires++;
			Nnoise_hits+=iNhits;
		}
	}

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL)continue;
		
		// If no FDC hits were produced by HDGeant, then we need to add
		// the branches needed to the HDDM tree
		if(hits->forwardDC == HDDM_NULL){
			hits->forwardDC = make_s_ForwardDC();
			hits->forwardDC->fdcChambers = (s_FdcChambers_t*)HDDM_NULL;
		}

		if(hits->forwardDC->fdcChambers == HDDM_NULL){
			hits->forwardDC->fdcChambers = make_s_FdcChambers(0);
			hits->forwardDC->fdcChambers->mult=0;
		}
		
		// Get existing hits
		s_FdcChambers_t *old_fdcchambers = hits->forwardDC->fdcChambers;
		unsigned int Nold = old_fdcchambers->mult;

		// If we were doing this "right" we'd conglomerate all of the noise
		// hits from the same chamber into the same s_FdcChamber_t structure.
		// That's a pain in the butt and really may only save a tiny bit of disk
		// space so we just add each noise hit back in as another chamber
		// structure.
		

		// Create FdcChambers structure that has enough slots for
		// both the real and noise hits.
		s_FdcChambers_t* fdcchambers = make_s_FdcChambers(Nwire_hits.size() + Nold);

		// Add real hits back in first
		fdcchambers->mult = 0;
		for(unsigned int j=0; j<Nold; j++){
			fdcchambers->in[fdcchambers->mult++] = old_fdcchambers->in[j];
			
			// We need to transfer ownership of the hits to the new fdcchambers
			// branch so they don't get deleted when old_fdcchambers is freed.
			s_FdcChamber_t *fdcchamber = &old_fdcchambers->in[j];
			fdcchamber->fdcAnodeWires = (s_FdcAnodeWires_t *)HDDM_NULL;
			fdcchamber->fdcCathodeStrips = (s_FdcCathodeStrips_t *)HDDM_NULL;
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_fdcchambers);
		hits->forwardDC->fdcChambers = fdcchambers;

		// Loop over wires with noise hits
		for(unsigned int j=0; j<Nwire_hits.size(); j++){
			s_FdcChamber_t *fdcchamber = &fdcchambers->in[fdcchambers->mult++];
			
			// Create structure for anode wires
			s_FdcAnodeWires_t *fdcAnodeWires = make_s_FdcAnodeWires(Nwire_hits[j]);
			fdcchamber->fdcAnodeWires = fdcAnodeWires;

			fdcAnodeWires->mult = 0;

			for(int k=0; k<Nwire_hits[j]; k++){
				// Get pointer to anode wire structure
				s_FdcAnodeWire_t *fdcAnodeWire = &fdcAnodeWires->in[fdcAnodeWires->mult++];

				// Create anode hits structure
				s_FdcAnodeHits_t *fdcanodehits = make_s_FdcAnodeHits(1);
				fdcAnodeWire->fdcAnodeHits = fdcanodehits;
				
				// Get pointer to anode hit structure
				fdcanodehits->mult = 1;
				s_FdcAnodeHit_t *fdcanodehit = &fdcanodehits->in[0];
				
				fdcanodehit->dE = 0.1; // what should this be?
				fdcanodehit->t = SampleRange(-FDC_TIME_WINDOW/2., +FDC_TIME_WINDOW/2.)*1.e9;
				fdcanodehit->d = 0.; // d=DOCA initialize to avoid any NAN
				
				fdcAnodeWire->wire = wire_number[j];
				
				fdcchamber->layer = (layer_number[j]-1)%3 + 1;
				fdcchamber->module = (layer_number[j]-1)/3 + 1;
			}
		}
	}
}

//-----------
// SmearFCAL
//-----------
void SmearFCAL(s_HDDM_t *hddm_s)
{
	/// Smear the FCAL hits using the nominal resolution of the individual blocks.
	/// The way this works is a little funny and warrants a little explanation.
	/// The information coming from hdgeant is truth information indexed by 
	/// row and column, but containing energy deposited and time. The mcsmear
	/// program will copy the truth information from the FcalTruthBlock to the
	/// FcalBlock branch, smearing the values with the appropriate detector
	/// resolution.
	///
	/// To access the "truth" values in DANA, get the DFCALHit objects using the
	/// "TRUTH" tag.
	
	// The code below is perhaps slightly odd in that it simultaneously creates
	// and fills the mirror (truth) branch while smearing the values in the
	// nominal hit branch.
	
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	if(!fcalGeom)fcalGeom = new DFCALGeometry();

	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			 hits->forwardEMcal == HDDM_NULL ||
			 hits->forwardEMcal->fcalBlocks == HDDM_NULL)continue;
		
		s_FcalBlocks_t *blocks = hits->forwardEMcal->fcalBlocks;
		for(unsigned int j=0; j<blocks->mult; j++){
			s_FcalBlock_t *block = &blocks->in[j];

			// Create FCAL hits structures to put smeared data into
			if(block->fcalHits!=HDDM_NULL)free(block->fcalHits);
			block->fcalHits = make_s_FcalHits(block->fcalTruthHits->mult);
			block->fcalHits->mult = block->fcalTruthHits->mult;

			for(unsigned int k=0; k<block->fcalTruthHits->mult; k++){
				s_FcalTruthHit_t *fcaltruthhit = &block->fcalTruthHits->in[k];
				s_FcalHit_t *fcalhit = &block->fcalHits->in[k];

				// Copy info from truth stream before doing anything else
				fcalhit->E = fcaltruthhit->E;
				fcalhit->t = fcaltruthhit->t;

				// Simulation simulates a grid of blocks for simplicity. 
				// Do not bother smearing inactive blocks. They will be
				// discarded in DEventSourceHDDM.cc while being read in
				// anyway.
				if(!fcalGeom->isBlockActive( block->row, block->column ))continue;

				// Smear the energy and timing of the hit
				double sigma = FCAL_PHOT_STAT_COEF/sqrt(fcalhit->E) ;
				fcalhit->E *= 1.0 + SampleGaussian(sigma);
				fcalhit->t += SampleGaussian(200.0E-3); // smear by 200 ps fixed for now 7/2/2009 DL
				
				// Apply a single block threshold. If the (smeared) energy is below this,
				// then set the energy and time to zero. 
				if(fcalhit->E < FCAL_BLOCK_THRESHOLD){fcalhit->E = fcalhit->t = 0.0;}

			} // k  (fcalhits)
		} // j  (blocks)
	} // i  (physicsEvents)

}

//-----------
// SmearCCAL
//-----------
void SmearCCAL(s_HDDM_t *hddm_s)
{
	/// Smear the CCAL hits using the same procedure as the FCAL above.
	/// See those comments for details.

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	if(!ccalGeom)ccalGeom = new DCCALGeometry();

	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			 hits->ComptonEMcal == HDDM_NULL ||
			 hits->ComptonEMcal->ccalBlocks == HDDM_NULL)continue;
		
		s_CcalBlocks_t *blocks = hits->ComptonEMcal->ccalBlocks;
		for(unsigned int j=0; j<blocks->mult; j++){
			s_CcalBlock_t *block = &blocks->in[j];

			// Create FCAL hits structures to put smeared data into
			if(block->ccalHits!=HDDM_NULL)free(block->ccalHits);
			block->ccalHits = make_s_CcalHits(block->ccalTruthHits->mult);
			block->ccalHits->mult = block->ccalTruthHits->mult;

			for(unsigned int k=0; k<block->ccalTruthHits->mult; k++){
				s_CcalTruthHit_t *ccaltruthhit = &block->ccalTruthHits->in[k];
				s_CcalHit_t *ccalhit = &block->ccalHits->in[k];

				// Copy info from truth stream before doing anything else
				ccalhit->E = ccaltruthhit->E;
				ccalhit->t = ccaltruthhit->t;

				// Simulation simulates a grid of blocks for simplicity. 
				// Do not bother smearing inactive blocks. They will be
				// discarded in DEventSourceHDDM.cc while being read in
				// anyway.
				if(!ccalGeom->isBlockActive( block->row, block->column ))continue;

				// Smear the energy and timing of the hit
				double sigma = CCAL_PHOT_STAT_COEF/sqrt(ccalhit->E) ;
				ccalhit->E *= 1.0 + SampleGaussian(sigma);
				ccalhit->t += SampleGaussian(200.0E-3); // smear by 200 ps fixed for now 7/2/2009 DL
				
				// Apply a single block threshold. If the (smeared) energy is below this,
				// then set the energy and time to zero. 
				if(ccalhit->E < CCAL_BLOCK_THRESHOLD){ccalhit->E = ccalhit->t = 0.0;}

			} // k  (fcalhits)
		} // j  (blocks)
	} // i  (physicsEvents)

}

//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{
  DBCALGeometry *bcalGeom = new DBCALGeometry();
  
  map< int, pair< int, int > > darkHits;
  map< int, pair< int, int > >::iterator darkHitItr;
  
  
  for( int m = 1; m <= 48; ++m ){
    for( int l = 1; l <= 4; ++l ){
      for( int s = 1; s <= 10; ++s ){
	
	pair< int, int > nHits( getDarkHits(), getDarkHits() );
	
	darkHits[DBCALGeometry::cellId( m, l, s )] = nHits;
      }
    }
  }
  
  double mevPerPE = 1 / 
    ( BCAL_PHOTONSPERSIDEPERMEV_INFIBER * BCAL_DEVICEPDE *
      BCAL_SAMPLING_FRACT );
  
  bcalInit(); //find inner cell threshold E
  
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return;
  
  for(unsigned int i=0; i<PE->mult; i++){ 
    
    float eUpStore[49][11][5] = {{{0}}}; //store energies [mod][lay][sect]
    float eDownStore[49][11][5]= {{{0}}};	  
    
    float eUpfADC[49][11][5]={{{0}}};
    float eDownfADC[49][11][5]={{{0}}};
    
    float tUpStore[49][11][5] = {{{0}}}; //store times [mod][lay][sect]
    float tDownStore[49][11][5]= {{{0}}};
    
    float tUpfADC[49][11][5]={{{0}}};
    float tDownfADC[49][11][5]={{{0}}};
    
    int fADCCellCount = 0;
    
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->barrelEMcal == HDDM_NULL ||
	hits->barrelEMcal->bcalCells == HDDM_NULL)continue;
    
    
    s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
    
    
    
    for(unsigned int j=0; j<cells->mult; j++){
      s_BcalCell_t *cell = &cells->in[j];
      
      // dcell is needed for dark hits cellID
      int dcell = DBCALGeometry::cellId( cell->module, cell->layer, cell->sector);
      darkHitItr = darkHits.find( dcell );
      
      
      //Create BCAL hits structure to put smeared data into
      
      if(cell->bcalSiPMUpHits!=HDDM_NULL)free(cell->bcalSiPMUpHits);
      cell->bcalSiPMUpHits = make_s_BcalSiPMUpHits(cell->bcalHits->mult);
      cell->bcalSiPMUpHits->mult = cell->bcalHits->mult;
      
      if(cell->bcalSiPMDownHits!=HDDM_NULL)free(cell->bcalSiPMDownHits);
      cell->bcalSiPMDownHits = make_s_BcalSiPMDownHits(cell->bcalHits->mult);
      cell->bcalSiPMDownHits->mult = cell->bcalHits->mult;			
      
      for(unsigned int k=0; k<cell->bcalHits->mult; k++){
	s_BcalHit_t *bcalhit = &cell->bcalHits->in[k];
	s_BcalSiPMUpHit_t *bcaluphit = &cell->bcalSiPMUpHits->in[k];
	s_BcalSiPMDownHit_t *bcaldownhit = &cell->bcalSiPMDownHits->in[k];
	
	float smearedE = bcalSamplingSmear( bcalhit->E );
	
	float upDist = ( bcalGeom->BCALFIBERLENGTH / 2 ) + bcalhit->zLocal;
	float downDist = ( bcalGeom->BCALFIBERLENGTH / 2 ) - bcalhit->zLocal;
	
	// sampling fluctuations are correlated between ends
	float upEnergy = smearedE * exp( -upDist / bcalGeom->ATTEN_LENGTH );
	float downEnergy = smearedE * exp( -downDist / bcalGeom->ATTEN_LENGTH );
	
	// independently smear time for both ends -- time smearing 
	// parameters come from data taken with beam at the center of 
	// the module so there is an implicit exp( ( -L / 2 ) / lambda ) 
	// that needs to be canceled out since we are working
	// at this stage with attenuated energies 
	float smearedtUp = 
	  bcalTimeSmear( bcalhit->t, 
			 upEnergy * exp( ( bcalGeom->BCALFIBERLENGTH / 2 ) / 
					 bcalGeom->ATTEN_LENGTH ) );
	float smearedtDown = 
	  bcalTimeSmear( bcalhit->t, 
			 downEnergy * exp( ( bcalGeom->BCALFIBERLENGTH / 2 ) / 
					   bcalGeom->ATTEN_LENGTH ) );
	
	darkHitItr = darkHits.find( dcell );
	if( darkHitItr != darkHits.end() ){
	  
	  upEnergy += ( darkHitItr->second.first * mevPerPE * k_MeV );
	  downEnergy += ( darkHitItr->second.second * mevPerPE * k_MeV );
	  
	  // now delete this from the map so we don't create
	  // additional hits later
	  darkHits.erase( darkHitItr );
	}
	
	// now offset times for propagation distance
	float upTime = smearedtUp + upDist / bcalGeom->C_EFFECTIVE;
	float downTime = smearedtDown + downDist / bcalGeom->C_EFFECTIVE;
	
	//If energy is smeared to negative, set to 0.
	if(upEnergy <= 0 || upTime!=upTime)
	  {
	    upEnergy = 0;
	    upTime = 0;
	  }
	if(downEnergy <= 0|| downTime!=downTime)
	  {
	    downEnergy = 0;
	    downTime = 0;
	  }
	
	bcaluphit->E = upEnergy;
	eUpStore[cell->module][cell->layer][cell->sector]+= upEnergy;
	bcaluphit->t = upTime;
	tUpStore[cell->module][cell->layer][cell->sector]= upTime;
	
	
	bcaldownhit->E = downEnergy;
	eDownStore[cell->module][cell->layer][cell->sector]+= downEnergy;
	bcaldownhit->t = downTime;
	tDownStore[cell->module][cell->layer][cell->sector]= downTime;
	
	
      } //k (bcal hits)
      
    } //j (cells)
    
    //Add in dark hits to empty cells
    for(int i = 1; i<=48;i++)
      {
	for(int j = 1; j<=10;j++)
	  {
	    for(int k = 1; k<=4; k++)
	      {
		int dcell = DBCALGeometry::cellId( i, j, k);
		darkHitItr = darkHits.find( dcell );
		if( darkHitItr != darkHits.end() ){
		  eUpStore[i][j][k] += ( darkHitItr->second.first * mevPerPE * k_MeV );
		  eDownStore[i][j][k] += ( darkHitItr->second.second * mevPerPE * k_MeV );
		  tUpStore[i][j][k] = SampleRange( -0.25 * BCAL_INTWINDOW_NS,
						   0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
		  tDownStore[i][j][k] = SampleRange( -0.25 * BCAL_INTWINDOW_NS,
						     0.75 * BCAL_INTWINDOW_NS ) * k_nsec;
		}
	      }
	  }
      }		     
    
    //Inner Cells, Summing
    for(int i=1;i<=48;i++)
      {
	for(int j=1;j<=bcalGeom->NBCALLAYS1;j++)
	  {
	    for(int k=1;k<=bcalGeom->NBCALSECS1;k++)
	      {
		for(int l=1;l<=(bcalGeom->BCALMID-1)/bcalGeom->NBCALLAYS1;l++)
		  {
		    for(int m=1;m<=4/bcalGeom->NBCALSECS1;m++)
		      {			 
			eUpfADC[i][j][k]+= eUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];//Sum energies
			eDownfADC[i][j][k]+= eDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			
			tUpfADC[i][j][k]+= eUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m] //Sum times weighted by Energy
			  *tUpStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			tDownfADC[i][j][k]+= eDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m]
			  *tDownStore[i][(j-1)*(bcalGeom->BCALMID-1)/(bcalGeom->NBCALLAYS1)+l][(k-1)*4/(bcalGeom->NBCALSECS1)+m];
			
		      }
		  }
		
		if(eUpfADC[i][j][k]!=0) //Divide time by Energy, to average weighted by energy, but make sure nonzero energy present.
		  tUpfADC[i][j][k] = tUpfADC[i][j][k]/eUpfADC[i][j][k];
		else
		  tUpfADC[i][j][k] = 0;
		
		if(eDownfADC[i][j][k]!=0)
		  tDownfADC[i][j][k] = tDownfADC[i][j][k]/eDownfADC[i][j][k];
		else
		  tDownfADC[i][j][k] = 0;
	      }
	  }
      }
    //Outer Cells, Summing
    for(int i=1;i<=48;i++)
      {
	for(int j=bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1;j++)
	  {
	    for(int k=1;k<=bcalGeom->NBCALSECS2;k++)
	      {
		for(int l=1;l<=(11-bcalGeom->BCALMID)/bcalGeom->NBCALLAYS2;l++)
		  {
		    for(int m=1;m<=4/bcalGeom->NBCALSECS2;m++)
		      {			 
			eUpfADC[i][j][k]+= eUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];//Sum energies
			eDownfADC[i][j][k]+= eDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];
			
			tUpfADC[i][j][k]+= eUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m] //Sum times weighted by energy
			  *tUpStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];
			tDownfADC[i][j][k]+= eDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m]
			  *tDownStore[i][(bcalGeom->BCALMID-1)+(j-bcalGeom->NBCALLAYS1-1)*(11-bcalGeom->BCALMID)/(bcalGeom->NBCALLAYS2)+l][(k-1)*4/(bcalGeom->NBCALSECS2)+m];			 
		      }
		  }
		
		if(eUpfADC[i][j][k]!=0) //Divide time by Energy, to average weighted by energy, but make sure nonzero energy present.
		  tUpfADC[i][j][k] = tUpfADC[i][j][k]/eUpfADC[i][j][k];
		else
		  tUpfADC[i][j][k] = 0;
		
		if(eDownfADC[i][j][k]!=0)
		  tDownfADC[i][j][k] = tDownfADC[i][j][k]/eDownfADC[i][j][k];
		else
		  tDownfADC[i][j][k] = 0;
	      }
	  }
      }
    
    //Naively multiply threshold by the square root of the number summed cells
    //in order to apply an effective 95% cut on fADC dark counts
    //passing.  Needs more thought.	
    float InnerThreshold = Bcal_CellInnerThreshold * 
      sqrt( ((bcalGeom->BCALMID-1)/bcalGeom->NBCALLAYS1) * (4/bcalGeom->NBCALSECS1) );
    float OuterThreshold = Bcal_CellInnerThreshold * 
      sqrt( ((11-bcalGeom->BCALMID)/bcalGeom->NBCALLAYS2) * (4/bcalGeom->NBCALSECS2) );	
    
    for(int i = 1;i<=48;i++)
      {
	for(int j = 1;j<=bcalGeom->NBCALLAYS1 ;j++)
	  {
	    for(int k = 1;k<=bcalGeom->NBCALSECS1;k++)
	      {
		if(eUpfADC[i][j][k]>InnerThreshold||eDownfADC[i][j][k]>InnerThreshold)
		  {		
		    fADCCellCount++;
		  }
	      }
	  }
      }
    for(int i = 1;i<=48;i++)
      {
	for(int j = bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1 ;j++)
	  {
	    for(int k = 1;k<=bcalGeom->NBCALSECS2 ;k++)
	      {
		if(eUpfADC[i][j][k]> OuterThreshold ||eDownfADC[i][j][k]>OuterThreshold)
		  {	
		    fADCCellCount++;
		  }
	      }
	  }
      }
    
    if(hits->barrelEMcal->bcalfADCCells!=HDDM_NULL)free(hits->barrelEMcal->bcalfADCCells);
    hits->barrelEMcal->bcalfADCCells = make_s_BcalfADCCells(fADCCellCount);
    hits->barrelEMcal->bcalfADCCells->mult = fADCCellCount;
    
    for(unsigned int j=0; j<hits->barrelEMcal->bcalfADCCells->mult; j++)
      {
	s_BcalfADCCell_t *fADCcell = &hits->barrelEMcal->bcalfADCCells->in[j];
	
	bool Etrigg = FALSE;		 
	
	//Inner Cells
	for(int i = 1;i<=48 && !Etrigg ;i++)
	  {
	    for(int j = 1;j<=bcalGeom->NBCALLAYS1 && !Etrigg ;j++)
	      {
		for(int k = 1;k<=bcalGeom->NBCALSECS1 && !Etrigg ;k++)
		  {
		    if(eUpfADC[i][j][k]>InnerThreshold||eDownfADC[i][j][k]>InnerThreshold)
		      {			       
			fADCcell->module = i;
			fADCcell->layer = j;
			fADCcell->sector = k;
			
			if(fADCcell->bcalfADCUpHits!=HDDM_NULL)free(fADCcell->bcalfADCUpHits);
			fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
			fADCcell->bcalfADCUpHits->mult = 1;
			
			if(fADCcell->bcalfADCDownHits!=HDDM_NULL)free(fADCcell->bcalfADCDownHits);
			fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
			fADCcell->bcalfADCDownHits->mult = 1;
			
			s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
			s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];				   
			
			if( eUpfADC[i][j][k] > InnerThreshold){
			  fadcuphit->E = eUpfADC[i][j][k];
			  fadcuphit->t = tUpfADC[i][j][k];
			}
			else {fadcuphit->E = 0; fadcuphit->t = 0;}
			
			if( eDownfADC[i][j][k] > InnerThreshold){
			  fadcdownhit->E = eDownfADC[i][j][k];
			  fadcdownhit->t = tDownfADC[i][j][k];
			}
			else {fadcdownhit->E = 0; fadcdownhit->t = 0;}			       
			
			Etrigg = TRUE;
			eUpfADC[i][j][k]=eDownfADC[i][j][k]=0;
		      }			        	   
		  }
	      }
	  }
	//Outer Cells
	for(int i = 1;i<=48 && !Etrigg ;i++)
	  {
	    for(int j = bcalGeom->NBCALLAYS1+1;j<=bcalGeom->NBCALLAYS2+bcalGeom->NBCALLAYS1 && !Etrigg ;j++)
	      {
		for(int k = 1;k<=bcalGeom->NBCALSECS2 && !Etrigg ;k++)
		  {
		    if(eUpfADC[i][j][k]> OuterThreshold ||eDownfADC[i][j][k]>OuterThreshold)
		      {			  
			fADCcell->module = i;
			fADCcell->layer = j;
			fADCcell->sector = k;
			
			if(fADCcell->bcalfADCUpHits!=HDDM_NULL)free(fADCcell->bcalfADCUpHits);
			fADCcell->bcalfADCUpHits = make_s_BcalfADCUpHits(1);
			fADCcell->bcalfADCUpHits->mult = 1;
			
			if(fADCcell->bcalfADCDownHits!=HDDM_NULL)free(fADCcell->bcalfADCDownHits);
			fADCcell->bcalfADCDownHits = make_s_BcalfADCDownHits(1);
			fADCcell->bcalfADCDownHits->mult = 1;
			
			s_BcalfADCUpHit_t *fadcuphit = &fADCcell->bcalfADCUpHits->in[0];
			s_BcalfADCDownHit_t *fadcdownhit = &fADCcell->bcalfADCDownHits->in[0];
			
			if( eUpfADC[i][j][k] > OuterThreshold){
			  fadcuphit->E = eUpfADC[i][j][k];
			  fadcuphit->t = tUpfADC[i][j][k];
			}
			else {fadcuphit->E = 0; fadcuphit->t = 0;}
			
			if( eDownfADC[i][j][k] > OuterThreshold){
			  fadcdownhit->E = eDownfADC[i][j][k];
			  fadcdownhit->t = tDownfADC[i][j][k];
			}
			else {fadcdownhit->E = 0; fadcdownhit->t = 0;}			       
			
			Etrigg = TRUE;
			eUpfADC[i][j][k]=eDownfADC[i][j][k]=0;
		      }			        	   
		  }
	      }
	  }
      }		                   		       
  } //i (PhysicsEvents)
  
}

void bcalInit()
{

// this the average number of (assumed single PE) 
  // dark pulses in a window

  double nAvg = BCAL_DARKRATE_GHZ * BCAL_INTWINDOW_NS;

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
  double darkTerm = exp( -nAvg );
  pdf[0] = darkTerm;

  for( int n_d = 1; n_d < 100; ++n_d ){

    darkTerm *= ( nAvg / n_d );

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

  double integral = 0;
  int nPEThresh = 0;
  while( integral < ( 1 - BCAL_MAXOCCUPANCY_FRACT ) ){

    // post increment includes zero and requires
    // one more PE than what breaks the loop
    integral += pdf[nPEThresh];
    ++nPEThresh;
  }

  // now get the photon theshold
  double photonThresh = nPEThresh / BCAL_DEVICEPDE;

  // now convert this into a energy threshold
  Bcal_CellInnerThreshold = 
    ( photonThresh / BCAL_PHOTONSPERSIDEPERMEV_INFIBER ) / 
    BCAL_SAMPLING_FRACT * k_MeV; 

}


float bcalSamplingSmear( float E )
{
    double sigmaSamp = BCAL_SAMPLINGCOEFA / sqrt( E ) + BCAL_SAMPLINGCOEFB;
    
    return( E * (1.0 + SampleGaussian(sigmaSamp)) );
}

float bcalTimeSmear( float t, float E )
{
  double sigmaT = BCAL_TIMEDIFFCOEFA / sqrt( E ) + BCAL_TIMEDIFFCOEFB;

  return( t + SampleGaussian(sigmaT) );
}

int getDarkHits()
{

  int darkPulse = SamplePoisson( BCAL_DARKRATE_GHZ* BCAL_INTWINDOW_NS );

  int xTalk = SamplePoisson( darkPulse * BCAL_XTALK_FRACT );

  return( xTalk + darkPulse );
}

//-----------
// SmearTOF
//-----------
void SmearTOF(s_HDDM_t *hddm_s)
{
  // Loop over Physics Events
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return;
  
  for(unsigned int i=0; i<PE->mult; i++){
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->forwardTOF == HDDM_NULL ||
	hits->forwardTOF->ftofCounters == HDDM_NULL)continue;
    
    
    s_FtofCounters_t* ftofCounters = hits->forwardTOF->ftofCounters;
    
    // Loop over counters
    
    for(unsigned int j=0;j<ftofCounters->mult; j++){
      s_FtofCounter_t *ftofCounter = &(ftofCounters->in[j]);		 
      
      // take care of North Hits
      s_FtofNorthTruthHits_t *ftofNorthTruthHits = ftofCounter->ftofNorthTruthHits;
      ftofCounter->ftofNorthHits = make_s_FtofNorthHits(ftofNorthTruthHits->mult);
      ftofCounter->ftofNorthHits->mult = ftofNorthTruthHits->mult;
      
      for (unsigned int m=0;m<ftofNorthTruthHits->mult;m++){
	s_FtofNorthTruthHit_t *ftofNorthTruthHit = &(ftofNorthTruthHits->in[m]);
	s_FtofNorthHit_t *ftofHit = &(ftofCounter->ftofNorthHits->in[m]);
	
	// Smear the time
	ftofHit->t = ftofNorthTruthHit->t + SampleGaussian(TOF_SIGMA);
	
	// Smear the energy
	double npe = (double)ftofNorthTruthHit->dE * 1000. *  TOF_PHOTONS_PERMEV;
	npe = npe +  SampleGaussian(sqrt(npe));
	float NewE = npe/TOF_PHOTONS_PERMEV/1000.;
	ftofHit->dE = NewE;
      } 
      
      // take care of South Hits
      s_FtofSouthTruthHits_t *ftofSouthTruthHits = ftofCounter->ftofSouthTruthHits;
      ftofCounter->ftofSouthHits = make_s_FtofSouthHits(ftofSouthTruthHits->mult);
      ftofCounter->ftofSouthHits->mult = ftofSouthTruthHits->mult;
      
      for (unsigned int m=0;m<ftofSouthTruthHits->mult;m++){
	s_FtofSouthTruthHit_t *ftofSouthTruthHit = &(ftofSouthTruthHits->in[m]);
	s_FtofSouthHit_t *ftofHit = &(ftofCounter->ftofSouthHits->in[m]);
	
	// Smear the time
	ftofHit->t = ftofSouthTruthHit->t + SampleGaussian(TOF_SIGMA);
	
	// Smear the energy
	double npe = (double)ftofSouthTruthHit->dE * 1000. *  TOF_PHOTONS_PERMEV;
	npe = npe +  SampleGaussian(sqrt(npe));
	float NewE = npe/TOF_PHOTONS_PERMEV/1000.;
	ftofHit->dE = NewE;

      }    
    } // end loop over all counters
  } 
}

//-----------
// SmearSTC - smear hits in the start counter
//-----------
void SmearSTC(s_HDDM_t *hddm_s)
{
  
 // Loop over Physics Events
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return;
  
  for(unsigned int i=0; i<PE->mult; i++){
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->startCntr == HDDM_NULL ||
	hits->startCntr->stcPaddles == HDDM_NULL)continue;
    s_StcPaddles_t *stcPaddles= hits->startCntr->stcPaddles;
    
    // Loop over counters
    for(unsigned int j=0;j<stcPaddles->mult; j++){
      s_StcPaddle_t *stcPaddle = &(stcPaddles->in[j]);
 
      s_StcTruthHits_t *stcTruthHits=stcPaddle->stcTruthHits;
      stcPaddle->stcHits=make_s_StcHits(stcTruthHits->mult);
      stcPaddle->stcHits->mult=stcTruthHits->mult;
      for (unsigned int m=0;m<stcTruthHits->mult;m++){
	s_StcTruthHit_t *stcTruthHit=&(stcTruthHits->in[m]);	
	s_StcHit_t *stcHit=&(stcPaddle->stcHits->in[m]);
	
	// smear the time
	stcHit->t=stcTruthHit->t + SampleGaussian(START_SIGMA);

	// Smear the energy
	double npe = (double)stcTruthHit->dE * 1000. *  START_PHOTONS_PERMEV;
	npe = npe +  SampleGaussian(sqrt(npe));
	float NewE = npe/START_PHOTONS_PERMEV/1000.;
	stcHit->dE = NewE;
      }
    }
  }
}

//-----------
// SmearCherenkov
//-----------
void SmearCherenkov(s_HDDM_t *hddm_s)
{



}

//-----------
// InitCDCGeometry
//-----------
void InitCDCGeometry(void)
{
  CDC_GEOMETRY_INITIALIZED = true;
  
  CDC_MAX_RINGS = 25;
  
  //-- This was cut and pasted from DCDCTrackHit_factory.cc on 10/11/2007 --
  
  float degrees0 = 0.0;
  float degrees6 = 6.0*M_PI/180.0;
  
  for(int ring=1; ring<=CDC_MAX_RINGS; ring++){
    int myNstraws=0;
    float radius = 0.0;
    float stereo=0.0;
    float phi_shift=0.0;
    float deltaX=0.0, deltaY=0.0;
    float rotX=0.0, rotY=0.0;
    switch(ring){
      // axial
    case  1:	myNstraws=  42;	radius= 10.7219;	stereo=  degrees0; phi_shift= 0.00000;	break;
    case  2:	myNstraws=  42;	radius= 12.097;	        stereo=  degrees0; phi_shift= 4.285714;	break;
    case  3:	myNstraws=  54;	radius= 13.7803;	stereo=  degrees0; phi_shift= 2.00000;	break;
    case  4:	myNstraws=  54;	radius= 15.1621;	stereo=  degrees0; phi_shift= 5.3333333;	break;
      
      // -stereo
    case  5:	myNstraws=  66;	radius= 16.9321;	stereo= -degrees6; phi_shift= 0.33333;	break;
    case  6:	myNstraws=  66;	phi_shift= 0.33333;	deltaX= 18.2948 ;	deltaY= 0.871486;	rotX=-6.47674;	rotY=-0.302853;	break;
    case  7:	myNstraws=  80;	radius= 20.5213;	stereo= -degrees6; phi_shift= -0.5000;	break;
    case  8:	myNstraws=  80;	phi_shift= -0.5000;	deltaX= 21.8912;	deltaY= 0.860106;	rotX=-6.39548;	rotY=-0.245615;	break;
      
      // +stereo
    case  9:	myNstraws=  93;	radius= 23.8544;	stereo= +degrees6; phi_shift= 1.1000;	break;
    case 10:	myNstraws=  93;	phi_shift= 1.1000;	deltaX= 25.229;	deltaY= 0.852573;	rotX=+6.34142;	rotY=+0.208647;	break;
    case 11:	myNstraws= 106;	radius= 27.1877;	stereo= +degrees6; phi_shift= -1.40;	break;
    case 12:	myNstraws= 106;	phi_shift= -1.400;	deltaX= 28.5658;	deltaY= 0.846871;	rotX=+6.30035;	rotY=+0.181146;	break;
      
      // axial
    case 13:	myNstraws= 123;	radius= 31.3799;	stereo=  degrees0; phi_shift= 0.5000000;	break;
    case 14:	myNstraws= 123;	radius= 32.7747;	stereo=  degrees0; phi_shift= 1.9634146;	break;
    case 15:	myNstraws= 135;	radius= 34.4343;	stereo=  degrees0; phi_shift= 1.0000000;	break;
    case 16:	myNstraws= 135;	radius= 35.8301;	stereo=  degrees0; phi_shift= 2.3333333;	break;
      
      // -stereo
    case 17:	myNstraws= 146;	radius= 37.4446;	stereo= -degrees6; phi_shift= 0.2;	break;
    case 18:	myNstraws= 146;	phi_shift= 0.2;	deltaX= 38.8295;	deltaY= 0.835653;	rotX=-6.21919 ;	rotY=-0.128247;	break;
    case 19:	myNstraws= 158;	radius= 40.5364;	stereo= -degrees6; phi_shift= 0.7;	break;
    case 20:	myNstraws= 158;	phi_shift= 0.7;	deltaX=41.9225 ;	deltaY= 0.833676;	rotX=-6.20274;	rotY=-0.118271;	break;
      
      // +stereo
    case 21:	myNstraws= 170;	radius= 43.6152;	stereo= +degrees6; phi_shift= 1.1000;	break;
    case 22:	myNstraws= 170;	phi_shift= 1.1000;	deltaX=45.0025  ;	deltaY= 0.831738;	rotX=+6.18859;	rotY=+0.109325;	break;
    case 23:	myNstraws= 182;	radius= 46.6849;	stereo= +degrees6; phi_shift= 1.40;	break;
    case 24:	myNstraws= 182;	phi_shift= 1.400;	deltaX= 48.0733;	deltaY= 0.829899;	rotX=+6.1763;	rotY=+0.101315;	break;
      
      // axial
    case 25:	myNstraws= 197;	radius= 50.37;	stereo=  degrees0; phi_shift= 0.200000000;	break;
    case 26:	myNstraws= 197;	radius= 51.77;	stereo=  degrees0; phi_shift= 1.113705000;	break;
    case 27:	myNstraws= 209;	radius= 53.363;	stereo=  degrees0; phi_shift= 0.800000000;	break;
    case 28:	myNstraws= 209;	radius= 54.76;	stereo=  degrees0; phi_shift= 1.661244;	break;
    default:
      cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-28 inclusive!"<<endl;
    }
    NCDC_STRAWS.push_back(myNstraws);
    CDC_RING_RADIUS.push_back(radius);
  }
  
  double Nstraws = 0;
  double alpha = 0.0;
  for(unsigned int i=0; i<NCDC_STRAWS.size(); i++){
    Nstraws += (double)NCDC_STRAWS[i];
    alpha += (double)NCDC_STRAWS[i]/CDC_RING_RADIUS[i];
  }
}


//-----------
// InitFDCGeometry
//-----------
void InitFDCGeometry(void)
{
  FDC_GEOMETRY_INITIALIZED = true;
  
  int FDC_NUM_LAYERS = 24;
  //int WIRES_PER_PLANE = 96;
  //int WIRE_SPACING = 1.116;
  
  for(int layer=1; layer<=FDC_NUM_LAYERS; layer++){
    
    float degrees00 = 0.0;
    float degrees60 = M_PI*60.0/180.0;
    
    float angle=0.0;
    float z_anode=212.0+95.5;
    switch(layer){
    case  1: z_anode+= -92.5-2.0;	angle=  degrees00; break;
    case  2: z_anode+= -92.5+0.0;	angle= +degrees60; break;
    case  3: z_anode+= -92.5+2.0;	angle= -degrees60; break;
    case  4: z_anode+= -86.5-2.0;	angle=  degrees00; break;
    case  5: z_anode+= -86.5+0.0;	angle= +degrees60; break;
    case  6: z_anode+= -86.5+2.0;	angle= -degrees60; break;
      
    case  7: z_anode+= -32.5-2.0;	angle=  degrees00; break;
    case  8: z_anode+= -32.5+0.0;	angle= +degrees60; break;
    case  9: z_anode+= -32.5+2.0;	angle= -degrees60; break;
    case 10: z_anode+= -26.5-2.0;	angle=  degrees00; break;
    case 11: z_anode+= -26.5+0.0;	angle= +degrees60; break;
    case 12: z_anode+= -26.5+2.0;	angle= -degrees60; break;
      
    case 13: z_anode+= +26.5-2.0;	angle=  degrees00; break;
    case 14: z_anode+= +26.5+0.0;	angle= +degrees60; break;
    case 15: z_anode+= +26.5+2.0;	angle= -degrees60; break;
    case 16: z_anode+= +32.5-2.0;	angle=  degrees00; break;
    case 17: z_anode+= +32.5+0.0;	angle= +degrees60; break;
    case 18: z_anode+= +32.5+2.0;	angle= -degrees60; break;
      
    case 19: z_anode+= +86.5-2.0;	angle=  degrees00; break;
    case 20: z_anode+= +86.5+0.0;	angle= +degrees60; break;
    case 21: z_anode+= +86.5+2.0;	angle= -degrees60; break;
    case 22: z_anode+= +92.5-2.0;	angle=  degrees00; break;
    case 23: z_anode+= +92.5+0.0;	angle= +degrees60; break;
    case 24: z_anode+= +92.5+2.0;	angle= -degrees60; break;
    }
    
    FDC_LAYER_Z.push_back(z_anode);
  }
  
  // Coefficient used to calculate FDCsingle wire rate. We calculate
  // it once here just to save calculating it for every wire in every event
  FDC_RATE_COEFFICIENT = exp(-log(4.0)/23.0)/2.0/log(24.0)*FDC_TIME_WINDOW/1000.0E-9;
  
  // Something is a little off in my calculation above so I scale it down via
  // an emprical factor:
  FDC_RATE_COEFFICIENT *= 0.353;
}

