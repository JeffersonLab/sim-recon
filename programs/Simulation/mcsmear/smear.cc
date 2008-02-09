// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.h"

float RANDOM_MAX = (float)(0x7FFFFFFF);
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl

void SmearCDC(s_HDDM_t *hddm_s);
void AddNoiseHitsCDC(s_HDDM_t *hddm_s);
void SmearFDC(s_HDDM_t *hddm_s);
void AddNoiseHitsFDC(s_HDDM_t *hddm_s);
void SmearFCAL(s_HDDM_t *hddm_s);
void SmearBCAL(s_HDDM_t *hddm_s);
void SmearTOF(s_HDDM_t *hddm_s);
void SmearUPV(s_HDDM_t *hddm_s);
void SmearCherenkov(s_HDDM_t *hddm_s);
void InitCDCGeometry(void);
void InitFDCGeometry(void);

bool CDC_GEOMETRY_INITIALIZED = false;
int CDC_MAX_RINGS=0;
vector<unsigned int> NCDC_STRAWS;
vector<double> CDC_RING_RADIUS;

bool FDC_GEOMETRY_INITIALIZED = false;
unsigned int NFDC_WIRES_PER_PLANE;
vector<double> FDC_LAYER_Z;
double FDC_RATE_COEFFICIENT;


double SampleGaussian(double sigma);
double SampleRange(double x1, double x2);

// Do we or do we not add noise hits
bool ADD_NOISE = true;

// Do we or do we not smear real hits
bool SMEAR_HITS = true;

// The error on the drift time in the CDC. The drift times
// for the actual CDC hits coming from the input file
// are smeared by a gaussian with this sigma.
double CDC_TDRIFT_SIGMA = 5.0E-9;	// in seconds

// The time window for which CDC hits are accumulated.
// This is used to determine the number of background
// hits in the CDC for a given event.
double CDC_TIME_WINDOW = 1000.0E-9; // in seconds
 
// The error on the drift time in the CDC. The drift times
// for the actual CDC hits coming from the input file
// are smeared by a gaussian with this sigma.
double FDC_TDRIFT_SIGMA = 2.0E-9;	// in seconds

// The error in the distance along the wire as measured by
// the cathodes. This should NOT include the Lorentz
// effect which is already included in hdgeant. It
// should include any fluctuations due to ion trail density
// etc.
double FDC_CATHODE_SIGMA = 150.0; // in microns 

// The FDC pedestal noise is used to smear the cathode ADC
// values such that the position along the wire has the resolution
// specified by FDC_CATHODE_SIGMA.
double FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// Time window for acceptance of FDC hits
double FDC_TIME_WINDOW = 1000.0E-9; // in seconds


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
	if(SMEAR_HITS)SmearBCAL(hddm_s);
	if(SMEAR_HITS)SmearTOF(hddm_s);
	if(SMEAR_HITS)SmearUPV(hddm_s);
	if(SMEAR_HITS)SmearCherenkov(hddm_s);
}

//-----------
// SmearCDC
//-----------
void SmearCDC(s_HDDM_t *hddm_s)
{
	/// Smear the drift times of all CDC hits.

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
			for(unsigned int j=0; j<cdcstraw->cdcStrawHits->mult; j++){
				s_CdcStrawHit_t *strawhit = &cdcstraw->cdcStrawHits->in[j];

				// Smear out the CDC drift time using the specified sigma.
				// This should include both timing resolution and ion trail
				// density effects.
				double delta_t = SampleGaussian(CDC_TDRIFT_SIGMA);
				strawhit->t += delta_t;
				
				// If the time is negative, reject this smear and try again
				if(strawhit->t<0)j--;
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
				strawhit->t = SampleRange(-CDC_TIME_WINDOW/2.0, +CDC_TIME_WINDOW/2.0);
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
			    s_FdcCathodeHits_t *hits=strip->fdcCathodeHits;
			    if (hits==HDDM_NULL)continue;
			    s_FdcCathodeHit_t *hit=hits->in;
			    for (unsigned int s=0;s<hits->mult;s++,hit++){
			      hit->q+=SampleGaussian(FDC_PED_NOISE);
			    }
			  }
			}

			// Add drift time varation to the anode data 
			s_FdcAnodeWires_t *wires=fdcChamber->fdcAnodeWires;
			
			if (wires!=HDDM_NULL){
			  s_FdcAnodeWire_t *wire=wires->in;
			  for (unsigned int k=0;k<wires->mult;k++,wire++){
			    s_FdcAnodeHits_t *hits=wire->fdcAnodeHits;
			    if (hits==HDDM_NULL)continue;
			    s_FdcAnodeHit_t *hit=hits->in;
			    for (unsigned int s=0;s<hits->mult;s++,hit++){
			      hit->t+=SampleGaussian(FDC_TDRIFT_SIGMA);
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
				fdcanodehit->t = SampleRange(0.0, FDC_TIME_WINDOW);
				
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



}
//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{



}

//-----------
// SmearTOF
//-----------
void SmearTOF(s_HDDM_t *hddm_s)
{



}

//-----------
// SmearUPV
//-----------
void SmearUPV(s_HDDM_t *hddm_s)
{



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
		switch(ring){
			case  1:	myNstraws=  43;	radius= 10.960;	stereo=  degrees0; break;
			case  2:	myNstraws=  50;	radius= 12.741;	stereo=  degrees0; break;
			case  3:	myNstraws=  57;	radius= 14.522;	stereo=  degrees0; break;
			case  4:	myNstraws=  64;	radius= 16.304;	stereo= -degrees6; break;
			case  5:	myNstraws=  71;	radius= 18.086;	stereo= -degrees6; break;
			case  6:	myNstraws=  78;	radius= 19.868;	stereo= +degrees6; break;
			case  7:	myNstraws=  85;	radius= 21.65 ;	stereo= +degrees6; break;
			case  8:	myNstraws=  99;	radius= 25.214;	stereo=  degrees0; break;
			case  9:	myNstraws= 106;	radius= 26.997;	stereo=  degrees0; break;
			case 10:	myNstraws= 113;	radius= 28.779;	stereo=  degrees0; break;
			case 11:	myNstraws= 120;	radius= 30.561;	stereo=  degrees0; break;
			case 12:	myNstraws= 127;	radius= 32.344;	stereo=  degrees0; break;
			case 13:	myNstraws= 134;	radius= 34.126;	stereo= -degrees6; break;
			case 14:	myNstraws= 141;	radius= 35.908;	stereo= -degrees6; break;
			case 15:	myNstraws= 148;	radius= 37.691;	stereo= +degrees6; break;
			case 16:	myNstraws= 155;	radius= 39.473;	stereo= +degrees6; break;
			case 17:	myNstraws= 166;	radius= 42.274;	stereo=  degrees0; break;
			case 18:	myNstraws= 173;	radius= 44.057;	stereo=  degrees0; break;
			case 19:	myNstraws= 180;	radius= 45.839;	stereo=  degrees0; break;
			case 20:	myNstraws= 187;	radius= 47.621;	stereo=  degrees0; break;
			case 21:	myNstraws= 194;	radius= 49.404;	stereo=  degrees0; break;
			case 22:	myNstraws= 201;	radius= 51.186;	stereo=  degrees0; break;
			case 23:	myNstraws= 208;	radius= 52.969;	stereo=  degrees0; break;
			case 24:	myNstraws= 215;	radius= 54.751;	stereo=  degrees0; break;
			case 25:	myNstraws= 222;	radius= 56.534;	stereo=  degrees0; break;
			default:
				cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-23 inclusive!"<<endl;
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
