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

// The straw diameter is 1.6cm
double CDC_R = 0.8;			// cm

// The stereo pitch is 6 degrees. For hit-based tracks, the
// hit could come anywhere in the 1.6cm diameter straw. Thus,
// the z-range corresponding to the straw diameter is
// 1.6cm * cot(6 degrees) = 15cm. 
double CDC_Z_SIGMA = 7.5;	// cm

// The occupancy of the CDC is used to determine how many noise
// hits to add per event. There are 3240 wires in the CDC. The
// Number of noise hits is sampled from a gaussian with a sigma
// equal to the sqrt(Nwires * CDCOccupancy).
double CDC_AVG_NOISE_HITS = 0.03*3240.0; // 0.03 = 3% occupancy
double CDC_OCCUPANCY = 0.03; // 0.03 = 3% occupancy
double CDC_TIME_WINDOW = 1000.0E-9; // in seconds

// For simulated data, we use a linear drift velocity function.
// This is the slope of that relation and should be aligned 
// with what is used by the simulation.
double CDC_DRIFT_VELOCITY = 55.0E-4; // Use number hardwired in simulation for now

// The position resolution of the CDC wires.
double CDC_POSITION_RESOLUTION = 150.0E-4; // in cm
 
// The wire spacing in the FDC is 0.5cm. In the actual data,
// multiple planes will be combined into psuedo points since
// each plane only measures in one direction (not counting z).
// At any rate, this should probably be an over estimate.
double FDC_R = 0.25;			// cm

// The z-coordinate of the wire and cathode planes are well
// known, but the hit will have some error on it. This should
// come primarily from the fact that the reconstructed x/y
// coordinates are at the DOCA point which will be off-plane.
// The actual error will come from the error on the track angle
// at the anode plane and the DOCA value measured by drift time.
// For here, we'll just use sigma=1mm.
double FDC_Z = 0.1;			// cm

// The occupancy of the FDC is used to determine how many noise
// hits to add per event. There are 2856 anode wires in the FDC.
// Each anode plane will have two cathode planes that will be
// used to resolve ambiguities and create "hits" that are passed
// on to the tracking code. We base the noise level on the
// number of anode wires.
double FDC_AVG_NOISE_HITS = 0.01*2856.0; // 0.01 = 1% occupancy

// Pedestal noise for FDC strips
double FDC_CATHODE_SIGMA= 200. ; // microns
double FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// Drift time variation for FDC anode wires
double FDC_DRIFT_SIGMA=200.0/55.0; // 200 microns/ (55 microns/ns)

// Drift time variation for CDC anode wires
double CDC_DRIFT_SIGMA=200.0/55.0; // 200 microns/ (55 microns/ns)

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

				// The resolution will depend on dE/dx which in turn depends on
				// the particle type and momentum. We don't have any of those
				// readily available here, but we do have the dE which could be
				// used to adjust the resolution of the hit.
				//
				// For now, we apply a simple 150 micron gaussian resolution
				// everywhere.
				double sigma_t = CDC_POSITION_RESOLUTION/CDC_DRIFT_VELOCITY;
				double delta_t = SampleGaussian(sigma_t);
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
			int iNhits = floor(Nhits);
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
			s_FdcTruthPoints_t *fdcTruthPoints = fdcChamber->fdcTruthPoints;
			if(fdcTruthPoints == HDDM_NULL)continue;
			
			s_FdcTruthPoint_t *truth = fdcTruthPoints->in;
			for(unsigned int k=0; k<fdcTruthPoints->mult; k++, truth++){
							
				// Take the simple approach here since
				// we don't know the orientation of the plane
				truth->x += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
				truth->y += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							
				// Z should be well known for the FDC, but smear it slightly anyway
				if(FDC_Z!=0.0)truth->z += SampleGaussian(FDC_Z);
			}
			
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
			      hit->t+=SampleGaussian(FDC_DRIFT_SIGMA);
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
	// The wire rates are obtained using a parameterization done
	// to calculate the event size for the August 29, 2007 online
	// meeting. This parameterization is almost already obsolete.
	// 11/8/2007 D. L.
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
			int iNhits = floor(Nhits);
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
				// Create single anode wire structure
				s_FdcAnodeWire_t *fdcAnodeWire = &fdcAnodeWires->in[fdcAnodeWires->mult++];

				// Create anode hits structure
				s_FdcAnodeHits_t *fdcanodehits = make_s_FdcAnodeHits(1);
				fdcAnodeWire->fdcAnodeHits = fdcanodehits;
				
				// Create single anode hit
				fdcanodehits->mult = 1;
				s_FdcAnodeHit_t *fdcanodehit = &fdcanodehits->in[0];
				
				fdcanodehit->dE = 0.1; // what should this be?
				fdcanodehit->t = SampleRange(0.0, 4.0);
				
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

	CDC_MAX_RINGS = 23;

	//-- This was cut and pasted from DCDCTrackHit_factory.cc on 10/11/2007 --

	float degrees0 = 0.0;
	float degrees6 = 6.0*M_PI/180.0;

	for(int ring=1; ring<=CDC_MAX_RINGS; ring++){
		int myNstraws=0;
		float radius = 0.0;
		float stereo=0.0;
		switch(ring){
			case  1:	myNstraws=  63;	radius= 16.049;	stereo=  degrees0; break;
			case  2:	myNstraws=  70;	radius= 17.831;	stereo=  degrees0; break;
			case  3:	myNstraws=  77;	radius= 19.613;	stereo=  degrees0; break;
			case  4:	myNstraws=  84;	radius= 21.395;	stereo=  degrees0; break;
			case  5:	myNstraws=  91;	radius= 23.178;	stereo= +degrees6; break;
			case  6:	myNstraws=  98;	radius= 24.960;	stereo= +degrees6; break;
			case  7:	myNstraws= 105;	radius= 26.742;	stereo= -degrees6; break;
			case  8:	myNstraws= 112;	radius= 28.524;	stereo= -degrees6; break;
			case  9:	myNstraws= 126;	radius= 32.089;	stereo=  degrees0; break;
			case 10:	myNstraws= 133;	radius= 33.871;	stereo=  degrees0; break;
			case 11:	myNstraws= 140;	radius= 35.654;	stereo=  degrees0; break;
			case 12:	myNstraws= 147;	radius= 37.435;	stereo=  degrees0; break;
			case 13:	myNstraws= 154;	radius= 39.218;	stereo=  degrees0; break;
			case 14:	myNstraws= 161;	radius= 41.001;	stereo= +degrees6; break;
			case 15:	myNstraws= 168;	radius= 42.783;	stereo= +degrees6; break;
			case 16:	myNstraws= 175;	radius= 44.566;	stereo= -degrees6; break;
			case 17:	myNstraws= 182;	radius= 46.348;	stereo= -degrees6; break;
			case 18:	myNstraws= 193;	radius= 49.149;	stereo=  degrees0; break;
			case 19:	myNstraws= 200;	radius= 50.932;	stereo=  degrees0; break;
			case 20:	myNstraws= 207;	radius= 52.714;	stereo=  degrees0; break;
			case 21:	myNstraws= 214;	radius= 54.497;	stereo=  degrees0; break;
			case 22:	myNstraws= 221;	radius= 56.279;	stereo=  degrees0; break;
			case 23:	myNstraws= 228;	radius= 58.062;	stereo=  degrees0; break;
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
	FDC_RATE_COEFFICIENT = exp(-log(4)/23.0)/2.0/log(24.0)*FDC_TIME_WINDOW/1000.0E-9;
	
	// Something is a little off in my calculation above so I scale it down via
	// an emprical factor:
	FDC_RATE_COEFFICIENT *= 0.353;
}
