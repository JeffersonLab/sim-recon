// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.h"

float RANDOM_MAX = (float)(0x7FFFFFFF);

void SmearCDC(s_HDDM_t *hddm_s);
void AddNoiseHitsCDC(s_HDDM_t *hddm_s);
void SmearFDC(s_HDDM_t *hddm_s);
void AddNoiseHitsFDC(s_HDDM_t *hddm_s);
void SmearFCAL(s_HDDM_t *hddm_s);
void SmearBCAL(s_HDDM_t *hddm_s);
void SmearTOF(s_HDDM_t *hddm_s);
void SmearUPV(s_HDDM_t *hddm_s);
void SmearCherenkov(s_HDDM_t *hddm_s);

double SampleGaussian(double sigma);
double SampleRange(double x1, double x2);

// Do we or do we not add noise hits
bool ADD_NOISE = false;

// Do we or do we not smear real hits
bool SMEAR_HITS = true;

// The straw diameter is 1.6cm
float CDC_R = 0.8;			// cm

// The stereo pitch is 6 degrees. For hit-based tracks, the
// hit could come anywhere in the 1.6cm diameter straw. Thus,
// the z-range corresponding to the straw diameter is
// 1.6cm * cot(6 degrees) = 15cm. 
float CDC_Z_SIGMA = 7.5;	// cm

// The occupancy of the CDC is used to determine how many noise
// hits to add per event. There are 3240 wires in the CDC. The
// Number of noise hits is sampled from a gaussian with a sigma
// equal to the sqrt(Nwires * CDCOccupancy).
float CDC_AVG_NOISE_HITS = 0.03*3240.0; // 0.03 = 3% occupancy

// The wire spacing in the FDC is 0.5cm. In the actual data,
// multiple planes will be combined into psuedo points since
// each plane only measures in one direction (not counting z).
// At any rate, this should probably be an over estimate.
float FDC_R = 0.25;			// cm

// The z-coordinate of the wire and cathode planes are well
// known, but the hit will have some error on it. This should
// come primarily from the fact that the reconstructed x/y
// coordinates are at the DOCA point which will be off-plane.
// The actual error will come from the error on the track angle
// at the anode plane and the DOCA value measured by drift time.
// For here, we'll just use sigma=1mm.
float FDC_Z = 0.1;			// cm

// The occupancy of the FDC is used to determine how many noise
// hits to add per event. There are 2856 anode wires in the FDC.
// Each anode plane will have two cathode planes that will be
// used to resolve ambiguities and create "hits" that are passed
// on to the tracking code. We base the noise level on the
// number of anode wires.
float FDC_AVG_NOISE_HITS = 0.01*2856.0; // 0.01 = 1% occupancy

// Pedestal noise for FDC strips
float FDC_CATHODE_SIGMA= 200. ; // microns
float FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// Drift time variation for FDC anode wires
float FDC_DRIFT_SIGMA=200.0/55.0; // 200 microns/ (55 microns/ns)

// Drift time variation for CDC anode wires
float CDC_DRIFT_SIGMA=200.0/55.0; // 200 microns/ (55 microns/ns)


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
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->centralDC == HDDM_NULL ||
			hits->centralDC->cdcTruthPoints == HDDM_NULL)continue;
		
		s_CdcTruthPoints_t *cdctruthpoints = hits->centralDC->cdcTruthPoints;
		s_CdcTruthPoint_t *cdctruthpoint = cdctruthpoints->in;
		for(unsigned int j=0; j<cdctruthpoints->mult; j++, cdctruthpoint++){

			// Here we want to move the point to a position
			// randomly sampled from a disc in the x/y
			// plane centered on the actual point.
			float r = cdctruthpoint->r;
			float phi = cdctruthpoint->phi;
			float z = cdctruthpoint->z;
			float x = r*cos(phi);
			float y = r*sin(phi);
								
			// since this needs to be evenly distributed over
			// the disc of radius CDC_R, we may try several times
			float deltaX, deltaY;
			do{
				float sx = (float)random()/RANDOM_MAX;
				deltaX = 2.0*(sx-0.5)*CDC_R;
				float sy = (float)random()/RANDOM_MAX;
				deltaY = 2.0*(sy-0.5)*CDC_R;
			}while(sqrt(deltaX*deltaX + deltaY*deltaY)>CDC_R);
			x+=deltaX;
			y+=deltaY;
			if(CDC_Z_SIGMA!=0.0)z+= SampleGaussian(CDC_Z_SIGMA);

			cdctruthpoint->r = sqrt(x*x + y*y);
			cdctruthpoint->phi = atan2(y,x);
			if(cdctruthpoint->phi<0.0)cdctruthpoint->phi += 2.0*M_PI;
			cdctruthpoint->z = z;
		}

		// Add drift time varation to the anode data 
		if (hits->centralDC->cdcStraws != HDDM_NULL){
			for(unsigned int k=0; k<hits->centralDC->cdcStraws->mult; k++){
				s_CdcStraw_t *cdcstraw = &hits->centralDC->cdcStraws->in[k];
				for(unsigned int j=0; j<cdcstraw->cdcStrawHits->mult; j++){
					s_CdcStrawHit_t *strawhit = &cdcstraw->cdcStrawHits->in[j];

					strawhit->t+=SampleGaussian(CDC_DRIFT_SIGMA);
				}
			}
		}
	}
}

//-----------
// AddNoiseHitsCDC
//-----------
void AddNoiseHitsCDC(s_HDDM_t *hddm_s)
{
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
		
		if(hits->centralDC->cdcTruthPoints == HDDM_NULL){
			hits->centralDC->cdcTruthPoints = make_s_CdcTruthPoints(0);
			hits->centralDC->cdcTruthPoints->mult=0;
		}

		// Get existing hits
		s_CdcTruthPoints_t *old_cdctruthpoints = hits->centralDC->cdcTruthPoints;
		unsigned int Nold = old_cdctruthpoints->mult;
		
		// How many noise hits to add
		int Nhits = (int)(CDC_AVG_NOISE_HITS + SampleGaussian(sqrt(CDC_AVG_NOISE_HITS)));
		if(Nhits<0)Nhits=0;
		s_CdcTruthPoints_t* cdctruthpoints = make_s_CdcTruthPoints((unsigned int)Nhits + Nold);

		// Add real hits back in first
		cdctruthpoints->mult = 0;
		for(unsigned int j=0; j<Nold; j++){
			cdctruthpoints->in[cdctruthpoints->mult++] = old_cdctruthpoints->in[j];
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_cdctruthpoints);
		hits->centralDC->cdcTruthPoints = cdctruthpoints;
		
		// Add noise hits. We add 1/3 of the hits evenly distributed
		// throughout the volume. The remaining 2/3 of the noise hits
		// will be distributed evenly in r and phi which gives a 
		// 1/r density distribution
		int j;
		for(j=0; j<Nhits/3; j++){
			s_CdcTruthPoint_t *cdc = &cdctruthpoints->in[cdctruthpoints->mult++];
			cdc->dEdx = 0.0;
			cdc->dradius = 0.0;
			
			// To get a uniform sampling, we have to sample evenly in X/Y
			// until we find a point in the CDC fiducial area
			double x,y,r = 0.0;
			while(r<13.0 || r>59.0){
				x = SampleRange(-60.0, 60.0);
				y = SampleRange(-60.0, 60.0);
				r = sqrt(x*x + y*y);
			}
			
			cdc->phi = M_PI + atan2(y,x);
			cdc->primary = 1;
			cdc->r = r;
			cdc->track = 0;
			cdc->z = SampleRange(17.0, 217.0);
		}
		for(; j<Nhits; j++){
			s_CdcTruthPoint_t *cdc = &cdctruthpoints->in[cdctruthpoints->mult++];
			cdc->dEdx = 0.0;
			cdc->dradius = 0.0;
			cdc->phi = SampleRange(0.0, 2.0*M_PI);
			cdc->primary = 1;
			cdc->r = SampleRange(0.0, 59.0);
			cdc->track = 0;
			cdc->z = SampleRange(17.0, 217.0);
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
			      hit->dE+=SampleGaussian(FDC_PED_NOISE);
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
			hits->forwardDC->fdcChambers = make_s_FdcChambers(1);
			hits->forwardDC->fdcChambers->mult=1;
			hits->forwardDC->fdcChambers->in[0].fdcAnodeWires = (s_FdcAnodeWires_t*)HDDM_NULL;
			hits->forwardDC->fdcChambers->in[0].fdcCathodeStrips = (s_FdcCathodeStrips_t*)HDDM_NULL;
			hits->forwardDC->fdcChambers->in[0].fdcTruthPoints = (s_FdcTruthPoints_t*)HDDM_NULL;
		}

		if(hits->forwardDC->fdcChambers->in[0].fdcTruthPoints == HDDM_NULL){
			hits->forwardDC->fdcChambers->in[0].fdcTruthPoints = make_s_FdcTruthPoints(0);
			hits->forwardDC->fdcChambers->in[0].fdcTruthPoints->mult=0;
		}

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		s_FdcChamber_t *fdcChamber = fdcChambers->in;
		for(unsigned int j=0; j<fdcChambers->mult; j++, fdcChamber++){
			s_FdcTruthPoints_t *old_fdcTruthPoints = fdcChamber->fdcTruthPoints;
			if(old_fdcTruthPoints == HDDM_NULL)continue;
			
			// Get existing hits
			unsigned int Nold = old_fdcTruthPoints->mult;
			
			// How many noise hits to add
			float Nchamber_noise_hits = FDC_AVG_NOISE_HITS/fdcChambers->mult;
			int Nhits = (int)(Nchamber_noise_hits + SampleGaussian(sqrt(Nchamber_noise_hits)));
			if(Nhits<0)Nhits = 0;
			s_FdcTruthPoints_t* fdcTruthPoints = make_s_FdcTruthPoints((unsigned int)Nhits + Nold);
			fdcTruthPoints->mult = 0;
	
			// Add real hits back in first
			s_FdcTruthPoint_t *truth = fdcTruthPoints->in;
			for(unsigned int k=0; k<Nold; k++, truth++){
				fdcTruthPoints->in[fdcTruthPoints->mult++] = old_fdcTruthPoints->in[k];
			}
		
			// Delete memory used for old hits structure and
			// replace pointer in HDDM tree with ours
			free(old_fdcTruthPoints);
			fdcChamber->fdcTruthPoints = fdcTruthPoints;

			// Add noise hits. We add 1/3 of the hits evenly distributed
			// throughout the volume. The remaining 2/3 of the noise hits
			// will be distributed evenly in r and phi which gives a 
			// 1/r density distribution
			int k;
			for(k=0; k<Nhits/3; k++){
				s_FdcTruthPoint_t *fdc = &fdcTruthPoints->in[fdcTruthPoints->mult++];
				fdc->dEdx = 0.0;
				fdc->dradius = 0.0;
				fdc->primary = 1;
				fdc->track = 0;
			
				// To get a uniform sampling, we have to sample evenly in X/Y
				// until we find a point in the FDC fiducial area
				double x,y,r = 0.0;
				while(r<3.5 || r>59.0){
					x = SampleRange(-60.0, 60.0);
					y = SampleRange(-60.0, 60.0);
					r = sqrt(x*x + y*y);
				}
				fdc->x = x;
				fdc->y = y;
				fdc->z = 240.0 + SampleRange(-6.0, +6.0);
				fdc->z += ((float)(random()%4)) * 52.0;
			}
			// 1/r distributed noise hits
			for(; k<Nhits; k++){
				s_FdcTruthPoint_t *fdc = &fdcTruthPoints->in[fdcTruthPoints->mult++];
				fdc->dEdx = 0.0;
				fdc->dradius = 0.0;
				fdc->primary = 1;
				fdc->track = 0;

				double phi = SampleRange(0.0, 2.0*M_PI);
				double r = SampleRange(0.0, 60.0);

				fdc->x = r*cos(phi);
				fdc->y = r*sin(phi);
				fdc->z = 240.0 + SampleRange(-6.0, +6.0);
				fdc->z += ((float)(random()%4)) * 52.0;
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
