// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "hddm_s.h"

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

// The occupancy of the FDC is used to determine how many noise
// hits to add per event. There are 2856 anode wires in the FDC.
// Each anode plane will have two cathode planes that will be
// used to resolve ambiguities and create "hits" that are passed
// on to the tracking code. We base the noise level on the
// number of anode wires.
float FDC_AVG_NOISE_HITS = 0.01*2856.0; // 0.01 = 1% occupancy

//-----------
// Smear
//-----------
void Smear(s_HDDM_t *hddm_s)
{
	SmearCDC(hddm_s);
	AddNoiseHitsCDC(hddm_s);
	SmearFDC(hddm_s);
	AddNoiseHitsFDC(hddm_s);
	SmearFCAL(hddm_s);
	SmearBCAL(hddm_s);
	SmearTOF(hddm_s);
	SmearUPV(hddm_s);
	SmearCherenkov(hddm_s);
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
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(unsigned int j=0;j<rings->mult;j++){
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(unsigned int k=0;k<straws->mult;k++){
						s_CdcPoints_t *cdcpoints = straws->in[k].cdcPoints;
						if(cdcpoints){
							for(unsigned int m=0;m<cdcpoints->mult;m++){
								// Here we want to move the point to a position
								// randomly sampled from a disc in the x/y
								// plane centered on the actual point.
								float r = cdcpoints->in[m].r;
								float phi = cdcpoints->in[m].phi;
								float z = cdcpoints->in[m].z;
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
								z+= SampleGaussian(CDC_Z_SIGMA);

								cdcpoints->in[m].r = sqrt(x*x + y*y);
								cdcpoints->in[m].phi = atan2(y,x);
								if(cdcpoints->in[m].phi<0.0)cdcpoints->in[m].phi += 2.0*M_PI;
								cdcpoints->in[m].z = z;
							}
						}
					}
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
	// Since the current pattern finding algorithm only uses values
	// in the s_CdcPoint_t structure, it doesn't really matter what
	// ring,straw we put the noise hits under. Thus, we put them under
	// the first straw of the first ring. If there is no ring/straw, then
	// no hits are added.

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(!rings)return;
		if(!rings->mult)return;
		s_Straws_t *straws = rings->in[0].straws;
		if(!straws)return;
		if(!straws->mult)return;
		
		// Get existing hits
		int Nold = 0;
		s_CdcPoints_t *old_cdcpoints = straws->in[0].cdcPoints;
		if(old_cdcpoints)Nold = old_cdcpoints->mult;
		
		// How many noise hits to add
		int Nhits = (int)(CDC_AVG_NOISE_HITS + SampleGaussian(sqrt(CDC_AVG_NOISE_HITS)));
		s_CdcPoints_t* cdcpoints = make_s_CdcPoints(Nhits + Nold);
		cdcpoints->mult = 0;
		
		// Add real hits back in first
		if(old_cdcpoints){
			for(unsigned int i=0; i<old_cdcpoints->mult; i++){
				cdcpoints->in[cdcpoints->mult++] = old_cdcpoints->in[i];
			}
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_cdcpoints);
		straws->in[0].cdcPoints = cdcpoints;
		
		// Add noise hits
		for(int i=0; i<Nhits; i++){
			s_CdcPoint_t *cdc = &cdcpoints->in[cdcpoints->mult++];
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
	}
}

//-----------
// SmearFDC
//-----------
void SmearFDC(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Chambers_t *chambers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		
		for(unsigned int j=0;j<chambers->mult;j++){

			s_AnodePlanes_t *anodeplanes = chambers->in[j].anodePlanes;
			if(anodeplanes){
			
				for(unsigned int k=0;k<anodeplanes->mult;k++){
					s_Wires_t *wires = anodeplanes->in[k].wires;
					if(!wires)continue;
				
					for(unsigned int m=0;m<wires->mult;m++){
						s_FdcPoints_t *fdcPoints = wires->in[m].fdcPoints;
						if(!fdcPoints)continue;
						for(unsigned int n=0;n<fdcPoints->mult;n++){
							float x = fdcPoints->in[n].x;
							float y = fdcPoints->in[n].y;
							
							// Take the simple approach here since
							// we don't know the orientation of the plane
							x += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							y += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							
							// Z should be well known for the FDC
							fdcPoints->in[n].x = x;
							fdcPoints->in[n].y = y;
						}
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
	// Since the current pattern finding algorithm only uses values
	// in the s_FdcPoint_t structure, it doesn't really matter what
	// plane,wire we put the noise hits under. Thus, we put them under
	// the first wire of the first plane. If there is no plane/wire, then
	// no hits are added.

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Chambers_t *chambers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		if(!chambers->mult)return;
		s_AnodePlanes_t *anodeplanes = chambers->in[0].anodePlanes;
		if(!anodeplanes)return;
		if(!anodeplanes->mult)return;
		s_Wires_t *wires = anodeplanes->in[0].wires;
		if(!wires)continue;
		if(!wires->mult)continue;
		
		// Get existing hits
		int Nold = 0;
		s_FdcPoints_t *old_fdcpoints = wires->in[0].fdcPoints;
		if(old_fdcpoints)Nold = old_fdcpoints->mult;
		
		// How many noise hits to add
		int Nhits = (int)(FDC_AVG_NOISE_HITS + SampleGaussian(sqrt(FDC_AVG_NOISE_HITS)));
		s_FdcPoints_t* fdcpoints = make_s_FdcPoints(Nhits + Nold);
		fdcpoints->mult = 0;
		
		// Add real hits back in first
		if(old_fdcpoints){
			for(unsigned int i=0; i<old_fdcpoints->mult; i++){
				fdcpoints->in[fdcpoints->mult++] = old_fdcpoints->in[i];
			}
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_fdcpoints);
		wires->in[0].fdcPoints = fdcpoints;
		
		// Add noise hits
		for(int i=0; i<Nhits; i++){
			s_FdcPoint_t *fdc = &fdcpoints->in[fdcpoints->mult++];
			fdc->dEdx = 0.0;
			fdc->dradius = 0.0;
			fdc->primary = 1;
			fdc->track = 0;
			
			// To get a uniform sampling, we have to sample evenly in X/Y
			// until we find a point in the CDC fiducial area
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
