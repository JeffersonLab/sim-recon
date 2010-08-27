// $Id$
//
//    File: DEventProcessor_trkres_tree.cc
// Created: Tue Apr  7 14:54:33 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include "DEventProcessor_trkres_tree.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
using namespace std;

#include <TROOT.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrajectoryPoint.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_trkres_tree());
}
} // "C"


//------------------
// DEventProcessor_track_hists
//------------------
DEventProcessor_trkres_tree::DEventProcessor_trkres_tree()
{
	trkres_ptr = &trkres;

	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_track_hists
//------------------
DEventProcessor_trkres_tree::~DEventProcessor_trkres_tree()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_trkres_tree::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("TRACKING");
	if(!dir)dir = new TDirectoryFile("TRACKING","TRACKING");
	dir->cd();

	ttrkres = new TTree("track","Track Res.");
	ttrkres->Branch("E","trackres",&trkres_ptr);

	dir->cd("../");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trkres_tree::brun(JEventLoop *loop, int runnumber)
{
	pthread_mutex_lock(&mutex);
	
	DApplication *dapp =  dynamic_cast<DApplication*>(loop->GetJApplication());
	if(dapp){
		bfield = dapp->GetBfield();
	}else{
		_DBG_<<"Unable to get DApplication pointer! (JApplication* = "<<loop->GetJApplication()<<")"<<endl;
		exit(-1);
	}
	
	// These are copied from DTrackFitterALT1.cc
	SIGMA_CDC = 0.0150;
	SIGMA_FDC_ANODE = 0.0200;
	SIGMA_FDC_CATHODE = 0.0200;

	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_CDC",						SIGMA_CDC);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_ANODE",				SIGMA_FDC_ANODE);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_CATHODE",			SIGMA_FDC_CATHODE);

	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trkres_tree::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrajectoryPoint*> trajpoints;
	vector<const DCDCTrackHit*> cdchits;
	vector<const DFDCPseudo*> fdchits;
	const DMCThrown *mcthrown;

	loop->Get(trajpoints);
	loop->Get(cdchits);
	loop->Get(fdchits);
	loop->GetSingle(mcthrown);

	// Assume all hits belong to this one thrown track
	// (this should only be used on data procuded with
	// NOSECONDARIES set to 1 !!)
	vector<meas_t> meas; // container to hold measurements

	// Loop over CDC hits, adding them to list
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCoordinateSystem *wire = cdchits[i]->wire;
		
		meas_t m;
		m.err = SIGMA_CDC;
		m.errc = 0.0;
		
		// Find trajectory point closest to this wire
		m.traj = FindTrajectoryPoint(wire, m.radlen, m.s, trajpoints);
		double Bx, By, Bz;
		bfield->GetField(m.traj->x, m.traj->y, m.traj->z, Bx, By, Bz);
		m.B = sqrt(Bx*Bx + By*By + Bz*Bz);
		
		meas.push_back(m);
	}

	// Loop over FDC hits, adding them to list
	for(unsigned int i=0; i<fdchits.size(); i++){
		const DCoordinateSystem *wire = fdchits[i]->wire;
		
		meas_t m;
		m.err = SIGMA_FDC_ANODE;
		m.errc = 0.0;
		
		// Find trajectory point closest to this wire
		m.traj = FindTrajectoryPoint(wire, m.radlen, m.s, trajpoints);
		double Bx, By, Bz;
		bfield->GetField(m.traj->x, m.traj->y, m.traj->z, Bx, By, Bz);
		m.B = sqrt(Bx*Bx + By*By + Bz*Bz);
		
		meas.push_back(m);
	}

	if(meas.size()<5)return NOERROR;

	double deltak, pt_res;
	GetPtRes(meas, deltak, pt_res);

	double theta_res;
	GetThetaRes(meas, theta_res);

	double pt = sqrt(pow((double)meas[0].traj->px,2.0) + pow((double)meas[0].traj->py,2.0));
	double p  = sqrt(pow((double)meas[0].traj->pz,2.0) + pow(pt,2.0));
	double theta = asin(pt/p);
	if(theta<0.0)theta+=2.0*M_PI;
	
	DVector3 dthrown = mcthrown->momentum();

	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	trkres.event = eventnumber;
	trkres.recon.SetXYZ(meas[0].traj->px, meas[0].traj->py, meas[0].traj->pz);
	trkres.thrown.SetXYZ(dthrown.X(), dthrown.Y(), dthrown.Z());
	trkres.deltak = deltak;
	trkres.pt_res = pt_res;
	trkres.p_res = (pt_res*pt + fabs(pt/tan(theta)*theta_res))/sin(theta)/p;
	trkres.theta_res = theta_res;
	trkres.phi_res = 0;

	ttrkres->Fill();

	// Unlock mutex
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// FindTrajectoryPoint
//------------------
const DMCTrajectoryPoint* DEventProcessor_trkres_tree::FindTrajectoryPoint(const DCoordinateSystem *wire, double &radlen, double &s, vector<const DMCTrajectoryPoint*> trajpoints)
{
	// Loop over all points, keeping track of the one closest to the wire
	const DMCTrajectoryPoint* best=NULL;
	double min_dist = 1.0E6;
	s = 0.0;
	radlen = 0.0;
	double best_radlen=0.0;
	double best_s=0.0;
	for(unsigned int i=0; i<trajpoints.size(); i++){
		const DMCTrajectoryPoint* traj = trajpoints[i];
		DVector3 d(traj->x-wire->origin.X(), traj->y-wire->origin.Y(), traj->z-wire->origin.Z());
		double d_dot_udir = d.Dot(wire->udir);
		double dist = sqrt(d.Mag2() - d_dot_udir*d_dot_udir);
		double dist_past_end = fabs(d_dot_udir) - wire->L/2.0;
		if(dist_past_end>0.0) dist = sqrt(dist*dist + dist_past_end*dist_past_end);
		
		s += traj->step;
		radlen += (double)traj->step/(double)traj->radlen;
		
		if(dist<min_dist){
			min_dist = dist;
			best = traj;
			best_s = s;
			best_radlen = radlen;
		}
	}
	
	// Copy integrated steps and radiation lengths for this step into references
	s = best_s;
	radlen = best_radlen;

	return best;
}

//------------------
// GetPtRes
//------------------
void DEventProcessor_trkres_tree::GetPtRes(vector<meas_t> &meas, double &deltak, double &pt_res)
{
	// Calculate using second method (eq. 28.47 PDG July 2008, pg. 309)
	double Vs1 = 0.0;
	double Vs2 = 0.0;
	double Vs3 = 0.0;
	double Vs4 = 0.0;
	double w = 0.0;
	double Bavg = 0.0;
	for(unsigned int i=0; i<meas.size(); i++){
		double s = meas[i].s - meas[0].s;
		double e = meas[i].err;
		
		w += 1.0/(e*e);
		Vs1 += s/(e*e);
		Vs2 += s*s/(e*e);
		Vs3 += s*s*s/(e*e);
		Vs4 += s*s*s*s/(e*e);
		Bavg += meas[i].B;
	}
	double Vss = (Vs2 - Vs1*Vs1)/w;
	double Vs2s2 = (Vs4 - Vs2*Vs2)/w;
	double Vss2 = (Vs3 - Vs2*Vs1)/w;
	
	double deltak_res = sqrt(4.0/w*Vss/(Vss*Vs2s2 - Vss2*Vss2));

	// Resolution due to multiple scattering
	double radlen = meas[meas.size()-1].radlen; // we want total integrated from vertex
	double L = meas[meas.size()-1].s - meas[0].s; // pathlength through detector
	double pt = sqrt(pow((double)meas[0].traj->px,2.0) + pow((double)meas[0].traj->py,2.0));
	double p  = sqrt(pow((double)meas[0].traj->pz,2.0) + pow(pt,2.0));
	double beta = p/sqrt(p*p + pow(0.13957,2.0)); // assume pion (could get this from meas[0].traj->part)
	double lambda = M_PI/2.0 - asin(pt/p);
	double deltak_ms = 0.016*sqrt(radlen)/(L*p*beta*pow(cos(lambda), 2.0));
	
	deltak = sqrt(deltak_res*deltak_res + deltak_ms*deltak_ms);

	// Use the average B-field to estimate the curvature based on the 
	// momentum at the first trajectory point.
	Bavg/=(double)meas.size();
	double k = 0.003*Bavg/pt;

	pt_res = deltak/k;
}

//------------------
// GetThetaRes
//------------------
void DEventProcessor_trkres_tree::GetThetaRes(vector<meas_t> &meas, double &theta_res)
{
	// Error due to position resolution.
	// This is based on the method outlined in Mark Ito's GlueX note 1015-v2, page 7.
	
	// This essentially is eq. 12 from Mark's paper, but in a different form. This
	// form I got off of Wikipedia and it includes individual errors. The wikipedia
	// one also divides by N-2 rather than N. (I believe this may be due to the 2
	// parameters of the linear fit).
	double s_avg = 0.0;
	double err2_avg = 0.0;
	for(unsigned int i=0; i<meas.size(); i++){
		double s = meas[i].s - meas[0].s;
		double e = meas[i].err;
		
		s_avg += s;
		err2_avg += e*e;
	}
	s_avg/=(double)meas.size();
	err2_avg/=(double)(meas.size()-2);
	
	double sum=0.0;
	for(unsigned int i=0; i<meas.size(); i++){
		double s = meas[i].s - meas[0].s;

		sum += pow(s - s_avg, 2.0);
	}

	// The difference from Mark's note and here is that in the note he had a factor
	// of 1/sec^2(theta) arising from the derivative of arctan(theta). However, in 
	// that example, the error bars are always in the y-direction which is not always
	// perpendicular to the "track". This leads to a zero contribution to the error
	// at pi/2 in his formula. Here, we take the result at theta=0 since that is where
	// the error direction is perpendicular to the track direction as it always is
	// in the case of real, helical tracks. Note that I'm tempted to take an arctan
	// of the dtheta_res here, but I think the derivative method is probably correct.
	// In the limit of the error being small compared to the track length, the small
	// angle approximation is valid so it doesn't really matter.
	double dtheta_res = sqrt(err2_avg/sum);
	
	// Error due to multiple scattering
	double radlen = meas[meas.size()-1].radlen; // we want total integrated from vertex
	double pt = sqrt(pow((double)meas[0].traj->px,2.0) + pow((double)meas[0].traj->py,2.0));
	double p  = sqrt(pow((double)meas[0].traj->pz,2.0) + pow(pt,2.0));
	double beta = p/sqrt(p*p + pow(0.13957,2.0)); // assume pion (could get this from meas[0].traj->part)
	double dtheta_ms = 0.0136*sqrt(radlen)/(beta*p)*(1.0+0.038*log(radlen))/sqrt(3.0);

	theta_res = sqrt(dtheta_res*dtheta_res + dtheta_ms*dtheta_ms);
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trkres_tree::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trkres_tree::fini(void)
{

	return NOERROR;
}

