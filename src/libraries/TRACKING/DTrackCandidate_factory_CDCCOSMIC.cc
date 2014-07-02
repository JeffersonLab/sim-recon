// $Id$
//
//    File: DTrackCandidate_factory_CDCCOSMIC.cc
// Created: Sat Jun 28 16:50:07 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.2.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <CDC/DCDCTrackHit.h>

#include "DTrackCandidate_factory_CDCCOSMIC.h"
using namespace jana;


// This simple track finder is made for fitting straight line cosmic tracks 
// in the CDC. It will automatically make a track from all CDC hits it finds. 
// It does this by doing a linear regression on the CDC axial hits
// in the X/Y plane.

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::init(void)
{
	bfield = new DMagneticFieldMapNoField(japp);
	rt = new DReferenceTrajectory(bfield);
	rt->Rmax_interior = 100.0; // (cm)  set larger swim volume so we can swim through the BCAL 
	rt->Rmax_exterior = 200.0; // (cm)  set larger swim volume so we can swim through the BCAL 

	// Optionally fill residual vs ring
	bool FILL_HIST = false;
	gPARMS->SetDefaultParameter("FILL_HIST", FILL_HIST);
	if(FILL_HIST){
		residual_vs_ring = new TH2D("residual_vs_ring", "CDC Cosmics straight line fitter", 250, 0.0, 100.0, 28, 0.5, 28.5);
		residual_vs_ring->SetXTitle("residual (mm)");
		residual_vs_ring->SetYTitle("ring (a.k.a. layer)");

		h_chisq = new TH1D("chisq", "#chi^2", 250, 0.0, 200.0);
		h_Ndof = new TH1D("Ndof", "Degrees of Freedom", 201, -0.5, 200.5);
		h_chisq_per_Ndof = new TH1D("chisq_per_Ndof", "#chi^2/Ndof", 250, 0.0, 10.0);
	}

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DCDCTrackHit*> cdchits;
	loop->Get(cdchits);
	
	// Sort CDC hits into axial and stereo
	vector<const DCDCTrackHit*> axial_hits;
	vector<const DCDCTrackHit*> stereo_hits;
	for(unsigned int i=0; i< cdchits.size(); i++){
		if( !cdchits[i]->is_stereo ){
			axial_hits.push_back(cdchits[i]);
		}else{
			stereo_hits.push_back(cdchits[i]);
		}
	}

	if( axial_hits.size() < 3 || stereo_hits.size() < 3) return NOERROR;
	
	// Cosmics are mostly virtical so parameterize x as a function of y
	double N=0.0, Sxy=0.0, Sx=0.0, Sy=0.0, Sy2=0.0;
	for(unsigned int i=0; i<axial_hits.size(); i++){
		const DCDCWire *wire = axial_hits[i]->wire;
		double x = wire->origin.X();
		double y = wire->origin.Y();
		
		N += 1.0;
		Sxy += x*y;
		Sx += x;
		Sy += y;
		Sy2 += y*y;

	}
	
	double m = (N*Sxy - Sx*Sy)/(N*Sy2 - Sy*Sy);
	double b = (Sx - m*Sy)/N;
	
	// alpha is angle representing direction of momentum vector in
	// lab coordinates. For cosmics, this will generally be pointing
	// downwards so alpha will be ~270 degrees.
	// adjust angle from x vs. y to pointing generally towards gravity in lab frame
	double alpha = 3.0*M_PI_2 - atan(m); 
	
	// Make momentum vector be unit vector in X/Y plane
	DVector3 mom(cos(alpha), sin(alpha), 0.0);

	// Find vertex location in X/Y by finding point were line crosses
	// R=90cm, choosing the value larger in Y
	double R = 90.0; // cm
	double A = 1.0 + m*m;
	double B = 2.0*m;
	double C = b*b - R*R;
	double sqroot = sqrt(B*B - 4.0*A*C);
	double y_vertex = (-B + sqroot)/(2.0*A);
	double x_vertex = y_vertex*m + b;

	// Loop over stereo wires, calculating intersection with track
	// Do this as a linear regression of z vs r
	N = 0.0;
	Sy = 0.0;
	Sy2 = 0.0;
	double Syz=0.0, Sz=0.0;
	vector<const DCDCTrackHit*> stereo_hits_used;
	for(unsigned int i=0; i<stereo_hits.size(); i++){
		const DCDCWire *wire = stereo_hits[i]->wire;

		// Need to find intersection in X/Y plane
		DVector2 wpos(wire->origin.X(), wire->origin.Y());
		DVector2 wdir(wire->udir.X(), wire->udir.Y());
		double wnorm = wdir.Mod();
		wdir /= wnorm;
		DVector2 tpos(x_vertex, y_vertex);
		DVector2 tdir(mom.X(), mom.Y());
		DVector2 pos_diff = wpos - tpos;
		
		// wpos, wdir, tpos, tdir are all 2D vectors while alpha and beta are 
		// scalers in the following:
		//
		// wpos + alpha*wdir = tpos + beta*tdir
		//
		// Solving for alpha gives:
		double alpha = (pos_diff.X()*tdir.Y() - pos_diff.Y()*tdir.X()) / (tdir.X()*wdir.Y() - tdir.Y()*wdir.X());
		
		// Scale alpha by 1/wnorm to get proper scaling factor for 3D direction vector
		alpha /= wnorm;
		
		// If the crossing point is beyond the end of the wire, skip this stereo hit
		if(fabs(alpha) > wire->L/2.0) continue;
		
		DVector3 poca_wire = wire->origin + alpha*wire->udir;
		double z = poca_wire.Z();
		double y = poca_wire.Y();

		N += 1.0;
		Syz += y*z;
		Sy += y;
		Sz += z;
		Sy2 += y*y;
		
		stereo_hits_used.push_back(stereo_hits[i]);
	}
	
	// It's possible not enough stereo straws survived for a fit
	if(N<3.0) return NOERROR;
	
	double m_yz = (N*Syz - Sy*Sz)/(N*Sy2 - Sy*Sy);
	double b_yz = (Sz - m_yz*Sy)/N;
	double theta = 3.0*M_PI_2 - atan(m_yz);
	
	double z_vertex = m_yz*y_vertex + b_yz;
	
	// Adjust momentum vector to use calculated theta
	mom *= fabs(sin(theta));
	mom.SetZ(cos(theta));

	// Make momentum vector be 100GeV so hdview2 will draw as stright line
	// (hdview2 creates its own rt using the field map)
	mom *= 100.0;
	
	DVector3 pos(x_vertex, y_vertex, z_vertex);
	
	DTrackCandidate *can = new DTrackCandidate();
	can->setMomentum(mom);
	can->setPosition(pos);
	can->setCharge(1.0); // make track draw as stright line
	can->setMass(0.139);
	can->setPID(PiPlus);
	can->rt = rt;

	rt->Swim(pos, mom, 1.0);
	
	CalcChisq(can, axial_hits, stereo_hits_used);
	
	_data.push_back(can);

	return NOERROR;
}

//------------------
// CalcChisq
//------------------
void DTrackCandidate_factory_CDCCOSMIC::CalcChisq(DTrackCandidate *can, vector<const DCDCTrackHit*> &axial_hits, vector<const DCDCTrackHit*> &stereo_hits)
{
	// Calculate the chi-sq and Ndof for the candidate based
	// on the wires in the two lists provided. Assume the
	// candiates rt member is valid and that all hits in the
	// two lists are used in the fit.
	//
	// Axial hits assume a resolution of the cell size over
	// sqrt(12). Stereo hits scle this up by 1/sin(6 degrees)
	
	// Ndof is easy
	can->Ndof = (int)(axial_hits.size() + stereo_hits.size()) - 4; // 4 is for fit parameters
	
	// Hit resolutions
	double sigma_axial = 0.8 / sqrt(12.0);
	double sigma_stereo = sigma_axial/cos(6.0/57.3); // 6 degree stereo angle
	
	can->chisq = 0.0;

	// Point on track and momentum direction
	DVector3 tpos = rt->swim_steps[0].origin;
	DVector3 tdir = rt->swim_steps[0].mom;
	tdir *= 1.0/tdir.Mag();

	// Axial
	for(unsigned int i=0; i<axial_hits.size(); i++){
		const DCDCWire *wire = axial_hits[i]->wire;
		
		// Distance between two skew lines is given by:
		//
		//  dist = pos_diff . n
		// i.e. dot product between distance between any two points on the two
		// lines and the normalized vector perpendicular to the two lines.
		DVector3 pos_diff = wire->origin - tpos;
		DVector3 n = wire->udir.Cross(tdir);
		double dist = pos_diff.Dot(n);
		
		//double dist = can->rt->DistToRT(wire);
		
		double dchi = dist/sigma_axial;
		can->chisq += dchi*dchi;
		
		if(residual_vs_ring)residual_vs_ring->Fill(fabs(dist)*10.0, wire->ring);
	}

	// Stereo
	for(unsigned int i=0; i<stereo_hits.size(); i++){
		const DCDCWire *wire = stereo_hits[i]->wire;
		
		DVector3 pos_diff = wire->origin - tpos;
		DVector3 n = wire->udir.Cross(tdir);
		double dist = pos_diff.Dot(n);

		//double dist = can->rt->DistToRT(wire);
		
		double dchi = dist/sigma_stereo;
		can->chisq += dchi*dchi;

		if(residual_vs_ring) residual_vs_ring->Fill(fabs(dist)*10.0, wire->ring);
	}

	if(h_chisq){
		h_chisq->Fill(can->chisq);
		h_Ndof->Fill((double)can->Ndof);
		h_chisq_per_Ndof->Fill(can->chisq/(double)can->Ndof);
	}
}

//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::fini(void)
{
	if(rt) delete rt;
	if(bfield) delete bfield;
	
	rt = NULL;
	bfield = NULL;

	return NOERROR;
}

