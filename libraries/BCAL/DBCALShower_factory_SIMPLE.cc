// $Id$
//
//    File: DBCALShower_factory_SIMPLE.cc
// Created: Tue Mar 20 14:44:12 EDT 2007
// Creator: davidl (on Darwin swire-b241.jlab.org 8.9.0 powerpc)
//

#include <cmath>
using namespace std;

#include <TROOT.h>

#include <DANA/DApplication.h>
#include "DBCALShower_factory_SIMPLE.h"


bool pshowerSort_C(const DBCALShower_factory_SIMPLE::pshower_t &hit1, const DBCALShower_factory_SIMPLE::pshower_t &hit2) {
	return hit1.E > hit2.E;
}


//------------------
// DBCALShower_factory_SIMPLE (Constructor)
//------------------
DBCALShower_factory_SIMPLE::DBCALShower_factory_SIMPLE()
{
	UP_DOWN_COINCIDENCE_WINDOW = 50.0; // in ns
	ENERGY_SCALE_FACTOR = 1.0; // scale factor for converting energy to GeV (1.1 is empirical)
	SIGNAL_VELOCITY = 16.75; // in cm/ns
	Z_CENTER = 17.0+0.5+390.0/2.0; // in cm (0.5 is empirical)
	MIN_CLUSTER_SPACING = 50.0; // in cm
	MIN_SHOWER_ENERGY = 0.001; // in GeV
	DEBUG_HISTS = true;

	gPARMS->SetDefaultParameter("BCAL:UP_DOWN_COINCIDENCE_WINDOW",	UP_DOWN_COINCIDENCE_WINDOW);
	gPARMS->SetDefaultParameter("BCAL:ENERGY_SCALE_FACTOR",	ENERGY_SCALE_FACTOR);
	gPARMS->SetDefaultParameter("BCAL:SIGNAL_VELOCITY",	SIGNAL_VELOCITY);
	gPARMS->SetDefaultParameter("BCAL:Z_CENTER",	Z_CENTER);
	gPARMS->SetDefaultParameter("BCAL:MIN_CLUSTER_SPACING",	MIN_CLUSTER_SPACING);
	gPARMS->SetDefaultParameter("BCAL:MIN_SHOWER_ENERGY",	MIN_SHOWER_ENERGY);
	gPARMS->SetDefaultParameter("BCAL:DEBUG_HISTS",	DEBUG_HISTS);

	hit_element_dist=NULL;
}

//------------------
// brun
//------------------
jerror_t DBCALShower_factory_SIMPLE::brun(JEventLoop *eventLoop, int runnumber)
{
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());

	if(DEBUG_HISTS){
		dapp->Lock();
		
		// Histograms may already exist. (Another thread may have created them)
		// Try and get pointers to the existing ones.
		hit_element_dist = (TH1F*)gROOT->FindObject("hit_element_dist");
		x_vs_y_vs_z = (TH3F*)gROOT->FindObject("x_vs_y_vs_z");
		r_shower = (TH1F*)gROOT->FindObject("r_shower");
		phi_shower = (TH1F*)gROOT->FindObject("phi_shower");
		z_shower = (TH1F*)gROOT->FindObject("z_shower");
		E_shower = (TH1F*)gROOT->FindObject("E_shower");
		r_element = (TH1F*)gROOT->FindObject("r_element");
		phi_element = (TH1F*)gROOT->FindObject("phi_element");
		z_element = (TH1F*)gROOT->FindObject("z_element");
		E_element = (TH1F*)gROOT->FindObject("E_element");
		E_upstream_vs_downstream = (TH2F*)gROOT->FindObject("E_upstream_vs_downstream");
		E_clustersum = (TH1F*)gROOT->FindObject("E_clustersum");
		E_elementsum = (TH1F*)gROOT->FindObject("E_elementsum");

		if(!hit_element_dist)hit_element_dist = new TH1F("hit_element_dist","Distance between hit elements(cm)",300, 0.0, 300.0);
		if(!x_vs_y_vs_z)x_vs_y_vs_z = new TH3F("x_vs_y_vs_z","X,Y,Z distribution of shower centers",100, -100.0, 100.0, 100, -100.0, 100.0, 100, -50.0, 550.0);
		if(!r_shower)r_shower = new TH1F("r_shower","radial distribution of shower centers",500, 0.0, 100.0);
		if(!phi_shower)phi_shower = new TH1F("phi_shower","#phi distribution of shower centers",360, 0.0, 360.0);
		if(!z_shower)z_shower = new TH1F("z_shower","z distribution of shower centers",700, -50.0, 650.0);
		if(!E_shower)E_shower = new TH1F("E_shower","Energy of shower (GeV)",1000, 0.0, 8.0);
		if(!r_element)r_element = new TH1F("r_element","radial distribution of element centers",500, 0.0, 100.0);
		if(!phi_element)phi_element = new TH1F("phi_element","#phi distribution of element centers",360, 0.0, 360.0);
		if(!z_element)z_element = new TH1F("z_element","z distribution of element centers",700, -50.0, 650.0);
		if(!E_element)E_element = new TH1F("E_element","Energy of element (GeV)",1000, 0.0, 2.0);
		if(!E_upstream_vs_downstream)E_upstream_vs_downstream = new TH2F("E_upstream_vs_downstream","Amplitude of upstream readout vs. downsteam readout for single elements",100, 0.0, 1.0, 100, 0.0, 1.0);
		if(!E_clustersum)E_clustersum = new TH1F("E_clustersum","Energy sum of all clusters (GeV)",1000, 0.0, 8.0);
		if(!E_elementsum)E_elementsum = new TH1F("E_elementsum","Energy sum of all elements (GeV)",1000, 0.0, 8.0);

		dapp->Unlock();
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALShower_factory_SIMPLE::evnt(JEventLoop *loop, int eventnumber)
{
	// Get individual hits
	vector<const DBCALHit*> bcalhits;
	loop->Get(bcalhits);
	
	// Sort hits into upstream and downstream hits
	vector<const DBCALHit*> upstream_bcalhits;
	vector<const DBCALHit*> downstream_bcalhits;
	for(unsigned int i=0; i<bcalhits.size(); i++){
		const DBCALHit *hit = bcalhits[i];
		if(hit->end == DBCALHit::UPSTREAM){upstream_bcalhits.push_back(hit); continue;}
		if(hit->end == DBCALHit::DOWNSTREAM){downstream_bcalhits.push_back(hit); continue;}
		
		_DBG_<<"Heyyy! How'd I get here?!?"<<endl;
	}

	// Loop over hits finding coincidences between upstream and downstream
	// hits. Put them in pshower_t objects.
	vector<bool> up_used(upstream_bcalhits.size(), false);
	vector<bool> down_used(downstream_bcalhits.size(), false);
	vector<pshower_t> pshowers;
	for(unsigned int i=0; i<upstream_bcalhits.size(); i++){
		const DBCALHit *up = upstream_bcalhits[i];
		for(unsigned int j=0; j<downstream_bcalhits.size(); j++){
			const DBCALHit *down = downstream_bcalhits[j];

			// Check that we're looking at the same module, layer, sector
			if(up->module != down->module)continue;
			if(up->layer  != down->layer)continue;
			if(up->sector != down->sector)continue;
			
			// Check timing
			double tdiff = up->t - down->t;
			if(fabs(tdiff)>UP_DOWN_COINCIDENCE_WINDOW)continue;
			
			up_used[i] = true;
			down_used[j] = true;

			// Looks like a match! Add to list of coincidences
			pshower_t pshower;
			pshower.upstream_hit = up;
			pshower.downstream_hit = down;
			pshower.E = sqrt(up->E * down->E)*ENERGY_SCALE_FACTOR;
			pshower.t = (up->t + down->t)/2.0;
			pshower.z = Z_CENTER + SIGNAL_VELOCITY*(up->t - down->t)/2.0;
			ModuleLayerSectorToPhiR(up->module, up->layer, up->sector, pshower.phi, pshower.R);
			pshower.used = false;

			pshowers.push_back(pshower);

			if(DEBUG_HISTS){
				r_element->Fill(pshower.R);
				phi_element->Fill(pshower.phi*180.0/M_PI);
				z_element->Fill(pshower.z);
				E_element->Fill(pshower.E);
				E_upstream_vs_downstream->Fill(down->E, up->E);
			}
		}
	}
	
	// Some hits are not used because a matching hit on the other end
	// was not found. This can represent a signifacnt amount of energy
	// in the cluster. Here, we loop over the hits and turn any un-used
	// single end hits into pshowers. Note that for now, we assume
	// knowledge of the attenutation length, signal propagation velocity,
	// and time offset. All of these should be known to some degree for
	// real data.
	for(unsigned int i=0; i<upstream_bcalhits.size(); i++){
		if(up_used[i])continue; // skip hits already used above
		const DBCALHit *up = upstream_bcalhits[i];
		
		pshower_t pshower;
		pshower.upstream_hit = up;
		pshower.downstream_hit = NULL;
		pshower.E = up->E*ENERGY_SCALE_FACTOR;
		pshower.t = up->t;
		pshower.z = Z_CENTER + SIGNAL_VELOCITY*up->t;
		ModuleLayerSectorToPhiR(up->module, up->layer, up->sector, pshower.phi, pshower.R);
		pshower.used = false;

		pshowers.push_back(pshower);
	}
	for(unsigned int i=0; i<downstream_bcalhits.size(); i++){
		if(down_used[i])continue; // skip hits already used above
		const DBCALHit *down = downstream_bcalhits[i];
		
		pshower_t pshower;
		pshower.upstream_hit = NULL;
		pshower.downstream_hit = down;
		pshower.E = down->E*ENERGY_SCALE_FACTOR;
		pshower.t = down->t;
		pshower.z = Z_CENTER - SIGNAL_VELOCITY*down->t;
		ModuleLayerSectorToPhiR(down->module, down->layer, down->sector, pshower.phi, pshower.R);
		pshower.used = false;

		pshowers.push_back(pshower);
	}

	// Sort pshowers by energy
	sort(pshowers.begin(), pshowers.end(), pshowerSort_C);

	// Sort pshowers into clusters. In this simple algorithm, the
	// largest energy pshower is assumed to be near a cluster center.
	// All pshowers within a certain distance of this are added to
	// the cluster and marked as used. The cycle starts again with
	// the largest energy, un-used pshower and continues until all
	// pshowers have been made into clusters
	for(unsigned int i=0; i<pshowers.size(); i++){
		if(pshowers[i].used)continue;
		pshowers[i].used = true;
		
		// An unused pshower, start a new cluster
		vector<const pshower_t*> cluster;
		cluster.clear();
		cluster.push_back(&pshowers[i]);
		for(unsigned int j=i+1; j<pshowers.size(); j++){
			if(pshowers[j].used)continue;
			
			double &R1 = pshowers[i].R;
			double &R2 = pshowers[j].R;
			double &phi1 = pshowers[i].phi;
			double &phi2 = pshowers[j].phi;
			double &z1 = pshowers[i].z;
			double &z2 = pshowers[j].z;
			double dist = sqrt(R1*R1 + R2*R2 - 2.0*R1*R2*cos(phi1-phi2) + pow(z1-z2, 2.0));
			if(DEBUG_HISTS)hit_element_dist->Fill(dist);
			if(dist>MIN_CLUSTER_SPACING)continue;
			
			pshowers[j].used = true;
			cluster.push_back(&pshowers[j]);
		}
		
		if(cluster.size()==0)break;
	
		// Loop over pshowers in this cluster and calculate logarithmically weighted means
		double x=0.0, y=0.0, z=0.0;
		double norm=0.0;
		double tmin = 1.0E6;
		double Etot = 0.0;
		for(unsigned int j=0; j<cluster.size(); j++){
			const pshower_t *pshower = cluster[j];
			double logE = log(pshower->E*1000.0);
			//double logE = pshower->E;
			norm += logE;
			x += pshower->R*cos(pshower->phi)*logE;
			y += pshower->R*sin(pshower->phi)*logE;
			z += pshower->z*logE;
			if(pshower->t<tmin)tmin = pshower->t;
			Etot += pshower->E;
		}
		if(Etot<MIN_SHOWER_ENERGY)continue;
		
		x /= norm;
		y /= norm;
		z /= norm;
		
		// Create a DBCALShower from this cluster
		DBCALShower *shower = new DBCALShower();
		shower->E = Etot;
		shower->Ecorr = Etot;
		shower->x = x;
		shower->y = y;
		shower->z = z;
		shower->t = tmin;
		
		_data.push_back(shower);
		
		if(DEBUG_HISTS){
			r_shower->Fill(sqrt(x*x + y*y));
			double phi = atan2(y,x)*180.0/M_PI;
			if(phi<0.0)phi+=360.0;
			phi_shower->Fill(phi);
			z_shower->Fill(z);
			x_vs_y_vs_z->Fill(x,y,z);
			E_shower->Fill(Etot);
		}

	}
	
	// A few more debugging histos
	if(DEBUG_HISTS){
		// Total energy in all clusters(showers)
		double Etot = 0.0;
		for(unsigned int i=0; i<_data.size(); i++)Etot+=_data[i]->E;
		E_clustersum->Fill(Etot);
		
		// Total energy in all pshowers
		Etot = 0.0;
		for(unsigned int i=0; i<pshowers.size(); i++)Etot+=pshowers[i].E;
		E_elementsum->Fill(Etot);
	}

	return NOERROR;
}

//------------------
// ModuleLayerSectorToPhiR
//------------------
void DBCALShower_factory_SIMPLE::ModuleLayerSectorToPhiR(int module, int layer, int sector, double &phi, double &R)
{
	// This is where the geometry of the BCAL readout is implemented.
	// At this point in time (it could change!) the geomeetry in HDDS implements
	// what is described in GlueX-doc-708. Namely:
	//
	// - A 4x5 division of the first 10cm (sectors x layers)
	// - a 3x4 division of the last 12.46 cm
	// - layers go radially out starting at R=65.0cm
	// - the bcal is made up of 48 identical modules (7.5 degrees each)
	// - the first module is centered at phi=3.75 degrees
	
	// Check that values are in range
	if(module<1 || module>48){
		_DBG_<<"Module out of range: "<<module<<endl;
		return;
	}
	if(layer<1 || layer>9){
		_DBG_<<"Layer out of range: "<<layer<<endl;
		return;
	}
	
	// Phi angle of center of module
	double phi0 = ((double)(module-1)*7.5 + 3.75)*M_PI/180.0;
	
	if(layer<=5){
		// in first group
		if(sector<1 || sector>4){
			_DBG_<<"Sector out of range: "<<sector<<endl;
			return;
		}
		double phi_sector = 7.5/4.0*M_PI/180.0;
		double dphi = ((double)sector-2.5)*phi_sector;
		phi = phi0 + dphi;
		R=65.0+2.0*((double)layer-0.5)/cos(dphi);
	}else{
		if(sector<1 || sector>3){
			_DBG_<<"Sector out of range: "<<sector<<endl;
			return;
		}
		double phi_sector = 7.5/3.0*M_PI/180.0;
		double dphi = ((double)sector-2.0)*phi_sector;
		phi = phi0 + dphi;
		R=75.0+(12.46/4)*((double)layer-5.5)/cos(dphi);
	}
}

//------------------
// toString
//------------------
const string DBCALShower_factory_SIMPLE::toString(void)
{

	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)
        return string(); // don't print anything if we have no data!

  printheader("row:      x:      y:      z:       t:     E_seen:     E_corr:");

  for(unsigned int i = 0; i < _data.size(); i++) {
		DBCALShower *s = _data[i];
    
		printnewrow();
		printcol("%d",	i);
		printcol("%5.2f",	s->x);
		printcol("%5.2f",	s->y);
		printcol("%5.2f",	s->z);
		printcol("%5.3f",	s->t);
		printcol("%5.3f",	s->E);
		printcol("%5.3f",	s->Ecorr);
		printrow();
	}

	return _table;

}
