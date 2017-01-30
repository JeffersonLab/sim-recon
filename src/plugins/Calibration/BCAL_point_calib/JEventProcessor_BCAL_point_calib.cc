// $Id$
//
//    File: JEventProcessor_BCAL_point_calib.cc
// Created: Mon Sep 26 09:38:23 EDT 2016
// Creator: gvasil (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_point_calib.h"

#include <TLorentzVector.h>
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "BCAL/DBCALGeometry.h"
#include <vector>
#include "TGraphErrors.h"
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace jana;
using namespace std;

		// histograms for storing matched shower info
		static TH1I* h1_Num_matched_showers = NULL;
                static TH1I* h1_E = NULL;
		static TH1I* h1_E_raw = NULL;
		static TH1I* h1_x = NULL;
		static TH1I* h1_y = NULL;
		static TH1I* h1_z = NULL;
		static TH1I* h1_t = NULL;
		static TH1I* h1_N_cell = NULL;

		// histograms for storing matched track info
                static TH1I* h1trk_Num_matched_tracks = NULL;
		static TH1I* h1trk_FOM = NULL;
		static TH1I* h1trk_pmag = NULL;
		static TH1I* h1trk_px = NULL;
		static TH1I* h1trk_py = NULL;
		static TH1I* h1trk_pz = NULL;			
		static TH1I* h1trk_x = NULL;
		static TH1I* h1trk_y = NULL;
		static TH1I* h1trk_z = NULL;
		static TH1I* h1trk_energy = NULL;

		static TH2I* h2_Evsenergy = NULL;
		static TH2I* h2_pmagvsenergy = NULL;

		// histograms for storing point info
		static TH1I* h1pt_Num_points = NULL;
		static TH1I* h1pt_module = NULL;
		static TH1I* h1pt_layer = NULL;
		static TH1I* h1pt_sector = NULL;
		static TH1I* h1pt_cell_id = NULL;
		static TH1I* h1pt_energy = NULL;
		static TH1I* h1pt_energy_US = NULL;
		static TH1I* h1pt_energy_DS = NULL;
		static TH1I* h1pt_z = NULL;
		static TH1I* h1pt_r = NULL;
		static TH1I* h1pt_r_size = NULL;
		static TH1I* h1pt_phi = NULL;
		static TH1I* h1pt_t = NULL;

		// histograms to check z_track and z_point correlation
		static TH1I* h1_ztrack_minus_zpoint = NULL;
		static TH1I* h1_abs_ztrack_minus_zpoint = NULL;
		static TH2I* h2_ztrack_minus_zpoint_vs_ztrack = NULL;

		static TH2I* h2_zpoint_vs_ztrack = NULL;
		static TH2I* h2_zpoint_vs_ztrack_thrown = NULL;

		// create a graph for each channel (cell)
		static TGraph* h2_tgraph[768] = {NULL};
		// graphs of thrown points for each channel (cell)
		static TGraph* h2_tgraph_thrown[768] = {NULL};

		// histogram for the percentage of thrown points per channel
		static TH1D* h1_thrown_per_channel = NULL;

		int event_count=0;

		// vectors for storing the z_track and z_point for each channel
		vector< vector< double > > z_track_graph(768, vector< double >() );
		vector< vector< double > > z_point_graph(768, vector< double >() );
		// vectors for storing the thrown points for each channel (if cuts are applied)
		vector< vector< double > > z_track_thrown(768, vector< double >() );
		vector< vector< double > > z_point_thrown(768, vector< double >() );

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_point_calib());
}
} // "C"


//------------------
// JEventProcessor_BCAL_point_calib (Constructor)
//------------------
JEventProcessor_BCAL_point_calib::JEventProcessor_BCAL_point_calib()
{

}

//------------------
// ~JEventProcessor_BCAL_point_calib (Destructor)
//------------------
JEventProcessor_BCAL_point_calib::~JEventProcessor_BCAL_point_calib()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_point_calib::init(void)
{
	// This is called once at program startup. 
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	
	// japp->RootFillLock(this); 
        DEBUG=false;
	VERBOSE=false;
        gPARMS->SetDefaultParameter("BCALPOINTCALIB:DEBUG", DEBUG, "Control the creation of extra histograms");
        gPARMS->SetDefaultParameter("BCALPOINTCALIB:DEBUG", VERBOSE, "Control verbose output");

	if(h1_Num_matched_showers != NULL){
		return NOERROR;
	}
	
 	// lock all root operations
 	japp->RootWriteLock();

	 //	TDirectory *dir = new TDirectoryFile("BCAL","BCAL");
	 //	dir->cd();
	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcal_point_calibs")->cd();


        int nbins = 100;
	int nbins2 = 500;


	h1_Num_matched_showers = new TH1I("h1_Num_matched_showers", "Number of matched showers per event", 20,0,20);
	h1_Num_matched_showers->SetXTitle("Number of showers");
	h1_Num_matched_showers->SetYTitle("counts");
	  
	h1trk_Num_matched_tracks = new TH1I("h1trk_Num_matched_tracks", "Number of matched tracks per event", 20,0,20);
	h1trk_Num_matched_tracks->SetXTitle("Number of tracks");
	h1trk_Num_matched_tracks->SetYTitle("counts");

	// thrown points per channel (if cuts are applied)
	h1_thrown_per_channel = new TH1D("h1_thrown_per_channel","Percentage of thrown points per Channel",800,0,800);
	h1_thrown_per_channel->SetXTitle("Channel #");
	h1_thrown_per_channel->SetYTitle("Percentage (%)");
	h1_thrown_per_channel->SetMarkerStyle(20);

	h1pt_Num_points = new TH1I("h1pt_Num_points", "pt points per shower",20,0,20);
	h1pt_Num_points->SetXTitle("Points per Shower");
	h1pt_Num_points->SetYTitle("counts");

	if(DEBUG) {
	  
	  // shower info
	  h1_E = new TH1I("h1_E", "Matched energy per shower",nbins,0,2);
	  h1_E->SetXTitle("Energy per shower (GeV)");
	  h1_E->SetYTitle("counts");
	  h1_E_raw = new TH1I("h1_E_raw", "Matched raw energy per shower",nbins,0,2);
	  h1_E_raw->SetXTitle("Raw Energy per shower (GeV)");
	  h1_E_raw->SetYTitle("counts");
	  h1_x = new TH1I("h1_x", "x pos of shower",nbins2,-100,100);
	  h1_x->SetXTitle("x (cm)");
	  h1_x->SetYTitle("counts");
	  h1_y = new TH1I("h1_y", "y pos of shower",nbins2,-100,100);
	  h1_y->SetXTitle("y (cm)");
	  h1_y->SetYTitle("counts");
	  h1_z = new TH1I("h1_z", "z pos of shower",nbins2,0,500);
	  h1_z->SetXTitle("z (cm)");
	  h1_z->SetYTitle("counts");
	  h1_t = new TH1I("h1_t", "time of shower",nbins2,-100,100);
	  h1_t->SetXTitle("time of shower (ns)");
	  h1_t->SetYTitle("counts");
	  h1_N_cell = new TH1I("h1_N_cell", "N_cell of shower",20,0,20);
	  h1_N_cell->SetXTitle("N_cell");
	  h1_N_cell->SetYTitle("counts");

	  // track info
	  h1trk_FOM = new TH1I("h1trk_FOM", "FOM for matched tracks", nbins,0,1);
	  h1trk_FOM->SetXTitle("FOM");
	  h1trk_FOM->SetYTitle("counts");
	  h1trk_pmag = new TH1I("h1trk_mag", "pmag for matched tracks", nbins,0,4);
	  h1trk_pmag->SetXTitle("pmag (GeV)");
	  h1trk_pmag->SetYTitle("counts");
	  h1trk_px = new TH1I("h1trk_px", "px for matched tracks", nbins,0,2);
	  h1trk_px->SetXTitle("px (GeV)");
	  h1trk_px->SetYTitle("counts");
	  h1trk_py = new TH1I("h1trk_py", "py for matched tracks", nbins,0,2);
	  h1trk_py->SetXTitle("Number of tracks");
	  h1trk_py->SetYTitle("counts");
	  h1trk_pz = new TH1I("h1trk_pz", "pz for matched tracks", nbins,0,4);
	  h1trk_pz->SetXTitle("pz (GeV)");
	  h1trk_pz->SetYTitle("counts");
	  h1trk_x = new TH1I("h1trk_x", "x for matched tracks", nbins2,-30,30);
	  h1trk_x->SetXTitle("x (cm)");
	  h1trk_x->SetYTitle("counts");
	  h1trk_y = new TH1I("h1trk_y", "y for matched tracks", nbins2,-30,30);
	  h1trk_y->SetXTitle("y (cm)");
	  h1trk_y->SetYTitle("counts");
	  h1trk_z = new TH1I("h1trk_z", "z for matched tracks", nbins2,-50,250);
	  h1trk_z->SetXTitle("z (cm)");
	  h1trk_z->SetYTitle("counts");
	  h1trk_energy = new TH1I("h1trk_energy", "energy for matched tracks", nbins,0,4);
	  h1trk_energy->SetXTitle("energy (GeV)");
	  h1trk_energy->SetYTitle("counts");
	  
	  h2_Evsenergy = new TH2I("h2_Evsenergy", "E vs energy matched tracks", nbins,0,4,nbins,0,4);
	  h2_Evsenergy->SetXTitle("track energy (GeV)");
	  h2_Evsenergy->SetYTitle("Eshower (GeV)");
	  h2_pmagvsenergy = new TH2I("h2_pmagvsenergy", "pmag vs energy matched tracks", nbins,0,4,nbins,0,4);
	  h2_pmagvsenergy->SetXTitle("track energy (GeV)");
	  h2_pmagvsenergy->SetYTitle("track pmag (GeV)");
	  
	  // point info
	  h1pt_module = new TH1I("h1pt_module", "pt module number",50,0,50);
	  h1pt_module->SetXTitle("Module");
	  h1pt_module->SetYTitle("counts");
	  h1pt_layer = new TH1I("h1pt_layer", "pt layer number",5,0,5);
	  h1pt_layer->SetXTitle("Layer");
	  h1pt_layer->SetYTitle("counts");
	  h1pt_sector = new TH1I("h1pt_sector", "pt sector number",5,0,5);
	  h1pt_sector->SetXTitle("Sector");
	  h1pt_sector->SetYTitle("counts");
	  h1pt_cell_id = new TH1I("h1pt_cell_id", "pt cell id number",800,0,800);
	  h1pt_cell_id->SetXTitle("Cell ID");
	  h1pt_cell_id->SetYTitle("counts");
	  h1pt_energy = new TH1I("h1pt_energy", "pt energy",nbins,0,1);
	  h1pt_energy->SetXTitle("Point Energy");
	  h1pt_energy->SetYTitle("counts");
	  h1pt_energy_US = new TH1I("h1pt_energy_US", "pt energy US",nbins,0,1);
	  h1pt_energy_US->SetXTitle("Point Energy US");
	  h1pt_energy_US->SetYTitle("counts");
	  h1pt_energy_DS = new TH1I("h1pt_energy_DS", "pt energy DS",nbins,0,1);
	  h1pt_energy_DS->SetXTitle("Point Energy DS");
	  h1pt_energy_DS->SetYTitle("counts");
	  h1pt_z = new TH1I("h1pt_z", "z for point",nbins2,0,500);
	  h1pt_z->SetXTitle("z (cm)");
	  h1pt_z->SetYTitle("counts");
	  h1pt_r = new TH1I("h1pt_r", "r for point",nbins,60,90);
	  h1pt_r->SetXTitle("r (cm)");
	  h1pt_r->SetYTitle("counts");
	  h1pt_r_size = new TH1I("h1pt_r_size", "r size of point", 20, 0, 20);
	  h1pt_r_size->SetXTitle("r size (cm)");
	  h1pt_r_size->SetYTitle("counts");
	  h1pt_phi = new TH1I("h1pt_phi", "phi of point", 700, 0, 7);
	  h1pt_phi->SetXTitle("phi (rad)");
	  h1pt_phi->SetYTitle("counts");
	  h1pt_t = new TH1I("h1pt_t", "t for point",nbins2,-100,100);
	  h1pt_t->SetXTitle("t");
	  h1pt_t->SetYTitle("counts");
	  
	  // z_track and z_point correlation
	  h1_ztrack_minus_zpoint = new TH1I("h1_ztrack_minus_zpoint", "z_track - z_point", 400, -200, 200);
	  h1_ztrack_minus_zpoint->SetXTitle("z_track - z_point (cm)");
	  h1_ztrack_minus_zpoint->SetYTitle("counts");

	  h1_abs_ztrack_minus_zpoint = new TH1I("h1_abs_ztrack_minus_zpoint", "|z_track - z_point|", 400, 0, 400);
	  h1_abs_ztrack_minus_zpoint->SetXTitle("|z_track - z_point| (cm)");
	  h1_abs_ztrack_minus_zpoint->SetYTitle("counts");
	  h2_ztrack_minus_zpoint_vs_ztrack = new TH2I("h2_ztrack_minus_zpoint_vs_ztrack", "z_track - z_point, vs z_track", 500, -100, 400, 400, -200, 200); 
	  h2_ztrack_minus_zpoint_vs_ztrack->SetXTitle("z_track (cm)");
	  h2_ztrack_minus_zpoint_vs_ztrack->SetYTitle("z_track - z_point (cm)");
	  h2_zpoint_vs_ztrack = new TH2I("h2_zpoint_vs_ztrack", "z point vs z track for the BCAL",1000,0,500,1000,0,500);
	  h2_zpoint_vs_ztrack->SetXTitle("z of track (cm)");
	  h2_zpoint_vs_ztrack->SetYTitle("z of point (cm)"); 
	  h2_zpoint_vs_ztrack_thrown = new TH2I("h2_zpoint_vs_ztrack_thrown", "thrown z point vs z track for the BCAL",1000,0,500,1000,0,500);
	  h2_zpoint_vs_ztrack_thrown->SetXTitle("z of track (cm)");
	  h2_zpoint_vs_ztrack_thrown->SetYTitle("z of point (cm)"); 
	  
	}

	// back to main dir
	main->cd();

	// japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	
	// unlock
	japp->RootUnLock();

	return NOERROR;

}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_point_calib::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	Run_Number = runnumber;
	//BCAL_Neutrals->Fill();
	//cout << " run number = " << RunNumber << endl;
	/*
	//Optional: Retrieve REST writer for writing out skims
	eventLoop->GetSingle(dEventWriterREST);
	*/

	//vector<const DTrackFinder *> finders;
	//eventLoop->Get(finders);
	//finder = const_cast<DTrackFinder*>(finders[0]);

	/*
	//Recommeded: Create output ROOT TTrees (nothing is done if already created)
	eventLoop->GetSingle(dEventWriterROOT);
	dEventWriterROOT->Create_DataTrees(eventLoop);
	*/

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_point_calib::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootFillLock(this);
	//  ... fill historgrams or trees ...
	// japp->RootFillUnLock(this);
	
	vector<const DBCALShower*> locBCALShowers;
	vector<const DBCALHit*> bcalhits;
	vector<const DBCALPoint*> locBCALPoints;

	loop->Get(locBCALShowers);
	loop->Get(bcalhits);
	loop->Get(locBCALPoints);

	vector<const DTrackTimeBased*> locTrackTimeBased;
	loop->Get(locTrackTimeBased);

	vector <const DBCALShower *> matchedShowers;
	vector < const DBCALShower *> matchedShowersneg;
	vector < const DBCALShower *> matchedShowerspos;
	vector <const DTrackTimeBased* > matchedTracks;
	DVector3 mypos(0.0,0.0,0.0);

	map< const DBCALShower*, vector<const DBCALPoint*> > matchedShowerPoints_cache;

	// the following two loops loop through a) the showers in each event and b) the tracks in each event and determine the z and phi of both. If the quantities dZ and dPhi are "small" enough, the relevant track and shower are appended to the end of the matchedTracks and matchedShowers vectors respectively

	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
	  for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	
	        double x = locBCALShowers[j]->x;
	        double y = locBCALShowers[j]->y;
		double z = locBCALShowers[j]->z;

		DVector3 pos_bcal(x,y,z);
		double R = pos_bcal.Perp();
		double phi = pos_bcal.Phi();

		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos); // takes R and returns mypos

		double dPhi = TMath::Abs(mypos.Phi()-phi);
		double dZ = TMath::Abs(mypos.Z() - z);

                // save showers matched to tracks
		if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) {
			matchedShowers.push_back(locBCALShowers[j]);
			matchedTracks.push_back(locTrackTimeBased[i]);
		}

	  } // end of shower loop
	} // end of track loop



	event_count++;
	if (VERBOSE && (event_count%100 == 0)) printf ("Event count=%d, EventNumber=%lu\n",event_count,eventnumber);

	// vectors to store the intersection of the matched tracks with the middle of each BCAL layer
	DVector3 mypos_1(0.0,0.0,0.0); // middle of Layer 1
	DVector3 mypos_2(0.0,0.0,0.0); // middle of Layer 2
	DVector3 mypos_3(0.0,0.0,0.0); // middle of Layer 3
	DVector3 mypos_4(0.0,0.0,0.0); // middle of Layer 4
	
	// get the inner radius of each layer (and the outer radius of layer 4)
	float* bcal_radii;
	bcal_radii = DBCALGeometry::GetBCAL_radii();
	
	// the radii at the middle of each layer
	double R1 = 0.5*(bcal_radii[0] + bcal_radii[1]); 
	double R2 = 0.5*(bcal_radii[1] + bcal_radii[2]); 
	double R3 = 0.5*(bcal_radii[2] + bcal_radii[3]); 
	double R4 = 0.5*(bcal_radii[3] + bcal_radii[4]); 

        // loop over matched showers
        Int_t numshowers_per_event = matchedShowers.size();
	if (numshowers_per_event > 0) 
	  h1_Num_matched_showers->Fill(numshowers_per_event);
        Int_t numtracks_per_event = matchedTracks.size();
        if (numtracks_per_event > 0) 
	  h1trk_Num_matched_tracks->Fill(numtracks_per_event);

	// cache points for each matched points to speed up the locked loop below
	for(int i=0; i<numshowers_per_event; i++)  {
		vector<const DBCALPoint*> points;
		matchedShowers[i]->Get(points);
		matchedShowerPoints_cache[ matchedShowers[i] ] = points;
	}

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	for(int i=0; i<numshowers_per_event; i++){
	        if(DEBUG) {
		  // fill histogram information for matched showers
		  Float_t E = matchedShowers[i]->E;
		  Float_t E_raw = matchedShowers[i]->E_raw;
		  Float_t x = matchedShowers[i]->x;
		  Float_t y = matchedShowers[i]->y;
		  Float_t z = matchedShowers[i]->z;
		  Float_t t = matchedShowers[i]->t;
		  Int_t N_cell = matchedShowers[i]->N_cell;
		  
		  h1_E->Fill(E);
		  h1_E_raw->Fill(E_raw);
		  h1_x->Fill(x);
		  h1_y->Fill(y);
		  h1_z->Fill(z);
		  h1_t->Fill(t);   
		  h1_N_cell->Fill(N_cell);

		  // fill histograms related to matched track
		  double FOM = matchedTracks[i]->FOM;
		  double pmag = matchedTracks[i]->pmag();
		  double px = matchedTracks[i]->px();
		  double py = matchedTracks[i]->py();
		  double pz = matchedTracks[i]->pz();
		  double x_track = matchedTracks[i]->x();
		  double y_track = matchedTracks[i]->y();
		  double z_track = matchedTracks[i]->z();

		  double energy = matchedTracks[i]->energy();
		  h1trk_FOM->Fill(FOM);
		  h1trk_pmag->Fill(pmag);
		  h1trk_px->Fill(px);
		  h1trk_py->Fill(py);
		  h1trk_pz->Fill(pz);
		  h1trk_x->Fill(x_track);
		  h1trk_y->Fill(y_track);
		  h1trk_z->Fill(z_track);
		  h1trk_energy->Fill(energy);

		  h2_Evsenergy->Fill(energy,E);
		  h2_pmagvsenergy->Fill(energy,pmag);
		}

		// get mypos for given radius R
		matchedTracks[i]->rt->GetIntersectionWithRadius(R1, mypos_1); // middle of layer 1
		matchedTracks[i]->rt->GetIntersectionWithRadius(R2, mypos_2); // middle of layer 2
		matchedTracks[i]->rt->GetIntersectionWithRadius(R3, mypos_3); // middle of layer 3
		matchedTracks[i]->rt->GetIntersectionWithRadius(R4, mypos_4); // middle of layer 4

		// create an array with the values of z_of_track for each layer
		double mypos_z_array[4] = {};
		mypos_z_array[0] = mypos_1.Z(); // Layer 1
		mypos_z_array[1] = mypos_2.Z(); // Layer 2
		mypos_z_array[2] = mypos_3.Z(); // Layer 3 
		mypos_z_array[3] = mypos_4.Z(); // Layer 4

		double target_center = 65.0;

		// get associated information in the shower

		//vector<const DBCALPoint*> points;
		//matchedShowers[i]->Get(points);
		vector<const DBCALPoint*> &points = matchedShowerPoints_cache[ matchedShowers[i] ];

		// loop over points in shower 
		Int_t numpoints_per_shower = points.size();
		h1pt_Num_points->Fill(numpoints_per_shower);

		for (int jj=0; jj<numpoints_per_shower; jj++) {
			int module = points[jj]->module(); 
			int layer = points[jj]->layer();
			int sector = points[jj]->sector();
			int cell_id = (module-1)*16 + (layer-1)*4 + sector;
			double energy = points[jj]->E();
			double energy_US = points[jj]->E_US();
			double energy_DS = points[jj]->E_DS();
			double z_point_wrt_target_center = points[jj]->z();
			double z_point = z_point_wrt_target_center + target_center;
			double r_point = points[jj]->r();
			double phi_point = points[jj]->phi();
			double t_point = points[jj]->t();


			if(DEBUG) { 
			  // fill histograms with point information
			  h1pt_module->Fill(module);
			  h1pt_layer->Fill(layer);
			  h1pt_sector->Fill(sector);
			  h1pt_cell_id->Fill(cell_id);
			  h1pt_energy->Fill(energy);
			  h1pt_energy_US->Fill(energy_US);
			  h1pt_energy_DS->Fill(energy_DS);
			  h1pt_z->Fill(z_point);  
			  h1pt_r->Fill(r_point);
			  h1pt_phi->Fill(phi_point);
			  h1pt_t->Fill(t_point);
			}

			// accept only positive z-coordinates
			if(mypos_z_array[layer-1] > 0 && z_point > 0){
      			        if(DEBUG) { 
				  // fill histograms for z_track and z_point correlation
				  h1_ztrack_minus_zpoint->Fill(mypos_z_array[layer-1] - z_point);
				  h1_abs_ztrack_minus_zpoint->Fill(TMath::Abs(mypos_z_array[layer-1] - z_point));
				  h2_ztrack_minus_zpoint_vs_ztrack->Fill(mypos_z_array[layer-1], mypos_z_array[layer-1] - z_point);
				  
				  h2_zpoint_vs_ztrack->Fill(mypos_z_array[layer-1],z_point);
				}

				// introduce |z_track - z_point| cuts  here, if needed
				if(TMath::Abs(mypos_z_array[layer-1] - z_point) < 400.0){
					// fill z_track vs z_point graphs
					z_track_graph[cell_id-1].push_back(mypos_z_array[layer-1]);
					z_point_graph[cell_id-1].push_back(z_point);
				}
				else
				{
					// fill graphs of rejected points, if cuts were applied
					z_track_thrown[cell_id-1].push_back(mypos_z_array[layer-1]);
					z_point_thrown[cell_id-1].push_back(z_point);
					h2_zpoint_vs_ztrack_thrown->Fill(mypos_z_array[layer-1],z_point);
				} 
			}

		}  // end loop over points in shower

	}   // end loop over matched showers

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	/*
	//Optional: Save event to output REST file. Use this to create skims.
	dEventWriterREST->Write_RESTEvent(loop, "BCAL_Shower"); //string is part of output file name
	*/


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_point_calib::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_point_calib::fini(void)
{
	// Called before program exit after event processing is finished.
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	// lock all root operations
	// japp->RootWriteLock();

	int channels = 768;

	// array to store the percentage of thrown points per channel
	double percent[768] = {0};
	double number_of_points = 0.0;
	double number_of_thrown_points = 0.0;

	// loop over channels (cells) and fill graphs
	for(int m = 0; m < channels; ++m){
		h2_tgraph[m] = new TGraph(z_track_graph[m].size(), &(z_track_graph[m][0]), &(z_point_graph[m][0]));
		h2_tgraph[m]->SetTitle(Form("z of point vs z of track for channel %i", m+1));
		h2_tgraph[m]->GetXaxis()->SetTitle("z of track (cm)");
		h2_tgraph[m]->GetYaxis()->SetTitle("z of point (cm)");
		h2_tgraph[m]->GetXaxis()->SetLimits(0.0, 500.0);
		h2_tgraph[m]->GetYaxis()->SetRangeUser(0.0, 500.0);
		h2_tgraph[m]->SetMarkerStyle(20);

		h2_tgraph[m]->Draw("A");
		h2_tgraph[m]->Write(Form("h2_tgraph[%i]", m+1));

		// graphs of the thrown points, if cuts were applied
		h2_tgraph_thrown[m] = new TGraph(z_track_thrown[m].size(), &(z_track_thrown[m][0]), &(z_point_thrown[m][0]));
		h2_tgraph_thrown[m]->SetTitle(Form("Thrown z of point vs z of track for channel %i", m+1));
		h2_tgraph_thrown[m]->GetXaxis()->SetTitle("z of track (cm)");
		h2_tgraph_thrown[m]->GetYaxis()->SetTitle("z of point (cm)");
		h2_tgraph_thrown[m]->GetXaxis()->SetLimits(0.0, 500.0);
		h2_tgraph_thrown[m]->GetYaxis()->SetRangeUser(0.0, 500.0);
		h2_tgraph_thrown[m]->SetMarkerStyle(20);
		h2_tgraph_thrown[m]->Draw("A");
		h2_tgraph_thrown[m]->Write(Form("h2_tgraph_thrown[%i]", m+1));

		// fill histos with percentages of thrown points per channel
		number_of_points = h2_tgraph[m]->GetN();
		number_of_thrown_points = h2_tgraph_thrown[m]->GetN();
		percent[m] = 100 * (number_of_thrown_points / number_of_points);
		if(percent[m] >= 0.0 && percent[m] <= 100.0){
			h1_thrown_per_channel->Fill(m+1, percent[m]);
		}
	} 

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	
 	// unlock
 	//japp->RootUnLock();

	return NOERROR;
}

