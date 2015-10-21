// $Id$
//
//    File: JEventProcessor_BCAL_Eff.cc  Copied from Will McGinley
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_Eff.h"

#include <TLorentzVector.h>
#include "TMath.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>



		static TH1I* h1_Num_matched_showers = NULL;
                /*static TH1I* h1_E = NULL;
		static TH1I* h1_E_raw = NULL;
		static TH1I* h1_x = NULL;
		static TH1I* h1_y = NULL;
		static TH1I* h1_z = NULL;
		static TH1I* h1_t = NULL;
		static TH1I* h1_xErr = NULL;
		static TH1I* h1_yErr = NULL;
		static TH1I* h1_zErr = NULL;
		static TH1I* h1_tErr = NULL;
		static TH1I* h1_N_cell = NULL;*/

                static TH1I* h1trk_Num_matched_tracks = NULL;
		static TH1I* h1trk_FOM = NULL;
		static TH1I* h1trk_pmag = NULL;
		static TH1I* h1trk_px = NULL;
		static TH1I* h1trk_py = NULL;
		static TH1I* h1trk_pz = NULL;
		static TH1I* h1trk_energy = NULL;

		static TH2I* h2_Evsenergy = NULL;
		static TH2I* h2_pmagvsenergy = NULL;

		static TH1I* h1clust_Num_matched_clusters = NULL;
		static TH1I* h1clust_nCells = NULL;

		static TH1I* h1pt_Num_points = NULL;
		static TH1I* h1pt_module = NULL;
		static TH1I* h1pt_layer = NULL;
		static TH1I* h1pt_sector = NULL;
		static TH1I* h1pt_cell_id = NULL;
		static TH1I* h1pt_energy = NULL;
		static TH1I* h1pt_energy_US = NULL;
		static TH1I* h1pt_energy_DS = NULL;

		static TH1I* h1eff_layertot = NULL;
		static TH1I* h1eff_layer = NULL;
		static TH1F* h1eff_eff = NULL;

		static TH1I* h1eff_cellidtot = NULL;
		static TH1I* h1eff_cellid = NULL;
		static TH1F* h1eff_cellideff = NULL;

		static TH1I* h1eff2_layertot = NULL;
		static TH1I* h1eff2_layer = NULL;
		static TH1F* h1eff2_eff2 = NULL;

		static TH1I* h1eff2_cellidtot = NULL;
		static TH1I* h1eff2_cellid = NULL;
		static TH1F* h1eff2_cellideff2 = NULL;

		static TH1I* h1E_US_layer1 = NULL;
		static TH1I* h1E_DS_layer1 = NULL;
		static TH2I* h2E_USvsDS_layer1 = NULL;

		int event_count=0;


// Routine used to create our JEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_BCAL_Eff()); //register this plugin
	}
} // "C"
   
//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_Eff::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	
	 japp->RootWriteLock();

	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	if(h1_Num_matched_showers != NULL){
		japp->RootUnLock();
		return NOERROR;
	}
	
	//  ... create historgrams or trees ...

	 //	TDirectory *dir = new TDirectoryFile("BCAL","BCAL");
	 //	dir->cd();
	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcal_eff")->cd();

       
	 // double pi0_mass = 0.1349766;
        int nbins=100;

	h1_Num_matched_showers = new TH1I("h1_Num_matched_showers", "Number of matched showers per event", 20,0,20);
	h1_Num_matched_showers->SetXTitle("Number of showers");
	h1_Num_matched_showers->SetYTitle("counts");

	/*h1_E = new TH1I("h1_E", "Matched energy per shower",nbins,0,2);
	h1_E->SetXTitle("Energy per shower (GeV)");
	h1_E->SetYTitle("counts");
	h1_E_raw = new TH1I("h1_E_raw", "Matched raw energy per shower",nbins,0,2);
	h1_E_raw->SetXTitle("Raw Energy per shower (GeV)");
	h1_E_raw->SetYTitle("counts");
	h1_x = new TH1I("h1_x", "x pos of shower",nbins,-100,100);
	h1_x->SetXTitle("x (cm)");
	h1_x->SetYTitle("counts");
	h1_y = new TH1I("h1_y", "y pos of shower",nbins,-100,100);
	h1_y->SetXTitle("y (cm)");
	h1_y->SetYTitle("counts");
	h1_z = new TH1I("h1_z", "z pos of shower",nbins,0,500);
	h1_z->SetXTitle("z (cm)");
	h1_z->SetYTitle("counts");
	h1_t = new TH1I("h1_t", "time of shower",nbins,-100,100);
	h1_t->SetXTitle("time of shower (ns)");
	h1_t->SetYTitle("counts");
	h1_xErr = new TH1I("h1_xErr", "x err of shower",nbins,0,20);
	h1_xErr->SetXTitle("x err (cm)");
	h1_xErr->SetYTitle("counts");
	h1_yErr = new TH1I("h1_yErr", "y err of shower",nbins,0,20);
	h1_yErr->SetXTitle("y err (cm)");
	h1_yErr->SetYTitle("counts");
	h1_zErr = new TH1I("h1_zErr", "z err of shower",nbins,0,20);
	h1_zErr->SetXTitle("z err (cm)");
	h1_zErr->SetYTitle("counts");
	h1_tErr = new TH1I("h1_tErr", "t err of shower",nbins,0,5);
	h1_tErr->SetXTitle("t err (ns)");
	h1_tErr->SetYTitle("counts");
	h1_N_cell = new TH1I("h1_N_cell", "N_cell of shower",20,0,20);
	h1_N_cell->SetXTitle("N_cell");
	h1_N_cell->SetYTitle("counts");*/

	h1trk_Num_matched_tracks = new TH1I("h1trk_Num_matched_tracks", "Number of matched tracks per event", 20,0,20);
	h1trk_Num_matched_tracks->SetXTitle("Number of tracks");
	h1trk_Num_matched_tracks->SetYTitle("counts");
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
	h1trk_energy = new TH1I("h1trk_energy", "energy for matched tracks", nbins,0,4);
	h1trk_energy->SetXTitle("energy (GeV)");
	h1trk_energy->SetYTitle("counts");

	h2_Evsenergy = new TH2I("h2_Evsenergy", "E vs energy matched tracks", nbins,0,4,nbins,0,4);
	h2_Evsenergy->SetXTitle("track energy (GeV)");
	h2_Evsenergy->SetYTitle("Eshower (GeV)");
	h2_pmagvsenergy = new TH2I("h2_pmagvsenergy", "pmag vs energy matched tracks", nbins,0,4,nbins,0,4);
	h2_pmagvsenergy->SetXTitle("track energy (GeV)");
	h2_pmagvsenergy->SetYTitle("track pmag (GeV)");

	h1clust_Num_matched_clusters = new TH1I("h1_Num_matched_clusters", "Number of matched clusters per event", 20,0,20);
	h1clust_Num_matched_clusters->SetXTitle("Number of clusters");
	h1clust_Num_matched_clusters->SetYTitle("counts");
	h1clust_nCells = new TH1I("h1clust_nCells", "Cluster nCells",20,0,20);
	h1clust_nCells->SetXTitle("Cluster nCells");
	h1clust_nCells->SetYTitle("counts");

	h1pt_Num_points = new TH1I("h1pt_Num_points", "pt points per cluster",20,0,20);
	h1pt_Num_points->SetXTitle("Points per cluster");
	h1pt_Num_points->SetYTitle("counts");
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

	h1eff_layertot = new TH1I("h1eff_layertot", "eff layertot",5,0,5);
	h1eff_layertot->SetXTitle("Layer");
	h1eff_layertot->SetYTitle("counts");
	//h1eff_layertot->Sumw2();
	h1eff_layer = new TH1I("h1eff_layer", "eff layer hit",5,0,5);
	h1eff_layer->SetXTitle("Layer hit");
	h1eff_layer->SetYTitle("counts");
	//h1eff_layer->Sumw2();
	h1eff_eff = new TH1F("h1eff_eff", "efficiency",5,0,5);
	h1eff_eff->SetXTitle("Layer");
	h1eff_eff->SetYTitle("Efficiency");
        h1eff_eff->Sumw2();

	h1eff_cellidtot = new TH1I("h1eff_cellidtot", "eff cellidtot",200,0,200);
	h1eff_cellidtot->SetXTitle("Cellid");
	h1eff_cellidtot->SetYTitle("counts");
	//h1eff_cellidtot->Sumw2();
	h1eff_cellid = new TH1I("h1eff_cellid", "eff cellid hit",200,0,200);
	h1eff_cellid->SetXTitle("Cellid hit");
	h1eff_cellid->SetYTitle("counts");
	//h1eff_cellid->Sumw2();
	h1eff_cellideff = new TH1F("h1eff_cellideff", "efficiency",200,0,200);
	h1eff_cellideff->SetXTitle("Cellid");
	h1eff_cellideff->SetYTitle("Efficiency");
        h1eff_cellideff->Sumw2();

	h1eff2_layertot = new TH1I("h1eff2_layertot", "eff2 layertot",5,0,5);
	h1eff2_layertot->SetXTitle("Layer");
	h1eff2_layertot->SetYTitle("counts");
	//h1eff2_layertot->Sumw2();
	h1eff2_layer = new TH1I("h1eff2_layer", "eff2 layer hit",5,0,5);
	h1eff2_layer->SetXTitle("Layer hit");
	h1eff2_layer->SetYTitle("counts");
	//h1eff2_layer->Sumw2();
	h1eff2_eff2 = new TH1F("h1eff2_eff2", "eff2 efficiency",5,0,5);
	h1eff2_eff2->SetXTitle("Layer");
	h1eff2_eff2->SetYTitle("Efficiency");
	h1eff2_eff2->Sumw2();

	h1eff2_cellidtot = new TH1I("h1eff2_cellidtot", "eff2 cellidtot",200,0,200);
	h1eff2_cellidtot->SetXTitle("Cellid");
	h1eff2_cellidtot->SetYTitle("counts");
	//h1eff2_cellidtot->Sumw2();
	h1eff2_cellid = new TH1I("h1eff2_cellid", "eff2 cellid hit",200,0,200);
	h1eff2_cellid->SetXTitle("Cellid hit");
	h1eff2_cellid->SetYTitle("counts");
	//h1eff2_cellid->Sumw2();
	h1eff2_cellideff2 = new TH1F("h1eff2_cellideff2", "eff2 efficiency",200,0,200);
	h1eff2_cellideff2->SetXTitle("Cellid");
	h1eff2_cellideff2->SetYTitle("Efficiency");
	h1eff2_cellideff2->Sumw2();

	h1E_US_layer1 = new TH1I("h1E_US_layer1", "E US layer1",nbins,0,1);
	h1E_US_layer1->SetXTitle("Upstream E (GeV)");
	h1E_US_layer1->SetYTitle("Counts");
	h1E_DS_layer1 = new TH1I("h1E_DS_layer1", "E DS layer1",nbins,0,1);
	h1E_DS_layer1->SetXTitle("Downstream E (GeV)");
	h1E_DS_layer1->SetYTitle("Counts");
	h2E_USvsDS_layer1 = new TH2I("h2E_USvsDS_layer1", "E US vs DS ",nbins,0,1,nbins,0,1);
	h2E_USvsDS_layer1->SetXTitle("Downstream E (GeV)");
	h2E_USvsDS_layer1->SetYTitle("Upstream E (GeV)");	


	
	// back to main dir
	main->cd();

        japp->RootUnLock();


	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes
	Run_Number = locRunNumber;
	//BCAL_Neutrals->Fill();
	//cout << " run number = " << RunNumber << endl;
	/*
	//Optional: Retrieve REST writer for writing out skims
	locEventLoop->GetSingle(dEventWriterREST);
	*/

	//vector<const DTrackFinder *> finders;
	//locEventLoop->Get(finders);
	//finder = const_cast<DTrackFinder*>(finders[0]);

	/*
	//Recommeded: Create output ROOT TTrees (nothing is done if already created)
	locEventLoop->GetSingle(dEventWriterROOT);
	dEventWriterROOT->Create_DataTrees(locEventLoop);
	*/

	return NOERROR;
}

//------------------
// evnt
//------------------


jerror_t JEventProcessor_BCAL_Eff::evnt(jana::JEventLoop* locEventLoop, int locEventNumber)
{

	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software

	vector<const DBCALShower*> locBCALShowers;
	vector<const DBCALHit*> bcalhits;
	vector<const DBCALCluster*> locBCALClusters;
	vector<const DBCALPoint*> locBCALPoints;
	vector<const DVertex*> kinfitVertex;
	//const DDetectorMatches* locDetectorMatches = NULL;
	//locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locBCALShowers);
	locEventLoop->Get(bcalhits);
	locEventLoop->Get(locBCALClusters);
	locEventLoop->Get(locBCALPoints);
	locEventLoop->Get(kinfitVertex);

	vector<const DTrackTimeBased*> locTrackTimeBased;
	locEventLoop->Get(locTrackTimeBased);

	vector <const DBCALShower *> matchedShowers;
	vector < const DBCALShower *> matchedShowersneg;
	vector < const DBCALShower *> matchedShowerspos;
	vector <const DTrackTimeBased* > matchedTracks;
	DVector3 mypos(0.0,0.0,0.0);

	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
	  for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	
	        double x = locBCALShowers[j]->x;
	        double y = locBCALShowers[j]->y;
		double z = locBCALShowers[j]->z;
		// double E = locBCALShowers[j]->E;
		DVector3 pos_bcal(x,y,z);
		double R = pos_bcal.Perp();
		double phi = pos_bcal.Phi();
		// double FOM = TMath::Prob(locTrackTimeBased[i]->chisq, locTrackTimeBased[i]->Ndof);

		locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);

		// double q = locTrackTimeBased[i]->rt->q;
		locTrackTimeBased[i]->momentum().Mag();
		double trkmass = locTrackTimeBased[i]->mass();
		double dPhi = TMath::Abs(mypos.Phi()-phi);
		double dZ = TMath::Abs(mypos.Z() - z);

                // save showers matched to tracks
		
		  // assume program is run with plugin -PTRKFIT:MASS_HYPOTHESES_NEGATIVE=0.13957 -PTRKFIT:MASS_HYPOTHESES_POSITIVE=0.13957
		// if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R) {
		  if(dZ < 30.0 && dPhi < 0.18 && mypos.Perp() == R && trkmass < 0.15) {   // select pion hypothesis only
		  matchedShowers.push_back(locBCALShowers[j]);
	          matchedTracks.push_back(locTrackTimeBased[i]);
		  // printf ("Matched event=%d, i=%d, j=%d, p=%f, Ztrk=%f Zshr=%f, Phitrk=%f Phishw=%f\n",locEventNumber,i,j,p,mypos.Z(),z,mypos.Phi(),pos_bcal.Phi());
		}
	  }
	}

	japp->RootWriteLock();
	event_count++;
//	if (event_count%100 == 0) printf ("Event count=%d, EventNumber=%d\n",event_count,locEventNumber);


        // Get kinematic-fitted vertex

#if 0  // disabling since it doesn't do anything other than cause compiler warnings  10/15/2015 DL
	double kinfitVertexX, kinfitVertexY, kinfitVertexZ, kinfitVertexT;
	for (unsigned int i = 0 ; i < kinfitVertex.size(); i++)
	{
		kinfitVertexX = kinfitVertex[i]->dSpacetimeVertex.X();
		kinfitVertexY = kinfitVertex[i]->dSpacetimeVertex.Y();
		kinfitVertexZ = kinfitVertex[i]->dSpacetimeVertex.Z();
		kinfitVertexT = kinfitVertex[i]->dSpacetimeVertex.T();
		//		goodVertexZ->Fill(kinfitVertexZ);
	}
#endif

        // now loop over matched showers

        Int_t numshowers_per_event = matchedShowers.size();
	if (numshowers_per_event > 0) h1_Num_matched_showers->Fill(numshowers_per_event);
        Int_t numtracks_per_event = matchedTracks.size();
        if (numtracks_per_event > 0) h1trk_Num_matched_tracks->Fill(numtracks_per_event);
        // if (numtracks_per_event > 0) printf ("Event=%d nshower=%d ntracks=%d\n\n",locEventNumber,numshowers_per_event,numtracks_per_event);
	
	for(int i=0; i<numshowers_per_event; i++){

	  // fill histogram information for matched showers
          Float_t E = matchedShowers[i]->E;
          /*Float_t E_raw = matchedShowers[i]->E_raw;
          Float_t x = matchedShowers[i]->x;
          Float_t y = matchedShowers[i]->y;
          Float_t z = matchedShowers[i]->z;
          Float_t t = matchedShowers[i]->t;
          Float_t xErr = matchedShowers[i]->xErr;
          Float_t yErr = matchedShowers[i]->yErr;
          Float_t zErr = matchedShowers[i]->zErr;
          Float_t tErr = matchedShowers[i]->tErr;
          Int_t N_cell = matchedShowers[i]->N_cell;

          h1_E->Fill(E);
          h1_E_raw->Fill(E_raw);
          h1_x->Fill(x);
          h1_y->Fill(y);
          h1_z->Fill(z);
          h1_t->Fill(t);   
          h1_xErr->Fill(xErr);
          h1_yErr->Fill(yErr);
          h1_zErr->Fill(zErr);
          h1_tErr->Fill(tErr);
          h1_N_cell->Fill(N_cell);*/

	  // fill histograms related to matched track

	  double FOM = matchedTracks[i]->FOM;
          double pmag = matchedTracks[i]->pmag();
          double px = matchedTracks[i]->px();
          double py = matchedTracks[i]->py();
          double pz = matchedTracks[i]->pz();
          double energy = matchedTracks[i]->energy();
          h1trk_FOM->Fill(FOM);
          h1trk_pmag->Fill(pmag);
          h1trk_px->Fill(px);
          h1trk_py->Fill(py);
          h1trk_pz->Fill(pz);
          h1trk_energy->Fill(energy);

          h2_Evsenergy->Fill(energy,E);
          h2_pmagvsenergy->Fill(energy,pmag);

	  // get associated information in the shower

	  vector<const DBCALCluster*> clusters;
	  matchedShowers[i]->Get(clusters);

	  vector <int> cluster_modules;  // vector with list of modules  in cluster

	  // loop over associated clusters in shower
	  Int_t numclusters_per_event = clusters.size();
	  h1clust_Num_matched_clusters->Fill(numclusters_per_event);
	  for(int j=0; j<numclusters_per_event; j++) {

	    Int_t clusters_nCells = clusters[j]->nCells();
	    h1clust_nCells->Fill(clusters_nCells);

	    // get associated information in the shower

	    vector<const DBCALPoint*> points;
	    clusters[j]->Get(points);
            cluster_modules.clear();

	    // loop over points in cluster
	    Int_t numpoints_per_cluster = points.size();
	    h1pt_Num_points->Fill(numpoints_per_cluster);

	    // initialize hit variables
	    Int_t hitlayer1=0;
	    Int_t hitlayer2=0;
	    Int_t hitlayer3=0;
	    Int_t hitlayer4=0;
	    Int_t module_ndx=0;    // use one of the modules in cluster as index
	    for (int jj=0; jj<numpoints_per_cluster; jj++) {

	       int module = points[jj]->module(); 
	       int layer = points[jj]->layer();
	       int sector = points[jj]->sector();
	       // int cell_id = DBCALGeometry::cellId(module, layer, sector);
	       int cell_id = (module-1)*16 + (layer-1)*4 + sector;
	       //	cout << " cell ID = " << cell_id << " module = " << module1 << " layer = " << layer1 << " sector = " << sector1 << endl;
	       double energy = points[jj]->E();
	       double energy_US = points[jj]->E_US();
	       double energy_DS = points[jj]->E_DS();

	       // fill histograms with point information
	       h1pt_module->Fill(module);
	       h1pt_layer->Fill(layer);
	       h1pt_sector->Fill(sector);
	       h1pt_cell_id->Fill(cell_id);
	       h1pt_energy->Fill(energy);
	       h1pt_energy_US->Fill(energy_US);
	       h1pt_energy_DS->Fill(energy_DS);

	       if (layer == 1) hitlayer1 = 1;
	       if (layer == 2) hitlayer2 = 2;
	       if (layer == 3) hitlayer3 = 3;
	       if (layer == 4) hitlayer4 = 4;
	       // printf ("event=%d, pt jj=%d, module=%d, layer=%d, sector=%d, cell_id=%d, hitlayer1=%d, hitlayer2=%d, hitlayer3=%d, hitlayer4=%d\n",locEventNumber,jj,module,layer, sector,cell_id,hitlayer1,hitlayer2,hitlayer3,hitlayer4);

               // fill vector with module number if not yet in list
	       bool module_in_list = false;
	       for (unsigned int jjj=0; jjj<cluster_modules.size(); jjj++) {
	         if (module == cluster_modules[jjj]) module_in_list = true; 
	       }
	       if (!module_in_list) cluster_modules.push_back(module);
	       module_ndx = module;

            }

	    // fill efficiency histograms, once per cluster

	    if (hitlayer2 + hitlayer3 + hitlayer4 > 0) {
	      h1eff_layertot->Fill(1);
	      h1eff_layer->Fill(hitlayer1);
	      int cellid = (module_ndx-1)*4 + 1;
	      h1eff_cellidtot->Fill(cellid);
	      if (hitlayer1>0) h1eff_cellid->Fill(cellid);
	    }

	    if (hitlayer3 + hitlayer4 > 0) {
	      h1eff_layertot->Fill(2);
	      h1eff_layer->Fill(hitlayer2);
	      int cellid = (module_ndx-1)*4 + 2;
	      h1eff_cellidtot->Fill(cellid);
	      if (hitlayer2>0) h1eff_cellid->Fill(cellid);
	    }
	    if (hitlayer4 > 0) {
	      h1eff_layertot->Fill(3);
	      h1eff_layer->Fill(hitlayer3);
	      int cellid = (module_ndx-1)*4 + 3;
	      h1eff_cellidtot->Fill(cellid);
	      if (hitlayer3>0) h1eff_cellid->Fill(cellid);
	    }

	    if (hitlayer1 + hitlayer2 + hitlayer3 > 3) {
	      h1eff_layertot->Fill(4);
	      h1eff_layer->Fill(hitlayer4);
	      int cellid = (module_ndx-1)*4 + 4;
	      h1eff_cellidtot->Fill(cellid);
	    if (hitlayer4>0) h1eff_cellid->Fill(cellid);
	    }
	    
	  }

	    int hitlayer1=0;
	    int hitlayer2=0;
	    int hitlayer3=0;
	    int hitlayer4=0;
	    Int_t module_ndx=0;    // use one of the modules in cluster as index
	    // loop over all points in BCAL
	    for (unsigned int jjj=0; jjj<bcalhits.size(); jjj++) {
	       int module = bcalhits[jjj]->module; 

	       // histogram all hits, whether or not in clustersssss

	       if (bcalhits[jjj]->end == 0) h1E_US_layer1->Fill(bcalhits[jjj]->E);
	       if (bcalhits[jjj]->end == 1) h1E_DS_layer1->Fill(bcalhits[jjj]->E);

	       // check to see if module is in cluster list
	       bool module_in_list = false;
	       for (unsigned int k=0; k<cluster_modules.size(); k++) {
	         if (module == cluster_modules[k]) module_in_list = true; 
	       }
	       if (!module_in_list) continue;

	       int layer = bcalhits[jjj]->layer;
	       // int sector = bcalhits[jjj]->sector;
	       // int end = bcalhits[jjj]->end;
	       // printf ("locPts event=%d, module=%d, layer=%d sector=%d end=%d\n",locEventNumber,module,layer,sector,end);

	       if (layer == 1) hitlayer1 = 1;
	       if (layer == 2) hitlayer2 = 2;
	       if (layer == 3) hitlayer3 = 3;
	       if (layer == 4) hitlayer4 = 4;
               module_ndx = module;

	    }

	    // fill efficiency histograms, once per cluster

	    if (hitlayer2 + hitlayer3 + hitlayer4 > 0) {
	      h1eff2_layertot->Fill(1);
	      h1eff2_layer->Fill(hitlayer1);;
	      int cellid = (module_ndx-1)*4 + 1;
	      h1eff2_cellidtot->Fill(cellid);
	      if (hitlayer1>0) h1eff2_cellid->Fill(cellid);
	    }

	    if (hitlayer3 + hitlayer4 > 0) {
	      h1eff2_layertot->Fill(2);
	      h1eff2_layer->Fill(hitlayer2);;
	      int cellid = (module_ndx-1)*4 + 2;
	      h1eff2_cellidtot->Fill(cellid);
	      if (hitlayer2>0) h1eff2_cellid->Fill(cellid);
	    }
	    if (hitlayer4 > 0) {
	      h1eff2_layertot->Fill(3);
	      h1eff2_layer->Fill(hitlayer3);;
	      int cellid = (module_ndx-1)*4 + 3;
	      h1eff2_cellidtot->Fill(cellid);
	      if (hitlayer3>0) h1eff2_cellid->Fill(cellid);
	    }

	    if (hitlayer1 + hitlayer2 + hitlayer3 > 3) {
	      h1eff2_layertot->Fill(4);
	      h1eff2_layer->Fill(hitlayer4);;
	      int cellid = (module_ndx-1)*4 + 4;
	      h1eff2_cellidtot->Fill(cellid);
	      if (hitlayer4>0) h1eff2_cellid->Fill(cellid);
	    }

	  

	}   // end loop over matched showers

        //UnlockState();	
	japp->RootUnLock();


	/*
	//Optional: Save event to output REST file. Use this to create skims.
	dEventWriterREST->Write_RESTEvent(locEventLoop, "BCAL_Shower"); //string is part of output file name
	*/

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

  japp->RootWriteLock();

  // h1eff_eff->Divide(h1eff_layer,h1eff_layertot);

  japp->RootUnLock();


	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

  japp->RootWriteLock();

  h1eff_eff->Divide(h1eff_layer,h1eff_layertot,1,1,"B");
  h1eff2_eff2->Divide(h1eff2_layer,h1eff2_layertot,1,1,"B");

  h1eff_cellideff->Divide(h1eff_cellid,h1eff_cellidtot,1,1,"B");
  h1eff2_cellideff2->Divide(h1eff2_cellid,h1eff2_cellidtot,1,1,"B");


  japp->RootUnLock();
	return NOERROR;
}

