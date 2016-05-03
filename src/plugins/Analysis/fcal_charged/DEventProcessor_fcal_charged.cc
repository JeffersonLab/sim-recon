// $Id$
//
//    File: DEventProcessor_fcal_charged.cc  
// Modified version of BCAL_Eff
//

#include "DEventProcessor_fcal_charged.h"

#include <TLorentzVector.h>
#include "TMath.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALCluster.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"
#include "GlueX.h"

#include "BCAL/DBCALGeometry.h"
#include "FCAL/DFCALGeometry.h"

#include "TOF/DTOFPoint.h"
#include "DVector3.h"
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

DFCALGeometry fcalgeom;


static TH1I* h1_stats = NULL;
static TH1I* h1_deltaX = NULL;
static TH1I* h1_deltaY = NULL;
static TH1I* h1_deltaX_tof = NULL;
static TH1I* h1_deltaY_tof = NULL;
static TH1I* h1_Efcal = NULL;
static TH1I* h1_Ebcal = NULL;
static TH1I* h1_tfcal = NULL;
static TH1I* h1_intOverPeak = NULL;

static TH1I* h1_N9 = NULL;
static TH1I* h1_E9 = NULL;
static TH1I* h1_t9 = NULL;
static TH1I* h1_t9sigma = NULL;
static TH1I* h1_N1 = NULL;
static TH1I* h1_E1 = NULL;
static TH1I* h1_t1 = NULL;
static TH1I* h1_t1sigma = NULL;

static TH2I* h2_dEdx_vsp = NULL;
static TH2I* h2_dEdx9_vsp = NULL;
static TH2I* h2_YvsX9 = NULL;
static TH2I* h2_YvsXcheck = NULL;
static TH2I* h2_YvsX1 = NULL;
static TH2I* h2_E9vsp = NULL;
static TH2I* h2_E1vsp = NULL;
static TH2I* h2_E2vsp = NULL;
static TH2I* h2_E1vspPos = NULL;
static TH2I* h2_dEdxvsE9 = NULL;
static TH2I* h2_intOverPeak_vsEfcal = NULL;
static TH2I* h2_E1_vsintOverPeak = NULL;
static TH2I* h2_E1peak_vsp = NULL;
static TH2I* h2_E1vsRlocal = NULL;
static TH2I* h2_E9_vsTheta = NULL;
static TH2I* h2_E9_vsPhi = NULL;
static TH2I* h2_E1_vsTheta = NULL;
static TH2I* h2_E1_vsPhi = NULL;
static TH2D* h2_E9_vsThetaPhiw = NULL;
static TH2D* h2_E1_vsThetaPhiw  = NULL;
static TH2D* h2_E1_vsThetaPhi  = NULL;




// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_fcal_charged()); //register this plugin
	}
} // "C"
   
//------------------
// init
//------------------
jerror_t DEventProcessor_fcal_charged::init(void)
{
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	if(h1_deltaX != NULL){
	  printf ("TEST>> DEventProcessor_fcal_charged::init - Other threads continue\n");
		return NOERROR;
	}

	
	//  ... create historgrams or trees ...

	 //	TDirectory *dir = new TDirectoryFile("FCAL","FCAL");
	 //	dir->cd();
	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("fcal_charged")->cd();

       
	 // double pi0_mass = 0.1349766;
        int const nbins=100;

        h1_stats = new TH1I("h1_stats", "Run Statistics", 20,0,20);

	h1_deltaX = new TH1I("h1_deltaX", "Fcal X - Track X", nbins,-25,25);
	h1_deltaX->SetXTitle("Fcal - Track X (cm)");
	h1_deltaX->SetYTitle("counts");
	h1_deltaY = new TH1I("h1_deltaY", "Fcal Y - Track X", nbins,-25,25);
	h1_deltaY->SetXTitle("Fcal - Track Y (cm)");
	h1_deltaY->SetYTitle("counts");

	h1_deltaX_tof = new TH1I("h1_deltaX_tof", "TOF X - Track X E1", nbins,-25,25);
	h1_deltaX_tof->SetXTitle("TOF X - Track X (cm)");
	h1_deltaX_tof->SetYTitle("counts");
	h1_deltaY_tof = new TH1I("h1_deltaY_tof", "TOF Y - Track Y E1", nbins,-25,25);
	h1_deltaY_tof->SetXTitle("TOF Y - Track Y (cm)");
	h1_deltaY_tof->SetYTitle("counts");

	h1_Ebcal = new TH1I("h1_Ebcal", "Sum of point energy", nbins,0,1);
	h1_Ebcal->SetXTitle("Bcal point energy (GeV)");
	h1_Ebcal->SetYTitle("counts");

	h1_Efcal = new TH1I("h1_Efcal", "Hit energy", nbins,0,4);
	h1_Efcal->SetXTitle("Fcal hit energy (GeV)");
	h1_Efcal->SetYTitle("counts");
	h1_tfcal = new TH1I("h1_tfcal", "Hit time", 250,-25,225);
	h1_tfcal->SetXTitle("Fcal hit time (ns)");
	h1_tfcal->SetYTitle("counts");
	h1_intOverPeak = new TH1I("h1_intOverPeak", "Hit intOverPeak",nbins,0,25);
	h1_intOverPeak->SetXTitle("Fcal hit intOverPeak");
	h1_intOverPeak->SetYTitle("counts");

	h2_dEdx_vsp = new TH2I("h2_dEdx_vsp", "Track dEdx vs p",nbins,0,4,nbins,0,10);
	h2_dEdx_vsp->SetXTitle("p (GeV/c)");
	h2_dEdx_vsp->SetYTitle("dEdx (keV/cm)");
	h2_dEdx9_vsp = new TH2I("h2_dEdx9_vsp", "Track dEdx9 vs p",nbins,0,4,nbins,0,10);
	h2_dEdx9_vsp->SetXTitle("p (GeV/c)");
	h2_dEdx9_vsp->SetYTitle("dEdx (keV/cm)");

	h1_N9 = new TH1I("h1_N9", "Hit N9",25,0,25);
	h1_N9->SetXTitle("Number of Hits N9");
	h1_N9->SetYTitle("counts");
	h1_E9 = new TH1I("h1_E9", "Energy E9",nbins,0,2);
	h1_E9->SetXTitle("Energy E9 (GeV)");
	h1_E9->SetYTitle("counts");
	h1_t9 = new TH1I("h1_t9", "Time t9",250,-25,225);
	h1_t9->SetXTitle("Time t9 (ns)");
	h1_t9->SetYTitle("counts");
	h1_t9sigma = new TH1I("h1_t9sigma", "Time t9sigma",nbins,0,10);
	h1_t9sigma->SetXTitle("Time spread t9sigma (ns)");
	h1_t9sigma->SetYTitle("counts");

	h2_YvsX9 = new TH2I("h2_YvsX9", "Hit position Y vs X, E9",240,-120,120,240,-120,120);
	h2_YvsX9->SetXTitle("X (cm)");
	h2_YvsX9->SetYTitle("Y (cm)");
	h2_YvsXcheck = new TH2I("h2_YvsXcheck", "Delta Hit Y vs X Checkered",nbins,-5,5,nbins,-5,5);
	h2_YvsXcheck->SetXTitle("X (cm)");
	h2_YvsXcheck->SetYTitle("Y (cm)");
	h2_E1vsRlocal = new TH2I("h2_E1vsRlocal", "E1 vs Rtrk rel to Block center (cm)",nbins,0,5,nbins,0,4);
	h2_E1vsRlocal->SetXTitle("Rlocal (cm)");
	h2_E1vsRlocal->SetYTitle("E1 (GeV)");
	h2_E9vsp = new TH2I("h2_E9vsp", "E9 vs p",nbins,0,4,nbins,0,4);
	h2_E9vsp->SetXTitle("p (GeV)");
	h2_E9vsp->SetYTitle("E9 (GeV)");
	h2_dEdxvsE9 = new TH2I("h2_dEdxvsE9", "dEdx vs E9 energy",nbins,0,4,nbins,0,4);
	h2_dEdxvsE9->SetXTitle("E9 (GeV)");
	h2_dEdxvsE9->SetYTitle("dEdx (keV/cm)");

	h1_N1 = new TH1I("h1_N1", "Hit N1",25,0,25);
	h1_N1->SetXTitle("Number of Hits N1");
	h1_N1->SetYTitle("counts");
	h1_E1 = new TH1I("h1_E1", "Energy E1",nbins,0,2);
	h1_E1->SetXTitle("Energy E1 (GeV)");
	h1_E1->SetYTitle("counts");
	h1_t1 = new TH1I("h1_t1", "Time t1",250,-25,225);
	h1_t1->SetXTitle("Time t1 (ns)");
	h1_t1->SetYTitle("counts");
	h1_t1sigma = new TH1I("h1_t1sigma", "Time t1sigma",nbins,0,10);
	h1_t1sigma->SetXTitle("Time spread t1sigma (ns)");
	h1_t1sigma->SetYTitle("counts");

	h2_YvsX1 = new TH2I("h2_YvsX1", "Hit position Y vs X, E1",240,-120,120,240,-120,120);
	h2_YvsX1->SetXTitle("X (cm)");
	h2_YvsX1->SetYTitle("Y (cm)");
	h2_E1vsp = new TH2I("h2_E1vsp", "E1 vs p",nbins,0,4,nbins,0,4);
	h2_E1vsp->SetXTitle("p (GeV)");
	h2_E1vsp->SetYTitle("E1 (GeV)");
	h2_E2vsp = new TH2I("h2_E2vsp", "E2 vs p",nbins,0,4,nbins,0,4);
	h2_E2vsp->SetXTitle("p (GeV)");
	h2_E2vsp->SetYTitle("E2 (GeV)");
	h2_E1vspPos = new TH2I("h2_E1vspPos", "E1 vs p Q>0",nbins,0,4,nbins,0,4);
	h2_E1vspPos->SetXTitle("p (GeV)");
	h2_E1vspPos->SetYTitle("E1 (GeV)");
	h2_E1peak_vsp = new TH2I("h2_E1peak_vsp", "E1peak*6 vs p",nbins,0,4,nbins,0,4);
	h2_E1peak_vsp->SetXTitle("p (GeV)");
	h2_E1peak_vsp->SetYTitle("E1 peak*6(GeV)");
	h2_intOverPeak_vsEfcal = new TH2I("h2_intOverPeak_vsEfcal", "intOverPeak vs Efcal",nbins,0,1,nbins,0,25);
	h2_intOverPeak_vsEfcal->SetXTitle("Efcal (GeV)");
	h2_intOverPeak_vsEfcal->SetYTitle("IntOverPeak");
	h2_E1_vsintOverPeak = new TH2I("h2_E1_vsintOverPeak", "E1 vs intOverPeak",nbins,0,25,nbins,0,4);
	h2_E1_vsintOverPeak->SetXTitle("IntOverPeak");
	h2_E1_vsintOverPeak->SetYTitle("E1 (GeV)");

	h2_E9_vsTheta = new TH2I("h2_E9_vsTheta", "E9 vs Theta",90,0,30,nbins,0,4);
	h2_E9_vsTheta->SetXTitle("Theta (deg)");
	h2_E9_vsTheta->SetYTitle("E9 (GeV)");
	h2_E1_vsTheta = new TH2I("h2_E1_vsTheta", "E1 vs Theta",90,0,30,nbins,0,4);
	h2_E1_vsTheta->SetXTitle("Theta (deg)");
	h2_E1_vsTheta->SetYTitle("E1 (GeV)");
	h2_E9_vsPhi = new TH2I("h2_E9_vsPhi", "E9 vs Phi",90,-180,180,nbins,0,4);
	h2_E9_vsPhi->SetXTitle("Phi (deg)");
	h2_E9_vsPhi->SetYTitle("E9 (GeV)");
	h2_E1_vsPhi = new TH2I("h2_E1_vsPhi", "E1 vs Phi",90,-180,180,nbins,0,4);
	h2_E1_vsPhi->SetXTitle("Phi (deg)");
	h2_E1_vsPhi->SetYTitle("E1 (GeV)");

	h2_E9_vsThetaPhiw = new TH2D("h2_E9_vsThetaPhiw", "E9 vs ThetaPhiw",90,-180,180,90,0,30);
	h2_E9_vsThetaPhiw->SetXTitle("Phi (deg)");
	h2_E9_vsThetaPhiw->SetYTitle("Theta (deg)");
	h2_E1_vsThetaPhi = new TH2D("h2_E1_vsThetaPhi", "Unity vs ThetaPhi",90,-180,180,90,0,30);
	h2_E1_vsThetaPhi->SetXTitle("Phi (deg)");
	h2_E1_vsThetaPhi->SetYTitle("Theta (deg)");
	h2_E1_vsThetaPhiw = new TH2D("h2_E1_vsThetaPhiw", "E1 vs ThetaPhiw",90,-180,180,90,0,30);
	h2_E1_vsThetaPhiw->SetXTitle("Phi (deg)");
	h2_E1_vsThetaPhiw->SetYTitle("Theta (deg)");

	// back to main dir
	printf ("TEST>> DEventProcessor_fcal_charged::init - First thread created histograms\n");
	main->cd();

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t DEventProcessor_fcal_charged::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes
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


jerror_t DEventProcessor_fcal_charged::evnt(jana::JEventLoop* locEventLoop, int locEventNumber)
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

	vector<const DFCALShower*> locFCALShowers;
	vector<const DBCALPoint*> bcalpoints;
	vector<const DTOFPoint*> tofpoints;
	vector<const DFCALHit*> fcalhits;
	vector<const DFCALCluster*> locFCALClusters;
	vector<const DVertex*> kinfitVertex;
	//const DDetectorMatches* locDetectorMatches = NULL;
	//locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locFCALShowers);
	locEventLoop->Get(bcalpoints);;
	locEventLoop->Get(tofpoints);
	locEventLoop->Get(fcalhits);
	locEventLoop->Get(locFCALClusters);
	locEventLoop->Get(kinfitVertex);

	vector<const DChargedTrack*> locChargedTrack;
	locEventLoop->Get(locChargedTrack);

	DVector3 trkpos(0.0,0.0,0.0);
	DVector3 proj_mom(0.0,0.0,0.0);
	DVector3 trkpos_tof(0.0,0.0,0.0);
	DVector3 proj_mom_tof(0.0,0.0,0.0);
	Double_t zfcal=638;
	Double_t ztof=618.8;
	DVector3 fcal_origin(0.0,0.0,zfcal);
	DVector3 fcal_normal(0.0,0.0,1.0);
	DVector3 tof_origin(0.0,0.0,ztof);
	DVector3 tof_normal(0.0,0.0,1.0);

	/* Double_t Lfcal=45.;
	Double_t zfcal_back=zfcal+Lfcal;
	DVector3 trkpos_fcal_back(0.0,0.0,0.0);;
	DVector3 proj_mom_back(0.0,0.0,0.0);
	DVector3 fcal_origin_back(0.0,0.0,zfcal_back);
	DVector3 fcal_normal_back(0.0,0.0,1.0);*/

	if (locEventNumber%100 == 0) printf ("EventNumber=%d\n",locEventNumber);

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	double p;
	double dEdx;
	double pmin=1;


        for (unsigned int i=0; i < locChargedTrack.size() ; ++i){
	  const DChargedTrack *ctrk = locChargedTrack[i];
	  const DChargedTrackHypothesis *bestctrack = ctrk->Get_BestFOM();
	  vector<const DTrackTimeBased*> locTrackTimeBased;
	  bestctrack->Get(locTrackTimeBased);     // locTrackTimeBased[0] contains the best FOM track

	  double FOM = TMath::Prob(locTrackTimeBased[0]->chisq, locTrackTimeBased[0]->Ndof);

	  // get intersection with fcal front face and backface;

	  if (locTrackTimeBased[0]->rt->GetIntersectionWithPlane(fcal_origin,fcal_normal,trkpos,proj_mom,NULL,NULL,NULL,SYS_FCAL)==NOERROR){
	  // if (locTrackTimeBased[0]->rt->GetIntersectionWithPlane(fcal_origin_back,fcal_normal_back,trkpos_fcal_back,proj_mom_back,NULL,NULL,NULL)!=NOERROR)
	  //  continue;       // get track projection to back of fcal
	    h1_stats->Fill(0);
	    if (locTrackTimeBased[0]->rt->GetIntersectionWithPlane(tof_origin,tof_normal,trkpos_tof,proj_mom_tof,NULL,NULL,NULL,SYS_TOF)!=NOERROR) {
	      continue;     // get track projection to tof plane
	    }
	    h1_stats->Fill(1);
	    
	    // sum up all energy and Bcal and keep only if likely BCAL trigger

	    float bcal_energy = 0;
	    for (unsigned int j=0; j < bcalpoints.size(); ++j) {
	      bcal_energy += bcalpoints[j]->E();
		}
	    if (bcal_energy < 0.15) {
	      continue;              // require that bcal be sufficient for trigger
	    }
	    h1_Ebcal->Fill(bcal_energy);
	    h1_stats->Fill(2);
	    

	    double q = locTrackTimeBased[0]->rt->q;
	    p = locTrackTimeBased[0]->momentum().Mag();
	    double trkmass = locTrackTimeBased[0]->mass();

	    double radius = sqrt(trkpos.X()*trkpos.X() + trkpos.Y()*trkpos.Y());  // distance of track from beamline
	    dEdx = locTrackTimeBased[0]->dEdx()*1e6;  // convert to keV
	    if(!(trkmass < 0.15 && FOM>0.01 && radius>20)) {
	      continue;   // select tracks of interest
	    }
	    h1_stats->Fill(3);

	    if (! (p>pmin)) {
	      continue;      // require minimum momentum to pion analyzis
	    }
	    h1_stats->Fill (4);
	    
	    h2_dEdx_vsp->Fill(p,dEdx);
	    // use position at the back of fcal  // asume zero field out here
	    // double xfcal_back = trkpos.X() + proj_mom.X()*Lfcal/proj_mom.Z();
	    // double yfcal_back = trkpos.Y() + proj_mom.Y()*Lfcal/proj_mom.Z();
	    double trkposX = trkpos.X();
	    double trkposY = trkpos.Y();
	    int trkrow = fcalgeom.row((float)trkposY);
	    int trkcol = fcalgeom.column((float)trkposX);
	    double dX_tof = 10000;
	    double dY_tof = 10000;
	    double dR_tof = 10000;

	    // get tof matches

	    for (unsigned int j=0; j < tofpoints.size(); ++j) {
	    double tofx = 10000;
	    double tofy = 10000;
	      if (tofpoints[j]->Is_XPositionWellDefined() && tofpoints[j]->Is_YPositionWellDefined()) {
		DVector3 pos_tof = tofpoints[j]->pos;
		tofx = pos_tof.X();
		tofy = pos_tof.Y();
		printf ("Event=%d tofx=%f, tofy=%f, trkx=%f, trky=%f \n",locEventNumber,tofx,tofy,trkposX,trkposY);
	        double dX_tof1 = tofx - trkpos_tof.X();
	        double dY_tof1 = tofy - trkpos_tof.Y();
		if (sqrt(dX_tof1*dX_tof1 + dY_tof1*dY_tof1) < dR_tof) {
		  dX_tof = dX_tof1;
		  dY_tof = dY_tof1;
		  dR_tof =  sqrt(dX_tof*dX_tof + dY_tof*dY_tof);		
		}
	      }
	    }

	    printf ("\n1 Event=%d dX_tof=%f, DY_tof=%f, trkx=%f, trky=%f \n",locEventNumber,dX_tof,dY_tof,trkposX,trkposY);
	    if ( !(abs(dX_tof)<10 && abs(dY_tof)<10) )  continue;

	    h1_stats->Fill(5);
	    printf ("2 Event=%d dX_tof=%f, DY_tof=%f, trkx=%f, trky=%f \n",locEventNumber,dX_tof,dY_tof,trkposX,trkposY);
	    printf ("Event=%d trkmass=%f radius=%f p=%f dEdx=%g,trkrow=%d trkcol=%d x=%f y=%f z=%f\n",locEventNumber,trkmass,radius,p,dEdx,trkrow,trkcol,trkposX,trkposY,trkpos.Z());

	    // get components of p at FCAL
	    // double theta = proj_mom.Theta()*TVector3::RAD2DEG;
	    double theta = proj_mom.Theta()*180./M_PI;
	    double phi = proj_mom.Phi()*180./M_PI;
	    printf ("Event=%d p=%f,theta=%f, phi=%f\n",locEventNumber,p,theta,phi);


            // loop over fcal hits

	    double E9=0;  // energy, x, y of 9 blocks surrounding track
	    double E9peak=0;
	    double x9=0;
	    double y9=0;
	    double t9=0;
	    double t9sq=0;
	    double t9sigma=0;
	    int N9=0;
	    int Delta_block=2;   // =1 for 3x3, =2 for 5x5
	    int row_E1=-1000;
	    int col_E1=-1000;
	    int drow_E1=1000;
	    int dcol_E1=1000;
	    double dX_E1=-1000;
	    double dY_E1=-1000;

	    int row_offset = 0;   // offset actual row to look for randoms
	    int col_offset = 0;
	    trkrow += row_offset;
	    trkcol += col_offset;

	    for (unsigned int j=0; j < fcalhits.size(); ++j) {
	      int row = fcalhits[j]->row;
	      int col = fcalhits[j]->column;
	      double x = fcalhits[j]->x;
	      double y = fcalhits[j]->y;
	      double Efcal = fcalhits[j]->E;
	      double tfcal= fcalhits[j]->t;
	      double intOverPeak = fcalhits[j]->intOverPeak;
	      // printf ("Event=%d row=%d col=%d x=%f y=%f Efcal=%f tfcal=%f intOverPeak=%f\n",locEventNumber,row,col,x,y,Efcal,tfcal,intOverPeak);

	      // fill histograms
	      int drow = abs(row - trkrow);
	      int dcol = abs(col - trkcol);

	      h1_deltaX->Fill(x - trkposX);
	      h1_deltaY->Fill(y - trkposY);
	      h1_Efcal->Fill(Efcal);
	      h1_tfcal->Fill(tfcal);
	      h1_intOverPeak->Fill(intOverPeak);
	      h2_intOverPeak_vsEfcal->Fill(Efcal,intOverPeak);

	      // select hits 
	      if (drow<=Delta_block && dcol<=Delta_block && (tfcal>-15 && tfcal<50) && (intOverPeak>2.5 && intOverPeak<9)) {
	      // if (drow<=Delta_block && dcol<=Delta_block && (tfcal>-15 && tfcal<50)) {
		  E9 += Efcal;
		  E9peak += Efcal*6/intOverPeak;   // factor of 6 so that E9peak ~ E9
		  x9 += x;
		  y9 += y;
		  t9 += tfcal;
		  t9sq += tfcal*tfcal;
		  N9 += 1;

		  //save for later
		  row_E1 = row;
		  col_E1 = col;
		  drow_E1 = drow;
		  dcol_E1 = dcol;
		  dX_E1 = x - trkposX;    
		  dY_E1 = y - trkposY;
		  printf ("Event=%d, trkposX=%f, trkposY=%f, dX_E1=%f, dY_E1=%f\n",locEventNumber,trkposX,trkposY,dX_E1,dY_E1);
	      }

	    } // end loop over fcal hits

	    x9 = N9>0? x9/N9 : 0;
	    y9 = N9>0? y9/N9 : 0;
	    t9 = N9>0? t9/N9 : 0;
	    t9sigma = N9>0? sqrt(t9sq/N9 - t9*t9): 0;
	    // printf ("\nEvent=%d N9=%d E9=%f x9=%f y9=%f t9=%f t9sigma=%f\n",locEventNumber,N9,E9,x9,y9,t9,t9sigma);


	    h1_stats->Fill(6);
	    if (E9 > 0.8 && theta<4) continue;
	    h1_stats->Fill(7);
	    if (N9>0) {
	      h1_N9->Fill(N9);
	      h1_E9->Fill(E9);
	      h1_t9->Fill(t9);
	      h1_t9sigma->Fill(t9sigma);
	      h2_YvsX9->Fill(trkposX,trkposY);
	      h2_dEdx9_vsp->Fill(p,dEdx);
	      h2_E9vsp->Fill(p,E9);
	      h2_dEdxvsE9->Fill(E9,dEdx);
	      h1_stats->Fill(8);
	      h2_E9_vsThetaPhiw->Fill(phi,theta,E9);
	      h2_E9_vsTheta->Fill(theta,E9);
	      h2_E9_vsPhi->Fill(phi,E9);
	    }
	    if (N9==1) {
	      h1_N1->Fill(N9);
	      h1_E1->Fill(E9);
	      h1_t1->Fill(t9);
	      h1_t1sigma->Fill(t9sigma);
	      h2_YvsX1->Fill(trkposX,trkposY);
	      h2_E1vsp->Fill(p,E9);
	      if (q > 0) h2_E1vspPos->Fill(p,E9);
	      h2_E1peak_vsp->Fill(p,E9peak);
	      h2_E1_vsintOverPeak->Fill(E9*6/E9peak,E9);   // note that this only works because a single cell fires
	      h1_deltaX_tof->Fill(dX_tof);
	      h1_deltaY_tof->Fill(dY_tof);
	      h1_stats->Fill(9);
	      h2_E1_vsThetaPhi->Fill(phi,theta);
	      h2_E1_vsThetaPhiw->Fill(phi,theta,E9);
	      h2_E1_vsTheta->Fill(theta,E9);
	      h2_E1_vsPhi->Fill(phi,E9);
	      double Rlocal = sqrt(dX_E1*dX_E1 + dY_E1*dY_E1);    
	      h2_E1vsRlocal->Fill(Rlocal,E9);
	      // if (((row_E1%2==0 && col_E1%2==0) || (row_E1%2==1 && col_E1%2==1)) ) h2_YvsXcheck->Fill(dX_E1,dY_E1);
	      h2_YvsXcheck->Fill(dX_E1,dY_E1);
	    }
	    if (N9==1 || N9==2) {
	      h2_E2vsp->Fill(p,E9);
	      h1_stats->Fill(10);
	    }



	  }   // end track intersection if statement
	}

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


	/*
	//Optional: Save event to output REST file. Use this to create skims.
	dEventWriterREST->Write_RESTEvent(locEventLoop, "FCAL_Shower"); //string is part of output file name
	*/

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_fcal_charged::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.


	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_fcal_charged::fini(void)
{
	// Called before program exit after event processing is finished.  

	return NOERROR;
}

