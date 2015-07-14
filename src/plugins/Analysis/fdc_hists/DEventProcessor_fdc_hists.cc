// $Id: DEventProcessor_fdc_hists.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_fdc_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>
#include <TMath.h>
#include <TROOT.h>

#include "DEventProcessor_fdc_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCSegment.h>
#include <FDC/DFDCCathodeCluster.h>

#include <BCAL/DBCALShower.h>
#include <DVector2.h>

#define EPS 1e-3
#define ITER_MAX 10
#define ADJACENT_MATCH_RADIUS 1.0
#define MATCH_RADIUS 2.0

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_fdc_hists());
}
} // "C"


bool dEdxSort(double a,double b){
  return (a<b);
}


//------------------
// DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::DEventProcessor_fdc_hists()
{
	fdc_ptr = &fdc;
	fdchit_ptr = &fdchit;

	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::~DEventProcessor_fdc_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_fdc_hists::init(void)
{
  cout << "initializing" <<endl;

	mT0=0.;
	myevt=0;

	DoAlign=false;
	//DoAlign=true;
	alignments.resize(24);
	if (DoAlign){
	  for (unsigned int i=0;i<24;i++){
	    alignments[i].A(kDx)=0.;
	    alignments[i].A(kDy)=0.;
	    alignments[i].A(kDPhi)=0.;
	    alignments[i].E(kDx,kDx)=1.0;
	    alignments[i].E(kDy,kDy)=1.0;
	    alignments[i].E(kDPhi,kDPhi)=1.0;
	  }
	}
	
	TDirectory *dir = new TDirectoryFile("FDC","FDC");
	dir->cd();

	// Create Tree
	fdctree = new TTree("fdc","FDC algnments");
	fdcbranch = fdctree->Branch("T","FDC_branch",&fdc_ptr);

	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_fdc_hists::brun(JEventLoop *loop, int runnumber)
{	
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  dgeom->GetFDCWires(fdcwires);

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_dz,endplate_rmin,endplate_rmax;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z+=0.5*endplate_dz;

  dapp->Lock();
    

  Hqratio_vs_wire= (TH2F *)gROOT->FindObject("Hqratio_vs_wire");
  if (!Hqratio_vs_wire)
    Hqratio_vs_wire= new TH2F("Hqratio_vs_wire","Charge ratio vs wire number",
			      2304,0.5,2304.5,100,-0.5,0.5);

  Hdelta_z_vs_wire= (TH2F*)gROOT->FindObject("Hdelta_z_vs_wire");
  if (!Hdelta_z_vs_wire)
    Hdelta_z_vs_wire= new TH2F("Hdelta_z_vs_wire","wire z offset vs wire number",
			      2304,0.5,2304.5,100,-0.1,0.1);

    Hprob = (TH1F*)gROOT->FindObject("Hprob");
    if (!Hprob) 
      Hprob=new TH1F("Hprob","Confidence level for wire-based fit",
			  100,0.0,1.); 

    Hreduced_chi2 = (TH1F*)gROOT->FindObject("Hreduced_chi2");
    if (!Hreduced_chi2) 
      Hreduced_chi2=new TH1F("Hreduced_chi2","chi^2/ndf",
			  1000,0,100); 
    Htime_prob = (TH1F*)gROOT->FindObject("Htime_prob");
    if (!Htime_prob) 
      Htime_prob=new TH1F("Htime_prob","Confidence level for time-based fit",
			  100,0,1);

    Hures_vs_layer=(TH2F*)gROOT->FindObject("Hures_vs_layer");
    if (!Hures_vs_layer){
      Hures_vs_layer=new TH2F("Hures_vs_layer","wire-based residuals",
				  24,0.5,24.5,200,-0.5,0.5);
    }  
    Hvres_vs_layer=(TH2F*)gROOT->FindObject("Hvres_vs_layer");
    if (!Hvres_vs_layer){
      Hvres_vs_layer=new TH2F("Hvres_vs_layer","residual for position along wire",
				  24,0.5,24.5,200,-0.5,0.5);
    } 
  
    Hcand_ty_vs_tx=(TH2F*)gROOT->FindObject("Hcand_ty_vs_tx");
    if (!Hcand_ty_vs_tx){
      Hcand_ty_vs_tx=new TH2F("Hcand_ty_vs_tx","candidate ty vs tx",100,-1,1,
			      100,-1,1);
    }    
    Hty_vs_tx=(TH2F*)gROOT->FindObject("Hty_vs_tx");
    if (!Hty_vs_tx){
      Hty_vs_tx=new TH2F("Hty_vs_tx","wire-based ty vs tx",100,-1,1,
			      100,-1,1);
    }  
    Htime_ty_vs_tx=(TH2F*)gROOT->FindObject("Htime_ty_vs_tx");
    if (!Htime_ty_vs_tx){
      Htime_ty_vs_tx=new TH2F("Htime_ty_vs_tx","time-based ty vs tx",100,-1,1,
			      100,-1,1);
    }
    Htime_y_vs_x=(TH3F*)gROOT->FindObject("Htime_y_vs_x");
    if (!Htime_y_vs_x){
      Htime_y_vs_x=new TH3F("Htime_y_vs_x","time-based y vs x",24,0.5,24.5,100,-48.5,48.5,
			    100,-48.5,48.5);
    }
    
    Hdrift_time=(TH2F*)gROOT->FindObject("Hdrift_time");
    if (!Hdrift_time){
      Hdrift_time=new TH2F("Hdrift_time",
			   "doca vs drift time",201,-21,381,100,0,1);
    }  
    Hdrift_integral=(TH1F*)gROOT->FindObject("Hdrift_integral");
    if (!Hdrift_integral){
      Hdrift_integral=new TH1F("Hdrift_integral",
			   "drift time integral",140,-20,260);
    } 
    Hz_target=(TH1F*)gROOT->FindObject("Hz_target");
    if (!Hz_target){
      Hz_target=new TH1F("Hz_target",
			   "Target position",260,0,130);
    }

    Hres_vs_drift_time=(TH2F*)gROOT->FindObject("Hres_vs_drift_time");
    if (!Hres_vs_drift_time){
      Hres_vs_drift_time=new TH2F("Hres_vs_drift_time","Residual vs drift time",320,-20,300,1000,-1,1);
    } 
    Hdv_vs_dE=(TH2F*)gROOT->FindObject("Hdv_vs_dE");
    if (!Hdv_vs_dE){
      Hdv_vs_dE=new TH2F("Hdv_vs_dE","dv vs energy dep",100,0,20e-6,200,-1,1);
    }

    Hxcand_prob = (TH1F *)gROOT->FindObject("Hscand_prob");
    if (!Hxcand_prob){
      Hxcand_prob=new TH1F("Hxcand_prob","x candidate prob",100,0,1);
    }  
    Hycand_prob = (TH1F *)gROOT->FindObject("Hycand_prob");
    if (!Hycand_prob){
      Hycand_prob=new TH1F("Hycand_prob","y candidate prob",100,0,1);
    } 
    
    Hfcal_match=(TH1F*)gROOT->FindObject("Hfcal_match");
    if (!Hfcal_match){
      Hfcal_match=new TH1F("Hfcal_match",
			   "match dr for FCAL",200,0,10);
    }  

    Hbcal_match=(TH1F*)gROOT->FindObject("Hbcal_match");
    if (!Hbcal_match){
      Hbcal_match=new TH1F("Hbcal_match",
			   "match dr for BCAL",200,0,10);
    }  
     

    Htheta=(TH1F*)gROOT->FindObject("Htheta");
    if (!Htheta){
      Htheta=new TH1F("Htheta",
			   "#theta",100,0,25);
    }  
    HdEdx=(TH1F*)gROOT->FindObject("HdEdx");
    if (!HdEdx){
      HdEdx=new TH1F("HdEdx","dEdx",100,0,10e-6);
    }  
  
  dapp->Unlock();

  JCalibration *jcalib = dapp->GetJCalibration(0);  // need run number here
  vector< map<string, float> > tvals;
  if (jcalib->Get("FDC/fdc_drift_Bzero", tvals)==false){
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];   
      map<string,float>::iterator iter = row.begin();
      fdc_drift_table[i] = iter->second;
      Hdrift_integral->Fill(2.*i-20,fdc_drift_table[i]/0.5);
    }
  }
  else{
    jerr << " FDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }



  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_fdc_hists::erun(void)
{
 

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_fdc_hists::fini(void)
{

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_fdc_hists::evnt(JEventLoop *loop, int eventnumber){
  myevt++;
  
  // Get BCAL showers, FCAL showers and FDC space points
  vector<const DFCALShower*>fcalshowers;
  loop->Get(fcalshowers);
  vector<const DBCALShower*>bcalshowers;
  //loop->Get(bcalshowers);
  vector<const DFDCPseudo*>pseudos;
  loop->Get(pseudos);

  for (unsigned int i=0;i<pseudos.size();i++){
    vector<const DFDCCathodeCluster*>cathode_clusters;
    pseudos[i]->GetT(cathode_clusters);
    
    unsigned int wire_number
	=96*(pseudos[i]->wire->layer-1)+pseudos[i]->wire->wire;
    double q1=0,q2=0;
    for (unsigned int j=0;j<cathode_clusters[0]->members.size();j++){
      double q=cathode_clusters[0]->members[j]->q;
      if (q>q1) q1=q;
    }   
    for (unsigned int j=0;j<cathode_clusters[1]->members.size();j++){
      double q=cathode_clusters[1]->members[j]->q;
      if (q>q2) q2=q;
    }
    
    double ratio=q2/q1;
    Hqratio_vs_wire->Fill(wire_number,ratio-1.);
    Hdelta_z_vs_wire->Fill(wire_number,0.5336*(1.-ratio)/(1.+ratio));
  }

  if (pseudos.size()>2 && (bcalshowers.size()>0 || fcalshowers.size()>0)){
    // Group FDC hits by package
    vector<const DFDCPseudo*>packages[4];
    for (unsigned int i=0;i<pseudos.size();i++){
      packages[(pseudos[i]->wire->layer-1)/6].push_back(pseudos[i]);
    }
    
    // Link hits in each package together into track segments
    vector<segment_t>segments[4];
    for (unsigned int i=0;i<4;i++){  
      FindSegments(packages[i],segments[i]);
    }
    
    // Link the segments together to form track candidadates
    vector<vector<const DFDCPseudo *> >LinkedSegments;
    LinkSegments(segments,LinkedSegments);
  
    // Loop over linked segments
    for (unsigned int k=0;k<LinkedSegments.size();k++){
      vector<const DFDCPseudo *>hits=LinkedSegments[k];
    
      // Perform preliminary line fits for current set of linked segments    
      double var_x,var_tx,cov_x_tx,var_y,var_ty,cov_y_ty,chi2x,chi2y;
      int myndf=hits.size()-2; 
      DMatrix4x1 S=FitLine(hits,var_x,cov_x_tx,var_tx,chi2x,var_y,cov_y_ty,
			 var_ty,chi2y);     
      double probx=TMath::Prob(chi2x,myndf);
      double proby=TMath::Prob(chi2y,myndf);
    
      Hxcand_prob->Fill(probx);
      Hycand_prob->Fill(proby);
      
      // Match to outer detectors
      bool got_match=false;
      // First match to FCAL 
      double dz=0.,outer_time=0.;
      double drmin=1000.;
      for (unsigned int i=0;i<fcalshowers.size();i++){
	double fcal_z=fcalshowers[i]->getPosition().z();
	double x=S(state_x)+fcal_z*S(state_tx);
	double y=S(state_y)+fcal_z*S(state_ty);
	double dx=fcalshowers[i]->getPosition().x()-x;
	double dy=fcalshowers[i]->getPosition().y()-y;
	double dr=sqrt(dx*dx+dy*dy);

	if (dr<drmin){
	  drmin=dr;
	  dz=fcal_z-endplate_z;
	  outer_time=fcalshowers[i]->getTime();
	}
	
      }
      Hfcal_match->Fill(drmin);
      if (drmin<4.){
	got_match=true;
	//	outer_time-=2.218; //empirical correction
      }
      else{
	// Match to BCAL
	drmin=1000.;
	for (unsigned int i=0;i<bcalshowers.size();i++){
	  double bcal_z=bcalshowers[i]->z; 
	  double R2=bcalshowers[i]->x*bcalshowers[i]->x
	    +bcalshowers[i]->y*bcalshowers[i]->y;	 
	  double x0=S(state_x);
	  double y0=S(state_y);
	  double tx=S(state_tx);
	  double ty=S(state_ty);
	  double myz=(-x0*tx-y0*ty
		      +sqrt(R2*(tx*tx+ty*ty)-(x0*ty-y0*tx)*(x0*ty-y0*tx)))
	    /(tx*tx+ty*ty);
	  double x=x0+myz*tx;
	  double y=y0+myz*ty;
	  double dx=bcalshowers[i]->x-x;
	  double dy=bcalshowers[i]->y-y;
	  double dr=sqrt(dx*dx+dy*dy);
  
	  if (dr<drmin){
	    drmin=dr;
	    dz=bcal_z-endplate_z;
	    outer_time=bcalshowers[i]->t;
	  }	
	}
	Hbcal_match->Fill(drmin);

	if (drmin<4.) 
	  got_match=true;
      }
      if (got_match){
	Hcand_ty_vs_tx->Fill(S(state_tx),S(state_ty));  
	
	//compute tangent of dip angle and related angular quantities
	double tanl=1./sqrt(S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty)); 
	double theta=90-180/M_PI*atan(tanl); 
	double sinl=sin(atan(tanl));
    
	Htheta->Fill(theta);

	// Compute dEdx
	double dEdx=0.;
	vector<double>dElist;
	for (unsigned int i=0;i<hits.size();i++){
	  dElist.push_back(hits[i]->dE);
	}
	sort(dElist.begin(),dElist.end(),dEdxSort);
	unsigned int N=dElist.size()/2;
	for (unsigned int i=0;i<N;i++){
	  dEdx+=dElist[i];
	}
	dEdx/=N*sinl;

	HdEdx->Fill(dEdx);

	// Estimate z-position of "vertex" in target
	double one_over_tx2=1./(S(state_tx)*S(state_tx));
	double varz_x
	  =one_over_tx2*(var_x+S(state_x)*S(state_x)*one_over_tx2*var_tx
			 -2.*S(state_x)/S(state_tx)*cov_x_tx);
	double one_over_ty2=1./(S(state_ty)*S(state_ty));
	double varz_y
	  =one_over_ty2*(var_y+S(state_y)*S(state_y)*one_over_ty2*var_ty
			 -2.*S(state_x)/S(state_ty)*cov_y_ty);
	// Weighted average from x and y line fits, ignoring x-y correlations
	double z_target
	  =(-S(state_x)/(S(state_tx)*varz_x)-S(state_y)/(S(state_ty)*varz_y))/
	  (1./varz_x+1/varz_y);

	Hz_target->Fill(z_target);
		
	// Estimate for t0 at the beginning of track assuming particle is 
	// moving at the speed of light
	mT0=outer_time-dz/(29.98*sinl);

	// Run the Kalman Filter algorithm
	DoFilter(S,hits);
      }
    }
  }
   
  return NOERROR;
}

// Steering routine for the kalman filter
jerror_t 
DEventProcessor_fdc_hists::DoFilter(DMatrix4x1 &S,
				    vector<const DFDCPseudo *>&hits){
  unsigned int num_hits=hits.size();
  //  vector<strip_update_t>strip_updates(num_hits);
  //vector<strip_update_t>smoothed_strip_updates(num_hits);
  vector<update_t>updates(num_hits);
  vector<update_t>smoothed_updates(num_hits);

  int NEVENTS=200000;
  double anneal_factor=1.;
  if (DoAlign){
    anneal_factor=pow(1e6,(double(NEVENTS-myevt))/(NEVENTS-1.));
    if (myevt>NEVENTS) anneal_factor=1.;
  }
  //anneal_factor=1.;
  //anneal_factor=10000.;

  // Best guess for state vector at "vertex"
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;
  // double start_z=hits[0]->wire->origin.z()-1.;
  S(state_x)+=endplate_z*S(state_tx);
  S(state_y)+=endplate_z*S(state_ty);
  SetReferenceTrajectory(endplate_z,S,trajectory,hits);
      
  // Intial guess for covariance matrix
  DMatrix4x4 C,C0,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=1.;
  C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;
  
  // Chi-squared and degrees of freedom
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0,ndof_old=0;
  unsigned iter=0;
  for(;;){
    iter++;
    chi2_old=chi2; 
    ndof_old=ndof;
    C=C0;
    if (KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof)
	!=NOERROR) break;
	
    //printf("=======chi2 %f\n",chi2);
    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1 || iter==ITER_MAX) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S;
    
    // run the smoother (opposite direction to filter)
    Smooth(S,C,trajectory,hits,updates,smoothed_updates);
  }
      
  if (iter>1){
    double reduced_chi2=chi2_old/double(ndof_old);
    double prob=TMath::Prob(chi2_old,ndof_old);
    Hprob->Fill(prob);
    Hreduced_chi2->Fill(reduced_chi2);
	
    if (prob>0.05){
      Hty_vs_tx->Fill(Sbest(state_tx),Sbest(state_ty));
      
      for (unsigned int i=0;i<smoothed_updates.size();i++){
	unsigned int layer=hits[i]->wire->layer;
	Hures_vs_layer->Fill(layer,smoothed_updates[i].ures);
	Hvres_vs_layer->Fill(layer,smoothed_updates[i].vres);
	if (layer==1){
	  Hres_vs_drift_time->Fill(smoothed_updates[i].drift_time,
				   smoothed_updates[i].ures);
	  Hdrift_time->Fill(smoothed_updates[i].drift_time,
			    smoothed_updates[i].doca);
 
	  Hdv_vs_dE->Fill(hits[i]->dE,smoothed_updates[i].vres);
	}
      }
    }
    if (DoAlign){
      FindOffsets(hits,smoothed_updates);
      
      for (unsigned int layer=0;layer<24;layer++){
	// Set up to fill tree
	double dxr=alignments[layer].A(kDx);
	double dyr=alignments[layer].A(kDy);  
	fdc.dPhi=alignments[layer].A(kDPhi);
	double cosdphi=cos(fdc.dPhi);
	double sindphi=sin(fdc.dPhi);
	double dx=dxr*cosdphi-dyr*sindphi;
	double dy=dxr*sindphi+dyr*cosdphi;
	fdc.dX=dx;
	fdc.dY=dy;
	fdc.layer=layer;
	fdc.N=myevt;
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	fdctree->Fill();
	
	// Unlock mutex
	pthread_mutex_unlock(&mutex);
      }
    }
   
    

    /*
    printf("-------Event %d\n",myevt);
    for (unsigned int i=0;i<24;i++) 
      printf("A %f %f %f\n",alignments[i].A(kDx),alignments[i].A(kDy),
	     alignments[i].A(kDPhi));
    */

    return NOERROR;
  }

  return VALUE_OUT_OF_RANGE;
}


// Link segments from package to package by doing straight-line projections
jerror_t
DEventProcessor_fdc_hists::LinkSegments(vector<segment_t>segments[4], 
					vector<vector<const DFDCPseudo *> >&LinkedSegments){
  vector<const DFDCPseudo *>myhits;
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<segments[i].size();j++){
      if (segments[i][j].matched==false){
	myhits.assign(segments[i][j].hits.begin(),segments[i][j].hits.end());
	
	unsigned i_plus_1=i+1; 
	if (i_plus_1<4){
	  double tx=segments[i][j].S(state_tx);
	  double ty=segments[i][j].S(state_ty);
	  double x0=segments[i][j].S(state_x);
	  double y0=segments[i][j].S(state_y);
	  
	  for (unsigned int k=0;k<segments[i_plus_1].size();k++){
	    if (segments[i_plus_1][k].matched==false){
	      double z=segments[i_plus_1][k].hits[0]->wire->origin.z();
	      DVector2 proj(x0+tx*z,y0+ty*z);
	      
	      if ((proj-segments[i_plus_1][k].hits[0]->xy).Mod()<MATCH_RADIUS){
		segments[i_plus_1][k].matched=true;
		myhits.insert(myhits.end(),segments[i_plus_1][k].hits.begin(),
				segments[i_plus_1][k].hits.end());
		
		unsigned int i_plus_2=i_plus_1+1;
		if (i_plus_2<4){
		  tx=segments[i_plus_1][k].S(state_tx);
		  ty=segments[i_plus_1][k].S(state_ty);
		  x0=segments[i_plus_1][k].S(state_x);
		  y0=segments[i_plus_1][k].S(state_y);
		  
		  for (unsigned int m=0;m<segments[i_plus_2].size();m++){
		    if (segments[i_plus_2][m].matched==false){
		      z=segments[i_plus_2][m].hits[0]->wire->origin.z();
		      proj.Set(x0+tx*z,y0+ty*z);
		      
		      if ((proj-segments[i_plus_2][m].hits[0]->xy).Mod()<MATCH_RADIUS){
			segments[i_plus_2][m].matched=true;
			myhits.insert(myhits.end(),segments[i_plus_2][m].hits.begin(),
				      segments[i_plus_2][m].hits.end());
			
			unsigned int i_plus_3=i_plus_2+1;
			if (i_plus_3<4){
			  tx=segments[i_plus_2][m].S(state_tx);
			  ty=segments[i_plus_2][m].S(state_ty);
			  x0=segments[i_plus_2][m].S(state_x);
			  y0=segments[i_plus_2][m].S(state_y);
			  
			  for (unsigned int n=0;n<segments[i_plus_3].size();n++){
			    if (segments[i_plus_3][n].matched==false){
			      z=segments[i_plus_3][n].hits[0]->wire->origin.z();
			      proj.Set(x0+tx*z,y0+ty*z);
			      
			      if ((proj-segments[i_plus_3][n].hits[0]->xy).Mod()<MATCH_RADIUS){
				segments[i_plus_3][n].matched=true;
				myhits.insert(myhits.end(),segments[i_plus_3][n].hits.begin(),
					      segments[i_plus_3][n].hits.end());
				
				break;
			      } // matched a segment
			    }
			  }  // loop over last set of segments
			} // if we have another package to loop over
			break;
		      } // matched a segment
		    }
		  } // loop over second-to-last set of segments
		}
		break;
	      } // matched a segment
	    }
	  } // loop over third-to-last set of segments
	}	
	LinkedSegments.push_back(myhits);
	myhits.clear();
      } // check if we have already used this segment
    } // loop over first set of segments
  } // loop over packages
  
  return NOERROR;
}

// Find segments by associating adjacent hits within a package together.
jerror_t DEventProcessor_fdc_hists::FindSegments(vector<const DFDCPseudo*>&points,
					vector<segment_t>&segments){
  if (points.size()==0) return RESOURCE_UNAVAILABLE;
  vector<int>used(points.size());

  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=points[0]->wire->origin.z();
  vector<unsigned int>x_list;
  x_list.push_back(0);
  for (unsigned int i=0;i<points.size();i++){
    used.push_back(false);
    if (points[i]->wire->origin.z()!=old_z){
      x_list.push_back(i);
    }
    old_z=points[i]->wire->origin.z();
  }
  x_list.push_back(points.size()); 

  unsigned int start=0;
  // loop over the start indices, starting with the first plane
  while (start<x_list.size()-1){
    // Now loop over the list of track segment start points
    for (unsigned int i=x_list[start];i<x_list[start+1];i++){
      if (used[i]==false){
	used[i]=true;
	
	// Point in the current plane in the package 
	DVector2 XY=points[i]->xy;
	
	// Create list of nearest neighbors
	vector<const DFDCPseudo*>neighbors;
	neighbors.push_back(points[i]);
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<x_list.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	    delta=(XY-points[m]->xy).Mod();
	    if (delta<delta_min && delta<MATCH_RADIUS){
	      delta_min=delta;
	      match=m;
	    }
	  }	
	  if (match!=0 
	      && used[match]==false
	      ){
	    XY=points[match]->xy;
	    used[match]=true;
	    neighbors.push_back(points[match]);	  
	  }
	}
	unsigned int num_neighbors=neighbors.size();

	bool do_sort=false;
	// Look for hits adjacent to the ones we have in our segment candidate
	for (unsigned int k=0;k<points.size();k++){
	  if (!used[k]){
	    for (unsigned int j=0;j<num_neighbors;j++){
	      delta=(points[k]->xy-neighbors[j]->xy).Mod();

	      if (delta<ADJACENT_MATCH_RADIUS && 
		  abs(neighbors[j]->wire->wire-points[k]->wire->wire)<=1
		  && neighbors[j]->wire->origin.z()==points[k]->wire->origin.z()){
		used[k]=true;
		neighbors.push_back(points[k]);
		do_sort=true;
	      }      
	    }
	  }
	} // loop looking for hits adjacent to hits on segment

	if (neighbors.size()>4){
	  segment_t mysegment;
	  mysegment.matched=false;
	  mysegment.S=FitLine(neighbors);
	  mysegment.hits=neighbors;
	  segments.push_back(mysegment);
	}
      }
    }// loop over start points within a plane
    
    // Look for a new plane to start looking for a segment
    while (start<x_list.size()-1){
      if (used[x_list[start]]==false) break;
      start++;
    }

  }

  return NOERROR;
}


DMatrix4x1 
DEventProcessor_fdc_hists::FitLine(vector<const DFDCPseudo*> &fdchits){
  double var_x,var_tx,cov_x_tx,chi2x;
  double var_y,var_ty,cov_y_ty,chi2y;
  
  return FitLine(fdchits,var_x,cov_x_tx,var_tx,chi2x,var_y,cov_y_ty,var_ty,
		 chi2y);
}


// Use linear regression on the hits to obtain a first guess for the state
// vector.  Method taken from Numerical Recipes in C.
DMatrix4x1 
DEventProcessor_fdc_hists::FitLine(vector<const DFDCPseudo*> &fdchits,
				   double &var_x,double &cov_x_tx,
				   double &var_tx,double &chi2x,
				   double &var_y,double &cov_y_ty,
				   double &var_ty,double &chi2y){
  double S1=0.;
  double S1z=0.;
  double S1y=0.;
  double S1zz=0.;
  double S1zy=0.;  
  double S2=0.;
  double S2z=0.;
  double S2x=0.;
  double S2zz=0.;
  double S2zx=0.;

  double sig2v=0.01; // rough guess;

  for (unsigned int i=0;i<fdchits.size();i++){
    double cosa=fdchits[i]->wire->udir.y();
    double sina=fdchits[i]->wire->udir.x();
    double x=fdchits[i]->xy.X();
    double y=fdchits[i]->xy.Y();
    double z=fdchits[i]->wire->origin.z();
    double sig2x=cosa*cosa/12+sina*sina*sig2v;
    double sig2y=sina*sina/12+cosa*cosa*sig2v;
    double one_over_var1=1/sig2y;
    double one_over_var2=1/sig2x;

    S1+=one_over_var1;
    S1z+=z*one_over_var1;
    S1y+=y*one_over_var1;
    S1zz+=z*z*one_over_var1;
    S1zy+=z*y*one_over_var1;    
    
    S2+=one_over_var2;
    S2z+=z*one_over_var2;
    S2x+=x*one_over_var2;
    S2zz+=z*z*one_over_var2;
    S2zx+=z*x*one_over_var2;
  }
  double D1=S1*S1zz-S1z*S1z;
  double y_intercept=(S1zz*S1y-S1z*S1zy)/D1;
  double y_slope=(S1*S1zy-S1z*S1y)/D1;
  double D2=S2*S2zz-S2z*S2z;
  double x_intercept=(S2zz*S2x-S2z*S2zx)/D2;
  double x_slope=(S2*S2zx-S2z*S2x)/D2;

  // Covariance matrix for x
  var_x=S2zz/D2;
  var_tx=S2/D2;
  cov_x_tx=-S2z/D2;

  // Covariance matrix for y
  var_y=S1zz/D1;
  var_ty=S1/D1;
  cov_y_ty=-S1z/D1;

  // Compute chi2 for the line fits, ignoring correlations between x and y
  chi2x=0;
  chi2y=0;
  for (unsigned int i=0;i<fdchits.size();i++){
    double cosa=fdchits[i]->wire->udir.y();
    double sina=fdchits[i]->wire->udir.x(); 
    double sig2x=cosa*cosa/12+sina*sina*sig2v;
    double sig2y=sina*sina/12+cosa*cosa*sig2v;
    double one_over_var1=1/sig2y;
    double one_over_var2=1/sig2x;

    double z=fdchits[i]->wire->origin.z();
    double dx=fdchits[i]->xy.X()-(x_intercept+x_slope*z);
    double dy=fdchits[i]->xy.Y()-(y_intercept+y_slope*z);

    chi2x+=dx*dx*one_over_var2;
    chi2y+=dy*dy*one_over_var1;
  }
  return DMatrix4x1(x_intercept,y_intercept,x_slope,y_slope);

}


// Kalman smoother 
jerror_t DEventProcessor_fdc_hists::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
				      deque<trajectory_t>&trajectory,
				      vector<const DFDCPseudo *>&hits,
				      vector<strip_update_t>updates,
				    vector<strip_update_t>&smoothed_updates
					   ){
  DMatrix4x1 S; 
  DMatrix4x4 C,dC;
  DMatrix4x4 JT,A;

  unsigned int max=trajectory.size()-1;
  S=(trajectory[max].Skk);
  C=(trajectory[max].Ckk);
  JT=(trajectory[max].J.Transpose());
  //Ss=S;
  //Cs=C;
  for (unsigned int m=max-1;m>0;m--){
    if (trajectory[m].h_id==0){
      A=trajectory[m].Ckk*JT*C.Invert();
      Ss=trajectory[m].Skk+A*(Ss-S);
      Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    else if (trajectory[m].h_id>0){
      unsigned int id=trajectory[m].h_id-1;
      A=updates[id].C*JT*C.Invert();
      dC=A*(Cs-C)*A.Transpose();
      Ss=updates[id].S+A*(Ss-S);
      Cs=updates[id].C+dC;
      /*
      printf("-------\n");
      updates[id].C;
      Cs.Print();
      */

      // Nominal rotation of wire planes
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      
      // State vector
      double x=Ss(state_x);
      double y=Ss(state_y);
      double tx=Ss(state_tx);
      double ty=Ss(state_ty);
 
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x3 E=alignments[layer].E;
      DMatrix3x1 A=alignments[layer].A;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));

      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;

      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);

      // Smoothed residuals
      double uwire=hits[id]->w;
      double v=hits[id]->s;  
      double d=(upred-uwire)*cosalpha;
      smoothed_updates[id].vres=v-vpred+tv*d*sinalpha;
      smoothed_updates[id].ures=(d>0?1.:-1.)*updates[id].drift-d; 
      
      smoothed_updates[id].id=trajectory[m].h_id;
      smoothed_updates[id].drift=updates[id].drift;
      smoothed_updates[id].drift_time=updates[id].drift_time;
      smoothed_updates[id].S=Ss;
      smoothed_updates[id].C=Cs;
      smoothed_updates[id].R=updates[id].R-updates[id].H*dC*updates[id].H_T;
    }
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  /*
    printf("-----end\n");
  Ss.Print();
  Cs.Print();
  printf("--ckk \n");
  C.Print();
  */

  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}

// Kalman smoother 
jerror_t DEventProcessor_fdc_hists::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
					   deque<trajectory_t>&trajectory,
					   vector<const DFDCPseudo *>&hits,
					   vector<update_t>updates,
					   vector<update_t>&smoothed_updates
					   ){
  DMatrix4x1 S; 
  DMatrix4x4 C,dC;
  DMatrix4x4 JT,A;

  unsigned int max=trajectory.size()-1;
  S=(trajectory[max].Skk);
  C=(trajectory[max].Ckk);
  JT=(trajectory[max].J.Transpose());
  //Ss=S;
  //Cs=C;
  for (unsigned int m=max-1;m>0;m--){
    if (trajectory[m].h_id==0){
      A=trajectory[m].Ckk*JT*C.Invert();
      Ss=trajectory[m].Skk+A*(Ss-S);
      Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    else if (trajectory[m].h_id>0){
      unsigned int id=trajectory[m].h_id-1;
      A=updates[id].C*JT*C.Invert();
      dC=A*(Cs-C)*A.Transpose();
      Ss=updates[id].S+A*(Ss-S);
      Cs=updates[id].C+dC;
      /*
      printf("-------\n");
      updates[id].C;
      Cs.Print();
      */
      
        // Nominal rotation of wire planes
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      
      // State vector
      double x=Ss(state_x);
      double y=Ss(state_y);
      double tx=Ss(state_tx);
      double ty=Ss(state_ty);
 
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x3 E=alignments[layer].E;
      DMatrix3x1 A=alignments[layer].A;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));

      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;

      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);

      // Smoothed residuals
      double uwire=hits[id]->w;
      double v=hits[id]->s;  
      double d=(upred-uwire)*cosalpha;
      smoothed_updates[id].vres=v-vpred+tv*d*sinalpha;
      smoothed_updates[id].ures=(d>0?1.:-1.)*updates[id].drift-d;
      smoothed_updates[id].doca=fabs(d);

      smoothed_updates[id].id=trajectory[m].h_id;
      smoothed_updates[id].drift=updates[id].drift;
      smoothed_updates[id].drift_time=updates[id].drift_time;
      smoothed_updates[id].S=Ss;
      smoothed_updates[id].C=Cs;
      smoothed_updates[id].R=updates[id].R-updates[id].H*dC*updates[id].H_T;
    }
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  /*
    printf("-----end\n");
  Ss.Print();
  Cs.Print();
  printf("--ckk \n");
  C.Print();
  */

  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}

// Perform Kalman Filter for the current trajectory
jerror_t 
DEventProcessor_fdc_hists::KalmanFilter(double anneal_factor,
					DMatrix4x1 &S,DMatrix4x4 &C,
			       vector<const DFDCPseudo *>&hits,
			       deque<trajectory_t>&trajectory,
			       vector<strip_update_t>&updates,
			       double &chi2,unsigned int &ndof){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  double V=0;  // Measurement covariance 
  double Mdiff=0.;
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  trajectory[0].Skk=S;
  trajectory[0].Ckk=C;
  for (unsigned int k=1;k<trajectory.size();k++){
    if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
      return UNRECOVERABLE_ERROR;
    
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();
    
    // Save the current state and covariance matrix 
    trajectory[k].Skk=S;
    trajectory[k].Ckk=C;

    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;

    // Correct S and C for the hit 
    if (trajectory[k].h_id>0){
      unsigned int id=trajectory[k].h_id-1;
      
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      
      // State vector
      double x=S(state_x);
      double y=S(state_y);
      double tx=S(state_tx);
      double ty=S(state_ty);
      if (isnan(x) || isnan(y)) return UNRECOVERABLE_ERROR;
 
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x3 E=alignments[layer].E;
      DMatrix3x1 A=alignments[layer].A;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));

      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;

      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);
       
      // Difference between measurement and projection
      for (int m=trajectory[k].num_hits-1;m>=0;m--){
	unsigned int my_id=id+m;
	double uwire=hits[my_id]->w;
	double v=hits[my_id]->s;

	//printf("cospsi %f sinpsi %f\n",cospsi,sinpsi);
	//printf("u %f %f v %f %f\n",uwire,upred,v,vpred);

	double drift_time=hits[my_id]->time-trajectory[k].t-mT0;
	double drift=GetDriftDistance(drift_time);
	updates[my_id].drift=drift;
	updates[my_id].drift_time=drift_time;

	double du=upred-uwire;
	double d=du*cosalpha;

	Mdiff=v-vpred+tv*d*sinalpha;
	
	//printf("tdrift %f d %f %f\n",drift_time,drift,d);
	//printf("dv %f ddv %f\n",v-vpred,tv*d*sinalpha);

	// Variance of measurement error
	double sigma=3.7e-8/hits[my_id]->dE+0.0062;
	V=anneal_factor*sigma*sigma;
	
	// Matrix for transforming from state-vector space to measurement space
	double sinalpha_cosalpha=sinalpha*cosalpha;
	H_T(state_x)=sinpsi-tv*cospsi*sinalpha_cosalpha;
	H_T(state_y)=cospsi+tv*sinpsi*sinalpha_cosalpha;

	double temp=tv*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha);
	H_T(state_tx)=-d*(sinpsi*sinalpha+cospsi*temp);
	H_T(state_ty)=-d*(cospsi*sinalpha-sinpsi*temp);

	// H-matrix transpose
	H(state_x)=H_T(state_x);
	H(state_y)=H_T(state_y);
	H(state_tx)=H_T(state_tx);
	H(state_ty)=H_T(state_ty);
		
	updates[my_id].H=H;
	updates[my_id].H_T=H_T;
		
	double Vtemp=V+H*C*H_T;

	if (DoAlign){	
	  // Matrices to rotate alignment error matrix into measurement space
	  DMatrix1x3 G;
	  DMatrix3x1 G_T;
	  
	  G_T(kDx)=-sina+tv*cosa*sinalpha_cosalpha;	
	  G_T(kDy)=-cosa-tv*sina*sinalpha_cosalpha;
	  
	  G_T(kDPhi)=-x*cospsi+y*sinpsi
	    +tu*d*sinalpha-tv*(x*sinpsi+y*cospsi)*sinalpha_cosalpha
	    -tu*tv*d*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha)
	    ;
	  
	  // G-matrix transpose
	  G(kDx)=G_T(kDx);
	  G(kDy)=G_T(kDy);
	  G(kDPhi)=G_T(kDPhi);
	  
	  Vtemp=Vtemp+G*E*G_T;
	}
	
	// Compute Kalman gain matrix
	K=(1./Vtemp)*(C*H_T);
	
	// Update the state vector 
	S+=Mdiff*K;
	updates[my_id].S=S;

	// Update state vector covariance matrix
	C=C-K*(H*C);    
	updates[my_id].C=C;
	
	// Update chi2 for this trajectory
	double scale=1-H*K;
	updates[my_id].R=scale*V;

	//printf("chi2 %f for %d\n",RC.Chi2(R),my_id);
	
	chi2+=scale*Mdiff*Mdiff/V;
	ndof++;
      }

    }

  }
  chi2*=anneal_factor;

  ndof-=4;

  return NOERROR;
}




// Perform Kalman Filter for the current trajectory
jerror_t 
DEventProcessor_fdc_hists::KalmanFilter(double anneal_factor,
					DMatrix4x1 &S,DMatrix4x4 &C,
			       vector<const DFDCPseudo *>&hits,
			       deque<trajectory_t>&trajectory,
			       vector<update_t>&updates,
			       double &chi2,unsigned int &ndof){
  DMatrix2x4 H;  // Track projection matrix
  DMatrix4x2 H_T; // Transpose of track projection matrix 
  DMatrix4x2 K;  // Kalman gain matrix
  DMatrix2x2 V(0.020833,0,0,0);  // Measurement covariance 
  DMatrix2x2 Vtemp;
  DMatrix2x1 Mdiff;
  DMatrix2x2 InvV; // Inverse of error matrix
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  trajectory[0].Skk=S;
  trajectory[0].Ckk=C;
  for (unsigned int k=1;k<trajectory.size();k++){
    if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
      return UNRECOVERABLE_ERROR;
    
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();
    
    // Save the current state and covariance matrix 
    trajectory[k].Skk=S;
    trajectory[k].Ckk=C;

    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;

    // Correct S and C for the hit 
    if (trajectory[k].h_id>0){
      unsigned int id=trajectory[k].h_id-1;
      
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      
      // State vector
      double x=S(state_x);
      double y=S(state_y);
      double tx=S(state_tx);
      double ty=S(state_ty);
      if (isnan(x) || isnan(y)) return UNRECOVERABLE_ERROR;
 
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x3 E=alignments[layer].E;
      DMatrix3x1 A=alignments[layer].A;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));

      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;

      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);
       
      // Difference between measurement and projection
      for (int m=trajectory[k].num_hits-1;m>=0;m--){
	unsigned int my_id=id+m;
	double uwire=hits[my_id]->w;
	double v=hits[my_id]->s;

	//printf("cospsi %f sinpsi %f\n",cospsi,sinpsi);
	//printf("u %f %f v %f %f\n",uwire,upred,v,vpred);
	

	// Find drift distance
	double drift_time=hits[my_id]->time-trajectory[k].t-mT0;
	double drift=GetDriftDistance(drift_time);
	updates[my_id].drift=drift;
	updates[my_id].drift_time=drift_time;

	double du=upred-uwire;
	double d=du*cosalpha;

	double sign=(du>0)?1.:-1.;
	
	// Difference between measured and predicted vectors
	//drift=0.25;
	Mdiff(0)=sign*drift-d;
	Mdiff(1)=v-vpred+tv*d*sinalpha;
	
	//printf("tdrift %f d %f %f\n",drift_time,drift,d);
	//printf("dv %f ddv %f\n",v-vpred,tv*d*sinalpha);

	// Variance of measurement error
	V(0,0)=anneal_factor*GetDriftVariance(drift_time);
	
	double sigma=3.7e-8/hits[my_id]->dE+0.0062;
	V(1,1)=anneal_factor*sigma*sigma;
	
	// Matrix for transforming from state-vector space to measurement space
	double sinalpha_cosalpha=sinalpha*cosalpha;
	H_T(state_x,0)=cospsi*cosalpha;
	H_T(state_y,0)=-sinpsi*cosalpha;

	H_T(state_x,1)=sinpsi-tv*cospsi*sinalpha_cosalpha;
	H_T(state_y,1)=cospsi+tv*sinpsi*sinalpha_cosalpha;
	
	double temp=d*sinalpha_cosalpha;
	H_T(state_tx,0)=-temp*cospsi;
	H_T(state_ty,0)=+temp*sinpsi;

	temp=tv*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha);
	H_T(state_tx,1)=-d*(sinpsi*sinalpha+cospsi*temp);
	H_T(state_ty,1)=-d*(cospsi*sinalpha-sinpsi*temp);

	// H-matrix transpose
	H(0,state_x)=H_T(state_x,0);
	H(0,state_y)=H_T(state_y,0);

	H(1,state_x)=H_T(state_x,1);
	H(1,state_y)=H_T(state_y,1);

	H(0,state_tx)=H_T(state_tx,0);
	H(0,state_ty)=H_T(state_ty,0);
	
	updates[my_id].H=H;
	updates[my_id].H_T=H_T;
		
	DMatrix2x2 Vtemp=V+H*C*H_T;

	if (DoAlign){	
	  // Matrices to rotate alignment error matrix into measurement space
	  DMatrix2x3 G;
	  DMatrix3x2 G_T;
	  
	  G_T(kDx,0)=-cosa*cosalpha;
	  G_T(kDy,0)=+sina*cosalpha;
	  
	  G_T(kDx,1)=-sina+tv*cosa*sinalpha_cosalpha;	
	  G_T(kDy,1)=-cosa-tv*sina*sinalpha_cosalpha;
	  
	  G_T(kDPhi,0)=(x*sinpsi+y*cospsi)*cosalpha;
	  G_T(kDPhi,1)=-x*cospsi+y*sinpsi
	    +tu*d*sinalpha-tv*(x*sinpsi+y*cospsi)*sinalpha_cosalpha
	    -tu*tv*d*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha)
	    ;
	  
	  // G-matrix transpose
	  G(0,kDx)=G_T(kDx,0);
	  G(0,kDy)=G_T(kDy,0);
	  G(1,kDx)=G_T(kDx,1);
	  G(1,kDy)=G_T(kDy,1);
	  G(0,kDPhi)=G_T(kDPhi,0);
	  G(1,kDPhi)=G_T(kDPhi,1);
	  
	  Vtemp=Vtemp+G*E*G_T;
	}

	// Variance for this hit
	InvV=Vtemp.Invert();
	
	// Compute Kalman gain matrix
	K=(C*H_T)*InvV;
	
	// Update the state vector 
	S+=K*Mdiff;
	updates[my_id].S=S;

	// Update state vector covariance matrix
	C=C-K*(H*C);    
	updates[my_id].C=C;
	
	// Update chi2 for this trajectory
	DMatrix2x1 R=Mdiff-H*K*Mdiff;
	DMatrix2x2 RC=V-H*K*V;
	updates[my_id].R=RC;

	//printf("chi2 %f for %d\n",RC.Chi2(R),my_id);
	
	chi2+=RC.Chi2(R);
	ndof+=2;
      }

    }

  }
  chi2*=anneal_factor;
  ndof-=4;

  return NOERROR;
}

// Reference trajectory for the track
jerror_t DEventProcessor_fdc_hists
::SetReferenceTrajectory(double z,DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 vector<const DFDCPseudo *>&pseudos){
  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);
    
  double dz=1.1;
  double t=0.;
  trajectory_t temp;
  temp.S=S;
  temp.J=J;
  temp.Skk=Zero4x1;
  temp.Ckk=Zero4x4;
  temp.h_id=0;	  
  temp.num_hits=0;
  temp.z=z;
  temp.t=0.;
  trajectory.push_front(temp);
  double zhit=z;
  double old_zhit=z;
  unsigned int itrajectory=0;
  for (unsigned int i=0;i<pseudos.size();i++){  
    zhit=pseudos[i]->wire->origin.z();

    if (fabs(zhit-old_zhit)<EPS){
      trajectory_t temp;
      temp.J=J;
      temp.S=trajectory[0].S;
      temp.t=trajectory[0].t;
      temp.h_id=i+1;
      temp.z=trajectory[0].z;
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;
      trajectory.push_front(temp); 

      continue;
    }
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      
      trajectory_t temp;
      temp.J=J;
      temp.J(state_x,state_tx)=dz;
      temp.J(state_y,state_ty)=dz;
      // Flight time: assume particle is moving at the speed of light
      temp.t=(t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
	      /29.98);
      //propagate the state to the next z position
      temp.S(state_x)=S(state_x)+S(state_tx)*dz;
      temp.S(state_y)=S(state_y)+S(state_ty)*dz;
      temp.S(state_tx)=S(state_tx);
      temp.S(state_ty)=S(state_ty);
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;  
      temp.h_id=0;	  
      temp.num_hits=0;
      if (new_z>zhit){
	new_z=zhit;
	temp.h_id=i+1;
	temp.num_hits=1;
	done=true;
      }
      temp.z=new_z;
      trajectory.push_front(temp); 
      S=temp.S;
      itrajectory++;
      z=new_z;
    }	   
    old_zhit=zhit;
  }
  temp.Skk=Zero4x1;
  temp.Ckk=Zero4x4;
  temp.h_id=0;	  
  temp.z=z+dz;
  temp.J=J;
  temp.J(state_x,state_tx)=-dz;
  temp.J(state_y,state_ty)=-dz;
  // Flight time: assume particle is moving at the speed of light
  temp.t=(t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
	  /29.98);
  //propagate the state to the next z position
  temp.S(state_x)=S(state_x)+S(state_tx)*dz;
  temp.S(state_y)=S(state_y)+S(state_ty)*dz;
  temp.S(state_tx)=S(state_tx);
  temp.S(state_ty)=S(state_ty);
  S=temp.S;
  trajectory.push_front(temp);

  if (false){
    printf("Trajectory:\n");
    for (unsigned int i=0;i<trajectory.size();i++){
    printf(" x %f y %f z %f hit %d\n",trajectory[i].S(state_x),
	   trajectory[i].S(state_y),trajectory[i].z,trajectory[i].h_id); 
    }
  }

  return NOERROR;
}

// Crude approximation for the variance in drift distance due to smearing
double DEventProcessor_fdc_hists::GetDriftVariance(double t){
  if (t<0) t=0;
  if (t>115)t=115;
  double tpar[5]={0.0162,-0.001064,3.59e-5,-5.04e-7,2.46e-9};
  double sigma=tpar[0]+0.008;
  for (int i=1;i<5;i++) sigma+=tpar[i]*pow(t,i);

  //double sigma=0.0159-0.000281*t+3.82e-6*t*t;

  return sigma*sigma;
}

#define FDC_T0_OFFSET 20.
// Interpolate on a table to convert time to distance for the fdc
double DEventProcessor_fdc_hists::GetDriftDistance(double t){
  int id=int((t+FDC_T0_OFFSET)/2.);
  if (id<0) id=0;
  if (id>138) id=138;
  double d=fdc_drift_table[id];

  if (id!=138){
    double frac=0.5*(t+FDC_T0_OFFSET-2.*double(id));
    double dd=fdc_drift_table[id+1]-fdc_drift_table[id];
    d+=frac*dd;
  }
  return d;
}

jerror_t 
DEventProcessor_fdc_hists::FindOffsets(vector<const DFDCPseudo *>&hits,
				      vector<update_t>smoothed_updates){
  DMatrix2x3 G;//matrix relating alignment vector to measurement coords
  DMatrix3x2 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();

  for (unsigned int i=0;i<num_hits;i++){
    if (smoothed_updates[i].id>0){
      unsigned int id=smoothed_updates[i].id-1;
      double x=smoothed_updates[i].S(state_x);
      double y=smoothed_updates[i].S(state_y);
      double tx=smoothed_updates[i].S(state_tx);
      double ty=smoothed_updates[i].S(state_ty);
      
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      double uwire=hits[id]->w;
      double v=hits[id]->s;
      
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x1 A=alignments[layer].A;
      DMatrix3x3 E=alignments[layer].E;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));
      
      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;
      
      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;
      double du=upred-uwire;
      double sign=(du>0)?1.:-1.;
      
      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);
      
      // Transform from alignment vector coords to measurement coords
      G_T(kDx,0)=-cosa*cosalpha;
      G_T(kDy,0)=+sina*cosalpha;
	    
      double sinalpha_cosalpha=sinalpha*cosalpha;
      G_T(kDx,1)=-sina+tv*cosa*sinalpha_cosalpha;	
      G_T(kDy,1)=-cosa-tv*sina*sinalpha_cosalpha;
      
      double d=du*cosalpha;
      G_T(kDPhi,0)=(x*sinpsi+y*cospsi)*cosalpha;
      G_T(kDPhi,1)=-x*cospsi+y*sinpsi
	  +tu*d*sinalpha-tv*(x*sinpsi+y*cospsi)*sinalpha_cosalpha
	-tu*tv*d*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha)
	;
      
      // G-matrix transpose
      G(0,kDx)=G_T(kDx,0);
      G(0,kDy)=G_T(kDy,0);
      G(1,kDx)=G_T(kDx,1);
      G(1,kDy)=G_T(kDy,1);
      G(0,kDPhi)=G_T(kDPhi,0);
      G(1,kDPhi)=G_T(kDPhi,1);
    
      // Inverse of error matrix
      DMatrix2x2 InvV=(smoothed_updates[i].R+G*E*G_T).Invert();
      
      // Difference between measurement and projection
      DMatrix2x1 Mdiff;
      Mdiff(0)=sign*smoothed_updates[i].drift-du*cosalpha;
      Mdiff(1)=v-vpred+tv*du*cosalpha*sinalpha;
      
      // update the alignment vector and covariance
      DMatrix3x2 Ka=(E*G_T)*InvV;
      DMatrix3x1 dA=Ka*Mdiff;
      DMatrix3x3 Etemp=E-Ka*G*E;
      if (Etemp(0,0)>0 && Etemp(1,1)>0 && Etemp(2,2)>0){
	alignments[layer].E=Etemp;
	alignments[layer].A=A+Ka*Mdiff;	  
      }
    }
  }

  return NOERROR;
}

jerror_t 
DEventProcessor_fdc_hists::FindOffsets(vector<const DFDCPseudo *>&hits,
			       vector<strip_update_t>smoothed_updates){
  DMatrix1x3 G;//matrix relating alignment vector to measurement coords
  DMatrix3x1 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();

  for (unsigned int i=0;i<num_hits;i++){
    if (smoothed_updates[i].id>0){
      unsigned int id=smoothed_updates[i].id-1;
      double x=smoothed_updates[i].S(state_x);
      double y=smoothed_updates[i].S(state_y);
      double tx=smoothed_updates[i].S(state_tx);
      double ty=smoothed_updates[i].S(state_ty);
      
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      double uwire=hits[id]->w;
      double v=hits[id]->s;
      
      // Get the aligment vector and error matrix for this layer
      unsigned int layer=hits[id]->wire->layer-1;
      DMatrix3x1 A=alignments[layer].A;
      DMatrix3x3 E=alignments[layer].E;
      double dx=A(kDx);
      double dy=A(kDy);
      double sindphi=sin(A(kDPhi));
      double cosdphi=cos(A(kDPhi));
      
      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;
      
      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi-dx*cosa+dy*sina;
      double vpred=x*sinpsi+y*cospsi-dx*sina-dy*cosa;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi-ty*cospsi;
      double du=upred-uwire;
      
      // Variables for angle of incidence with respect to the z-direction
      // in the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);
      
      // Transform from alignment vector coords to measurement coords
      double sinalpha_cosalpha=sinalpha*cosalpha;
      G_T(kDx)=-sina+tv*cosa*sinalpha_cosalpha;	
      G_T(kDy)=-cosa-tv*sina*sinalpha_cosalpha;
      
      double d=du*cosalpha;
      G_T(kDPhi)=-x*cospsi+y*sinpsi
	+tu*d*sinalpha-tv*(x*sinpsi+y*cospsi)*sinalpha_cosalpha
	-tu*tv*d*cosalpha*(cosalpha*cosalpha-sinalpha*sinalpha)
	;
      
      // G-matrix transpose
      G(kDx)=G_T(kDx);
      G(kDy)=G_T(kDy);
      G(kDPhi)=G_T(kDPhi);
   
      // Inverse of variance
      double InvV=1./(smoothed_updates[i].R+G*E*G_T);
      
      // Difference between measurement and projection
      double Mdiff=v-vpred+tv*d*sinalpha;
      
      // update the alignment vector and covariance
      DMatrix3x1 Ka=InvV*(E*G_T);
      DMatrix3x1 dA=Mdiff*Ka;
      DMatrix3x3 Etemp=E-Ka*G*E;
      if (Etemp(0,0)>0 && Etemp(1,1)>0 && Etemp(2,2)>0){
	alignments[layer].E=Etemp;
	alignments[layer].A=A+dA;	  
      }
    }
  }

  return NOERROR;
}
