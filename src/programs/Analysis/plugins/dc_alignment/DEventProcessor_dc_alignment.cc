// $Id$
//
//    File: DEventProcessor_dc_alignment.cc
// Created: Thu Oct 18 17:15:41 EDT 2012
// Creator: staylor (on Linux ifarm1102 2.6.18-274.3.1.el5 x86_64)
//

#include "DEventProcessor_dc_alignment.h"
using namespace jana;

#include <TROOT.h>

// Routine used to create our DEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_dc_alignment());
}
} // "C"


bool cdc_hit_cmp(const DCDCTrackHit *a,const DCDCTrackHit *b){
  
  return(a->wire->origin.Y()>b->wire->origin.Y());
}


bool bcal_cmp(const bcal_match_t &a,const bcal_match_t &b){
  return (a.match->y>b.match->y);
}


//------------------
// DEventProcessor_dc_alignment (Constructor)
//------------------
DEventProcessor_dc_alignment::DEventProcessor_dc_alignment()
{
  fdc_ptr = &fdc;
  
  pthread_mutex_init(&mutex, NULL);

}

//------------------
// ~DEventProcessor_dc_alignment (Destructor)
//------------------
DEventProcessor_dc_alignment::~DEventProcessor_dc_alignment()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_dc_alignment::init(void)
{
  mT0=0.;
  myevt=0;
 
  printf("Initializing..........\n");
 
  alignments.resize(24);
  for (unsigned int i=0;i<24;i++){
    alignments[i].A(kDx)=0.;
    alignments[i].A(kDy)=0.;
    alignments[i].A(kDPhi)=0.;
    alignments[i].E(kDx,kDx)=1.0;
    alignments[i].E(kDy,kDy)=1.0;
    alignments[i].E(kDPhi,kDPhi)=1.0;
  }
 	
  unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
			      135,135,146,146,158,158,170,170,182,182,197,197,
			      209,209};
  for (unsigned int i=0;i<28;i++){
    vector<cdc_align_t>tempvec;
    for (unsigned int j=0;j<numstraws[i];j++){
      cdc_align_t temp;
      temp.A=DMatrix4x1();
      temp.E=DMatrix4x4(1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.);
      tempvec.push_back(temp);
    }
    cdc_alignments.push_back(tempvec);
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
jerror_t DEventProcessor_dc_alignment::brun(JEventLoop *loop, int runnumber)
{	
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  dgeom  = dapp->GetDGeometry(runnumber);
  dgeom->GetFDCWires(fdcwires);

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_dz,endplate_rmin,endplate_rmax;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z+=0.5*endplate_dz;

  dapp->Lock();
    
  Hprob = (TH1F*)gROOT->FindObject("Hprob");
  if (!Hprob){
    Hprob=new TH1F("Hprob","Confidence level for time-based fit",100,0.0,1.); 
  } 
  Hprelimprob = (TH1F*)gROOT->FindObject("Hprelimprob");
  if (!Hprelimprob){
    Hprelimprob=new TH1F("Hprelimprob","Confidence level for prelimary fit",100,0.0,1.); 
  } 
  Hcdc_prob = (TH1F*)gROOT->FindObject("Hcdc_prob");
  if (!Hcdc_prob){
    Hcdc_prob=new TH1F("Hcdc_prob","Confidence level for time-based fit",100,0.0,1.); 
  } 
  Hcdc_prelimprob = (TH1F*)gROOT->FindObject("Hcdc_prelimprob");
  if (!Hcdc_prelimprob){
    Hcdc_prelimprob=new TH1F("Hcdc_prelimprob","Confidence level for prelimary fit",100,0.0,1.); 
  } 
  

  Hmatch = (TH1F*)gROOT->FindObject("Hmatch");
  if (!Hmatch){
    Hmatch=new TH1F("Hmatch","Segment matching distance",100,0.0,25.); 
  }
  Hbeta = (TH1F*)gROOT->FindObject("Hbeta");
  if (!Hbeta){
    Hbeta=new TH1F("Hbeta","Estimate for #beta",100,0.0,1.5); 
  }  
  HdEdx = (TH1F*)gROOT->FindObject("HdEdx");
  if (!HdEdx){
    HdEdx=new TH1F("HdEdx","Estimate for dE/dx",100,0.0,1e-5); 
  }
  HdEdx_vs_beta = (TH2F*)gROOT->FindObject("HdEdx_vs_beta");
  if (!HdEdx_vs_beta){
    HdEdx_vs_beta=new TH2F("HdEdx_vs_beta","dE/dx vs p/M",100,0,1.5,100,0.0,1e-5); 
  }

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
  Hcdcdrift_time=(TH2F*)gROOT->FindObject("Hcdcdrift_time");
  if (!Hcdcdrift_time){
    Hcdcdrift_time=new TH2F("Hcdcdrift_time",
			 "cdc doca vs drift time",801,-21,781,100,0,1);
  }  
  Hcdcres_vs_drift_time=(TH2F*)gROOT->FindObject("Hcdcres_vs_drift_time");
  if (!Hcdcres_vs_drift_time){
    Hcdcres_vs_drift_time=new TH2F("Hcdcres_vs_drift_time","cdc Residual vs drift time",800,-20,780,200,-0.1,0.1);
  } 
  Hdrift_time=(TH2F*)gROOT->FindObject("Hdrift_time");
  if (!Hdrift_time){
    Hdrift_time=new TH2F("Hdrift_time",
			 "doca vs drift time",201,-21,381,100,0,1);
  }  
  Hres_vs_drift_time=(TH2F*)gROOT->FindObject("Hres_vs_drift_time");
  if (!Hres_vs_drift_time){
    Hres_vs_drift_time=new TH2F("Hres_vs_drift_time","Residual vs drift time",320,-20,300,1000,-1,1);
  } 
  Hdv_vs_dE=(TH2F*)gROOT->FindObject("Hdv_vs_dE");
  if (!Hdv_vs_dE){
    Hdv_vs_dE=new TH2F("Hdv_vs_dE","dv vs energy dep",100,0,20e-6,200,-1,1);
  }

  Hbcalmatch=(TH2F*)gROOT->FindObject("Hbcalmatch");
  if (!Hbcalmatch){
    Hbcalmatch=new TH2F("Hbcalmatch","BCAL #deltar vs #deltaz",100,-25.,25.,
			100,0.,10.);
  }
  
  dapp->Unlock();

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_dc_alignment::erun(void)
{
 

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_dc_alignment::fini(void)
{

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_dc_alignment::evnt(JEventLoop *loop, int eventnumber){
  myevt++;
  
  // Get BCAL showers, FCAL showers and FDC space points
  vector<const DFCALShower*>fcalshowers;
  loop->Get(fcalshowers);
  vector<const DBCALShower*>bcalshowers;
  loop->Get(bcalshowers);
  vector<const DFDCPseudo*>pseudos;
  loop->Get(pseudos);
  vector<const DCDCTrackHit*>cdcs;
  loop->Get(cdcs);

  if (cdcs.size()>4 && bcalshowers.size()>0){
    vector<const DCDCTrackHit*>superlayers[5];
    for (unsigned int i=0;i<cdcs.size();i++){
      int ring=cdcs[i]->wire->ring;
      if (ring<=4) superlayers[0].push_back(cdcs[i]);
      else if (ring<=12) superlayers[1].push_back(cdcs[i]);
      else if (ring<=16) superlayers[2].push_back(cdcs[i]);
      else if (ring<=24) superlayers[3].push_back(cdcs[i]);
      else superlayers[4].push_back(cdcs[i]);
    }
    
    mMinTime=10000.;

    vector<cdc_segment_t>axial_segments;
    vector<cdc_segment_t>stereo_segments;
    FindSegments(superlayers[0],axial_segments);
    FindSegments(superlayers[1],stereo_segments);
    FindSegments(superlayers[2],axial_segments);
    FindSegments(superlayers[3],stereo_segments);
    FindSegments(superlayers[4],axial_segments);

    if (axial_segments.size()>0 && stereo_segments.size()>0){
      vector<cdc_track_t>tracks;
      LinkSegments(axial_segments,stereo_segments,tracks);

      for (unsigned int i=0;i<tracks.size();i++){
	double chi2x=0.,chi2y=0.;
	DMatrix4x1 S=GuessForStateVector(tracks[i],chi2x,chi2y);
	DMatrix4x1 Sbest=S;
	
	// Match to outer detectors 
	if (MatchOuterDetectors(fcalshowers,bcalshowers,S)){
	  // move S to the z-position where we found the match
	  S(state_x)+=mOuterZ*S(state_tx);
	  S(state_y)+=mOuterZ*S(state_ty);

	  // Add lists of stereo and axial hits associated with this track 
	  // and sort
	  vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
	  hits.insert(hits.end(),tracks[i].stereo_hits.begin(),tracks[i].stereo_hits.end());
	  sort(hits.begin(),hits.end(),cdc_hit_cmp);

	  // Run the Kalman Filter algorithm
	  DoFilter(S,hits);
	
	} // match outer detectors
      }
    }
  }

  if (pseudos.size()>4 && (fcalshowers.size()>0 || bcalshowers.size()>0)){
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
    
      mMinTime=1000.;
      mOuterTime=1000.;
      mOuterZ=0.;
      for (unsigned int m=0;m<hits.size();m++){
	if (hits[m]->time<mMinTime){
	  mMinTime=hits[m]->time;
	  mMinTimeID=m;
	}
      }

      // Perform preliminary line fits for current set of linked segments    
      double var_x,var_tx,cov_x_tx,var_y,var_ty,cov_y_ty,chi2x,chi2y;
      DMatrix4x1 S=FitLine(hits,var_x,cov_x_tx,var_tx,chi2x,var_y,cov_y_ty,
			 var_ty,chi2y); 
      double prelim_prob=TMath::Prob(chi2x,hits.size()-2);
      Hprelimprob->Fill(prelim_prob);

      vector<double>dE(hits.size());
      for (unsigned int j=0;j<dE.size();j++) dE[j]=hits[j]->dE;
      std::sort(dE.begin(),dE.end());
      double dEsum=0.;
      double numE=0.;
      for (unsigned int j=0;j<dE.size()/2;j++){
	dEsum+=dE[j];
	numE+=1.0;
      }
      dEsum/=numE*1.0*sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
      HdEdx->Fill(dEsum);

      // Match to outer detectors 
      if (MatchOuterDetectors(fcalshowers,bcalshowers,S)){
                
	// Run the Kalman Filter algorithm
	if (DoFilter(S,hits)==NOERROR){
	  HdEdx_vs_beta->Fill(mBeta,dEsum);
	}
      }	
    }
  }
   
  return NOERROR;
}

// Steering routine for the kalman filter
jerror_t 
DEventProcessor_dc_alignment::DoFilter(DMatrix4x1 &S,
				       vector<const DCDCTrackHit *>&hits){
  unsigned int numhits=hits.size();
  unsigned int maxindex=numhits-1;

  int NEVENTS=200000;
  double anneal_factor=pow(1e6,(double(NEVENTS-myevt))/(NEVENTS-1.));
  if (myevt>NEVENTS) anneal_factor=1.;  
  //anneal_factor=1.;

  // Create an initial reference trajectory
  deque<trajectory_t>trajectory;
  deque<trajectory_t>best_traj;
  // if (SetReferenceTrajectory(mOuterZ,S,trajectory,hits[maxindex])==NOERROR)
    {
    // State vector to store "best" values
    DMatrix4x1 Sbest;

    // Covariance matrix
    DMatrix4x4 C0,C,Cbest;
    C0(state_x,state_x)=C0(state_y,state_y)=1.;
    C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;
    
    vector<cdc_update_t>updates(hits.size());
    vector<cdc_update_t>best_updates;
    double chi2=1e16,chi2_old=1e16;
    unsigned int ndof=0,ndof_old=0;
    unsigned int iter=0;
    
    //printf("wirebased-----------\n");

    // Perform a wire-based pass
    for(iter=0;iter<20;iter++){
      chi2_old=chi2; 
      ndof_old=ndof;

      trajectory.clear();
      if (SetReferenceTrajectory(mOuterZ,S,trajectory,
				 hits[maxindex])!=NOERROR) break;

      C=C0;
      if (KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof)!=NOERROR)
	break;
	      
      //printf(">>>>>>chi2 %f ndof %d\n",chi2,ndof);

      if (fabs(chi2_old-chi2)<0.1 || chi2>chi2_old+1.0) break;  
      
      // Save the current state and covariance matrixes
      Cbest=C;
      Sbest=S;
      
      // run the smoother (opposite direction to filter)
      //Smooth(S,C,trajectory,updates);	      
    }
    if (iter>0){
      double prelimprob=TMath::Prob(chi2_old,ndof_old);
      Hcdc_prelimprob->Fill(prelimprob);

      //printf("cdc prob %f\n",prelimprob);

      if (prelimprob>0.01){ 
	// Perform a time-based pass
	S=Sbest;
	chi2=1e16;

	//printf("Timebased-----------\n");
	for (iter=0;iter<20;iter++){
	  chi2_old=chi2; 
	  ndof_old=ndof;
	  
	  trajectory.clear();
	  if (SetReferenceTrajectory(mOuterZ,S,trajectory,hits[maxindex])
	      ==NOERROR){
	    C=C0;
	    KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof,true);
  
	    //printf(">>>>>>chi2 %f ndof %d\n",chi2,ndof);
	    if (fabs(chi2-chi2_old)<0.1 || chi2>chi2_old+1.0 || ndof!=ndof_old) break;
	    
	    Sbest=S;
	    Cbest=C;
	    best_updates.assign(updates.begin(),updates.end());
	    best_traj.assign(trajectory.begin(),trajectory.end());
	  }
	  else break;
	}
	if (iter>0){
	  double prob=TMath::Prob(chi2_old,ndof_old);
	  Hcdc_prob->Fill(prob);

	  if (prob>0.01){
	    // run the smoother (opposite direction to filter)
	    vector<cdc_update_t>smoothed_updates(updates.size());
	    Smooth(Sbest,Cbest,best_traj,hits,best_updates,smoothed_updates);

	    for (unsigned int k=0;k<smoothed_updates.size();k++){
	      double tdrift=smoothed_updates[k].drift_time;
	      double d=smoothed_updates[k].doca;
	      double res=smoothed_updates[k].res;
	      Hcdcres_vs_drift_time->Fill(tdrift,res);
	      Hcdcdrift_time->Fill(tdrift,d);
	    }
	  }
	}
      }
    }
  }




  return NOERROR;
}


// Steering routine for the kalman filter
jerror_t 
DEventProcessor_dc_alignment::DoFilter(DMatrix4x1 &S,
				    vector<const DFDCPseudo *>&hits){
  unsigned int num_hits=hits.size();
  //  vector<strip_update_t>strip_updates(num_hits);
  //vector<strip_update_t>smoothed_strip_updates(num_hits);
  vector<update_t>updates(num_hits);
  vector<update_t>best_updates;
  vector<update_t>smoothed_updates(num_hits);

  int NEVENTS=200000;
  double anneal_factor=pow(1e6,(double(NEVENTS-myevt))/(NEVENTS-1.));
  if (myevt>NEVENTS) anneal_factor=1.;  
  //anneal_factor=1.;
  //anneal_factor=1e3;

  // Best guess for state vector at "vertex"
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;
  deque<trajectory_t>best_traj;
  // double start_z=hits[0]->wire->origin.z()-1.;
  S(state_x)+=endplate_z*S(state_tx);
  S(state_y)+=endplate_z*S(state_ty);
  //  SetReferenceTrajectory(endplate_z,S,trajectory,hits);
      
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

    trajectory.clear();
    if (SetReferenceTrajectory(endplate_z,S,trajectory,hits)!=NOERROR) break;
    C=C0;
    if (KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof)
	!=NOERROR) break;

    //printf("== event %d == iter %d =====chi2 %f ndof %d \n",myevt,iter,chi2,ndof);
    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1 || iter==ITER_MAX) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S;
    best_updates.assign(updates.begin(),updates.end());
    best_traj.assign(trajectory.begin(),trajectory.end());
    // run the smoother (opposite direction to filter)
    //Smooth(S,C,trajectory,hits,updates,smoothed_updates);
  }
      
  if (iter>1){
    double prob=TMath::Prob(chi2_old,ndof_old);
    Hprob->Fill(prob);
  
    //printf("prob %f\n",prob);

    if (prob>0.1){
      // run the smoother (opposite direction to filter)
      Smooth(Sbest,Cbest,best_traj,hits,best_updates,smoothed_updates);

      Hbeta->Fill(mBeta);
      for (unsigned int i=0;i<smoothed_updates.size();i++){
	unsigned int layer=hits[i]->wire->layer;
	Hures_vs_layer->Fill(layer,smoothed_updates[i].res(0));
	Hvres_vs_layer->Fill(layer,smoothed_updates[i].res(1));
	if (layer==1){
	  Hres_vs_drift_time->Fill(smoothed_updates[i].drift_time,
				   smoothed_updates[i].res(0));
	  Hdrift_time->Fill(smoothed_updates[i].drift_time,
			    smoothed_updates[i].doca);
 
	  Hdv_vs_dE->Fill(hits[i]->dE,smoothed_updates[i].res(1));
	}
      }

      FindOffsets(hits,smoothed_updates);
      
      for (unsigned int layer=0;layer<24;layer++){
	// Set up to fill tree
	double dxr=alignments[layer].A(kDx);
	double dyr=alignments[layer].A(kDy);  
	fdc.dPhi=180./M_PI*alignments[layer].A(kDPhi);
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
      return NOERROR;
    }
  }
   
  
  return VALUE_OUT_OF_RANGE;
}


// Link segments from package to package by doing straight-line projections
jerror_t
DEventProcessor_dc_alignment::LinkSegments(vector<segment_t>segments[4], 
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

// Find segments in cdc axial layers
jerror_t DEventProcessor_dc_alignment::FindSegments(vector<const DCDCTrackHit*>&hits,
						    vector<cdc_segment_t>&segments){

  vector<unsigned int>ring_boundaries;
  vector<bool>used_in_segment(hits.size());
  int last_ring=-1;
  for (unsigned int i=0;i<hits.size();i++){
    int ring=hits[i]->wire->ring;
    if (ring!=last_ring){
      ring_boundaries.push_back(i);
    }
    last_ring=ring;  
  }
  ring_boundaries.push_back(hits.size());

  unsigned int start=0;
  while (start<ring_boundaries.size()-1){
    for (unsigned int i=ring_boundaries[start];i<ring_boundaries[start+1];i++){
      if (used_in_segment[i]==false){
	used_in_segment[i]=true;
	
	// Current wire position 
	DVector3 pos=hits[i]->wire->origin;

	// Create list of nearest neighbors
	vector<const DCDCTrackHit *>neighbors;
	neighbors.push_back(hits[i]);	  
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<ring_boundaries.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=ring_boundaries[k];m<ring_boundaries[k+1];m++){
	    delta=(pos-hits[m]->wire->origin).Perp();
	    if (delta<delta_min && delta<CDC_MATCH_RADIUS){
	      delta_min=delta;
	      match=m;
	    }
	  }
	  // Hcdc_match->Fill(delta_min);
	  if (//match!=0 
		//&& 
	      used_in_segment[match]==false
	      ){
	    pos=hits[match]->wire->origin;
	    used_in_segment[match]=true;
	    neighbors.push_back(hits[match]);
	  }
	}
	
	if (neighbors.size()>1){
	  cdc_segment_t mysegment;
	  mysegment.matched=false;
	  mysegment.dir=neighbors[neighbors.size()-1]->wire->origin
	    -neighbors[0]->wire->origin;
	  mysegment.dir.SetMag(1.);
	  mysegment.hits=neighbors;
	  segments.push_back(mysegment);
	}
      }
    } // loop over start points in a ring
    
    // Look for a new ring to start looking for a segment
    while (start<ring_boundaries.size()-1){
      if (used_in_segment[ring_boundaries[start]]==false) break;
	start++;
    }
  }
  
  return NOERROR;
}


// Find segments by associating adjacent hits within a package together.
jerror_t DEventProcessor_dc_alignment::FindSegments(vector<const DFDCPseudo*>&points,
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
	  Hmatch->Fill(delta_min);
	  if (//match!=0 
	      //&& 
	      used[match]==false
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
DEventProcessor_dc_alignment::FitLine(vector<const DFDCPseudo*> &fdchits){
  double var_x,var_tx,cov_x_tx,chi2x;
  double var_y,var_ty,cov_y_ty,chi2y;
  
  return FitLine(fdchits,var_x,cov_x_tx,var_tx,chi2x,var_y,cov_y_ty,var_ty,
		 chi2y);
}


// Use linear regression on the hits to obtain a first guess for the state
// vector.  Method taken from Numerical Recipes in C.
DMatrix4x1 
DEventProcessor_dc_alignment::FitLine(vector<const DFDCPseudo*> &fdchits,
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

  double sig2v=0.04; // rough guess;

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
jerror_t DEventProcessor_dc_alignment::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
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
      smoothed_updates[id].t=updates[id].t;
      smoothed_updates[id].S=Ss;
      smoothed_updates[id].C=Cs;
      smoothed_updates[id].R=updates[id].R-updates[id].H*dC*updates[id].H_T;
    }
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }

  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}

// Kalman smoother 
jerror_t DEventProcessor_dc_alignment::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
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
      unsigned int first_id=trajectory[m].h_id-1;
      for (int k=trajectory[m].num_hits-1;k>=0;k--){
	unsigned int id=first_id+k;
	A=updates[id].C*JT*C.Invert();
	dC=A*(Cs-C)*A.Transpose();
	Ss=updates[id].S+A*(Ss-S);
	Cs=updates[id].C+dC;
      
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
	smoothed_updates[id].res(1)=v-vpred+tv*d*sinalpha;
	smoothed_updates[id].res(0)=(d>0?1.:-1.)*updates[id].drift-d;
	smoothed_updates[id].doca=fabs(d);
	
	smoothed_updates[id].drift=updates[id].drift;
	smoothed_updates[id].drift_time=updates[id].drift_time;
	smoothed_updates[id].S=Ss;
	smoothed_updates[id].C=Cs;
	smoothed_updates[id].R=updates[id].R-updates[id].H*dC*updates[id].H_T;
      }
    }
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  
  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}

// Kalman smoother 
jerror_t 
DEventProcessor_dc_alignment::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
				     deque<trajectory_t>&trajectory,
				     vector<const DCDCTrackHit *>&hits,
				     vector<cdc_update_t>&updates,
				     vector<cdc_update_t>&smoothed_updates
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
    else{
      unsigned int id=trajectory[m].h_id-1;
      A=updates[id].C*JT*C.Invert();
      dC=A*(Cs-C)*A.Transpose();
      Ss=updates[id].S+A*(Ss-S);
      Cs=updates[id].C+dC;
         
      // CDC index and wire position variables
      DVector3 origin=hits[id]->wire->origin;
      unsigned int ring=hits[id]->wire->ring-1;
      unsigned int straw=hits[id]->wire->straw-1;
      DVector3 dOrigin(cdc_alignments[ring][straw].A(dX),
		       cdc_alignments[ring][straw].A(dY),0.);
      origin+=dOrigin;
      DVector3 vhat=hits[id]->wire->udir;
      DVector3 dV(cdc_alignments[ring][straw].A(dVx),
		  cdc_alignments[ring][straw].A(dVy),0.);
      //dV.Print();
      vhat+=dV;
      vhat.SetZ((vhat.z()>0?1.:-1.)*sqrt(1.-vhat.Perp2()));

      // doca using smoothed state vector
      double d=FindDoca(trajectory[m].z,Ss,vhat,origin);
      smoothed_updates[id].doca=d;
      smoothed_updates[id].res=updates[id].drift-d;
 
      smoothed_updates[id].drift=updates[id].drift;
      smoothed_updates[id].drift_time=updates[id].drift_time;
      smoothed_updates[id].S=Ss;
      smoothed_updates[id].C=Cs;
      smoothed_updates[id].V=updates[id].V-updates[id].H*dC*updates[id].H_T;
            

      // Reset h_id for this position along the reference trajectory
      trajectory[m].h_id=0;
    }

    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  
  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}



// Perform Kalman Filter for the current trajectory
jerror_t 
DEventProcessor_dc_alignment::KalmanFilterCathodesOnly(double anneal_factor,
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
      double tv=tx*sinpsi+ty*cospsi;

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
	updates[my_id].t=trajectory[k].t;

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
  //chi2*=anneal_factor;

  ndof-=4;

  return NOERROR;
}

// Perform the Kalman Filter for the current set of cdc hits
jerror_t 
DEventProcessor_dc_alignment::KalmanFilter(double anneal_factor,
					   DMatrix4x1 &S,DMatrix4x4 &C,
					   vector<const DCDCTrackHit *>&hits,
					   deque<trajectory_t>&trajectory,
					   vector<cdc_update_t>&updates,
					   double &chi2,unsigned int &ndof,
					   bool timebased){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

  double doca2=0.;
 
  // CDC index and wire position variables
  unsigned int cdc_index=hits.size()-1;
  bool more_hits=true;
  DVector3 origin=hits[cdc_index]->wire->origin;
  unsigned int ring=hits[cdc_index]->wire->ring-1;
  unsigned int straw=hits[cdc_index]->wire->straw-1;
  DVector3 dOrigin(cdc_alignments[ring][straw].A(dX),
		   cdc_alignments[ring][straw].A(dY),0.);
  origin+=dOrigin;
  DVector3 vhat=hits[cdc_index]->wire->udir;
  DVector3 dV(cdc_alignments[ring][straw].A(dVx),
	      cdc_alignments[ring][straw].A(dVy),0.);
  //dV.Print();
  vhat+=dV;
  vhat.SetZ((vhat.z()>0?1.:-1.)*sqrt(1.-vhat.Perp2()));
  double z0=origin.z();
  double vz=vhat.z();
  DVector3 wirepos=origin+((trajectory[0].z-z0)/vz)*vhat;

  /// compute initial doca^2 to first wire
  double dx=S(state_x)-wirepos.x();
  double dy=S(state_y)-wirepos.y();
  double old_doca2=dx*dx+dy*dy;
  
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
    
    // Position along wire
    wirepos=origin+((trajectory[k].z-z0)/vz)*vhat;

    // New doca^2
    dx=S(state_x)-wirepos.x();
    dy=S(state_y)-wirepos.y();
    doca2=dx*dx+dy*dy;

    if (doca2>old_doca2 && more_hits){
      // zero-position and direction of line describing particle trajectory
      double tx=S(state_tx),ty=S(state_ty);
      DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
      DVector3 uhat(tx,ty,1.);
      uhat.SetMag(1.); 
      double ux=uhat.x(),uy=uhat.y(),uz=uhat.z();

      // Find the true doca to the wire
      DVector3 diff=pos0-origin;
      double dx0=diff.x(),dy0=diff.y(),dz0=diff.z();
      double vhat_dot_diff=diff.Dot(vhat);
      double uhat_dot_diff=diff.Dot(uhat);
      double uhat_dot_vhat=uhat.Dot(vhat);
      double D=1.-uhat_dot_vhat*uhat_dot_vhat;
      double N=uhat_dot_vhat*vhat_dot_diff-uhat_dot_diff;
      double N1=vhat_dot_diff-uhat_dot_vhat*uhat_dot_diff;
      double scale=1./D;
      double s=scale*N;
      double t=scale*N1;
      diff+=s*uhat-t*vhat;
      double d=diff.Mag();

      // The next measurement and its variance
      double tdrift=hits[cdc_index]->tdrift-mT0-trajectory[k].t;
      double dmeas=0.39;
      double V=1.0*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor
      if (timebased){
	dmeas=cdc_drift_distance(tdrift);
	V=anneal_factor*cdc_variance(tdrift);
      }

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      H(state_x)=H_T(state_x)=diff.x()*one_over_d;
      H(state_y)=H_T(state_y)=diff.y()*one_over_d;

      double one_plus_tx2_plus_ty2=1.+tx*tx+ty*ty;
      double scale2=1./(one_plus_tx2_plus_ty2*sqrt(one_plus_tx2_plus_ty2));

      double duxdtx=scale2*(1.+ty*ty);
      double duydtx=-scale2*tx*ty;
      double duzdtx=-scale2*tx;
      double duxdty=duydtx;
      double duydty=scale2*(1.+tx*tx);
      double duzdty=-scale2*ty;

      double vx=vhat.x(),vy=vhat.y();    
      double dDdux=-2.*uhat_dot_vhat*vx;
      double dDduy=-2.*uhat_dot_vhat*vy;
      double dDduz=-2.*uhat_dot_vhat*vz;

      double dNdux=(vx*vx-1.)*dx0+vx*vy*dy0+vx*vz*dz0;
      double dNduy=vy*vx*dx0+(vy*vy-1.)*dy0+vy*vz*dz0;
      double dNduz=vz*vx*dx0+vz*vy*dy0+(vz*vz-1.)*dz0;
 
      double ratio=N*scale;
      double dsdux=scale*(dNdux-ratio*dDdux);
      double dsduy=scale*(dNduy-ratio*dDduy);
      double dsduz=scale*(dNduz-ratio*dDduz);
      double dsdtx=dsdux*duxdtx+dsduy*duydtx+dsduz*duzdtx;
      double dsdty=dsdux*duxdty+dsduy*duydty+dsduz*duzdty;

      double dN1dux=-(2.*ux*vx+uy*vy+uz*vz)*dx0-vx*uy*dy0-vx*uz*dz0;
      double dN1duy=-vy*ux*dx0-(ux*vx+2.*uy*vy+uz*vz)*dy0-vy*uz*dz0;
      double dN1duz=-vz*ux*dx0-uz*uy*dy0-(2.*vz*uz+vy*uy+vx*ux)*dz0;

      double ratio1=N1*scale;
      double dtdux=scale*(dN1dux-ratio1*dDdux); 
      double dtduy=scale*(dN1duy-ratio1*dDduy);
      double dtduz=scale*(dN1duz-ratio1*dDduz);
      double dtdtx=dtdux*duxdtx+dtduy*duydtx+dtduz*duzdtx;  
      double dtdty=dtdux*duxdty+dtduy*duydty+dtduz*duzdty;
  
      double dx1=diff.x(),dy1=diff.y(),dz1=diff.z();
      H(state_tx)=H_T(state_tx)=one_over_d*(dx1*(ux*dsdtx+duxdtx*s-vx*dtdtx)
					    +dy1*(uy*dsdtx+duydtx*s-vy*dtdtx)
					    +dz1*(uz*dsdtx+duzdtx*s-vz*dtdtx));
      H(state_ty)=H_T(state_ty)=one_over_d*(dx1*(ux*dsdty+duxdty*s-vx*dtdty)
					    +dy1*(uy*dsdty+duydty*s-vy*dtdty)
					    +dz1*(uz*dsdty+duzdty*s-vz*dtdty));

      // Matrices to rotate alignment error matrix into measurement space
      DMatrix1x4 G;
      DMatrix4x1 G_T;

      double dNdvx=ux*vhat_dot_diff+uhat_dot_vhat*dx0;
      double dDdvx=-2.*uhat_dot_vhat*ux;
      double dsdvx=scale*(dNdvx-ratio*dDdvx);

      double dNdvy=uy*vhat_dot_diff+uhat_dot_vhat*dy0;
      double dDdvy=-2.*uhat_dot_vhat*uy;
      double dsdvy=scale*(dNdvy-ratio*dDdvy);

      double dNdvz=uz*vhat_dot_diff+uhat_dot_vhat*dz0;
      double dDdvz=-2.*uhat_dot_vhat*uz;
      double dsdvz=scale*(dNdvz-ratio*dDdvz);

      double dsddx=scale*(ux-vx*uhat_dot_vhat);
      double dsddvx=dsdvx-(vx/vz)*dsdvz;
      double dsddy=scale*(uy-vy*uhat_dot_vhat);
      double dsddvy=dsdvy-(vy/vz)*dsdvz;

      double dN1dvx=dx0-ux*uhat_dot_diff;
      double dtdvx=scale*(dN1dvx-ratio1*dDdvx);

      double dN1dvy=dy0-uy*uhat_dot_diff;
      double dtdvy=scale*(dN1dvy-ratio1*dDdvy);

      double dN1dvz=dz0-uz*uhat_dot_diff;
      double dtdvz=scale*(dN1dvz-ratio1*dDdvz);

      double dtddx=scale*(uhat_dot_vhat*ux-vx);
      double dtddvx=dtdvx-(vx/vz)*dtdvz;
      double dtddy=scale*(uhat_dot_vhat*uy-vy);
      double dtddvy=dtdvy-(vy/vz)*dtdvz;

      G_T(dX)=one_over_d*(dx1*(-1.+ux*dsddx-vx*dtddx)+dy1*(uy*dsddx-vy*dtddx)
			  +dz1*(uz*dsddx-vz*dtddx));
      G_T(dY)=one_over_d*(dx1*(ux*dsddy-vx*dtddy)+dy1*(-1.+uy*dsddy-vy*dtddy)
			  +dz1*(uz*dsddy-vz*dtddy));
      G_T(dVx)=one_over_d*(dx1*(ux*dsddvx-t-vx*dtddvx)+dy1*(uy*dsddvx-vy*dtddvx)
			   +dz1*(uz*dsddvx-(vx/vz)*t-vz*dtddvx));
      G_T(dVy)=one_over_d*(dx1*(ux*dsddvy-vx*dtddvy)+dy1*(uy*dsddvy-t-vy*dtddvy)
			   +dz1*(uz*dsddvy-(vy/vz)*t-vz*dtddvy));
      G(dX)=G_T(dX);
      G(dY)=G_T(dY);
      G(dVx)=G_T(dVx);
      G(dVy)=G_T(dVy);			   
      	
      // inverse of variance including prediction
      DMatrix4x4 E=cdc_alignments[ring][straw].E;
      V+=G*E*G_T;
      double InvV=1./(V+H*C*H_T);

      //printf("dV %f\n",(G*E*G_T));
      
      // Compute Kalman gain matrix
      K=InvV*(C*H_T);

      // Update state vector covariance matrix
      DMatrix4x4 Ctest=C-K*(H*C);

      // Check that Ctest is positive definite
      if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0){
	C=Ctest;

	// Update the state vector 
	//S=S+res*K;
	S+=res*K;

	// Compute new residual 
	d=FindDoca(trajectory[k].z,S,vhat,origin);
	res=dmeas-d;
	
	// Update chi2 for this segment
	chi2+=res*res/(V-H*C*H_T);
	ndof++;	

      }
      else return VALUE_OUT_OF_RANGE;

      updates[cdc_index].S=S;
      updates[cdc_index].C=C;
      updates[cdc_index].drift=dmeas;
      updates[cdc_index].drift_time=tdrift;
      updates[cdc_index].doca=d;
      updates[cdc_index].res=res;
      updates[cdc_index].V=V;
      updates[cdc_index].H_T=H_T;
      updates[cdc_index].H=H;

      trajectory[k].h_id=cdc_index+1;

      // move to next cdc hit
      if (cdc_index>0){
	cdc_index--;

	//New wire position
	origin=hits[cdc_index]->wire->origin;
	vhat=hits[cdc_index]->wire->udir;
	ring=hits[cdc_index]->wire->ring-1;
	straw=hits[cdc_index]->wire->straw-1;
	dOrigin.SetXYZ(cdc_alignments[ring][straw].A(dX),
		       cdc_alignments[ring][straw].A(dY),0.);
	origin+=dOrigin;
	dV.SetXYZ(cdc_alignments[ring][straw].A(dVx),
		  cdc_alignments[ring][straw].A(dVy),0.);
	vhat+=dV;
	vhat.SetZ((vhat.z()>0?1.:-1.)*sqrt(1.-vhat.Perp2()));

	z0=origin.z();
	vz=vhat.z();	
	wirepos=origin+((trajectory[k].z-z0)/vz)*vhat;
	
	// New doca^2
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca2=dx*dx+dy*dy;	
	
      }
      else more_hits=false;
    }
  
    old_doca2=doca2;
  }

  ndof-=4;

  return NOERROR;
}

// Perform Kalman Filter for the current trajectory
jerror_t 
DEventProcessor_dc_alignment::KalmanFilter(double anneal_factor,
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
      
      // Get the alignment vector and error matrix for this layer
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
      double tv=tx*sinpsi+ty*cospsi;

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
	updates[my_id].t=trajectory[k].t;

	double du=upred-uwire;
	double d=du*cosalpha;

	if (my_id==mMinTimeID){ 
	  double par[4]={0.23682,6.85469,1141.81,-859.02};
	  double d_sq=d*d;    
	  double t_from_d=par[0]+par[1]*fabs(d)+par[2]*d_sq+par[3]*fabs(d)*d_sq;
	  mBeta=(mOuterZ-trajectory[k].z)/29.98*sqrt(1.+tx*tx+ty*ty)
	    /(mOuterTime-(hits[my_id]->time-t_from_d));
	}

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

	DMatrix2x2 Vtemp=V+G*E*G_T;

	// Variance for this hit
	InvV=(Vtemp+H*C*H_T).Invert();
	
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
	DMatrix2x2 RC=Vtemp-H*C*H_T;
;
	updates[my_id].R=RC;
	
	chi2+=RC.Chi2(R);
	ndof+=2;
      }

    }

  }
  //  chi2*=anneal_factor;
  ndof-=4;

  return NOERROR;
}

//Reference trajectory for the track for cdc tracks
jerror_t DEventProcessor_dc_alignment
::SetReferenceTrajectory(double z,DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 const DCDCTrackHit *last_cdc){ 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double ds=1.0;
  double dz=(S(state_ty)>0.?-1.:1.)*ds/sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  double t=0.;
  trajectory_t temp;

  //y-position after which we cut off the loop
  double min_y=last_cdc->wire->origin.y()-5.;
  do{
    double newz=z+dz;
    temp.Skk=Zero4x1;
    temp.Ckk=Zero4x4;
    temp.h_id=0;	  
    temp.z=newz;
    temp.J=J;
    temp.J(state_x,state_tx)=-dz;
    temp.J(state_y,state_ty)=-dz;
    // Flight time: assume particle is moving at the speed of light
    temp.t=(t+=ds/29.98);
    //propagate the state to the next z position
    temp.S(state_x)=S(state_x)+S(state_tx)*dz;
    temp.S(state_y)=S(state_y)+S(state_ty)*dz;
    temp.S(state_tx)=S(state_tx);
    temp.S(state_ty)=S(state_ty);
    S=temp.S;
    trajectory.push_front(temp);
    
    z=newz;
  }while (S(state_y)>min_y);

  if (trajectory.size()<2) return UNRECOVERABLE_ERROR;
  if (false)
    {
    printf("Trajectory:\n");
    for (unsigned int i=0;i<trajectory.size();i++){
      printf(" x %f y %f z %f\n",trajectory[i].S(state_x),
	     trajectory[i].S(state_y),trajectory[i].z); 
    }
  }




  return NOERROR;
}


// Reference trajectory for the track
jerror_t DEventProcessor_dc_alignment
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
      /*
      trajectory_t temp;
      temp.J=J;
      temp.S=trajectory[0].S;
      temp.t=trajectory[0].t;
      temp.h_id=i+1;
      temp.z=trajectory[0].z;
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;
      trajectory.push_front(temp); 
      */
      trajectory[0].num_hits++;
      continue;
    }
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      
      trajectory_t temp;
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
    printf(" x %f y %f z %f first hit %d num in layer %d\n",trajectory[i].S(state_x),
	   trajectory[i].S(state_y),trajectory[i].z,trajectory[i].h_id,
	   trajectory[i].num_hits); 
    }
  }

  return NOERROR;
}

// Crude approximation for the variance in drift distance due to smearing
double DEventProcessor_dc_alignment::GetDriftVariance(double t){
  if (t<0) t=0;
  else if (t>110.) t=110.;
  double sigma=0.01639/sqrt(t+1.)+5.405e-3+4.936e-4*exp(0.09654*(t-66.86)); 
  return sigma*sigma;
}

#define FDC_T0_OFFSET 20.
// convert time to distance for the fdc
double DEventProcessor_dc_alignment::GetDriftDistance(double t){
  if (t<0.) return 0.;
  double d=0.0268*sqrt(t)/*-3.051e-4*/+7.438e-4*t;
  if (d>0.5) d=0.5;
  return d;
}

jerror_t 
DEventProcessor_dc_alignment::FindOffsets(vector<const DFDCPseudo *>&hits,
				      vector<update_t>smoothed_updates){
  DMatrix2x3 G;//matrix relating alignment vector to measurement coords
  DMatrix3x2 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();


  for (unsigned int i=0;i<num_hits;i++){ 
    double x=smoothed_updates[i].S(state_x);
    double y=smoothed_updates[i].S(state_y);
    double tx=smoothed_updates[i].S(state_tx);
    double ty=smoothed_updates[i].S(state_ty);
    
    double cosa=hits[i]->wire->udir.y();
    double sina=hits[i]->wire->udir.x();
    double uwire=hits[i]->w;
    //double v=hits[i]->s;
      
    // Get the aligment vector and error matrix for this layer
    unsigned int layer=hits[i]->wire->layer-1;
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
    double tu=tx*cospsi-ty*sinpsi;
    double tv=tx*sinpsi-ty*cospsi;
    double du=upred-uwire;
    
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
    
    // update the alignment vector and covariance
    DMatrix3x2 Ka=(E*G_T)*InvV;
    DMatrix3x1 dA=Ka*smoothed_updates[i].res;
    DMatrix3x3 Etemp=E-Ka*G*E;
    if (Etemp(0,0)>0 && Etemp(1,1)>0 && Etemp(2,2)>0){
      alignments[layer].E=Etemp;
      alignments[layer].A=A+dA;	  
    }
    else {
      printf("-------t= %f\n",smoothed_updates[i].drift_time);
      E.Print();
      Etemp.Print();
      smoothed_updates[i].R.Print();
    }
  }

  return NOERROR;
}

jerror_t 
DEventProcessor_dc_alignment::FindOffsets(vector<const DFDCPseudo *>&hits,
			       vector<strip_update_t>smoothed_updates){
  DMatrix1x3 G;//matrix relating alignment vector to measurement coords
  DMatrix3x1 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();

  for (unsigned int i=0;i<num_hits;i++){
    double x=smoothed_updates[i].S(state_x);
    double y=smoothed_updates[i].S(state_y);
    double tx=smoothed_updates[i].S(state_tx);
    double ty=smoothed_updates[i].S(state_ty);
    
    double cosa=hits[i]->wire->udir.y();
    double sina=hits[i]->wire->udir.x();
    double uwire=hits[i]->w;
    double v=hits[i]->s;
    
    // Get the aligment vector and error matrix for this layer
    unsigned int layer=hits[i]->wire->layer-1;
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

  return NOERROR;
}

// Routine to match preliminary track to outer detectors
bool DEventProcessor_dc_alignment::MatchOuterDetectors(vector<const DFCALShower *>&fcalshowers,
		      vector<const DBCALShower *>&bcalshowers,
		      const DMatrix4x1 &S){
    // First match to FCAL 
  double dz=0.;
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
      mOuterTime=fcalshowers[i]->getTime();
      mOuterZ=fcal_z;
    }
    
  }
  
  if (drmin<4.){ 
    //compute tangent of dip angle and related angular quantities
    double tanl=1./sqrt(S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));  
    double sinl=sin(atan(tanl));
    
    // Estimate for t0 at the beginning of track assuming particle is 
    // moving at the speed of light
    mT0=mOuterTime-dz/(29.98*sinl);

    //printf("t %f T0 %f\n",mOuterTime,mT0);

    return true;
  }
  else{
    // Match to BCAL
 
    // zero-position and direction of line describing particle trajectory
    DVector3 pos0(S(state_x),S(state_y),0.);
    DVector3 vhat(S(state_tx),S(state_ty),1.);
    vhat.SetMag(1.);
    
    // Keep list of matches
    vector<bcal_match_t>matching_bcals;

    // loop over the showers
    for (unsigned int i=0;i<bcalshowers.size();i++){
      // Match in x and y
      DVector3 pos1(bcalshowers[i]->x,bcalshowers[i]->y,bcalshowers[i]->z);
      DVector3 uhat(0.,0.,1.);
      DVector3 diff=pos1-pos0;
      double vhat_dot_uhat=vhat.Dot(uhat);
      double scale=1./(1.-vhat_dot_uhat*vhat_dot_uhat);
      double s=scale*(vhat_dot_uhat*diff.Dot(vhat)-diff.Dot(uhat));
      double t=scale*(diff.Dot(vhat)-vhat_dot_uhat*diff.Dot(uhat));
      double dr=(diff+s*uhat-t*vhat).Mag();
      
      Hbcalmatch->Fill(s,dr);
      if (dr<3.0){
	bcal_match_t temp;
	temp.ztrack=bcalshowers[i]->z+s;
	temp.match=bcalshowers[i];
	matching_bcals.push_back(temp);
      }	
    }
    if (matching_bcals.size()>0){
      sort(matching_bcals.begin(),matching_bcals.end(),bcal_cmp);
 
      mT0=matching_bcals[0].match->t;
      mOuterZ=matching_bcals[0].ztrack;

      return true;
    }
  }
 
  return false;
}

jerror_t 
DEventProcessor_dc_alignment::LinkSegments(vector<cdc_segment_t>&axial_segments,
					   vector<cdc_segment_t>&stereo_segments,
					   vector<cdc_track_t>&LinkedSegments){

  unsigned int num_axial=axial_segments.size();
  for (unsigned int i=0;i<num_axial-1;i++){  
    if (axial_segments[i].matched==false){
      cdc_track_t mytrack;
      mytrack.axial_hits=axial_segments[i].hits;

      DVector3 pos0=axial_segments[i].hits[0]->wire->origin;
      DVector3 vhat=axial_segments[i].dir;

      for (unsigned int j=i+1;j<num_axial;j++){
	if (axial_segments[j].matched==false){
	  DVector3 pos1=axial_segments[j].hits[0]->wire->origin;
	  DVector3 dir1=axial_segments[j].hits[0]->wire->udir;
	  DVector3 diff=pos1-pos0;
	  double s=diff.Dot(vhat);
	  double d=(diff-s*vhat).Mag();
	  if (d<CDC_MATCH_RADIUS){
	    axial_segments[j].matched=true;	   
	    mytrack.axial_hits.insert(mytrack.axial_hits.end(),
				  axial_segments[j].hits.begin(),
				  axial_segments[j].hits.end());
	    sort(mytrack.axial_hits.begin(),mytrack.axial_hits.end(),
		 cdc_hit_cmp);
	    
	    vhat=mytrack.axial_hits[mytrack.axial_hits.size()-1]->wire->origin
	      -mytrack.axial_hits[0]->wire->origin;
	    vhat.SetMag(1.);
	  }
	}
      }
    
      // Now try to associate stereo hits with this track
      unsigned int num_stereo=0;
      pos0=mytrack.axial_hits[0]->wire->origin;
      for (unsigned int j=0;j<stereo_segments.size();j++){
	if (stereo_segments[j].matched==false){
	  DVector3 pos1=stereo_segments[j].hits[0]->wire->origin;
	  DVector3 uhat=stereo_segments[j].hits[0]->wire->udir;
	  DVector3 diff=pos1-pos0;
	  double vhat_dot_uhat=vhat.Dot(uhat);
	  double scale=1./(1.-vhat_dot_uhat*vhat_dot_uhat);
	  double s=scale*(vhat_dot_uhat*diff.Dot(vhat)-diff.Dot(uhat));
	  double t=scale*(diff.Dot(vhat)-vhat_dot_uhat*diff.Dot(uhat));
	  double d=(diff+s*uhat-t*vhat).Mag();
	  if (d<CDC_MATCH_RADIUS){
	    num_stereo++;
	    stereo_segments[j].matched=true;
	    mytrack.stereo_hits.insert(mytrack.stereo_hits.end(),
				       stereo_segments[j].hits.begin(),
				       stereo_segments[j].hits.end());
	    sort(mytrack.stereo_hits.begin(),mytrack.stereo_hits.end(),
		 cdc_hit_cmp);
	  }
	}
      }
      if (num_stereo>1){
	mytrack.dir=vhat;
	LinkedSegments.push_back(mytrack);
      }
    }
  }

  return NOERROR;
}


// Compute initial guess for state vector (x,y,tx,ty) for a track in the CDC
// using two stereo wires
DMatrix4x1 
DEventProcessor_dc_alignment::GuessForStateVector(cdc_track_t &track,
						  double &chi2x,double &chi2y){

  // Parameters for line in x-y plane
  double vx=track.dir.x();
  double vy=track.dir.y();
  DVector3 pos0=track.axial_hits[0]->wire->origin;

  // Intersection of line in xy-plane with first stereo straw
  DVector3 origin_s=track.stereo_hits[0]->wire->origin;
  DVector3 dir_s=track.stereo_hits[0]->wire->udir;
  double ux_s=dir_s.x();
  double uy_s=dir_s.y();
  double dx=pos0.x()-origin_s.x();
  double dy=pos0.y()-origin_s.y();
  double s=(dx*vy-dy*vx)/(ux_s*vy-uy_s*vx);
  DVector3 pos1=origin_s+s*dir_s;

  // Intersection of line in xy-plane with last stereo straw
  unsigned int last_index=track.stereo_hits.size()-1;
  origin_s=track.stereo_hits[last_index]->wire->origin;
  dir_s=track.stereo_hits[last_index]->wire->udir;
  ux_s=dir_s.x();
  uy_s=dir_s.y();
  dx=pos0.x()-origin_s.x();
  dy=pos0.y()-origin_s.y();
  s=(dx*vy-dy*vx)/(ux_s*vy-uy_s*vx);
  DVector3 pos2=origin_s+s*dir_s;

  // Some code related to the drift time that may or may not be useful later
  //  double dt=hits[k]->tdrift-min_time;
  //  double d=0.019+0.0279*sqrt(dt);
  //  double scale=sqrt(vx*vx+vy*vy)/(ux_s*vy-uy_s*vx);
  //  double splus=s+d*scale;
  //  double sminus=s-d*scale;
  //  DVector3 pos_plus=origin_s+splus*dir_s;
  //  DVector3 pos_minus=origin_s+sminus*dir_s;

  // Estimate slopes and intercepts in xz and yz planes
  double z1=pos1.Z();
  double z2=pos2.Z();
  double delta_z=z1-z2;
  double x1=pos1.x();
  double x2=pos2.x();
  double x_slope=(x1-x2)/delta_z;
  double x_intercept=(x2*z1-x1*z2)/delta_z;  
  double y1=pos1.y();
  double y2=pos2.y();
  double y_slope=(y1-y2)/delta_z;
  double y_intercept=(y2*z1-y1*z2)/delta_z;

  //double norm_scale=1./sqrt(1.+x_slope*x_slope+y_slope*y_slope);
  
  return DMatrix4x1(x_intercept,y_intercept,x_slope,y_slope);
}

// Compute distance of closest approach between two lines
double DEventProcessor_dc_alignment::FindDoca(double z,const DMatrix4x1 &S,
					      const DVector3 &vhat,
					      const DVector3 &origin){
  DVector3 pos(S(state_x),S(state_y),z);
  DVector3 diff=pos-origin;
  
  DVector3 uhat(S(state_tx),S(state_ty),1.);
  uhat.SetMag(1.); 

  double vhat_dot_diff=diff.Dot(vhat);
  double uhat_dot_diff=diff.Dot(uhat);
  double uhat_dot_vhat=uhat.Dot(vhat);
  double D=1.-uhat_dot_vhat*uhat_dot_vhat;
  double N=uhat_dot_vhat*vhat_dot_diff-uhat_dot_diff;
  double N1=vhat_dot_diff-uhat_dot_vhat*uhat_dot_diff;
  double scale=1./D;
  double s=scale*N;
  double t=scale*N1;
  
  diff+=s*uhat-t*vhat;
  return diff.Mag();
}
