// $Id$
//
//    File: DEventProcessor_dc_alignment.cc
// Created: Thu Oct 18 17:15:41 EDT 2012
// Creator: staylor (on Linux ifarm1102 2.6.18-274.3.1.el5 x86_64)
//

#include "DEventProcessor_dc_alignment.h"
using namespace jana;

#include <TROOT.h>
#include <TCanvas.h>
#include <TPolyLine.h>

#define MAX_STEPS 1000

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

bool fdc_pseudo_cmp(const DFDCPseudo *a,const DFDCPseudo *b){
  return (a->wire->origin.z()<b->wire->origin.z());
}


bool bcal_cmp(const bcal_match_t &a,const bcal_match_t &b){
  return (a.match->y>b.match->y);
}

// Locate a position in vector xx given x
unsigned int DEventProcessor_dc_alignment::locate(vector<double>&xx,double x){
  int n=xx.size();
  if (x==xx[0]) return 0;
  else if (x==xx[n-1]) return n-2;

  int jl=-1;
  int ju=n;
  int ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    int jm=(ju+jl)>>1;
    if ( (x>=xx[jm])==ascnd)
      jl=jm;
    else
      ju=jm;
  } 
  return jl;
}


// Convert time to distance for the cdc
double DEventProcessor_dc_alignment::cdc_drift_distance(double t){
  double d=0.;
  if (t>cdc_drift_table[cdc_drift_table.size()-1]) return 0.78;
  if (t>0){
    unsigned int index=0;
    index=locate(cdc_drift_table,t);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(t-cdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 
  }
  return d;
}


// Interpolate on a table to convert time to distance for the fdc
double DEventProcessor_dc_alignment::fdc_drift_distance(double t){
  double d=0.;
  if (t>fdc_drift_table[fdc_drift_table.size()-1]) return 0.5;
  if (t>0){
    unsigned int index=0;
    index=locate(fdc_drift_table,t);
    double dt=fdc_drift_table[index+1]-fdc_drift_table[index];
    double frac=(t-fdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 
  }
  return d;
}



//------------------
// DEventProcessor_dc_alignment (Constructor)
//------------------
DEventProcessor_dc_alignment::DEventProcessor_dc_alignment()
{
  fdc_ptr = &fdc;
  fdc_c_ptr= &fdc_c;
  cdc_ptr = &cdc;
  
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
  one_over_zrange=1./150.;
 
  printf("Initializing..........\n");
 
  RUN_BENCHMARK=false;
  gPARMS->SetDefaultParameter("DCALIGN:RUN_BENCHMARK",RUN_BENCHMARK);
  USE_BCAL=false; 
  gPARMS->SetDefaultParameter("DCALIGN:USE_BCAL", USE_BCAL); 
  USE_FCAL=false; 
  gPARMS->SetDefaultParameter("DCALIGN:USE_FCAL", USE_FCAL);
  COSMICS=false; 
  gPARMS->SetDefaultParameter("DCALIGN:COSMICS", COSMICS);
  USE_DRIFT_TIMES=false;
  gPARMS->SetDefaultParameter("DCALIGN:USE_DRIFT_TIMES",USE_DRIFT_TIMES);
  READ_ANODE_FILE=false;
  gPARMS->SetDefaultParameter("DCALIGN:READ_ANODE_FILE",READ_ANODE_FILE);
  READ_CATHODE_FILE=false;
  gPARMS->SetDefaultParameter("DCALIGN:READ_CATHODE_FILE",READ_CATHODE_FILE);
  ALIGN_WIRE_PLANES=true;
  gPARMS->SetDefaultParameter("DCALIGN:ALIGN_WIRE_PLANES",ALIGN_WIRE_PLANES);
  FILL_TREE=false;
  gPARMS->SetDefaultParameter("DCALIGN:FILL_TREE",FILL_TREE);
  MIN_PSEUDOS=12;
  gPARMS->SetDefaultParameter("DCALIGN:MIN_PSEUDOS",MIN_PSEUDOS);
  MIN_INTERSECTIONS=10;
  gPARMS->SetDefaultParameter("DCALIGN:MIN_INTERSECTIONS",MIN_INTERSECTIONS);
 
  fdc_alignments.resize(24);
  for (unsigned int i=0;i<24;i++){
    fdc_alignments[i].A=DMatrix2x1();
    if (RUN_BENCHMARK==false){
      fdc_alignments[i].E=DMatrix2x2(0.000001,0.,0.,0.0001);
    }
    else{
      fdc_alignments[i].E=DMatrix2x2();
    }
  }
  fdc_cathode_alignments.resize(24);
  for (unsigned int i=0;i<24;i++){
    double var=0.0001;
    double var_phi=0.000001;
    if (RUN_BENCHMARK==false){
      fdc_cathode_alignments[i].E=DMatrix4x4(var_phi,0.,0.,0., 0.,var,0.,0., 
					   0.,0.,var_phi,0., 0.,0.,0.,var);
    }
    else{
      fdc_cathode_alignments[i].E=DMatrix4x4();
    }
    fdc_cathode_alignments[i].A=DMatrix4x1();
  }
  

  if (READ_ANODE_FILE){
    ifstream fdcfile("fdc_alignment.dat");
    // Skip first line, used to identify columns in file
    //char sdummy[40];
    //    fdcfile.getline(sdummy,40);
    // loop over remaining entries
    for (unsigned int i=0;i<24;i++){ 
      double du,dphi,dz;

      fdcfile >> dphi;
      fdcfile >> du;
      fdcfile >> dz;
      
      fdc_alignments[i].A(kU)=du;
      fdc_alignments[i].A(kPhiU)=dphi;
    }
    fdcfile.close();
  }
  
  if (READ_CATHODE_FILE){
    ifstream fdcfile("fdc_cathode_alignment.dat");
    // Skip first line, used to identify columns in file
    //char sdummy[40];
    //    fdcfile.getline(sdummy,40);
    // loop over remaining entries
    for (unsigned int i=0;i<24;i++){ 
      double du,dphiu,dv,dphiv;

      fdcfile >> dphiu;
      fdcfile >> du;
      fdcfile >> dphiv;
      fdcfile >> dv;
      
      fdc_cathode_alignments[i].A(kU)=du;
      fdc_cathode_alignments[i].A(kPhiU)=dphiu;  
      fdc_cathode_alignments[i].A(kV)=dv;
      fdc_cathode_alignments[i].A(kPhiV)=dphiv;

    }
    fdcfile.close();
  }

  fdc_drift_parms(0)=0.;
  fdc_drift_parms(1)=0.;
  fdc_drift_parms(2)=0.03;
 	
  unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
			      135,135,146,146,158,158,170,170,182,182,197,197,
			      209,209};
  for (unsigned int i=0;i<28;i++){
    vector<cdc_align_t>tempvec;
    for (unsigned int j=0;j<numstraws[i];j++){
      cdc_align_t temp;
      temp.A=DMatrix4x1(0.,0.,0.,0.);
      double var=0.01;
      if (RUN_BENCHMARK==false){
	temp.E=DMatrix4x4(var,0.,0.,0., 0.,var,0.,0., 0.,0.,var,0., 0.,0.,0.,var);
      }
      else {
	temp.E=DMatrix4x4();
      }
      tempvec.push_back(temp);
    }
    cdc_alignments.push_back(tempvec);
  }

  if (FILL_TREE){
    // Create Tree
    fdctree = new TTree("fdc","FDC alignments");
    fdcbranch = fdctree->Branch("T","FDC_branch",&fdc_ptr);
    
    // Create Tree
    fdcCtree = new TTree("fdc_c","FDC alignments");
    fdcCbranch = fdcCtree->Branch("T","FDC_c_branch",&fdc_c_ptr);
    
    // Create Tree
    cdctree = new TTree("cdc","CDC alignments");
    cdcbranch = cdctree->Branch("T","CDC_branch",&cdc_ptr);
  }
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_dc_alignment::brun(JEventLoop *loop, int runnumber)
{	
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  dgeom  = dapp->GetDGeometry(runnumber);
  //dgeom->GetFDCWires(fdcwires);

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_dz,endplate_rmin,endplate_rmax;
  dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z+=0.5*endplate_dz;

   
  JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());
  vector< map<string, double> > tvals;
  cdc_drift_table.clear();
  if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, double> &row = tvals[i];
      cdc_drift_table.push_back(1000.*row["t"]);
    }
  }
  else{
    jerr << " CDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }

  map<string, double> cdc_res_parms;
  jcalib->Get("CDC/cdc_resolution_parms", cdc_res_parms);
  CDC_RES_PAR1 = cdc_res_parms["res_par1"];
  CDC_RES_PAR2 = cdc_res_parms["res_par2"];

  fdc_drift_table.clear();
  if (jcalib->Get("FDC/fdc_drift_table", tvals)==false){    
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, double> &row = tvals[i];
      fdc_drift_table.push_back(1000.*row["t"]);
    }
  }
  else{
    jerr << " FDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }

  dapp->Lock();

  unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
			      135,135,146,146,158,158,170,170,182,182,197,197,
			      209,209};
  for (int i=0;i<28;i++){
    char title[40];
    sprintf(title,"cdc_residual_ring%d",i+1);
    Hcdc_ring_res[i]=(TH2F*)gROOT->FindObject(title);
    if (!Hcdc_ring_res[i]){
      Hcdc_ring_res[i]=new TH2F(title,title,numstraws[i],0.5,numstraws[i]+0.5,
				100,-1,1);
    }
  }
   for (int i=0;i<28;i++){
    char title[40];
    sprintf(title,"cdc_drift_time_ring%d",i+1);
    Hcdc_ring_time[i]=(TH2F*)gROOT->FindObject(title);
    if (!Hcdc_ring_time[i]){
      Hcdc_ring_time[i]=new TH2F(title,title,numstraws[i],0.5,numstraws[i]+0.5,
				 900,-100,800);
    }
  }
    
  Hprob = (TH1F*)gROOT->FindObject("Hprob");
  if (!Hprob){
    Hprob=new TH1F("Hprob","Confidence level for final fit",100,0.0,1.); 
  } 
  Hpseudo_prob = (TH1F*)gROOT->FindObject("Hpseudo_prob");
  if (!Hpseudo_prob){
    Hpseudo_prob=new TH1F("Hpseudo_prob","Confidence level for final fit",100,0.0,1.); 
  } 

  Hcdc_prob = (TH1F*)gROOT->FindObject("Hcdc_prob");
  if (!Hcdc_prob){
    Hcdc_prob=new TH1F("Hcdc_prob","Confidence level for time-based fit",100,0.0,1.); 
  } 
  Hcdc_prelimprob = (TH1F*)gROOT->FindObject("Hcdc_prelimprob");
  if (!Hcdc_prelimprob){
    Hcdc_prelimprob=new TH1F("Hcdc_prelimprob","Confidence level for prelimary fit",100,0.0,1.); 
  } 
  Hintersection_match = (TH1F*)gROOT->FindObject("Hintersection_match");
  if (!Hintersection_match){
    Hintersection_match=new TH1F("Hintersection_match","Intersection matching distance",100,0.0,25.); 
  }
  Hintersection_link_match = (TH1F*)gROOT->FindObject("Hintersection_link_match");
  if (!Hintersection_link_match){
    Hintersection_link_match=new TH1F("Hintersection_link_match","Segment matching distance",100,0.0,25.); 
  }

  Hcdcmatch = (TH1F*)gROOT->FindObject("Hcdcmatch");
  if (!Hcdcmatch){
    Hcdcmatch=new TH1F("Hcdcmatch","CDC hit matching distance",1000,0.0,50.); 
  }
  Hcdcmatch_stereo = (TH1F*)gROOT->FindObject("Hcdcmatch_stereo");
  if (!Hcdcmatch_stereo){
    Hcdcmatch_stereo=new TH1F("Hcdcmatch_stereo","CDC stereo hit matching distance",1000,0.0,50.); 
  }
  
  Hmatch = (TH1F*)gROOT->FindObject("Hmatch");
  if (!Hmatch){
    Hmatch=new TH1F("Hmatch","Segment matching distance",100,0.0,25.); 
  }
  Hlink_match = (TH1F*)gROOT->FindObject("Hlink_match");
  if (!Hlink_match){
    Hlink_match=new TH1F("link_match","Segment matching distance",100,0.0,25.); 
  }
  

  Hbeta = (TH1F*)gROOT->FindObject("Hbeta");
  if (!Hbeta){
    Hbeta=new TH1F("Hbeta","Estimate for #beta",100,0.0,1.5); 
    Hbeta->SetXTitle("#beta");
  }  

  Hztarg = (TH1F*)gROOT->FindObject("Hztarg");
  if (!Hztarg){
    Hztarg=new TH1F("Hztarg","Estimate for target z",1200,-300.0,300.0); 
  }
  Hures_vs_layer=(TH2F*)gROOT->FindObject("Hures_vs_layer");
  if (!Hures_vs_layer){
    Hures_vs_layer=new TH2F("Hures_vs_layer","Cathode u-view residuals",
			    24,0.5,24.5,200,-0.5,0.5);
  }  
  Hres_vs_layer=(TH2F*)gROOT->FindObject("Hres_vs_layer");
  if (!Hres_vs_layer){
    Hres_vs_layer=new TH2F("Hres_vs_layer","wire-based residuals",
			    24,0.5,24.5,200,-0.5,0.5);
  }  
  Hvres_vs_layer=(TH2F*)gROOT->FindObject("Hvres_vs_layer");
  if (!Hvres_vs_layer){
    Hvres_vs_layer=new TH2F("Hvres_vs_layer","Cathode v-view residuals",
			    24,0.5,24.5,200,-0.5,0.5);
  }  
  Hcdc_time_vs_d=(TH2F*)gROOT->FindObject("Hcdc_time_vs_d");
  if (!Hcdc_time_vs_d){
    Hcdc_time_vs_d=new TH2F("Hcdc_time_vs_d",
			    "cdc drift time vs doca",80,0,0.8,400,-20,780);
  } 
  Hcdcdrift_time=(TH2F*)gROOT->FindObject("Hcdcdrift_time");
  if (!Hcdcdrift_time){
    Hcdcdrift_time=new TH2F("Hcdcdrift_time",
			 "cdc doca vs drift time",801,-21,781,100,0,1);
  }  
  Hcdcres_vs_drift_time=(TH2F*)gROOT->FindObject("Hcdcres_vs_drift_time");
  if (!Hcdcres_vs_drift_time){
    Hcdcres_vs_drift_time=new TH2F("Hcdcres_vs_drift_time","cdc Residual vs drift time",400,-20,780,100,-1.,1.);
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
    Hbcalmatch=new TH2F("Hbcalmatch","BCAL #deltar vs #deltaz",100,-50.,50.,
			100,0.,10.);
  }
  Hbcalmatchxy=(TH2F*)gROOT->FindObject("Hbcalmatchxy");
  if (!Hbcalmatchxy){
    Hbcalmatchxy=new TH2F("Hbcalmatchxy","BCAL #deltay vs #deltax",400,-50.,50.,
			400,-50.,50.);
  }
  Hfcalmatch=(TH1F*)gROOT->FindObject("Hfcalmatch");
  if (!Hfcalmatch){
    Hfcalmatch=new TH1F("Hfcalmatch","FCAL #deltar",400,0.,50.);
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

  printf("Events processed = %d\n",myevt);

  if (RUN_BENCHMARK==false){
    ofstream cdcfile("cdc_alignment.dat");
    cdcfile << "Ring straw dXu dYu dXd dYd" << endl;
    for (unsigned int ring=0;ring<cdc_alignments.size();ring++){
      for (unsigned int straw=0;straw<cdc_alignments[ring].size();
	   straw++){
	cdcfile << ring+1 << " " << straw+1 << " " 
		<< cdc_alignments[ring][straw].A(k_dXu) << " " 
		<< cdc_alignments[ring][straw].A(k_dYu) << " "
		<< cdc_alignments[ring][straw].A(k_dXd) << " "
	      << cdc_alignments[ring][straw].A(k_dYd) << endl;
      }
    }
    cdcfile.close();
    
    if (ALIGN_WIRE_PLANES){
      ofstream fdcfile("fdc_alignment.dat");
      //fdcfile << "dPhi dU sig(dU)" <<endl;
      for (unsigned int layer=0;layer<24;layer++){
	double du=fdc_alignments[layer].A(kU);
	double dphi=fdc_alignments[layer].A(kPhiU);
	
      fdcfile <<  dphi <<" " <<" " << du << " " << "0." << endl;
      }
      fdcfile.close(); 
    }
    else{
      ofstream fdcfile("fdc_cathode_alignment.dat");
      for (unsigned int layer=0;layer<24;layer++){
	fdcfile << fdc_cathode_alignments[layer].A(kPhiU)
		<< " " << fdc_cathode_alignments[layer].A(kU)
		<< " " << fdc_cathode_alignments[layer].A(kPhiV)
		<< " " << fdc_cathode_alignments[layer].A(kV) << endl;
      }
      fdcfile.close(); 
    }
  }
  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_dc_alignment::evnt(JEventLoop *loop, int eventnumber){
  myevt++;
  
  // Get BCAL showers, FCAL showers and FDC space points
  vector<const DFCALShower*>fcalshowers;
  if (USE_FCAL) loop->Get(fcalshowers);
  vector<const DBCALShower*>bcalshowers;
  if (USE_BCAL)loop->Get(bcalshowers);

  vector<const DFDCPseudo*>pseudos;
  //loop->Get(pseudos);
  vector<const DCDCTrackHit*>cdcs;
  if (COSMICS) loop->Get(cdcs);
  vector<const DFDCIntersection*>intersections;
  if (ALIGN_WIRE_PLANES)loop->Get(intersections);
  else loop->Get(pseudos);

  if (cdcs.size()>20 /* && cdcs.size()<60*/){
    vector<const DCDCTrackHit*>axialhits;
    vector<const DCDCTrackHit*>stereohits;
    for (unsigned int i=0;i<cdcs.size();i++){
      int ring=cdcs[i]->wire->ring;
      if (ring<=4) axialhits.push_back(cdcs[i]);
      else if (ring<=12) stereohits.push_back(cdcs[i]);
      else if (ring<=16) axialhits.push_back(cdcs[i]);
      else if (ring<=24) stereohits.push_back(cdcs[i]);
      else axialhits.push_back(cdcs[i]);
    }

    // Here we group axial hits into segments
    vector<cdc_segment_t>axial_segments;
    vector<bool>used_axial(axialhits.size());
    FindSegments(axialhits,axial_segments,used_axial);

    if (axial_segments.size()>0 && stereohits.size()>0){
      // Here we link axial segments into tracks and associate stereo hits 
      // with each track
      vector<cdc_track_t>tracks;
      LinkSegments(axial_segments,used_axial,axialhits,stereohits,tracks);

      for (unsigned int i=0;i<tracks.size();i++){
	// Add lists of stereo and axial hits associated with this track 
	  // and sort
	vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
	hits.insert(hits.end(),tracks[i].stereo_hits.begin(),tracks[i].stereo_hits.end());
	sort(hits.begin(),hits.end(),cdc_hit_cmp);
	
	DMatrix4x1 S;
	if (USE_BCAL==false){
	  // Use earliest cdc time to estimate t0
	  double minT=1e6;
	  for (unsigned int j=0;j<hits.size();j++){
	    double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
	    double t_test=hits[j]->tdrift-L/29.98;
	    if (t_test<minT) minT=t_test;
	  }
	  mT0=minT;
	  //mT0=0.;

	  // Initial guess for state vector
	  if (GuessForStateVector(tracks[i],S)==NOERROR){
	    // Run the Kalman Filter algorithm
	    DoFilter(S,hits);	 
	  }
	}
	else if (MatchOuterDetectors(tracks[i],fcalshowers,bcalshowers,S)){
	  // Run the Kalman Filter algorithm
	  DoFilter(S,hits);
	  
	} // match outer detectors
      }
    }
  }
  //-------------------------------------------------------------------------
  // FDC alignment 
  //-------------------------------------------------------------------------
  if (ALIGN_WIRE_PLANES){ // Use intersection points to align wire planes
    if (intersections.size()>MIN_INTERSECTIONS
	//&&((fcalshowers.size()>0&&fcalshowers.size()<3)
	//   || (bcalshowers.size()>0&&bcalshowers.size()<3))
	){

      // Group FDC hits by package
      vector<const DFDCIntersection*>packages[4];
      for (unsigned int i=0;i<intersections.size();i++){
	packages[(intersections[i]->wire1->layer-1)/6].push_back(intersections[i]);
      }
      
      // Link hits in each package together into track segments
      vector<intersection_segment_t>segments[4];
      for (unsigned int i=0;i<4;i++){  
	FindSegments(packages[i],segments[i]);
      }
      // Link the segments together to form track candidadates
      vector<intersection_segment_t>LinkedSegments;
      LinkSegments(segments,LinkedSegments);
      
      // Loop over linked segments
      unsigned int num_linked_segments=LinkedSegments.size();
      if (num_linked_segments>0 && num_linked_segments<3){
	for (unsigned int k=0;k<LinkedSegments.size();k++){
	  vector<const DFDCIntersection *>intersections=LinkedSegments[k].hits;
	  
	  if (intersections.size()>MIN_INTERSECTIONS){
	    // Initial guess for state vector    
	    DMatrix4x1 S=LinkedSegments[k].S;
	    
	    // Move x and y to just before the first hit
	    double my_z=intersections[0]->pos.z()-2.;
	    S(state_x)+=my_z*S(state_tx);
	    S(state_y)+=my_z*S(state_ty);
	    
	    if (fabs(S(state_x))>48. || fabs(S(state_y))>48.) continue;
	    
	    DoFilter(my_z,S,intersections);
	  }
	}
      }
    }
  }
  else{  // Use pseudopoints to align cathode planes
    if (pseudos.size()>MIN_PSEUDOS
	//&&((fcalshowers.size()>0&&fcalshowers.size()<3)
	//   || (bcalshowers.size()>0&&bcalshowers.size()<3))
	){
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
      vector<segment_t>LinkedSegments;
      LinkSegments(segments,LinkedSegments);
      
      // Loop over linked segments
      size_t num_linked_segments=LinkedSegments.size();
      if (num_linked_segments>0 && num_linked_segments<3){
	for (unsigned int k=0;k<LinkedSegments.size();k++){
	  vector<const DFDCPseudo *>hits=LinkedSegments[k].hits;
	  
	  if (hits.size()>MIN_PSEUDOS){
	    sort(hits.begin(),hits.end(),fdc_pseudo_cmp);
	    
	    // Initial guess for state vector    
	    DMatrix4x1 S=LinkedSegments[k].S;
	  
	    // Move x and y to just before the first hit
	    double my_z=hits[0]->wire->origin.z()-1.;
	    S(state_x)+=my_z*S(state_tx);
	    S(state_y)+=my_z*S(state_ty);
	    
	    if (fabs(S(state_x))>48. || fabs(S(state_y))>48.) continue;
	    
	    mMinTime=1e6;
	    mOuterTime=1e6;
	    mOuterZ=0.;
	    double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
	    for (unsigned int m=0;m<hits.size();m++){
	      if (hits[m]->time<mMinTime){
		double L=(hits[m]->wire->origin.z()-my_z)*dsdz;
		mMinTime=hits[m]->time-L/29.98; // assume moving at speed of light
		mMinTimeID=m;
	      }
	    }
	    mT0=mMinTime;
	    
	    // Run the Kalman Filter algorithm
	    if (DoFilter(my_z,S,hits)==NOERROR){
	      
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
				       vector<const DCDCTrackHit *>&hits){
  unsigned int numhits=hits.size();
  unsigned int maxindex=numhits-1;

  int NEVENTS=300000;
  double anneal_factor=pow(1e4,(double(NEVENTS-myevt))/(NEVENTS-1.));
  if (myevt>NEVENTS) anneal_factor=1.;  
  anneal_factor=1.;
  if (RUN_BENCHMARK) anneal_factor=1.;

  // deques to store reference trajectories
  deque<trajectory_t>trajectory;
  deque<trajectory_t>best_traj;
  // if (SetReferenceTrajectory(mOuterZ,S,trajectory,hits[maxindex])==NOERROR)
    {
    // State vector to store "best" values
    DMatrix4x1 Sbest;

    // Covariance matrix
    DMatrix4x4 C0,C,Cbest;
    C0(state_x,state_x)=C0(state_y,state_y)=1.0;     
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

      if (fabs(chi2_old-chi2)<0.1 || chi2>chi2_old) break;  
      
      // Save the current state and covariance matrixes
      Cbest=C;
      Sbest=S; 
      best_updates.assign(updates.begin(),updates.end());
      best_traj.assign(trajectory.begin(),trajectory.end());

      // run the smoother (opposite direction to filter)
      //Smooth(S,C,trajectory,updates);	      
    }
    if (iter>0){
      double prelimprob=TMath::Prob(chi2_old,ndof_old);
      Hcdc_prelimprob->Fill(prelimprob);

      if (prelimprob>0.001){
 
	// Perform a time-based pass
	S=Sbest;
	chi2=1e16;
	/*
        printf("theta_x %f theta_y %f\n",180./M_PI*atan(S(state_tx)),
	       180./M_PI*atan(S(state_ty)));
	*/
	//printf("xyz %f %f %f\n",S(state_x),S(state_y),best_traj[0].z);

	// Match to outer detectors
	/*
	if (MatchOuterDetectors(fcalshowers,bcalshowers,S)){
	  // move S to the z-position where we found the match
	  S(state_x)+=mOuterZ*S(state_tx);
	  S(state_y)+=mOuterZ*S(state_ty);
	*/

	//printf("Timebased-----------\n");
	//if (false)
	for (iter=0;iter<20;iter++){
	  chi2_old=chi2; 
	  ndof_old=ndof;
	  
	  trajectory.clear();
	  if (SetReferenceTrajectory(mOuterZ,S,trajectory,hits[maxindex])
	      ==NOERROR){
	    C=C0;
	    KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof,true);
  
	    //printf(">>>>>>chi2 %f ndof %d\n",chi2,ndof);
	    if (fabs(chi2-chi2_old)<0.1  
		|| TMath::Prob(chi2,ndof)<TMath::Prob(chi2_old,ndof_old)) break;
	    
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

	  PlotLines(trajectory);

	  if (prob>1e-6)
	    {
	    // run the smoother (opposite direction to filter)
	    vector<cdc_update_t>smoothed_updates(updates.size());
	    for (unsigned int k=0;k<smoothed_updates.size();k++){
	      smoothed_updates[k].used_in_fit=false;
	    }
	    Smooth(Sbest,Cbest,best_traj,hits,best_updates,smoothed_updates);

	    for (unsigned int k=0;k<smoothed_updates.size();k++){
	      if (smoothed_updates[k].used_in_fit==true){
		double tdrift=smoothed_updates[k].drift_time;
		double d=smoothed_updates[k].doca;
		double res=smoothed_updates[k].res;
		int ring_id=smoothed_updates[k].ring_id;
		int straw_id=smoothed_updates[k].straw_id;
		Hcdcres_vs_drift_time->Fill(tdrift,res);
		Hcdcdrift_time->Fill(tdrift,d);
		Hcdc_time_vs_d->Fill(d,tdrift);
	        Hcdc_ring_res[ring_id]->Fill(straw_id+1,res); 
		Hcdc_ring_time[ring_id]->Fill(straw_id+1,tdrift);
	      }
	    }
	    
	    if (prob>0.001 && RUN_BENCHMARK==false){
	      FindOffsets(hits,smoothed_updates);
	    	
	      if (FILL_TREE){
		for (unsigned int ring=0;ring<cdc_alignments.size();ring++){
		  for (unsigned int straw=0;straw<cdc_alignments[ring].size();
		       straw++){
		    // Set up to fill tree
		    cdc.dXu=cdc_alignments[ring][straw].A(k_dXu);
		    cdc.dYu=cdc_alignments[ring][straw].A(k_dYu);
		    cdc.dXd=cdc_alignments[ring][straw].A(k_dXd);
		    cdc.dYd=cdc_alignments[ring][straw].A(k_dYd);
		    cdc.straw=straw+1;
		    cdc.ring=ring+1;
		    cdc.N=myevt;
		  
		
		    // Lock mutex
		    pthread_mutex_lock(&mutex);
		    
		    cdctree->Fill();
		    
		    // Unlock mutex
		    pthread_mutex_unlock(&mutex);
		    
		  }
		}
	      }
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
DEventProcessor_dc_alignment::DoFilter(double start_z,DMatrix4x1 &S,
				       vector<const DFDCPseudo *>&hits){
  unsigned int num_hits=hits.size();
  vector<update_t>updates(num_hits);
  vector<update_t>best_updates;
  vector<update_t>smoothed_updates(num_hits);  

  int NEVENTS=100000;
  double anneal_factor=pow(1e3,(double(NEVENTS-myevt))/(NEVENTS-1.));
  if (myevt>NEVENTS) anneal_factor=1.;  
  //anneal_factor=10.;
  if (RUN_BENCHMARK) anneal_factor=1.;
  //anneal_factor=1e3;

  // Best guess for state vector at the beginning of the trajectory
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;
  deque<trajectory_t>best_traj;
      
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
    if (SetReferenceTrajectory(start_z,S,trajectory,hits)!=NOERROR) break;
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
    Hpseudo_prob->Fill(prob);
    
    // printf("prob %f\n",prob);

    PlotLines(trajectory);

    if (prob>0.00001)
      {
      // run the smoother (opposite direction to filter)
      Smooth(Sbest,Cbest,best_traj,hits,best_updates,smoothed_updates);

      //Hbeta->Fill(mBeta);
      for (unsigned int i=0;i<smoothed_updates.size();i++){
	unsigned int layer=hits[i]->wire->layer;
     
	Hures_vs_layer->Fill(layer,smoothed_updates[i].res(0));
	Hvres_vs_layer->Fill(layer,smoothed_updates[i].res(1));
	Hdv_vs_dE->Fill(hits[i]->dE,smoothed_updates[i].res(1));

	Hdrift_time->Fill(smoothed_updates[i].drift_time,
			  smoothed_updates[i].doca);
      }

      if (prob>0.001 && RUN_BENCHMARK==false){
	FindOffsets(hits,smoothed_updates);
	
	if (FILL_TREE){
	  for (unsigned int layer=0;layer<24;layer++){
	    fdc_c.dPhiU=fdc_cathode_alignments[layer].A(kPhiU);
	    fdc_c.dU=fdc_cathode_alignments[layer].A(kU);
	    fdc_c.dPhiV=fdc_cathode_alignments[layer].A(kPhiV);
	    fdc_c.dV=fdc_cathode_alignments[layer].A(kV);
	    
	    fdc_c.layer=layer+1;
	    fdc_c.N=myevt;
	    
	    // Lock mutex
	    pthread_mutex_lock(&mutex);
	    
	    fdcCtree->Fill();
	  
	    // Unlock mutex
	    pthread_mutex_unlock(&mutex);
	  }
	}
      }
      return NOERROR;
    }
  }
   
  
  return VALUE_OUT_OF_RANGE;
}

// Steering routine for the kalman filter
jerror_t 
DEventProcessor_dc_alignment::DoFilter(double start_z,DMatrix4x1 &S,
				       vector<const DFDCIntersection *>&intersections){
  vector<intersection_hit_t>hits;
  unsigned int max_i=intersections.size()-1; 
  intersection_hit_t temp;
  for (unsigned int i=0;i<max_i;i++){    
    temp.wire=intersections[i]->wire1;
    temp.hit=intersections[i]->hit1;
    hits.push_back(temp);
    if (intersections[i]->wire2->layer!=intersections[i+1]->wire1->layer){
      temp.wire=intersections[i]->wire2;
      temp.hit=intersections[i]->hit2;
      hits.push_back(temp);
    }
  }
  temp.wire=intersections[max_i]->wire1;
  temp.hit=intersections[max_i]->hit1;
  hits.push_back(temp);
  temp.wire=intersections[max_i]->wire2;
  temp.hit=intersections[max_i]->hit2;
  hits.push_back(temp);

  unsigned int num_hits=hits.size();
  if (num_hits<10) return VALUE_OUT_OF_RANGE;

  mMinTime=1e6;
  mOuterTime=1e6;
  mOuterZ=0.;
  double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  for (unsigned int m=0;m<hits.size();m++){
    if (hits[m].hit->t<mMinTime){
      double L=(hits[m].wire->origin.z()-start_z)*dsdz;
      mMinTime=hits[m].hit->t-L/29.98; // assume moving at speed of light
      mMinTimeID=m;
    }
  }
  mT0=mMinTime;

  vector<wire_update_t>updates(num_hits);
  vector<wire_update_t>best_updates;
  vector<wire_update_t>smoothed_updates(num_hits);  

  int NEVENTS=75000;
  double anneal_factor=1.;
  if (USE_DRIFT_TIMES){
    anneal_factor=pow(1000.,(double(NEVENTS-myevt))/(NEVENTS-1.));
    if (myevt>NEVENTS) anneal_factor=1.;  
  }
  if (RUN_BENCHMARK) anneal_factor=1.;
  //anneal_factor=1e3;

  // Best guess for state vector at "vertex"
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;
  deque<trajectory_t>best_traj;
 
  // Intial guess for covariance matrix
  DMatrix4x4 C,C0,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=1.;
  C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.001;
  
  // Chi-squared and degrees of freedom
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0,ndof_old=0;
  unsigned iter=0;
  for(;;){
    iter++;
    chi2_old=chi2; 
    ndof_old=ndof;

    trajectory.clear();
    if (SetReferenceTrajectory(start_z,S,trajectory,hits)!=NOERROR) break;
    C=C0;
    if (KalmanFilter(anneal_factor,S,C,hits,trajectory,updates,chi2,ndof)
	!=NOERROR) break;

    //printf("== event %d == iter %d =====chi2 %f ndof %d \n",myevt,iter,chi2,ndof);
    if (chi2>chi2_old || iter==ITER_MAX) break;  
    
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
    
    PlotLines(trajectory);

    if (prob>0.001)
      {
      // run the smoother (opposite direction to filter)
      Smooth(Sbest,Cbest,best_traj,hits,best_updates,smoothed_updates);   

      //Hbeta->Fill(mBeta);
      for (unsigned int i=0;i<smoothed_updates.size();i++){
	unsigned int layer=hits[i].wire->layer;
     
	Hres_vs_layer->Fill(layer,smoothed_updates[i].ures);
	if (prob>0.1/*&&layer==smoothed_updates.size()/2*/){
	  Hdrift_time->Fill(smoothed_updates[i].drift_time,
			    smoothed_updates[i].doca);
	  Hres_vs_drift_time->Fill(smoothed_updates[i].drift_time,
				   smoothed_updates[i].ures);
	 
	}
	
      }
      
      if (RUN_BENCHMARK==false){
	FindOffsets(hits,smoothed_updates);
	
	if (FILL_TREE){
	  for (unsigned int layer=0;layer<24;layer++){
	    fdc.dPhi=fdc_alignments[layer].A(kPhiU);
	    fdc.dX=fdc_alignments[layer].A(kU);
	    
	    fdc.layer=layer+1;
	    fdc.N=myevt;
	    
	    // Lock mutex
	    pthread_mutex_lock(&mutex);
	    
	    fdctree->Fill();
	    
	    // Unlock mutex
	    pthread_mutex_unlock(&mutex);
	  }
	}
      }
      return NOERROR;

    }
  }
   
  
  return VALUE_OUT_OF_RANGE;
}




// Link segments from package to package by doing straight-line projections
jerror_t
DEventProcessor_dc_alignment::LinkSegments(vector<segment_t>segments[4], 
					   vector<segment_t>&LinkedSegments){
  vector<const DFDCPseudo *>myhits;
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<segments[i].size();j++){
      if (segments[i][j].matched==false){
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
	      
	      Hlink_match->Fill((proj-segments[i_plus_1][k].hits[0]->xy).Mod());

	      if ((proj-segments[i_plus_1][k].hits[0]->xy).Mod()<PSEUDO_LINK_MATCH_RADIUS){
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
		      
		      if ((proj-segments[i_plus_2][m].hits[0]->xy).Mod()<PSEUDO_LINK_MATCH_RADIUS){
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
			      
			      if ((proj-segments[i_plus_3][n].hits[0]->xy).Mod()<PSEUDO_LINK_MATCH_RADIUS){
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
	if (myhits.size()>0){ 
	  segments[i][j].matched=true;
	  myhits.insert(myhits.begin(),segments[i][j].hits.begin(),segments[i][j].hits.end());	
	  segment_t mysegments;
	  mysegments.hits.assign(myhits.begin(),myhits.end());
	  mysegments.S=FitLine(myhits);
	  LinkedSegments.push_back(mysegments);
	}	  
	myhits.clear();
      } // check if we have already used this segment
    } // loop over first set of segments
  } // loop over packages

  // Try to link tracks together
  if (LinkedSegments.size()>1){
    for (unsigned int i=0;i<LinkedSegments.size()-1;i++){
      DMatrix4x1 S=LinkedSegments[i].S;
      size_t last_index_1=LinkedSegments[i].hits.size()-1;
      int first_pack_1=(LinkedSegments[i].hits[0]->wire->layer-1)/6;
      int last_pack_1=(LinkedSegments[i].hits[last_index_1]->wire->layer-1)/6;
      for (unsigned int j=i+1;j<LinkedSegments.size();j++){
	size_t last_index_2=LinkedSegments[j].hits.size()-1;
	int first_pack_2=(LinkedSegments[j].hits[0]->wire->layer-1)/6;
	int last_pack_2=(LinkedSegments[j].hits[last_index_2]->wire->layer-1)/6;

	if (last_pack_1<first_pack_2 || first_pack_1 > last_pack_2){
	  double z=LinkedSegments[j].hits[0]->wire->origin.z();
	  DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	  double diff=(LinkedSegments[j].hits[0]->xy-proj).Mod();

	  if (diff<MATCH_RADIUS){
	    // Combine the hits from the two tracks and recompute the state 
	    // vector S
	    if (last_pack_1<first_pack_2){
	      LinkedSegments[i].hits.insert(LinkedSegments[i].hits.end(),
					    LinkedSegments[j].hits.begin(),
					    LinkedSegments[j].hits.end());
	    }
	    else{
	      LinkedSegments[i].hits.insert(LinkedSegments[i].hits.begin(),
					    LinkedSegments[j].hits.begin(),
					    LinkedSegments[j].hits.end());
	    }
	    LinkedSegments[i].S=FitLine(LinkedSegments[i].hits);

	    // Drop the second track from the list 
	    LinkedSegments.erase(LinkedSegments.begin()+j);
	    break;
	  }
	}
      }
    }
  }

  // Try to attach unmatched segments to tracks
  for (unsigned int i=0;i<LinkedSegments.size();i++){
    DMatrix4x1 S=LinkedSegments[i].S;
    size_t last_index_1=LinkedSegments[i].hits.size()-1;
    int first_pack_1=(LinkedSegments[i].hits[0]->wire->layer-1)/6;
    int last_pack_1=(LinkedSegments[i].hits[last_index_1]->wire->layer-1)/6; 
    for (unsigned int j=0;j<4;j++){
      for (unsigned int k=0;k<segments[j].size();k++){
	if (segments[j][k].matched==false){
	  int pack_2=(segments[j][k].hits[0]->wire->layer-1)/6;
	  if (pack_2<first_pack_1 || pack_2>last_pack_1){
	    double z=segments[j][k].hits[0]->wire->origin.z();
	    DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	    double diff=(segments[j][k].hits[0]->xy-proj).Mod();

	    if (diff<MATCH_RADIUS){
	      segments[j][k].matched=true;
	      
	      // Add hits and recompute S
	      if (pack_2<first_pack_1){
		LinkedSegments[i].hits.insert(LinkedSegments[i].hits.begin(),
					      segments[j][k].hits.begin(),
					      segments[j][k].hits.end());
	      }
	      else {
		LinkedSegments[i].hits.insert(LinkedSegments[i].hits.end(),
					      segments[j][k].hits.begin(),
					      segments[j][k].hits.end());
	      
	      }
	      LinkedSegments[i].S=FitLine(LinkedSegments[i].hits);
	    }
	  }
	}
      }
    }
  }


  
  return NOERROR;
}



// Link segments from package to package by doing straight-line projections
jerror_t
DEventProcessor_dc_alignment::LinkSegments(vector<intersection_segment_t>segments[4], 
					   vector<intersection_segment_t>&LinkedSegments){
  vector<const DFDCIntersection *>myhits;
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<segments[i].size();j++){
      if (segments[i][j].matched==false){
	unsigned i_plus_1=i+1; 
	if (i_plus_1<4){
	  double tx=segments[i][j].S(state_tx);
	  double ty=segments[i][j].S(state_ty);
	  double x0=segments[i][j].S(state_x);
	  double y0=segments[i][j].S(state_y);
	  
	  for (unsigned int k=0;k<segments[i_plus_1].size();k++){
	    if (segments[i_plus_1][k].matched==false){
	      double z=segments[i_plus_1][k].hits[0]->pos.z();
	      DVector2 proj(x0+tx*z,y0+ty*z);
	      DVector2 XY(segments[i_plus_1][k].hits[0]->pos.x(),
			  segments[i_plus_1][k].hits[0]->pos.y());
	     
	      Hintersection_link_match->Fill((proj-XY).Mod());

	      if ((proj-XY).Mod()<INTERSECTION_LINK_MATCH_RADIUS){
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
		      z=segments[i_plus_2][m].hits[0]->pos.z();
		      proj.Set(x0+tx*z,y0+ty*z);
		      XY.Set(segments[i_plus_2][m].hits[0]->pos.x(),
			     segments[i_plus_2][m].hits[0]->pos.y());
		      if ((proj-XY).Mod()<INTERSECTION_LINK_MATCH_RADIUS){
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
			      z=segments[i_plus_3][n].hits[0]->pos.z();
			      proj.Set(x0+tx*z,y0+ty*z);
			      XY.Set(segments[i_plus_3][n].hits[0]->pos.x(),
				     segments[i_plus_3][n].hits[0]->pos.y());
			      if ((proj-XY).Mod()<INTERSECTION_LINK_MATCH_RADIUS){
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
	if (myhits.size()>0){
	  segments[i][j].matched=true;
	  myhits.insert(myhits.begin(),segments[i][j].hits.begin(),segments[i][j].hits.end());	
	  intersection_segment_t mysegments;
	  mysegments.hits.assign(myhits.begin(),myhits.end());
	  mysegments.S=FitLine(myhits);
	  LinkedSegments.push_back(mysegments);
	}
	myhits.clear();
      } // check if we have already used this segment
    } // loop over first set of segments
  } // loop over packages

  
  // Try to link tracks together
  if (LinkedSegments.size()>1){
    for (unsigned int i=0;i<LinkedSegments.size()-1;i++){
      DMatrix4x1 S=LinkedSegments[i].S;
      size_t last_index_1=LinkedSegments[i].hits.size()-1;
      int first_pack_1=LinkedSegments[i].hits[0]->wire1->layer/6;
      int last_pack_1=LinkedSegments[i].hits[last_index_1]->wire1->layer/6;
      for (unsigned int j=i+1;j<LinkedSegments.size();j++){
	size_t last_index_2=LinkedSegments[j].hits.size()-1;
	int first_pack_2=LinkedSegments[j].hits[0]->wire1->layer/6;
	int last_pack_2=LinkedSegments[j].hits[last_index_2]->wire1->layer/6;
	if (last_pack_1<first_pack_2 || first_pack_1 > last_pack_2){
	  double z=LinkedSegments[j].hits[0]->pos.z();
	  DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	  DVector2 XY(LinkedSegments[j].hits[0]->pos.x(),LinkedSegments[j].hits[0]->pos.y());
	  double diff=(XY-proj).Mod();
	  if (diff<MATCH_RADIUS){
	    // Combine the hits from the two tracks and recompute the state 
	    // vector S
	    if (last_pack_1<first_pack_2){
	      LinkedSegments[i].hits.insert(LinkedSegments[i].hits.end(),
					    LinkedSegments[j].hits.begin(),
					    LinkedSegments[j].hits.end());
	    }
	    else{
	      LinkedSegments[i].hits.insert(LinkedSegments[i].hits.begin(),
					    LinkedSegments[j].hits.begin(),
					    LinkedSegments[j].hits.end());
	    }
	    LinkedSegments[i].S=FitLine(LinkedSegments[i].hits);

	    // Drop the second track from the list 
	    LinkedSegments.erase(LinkedSegments.begin()+j);
	    break;
	  }
	}
      }
    }
  }

  // Try to attach unmatched segments to tracks
  for (unsigned int i=0;i<LinkedSegments.size();i++){
    DMatrix4x1 S=LinkedSegments[i].S;
    size_t last_index_1=LinkedSegments[i].hits.size()-1;
    int first_pack_1=LinkedSegments[i].hits[0]->wire1->layer/6;
    int last_pack_1=LinkedSegments[i].hits[last_index_1]->wire1->layer/6; 
    for (unsigned int j=0;j<4;j++){
      for (unsigned int k=0;k<segments[j].size();k++){
	if (segments[j][k].matched==false){
	  int pack_2=segments[j][k].hits[0]->wire1->layer/6;
	  if (pack_2<first_pack_1 || pack_2>last_pack_1){
	    double z=segments[j][k].hits[0]->pos.z();
	    DVector2 proj(S(state_x)+z*S(state_tx),S(state_y)+z*S(state_ty));
	    DVector2 XY(segments[j][k].hits[0]->pos.x(),
		      segments[j][k].hits[0]->pos.y());
	    double diff=(XY-proj).Mod();

	    if (diff<MATCH_RADIUS){
	      segments[j][k].matched=true;
	      
	      // Add hits and recompute S
	      if (pack_2<first_pack_1){
		LinkedSegments[i].hits.insert(LinkedSegments[i].hits.begin(),
					      segments[j][k].hits.begin(),
					    segments[j][k].hits.end());
	      }
	      else {
		LinkedSegments[i].hits.insert(LinkedSegments[i].hits.end(),
					      segments[j][k].hits.begin(),
					    segments[j][k].hits.end());
	      
	      }
	      LinkedSegments[i].S=FitLine(LinkedSegments[i].hits);
	    }
	  }
	}
      }
    }
  }
  
  return NOERROR;
}


// Find segments in cdc axial layers
jerror_t 
DEventProcessor_dc_alignment::FindSegments(vector<const DCDCTrackHit*>&hits,
					   vector<cdc_segment_t>&segments,
					   vector<bool>&used_in_segment
					   ){
  if (hits.size()==0) return RESOURCE_UNAVAILABLE;

  // Group adjacent hits into pairs
  vector<pair<unsigned int,unsigned int> > pairs;
  for (unsigned int i=0;i<hits.size()-1;i++){
    for (unsigned int j=i+1;j<hits.size();j++){
      int r1=hits[i]->wire->ring;
      int r2=hits[j]->wire->ring;
      int s1=hits[i]->wire->straw;
      int s2=hits[j]->wire->straw;
      double d=(hits[i]->wire->origin-hits[j]->wire->origin).Perp();
      
      {
	pthread_mutex_lock(&mutex);
	Hcdcmatch->Fill(d);
	pthread_mutex_unlock(&mutex);
      }

      if ((abs(r1-r2)<=2 && d<CDC_MATCH_RADIUS) 
	  || (abs(r1-r2)==0 && abs(s1-s2)==1)){
	pair <unsigned int,unsigned int> mypair(i,j);
	pairs.push_back(mypair);
      }
    }
  }
  // Link pairs of hits together into segments
  for (unsigned int i=0;i<pairs.size();i++){
    if (used_in_segment[pairs[i].first]==false 
	&& used_in_segment[pairs[i].second]==false){
      vector<const DCDCTrackHit *>neighbors;
      unsigned int old=i;
      unsigned int old_first=pairs[old].first;
      unsigned int old_second=pairs[old].second;
      used_in_segment[old_first]=true;
      used_in_segment[old_second]=true;
      neighbors.push_back(hits[old_first]);
      neighbors.push_back(hits[old_second]);
      for (unsigned int j=i+1;j<pairs.size();j++){
	unsigned int first=pairs[j].first;
	unsigned int second=pairs[j].second;
	old_first=pairs[old].first;
	old_second=pairs[old].second;
	if ((used_in_segment[old_first] || used_in_segment[old_second])
	    && (first==old_first || first==old_second || second==old_second
		|| second==old_first)){
	  if (used_in_segment[first]==false){
	    used_in_segment[first]=true;
	    neighbors.push_back(hits[first]);
	  }
	  if (used_in_segment[second]==false){
	    used_in_segment[second]=true;
	    neighbors.push_back(hits[second]);
	  }  
	  if (used_in_segment[old_first]==false){
	    used_in_segment[old_first]=true;
	    neighbors.push_back(hits[old_first]);
	  }
	  if (used_in_segment[old_second]==false){
	    used_in_segment[old_second]=true;
	    neighbors.push_back(hits[old_second]);
	  }
	}
	old=j;
      }

      cdc_segment_t mysegment; 
      sort(neighbors.begin(),neighbors.end(),cdc_hit_cmp);
      mysegment.dir=neighbors[neighbors.size()-1]->wire->origin
	-neighbors[0]->wire->origin;
      mysegment.dir.SetMag(1.);
      mysegment.hits=neighbors;
      mysegment.matched=false;
      segments.push_back(mysegment);
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
	    if (delta<delta_min){
	      delta_min=delta;
	      if (delta<MATCH_RADIUS) match=m;
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

	if (neighbors.size()>3){
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

// Find segments by associating adjacent hits within a package together.
jerror_t DEventProcessor_dc_alignment::FindSegments(vector<const DFDCIntersection*>&points,
					vector<intersection_segment_t>&segments){
  if (points.size()==0) return RESOURCE_UNAVAILABLE;
  vector<int>used(points.size());

  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=points[0]->pos.z();
  vector<unsigned int>x_list;
  x_list.push_back(0);
  for (unsigned int i=0;i<points.size();i++){
    used.push_back(false);
    if (points[i]->pos.z()!=old_z){
      x_list.push_back(i);
    }
    old_z=points[i]->pos.z();
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
	DVector2 XY(points[i]->pos.x(),points[i]->pos.y());
	
	// Create list of nearest neighbors
	vector<const DFDCIntersection*>neighbors;
	neighbors.push_back(points[i]);
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<x_list.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	    DVector2 XY2(points[m]->pos.x(),points[m]->pos.y());
	    delta=(XY-XY2).Mod();
	    if (delta<delta_min){
	      delta_min=delta;
	      if (delta<MATCH_RADIUS) match=m;
	    }
	  }
	  Hintersection_match->Fill(delta_min);
	  if (//match!=0 
	      //&& 
	      used[match]==false
	      ){
	    XY.Set(points[match]->pos.x(),points[match]->pos.y());
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
	      DVector2 XY(points[k]->pos.x(),points[k]->pos.y());
	      DVector2 XY2(neighbors[j]->pos.x(),neighbors[j]->pos.y());
	      delta=(XY-XY2).Mod();

	      if (delta<ADJACENT_MATCH_RADIUS && 
		  neighbors[j]->pos.z()==points[k]->pos.z()){
		used[k]=true;
		neighbors.push_back(points[k]);
		do_sort=true;
	      }      
	    }
	  }
	} // loop looking for hits adjacent to hits on segment

	if (neighbors.size()>3){
	  intersection_segment_t mysegment;
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

  double sig2v=0.25; // rough guess;

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

// Use linear regression on the hits to obtain a first guess for the state
// vector.  Method taken from Numerical Recipes in C.
DMatrix4x1 
DEventProcessor_dc_alignment::FitLine(vector<const DFDCIntersection*> &fdchits){
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

  for (unsigned int i=0;i<fdchits.size();i++){
    double x=fdchits[i]->pos.X();
    double y=fdchits[i]->pos.Y();
    double z=fdchits[i]->pos.Z();

    S1+=1.0;  // assume all errors are the same
    S1z+=z;
    S1y+=y;
    S1zz+=z*z;
    S1zy+=z*y;    
    
    S2+=1.0;
    S2z+=z;
    S2x+=x;
    S2zz+=z*z;
    S2zx+=z*x;
  }
  double D1=S1*S1zz-S1z*S1z;
  double y_intercept=(S1zz*S1y-S1z*S1zy)/D1;
  double y_slope=(S1*S1zy-S1z*S1y)/D1;
  double D2=S2*S2zz-S2z*S2z;
  double x_intercept=(S2zz*S2x-S2z*S2zx)/D2;
  double x_slope=(S2*S2zx-S2z*S2x)/D2;

  return DMatrix4x1(x_intercept,y_intercept,x_slope,y_slope);

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
  DMatrix2x1 Mdiff;

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
	
	// Get the aligment parameters for this layer
	unsigned int layer=hits[id]->wire->layer-1;
	DMatrix4x1 A=fdc_cathode_alignments[layer].A;
	DMatrix2x1 Aw=fdc_alignments[layer].A;
	double delta_u=Aw(kU);
	double sindphi=sin(Aw(kPhiU));
	double cosdphi=cos(Aw(kPhiU));
	
	// Components of rotation matrix for converting global to local coords.
	double cospsi=cosa*cosdphi+sina*sindphi;
	double sinpsi=sina*cosdphi-cosa*sindphi;
	
	// x,y and tx,ty in local coordinate system	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	double upred_wire_plane=x*cospsi-y*sinpsi;
	double vpred_wire_plane=x*sinpsi+y*cospsi;
	double tu=tx*cospsi-ty*sinpsi;
	double tv=tx*sinpsi+ty*cospsi;
	
	// Variables for angle of incidence with respect to the z-direction in
	// the u-z plane
	double alpha=atan(tu);
	double cosalpha=cos(alpha);
	double sinalpha=sin(alpha);
	
	// Doca from wire
	double uwire=hits[id]->wire->u+delta_u;
	double d=(upred_wire_plane-uwire)*cosalpha;

	// Predicted avalanche position along the wire
	double vpred=vpred_wire_plane-tv*sinalpha*d;

	// predicted positions in two cathode planes' coordinate systems
	double phi_u=hits[id]->phi_u+A(kPhiU);
	double phi_v=hits[id]->phi_v+A(kPhiV);
	double cosphi_u=cos(phi_u);
	double sinphi_u=sin(phi_u);
	double cosphi_v=cos(phi_v);
	double sinphi_v=sin(phi_v);
	double vv=-vpred*sinphi_v+uwire*cosphi_v+A(kV);
	double vu=-vpred*sinphi_u+uwire*cosphi_u+A(kU);

	// Difference between measurements and predictions
	Mdiff(0)=hits[id]->u-vu;
	Mdiff(1)=hits[id]->v-vv;

	smoothed_updates[id].res=Mdiff;
	smoothed_updates[id].doca=fabs(d);
	
	smoothed_updates[id].drift=updates[id].drift;
	smoothed_updates[id].drift_time=updates[id].drift_time;
	smoothed_updates[id].S=Ss;
	smoothed_updates[id].C=Cs;
	smoothed_updates[id].V=updates[id].V-updates[id].H*dC*updates[id].H_T;
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
jerror_t DEventProcessor_dc_alignment::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
					   deque<trajectory_t>&trajectory,
					   vector<intersection_hit_t>&hits,
					   vector<wire_update_t>updates,
					   vector<wire_update_t>&smoothed_updates
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
	double cosa=hits[id].wire->udir.y();
	double sina=hits[id].wire->udir.x();
	
	// State vector
	double x=Ss(state_x);
	double y=Ss(state_y);
	double tx=Ss(state_tx);
	double ty=Ss(state_ty);
	
	// Get the aligment vector and error matrix for this layer
	unsigned int layer=hits[id].wire->layer-1;
	DMatrix2x2 E=fdc_alignments[layer].E;
	DMatrix2x1 A=fdc_alignments[layer].A;
	double delta_u=A(kU);
	double sindphi=sin(A(kPhiU));
	double cosdphi=cos(A(kPhiU));
	
	// Components of rotation matrix for converting global to local coords.
	double cospsi=cosa*cosdphi+sina*sindphi;
	double sinpsi=sina*cosdphi-cosa*sindphi;
	
	// x,y and tx,ty in local coordinate system	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	// (without alignment offsets)
	double upred=x*cospsi-y*sinpsi;
	double tu=tx*cospsi-ty*sinpsi;
	
	// Variables for angle of incidence with respect to the z-direction in
	// the u-z plane
	double alpha=atan(tu);
	double cosalpha=cos(alpha);
	
	// Smoothed residuals
	double uwire=hits[id].wire->u+delta_u;
	double d=(upred-uwire)*cosalpha;
	smoothed_updates[id].ures=(d>0?1.:-1.)*updates[id].drift-d;
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
  //printf("--------\n");
  for (unsigned int m=max-1;m>0;m--){
    if (trajectory[m].h_id==0){
      A=trajectory[m].Ckk*JT*C.Invert();
      Ss=trajectory[m].Skk+A*(Ss-S);
      Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    else{
      unsigned int id=trajectory[m].h_id-1;
      smoothed_updates[id].used_in_fit=false;
      //printf("%d:%d used ? %d\n",m,id,updates[id].used_in_fit);
      if (updates[id].used_in_fit){
	smoothed_updates[id].used_in_fit=true;

	A=updates[id].C*JT*C.Invert();
	dC=A*(Cs-C)*A.Transpose();
	Ss=updates[id].S+A*(Ss-S);
	Cs=updates[id].C+dC;
	
	// CDC index and wire position variables
	const DCDCWire *wire=hits[id]->wire;
	DVector3 origin=wire->origin;
	DVector3 wdir=wire->udir;
	
	unsigned int ring=hits[id]->wire->ring-1;
	unsigned int straw=hits[id]->wire->straw-1;
	UpdateWireOriginAndDir(ring,straw,origin,wdir);
	
	// doca using smoothed state vector
	double d=FindDoca(trajectory[m].z,Ss,wdir,origin);
	smoothed_updates[id].ring_id=ring;
	smoothed_updates[id].straw_id=straw;
	smoothed_updates[id].doca=d;
	smoothed_updates[id].res=updates[id].drift-d;	
	smoothed_updates[id].drift=updates[id].drift;
	smoothed_updates[id].drift_time=updates[id].drift_time;
	smoothed_updates[id].S=Ss;
	smoothed_updates[id].C=Cs;
	smoothed_updates[id].V=updates[id].V-updates[id].H*dC*updates[id].H_T;
	smoothed_updates[id].z=updates[id].z;

	// Reset h_id for this position along the reference trajectory
	trajectory[m].h_id=0;
      }
      else{
	A=trajectory[m].Ckk*JT*C.Invert();
	Ss=trajectory[m].Skk+A*(Ss-S);
	Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
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
  double V=1.15*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor

  for (unsigned int i=0;i<updates.size();i++){
    updates[i].used_in_fit=false;
  }

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

  double doca2=0.;
 
  // CDC index and wire position variables
  unsigned int cdc_index=hits.size()-1;
  bool more_hits=true;
  const DCDCWire *wire=hits[cdc_index]->wire;
  DVector3 origin=wire->origin;
  double z0=origin.z();
  double vz=wire->udir.z();
  DVector3 wdir=(1./vz)*wire->udir;

  // Wire offsets
  unsigned int ring=wire->ring-1;
  unsigned int straw=wire->straw-1;
  UpdateWireOriginAndDir(ring,straw,origin,wdir);

  DVector3 wirepos=origin+(trajectory[0].z-z0)*wdir;

  /// compute initial doca^2 to first wire
  double dx=S(state_x)-wirepos.X();
  double dy=S(state_y)-wirepos.Y();
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
    wirepos=origin+(trajectory[k].z-z0)*wdir;

    // New doca^2
    dx=S(state_x)-wirepos.X();
    dy=S(state_y)-wirepos.Y();
    doca2=dx*dx+dy*dy;

    if (doca2>old_doca2 && more_hits){

      // zero-position and direction of line describing particle trajectory
      double tx=S(state_tx),ty=S(state_ty);
      DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
      DVector3 tdir(tx,ty,1.);

      // Find the true doca to the wire
      DVector3 diff=pos0-origin;
      double dx0=diff.x(),dy0=diff.y();
      double wdir_dot_diff=diff.Dot(wdir);
      double tdir_dot_diff=diff.Dot(tdir);
      double tdir_dot_wdir=tdir.Dot(wdir);
      double tdir2=tdir.Mag2();
      double wdir2=wdir.Mag2();
      double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
      double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
      double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
      double scale=1./D;
      double s=scale*N;
      double t=scale*N1;
      diff+=s*tdir-t*wdir;
      double d=diff.Mag();

      // The next measurement and its variance
      double tdrift=hits[cdc_index]->tdrift-mT0-trajectory[k].t;
      double dmeas=0.39; 
      if (timebased){
	double drift_var=cdc_variance(tdrift);
	dmeas=cdc_drift_distance(tdrift);
	V=anneal_factor*drift_var;

	//printf("t0 %f t %f d %f %f V %f\n",mT0,hits[cdc_index]->tdrift,dmeas,d,V);
      }

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
 
      H(state_x)=H_T(state_x)=diffx*one_over_d;
      H(state_y)=H_T(state_y)=diffy*one_over_d;

      double wx=wdir.x(),wy=wdir.y();

      double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
      double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
      double dtdtx=scale*(dN1dtx-t*dDdtx);

      double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
      double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
      double dtdty=scale*(dN1dty-t*dDdty);

      double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
      double dsdtx=scale*(dNdtx-s*dDdtx);

      double dNdty=wy*wdir_dot_diff-wdir2*dy0;
      double dsdty=scale*(dNdty-s*dDdty);
      
      H(state_tx)=H_T(state_tx)
	=one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
		     +diffz*(dsdtx-dtdtx));
      H(state_ty)=H_T(state_ty)
	=one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
		     +diffz*(dsdty-dtdty));

      // Matrices to rotate alignment error matrix into measurement space
      DMatrix1x4 G;
      DMatrix4x1 G_T;      
      ComputeGMatrices(s,t,scale,tx,ty,tdir2,one_over_d,wx,wy,wdir2,tdir_dot_wdir,
		       tdir_dot_diff,wdir_dot_diff,dx0,dy0,diffx,diffy,diffz,
		       G,G_T);      
      
      // inverse of variance including prediction
      DMatrix4x4 E=cdc_alignments[ring][straw].E;
      double Vtemp=V+G*E*G_T;
      double InvV=1./(Vtemp+H*C*H_T);
      
      // Compute Kalman gain matrix
      K=InvV*(C*H_T);

      // Update state vector covariance matrix
      DMatrix4x4 Ctest=C-K*(H*C);

      //C.Print();
      //K.Print();
      //Ctest.Print();

      // Check that Ctest is positive definite
      if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0)
	{
	C=Ctest;

	// Update the state vector 
	//S=S+res*K;
	S+=res*K;

	// Compute new residual 
	d=FindDoca(trajectory[k].z,S,wdir,origin);
	res=dmeas-d;
	
	// Update chi2 for this segment
	Vtemp-=H*C*H_T;
	chi2+=res*res/Vtemp;
	ndof++;	
      }
      else{
	//	_DBG_ << "Bad C!" << endl;
	return VALUE_OUT_OF_RANGE;
      }

      updates[cdc_index].S=S;
      updates[cdc_index].C=C;
      updates[cdc_index].drift=dmeas;
      updates[cdc_index].drift_time=tdrift;
      updates[cdc_index].doca=d;
      updates[cdc_index].res=res;
      updates[cdc_index].V=Vtemp;
      updates[cdc_index].H_T=H_T;
      updates[cdc_index].H=H;
      updates[cdc_index].z=trajectory[k].z;
      updates[cdc_index].used_in_fit=true;

      trajectory[k].h_id=cdc_index+1;

      // move to next cdc hit
      if (cdc_index>0){
	cdc_index--;

	//New wire position
	wire=hits[cdc_index]->wire;
	origin=wire->origin;
	vz=wire->udir.z();
	wdir=(1./vz)*wire->udir;

	ring=hits[cdc_index]->wire->ring-1;
	straw=hits[cdc_index]->wire->straw-1;
	UpdateWireOriginAndDir(ring,straw,origin,wdir);

	wirepos=origin+((trajectory[k].z-z0))*wdir;
	
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
  DMatrix2x2 V(0.0075*anneal_factor,0.,0.,0.0075*anneal_factor);  // Measurement variance 
  DMatrix2x2 Vtemp,InvV;
  DMatrix2x1 Mdiff;
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
      DMatrix2x1 Aw=fdc_alignments[layer].A;
      double delta_u=Aw(kU);
      double sindphi=sin(Aw(kPhiU));
      double cosdphi=cos(Aw(kPhiU));
      
      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;
      
      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double vpred_wire_plane=y*cospsi+x*sinpsi;
      double upred_wire_plane=x*cospsi-y*sinpsi;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi+ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double cos2_alpha=cosalpha*cosalpha;
      double sinalpha=sin(alpha);
      double sin2_alpha=sinalpha*sinalpha;
      
      // Alignment parameters for cathode planes
      DMatrix4x4 E=fdc_cathode_alignments[layer].E;
      DMatrix4x1 A=fdc_cathode_alignments[layer].A;

      // Difference between measurement and projection
      for (int m=trajectory[k].num_hits-1;m>=0;m--){
	unsigned int my_id=id+m;
	double uwire=hits[my_id]->wire->u+delta_u;
	// (signed) distance of closest approach to wire
        double doca=(upred_wire_plane-uwire)*cosalpha;

	// Predicted avalanche position along the wire
	double vpred=vpred_wire_plane-tv*sinalpha*doca;

	// predicted positions in two cathode planes' coordinate systems
	double phi_u=hits[my_id]->phi_u+A(kPhiU);
	double phi_v=hits[my_id]->phi_v+A(kPhiV);
	double cosphi_u=cos(phi_u);
	double sinphi_u=sin(phi_u);
	double cosphi_v=cos(phi_v);
	double sinphi_v=sin(phi_v);
	double vv=-vpred*sinphi_v+uwire*cosphi_v+A(kV);
	double vu=-vpred*sinphi_u+uwire*cosphi_u+A(kU);

	// Difference between measurements and predictions
	Mdiff(0)=hits[my_id]->u-vu;
	Mdiff(1)=hits[my_id]->v-vv;

	// Start filling the update vector
	updates[my_id].drift_time=hits[my_id]->time-trajectory[k].t-mT0;
	
	// Matrix for transforming from state-vector space to measurement space
	double temp2=tv*sinalpha*cosalpha;
	double dvdy=cospsi+sinpsi*temp2;
	double dvdx=sinpsi-cospsi*temp2;
	
        H_T(state_x,0)=-dvdx*sinphi_u;
	H_T(state_y,0)=-dvdy*sinphi_u;
	H_T(state_x,1)=-dvdx*sinphi_v;
	H_T(state_y,1)=-dvdy*sinphi_v;
       
        double cos2_minus_sin2=cos2_alpha-sin2_alpha;
        double doca_cosalpha=doca*cosalpha;
	double dvdtx=-doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2);
	double dvdty=-doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2);

        H_T(state_tx,0)=-dvdtx*sinphi_u;
        H_T(state_ty,0)=-dvdty*sinphi_u;
	H_T(state_tx,1)=-dvdtx*sinphi_v;
        H_T(state_ty,1)=-dvdty*sinphi_v;

        // Matrix transpose H_T -> H
        H(0,state_x)=H_T(state_x,0);
        H(0,state_y)=H_T(state_y,0);
        H(0,state_tx)=H_T(state_tx,0);
        H(0,state_ty)=H_T(state_ty,0);
	H(1,state_x)=H_T(state_x,1);
        H(1,state_y)=H_T(state_y,1);
        H(1,state_tx)=H_T(state_tx,1);
        H(1,state_ty)=H_T(state_ty,1);

	updates[my_id].H=H;
	updates[my_id].H_T=H_T;
		
	// Matrices to rotate alignment error matrix into measurement space
	DMatrix2x4 G;
	DMatrix4x2 G_T;

	G_T(kU,0)=1.;
	G_T(kPhiU,0)=-vpred*cosphi_u-uwire*sinphi_u;
	G_T(kV,1)=1.;
	G_T(kPhiV,1)=-vpred*cosphi_v-uwire*sinphi_v;
	
	// G-matrix transpose
	G(0,kU)=G_T(kU,0);
	G(0,kPhiU)=G_T(kPhiU,0);	
	G(1,kV)=G_T(kV,1);
	G(1,kPhiV)=G_T(kPhiV,1);

	Vtemp=V+G*E*G_T;

	// Variance for this hit
	InvV=(Vtemp+H*C*H_T).Invert();
	
	// Compute Kalman gain matrix
	K=(C*H_T)*InvV;
	
	// Update the state vector 
	S+=K*Mdiff;

	// Update state vector covariance matrix
	C=C-K*(H*C);    

	// Update the filtered measurement covariane matrix and put results in 
	// update vector
	DMatrix2x2 RC=Vtemp-H*C*H_T;
	updates[my_id].res=Mdiff-H*K*Mdiff;
	updates[my_id].V=RC;		
	updates[my_id].S=S;
	updates[my_id].C=C;
	
	chi2+=RC.Chi2(updates[my_id].res);
	ndof+=2;
      }

    }

  }
  //  chi2*=anneal_factor;
  ndof-=4;

  return NOERROR;
}



// Perform Kalman Filter for the current trajectory
jerror_t 
DEventProcessor_dc_alignment::KalmanFilter(double anneal_factor,
					DMatrix4x1 &S,DMatrix4x4 &C,
			       vector<intersection_hit_t>&hits,
			       deque<trajectory_t>&trajectory,
			       vector<wire_update_t>&updates,
			       double &chi2,unsigned int &ndof){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  double V=0.020833;  // Measurement variance 
  double Vtemp,Mdiff,InvV;
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
      
      double cosa=hits[id].wire->udir.y();
      double sina=hits[id].wire->udir.x();
      
      // State vector
      double x=S(state_x);
      double y=S(state_y);
      double tx=S(state_tx);
      double ty=S(state_ty);
      if (isnan(x) || isnan(y)) return UNRECOVERABLE_ERROR;
      
      // Get the alignment vector and error matrix for this layer
      unsigned int layer=hits[id].wire->layer-1;
      DMatrix2x2 E=fdc_alignments[layer].E;
      DMatrix2x1 A=fdc_alignments[layer].A;
      double delta_u=A(kU);
      double sindphi=sin(A(kPhiU));
      double cosdphi=cos(A(kPhiU));
      
      // Components of rotation matrix for converting global to local coords.
      double cospsi=cosa*cosdphi+sina*sindphi;
      double sinpsi=sina*cosdphi-cosa*sindphi;
      
      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      // (without alignment offsets)
      double upred=x*cospsi-y*sinpsi;
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
	double uwire=hits[my_id].wire->u+delta_u;

	// Find drift distance
	double drift_time=hits[my_id].hit->t-trajectory[k].t-mT0;
	updates[my_id].drift_time=drift_time;
	updates[my_id].t=trajectory[k].t;

	double du=upred-uwire;
	double d=du*cosalpha;
	double sign=(du>0)?1.:-1.;
	
	// Difference between measured and predicted vectors
	// assume the track passes through the center of the cell
	double drift=0.25;
	if (USE_DRIFT_TIMES){
	  drift=0.;
	  if (drift_time>0){
	    drift=fdc_drift_distance(drift_time);

	    //V=0.0004+0.020433*(anneal_factor/1000.);
	    double sigma=0.0135-3.98e-4*drift_time+5.62e-6*drift_time*drift_time;
	    V=anneal_factor*sigma*sigma;
	  }
	}
	Mdiff=sign*drift-d;
	updates[my_id].drift=drift;

	// Matrix for transforming from state-vector space to measurement space
	double sinalpha_cosalpha=sinalpha*cosalpha;
	H_T(state_x)=cospsi*cosalpha;
	H_T(state_y)=-sinpsi*cosalpha;
	
	double temp=d*sinalpha_cosalpha;
	H_T(state_tx)=-temp*cospsi;
	H_T(state_ty)=+temp*sinpsi;

	// H-matrix transpose
	H(state_x)=H_T(state_x);
	H(state_y)=H_T(state_y);

	H(state_tx)=H_T(state_tx);
	H(state_ty)=H_T(state_ty);
	
	updates[my_id].H=H;
	updates[my_id].H_T=H_T;
		
	// Matrices to rotate alignment error matrix into measurement space
	DMatrix1x2 G;
	DMatrix2x1 G_T;
	
	G_T(kU)=-cosalpha;
	G_T(kPhiU)=cosalpha*(x*sinpsi+y*cospsi-tv*d);
	
	// G-matrix transpose
	G(kU)=G_T(kU);
	G(kPhiU)=G_T(kPhiU);

	Vtemp=V+G*E*G_T;

	// Variance for this hit
	InvV=1./(Vtemp+H*C*H_T);

	// Compute Kalman gain matrix
	K=InvV*(C*H_T);
	  
	// Update the state vector 
	S+=Mdiff*K;
	updates[my_id].S=S;
	
	// Update state vector covariance matrix
	C=C-K*(H*C);    
	updates[my_id].C=C;
	
	// Update chi2 for this trajectory 
	x=S(state_x);
	y=S(state_y);
	tx=S(state_tx);
	ty=S(state_ty);
	upred=x*cospsi-y*sinpsi;
	tu=tx*cospsi-ty*sinpsi;
	
	// Variables for angle of incidence with respect to the z-direction in
	// the u-z plane
	alpha=atan(tu);
	cosalpha=cos(alpha);
	du=upred-uwire;
	d=du*cosalpha;
	sinalpha=sin(alpha);
	
	sign=(du>0)?1.:-1.;
	Mdiff=sign*drift-d;
	
	double RC=Vtemp-H*C*H_T;
	updates[my_id].ures=Mdiff;
	updates[my_id].R=RC;
	
	chi2+=Mdiff*Mdiff/RC;
	ndof++;
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
  unsigned int numsteps=0;
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
    numsteps++;
  }while (S(state_y)>min_y && numsteps<MAX_STEPS);

  if (trajectory.size()<2) return UNRECOVERABLE_ERROR;
  //if (true)
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
    dz=1.1;
    
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
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;  
      temp.h_id=0;	  
      temp.num_hits=0;
      if (new_z>zhit){
	dz=zhit-z;
	new_z=zhit;
	temp.h_id=i+1;
	temp.num_hits=1;
	done=true;
      }
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
  dz=1.1;
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

// Reference trajectory for the track
jerror_t DEventProcessor_dc_alignment
::SetReferenceTrajectory(double z,DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 vector<intersection_hit_t>&hits){
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
  for (unsigned int i=0;i<hits.size();i++){  
    zhit=hits[i].wire->origin.z();
    dz=1.1;

    if (fabs(zhit-old_zhit)<EPS){
      trajectory[0].num_hits++;
      continue;
    }
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      
      trajectory_t temp;   
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;  
      temp.h_id=0;	  
      temp.num_hits=0;
      if (new_z>zhit){
	new_z=zhit;
	dz=zhit-z;
	temp.h_id=i+1;
	temp.num_hits=1;
	done=true;
      }
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
  dz=1.1;
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

void DEventProcessor_dc_alignment::UpdateWireOriginAndDir(unsigned int ring,
							  unsigned int straw,
							  DVector3 &origin,
							  DVector3 &wdir){
  double zscale=75.0/wdir.z();
  DVector3 upstream=origin-zscale*wdir;
  DVector3 downstream=origin+zscale*wdir;
  DVector3 du(cdc_alignments[ring][straw].A(k_dXu),
	      cdc_alignments[ring][straw].A(k_dYu),0.);
  DVector3 dd(cdc_alignments[ring][straw].A(k_dXd),
	      cdc_alignments[ring][straw].A(k_dYd),0.);
  upstream+=du;
  downstream+=dd;
  
  origin=0.5*(upstream+downstream);
  wdir=downstream-upstream;
  wdir.SetMag(1.);
}
							  


jerror_t 
DEventProcessor_dc_alignment::FindOffsets(vector<const DCDCTrackHit*>&hits,
					  vector<cdc_update_t>&updates){
  for (unsigned int i=0;i<updates.size();i++){
    if (updates[i].used_in_fit==true){
      // wire data
      const DCDCWire *wire=hits[i]->wire;
      DVector3 origin=wire->origin;
      DVector3 wdir=wire->udir;
      
      unsigned int ring=wire->ring-1;
      unsigned int straw=wire->straw-1;
      UpdateWireOriginAndDir(ring,straw,origin,wdir);

      // zero-position and direction of line describing particle trajectory
      double tx=updates[i].S(state_tx),ty=updates[i].S(state_ty);
      DVector3 pos0(updates[i].S(state_x),updates[i].S(state_y),updates[i].z);
      DVector3 diff=pos0-origin;
      double dx0=diff.x(),dy0=diff.y();
      DVector3 tdir(tx,ty,1.);
      double wdir_dot_diff=diff.Dot(wdir);
      double tdir_dot_diff=diff.Dot(tdir);
      double tdir_dot_wdir=tdir.Dot(wdir);
      double tdir2=tdir.Mag2();
      double wdir2=wdir.Mag2();
      double wx=wdir.x(),wy=wdir.y();
      double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
      double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
      double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
      double scale=1./D;
      double s=scale*N;
      double t=scale*N1;
      diff+=s*tdir-t*wdir;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
      double one_over_d=1./diff.Mag();
      
      // Matrices to rotate alignment error matrix into measurement space
      DMatrix1x4 G;
      DMatrix4x1 G_T;    
      ComputeGMatrices(s,t,scale,tx,ty,tdir2,one_over_d,wx,wy,wdir2,tdir_dot_wdir,
		       tdir_dot_diff,wdir_dot_diff,dx0,dy0,diffx,diffy,diffz,
		       G,G_T);      
      
      // Offset error matrix
      DMatrix4x4 E=cdc_alignments[ring][straw].E;
      
      // Inverse error
      double InvV=1./updates[i].V;
      
      // update the alignment vector and covariance
      DMatrix4x1 Ka=InvV*(E*G_T);
      DMatrix4x1 dA=updates[i].res*Ka;
      DMatrix4x4 Etemp=E-Ka*G*E; 
      //dA.Print();
      //Etemp.Print();
      if (Etemp(0,0)>0 && Etemp(1,1)>0 && Etemp(2,2)>0&&Etemp(3,3)>0.){  
	//cdc_alignments[ring][straw].A.Print();
	//dA.Print();
	//Etemp.Print();
	
	cdc_alignments[ring][straw].E=Etemp;
	cdc_alignments[ring][straw].A+=dA;	  
      }
    }
  }

  return NOERROR;
}


jerror_t 
DEventProcessor_dc_alignment::FindOffsets(vector<const DFDCPseudo *>&hits,
				      vector<update_t>&smoothed_updates){
  DMatrix2x4 G;//matrix relating alignment vector to measurement coords
  DMatrix4x2 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();

  for (unsigned int i=0;i<num_hits;i++){ 
    // Get the cathode planes aligment vector and error matrix for this layer
    unsigned int layer=hits[i]->wire->layer-1;
    DMatrix4x1 A=fdc_cathode_alignments[layer].A;
    DMatrix4x4 E=fdc_cathode_alignments[layer].E;

    // Rotation of wire planes
    double cosa=hits[i]->wire->udir.y();
    double sina=hits[i]->wire->udir.x();
      
    // State vector
    DMatrix4x1 S=smoothed_updates[i].S;
    double x=S(state_x);
    double y=S(state_y);
    double tx=S(state_tx);
    double ty=S(state_ty);
    if (isnan(x) || isnan(y)) return UNRECOVERABLE_ERROR;
      
    // Get the wire plane alignment vector and error matrix for this layer
    DMatrix2x1 Aw=fdc_alignments[layer].A;
    double delta_u=Aw(kU);
    double sindphi=sin(Aw(kPhiU));
    double cosdphi=cos(Aw(kPhiU));
    
    // Components of rotation matrix for converting global to local coords.
    double cospsi=cosa*cosdphi+sina*sindphi;
    double sinpsi=sina*cosdphi-cosa*sindphi;
      
    // x,y and tx,ty in local coordinate system	
    // To transform from (x,y) to (u,v), need to do a rotation:
    //   u = x*cosa-y*sina
    //   v = y*cosa+x*sina
    // (without alignment offsets)
    double vpred_wire_plane=y*cospsi+x*sinpsi;
    double upred_wire_plane=x*cospsi-y*sinpsi;
    double tu=tx*cospsi-ty*sinpsi;
    double tv=tx*sinpsi+ty*cospsi;
    double alpha=atan(tu);
    double cosalpha=cos(alpha);
    double sinalpha=sin(alpha);

    // Wire position in wire-plane local coordinate system
    double uwire=hits[i]->wire->u+delta_u;
    // (signed) distance of closest approach to wire
    double doca=(upred_wire_plane-uwire)*cosalpha;

    // Predicted avalanche position along the wire
    double vpred=vpred_wire_plane-tv*sinalpha*doca;

    // Matrices to rotate alignment error matrix into measurement space
    DMatrix2x4 G;
    DMatrix4x2 G_T;
    
    double phi_u=hits[i]->phi_u+A(kPhiU);
    double phi_v=hits[i]->phi_v+A(kPhiV);
     
    G_T(kU,0)=1.;
    G_T(kPhiU,0)=-vpred*cos(phi_u)-uwire*sin(phi_u);
    G_T(kV,1)=1.;
    G_T(kPhiV,1)=-vpred*cos(phi_v)-uwire*sin(phi_v);
    
    // update the alignment vector and covariance
    DMatrix4x2 Ka=(E*G_T)*smoothed_updates[i].V.Invert();
    DMatrix4x1 dA=Ka*smoothed_updates[i].res;
    DMatrix4x4 Etemp=E-Ka*G*E;
    if (Etemp(0,0)>0 && Etemp(1,1)>0 && Etemp(2,2)>0 && Etemp(3,3)>0){
      fdc_cathode_alignments[layer].E=Etemp;
      fdc_cathode_alignments[layer].A=A+dA;	  
    }
    else {
      /*
      printf("-------t= %f\n",smoothed_updates[i].drift_time);
      E.Print();
      Etemp.Print();
      */
    }
  }

  return NOERROR;
}

jerror_t 
DEventProcessor_dc_alignment::FindOffsets(vector<intersection_hit_t>&hits,
				      vector<wire_update_t>&smoothed_updates){
  DMatrix1x2 G;//matrix relating alignment vector to measurement coords
  DMatrix2x1 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();


  for (unsigned int i=0;i<num_hits;i++){ 
    double x=smoothed_updates[i].S(state_x);
    double y=smoothed_updates[i].S(state_y);
    double tx=smoothed_updates[i].S(state_tx);
    double ty=smoothed_updates[i].S(state_ty);
    
    double cosa=hits[i].wire->udir.y();
    double sina=hits[i].wire->udir.x();
      
    // Get the aligment vector and error matrix for this layer
    unsigned int layer=hits[i].wire->layer-1;
    DMatrix2x1 A=fdc_alignments[layer].A;
    DMatrix2x2 E=fdc_alignments[layer].E;
    double delta_u=A(kU);
    double sindphi=sin(A(kPhiU));
    double cosdphi=cos(A(kPhiU));
    
    // Components of rotation matrix for converting global to local coords.
    double cospsi=cosa*cosdphi+sina*sindphi;
    double sinpsi=sina*cosdphi-cosa*sindphi;
    
    // x,y and tx,ty in local coordinate system	
    // To transform from (x,y) to (u,v), need to do a rotation:
    //   u = x*cosa-y*sina
    //   v = y*cosa+x*sina
    // (without alignment offsets)
    double uwire=hits[i].wire->u+delta_u;
    double upred=x*cospsi-y*sinpsi;
    double tu=tx*cospsi-ty*sinpsi;
    double tv=tx*sinpsi+ty*cospsi;
    double du=upred-uwire;
    
    // Variables for angle of incidence with respect to the z-direction in
    // the u-z plane
    double alpha=atan(tu);
    double cosalpha=cos(alpha);
    
    // Transform from alignment vector coords to measurement coords
    G_T(kU)=-cosalpha;

    double d=du*cosalpha;
    G_T(kPhiU)=cosalpha*(x*sinpsi+y*cospsi-tv*d);
    
    // G-matrix transpose
    G(kU)=G_T(kU);
    G(kPhiU)=G_T(kPhiU);
    
    // Inverse of error "matrix"
    double InvV=1./smoothed_updates[i].R;
    
    // update the alignment vector and covariance
    DMatrix2x1 Ka=InvV*(E*G_T);
    DMatrix2x1 dA=smoothed_updates[i].ures*Ka;
    DMatrix2x2 Etemp=E-Ka*G*E;
    if (Etemp(0,0)>0 && Etemp(1,1)>0){
      fdc_alignments[layer].E=Etemp;
      fdc_alignments[layer].A=A+dA;
    }
    else {
      /*
      printf("-------t= %f\n",smoothed_updates[i].drift_time);
      E.Print();
      Etemp.Print();
      */
    }
  }

  return NOERROR;
}

// Match tracks in the cdc to the outer detectors
bool DEventProcessor_dc_alignment::MatchOuterDetectors(const cdc_track_t &track,
				        vector<const DFCALShower *>&fcalshowers,
			        	vector<const DBCALShower *>&bcalshowers,
						       DMatrix4x1 &S){
  
  double ux=track.dir.x();
  double uy=track.dir.y();
  DVector3 pos0=track.axial_hits[0]->wire->origin;
  double x0=pos0.x();
  double y0=pos0.y();

  for (unsigned int i=0;i<fcalshowers.size();i++){


  }

  // Keep list of matches
  vector<bcal_match_t>matching_bcals;

  for (unsigned int i=0;i<bcalshowers.size();i++){
    double x=bcalshowers[i]->x,y=bcalshowers[i]->y;
    double s=(x-x0)*ux+(y-y0)*uy;
    double x1=x0+s*ux;
    double y1=y0+s*uy;
    double dx=x1-x;
    double dy=y1-y;

    Hbcalmatchxy->Fill(dx,dy);

    if (fabs(dx)<2.7 && fabs(dy)<0.6){
      bcal_match_t temp;
      temp.xtrack=x1;
      temp.ytrack=y1;
      temp.match=bcalshowers[i];
      matching_bcals.push_back(temp);
    }
  }
  if (matching_bcals.size()>0){
    sort(matching_bcals.begin(),matching_bcals.end(),bcal_cmp);
    
    mT0=matching_bcals[0].match->t;
    mOuterZ=matching_bcals[0].match->z;

    S(state_x)=matching_bcals[0].match->x;
    S(state_y)=matching_bcals[0].match->y;
      
    if (COSMICS){
      if (matching_bcals.size()!=2 || 
	  matching_bcals[0].match->y*matching_bcals[1].match->y>0){ 
	if (GuessForStateVector(track,S)==NOERROR)return true;

	return false;
      }
      
      // Estimate for beta
      double dx=matching_bcals[0].match->x-matching_bcals[1].match->x;
      double dy=matching_bcals[0].match->y-matching_bcals[1].match->y;
      double dz=matching_bcals[0].match->z-matching_bcals[1].match->z;
      double beta=sqrt(dx*dx+dy*dy+dz*dz)
	/(29.98*(matching_bcals[1].match->t-matching_bcals[0].match->t));
      
      Hbeta->Fill(beta);
      
      // Use bcal hits to estimate slopes
      S(state_tx)=dx/dz;
      S(state_ty)=dy/dz;
    }
    else{
      if (GuessForStateVector(track,S)==NOERROR) return true;
      return false;
    }
    return true;
  }
  return false;
}




// Routine to match preliminary track to outer detectors
bool DEventProcessor_dc_alignment::MatchOuterDetectors(vector<const DFCALShower *>&fcalshowers,
		      vector<const DBCALShower *>&bcalshowers,
		      const DMatrix4x1 &S){
  //compute tangent of dip angle and related angular quantities
  double tanl=1./sqrt(S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));  
  double sinl=sin(atan(tanl));

  // First match to FCAL 
  double dz=0.;
  double drmin=1000.;
  for (unsigned int i=0;i<fcalshowers.size();i++){
    double fcal_z=fcalshowers[i]->getPosition().z();
    double mydz=fcal_z-endplate_z;
    double x=S(state_x)+mydz*S(state_tx);
    double y=S(state_y)+mydz*S(state_ty);
    double dx=fcalshowers[i]->getPosition().x()-x;
    double dy=fcalshowers[i]->getPosition().y()-y;
    double dr=sqrt(dx*dx+dy*dy);
      
    if (dr<drmin){
      drmin=dr;
      dz=mydz;
      mOuterTime=fcalshowers[i]->getTime();
      mOuterZ=fcal_z;
    }
    
  }
  Hfcalmatch->Fill(drmin);
  
  if (drmin<4.){   
    // Estimate for t0 at the beginning of track assuming particle is 
    // moving at the speed of light
    //mT0=mOuterTime-dz/(29.98*sinl);
    mT0=mOuterTime-(mOuterZ-endplate_z)/(29.98*sinl);

    return true;
  }
  // Disable matching to BCAL by setting "false" below
  else if (true/*false*/){
    // Match to BCAL
 
    // zero-position and direction of line describing particle trajectory
    DVector3 pos0(S(state_x),S(state_y),endplate_z);
    DVector3 vhat(S(state_tx),S(state_ty),1.);
    vhat.SetMag(1.);
    
    // Keep list of matches
    vector<bcal_match_t>matching_bcals;

    // loop over the showers
    for (unsigned int i=0;i<bcalshowers.size();i++){
      // Match in x and y
      DVector3 pos1(bcalshowers[i]->x,bcalshowers[i]->y,bcalshowers[i]->z);
      double s=(pos1-pos0).Dot(vhat);
      pos0+=s*vhat;
      DVector3 diff=pos1-pos0;

      double dr=diff.Perp();
      double dz=diff.z();
      Hbcalmatch->Fill(dz,dr);
      if (dr<2.0 && fabs(dz)<10.0){
	bcal_match_t temp;
	temp.ztrack=pos0.z();
	temp.match=bcalshowers[i];
	matching_bcals.push_back(temp);
      }	
    }
    if (matching_bcals.size()>0){
      sort(matching_bcals.begin(),matching_bcals.end(),bcal_cmp);

      mT0=matching_bcals[0].match->t;
      mOuterZ=matching_bcals[0].ztrack;

      
      if (COSMICS){
	if (matching_bcals.size()!=2) return false;
	if (matching_bcals[0].match->y*matching_bcals[1].match->y>0) return false;
	
	// Estimate for beta
	double dx=matching_bcals[0].match->x-matching_bcals[1].match->x;
	double dy=matching_bcals[0].match->y-matching_bcals[1].match->y;
	double dz=matching_bcals[0].match->z-matching_bcals[1].match->z;
	double beta=sqrt(dx*dx+dy*dy+dz*dz)
	  /(29.98*(matching_bcals[1].match->t-matching_bcals[0].match->t));

	Hbeta->Fill(beta);

      }
      else{
	// assume particle moving at speed of light
	//mT0-=(mOuterZ-ztarg)/(29.98*sinl); 
	mT0-=(mOuterZ-endplate_z)/(29.98*sinl); 
      }

      return true;
    }
  }
 
  return false;
}



// Match a CDC hit with a line starting at pos0 going in the vhat direction
bool DEventProcessor_dc_alignment::MatchCDCHit(const DVector3 &vhat,
					     const DVector3 &pos0,
					     const DCDCTrackHit *hit){  
  DVector3 pos1=hit->wire->origin;
  DVector3 uhat=hit->wire->udir;
  DVector3 diff=pos1-pos0;
  double vhat_dot_uhat=vhat.Dot(uhat);
  double scale=1./(1.-vhat_dot_uhat*vhat_dot_uhat);
  double s=scale*(vhat_dot_uhat*diff.Dot(vhat)-diff.Dot(uhat));
  double t=scale*(diff.Dot(vhat)-vhat_dot_uhat*diff.Dot(uhat));
  double d=(diff+s*uhat-t*vhat).Mag();

  if (d<CDC_MATCH_RADIUS) return true;

  return false;
}


// Link axial segments together to form track candidates and match to stereo 
// hits
jerror_t 
DEventProcessor_dc_alignment::LinkSegments(vector<cdc_segment_t>&axial_segments,
					   vector<bool>&used_axial,
					   vector<const DCDCTrackHit *>&axial_hits,
					   vector<const DCDCTrackHit *>&stereo_hits,
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

	  {
	    pthread_mutex_lock(&mutex);
	    Hcdcmatch_stereo->Fill(d);
	    pthread_mutex_unlock(&mutex);
	  }

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
      //  Position of the first axial wire in the track  
      pos0=mytrack.axial_hits[0]->wire->origin;

      // Grab axial hits not associated with segments
      bool got_match=false;
      for (unsigned int j=0;j<axial_hits.size();j++){
	if (used_axial[j]==false){
	  if (MatchCDCHit(vhat,pos0,axial_hits[j])){
	    used_axial[j]=true;
	    mytrack.axial_hits.push_back(axial_hits[j]);
	    got_match=true;
	  }
	}
      }
      // Resort if we added axial hits and recompute direction vector
      if (got_match){
	sort(mytrack.axial_hits.begin(),mytrack.axial_hits.end(),cdc_hit_cmp);
	
	vhat=mytrack.axial_hits[mytrack.axial_hits.size()-1]->wire->origin
	  -mytrack.axial_hits[0]->wire->origin;
	vhat.SetMag(1.);   
	pos0=mytrack.axial_hits[0]->wire->origin;
      }
      
    
      // Now try to associate stereo hits with this track
      vector<unsigned int>used_stereo(stereo_hits.size());
      for (unsigned int j=0;j<stereo_hits.size();j++){
	if (used_stereo[j]==false){
	  if (MatchCDCHit(vhat,pos0,stereo_hits[j])){
	    used_stereo[j]=true;
	    mytrack.stereo_hits.push_back(stereo_hits[j]);
	  }
	}
      }
      size_t num_stereo=mytrack.stereo_hits.size();
      size_t num_axial=mytrack.axial_hits.size();
      if (num_stereo>0 && num_stereo+num_axial>4){
	mytrack.dir=vhat;
	LinkedSegments.push_back(mytrack);
      }
    }
  }

  return NOERROR;
}
// Compute initial guess for state vector (x,y,tx,ty) for a track in the CDC
// by fitting a line to the intersections between the line in the xy plane and 
// the stereo wires.
jerror_t
DEventProcessor_dc_alignment::GuessForStateVector(const cdc_track_t &track,
						DMatrix4x1 &S){
  // Parameters for line in x-y plane
  double vx=track.dir.x();
  double vy=track.dir.y();
  DVector3 pos0=track.axial_hits[0]->wire->origin;
  double xa=pos0.x();
  double ya=pos0.y();

  double sumv=0,sumx=0,sumy=0,sumz=0,sumxx=0,sumyy=0,sumxz=0,sumyz=0;
  for (unsigned int i=0;i<track.stereo_hits.size();i++){
    // Intersection of line in xy-plane with this stereo straw
    DVector3 origin_s=track.stereo_hits[i]->wire->origin;
    DVector3 dir_s=track.stereo_hits[i]->wire->udir;
    double ux_s=dir_s.x();
    double uy_s=dir_s.y();
    double dx=xa-origin_s.x();
    double dy=ya-origin_s.y();
    double s=(dx*vy-dy*vx)/(ux_s*vy-uy_s*vx);
    DVector3 pos1=origin_s+s*dir_s;
    double x=pos1.x(),y=pos1.y(),z=pos1.z();
    
    if (z>17.0 && z<167.0){ // Check for CDC dimensions
      sumv+=1.;
      sumx+=x;
      sumxx+=x*x;
      sumy+=y;
      sumyy+=y*y;
      sumz+=z;
      sumxz+=x*z;
      sumyz+=y*z;
    }
  }
  double xdenom=sumv*sumxz-sumx*sumz;
  if (fabs(xdenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double ydenom=sumv*sumyz-sumy*sumz;
  if (fabs(ydenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double xtemp=sumv*sumxx-sumx*sumx;
  double xslope=xtemp/xdenom;
  double ytemp=sumv*sumyy-sumy*sumy;
  double yslope=ytemp/ydenom;

  //  double z0x=(sumxx*sumz-sumx*sumxz)/xtemp;
  double z0y=(sumyy*sumz-sumy*sumyz)/ytemp;
  
  // Increment just beyond point largest in y
  double delta_z=(yslope>0)?0.5:-0.5;

  //Starting z position
  mOuterZ=z0y+ya/yslope+delta_z;

  S(state_x)=xa+xslope*delta_z;
  S(state_y)=ya+yslope*delta_z;
  S(state_tx)=xslope;
  S(state_ty)=yslope;

  return NOERROR;
}

// Compute distance of closest approach between two lines
double DEventProcessor_dc_alignment::FindDoca(double z,const DMatrix4x1 &S,
					      const DVector3 &wdir,
					      const DVector3 &origin){
  DVector3 pos(S(state_x),S(state_y),z);
  DVector3 diff=pos-origin;
  
  DVector3 uhat(S(state_tx),S(state_ty),1.);
  uhat.SetMag(1.); 
  DVector3 vhat=wdir;
  //  vhat.SetMag(1.);

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


// Compute matrices for rotating the aligment error matrix into the measurement
// space
void 
DEventProcessor_dc_alignment::ComputeGMatrices(double s,double t,double scale,
					       double tx,double ty,double tdir2,
					       double one_over_d,
					       double wx,double wy,double wdir2,
					       double tdir_dot_wdir,
					       double tdir_dot_diff,
					       double wdir_dot_diff,
					       double dx0,double dy0,
					       double diffx,double diffy,
					       double diffz,
					       DMatrix1x4 &G,DMatrix4x1 &G_T){
  double dsdDx=scale*(tdir_dot_wdir*wx-wdir2*tx);
  double dsdDy=scale*(tdir_dot_wdir*wy-wdir2*ty);
  
  double dNdvx=tx*wdir_dot_diff+tdir_dot_wdir*dx0-2.*wx*tdir_dot_diff;
  double dDdvx=2.*wx*tdir2-2.*tdir_dot_wdir*tx;
  double dsdvx=scale*(dNdvx-s*dDdvx);
  
  double dNdvy=ty*wdir_dot_diff+tdir_dot_wdir*dy0-2.*wy*tdir_dot_diff;
  double dDdvy=2.*wy*tdir2-2.*tdir_dot_wdir*ty;;
  double dsdvy=scale*(dNdvy-s*dDdvy);
  
  double dsddxu=-0.5*dsdDx-one_over_zrange*dsdvx;
  double dsddxd=-0.5*dsdDx+one_over_zrange*dsdvx; 
  double dsddyu=-0.5*dsdDy-one_over_zrange*dsdvy;
  double dsddyd=-0.5*dsdDy+one_over_zrange*dsdvy;

  double dtdDx=scale*(tdir2*wx-tdir_dot_wdir*tx);
  double dtdDy=scale*(tdir2*wy-tdir_dot_wdir*ty);
  
  double dN1dvx=tdir2*dx0-tdir_dot_diff*tx;
  double dtdvx=scale*(dN1dvx-t*dDdvx);

  double dN1dvy=tdir2*dy0-tdir_dot_diff*ty;
  double dtdvy=scale*(dN1dvy-t*dDdvy);
      
  double dtddxu=-0.5*dtdDx-one_over_zrange*dtdvx;
  double dtddxd=-0.5*dtdDx+one_over_zrange*dtdvx; 
  double dtddyu=-0.5*dtdDy-one_over_zrange*dtdvy;
  double dtddyd=-0.5*dtdDy+one_over_zrange*dtdvy;
  
  double t_over_zrange=one_over_zrange*t;
  G(k_dXu)=one_over_d*(diffx*(-0.5+tx*dsddxu+t_over_zrange-wx*dtddxu)
		       +diffy*(ty*dsddxu-wy*dtddxu)+diffz*(dsddxu-dtddxu));
  G(k_dXd)=one_over_d*(diffx*(-0.5+tx*dsddxd-t_over_zrange-wx*dtddxd)
		       +diffy*(ty*dsddxd-wy*dtddxd)+diffz*(dsddxd-dtddxd));
  G(k_dYu)=one_over_d*(diffx*(tx*dsddyu-wx*dtddyu)+diffz*(dsddyu-dtddyu)
		       +diffy*(-0.5+ty*dsddyu+t_over_zrange-wy*dtddyu));  
  G(k_dYd)=one_over_d*(diffx*(tx*dsddyd-wx*dtddyd)+diffz*(dsddyd-dtddyd)
		       +diffy*(-0.5+ty*dsddyd-t_over_zrange-wy*dtddyd));
  G_T(k_dXu)=G(k_dXu);
  G_T(k_dXd)=G(k_dXd);
  G_T(k_dYu)=G(k_dYu);
  G_T(k_dYd)=G(k_dYd);
}

// If the event viewer is available, grab parts of the hdview2 display and 
// overlay the results of the line fit on the tracking views.
void DEventProcessor_dc_alignment::PlotLines(deque<trajectory_t>&traj){
  unsigned int last_index=traj.size()-1;

  TCanvas *c1=dynamic_cast<TCanvas *>(gROOT->FindObject("endviewA Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
    
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].S(state_x),traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].S(state_x),traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }

  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("endviewA Large Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].S(state_x),traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].S(state_x),traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }
  
  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("sideviewA Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].z,traj[last_index].S(state_x));
    line->SetNextPoint(traj[0].z,traj[0].S(state_x));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }

  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("sideviewB Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].z,traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].z,traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    delete line;
  }
  // end of drawing code
  
}
