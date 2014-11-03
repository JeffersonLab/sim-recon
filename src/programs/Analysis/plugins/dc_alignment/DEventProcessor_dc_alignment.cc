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
  READ_CDC_FILE=false;
  gPARMS->SetDefaultParameter("DCALIGN:READ_CDC_FILE",READ_CDC_FILE);
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

   if (READ_CDC_FILE){
    ifstream cdcfile("cdc_alignment.dat");  
    for (unsigned int ring=0;ring<cdc_alignments.size();ring++){
      for (unsigned int straw=0;straw<cdc_alignments[ring].size();
	   straw++){
	double dxu,dyu,dxd,dyd;

	cdcfile >> dxu;
	cdcfile >> dyu;
	cdcfile >> dxd;
	cdcfile >> dyd;

	cdc_alignments[ring][straw].A(k_dXu)=dxu;
	cdc_alignments[ring][straw].A(k_dYu)=dyu;
	cdc_alignments[ring][straw].A(k_dXd)=dxd;
	cdc_alignments[ring][straw].A(k_dYd)=dyd;
      }
    }
    cdcfile.close();
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
  if (jcalib->Get("CDC/cdc_drift_table::NoBField", tvals)==false){    
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

   // Get pointer to TrackFinder object 
  vector<const DTrackFinder *> finders;
  loop->Get(finders);

  if(finders.size()<1){
    _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }

  // Drop the const qualifier from the DTrackFinder pointer
  finder = const_cast<DTrackFinder*>(finders[0]);

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
    //cdcfile << "Ring straw dXu dYu dXd dYd" << endl;
    for (unsigned int ring=0;ring<cdc_alignments.size();ring++){
      for (unsigned int straw=0;straw<cdc_alignments[ring].size();
	   straw++){
	//	cdcfile << ring+1 << " " << straw+1 << " " 
	cdcfile << cdc_alignments[ring][straw].A(k_dXu) << " " 
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

  // Reset the track finder
  finder->Reset();
  
  // Get BCAL showers, FCAL showers and FDC space points
  vector<const DFCALShower*>fcalshowers;
  if (USE_FCAL) loop->Get(fcalshowers);
  vector<const DBCALShower*>bcalshowers;
  if (USE_BCAL)loop->Get(bcalshowers);

  vector<const DFDCPseudo*>pseudos;
  loop->Get(pseudos);
  vector<const DCDCTrackHit*>cdcs;
  if (COSMICS) loop->Get(cdcs);

  if (cdcs.size()>20 /* && cdcs.size()<60*/){
    // Add the hits to the finder helper class, link axial hits into segments
    // then link axial hits and stereo hits together to form track candidates
    for (size_t i=0;i<cdcs.size();i++) finder->AddHit(cdcs[i]);
    finder->FindAxialSegments();
    finder->LinkCDCSegments();

    // Get the list of linked segments and fit the hits to lines
    const vector<DTrackFinder::cdc_track_t>tracks=finder->GetCDCTracks();
    for (unsigned int i=0;i<tracks.size();i++){
      // Add lists of stereo and axial hits associated with this track 
      // and sort
      vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
      hits.insert(hits.end(),tracks[i].stereo_hits.begin(),tracks[i].stereo_hits.end());
      sort(hits.begin(),hits.end(),cdc_hit_cmp);
      
      // Use earliest cdc time to estimate t0
      double t0=1e6;
      for (unsigned int j=0;j<hits.size();j++){
	double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
	double t_test=hits[j]->tdrift-L/29.98;
	if (t_test<t0) t0=t_test;
	}
      
      // Initial guess for state vector
      DMatrix4x1 S(tracks[i].S);
      
      // Run the Kalman Filter algorithm
      DoFilter(t0,tracks[i].z,S,hits);	 
    }
  }
  
  //-------------------------------------------------------------------------
  // FDC alignment 
  //-------------------------------------------------------------------------
  if (pseudos.size()>MIN_PSEUDOS
	//&&((fcalshowers.size()>0&&fcalshowers.size()<3)
	//   || (bcalshowers.size()>0&&bcalshowers.size()<3))
      ){
    // Add hits to the track finder helper class, link hits into segments 
    // then link segments together to form track candidates
    for (size_t i=0;i<pseudos.size();i++) finder->AddHit(pseudos[i]);
    finder->FindFDCSegments();
    finder->LinkFDCSegments();
    
    // Get the list of linked segments
    const vector<DTrackFinder::fdc_segment_t>tracks=finder->GetFDCTracks();
    
    // Loop over linked segments
    for (unsigned int k=0;k<tracks.size();k++){
      vector<const DFDCPseudo *>hits=tracks[k].hits;
      
      if (hits.size()>MIN_PSEUDOS){
	sort(hits.begin(),hits.end(),fdc_pseudo_cmp);
	
	// Initial guess for state vector    
	DMatrix4x1 S(tracks[k].S);
	
	// Move x and y to just before the first hit
	double my_z=hits[0]->wire->origin.z()-1.;
	S(state_x)+=my_z*S(state_tx);
	S(state_y)+=my_z*S(state_ty);
	
	// Use earliest fdc time to estimate t0
	double t0=1e6;
	double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
	for (unsigned int m=0;m<hits.size();m++){
	  if (hits[m]->time<t0){
	    double L=(hits[m]->wire->origin.z()-my_z)*dsdz;
	    t0=hits[m]->time-L/29.98; // assume moving at speed of light
	  }
	}
	
	// Run the Kalman Filter algorithm
	if (ALIGN_WIRE_PLANES) DoFilterAnodePlanes(t0,my_z,S,hits);
	else DoFilterCathodePlanes(t0,my_z,S,hits);
      }
    } //loop over tracks
  } // minimimum number of pseudopoints?
   
  return NOERROR;
}

// Steering routine for the kalman filter
jerror_t 
DEventProcessor_dc_alignment::DoFilter(double t0,double OuterZ,DMatrix4x1 &S,
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
    if (SetReferenceTrajectory(t0,OuterZ,S,trajectory,
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
      
      //printf("Timebased-----------\n");
      //if (false)
      for (iter=0;iter<20;iter++){
	chi2_old=chi2; 
	ndof_old=ndof;
	
	trajectory.clear();
	if (SetReferenceTrajectory(t0,OuterZ,S,trajectory,hits[maxindex])
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
	    
	  }  // check on final fit CL
      } // at least one time-based fit worked?
    } // check on preliminary fit CL
  } // at least one iteration worked?

  return NOERROR;
}


// Steering routine for the kalman filter
jerror_t 
DEventProcessor_dc_alignment::DoFilterCathodePlanes(double t0,double start_z,
						    DMatrix4x1 &S,
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
    if (SetReferenceTrajectory(t0,start_z,S,trajectory,hits)!=NOERROR) break;
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
DEventProcessor_dc_alignment::DoFilterAnodePlanes(double t0,double start_z,
						  DMatrix4x1 &S,
				       vector<const DFDCPseudo *>&hits){
  unsigned int num_hits=hits.size();
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
    if (SetReferenceTrajectory(t0,start_z,S,trajectory,hits)!=NOERROR) break;
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
	unsigned int layer=hits[i]->wire->layer;
     
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
					   vector<const DFDCPseudo *>&hits,
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
	double cosa=hits[id]->wire->udir.y();
	double sina=hits[id]->wire->udir.x();
	
	// State vector
	double x=Ss(state_x);
	double y=Ss(state_y);
	double tx=Ss(state_tx);
	double ty=Ss(state_ty);
	
	// Get the aligment vector and error matrix for this layer
	unsigned int layer=hits[id]->wire->layer-1;
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
	double uwire=hits[id]->wire->u+delta_u;
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
	double d=finder->FindDoca(trajectory[m].z,Ss,wdir,origin);
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
  const double d_EPS=1e-8;
 
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
      double d=diff.Mag()+d_EPS; // prevent division by zero

      // The next measurement and its variance
      double tdrift=hits[cdc_index]->tdrift-trajectory[k].t;
      double dmeas=0.39; 
      if (timebased){
	double drift_var=cdc_variance(tdrift);
	dmeas=cdc_drift_distance(tdrift);
	V=anneal_factor*drift_var;

	//printf("t %f d %f %f V %f\n",hits[cdc_index]->tdrift,dmeas,d,V);
      }

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
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

      double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
      double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
      double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
      double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);
      
      H(state_x)=H_T(state_x)
	=one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
		     +diffz*(dsdx-dtdx));
      H(state_y)=H_T(state_y)
	=one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
		     +diffz*(dsdy-dtdy));

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
	d=finder->FindDoca(trajectory[k].z,S,wdir,origin);
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
  DMatrix2x2 V(0.0008*anneal_factor,0.,0.,0.0008*anneal_factor);  // Measurement variance 
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
	double vv=-vpred*sinphi_v-uwire*cosphi_v+A(kV);
	double vu=-vpred*sinphi_u-uwire*cosphi_u+A(kU);

	// Difference between measurements and predictions
	Mdiff(0)=hits[my_id]->u-vu;
	Mdiff(1)=hits[my_id]->v-vv;

	// Start filling the update vector
	updates[my_id].drift_time=hits[my_id]->time-trajectory[k].t;
	
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
			       vector<const DFDCPseudo *>&hits,
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
	double uwire=hits[my_id]->wire->u+delta_u;

	// Find drift distance
	double drift_time=hits[my_id]->time-trajectory[k].t;
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
::SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
			 deque<trajectory_t>&trajectory,
			 const DCDCTrackHit *last_cdc){ 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double ds=1.0;
  double dz=(S(state_ty)>0.?-1.:1.)*ds/sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  double t=t0;

  //y-position after which we cut off the loop
  double min_y=last_cdc->wire->origin.y()-5.;
  unsigned int numsteps=0;
  do{
    z+=dz;
    J(state_x,state_tx)=-dz;
    J(state_y,state_ty)=-dz;
    // Flight time: assume particle is moving at the speed of light
    t+=ds/29.98;
    //propagate the state to the next z position
    S(state_x)+=S(state_tx)*dz;
    S(state_y)+=S(state_ty)*dz;
    trajectory.push_front(trajectory_t(z,t0,S,J,Zero4x1,Zero4x4));

    numsteps++;
  }while (S(state_y)>min_y && numsteps<MAX_STEPS);

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
::SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
			 deque<trajectory_t>&trajectory,
			 vector<const DFDCPseudo *>&pseudos){
  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double dz=1.1;
  double t=t0;
  trajectory.push_front(trajectory_t(z,t0,S,J,Zero4x1,Zero4x4));

  double zhit=z;
  double old_zhit=z;
  for (unsigned int i=0;i<pseudos.size();i++){  
    zhit=pseudos[i]->wire->origin.z();
    dz=1.1;
    
    if (fabs(zhit-old_zhit)<EPS){
      trajectory[0].num_hits++;
      continue;
    }
    bool done=false;
    while (!done){
      double new_z=z+dz;	      
	  
      if (new_z>zhit){
	dz=zhit-z;
	new_z=zhit;
	done=true;
      }
      J(state_x,state_tx)=-dz;
      J(state_y,state_ty)=-dz;
      // Flight time: assume particle is moving at the speed of light
      t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;
      //propagate the state to the next z position
      S(state_x)+=S(state_tx)*dz;
      S(state_y)+=S(state_ty)*dz;
      
      trajectory.push_front(trajectory_t(new_z,t,S,J,Zero4x1,Zero4x4)); 
      if (done){
	trajectory[0].h_id=i+1;
	trajectory[0].num_hits=1;
      }
      
      z=new_z;
    }	   
    old_zhit=zhit;
  }
  // One last step
  dz=1.1;
  J(state_x,state_tx)=-dz;
  J(state_y,state_ty)=-dz;

  // Flight time: assume particle is moving at the speed of light
  t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;

  //propagate the state to the next z position
  S(state_x)+=S(state_tx)*dz;
  S(state_y)+=S(state_ty)*dz;
  trajectory.push_front(trajectory_t(z+dz,t,S,J,Zero4x1,Zero4x4));

  if (false)
    {
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
	DMatrix4x1 A=cdc_alignments[ring][straw].A+dA;
	// Restrict offsets to less than 2 mm 
	if (fabs(A(k_dXu))<0.2 && fabs(A(k_dXd))<0.2 && fabs(A(k_dYu))<0.2 
	    && fabs(A(k_dYd))<0.2){
	  cdc_alignments[ring][straw].E=Etemp;
	  cdc_alignments[ring][straw].A=A;	  
	}
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
DEventProcessor_dc_alignment::FindOffsets(vector<const DFDCPseudo *>&hits,
				      vector<wire_update_t>&smoothed_updates){
  DMatrix1x2 G;//matrix relating alignment vector to measurement coords
  DMatrix2x1 G_T; // .. and its transpose
  
  unsigned int num_hits=hits.size();


  for (unsigned int i=0;i<num_hits;i++){ 
    double x=smoothed_updates[i].S(state_x);
    double y=smoothed_updates[i].S(state_y);
    double tx=smoothed_updates[i].S(state_tx);
    double ty=smoothed_updates[i].S(state_ty);
    
    double cosa=hits[i]->wire->udir.y();
    double sina=hits[i]->wire->udir.x();
      
    // Get the aligment vector and error matrix for this layer
    unsigned int layer=hits[i]->wire->layer-1;
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
    double uwire=hits[i]->wire->u+delta_u;
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
