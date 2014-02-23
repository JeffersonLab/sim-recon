// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <stdio.h>
#include <unistd.h>

#include "MyProcessor.h"
#include "PID/DPhoton.h"
#include "PID/DKinFit.h"
#include "PID/DPid.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DMagneticFieldStepper.h"


int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;
float MAX_EVENTS = 1e20;
int COUNT = 0;
//TFile *fout = NULL;
//TH1F *hpi0[4];

vector<string> toprint;

#define ansi_escape		((char)0x1b)
#define ansi_bold 		ansi_escape<<"[1m"
#define ansi_black		ansi_escape<<"[30m"
#define ansi_red			ansi_escape<<"[31m"
#define ansi_green		ansi_escape<<"[32m"
#define ansi_blue			ansi_escape<<"[34m"
#define ansi_normal		ansi_escape<<"[0m"
#define ansi_up(A)		ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)		ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)	ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)		ansi_escape<<"["<<(A)<<"D"

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
  char name[256];
  ///< Called once at program start.
  fout = new TFile(OUTNAME, "RECREATE","Output file");
  cout << "Opened ROOT file " << OUTNAME <<endl;
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<4;j++)
    {
      for(int k=0;k<5;k++)
      {
        sprintf(name,"hmm%d_%d_%d",i,j,k);
        hmm[i][j][k] = new TH1F(name,name,100,-0.01,0.01);

        sprintf(name,"hcl%d_%d_%d",i,j,k);
        hcl[i][j][k] = new TH1F(name,name,100,0.0,1.0);

        for(int p=0;p<4;p++)
        {
          sprintf(name,"hinvmass%d_%d_%d_%d",i,j,k,p);
          hinvmass[i][j][k][p] = new TH1F(name,name,50,0.0,1.0);
        }

        for(int p=0;p<24;p++)
        {
          sprintf(name,"hpulls%d_%d_%d_%d",i,j,k,p);
          hpulls[i][j][k][p] = new TH1F(name,name,50,-5.0,5.0);
        }
      }
    }
  }

  // Create Magnetic Field Map
  new JParameterManager(); // normally, this is created by JApplication!
  bfield = new DMagneticFieldMapGlueX();

  // Get the field map
  const DMagneticFieldMapGlueX::DBfieldPoint_t* Bmap;
  int Npoints;
  bfield->GetTable(Bmap, Npoints);

  return NOERROR;
}

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{
  // Read in Magnetic field map
  JApplication* japp = eventLoop->GetJApplication();

  return NOERROR;
}

//------------------------------------------------------------------
// evnt
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{

  if(VERBOSE) cerr << "First thing in evnt....." << endl;

  vector<const DPhoton*> photons;
  vector<const DTrack*> tracks;
  vector<DKinematicData> kd_initialState;
  vector<DKinematicData> kd_finalState;
  vector<DKinematicData> kd_initialState_post;
  vector<DKinematicData> kd_finalState_post;
  vector<const DMCThrown*> mcthrown;
  DLorentzVector beam, targ, p, pip, pim;
  DLorentzVector beamfit, targfit, pfit, pipfit, pimfit;
  double mm2, mm2fit;
  double cl;
  double pulls[32];

  const DVector3 tofPlane_origin(0.0, 0.0, 616.0);
  const DVector3 tofPlane_normal(0.0, 0.0, 1.0);

  DMagneticFieldStepper dmagstep(bfield);
  dmagstep.SetStepSize(1.00);

  DKinematicData kd_beam = DKinematicData();
  DKinematicData kd_targ = DKinematicData();
  // Set the beam
  kd_beam.setMass(0.0);
  kd_beam.setCharge(0);
  kd_beam.setMassFixed();
  kd_beam.setMomentum( TVector3(0.0, 0.0, 9.0) );
  kd_beam.clearErrorMatrix();
  kd_initialState.push_back(kd_beam);
  // Set the target
  kd_targ.setMass(0.93827);
  kd_targ.setCharge(1);
  kd_targ.setMassFixed();
  kd_targ.setMomentum( TVector3(0.0, 0.0, 0.0) );
  kd_targ.clearErrorMatrix();
  kd_initialState.push_back(kd_targ);

  DPid *dpid = new DPid();
  dpid->SetVerbose(VERBOSE);

  float numclGood = 0;
  float clCut[5] = {-1.0, 0.01, 0.1, 0.2, 0.5};

  if(eventLoop->Get(tracks,"THROWN"))
  { 
    kd_finalState.clear();
    for(int i=0;i<(int)tracks.size();i++)
    {
      kd_finalState.push_back(*(tracks[i]));
      if(COUNT%100==0) cerr << COUNT << " " << MAX_EVENTS << endl;
      if(COUNT>=MAX_EVENTS) { fini(); exit(1); }
    }
    kd_finalState[0].setMass(0.13957);
    kd_finalState[1].setMass(0.13957);
    kd_finalState[2].setMass(0.93827);

    kd_finalState[0].setMassFixed();
    kd_finalState[1].setMassFixed();
    kd_finalState[2].setMassFixed();

    eventLoop->Get(mcthrown);

    beam = kd_initialState[0].lorentzMomentum();
    targ = kd_initialState[1].lorentzMomentum();
    pip = kd_finalState[0].lorentzMomentum();
    pim = kd_finalState[1].lorentzMomentum();
    p = kd_finalState[2].lorentzMomentum();

    mm2 = (beam + targ - p - pip - pim).M2();

    // Set up and run the fit
    dpid->SetFinal(kd_finalState);
    dpid->SetInitial(kd_initialState);
    //kfit->SetMissingParticle(0.134);
    //kfit->Fit();


    if(VERBOSE)
    {
      //cerr << "prob after fit: " << kfit->Prob() << endl;
      //cerr << "Chi2 after fit: " << kfit->Chi2() << endl;
    }


    COUNT++;
  }

  return NOERROR;
}

//------------------------------------------------------------------
// fini
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
  fout->cd();
  fout->Write();
  fout->cd();
  fout->Close();
  //  delete fout;
  cout<<endl<<"Closed ROOT file"<<endl;

  ///< Called after last event of last event source has been processed.
  return NOERROR;
}				
