// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <stdio.h>
#include <unistd.h>

#include "frames_TLV.h"

#include "MyProcessor.h"
#include "PID/DPhoton.h"
#include "PID/DKinFit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "GlueX.h"


int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;
float MAX_EVENTS = 1e20;
int COUNT = 0;
int HIT_DETECTOR[3][5];
DRandom rnd;
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

#define SPEED_OF_LIGHT 29.9792

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
  char name[256];
  ///< Called once at program start.
  fout = new TFile(OUTNAME, "RECREATE","Output file");
  cout << "Opened ROOT file " << OUTNAME <<endl;
  for(int i=0;i<15;i++)
  {
    sprintf(name,"hdetectorhits%d",i);
    hdetectorhits[i] = new TH3F(name,name,65,-130.0, 130.0, 65,-130.0, 130.0, 50, 0.0, 650);

    for(int j=0;j<5;j++)
    {
      sprintf(name,"hcl%d_%d",i,j);
      hcl[i][j] = new TH1F(name,name,100,0.0,1.0);

      sprintf(name,"hnumpass%d_%d",i,j);
      hnumpass[i][j] = new TH1F(name, name, 15, -0.5, 14.5);

      for(int p=0;p<16;p++)
      {
        sprintf(name,"hinvmass%d_%d_%d",i,j,p);
        if(p==0) hinvmass[i][j][p] = new TH1F(name,name,100,0.0,0.4);
        else     hinvmass[i][j][p] = new TH1F(name,name,100,0.5,1.0);

        sprintf(name,"hpvtheta%d_%d_%d",i,j,p);
        hpvtheta[i][j][p] = new TH2F(name,name,80,0.0,40.0, 80,0.0,8.0);

        sprintf(name,"hpvcosCM%d_%d_%d",i,j,p);
        hpvcosCM[i][j][p] = new TH2F(name,name,100,-1.0,1.0, 100,0.0,2.0);

        sprintf(name,"hphivcosGJ%d_%d_%d",i,j,p);
        hphivcosGJ[i][j][p] = new TH2F(name,name,100,-1.0,1.0, 100,-1.0,1.0);

      }

    }
  }

  rnd.SetSeed((unsigned)time(NULL));

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
  DLorentzVector meson[4];
  DLorentzVector meson01[15], meson23[15], meson0123[15];
  DLorentzVector cmFrame[15];
  DLorentzVector gjFrame[15];
  double mm2, mm2fit;
  double cl[15];
  double ndf[15];
  double pulls[32];
  float mm01[15], mm23[15], mm0123[15];

  const DVector3 tofPlane_origin(0.0, 0.0, 619.0);
  const DVector3 tofPlane_normal(0.0, 0.0, 1.0);
  double bcal_upstream_z = 17;
  double bcal_downstream_z = 407;
  bool hit_BCAL[5];
  bool hit_TOF[5];

  for(int j=0;j<3;j++)
    for(int i=0;i<5;i++)
      HIT_DETECTOR[j][i] = 0;

  double tof_x[2] = {-120.0, 120.0};
  double tof_y[2] = {-120.0, 120.0};
  double tof_inner_x[2] = {-6.0, 6.0};
  double tof_inner_y[2] = {-6.0, 6.0};
  double bcal_radius = 65.0;
  double path_length;
  DVector3 pos, mom;
  const DVector3 *dum3vec = new TVector3(1.0, 1.0, 1.0);
  double charge;

  DVector3 testPos(0.0, 0.0, 0.0);

  DMagneticFieldStepper dmagstep(bfield);
  dmagstep.SetStepSize(1.00);

  float PIMASS = 0.13957;
  float PI0MASS = 0.13498;
  float KMASS = 0.494;
  float PROTONMASS = 0.93827;

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
  kd_targ.setMass(PROTONMASS);
  kd_targ.setCharge(1);
  kd_targ.setMassFixed();
  kd_targ.setMomentum( TVector3(0.0, 0.0, 0.0) );
  kd_targ.clearErrorMatrix();
  kd_initialState.push_back(kd_targ);

  DKinFit *kfit = new DKinFit();
  kfit->SetVerbose(VERBOSE);

  float numclGood = 0;
  float clCut[5] = {-1.0, 0.01, 0.1, 0.2, 0.5};

  int k;
  int massHypIndex[18][5];
  //  pi+=0  pi-=1  pi+=2  pi-=3  p=4  
  k=0; massHypIndex[k][0]=0; massHypIndex[k][1]=1; massHypIndex[k][2]=2; massHypIndex[k][3]=3; massHypIndex[k][4]=4;

  float massHyp[5];

  int numpass[5] = {0,0,0,0,0};
  int numfail[5] = {0,0,0,0,0};

  // Get the MC truth stuff.
  eventLoop->Get(mcthrown);

  if(eventLoop->Get(tracks,"THROWN"))
  { 
    if(COUNT%100==0) cerr << COUNT << " " << MAX_EVENTS << endl;
    if(COUNT>=MAX_EVENTS) { fini(); exit(1); }

    kd_finalState.clear();
    for(int i=0;i<(int)tracks.size();i++)
    {
      kd_finalState.push_back(*(tracks[i]));
    }

    int i = 0;
    if(eventLoop->Get(photons))
    { 
      int numPhotons = photons.size();
      if(numPhotons == 2)
      {
        kd_finalState[2] = (*photons[0]);
        kd_finalState[3] = (*photons[1]);

          if(1) {massHyp[0]=PIMASS;massHyp[1]=PIMASS;massHyp[2]=0.0;massHyp[3]=0.0;massHyp[4]=PROTONMASS;}

          for(int j=0;j<5;j++)
          {
            kd_finalState[ massHypIndex[i][j] ].setMass(massHyp[j]);
            kd_finalState[j].setMassFixed();
          }

          // Set up and run the fit
          kfit->ResetForNewFit();
          kfit->SetFinal(kd_finalState);
          kfit->SetInitial(kd_initialState);
          std::vector<int> constraintParticles;
          constraintParticles.push_back(2); // photon
          constraintParticles.push_back(3); // photon
          kfit->SetMassConstraint(PI0MASS,constraintParticles);
          kfit->Fit();

          cl[i]=kfit->Prob();
          ndf[i]=kfit->Ndf();

          kd_finalState_post = kfit->GetFinal_out();

          meson[0] = kd_finalState_post[ massHypIndex[i][0] ].lorentzMomentum(); // pi+
          meson[1] = kd_finalState_post[ massHypIndex[i][1] ].lorentzMomentum(); // pi-
          meson[2] = kd_finalState_post[ massHypIndex[i][2] ].lorentzMomentum(); // gamma
          meson[3] = kd_finalState_post[ massHypIndex[i][3] ].lorentzMomentum(); // gamma
          p        = kd_finalState_post[ massHypIndex[i][4] ].lorentzMomentum(); // pi-

          meson01[i] = meson[0]+meson[1];
          meson23[i] = meson[2]+meson[3];
          meson0123[i] = meson01[i]+meson23[i];
          mm01[i] = meson01[i].M();
          mm23[i] = meson23[i].M();
          mm0123[i] = meson0123[i].M();

          for(int k=0;k<5;k++)
          {
            if(cl[i]>clCut[k] && ndf[i]>0 )
            {
              hcl[i][k]->Fill(cl[i]);
              hinvmass[i][k][0]->Fill(mm23[i]);
              hinvmass[i][k][1]->Fill(mm01[i]);
              hinvmass[i][k][2]->Fill(mm0123[i]);
              numpass[k]++;
              hnumpass[1][k]->Fill(i);
            }
            if(cl[i]<clCut[k])
            {
              numfail[k]++;
            }
          }

          if(VERBOSE)
          {
            cerr << "prob after fit: " << kfit->Prob() << endl;
            cerr << "Chi2 after fit: " << kfit->Chi2() << endl;
          }

        }

    }


    for(int k=0;k<5;k++)
    {
      hnumpass[0][k]->Fill(numpass[k]);
    }

    COUNT++;
  }

  delete kfit;
  delete dum3vec;


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
