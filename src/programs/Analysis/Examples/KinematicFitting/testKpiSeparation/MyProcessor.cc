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
        hinvmass[i][j][p] = new TH1F(name,name,160,0.2,1.8);

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
  DLorentzVector meson01[15], meson23[15];
  DLorentzVector cmFrame[15];
  DLorentzVector gjFrame[15];
  double mm2, mm2fit;
  double cl[15];
  double ndf[15];
  double pulls[32];
  float mm01[15], mm23[15];

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

  float PIMASS = 0.140;
  float KMASS = 0.494;
  float PROTONMASS = 0.938;

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
  k=1; massHypIndex[k][0]=0; massHypIndex[k][1]=1; massHypIndex[k][2]=4; massHypIndex[k][3]=3; massHypIndex[k][4]=2;
  k=2; massHypIndex[k][0]=4; massHypIndex[k][1]=1; massHypIndex[k][2]=2; massHypIndex[k][3]=3; massHypIndex[k][4]=0;
  //  K+=0  K-=1  pi+=2  pi-=3  p=4  
  k=3; massHypIndex[k][0]=0; massHypIndex[k][1]=1; massHypIndex[k][2]=4; massHypIndex[k][3]=3; massHypIndex[k][4]=2;
  k=4; massHypIndex[k][0]=0; massHypIndex[k][1]=1; massHypIndex[k][2]=2; massHypIndex[k][3]=3; massHypIndex[k][4]=4;
  k=5; massHypIndex[k][0]=2; massHypIndex[k][1]=1; massHypIndex[k][2]=0; massHypIndex[k][3]=3; massHypIndex[k][4]=4;
  k=6; massHypIndex[k][0]=2; massHypIndex[k][1]=1; massHypIndex[k][2]=4; massHypIndex[k][3]=3; massHypIndex[k][4]=0;
  k=7; massHypIndex[k][0]=4; massHypIndex[k][1]=1; massHypIndex[k][2]=2; massHypIndex[k][3]=3; massHypIndex[k][4]=0;
  k=8; massHypIndex[k][0]=4; massHypIndex[k][1]=1; massHypIndex[k][2]=0; massHypIndex[k][3]=3; massHypIndex[k][4]=2;

  k=9; massHypIndex[k][0]=0; massHypIndex[k][1]=3; massHypIndex[k][2]=4; massHypIndex[k][3]=1; massHypIndex[k][4]=2;
  k=10;massHypIndex[k][0]=0; massHypIndex[k][1]=3; massHypIndex[k][2]=2; massHypIndex[k][3]=1; massHypIndex[k][4]=4;
  k=11;massHypIndex[k][0]=2; massHypIndex[k][1]=3; massHypIndex[k][2]=0; massHypIndex[k][3]=1; massHypIndex[k][4]=4;
  k=12;massHypIndex[k][0]=2; massHypIndex[k][1]=3; massHypIndex[k][2]=4; massHypIndex[k][3]=1; massHypIndex[k][4]=0;
  k=13;massHypIndex[k][0]=4; massHypIndex[k][1]=3; massHypIndex[k][2]=2; massHypIndex[k][3]=1; massHypIndex[k][4]=0;
  k=14;massHypIndex[k][0]=4; massHypIndex[k][1]=3; massHypIndex[k][2]=0; massHypIndex[k][3]=1; massHypIndex[k][4]=2;

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
      pos = tracks[i]->position();
      mom = tracks[i]->momentum();
      charge = tracks[i]->charge();

      dmagstep.SetStartingParams(charge, dum3vec, dum3vec);
      if(VERBOSE) cerr << "Original pos: " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
      if(VERBOSE) cerr << "Original mom: " << mom.X() << " " << mom.Y() << " " << mom.Z() << endl;
      if(VERBOSE) cerr << "Original mom: " << mom.Mag() << " " << mom.CosTheta() << " " << mom.Phi() << endl;
      dmagstep.SwimToRadius(pos, mom, bcal_radius, &path_length);
      if(VERBOSE) cerr << "Orig BCAL Radius: " << path_length << " " << pos.Perp() << " " << pos.Z() << endl;

      hit_BCAL[i] = false;
      double z = pos.Z();
      double true_pathLength = -999.0;
      double true_t0 = -999.0;
      double true_t1 = -999.0;
      DetectorSystem_t detector = DetectorSystem_t(SYS_NULL);
      if( z>bcal_upstream_z && z<bcal_downstream_z )
      {
        hit_BCAL[i] = true;
        HIT_DETECTOR[0][i]++;
      }

      if(hit_BCAL[i])
      {
        hdetectorhits[0]->Fill(pos.X(), pos.Y(), pos.Z() );
        true_pathLength = path_length;
        detector = SYS_BCAL;
      }

      hit_TOF[i] = false;
      if ( !hit_BCAL[i] )
      {
        pos = tracks[i]->position();
        mom = tracks[i]->momentum();

        dmagstep.SetStartingParams(charge, dum3vec, dum3vec);
        dmagstep.SwimToPlane(pos, mom, tofPlane_origin, tofPlane_normal, &path_length);
        double x = pos.X();
        double y = pos.Y();

        if ( !(x<tof_x[0] && x>tof_x[1] && y<tof_y[0] && y>tof_y[1]) && // outside the TOF
             !(x>tof_inner_x[0] && x<tof_inner_x[1] && y>tof_inner_y[0] && y<tof_inner_y[1]) ) // Hole in center
        {
          hit_TOF[i] = true;
          HIT_DETECTOR[1][i]++;
        }

        if(hit_TOF[i])
        {
          hdetectorhits[1]->Fill(pos.X(), pos.Y(), pos.Z() );
          true_pathLength = path_length;
          //true_pathLength = 500;
          detector = SYS_TOF;
        }

      }
      if( hit_BCAL[i] || hit_TOF[i])
      {
        HIT_DETECTOR[2][i]++;
        double pl_res = 1.0; // 1cm for path length
        double t0_res = 0.050; // 50 ps
        double t1_res = 0.080; // 80 ps TOF
        if (hit_BCAL[i]) t1_res = 0.420; // 420 ps BCAL
        if(VERBOSE) cerr << "True beta: " << mcthrown[i]->lorentzMomentum().Beta() << endl;
        true_t1 = true_pathLength/(SPEED_OF_LIGHT*mcthrown[i]->lorentzMomentum().Beta());
        double smear_t0 = rnd.Gaus(0.0, t0_res);
        double smear_t1 = rnd.Gaus(0.0, t1_res);
        double smear_pathLength = rnd.Gaus(0.0, pl_res);
        double meas_t0 = 0.0 + smear_t0;
        double meas_t1 = true_t1 + smear_t1;
        double meas_pathLength = true_pathLength + smear_pathLength;
        kd_finalState[i].setT0(meas_t0, t0_res, SYS_TAGGER);
        kd_finalState[i].setT1(meas_t1, t1_res, detector);
        kd_finalState[i].setPathLength(meas_pathLength, pl_res);
        if(VERBOSE) cerr << "final t0: " << i << " " << kd_finalState[i].t0() << " " << kd_finalState[i].t0_err() << endl;
        if(VERBOSE) cerr << "final t1: " << i << " " << kd_finalState[i].t1() << " " << kd_finalState[i].t1_err() << endl;
        if(VERBOSE) cerr << "final pl: " << i << " " << kd_finalState[i].pathLength() << " " << kd_finalState[i].pathLength_err() << endl;
        if(VERBOSE) cerr << "final beta: " << i << " " << kd_finalState[i].measuredBeta() << " " << kd_finalState[i].lorentzMomentum().Beta() << endl;
        if(VERBOSE) cerr << "final lm: " << i << " " << kd_finalState[i].lorentzMomentum().E() << " " <<  kd_finalState[i].lorentzMomentum().Rho() << " " << kd_finalState[i].lorentzMomentum().E()/kd_finalState[i].lorentzMomentum().Rho() << endl;
      }
      if(VERBOSE) cerr << "BCAL/TOF: " << hit_BCAL[i] << " " << hit_TOF[i] << endl;
    }

    for(int i=0;i<15;i++)
    //for(int i=4;i<5;i++)
    {

      if(i<3) {massHyp[0]=PIMASS;massHyp[1]=PIMASS;massHyp[2]=PIMASS;massHyp[3]=PIMASS;massHyp[4]=PROTONMASS;}
      else    {massHyp[0]=KMASS;massHyp[1]=KMASS;massHyp[2]=PIMASS;massHyp[3]=PIMASS;massHyp[4]=PROTONMASS;}

      for(int j=0;j<5;j++)
      {
        kd_finalState[ massHypIndex[i][j] ].setMass(massHyp[j]);
        kd_finalState[j].setMassFixed();
      }

      // Set up and run the fit
      kfit->ResetForNewFit();
      if(USE_TIMING) kfit->UseTimingInformation();
      kfit->SetFinal(kd_finalState);
      kfit->SetInitial(kd_initialState);
      kfit->Fit();

      cl[i]=kfit->Prob();
      ndf[i]=kfit->Ndf();

      kd_finalState_post = kfit->GetFinal_out();

      meson[0] = kd_finalState_post[ massHypIndex[i][0] ].lorentzMomentum(); // K+ or pi+
      meson[1] = kd_finalState_post[ massHypIndex[i][1] ].lorentzMomentum(); // K+ or pi-
      meson[2] = kd_finalState_post[ massHypIndex[i][2] ].lorentzMomentum(); // pi+
      meson[3] = kd_finalState_post[ massHypIndex[i][3] ].lorentzMomentum(); // pi-
      p        = kd_finalState_post[ massHypIndex[i][4] ].lorentzMomentum(); // pi-

      meson01[i] = meson[0]+meson[1];
      meson23[i] = meson[2]+meson[3];
      mm01[i] = meson01[i].M();
      mm23[i] = meson23[i].M();

      // (orgVector, GJ system, recoil, target)
      cmFrame[i]=meson01[i];
      cmFrame[i].Boost(-(kd_targ.lorentzMomentum()+kd_beam.lorentzMomentum()).BoostVector());
      gjFrame[i]=GJFrame_TLV(meson[0], (p + meson23[i]), meson01[i], kd_targ.lorentzMomentum());
      //

      for(int k=0;k<5;k++)
      {
        if(cl[i]>clCut[k] && ndf[i]>0 )
        {
          hcl[i][k]->Fill(cl[i]);
          hinvmass[i][k][0]->Fill(mm01[i]);
          hinvmass[i][k][1]->Fill(mm23[i]);
          hinvmass[i][k][2]->Fill(mm01[i]);
          hinvmass[i][k][2]->Fill(mm23[i]);
          hpvtheta[i][k][0]->Fill(180.0*meson01[i].Theta()/TMath::Pi(),meson01[i].Rho()); 
          hpvcosCM[i][k][0]->Fill(cmFrame[i].CosTheta(),cmFrame[i].Rho()); 
          hphivcosGJ[i][k][0]->Fill(gjFrame[i].CosTheta(),gjFrame[i].Phi()/TMath::Pi()); 
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

    int goodPIndex[5] = {0, 4, 5, 10, 11};
    int acI[4] = {4, 5, 10, 11}; // anti-cut Index
    for(int i=0;i<15;i++)
    {
      if(i==0)       {acI[0]=4;acI[1]=5;acI[2]=10;acI[3]=11;}
      else if(i==4)  {acI[0]=0;acI[1]=5;acI[2]=10;acI[3]=11;}
      else if(i==5)  {acI[0]=0;acI[1]=4;acI[2]=10;acI[3]=11;}
      else if(i==10) {acI[0]=0;acI[1]=4;acI[2]=5;acI[3]=11;}
      else if(i==11) {acI[0]=0;acI[1]=4;acI[2]=5;acI[3]=10;}
      else           {acI[0]=4;acI[1]=5;acI[2]=10;acI[3]=11;}
      for(int k=0;k<5;k++)
      {
        if(numpass[k]==1 && cl[i]>clCut[k])
        {
          hinvmass[i][k][4]->Fill(mm01[i]);
          hinvmass[i][k][5]->Fill(mm23[i]);
          hinvmass[i][k][6]->Fill(mm01[i]);
          hinvmass[i][k][6]->Fill(mm23[i]);
          hpvtheta[i][k][4]->Fill(180.0*meson01[i].Theta()/TMath::Pi(),meson01[i].Rho()); 
          hpvcosCM[i][k][4]->Fill(cmFrame[i].CosTheta(),cmFrame[i].Rho()); 
          hphivcosGJ[i][k][4]->Fill(gjFrame[i].CosTheta(),gjFrame[i].Phi()/TMath::Pi()); 
        }
        if(numfail[1]==14 && numpass[k]==1 && cl[i]>clCut[k]) // Anti cut on 14 other hyp at 1% CL
        {
          hinvmass[i][k][8]->Fill(mm01[i]);
          hinvmass[i][k][9]->Fill(mm23[i]);
          hinvmass[i][k][10]->Fill(mm01[i]);
          hinvmass[i][k][10]->Fill(mm23[i]);
          hpvtheta[i][k][8]->Fill(180.0*meson01[i].Theta()/TMath::Pi(),meson01[i].Rho()); 
          hpvcosCM[i][k][8]->Fill(cmFrame[i].CosTheta(),cmFrame[i].Rho()); 
          hphivcosGJ[i][k][8]->Fill(gjFrame[i].CosTheta(),gjFrame[i].Phi()/TMath::Pi()); 
        }
        if(numpass[k]==1 && cl[i]>clCut[k] &&
            cl[acI[0]]<0.01 &&
            cl[acI[1]]<0.01 &&
            cl[acI[2]]<0.01 &&
            cl[acI[3]]<0.01 ) // Anti cut on only other proton hyp at 1% CL
        {
          hinvmass[i][k][12]->Fill(mm01[i]);
          hinvmass[i][k][13]->Fill(mm23[i]);
          hinvmass[i][k][14]->Fill(mm01[i]);
          hinvmass[i][k][14]->Fill(mm23[i]);
          hpvtheta[i][k][12]->Fill(180.0*meson01[i].Theta()/TMath::Pi(),meson01[i].Rho()); 
          hpvcosCM[i][k][12]->Fill(cmFrame[i].CosTheta(),cmFrame[i].Rho()); 
          hphivcosGJ[i][k][12]->Fill(gjFrame[i].CosTheta(),gjFrame[i].Phi()/TMath::Pi()); 
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
