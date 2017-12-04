#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <TRandom2.h>
#include <TRandom3.h>
#include "TTCS_crs.hh"
#include "TTCS_kine.hh"
//#include "kin_funcs.cc"
#include "kin_funcs.h"
#include <TLorentzVector.h>

// #include "CobremsGenerator.hh"

#define USEHDDM
#ifdef USEHDDM 
#include "HddmOut.h"
#endif

#include <TCanvas.h>

using namespace std;

int main(int argc, char **argv)
{
  TCanvas *c1 = new TCanvas();

  // output configuration
  int run = 10000;
  int Nsim = 25000;
  int seedVal = 0;

  //COMMAND LINE PARSING
  char *argptr;
  for (int i=1; i<argc; i++) {
      argptr = argv[i];
      if (*argptr == '-') {
          argptr++;
          switch (*argptr) {
          case 'N':
              run = atoi(++argptr);
              break;
          case 'n':
              Nsim = atoi(++argptr);
              break;
          case 'r':
              seedVal= atoi(++argptr);
              break;
          default:
              fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
              break;
          }
      }
  }

  
  cout << " run = " << run << endl;

  const double PI = 3.14159265358979312;
  const double radian = 57.2957795130823229;
  const double Mp = 0.9383;
  const double Me = 0.00051;
  const double t_lim = -1.5; // limit of t distribution Max(|t|)

  //GET THE HISTOGRAM FOR COHERENT BREMSTRAHLUNG SPECTRUM
  const double eLower = 8.;
  const double eUpper = 11.8;

  TH1D* hGvsE;
  string photonSpectrumFile="cobrems.root";
  TFile *inCoBrem=new TFile(photonSpectrumFile.c_str()); //Using spectrum from Richard
  hGvsE=(TH1D*)inCoBrem->Get("cobrem_vs_E");
  TH1D* hGvsEout = (TH1D*)hGvsE->Clone("hGvsEout");
  hGvsEout->Reset();
  hGvsEout->Rebin(30);
  int eBinLow = hGvsE->GetXaxis()->FindBin(eLower);
  int eBinHigh = hGvsE->GetXaxis()->FindBin(eUpper);
  hGvsE->GetXaxis()->SetRange(eBinLow,eBinHigh);
  double gMax = hGvsE->GetMaximum();

  //const double Eb = 12.;
  //const double Eg = Eb;

  
  
  //const double Q2min = 2*Mp*Eg + t_lim - (Eg/Mp)*( 2*Mp*Mp - t_lim - sqrt(t_lim*t_lim - 4*Mp*Mp*t_lim) );
  const double Q2min = 2*Me;
  const double Minv_min = sqrt(Q2min);

  TRandom3 rand;
  rand.SetSeed(seedVal);    // need to set this

  //TTCS_kine tcs_kin1(Mp, Eb);
  TTCS_crs crs_lmlp;

  TLorentzVector target(0., 0., 0., Mp);
  TLorentzVector Lcm;
  TLorentzVector Lbeam;

  TFile *file_out = new TFile("tcs_gen.root", "Recreate");

  TH2D *h_ph_h_ph_cm1 = new TH2D("h_ph_h_ph_cm1", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_th_g_th_cm1 = new TH2D("h_th_g_th_cm1", "", 200, 0., 180., 200, 0., 180.);
  TH1D *h_mee = new TH1D("mee", "", 200, 0, 4.);

  //================= Definition of Tree Variables =================
  double Minv, t, Q2;
  double psf, crs_BH, crs_INT, crs_int;
  double psf_flux, flux_factor;
  TLorentzVector L_em, L_ep, L_prot;
  TLorentzVector L_gprime;
  
  TTree *tr1 = new TTree("tr1", "TCS MC events");
  tr1->Branch("L_em", "TLorentzVector", &L_em, 3200, 99);
  tr1->Branch("L_ep", "TLorentzVector", &L_ep, 3200, 99);
  tr1->Branch("L_prot", "TLorentzVector", &L_prot, 3200, 99);
  tr1->Branch("Q2", &Q2, "Q2/D");
  tr1->Branch("t", &t, "t/D");
  tr1->Branch("psf", &psf, "psf/D");
  tr1->Branch("crs_BH", &crs_BH, "crs_BH/D");

#ifdef USEHDDM 
  //============Initialize HDDM output ===============
  char hddmFile[80];
  sprintf(hddmFile,"genOut.hddm");

  HddmOut hddmGo(hddmFile);
  int evtNumber = 0;
#endif 

  ofstream outf("ee.ascii");

  for( int i = 0; i < Nsim; i++ ) {

      if( i%50000 == 0)
        {
          cout.flush()<<"Processed "<<i<<" events, approximetely "<<double(100.*i/double(Nsim))<<"%\r";
        }
      
      // throw beam photon
      double eGamma = gMax;
      double yMax = gMax*1.02;
      double yVal = 0.0;
      double testValY = yMax + 10.0;
      int eBin;
      while(testValY>yVal){//Monte Carlo the event to get brem spectrum
          eGamma = rand.Uniform(eLower,eUpper); //Grab a photon energy
          testValY = rand.Uniform(0.0,yMax); //Grab a test value                  
          eBin = hGvsE->GetXaxis()->FindBin(eGamma);
          yVal = hGvsE->GetBinContent(eBin);
      }
      hGvsEout->Fill(eGamma);

      double Eb = eGamma;
      double Eg = Eb;
      TTCS_kine tcs_kin1(Mp, Eb);

      double s = Mp*Mp + 2*Mp*Eg;
      double t_min = T_min(0., Mp*Mp, Q2min, Mp*Mp, s);
      double t_max = T_max(0., Mp*Mp, Q2min, Mp*Mp, s);
      double psf_t = t_min - TMath::Max(t_max, t_lim);

      if( t_min > t_lim )
	{
	  t = rand.Uniform( t_min - psf_t,  t_min);
	  double Q2max = 2*Mp*Eg + t - (Eg/Mp)*( 2*Mp*Mp - t - sqrt(t*t - 4*Mp*Mp*t) ); // Page 182 of my notebook. Derived using "Q2max = s + t - 2Mp**2 + u_max" relation
	  
	  double psf_Q2 = Q2max - Q2min;

	  Q2 = rand.Uniform(Q2min, Q2min + psf_Q2);
	 
	  double u = 2*Mp*Mp + Q2 - s - t;
	  double th_qprime = acos((s*(t - u) - Mp*Mp*(Q2 - Mp*Mp))/sqrt(LambdaFunc(s, 0, Mp*Mp)*LambdaFunc(s, Q2, Mp*Mp))); //Byukling Kayanti (4.9)
	  double th_pprime = PI + th_qprime;
          
	  double Pprime = 0.5*sqrt(LambdaFunc(s, Q2, Mp*Mp)/s); // Momentum in c.m. it is the same for q_pr and p_pr
	  
	  Lbeam.SetPxPyPzE(0., 0., Eg, Eg);
	  Lcm.SetPxPyPzE(0., 0., Eg, Mp + Eg);
	  L_prot.SetPxPyPzE(Pprime*sin(th_pprime), 0., Pprime*cos(th_pprime), sqrt(Pprime*Pprime + Mp*Mp) );
	  L_gprime.SetPxPyPzE( Pprime*sin(th_qprime), 0., Pprime*cos(th_qprime), sqrt(Pprime*Pprime + Q2) );
	  
	  double psf_cos_th = 2.; // cos(th):(-1 : 1)
	  double psf_phi_cm = 2*PI;
	  
	  double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
	  double sin_th = sqrt( 1 - cos_th*cos_th );
	  double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);
	  
	  double El = sqrt(Q2)/2.; // Energy of lepton in the rest frame of qprime
	  double Pl = sqrt( El*El - Me*Me );
	  
	  L_em.SetPxPyPzE( Pl*sin_th*cos(phi_cm), Pl*sin_th*sin(phi_cm), Pl*cos_th, El );
	  L_ep.SetPxPyPzE( -Pl*sin_th*cos(phi_cm), -Pl*sin_th*sin(phi_cm), -Pl*cos_th, El );
	  
	  L_em.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
	  L_ep.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame

	  L_em.Boost(L_gprime.BoostVector()); // Move to the CM Frame
	  L_ep.Boost(L_gprime.BoostVector()); // Move to the CM Frame
	  
	  L_em.Boost(Lcm.BoostVector());  // Move to the Lab Frame
	  L_ep.Boost(Lcm.BoostVector());  // Move to the Lab Frame
	  
	  
	  L_gprime.Boost(Lcm.BoostVector());
	  L_prot.Boost(Lcm.BoostVector());	  
	  
	  double psf_phi_lab = 2*PI;
	  double phi_rot = rand.Uniform(0., psf_phi_lab);

	  L_prot.RotateZ(phi_rot);
	  L_gprime.RotateZ(phi_rot);
	  L_em.RotateZ(phi_rot);
	  L_ep.RotateZ(phi_rot);
	  tcs_kin1.SetLemLepLp(L_em, L_ep, L_prot);

	  psf = psf_t*psf_Q2*psf_phi_lab*psf_cos_th*psf_phi_cm;
	  
	  crs_BH = crs_lmlp.Eval_BH(s, Q2, t, -1, tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm());  // -1: cros section is not weighted by L/L0
	  
	  tr1->Fill();


      // monitoring histograms
	  h_ph_h_ph_cm1->Fill(phi_cm*TMath::RadToDeg(), tcs_kin1.GetPhi_cm());
	  h_th_g_th_cm1->Fill(acos(cos_th)*TMath::RadToDeg(), tcs_kin1.GetTheta_cm());
      if( (L_em.Theta()*180./TMath::Pi() > 1.) && (L_ep.Theta()*180./TMath::Pi() > 1.) )
          h_mee->Fill((L_em+L_ep).M(), psf*crs_BH);
      
#ifdef USEHDDM 
      // ======= HDDM output =========
      tmpEvt_t tmpEvt;
      tmpEvt.beam = Lbeam;
      tmpEvt.target = target;
      tmpEvt.q1 = L_em;
      tmpEvt.q2 = L_ep;
      tmpEvt.recoil = L_prot;
      tmpEvt.nGen = 3;
      tmpEvt.rxn = 2;
      //tmpEvt.weight = psf*crs_BH;
      //tmpEvt.weight = crs_BH;
      tmpEvt.weight = psf;
      evtNumber++;
      hddmGo.write(tmpEvt,run,evtNumber);
#endif

	}
      else
	{
	  cout<<" |t_min| > |t_lim|"<<endl;
	  cout<<" t_min =  "<<t_min<<"   t_lim = "<<t_lim<<"  Eg = "<<Eg<<"   Q2min = "<<Q2min<<endl;
	}
      
    }

  outf.close();
 
  tr1->Write();
  h_ph_h_ph_cm1->Write();
  h_th_g_th_cm1->Write();
  h_mee->Write();
  hGvsEout->Write(); 
  file_out->Close();
}
