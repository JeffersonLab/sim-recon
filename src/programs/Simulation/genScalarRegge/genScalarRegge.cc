// Main program for generating scalar events. 
#include "HDDM/hddm_s.h"
#include "particleType.h"

#include <TMath.h>
#include <TRandom3.h>
#include <TGenPhaseSpace.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
using namespace std;

extern "C"{
  void cobrems_(float *Ee,float *Epeak, float *emitmr, float *radt,
		float *dist,float *diameter,int *doPolFlux);
  float dntdx_(float* x);
  float dnidx_(float* x);
}

// Masses
const double m_p=0.93827; // GeV
const double m_p_sq=m_p*m_p;
double m_S=0.98; // GeV
double m_S_sq_R=m_S*m_S;
// Width
double width=0.;
// Coupling constants 
double gsq_rho_S_gamma=0.02537;
double gsq_omega_S_gamma=0.2283;

double Emin=3.,Emax=12.0; // GeV
double zmin=50.0,zmax=80.0; // cm, target extent
int Nevents=10000;
int runNo=9000;
bool debug=false;

// Diagnostic histograms
TH1D *thrown_t;
TH1D *thrown_dalitzZ;
TH1D *thrown_Egamma;
TH2D *thrown_dalitzXY;  
TH2D *thrown_theta_vs_p;
TH1D *cobrems_vs_E;

char input_file_name[50]="scalar.in";
char output_file_name[50]="scalar_gen.hddm";

void Usage(void){
  printf("genEtaRegge: generator for eta production based on Regge trajectory formalism.\n");
  printf(" Usage:  genEtaRegge <options>\n");
  printf("   Options:  -N<number of events> (number of events to generate)\n");
  printf("             -O<output.hddm>   (default: eta_gen.hddm)\n");
  printf("             -I<input.in>      (default: eta548.in)\n");
  printf("             -R<run number>    (default: 9000)\n");
  printf("             -h                (Print this message and exit.)\n");
  printf("Coupling constants, photon beam energy range, and eta decay products are\n");
  printf("specified in the <input.in> file.\n");

  exit(0);
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
  int seed=0;
  if (narg==1){
    Usage();
  }
  for(int i=1; i<narg; i++){
    char *ptr = argv[i];
    
    if(ptr[0] == '-'){
      switch(ptr[1]){
      case 'h': Usage(); break;
      case 'I':
	sscanf(&ptr[2],"%s",input_file_name);
	break;
      case 'O':
	sscanf(&ptr[2],"%s",output_file_name);
	break;
      case 'N':
	sscanf(&ptr[2],"%d",&Nevents);
	break;
      case 'R':
	sscanf(&ptr[2],"%d",&runNo);
	break;
      case 'S':
	sscanf(&ptr[2],"%d",&seed);
	break;
      case 'd':
	debug=true;
	break;
      default:
	break;
      }
    }
  }
}

// Non-resonant pi-pi or pi-eta background following Donnachie and Kalashnikova,
// arXiv:0806.3698v1
double BackGroundCrossSection(TLorentzVector &q /* beam */,
			      vector<Particle_t>&particle_types,
			      vector<TLorentzVector>&particles){
  
  int two_particles=particle_types[0]+particle_types[1];
  TLorentzVector p1(0,0,0.,ParticleMass(Proton));
  TLorentzVector p2=particles[2];
  TLorentzVector p=p1-p2;
  TLorentzVector v1=particles[0]-q;
  TLorentzVector v2=particles[1]-q;

  double p1_dot_p2=p1.Dot(p2);
  double q_dot_p=q.Dot(p);
  double q_dot_v1=q.Dot(v1);
  double p_dot_v1=p.Dot(v1);
  double q_dot_v2=q.Dot(v2);
  double p_dot_v2=p.Dot(v2);
  double v1sq=v1.M2();
  double v2sq=v2.M2();
  double v1sq_minus_v2sq=v1sq-v2sq;
  double v1_dot_v2=v1.Dot(v2);
  double psq=p.M2();
  double b1=q_dot_p*v1sq-q_dot_v1*p_dot_v1;
  double b2=q_dot_p*v2sq-q_dot_v2*p_dot_v2;
  TLorentzVector c1=p_dot_v1*q-q_dot_p*v1;
  TLorentzVector c2=p_dot_v2*q-q_dot_p*v2;
  TLorentzVector d1=q_dot_v1*v1-v1sq*q;
  TLorentzVector d2=q_dot_v2*v2-v2sq*q;
  TLorentzVector N1=b1*p1+p1.Dot(c1)*v1+p1.Dot(d1)*p;
  TLorentzVector N2=b2*p1+p1.Dot(c2)*v1+p1.Dot(d2)*p;
  double N1_N1=N1.Dot(N1);
  double N2_N2=N2.Dot(N2);
  double N1_N2=N1.Dot(N2);
  double M1_M2=p_dot_v1*p_dot_v1*q_dot_v2*q_dot_v2 
    - 2.*v1_dot_v2*q_dot_p*(p_dot_v1*q_dot_v2+p_dot_v2*q_dot_v1)
    + v1_dot_v2*v1_dot_v2*q_dot_p*q_dot_p + p_dot_v2*p_dot_v2*q_dot_v1*q_dot_v1
    + v2sq*p_dot_v1*q_dot_v1*q_dot_p + v1sq*p_dot_v2*q_dot_v2*q_dot_p
    + psq*v1_dot_v2*q_dot_v1*q_dot_v2 - psq*v2sq*q_dot_v1*q_dot_v1
    - psq*v1sq*q_dot_v2*q_dot_v2;
  double M1_M1=2.*p_dot_v1*p_dot_v1*q_dot_v1*q_dot_v1
    - 2.*v1sq*p_dot_v1*q_dot_v1*q_dot_p + v1sq*v1sq*q_dot_p*q_dot_p
    - psq*v1sq*q_dot_v1*q_dot_v1; 
  double M2_M2=2.*p_dot_v2*p_dot_v2*q_dot_v2*q_dot_v2
    - 2.*v2sq*p_dot_v2*q_dot_v2*q_dot_p + v2sq*v2sq*q_dot_p*q_dot_p
    - psq*v2sq*q_dot_v2*q_dot_v2;

  double s1=(v1+p1).M2();
  double s2=(v2+p1).M2();
  double t1=(v1-particles[1]).M2();
  double t2=(v2-particles[0]).M2();

  // Rho propagator for top exchange
  double m_rho=0.7685;
  double Gamma_rho=0.1462;
  double m_rhosq_minus_v1sq=m_rho*m_rho-v1sq;
  double Pi_rho_1_sq=1./(m_rhosq_minus_v1sq*m_rhosq_minus_v1sq
			 +m_rho*m_rho*Gamma_rho*Gamma_rho);
  double m_rhosq_minus_v2sq=m_rho*m_rho-v2sq;
  double Pi_rho_2_sq=1./(m_rhosq_minus_v2sq*m_rhosq_minus_v2sq
			 +m_rho*m_rho*Gamma_rho*Gamma_rho);

  double s0=1.;
  // Regge trajectory for rho
  double a_rho_1=0.55+0.8*t1;
  double a_rho_prime=0.8;
  double cos_rho_1=cos(M_PI*a_rho_1);
  double sin_rho_1=sin(M_PI*a_rho_1);
  double regge_rho_1=pow(s1/s0,a_rho_1-1.)*M_PI*a_rho_prime/(sin_rho_1*TMath::Gamma(a_rho_1)); // excluding phase factor
  double regge_rho_1_sq=regge_rho_1*regge_rho_1*0.5*(1.-cos_rho_1); 
  double a_rho_2=0.55+0.8*t2;
  double cos_rho_2=cos(M_PI*a_rho_2);
  double sin_rho_2=sin(M_PI*a_rho_2);
  double regge_rho_2=pow(s2/s0,a_rho_2-1.)*M_PI*a_rho_prime/(sin_rho_2*TMath::Gamma(a_rho_2)); // excluding phase factor
  double regge_rho_2_sq=regge_rho_2*regge_rho_2*0.5*(1.-cos_rho_2);
  double cos_rho_1_rho_2=cos(M_PI*(a_rho_2-a_rho_1));
  double sin_rho_1_rho_2=sin(M_PI*(a_rho_1-a_rho_2)); 

  // omega propagator for top exchange
  double m_omega=0.78265;
  double Gamma_omega=0.00849;
  double m_omegasq_minus_v1sq=m_omega*m_omega-v1sq;
  double Pi_omega_1_sq=1./(m_omegasq_minus_v1sq*m_omegasq_minus_v1sq
			   +m_omega*m_omega*Gamma_omega*Gamma_omega);
  double m_omegasq_minus_v2sq=m_omega*m_omega-v2sq;
  double Pi_omega_2_sq=1./(m_omegasq_minus_v2sq*m_omegasq_minus_v2sq
			   +m_omega*m_omega*Gamma_omega*Gamma_omega);

  // Regge trajectory for omega
  double a_omega_1=0.44+0.9*t1;
  double a_omega_prime=0.9; 
  double cos_omega_1=cos(M_PI*a_omega_1);
  double sin_omega_1=sin(M_PI*a_omega_1);
  double regge_omega_1=pow(s1/s0,a_omega_1-1.)*M_PI*a_omega_prime/(sin_omega_1*TMath::Gamma(a_omega_1)); // excluding phase factor
  double regge_omega_1_sq
    =regge_omega_1*regge_omega_1*0.5*(1.-cos_omega_1);
  double a_omega_2=0.44+0.9*t2;
  double cos_omega_2=cos(M_PI*a_omega_2);
  double sin_omega_2=sin(M_PI*a_omega_2);
  double regge_omega_2=pow(s2/s0,a_omega_2-1.)*M_PI*a_omega_prime/(sin_omega_2*TMath::Gamma(a_omega_2)); // excluding phase factor
  double regge_omega_2_sq
    =regge_omega_2*regge_omega_2*0.5*(1.-cos_omega_2); 
  double cos_omega_1_omega_2=cos(M_PI*(a_omega_2-a_omega_1));
  double sin_omega_1_omega_2=sin(M_PI*(a_omega_2-a_omega_1));

  // rho-omega interference 
  double cos_omega_1_rho_1=cos(M_PI*(a_omega_1-a_rho_1));
  double sin_omega_1_rho_1=sin(M_PI*(a_omega_1-a_rho_1));
  double sin_rho_1_omega_1=-sin_omega_1_rho_1;
  double cos_omega_2_rho_2=cos(M_PI*(a_omega_2-a_rho_2));
  double sin_omega_2_rho_2=sin(M_PI*(a_omega_2-a_rho_2));
  double sin_rho_2_omega_2=-sin_omega_2_rho_2;
  double regge_rho_1_omega_1=regge_rho_1*regge_omega_1;
  double regge_rho_2_omega_2=regge_rho_2*regge_omega_2; 
  double cos_rho_1_omega_2=cos(M_PI*(a_omega_2-a_rho_1));
  double sin_omega_2_rho_1=sin(M_PI*(a_omega_2-a_rho_1));
  double sin_rho_1_omega_2=-sin_omega_2_rho_1;
  double cos_rho_2_omega_1=cos(M_PI*(a_omega_1-a_rho_2));
  double sin_rho_2_omega_1=sin(M_PI*(a_omega_1-a_rho_2));
  double cos_omega_1_rho_1_sum=cos_omega_1_rho_1-cos_rho_1-cos_omega_1+1.; 
  double cos_omega_2_rho_2_sum=cos_omega_2_rho_2-cos_rho_2-cos_omega_2+1.;

  // Coupling constants  
  double g_omega_V=15.;
  double gsq_omega_V=g_omega_V*g_omega_V;
  double g_rho_V=3.4;
  double g_rho_T=11.0; // GeV^-1
  double gsq_rho_T=g_rho_T*g_rho_T;
  double g_rho_V_and_T=g_rho_V+2.*m_p*g_rho_T;
  double gsq_rho_V_and_T=g_rho_V_and_T*g_rho_V_and_T;

  // terms involving complex conjugates of Regge propagators and rho/omega propagtors
  double a1_a1=0.,b1_b1=0.;
  double a2_a2=0.,b2_b2=0.;
  double a1_a2=0.,b1_b2=0.;
  double b1_a1=0,b2_a2=0.,b1_a2_b2_a1=0.;
  double Csq=1.,zetasq=1./2.;
  // Compute square of amplitude    
  if (two_particles==(7+7)){ // pi0 pi0
    double a1_a1_rho_V_and_T_sq=gsq_rho_V_and_T*regge_rho_1_sq*Pi_omega_1_sq;
    double a1_a1_omega_V_sq=gsq_omega_V*regge_omega_1_sq*Pi_rho_1_sq;
    double a1_a1_omega_V_rho_V_and_T
      =g_rho_V_and_T*g_omega_V*regge_rho_1_omega_1
      *Pi_rho_1_sq*Pi_omega_1_sq*( (cos_omega_1_rho_1_sum)
				   *(m_rhosq_minus_v1sq*m_omegasq_minus_v1sq
				     +m_rho*m_omega*Gamma_rho*Gamma_omega)
				   +(sin_omega_1_rho_1-sin_rho_1+sin_omega_1)
				   *(m_rho*Gamma_rho*m_omegasq_minus_v1sq
				   -m_omega*Gamma_omega*m_rhosq_minus_v1sq));
    double a2_a2_rho_V_and_T_sq=gsq_rho_V_and_T*regge_rho_2_sq*Pi_omega_2_sq;
    double a2_a2_omega_V_sq=gsq_omega_V*regge_omega_2_sq*Pi_rho_2_sq;
    double a2_a2_omega_V_rho_V_and_T=g_rho_V_and_T*g_omega_V*regge_rho_2_omega_2
      *Pi_rho_2_sq*Pi_omega_2_sq*( (cos_omega_2_rho_2-cos_rho_2-cos_omega_2+1.)
				   *(m_rhosq_minus_v2sq*m_omegasq_minus_v2sq
				     +m_rho*m_omega*Gamma_rho*Gamma_omega)
				   +(sin_omega_2_rho_2-sin_rho_2+sin_omega_2)
				   *(m_rho*Gamma_rho*m_omegasq_minus_v2sq
				   -m_omega*Gamma_omega*m_rhosq_minus_v2sq));
    double a1_a2_rho_V_and_T_sq
      =gsq_rho_V_and_T*regge_rho_1*regge_rho_2*Pi_omega_1_sq*Pi_omega_2_sq
      *((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
	*(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq 
	  + m_omega*Gamma_omega*m_omega*Gamma_omega)
	+ (sin_rho_1_rho_2-sin_rho_2+sin_rho_1)
	*m_omega*Gamma_omega*(m_omegasq_minus_v2sq- m_omegasq_minus_v1sq));
    double a1_a2_omega_V_sq
      =gsq_omega_V*regge_omega_1*regge_omega_2*Pi_rho_1_sq*Pi_rho_2_sq
      *((cos_omega_1_omega_2-cos_omega_1-cos_omega_2+1.)
	*(m_rhosq_minus_v1sq*m_rhosq_minus_v2sq
	  +m_rho*m_rho*Gamma_rho*Gamma_rho)
	+ (sin_omega_1_omega_2+sin_omega_2-sin_omega_1)
	*m_rho*Gamma_rho*(m_rhosq_minus_v2sq-m_rhosq_minus_v1sq));
    double a1_a2_omega_V_rho_V_and_T_term1=g_rho_V_and_T*g_omega_V
      *Pi_omega_1_sq*Pi_rho_2_sq*regge_rho_1*regge_omega_2
      *((cos_rho_1_omega_2-cos_rho_1-cos_omega_2+1.)
	*(m_omegasq_minus_v1sq*m_rhosq_minus_v2sq
	  +m_rho*Gamma_rho*m_omega*Gamma_omega)
	+(sin_rho_1_omega_2+sin_omega_2-sin_rho_1)
	*(m_omega*Gamma_omega*m_rhosq_minus_v2sq
	  -m_rho*Gamma_rho*m_omegasq_minus_v1sq));
     double a1_a2_omega_V_rho_V_and_T_term2=g_rho_V_and_T*g_omega_V
      *Pi_omega_1_sq*Pi_omega_2_sq*regge_rho_2*regge_omega_1
      *((cos_rho_2_omega_1-cos_omega_1-cos_rho_2+1.)
	*(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq
	  +m_omega*Gamma_omega*m_omega*Gamma_omega)
	+(sin_rho_2-sin_omega_1-sin_rho_2_omega_1)
	*m_omega*Gamma_omega*v1sq_minus_v2sq);
     double b1_a1_rho_T_rho_V_and_T
       =-4.*g_rho_T*g_rho_V_and_T*regge_rho_1_sq*Pi_omega_1_sq;
     double b1_a1_omega_V_rho_V_and_T
       =-g_rho_T*g_omega_V*Pi_rho_1_sq*Pi_omega_1_sq*regge_omega_1*regge_rho_1
       *((cos_omega_1_rho_1_sum)*(m_rhosq_minus_v1sq*m_omegasq_minus_v1sq
				  +m_rho*Gamma_rho*m_omega*Gamma_omega)
	 +(sin_rho_1_omega_1-sin_rho_1+sin_omega_1)
		 *(m_omega*Gamma_omega*m_rhosq_minus_v1sq
		   -m_rho*Gamma_rho*m_omegasq_minus_v1sq));
     double b2_a2_rho_T_rho_V_and_T
       =-4.*g_rho_T*g_rho_V_and_T*regge_rho_2_sq*Pi_omega_2_sq;
     double b2_a2_omega_V_rho_V_and_T
      =-g_rho_T*g_omega_V*Pi_rho_2_sq*Pi_omega_2_sq*regge_rho_2*regge_omega_2
      *((cos_omega_2_rho_2_sum)*(m_rhosq_minus_v2sq*m_omegasq_minus_v2sq
				 +m_rho*Gamma_rho*m_omega*Gamma_omega)
	+(sin_rho_2_omega_2-sin_rho_2+sin_omega_2)
	*(m_omega*Gamma_omega*m_rhosq_minus_v2sq
	  -m_rho*Gamma_rho*m_omegasq_minus_v2sq));
     double b1_a2_rho_T_rho_V_and_T
       =-regge_rho_1*regge_rho_2*g_rho_T*g_rho_V_and_T
       *Pi_omega_1_sq*Pi_omega_2_sq
       *((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
	 *(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq
	   +m_omega*Gamma_omega*m_omega*Gamma_omega)
	 +(sin_rho_1_rho_2-sin_rho_1+sin_rho_2)
	 *m_omega*Gamma_omega*(m_omegasq_minus_v2sq-m_omegasq_minus_v1sq));
     double b1_a2_omega_V_rho_T=-g_rho_T*g_omega_V*Pi_omega_1_sq*Pi_rho_2_sq
       *regge_rho_1*regge_omega_2*((cos_rho_1_omega_2-cos_rho_1-cos_omega_2+1.)
				 *(m_omegasq_minus_v1sq*m_rhosq_minus_v2sq
				   +m_rho*m_omega*Gamma_omega*Gamma_rho)
				  +(sin_rho_1_omega_2-sin_rho_1+sin_omega_2)
				   *(m_omega*Gamma_omega*m_rhosq_minus_v2sq
				     -m_rho*Gamma_rho*m_omegasq_minus_v1sq));
     double b2_a1_rho_T_rho_V_and_T=b1_a2_rho_T_rho_V_and_T; 
     double b2_a1_omega_V_rho_T=-g_rho_T*g_omega_V*Pi_rho_1_sq*Pi_omega_2_sq
       *regge_rho_2*regge_omega_1*((cos_rho_2_omega_1-cos_rho_2-cos_omega_1+1.)
				 *(m_rhosq_minus_v1sq*m_omegasq_minus_v2sq
				   +m_rho*m_omega*Gamma_rho*Gamma_omega)
				  +(sin_rho_2_omega_1-sin_rho_2+sin_omega_1)
				   *(m_omega*Gamma_omega*m_rhosq_minus_v1sq
				     -m_rho*Gamma_rho*m_omegasq_minus_v2sq));
     
     b1_b1=4.*gsq_rho_T*regge_rho_1_sq*Pi_omega_1_sq; 
     b2_b2=4.*gsq_rho_T*regge_rho_2_sq*Pi_omega_2_sq;
     b1_b2=2.*gsq_rho_T*regge_rho_1*regge_rho_2
       *Pi_omega_1_sq*Pi_omega_2_sq
       *((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
	 *(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq
	   +m_omega*m_rho*Gamma_rho*Gamma_rho)
	 +(sin_rho_1_rho_2+sin_rho_2-sin_rho_1)
	 *m_omega*Gamma_omega*(m_omegasq_minus_v2sq-m_omegasq_minus_v1sq));
     
     a1_a1=a1_a1_rho_V_and_T_sq + (1./9.)*a1_a1_omega_V_sq
      + (1./6.)*a1_a1_omega_V_rho_V_and_T; 
     a2_a2=a2_a2_rho_V_and_T_sq + (1./9.)*a2_a2_omega_V_sq
       + (1./6.)*a2_a2_omega_V_rho_V_and_T;  
     a1_a2=(1./2.)*a1_a2_rho_V_and_T_sq + (1./18.)*a1_a2_omega_V_sq	      
       +(1./6.)*(a1_a2_omega_V_rho_V_and_T_term1
		 +a1_a2_omega_V_rho_V_and_T_term2);

     b1_a1=b1_a1_rho_T_rho_V_and_T+(1./3.)*b1_a1_omega_V_rho_V_and_T; 
     b2_a2=b2_a2_rho_T_rho_V_and_T+(1./3.)*b2_a2_omega_V_rho_V_and_T;
     b1_a2_b2_a1=b1_a2_rho_T_rho_V_and_T+b2_a1_rho_T_rho_V_and_T
       +(1./3.)*(b1_a2_omega_V_rho_T+b2_a1_omega_V_rho_T);
   
  }
  else if (two_particles==(7+17)){ // pi0 eta
    Csq=2./3.;
    zetasq=1.;

    double a1_a1_rho_V_and_T_sq=gsq_rho_V_and_T*regge_rho_1_sq*Pi_rho_1_sq;
    double a1_a1_omega_V_sq=gsq_omega_V*regge_omega_1_sq*Pi_omega_1_sq;
    double a1_a1_omega_V_rho_V_and_T=g_rho_V_and_T*g_omega_V*regge_rho_1_omega_1
      *Pi_rho_1_sq*Pi_omega_1_sq*( (cos_omega_1_rho_1_sum)
				   *(m_rhosq_minus_v1sq*m_omegasq_minus_v1sq
				     +m_rho*m_omega*Gamma_rho*Gamma_omega)
				   +(sin_omega_1_rho_1-sin_rho_1+sin_omega_1)
				   *(m_rho*Gamma_rho*m_omegasq_minus_v1sq
				   -m_omega*Gamma_omega*m_rhosq_minus_v1sq));
    double a2_a2_rho_V_and_T_sq=gsq_rho_V_and_T*regge_rho_2_sq*Pi_omega_2_sq;
    double a2_a2_omega_V_sq=gsq_omega_V*regge_omega_2_sq*Pi_rho_2_sq;
    double a2_a2_omega_V_rho_V_and_T=g_rho_V_and_T*g_omega_V*regge_rho_2_omega_2
      *Pi_rho_2_sq*Pi_omega_2_sq*( (cos_omega_2_rho_2_sum)
				 *(m_rhosq_minus_v2sq*m_omegasq_minus_v2sq
				   +m_rho*m_omega*Gamma_rho*Gamma_omega)
				  + (sin_omega_2_rho_2-sin_rho_2+sin_omega_2)
				   *(m_rho*Gamma_rho*m_omegasq_minus_v2sq
				     -m_omega*Gamma_omega*m_rhosq_minus_v2sq));
    double a1_a2_rho_V_and_T_sq
      =gsq_rho_V_and_T*regge_rho_1*regge_rho_2*Pi_rho_1_sq*Pi_omega_2_sq
      *((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
	*(m_rhosq_minus_v1sq*m_omegasq_minus_v2sq 
	  + m_rho*Gamma_rho*m_omega*Gamma_omega)
	+ (sin_rho_1_rho_2-sin_rho_2+sin_rho_1)
	*(m_rho*Gamma_rho*m_omegasq_minus_v2sq
	  - m_omega*Gamma_omega*m_rhosq_minus_v1sq));
    double a1_a2_omega_V_rho_V_and_T_term1=g_rho_V_and_T*g_omega_V
      *Pi_rho_1_sq*Pi_rho_2_sq*regge_rho_1*regge_omega_2
      *((cos_rho_1_omega_2-cos_rho_1-cos_omega_2+1.)
	*(m_rhosq_minus_v1sq*m_rhosq_minus_v2sq+m_rho*Gamma_rho*m_rho*Gamma_rho)
	+(sin_rho_1_omega_2+sin_omega_2-sin_rho_1)*m_rho*Gamma_rho*v1sq_minus_v2sq);
    double a1_a2_omega_V_rho_V_and_T_term2=g_rho_V_and_T*g_omega_V
      *Pi_omega_1_sq*Pi_omega_2_sq*regge_rho_2*regge_omega_1
      *((cos_rho_2_omega_1-cos_omega_1-cos_rho_2+1.)
	*(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq
	  +m_omega*Gamma_omega*m_omega*Gamma_omega)
	+(sin_rho_2-sin_omega_1-sin_rho_2_omega_1)
	*m_omega*Gamma_omega*v1sq_minus_v2sq);
    double a1_a2_omega_V_sq
      =gsq_omega_V*regge_omega_1*regge_omega_2*Pi_omega_1_sq*Pi_rho_2_sq
      *((cos_omega_1_omega_2-cos_omega_1-cos_omega_2+1.)
	*(m_omegasq_minus_v1sq*m_rhosq_minus_v2sq
	  +m_rho*m_omega*Gamma_rho*Gamma_omega)
	+ (sin_omega_1_omega_2+sin_omega_2-sin_omega_1)
	*(m_omega*Gamma_omega*m_rhosq_minus_v2sq
	  -m_rho*Gamma_rho*m_omegasq_minus_v1sq));
    double b1_a1_rho_T_rho_V_and_T
      =-4.*g_rho_T*g_rho_V_and_T*regge_rho_1_sq*Pi_rho_1_sq;
    double b1_a1_omega_V_rho_V_and_T
      =-g_rho_T*g_omega_V*Pi_rho_1_sq*Pi_omega_1_sq*regge_omega_1*regge_rho_1
      *((cos_omega_1_rho_1_sum)*(m_rhosq_minus_v1sq*m_omegasq_minus_v1sq
			       +m_rho*Gamma_rho*m_omega*Gamma_omega)
      +(sin_rho_1_omega_1-sin_rho_1+sin_omega_1)
		 *(m_rho*Gamma_rho*m_omegasq_minus_v1sq
		   -m_omega*Gamma_omega*m_rhosq_minus_v1sq));
    double b2_a2_rho_T_rho_V_and_T
      =-4.*g_rho_T*g_rho_V_and_T*regge_rho_2_sq*Pi_omega_2_sq;
    double b2_a2_omega_V_rho_V_and_T
      =-g_rho_T*g_omega_V*Pi_rho_2_sq*Pi_omega_2_sq*regge_rho_2*regge_omega_2
      *((cos_omega_2_rho_2_sum)*(m_rhosq_minus_v2sq*m_omegasq_minus_v2sq
				 +m_rho*Gamma_rho*m_omega*Gamma_omega)
	+(sin_rho_2_omega_2-sin_rho_2+sin_omega_2)
	*(m_omega*Gamma_omega*m_rhosq_minus_v2sq
	  -m_rho*Gamma_rho*m_omegasq_minus_v2sq));
    double b1_a2_rho_T_rho_V_and_T
      =-regge_rho_1*regge_rho_2*g_rho_T*g_rho_V_and_T
      *Pi_rho_1_sq*Pi_omega_2_sq*((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
				*(m_rhosq_minus_v1sq*m_omegasq_minus_v2sq
				  +m_omega*Gamma_omega*m_rho*Gamma_rho)
				  +(sin_rho_1_rho_2+sin_rho_1-sin_rho_2)
				  *(m_rho*Gamma_rho*m_omegasq_minus_v2sq
				    -m_omega*Gamma_omega*m_rhosq_minus_v1sq));
    double b1_a2_omega_V_rho_T=-g_rho_T*g_omega_V*Pi_rho_1_sq*Pi_rho_2_sq
      *regge_rho_1*regge_omega_2*((cos_rho_1_omega_2-cos_rho_1-cos_omega_2+1.)
				 *(m_rhosq_minus_v1sq*m_rhosq_minus_v2sq
				   +m_rho*m_rho*Gamma_rho*Gamma_rho)
				  +(sin_rho_1_omega_2+sin_rho_1-sin_omega_2)
				  *m_rho*Gamma_rho*(m_rhosq_minus_v2sq
						    -m_rhosq_minus_v1sq));
    double b2_a1_rho_T_rho_V_and_T=b1_a2_rho_T_rho_V_and_T; 
    double b2_a1_omega_V_rho_T=-g_rho_T*g_omega_V*Pi_omega_1_sq*Pi_omega_2_sq
      *regge_rho_1*regge_omega_2*((cos_rho_2_omega_1-cos_rho_2-cos_omega_1+1.)
				 *(m_omegasq_minus_v1sq*m_omegasq_minus_v2sq
				   +m_omega*m_omega*Gamma_omega*Gamma_omega)
				  +(sin_rho_2_omega_1+sin_rho_2-sin_omega_1)
				  *m_omega*Gamma_omega*(m_omegasq_minus_v1sq
							-m_omegasq_minus_v2sq));
     
    b1_b1=(4./9.)*gsq_rho_T*regge_rho_1_sq*Pi_rho_1_sq; 
    b2_b2=(4./9.)*gsq_rho_T*regge_rho_2_sq*Pi_omega_2_sq;
    b1_b2=(2./9.)*gsq_rho_T*regge_rho_1*regge_rho_2
      *Pi_rho_1_sq*Pi_omega_2_sq*((cos_rho_1_rho_2-cos_rho_1-cos_rho_2+1.)
				  *(m_rhosq_minus_v1sq*m_rhosq_minus_v2sq
				    +m_rho*m_rho*Gamma_rho*Gamma_rho)
				  - (sin_rho_1_rho_2-sin_rho_2+sin_rho_1)
				  *(m_omega*Gamma_omega*m_rhosq_minus_v1sq
				    - m_rho*Gamma_rho*m_omegasq_minus_v2sq));
   
    a1_a1=(1./9.)*a1_a1_rho_V_and_T_sq + a1_a1_omega_V_sq 
      + (1./6.)*a1_a1_omega_V_rho_V_and_T; 
    a2_a2=(1./9.)*a2_a2_rho_V_and_T_sq + a2_a2_omega_V_sq 
      + (1./6.)*a2_a2_omega_V_rho_V_and_T; 
    a1_a2=(1./18.)*a1_a2_rho_V_and_T_sq + (1/2.)*a1_a2_omega_V_sq	      
      +(1./6.)*(a1_a2_omega_V_rho_V_and_T_term1
		+a1_a2_omega_V_rho_V_and_T_term2);
          
    b1_a1=(1./9.)*b1_a1_rho_T_rho_V_and_T+(1./3.)*b1_a1_omega_V_rho_V_and_T;
    b2_a2=(1./9.)*b2_a2_rho_T_rho_V_and_T+(1./3.)*b2_a2_omega_V_rho_V_and_T;
    b1_a2_b2_a1=(1./9.)*(b1_a2_rho_T_rho_V_and_T+b2_a1_rho_T_rho_V_and_T)
      +(1./3.)*(b1_a2_omega_V_rho_T+b2_a1_omega_V_rho_T);
  }
  double T=Csq*((m_p_sq-p1_dot_p2)*(a1_a1*M1_M1 + a1_a2*M1_M2 + a2_a2*M2_M2)
		+ 2.*(a1_a1*N1_N1 + a1_a2*N1_N2 + a2_a2*N2_N2)
		+ 2.*m_p*(b1_a1*N1_N1 + b1_a2_b2_a1*N1_N2 + b2_a2*N2_N2)
		+ (m_p_sq+p1_dot_p2)*(b1_b1*N1_N1 + b1_b2*N1_N2 + b2_b2*N2_N2));
    
  double Msq=(particles[0]+particles[1]).M2();
  double s=(q+p1).M2();
  double s_minus_mp_sq=s-m_p_sq;
  double g0_sq=8.8; // GeV^-2
  double m1=ParticleMass(particle_types[0]);
  double m2=ParticleMass(particle_types[1]);
  double m1_plus_m2=m1+m2;
  double m1_minus_m2=m1-m2;
  double k=sqrt((Msq-m1_plus_m2*m1_plus_m2)*(Msq-m1_minus_m2*m1_minus_m2))
    /(2.*sqrt(Msq));
  double hbarc_sq=389.; // Convert to micro-barns
  double xsec=-zetasq*g0_sq*k*T*hbarc_sq
    /(256*M_PI*M_PI*M_PI*M_PI*s_minus_mp_sq*s_minus_mp_sq);
  
  return xsec;				 
}
/*
double InterferenceCrossSection(TLorentzVector &q /* beam,
				vector<Particle_t>&particle_types,
				vector<TLorentzVector>&particles){ 
  int two_particles=particle_types[0]+particle_types[1];
  
  // Four vectors
  TLorentzVector p1(0,0,0.,ParticleMass(Proton));
  TLorentzVector p2=particles[2];
  TLorentzVector p=p1-p2;
  TLorentzVector v1=particles[0]-q;
  TLorentzVector v2=particles[1]-q;

  // Mandelstam variables
  double s=(q+p1).M2();
  double t=(p1-p2).M2();
  double s1=(v1+p1).M2();
  double s2=(v2+p1).M2();
  double t1=(v1-particles[1]).M2();
  double t2=(v2-particles[0]).M2();

  // dot products 
  double p1_dot_p2=p1.Dot(p2);
  double q_dot_p=q.Dot(p);
  double q_dot_v1=q.Dot(v1);
  double p_dot_v1=p.Dot(v1);
  double q_dot_v2=q.Dot(v2);
  double p_dot_v2=p.Dot(v2);
  double v1sq=v1.M2();
  double v2sq=v2.M2();
  double v1sq_minus_v2sq=v1sq-v2sq;
  double v1_dot_v2=v1.Dot(v2);
  double psq=p.M2();
  double b1=q_dot_p*v1sq-q_dot_v1*p_dot_v1;
  double b2=q_dot_p*v2sq-q_dot_v2*p_dot_v2;
  TLorentzVector c1=p_dot_v1*q-q_dot_p*v1;
  TLorentzVector c2=p_dot_v2*q-q_dot_p*v2;
  TLorentzVector d1=q_dot_v1*v1-v1sq*q;
  TLorentzVector d2=q_dot_v2*v2-v2sq*q;
  TLorentzVector N1=b1*p1+p1.Dot(c1)*v1+p1.Dot(d1)*p;
  TLorentzVector N2=b2*p1+p1.Dot(c2)*v1+p1.Dot(d2)*p;
  TLorentzVector N=(q_dot_p)*p1-q.Dot(p1)*p;
  double N_N1=N.Dot(N1);
  double N_N2=N.Dot(N2);
  double M_M1=3*b1*q_dot_p+v1.Dot(c1)*q_dot_p+p_dot_v1*q.Dot(c1)
    +p.Dot(d1)*q_dot_p-psq*q.Dot(d1);  
  double M_M2=3*b2*q_dot_p+v2.Dot(c2)*q_dot_p+p_dot_v2*q.Dot(c2)
    +p.Dot(d2)*q_dot_p-psq*q.Dot(d2);


  // Rho propagator for top exchange for double-exchange diagrams
  double m_rho=0.7685;
  double Gamma_rho=0.1462;
  double m_rhosq_minus_v1sq=m_rho*m_rho-v1sq;
  double Pi_rho_1_sq=1./(m_rhosq_minus_v1sq*m_rhosq_minus_v1sq
			 +m_rho*m_rho*Gamma_rho*Gamma_rho);
  double m_rhosq_minus_v2sq=m_rho*m_rho-v2sq;
  double Pi_rho_2_sq=1./(m_rhosq_minus_v2sq*m_rhosq_minus_v2sq
			 +m_rho*m_rho*Gamma_rho*Gamma_rho);
  
  double s0=1.;
  // Regge trajectory for rho
  double a_rho_1=0.55+0.8*t1;
  double a_rho=0.55+0.8*t;
  double a_rho_prime=0.8;
  double cos_rho=cos(M_PI*a_rho);
  double sin_rho=sin(M_PI*a_rho);
  double cos_rho_1=cos(M_PI*a_rho_1);
  double sin_rho_1=sin(M_PI*a_rho_1);
  double regge_rho=pow(s/s0,a_rho-1.)*M_PI*a_rho_prime/(sin_rho*TMath::Gamma(a_rho)); // excluding phase factor
  double regge_rho_1=pow(s1/s0,a_rho_1-1.)*M_PI*a_rho_prime/(sin_rho_1*TMath::Gamma(a_rho_1)); // excluding phase factor
  double regge_rho_1_sq=regge_rho_1*regge_rho_1*0.5*(1.-cos_rho_1); 
  double a_rho_2=0.55+0.8*t2;
  double cos_rho_2=cos(M_PI*a_rho_2);
  double sin_rho_2=sin(M_PI*a_rho_2);
  double regge_rho_2=pow(s2/s0,a_rho_2-1.)*M_PI*a_rho_prime/(sin_rho_2*TMath::Gamma(a_rho_2)); // excluding phase factor
  double regge_rho_2_sq=regge_rho_2*regge_rho_2*0.5*(1.-cos_rho_2);
  double cos_rho_1_rho_2=cos(M_PI*(a_rho_2-a_rho_1));
  double sin_rho_1_rho_2=sin(M_PI*(a_rho_1-a_rho_2)); 

 // omega propagator for top exchange
  double m_omega=0.78265;
  double Gamma_omega=0.00849;
  double m_omegasq_minus_v1sq=m_omega*m_omega-v1sq;
  double Pi_omega_1_sq=1./(m_omegasq_minus_v1sq*m_omegasq_minus_v1sq
			   +m_omega*m_omega*Gamma_omega*Gamma_omega);
  double m_omegasq_minus_v2sq=m_omega*m_omega-v2sq;
  double Pi_omega_2_sq=1./(m_omegasq_minus_v2sq*m_omegasq_minus_v2sq
			   +m_omega*m_omega*Gamma_omega*Gamma_omega);

  // Regge trajectory for omega
  double a_omega=0.44+0.9*t;
  double a_omega_1=0.44+0.9*t1;
  double a_omega_prime=0.9; 
  double cos_omega=cos(M_PI*a_omega);
  double sin_omega=sin(M_PI*a_omega);
  double cos_omega_1=cos(M_PI*a_omega_1);
  double sin_omega_1=sin(M_PI*a_omega_1); 
  double regge_omega=pow(s/s0,a_omega-1.)*M_PI*a_omega_prime/(sin_omega*TMath::Gamma(a_omega)); // excluding phase factor
  double regge_omega_1=pow(s1/s0,a_omega_1-1.)*M_PI*a_omega_prime/(sin_omega_1*TMath::Gamma(a_omega_1)); // excluding phase factor
  double regge_omega_1_sq
    =regge_omega_1*regge_omega_1*0.5*(1.-cos_omega_1);
  double a_omega_2=0.44+0.9*t2;
  double cos_omega_2=cos(M_PI*a_omega_2);
  double sin_omega_2=sin(M_PI*a_omega_2);
  double regge_omega_2=pow(s2/s0,a_omega_2-1.)*M_PI*a_omega_prime/(sin_omega_2*TMath::Gamma(a_omega_2)); // excluding phase factor
  double regge_omega_2_sq
    =regge_omega_2*regge_omega_2*0.5*(1.-cos_omega_2); 
  double cos_omega_1_omega_2=cos(M_PI*(a_omega_2-a_omega_1));
  double sin_omega_1_omega_2=sin(M_PI*(a_omega_2-a_omega_1));

  // Breit-Wigner real and imaginary terms
  double ReB=0, ImB=0;
  // phase between background and resonant signal
  double phase=0.;

  // terms involving complex conjugates of Regge propagators and rho/omega propagtors
  double Re_a1_aS=0, Im_a1_aS=0, Re_a2_aS=0., Im_a2_aS=0.;
  double Re_b1_aS=0, Im_b1_aS=0, Re_b2_aS=0., Im_b2_aS=0.;
  double Re_a1_bS=0, Im_a1_bS=0., Re_b1_bS=0., Im_b2_bS=0.;

  // Compute square of amplitude    
  if (two_particles==(7+7)){ // pi0 pi0
    double Re_a1_aS_rho_V_and_T_sq
      =(1./2.)*gsq_rho_V_and_T*g_rho_S*regge_rho*regge_rho_1*Pi_omega_1_sq
      *((cos_rho_rho_1-cos_rho-cos_rho_1+1)
	*(m_omegasq_minus_v1sq*ReB+m_omega*Gamma_omega*ImB)
	+(sin_rho_rho_1-sin_rho_1+sin_rho)
	*(m_omega*Gamma_omega*ReB-m_omegasq_minus_v1sq*ImB));
    double Im_a1_aS_rho_V_and_T_sq
      =-(1./2)*gsq_rho_V_and_T*g_rho_S*regge_rho*regge_rho_1*Pi_omega_1_sq
      *((sin_rho_1-sin_rho-sin_rho_rho_1)
	*(m_omegasq_minus_v1sq*ReB+m_omega*Gamma_omega*ImB)
	+(cos_rho_rho_1-cos_rho-cos_rho_1+1)
	*(m_omega*Gamma_omega*ReB-m_omegasq_minus_v1sq*ImB));
    
    Re_a1_aS=Re_a1_aS_rho_V_and_T_sq;
    Im_a1_aS=Im_a1_aS_rho_V_and_T_sq;
  }

  // Real and imaginary parts of the trace term
  double ReT=(m_p_sq-p1_dot_p2)*(Re_a1_aS*M_M1+Re_a2_aS*M_M2
				 +2.*Re_a1_aS*N_N1+2.*Re_a2_aS*N_N2)
    +2.*m_p*((Re_b1_aS+Re_a1_bS)*N_N1+(Re_b2_aS+Re_a2_bS)*N_N2)
    +(m_p_sq+p1_dot_p2)*(Re_b1_bS*N_N1+Re_b2_bS*N_N2);
  double ImT=(m_p_sq-p1_dot_p2)*(Im_a1_aS*M_M1+Im_a2_aS*M_M2
				 +2.*Im_a1_aS*N_N1+2.*Im_a2_aS*N_N2)
    +2.*m_p*((Im_b1_aS+Im_a1_bS)*N_N1+(Im_b2_aS+Im_a2_bS)*N_N2)
    +(m_p_sq+p1_dot_p2)*(Im_b1_bS*N_N1+Im_b2_bS*N_N2);
  
  double Msq=(particles[0]+particles[1]).M2();
  double g0=sqrt(1.514e4/137); 
  double s_minus_mp_sq=s-m_p_sq;
  double m1=ParticleMass(particle_types[0]);
  double m2=ParticleMass(particle_types[1]);
  double m1_plus_m2=m1+m2;
  double m1_minus_m2=m1-m2;
  double k=sqrt((Msq-m1_plus_m2*m1_plus_m2)*(Msq-m1_minus_m2*m1_minus_m2))
    /(2.*sqrt(Msq));
  double hbarc_sq=389.; // Convert to micro-barns
  double xsec=-zeta*gR*g0*k*hbarc_sq*(cos(phase)*ReT+sin(phase)*ImT)
    /(256*M_PI*M_PI*M_PI*M_PI*s_minus_mp_sq*s_minus_mp_sq*(ReB*ReB+ImB*ImB));
}
*/

// Get parameters for Breit-Wigner distribution for resonance shape
void GetResonanceParameters(double m1,double m2, double m_S_sq,double &ReB,
			    double &ImB){
  double m1_plus_m2=m1+m2;
  double m1_minus_m2=m1-m2;  
  bool got_pipi=(fabs(m1_minus_m2)>0.01)?false:true;

  // "No structure" model for a0(980)/f0(980)
  // masses
  double M_S=sqrt(m_S_sq);
  double MK0=ParticleMass(KShort);
  double MKPlus=ParticleMass(KPlus);

  // coupling constants
  double gK=3.05; 
  double g_m1m2=2.82;
  if (got_pipi){
    gK=5.006;
    g_m1m2=1.705;
  }
  double gKsq=gK*gK;    
  double g_m1m2_sq=g_m1m2*g_m1m2;
    
  // kinematic factors
  double rhoK0sq=1.-4.*MK0*MK0/m_S_sq;
  double rhoKPlussq=1.-4.*MKPlus*MKPlus/m_S_sq;
  double rho_m1m2
   =sqrt((1.-m1_plus_m2*m1_plus_m2/m_S_sq)*(1-m1_minus_m2*m1_minus_m2/m_S_sq));

      // Real and imaginary parts of BW amplitude
  ReB=m_S_sq_R-m_S_sq;
  if (M_S<2.*MKPlus){
    ReB+=gKsq/(32.*M_PI)*sqrt(-rhoKPlussq);
  }
  if (M_S<2.*MK0){
    ReB+=gKsq/(32.*M_PI)*sqrt(-rhoK0sq);
  }
  ImB=g_m1m2_sq/(16*M_PI)*rho_m1m2;
  if (M_S>2.*MKPlus){
    ImB+=gKsq/(32.*M_PI)*sqrt(rhoKPlussq);
  }
  if (M_S>2.*MK0){
    ImB+=gKsq/(32.*M_PI)*sqrt(rhoK0sq);
  }
}


// Cross section dsigma/(dt/dM/dOmega) from Donnachie and Kalashnikova
double CrossSection(double m1,double m2, double ms_sq, double s, double t,
		    double gR,double ReB,double ImB){
  // Kinematic factors
  double mp_sq_minus_s=m_p_sq-s;
  double mp_sq_plus_s=m_p_sq+s;
  double mp_sq_minus_s_sq=mp_sq_minus_s*mp_sq_minus_s;
  double temp1=ms_sq*mp_sq_plus_s-mp_sq_minus_s_sq;
  double temp2=mp_sq_minus_s*sqrt(mp_sq_minus_s_sq-2.*ms_sq*mp_sq_plus_s
				  +ms_sq*ms_sq);
  double t1=(temp1+temp2)/(2.*s);
  double t2=(temp1-temp2)/(2.*s);
  double t_minus_ms_sq=t-ms_sq;
  double kin1=s*(t-t1)*(t-t2);
  double kin_aS_aS=0.5*kin1+0.25*t*t_minus_ms_sq*t_minus_ms_sq;
  double kin_aS_bS=0.5*m_p*kin1;
  double kin_bS_bS=0.125*(4.*m_p_sq-t)*kin1;

  // Coupling constants 
  double g_omega_V=15.;
  double gsq_omega_V=g_omega_V*g_omega_V;
  double g_rho_V=3.4;
  double gsq_rho_V=g_rho_V*g_rho_V;
  double g_rho_T=11.0; // GeV^-1
  double g_rho_V_and_T=g_rho_V+2.*m_p*g_rho_T;
  double gsq_rho_V_and_T=g_rho_V_and_T*g_rho_V_and_T;
  double g_rho_S_gamma=sqrt(gsq_rho_S_gamma);
  double g_omega_S_gamma=sqrt(gsq_omega_S_gamma);
  double gRsq=gR*gR;

  // s scale for regge trajectories
  double s0=1.;

  // Regge trajectory for omega
  double a_omega=0.44+0.9*t;
  double a_omega_prime=0.9;
  double cos_omega=cos(M_PI*a_omega);
  double regge_omega=pow(s/s0,a_omega-1.)*M_PI*a_omega_prime/(sin(M_PI*a_omega)*TMath::Gamma(a_omega)); // excluding phase factor
  double regge_omega_sq=regge_omega*regge_omega*0.5*(1.-cos_omega);

  // Regge cuts for omega
  double a_omega_P=0.52+0.196*t; // Pomeron
  double a_omega_f2=0.112+0.428*t;
  double dc=1.36;
  double regge_omega_P_cut=exp(dc*t)*pow(s/s0,a_omega_P-1.);
  double regge_omega_f2_cut=exp(dc*t)*pow(s/s0,a_omega_f2-1.);
  double C_omega_P_cut=0.271;
  double C_omega_f2_cut=0.926;
  double Csq_omega_P_cut=C_omega_P_cut*C_omega_P_cut;
  double Csq_omega_f2_cut=C_omega_f2_cut*C_omega_f2_cut;
  double C_omega_f2_C_omega_P=C_omega_f2_cut*C_omega_P_cut;

  // Regge trajectory for rho
  double a_rho=0.55+0.8*t;
  double a_rho_prime=0.8;
  double cos_rho=cos(M_PI*a_rho);
  double regge_rho=pow(s/s0,a_rho-1.)*M_PI*a_rho_prime/(sin(M_PI*a_rho)*TMath::Gamma(a_rho)); // excluding phase factor
  double regge_rho_sq=regge_rho*regge_rho*0.5*(1.-cos_rho);

  // Regge cuts for rho
  double a_rho_P=0.64+0.16*t; // Pomeron
  double a_rho_f2=0.222+0.404*t;
  double regge_rho_P_cut=exp(dc*t)*pow(s/s0,a_rho_P-1.);
  double regge_rho_f2_cut=exp(dc*t)*pow(s/s0,a_rho_f2-1.);
  double C_rho_P_cut=-2.86;
  double C_rho_f2_cut=-4.64;
  double Csq_rho_P_cut=C_rho_P_cut*C_rho_P_cut;
  double Csq_rho_f2_cut=C_rho_f2_cut*C_rho_f2_cut;

  // rho-omega interference
  double cos_rho_omega_sum=cos(M_PI*(a_rho-a_omega))-cos(M_PI*a_rho)
    -cos(M_PI*a_omega)+1.;
  // Regge propagator products
  double regge_omega_omega_P_cut=regge_omega*regge_omega_P_cut
    *(cos(M_PI*(a_omega-0.5*a_omega_P))-cos(M_PI_2*a_omega_P)) ;
  double regge_omega_omega_f2_cut=regge_omega*regge_omega_f2_cut
    *(cos(M_PI*(a_omega-0.5*a_omega_f2))-cos(M_PI_2*a_omega_f2));
  double regge_rho_omega_P_cut=regge_rho*regge_omega_P_cut
    *(cos(M_PI*(a_rho-0.5*a_omega_P))-cos(M_PI_2*a_omega_P)); 
  double regge_rho_omega_f2_cut=regge_rho*regge_omega_f2_cut
    *(cos(M_PI*(a_rho-0.5*a_omega_f2))-cos(M_PI_2*a_omega_f2));
  double regge_omega_f2_omega_P_cuts
    =regge_omega_f2_cut*regge_omega_P_cut*cos(M_PI_2*(a_omega_f2-a_omega_P));
  double regge_omega_P_cut_sq=regge_omega_P_cut*regge_omega_P_cut;
  double regge_omega_f2_cut_sq=regge_omega_f2_cut*regge_omega_f2_cut;
  double regge_rho_P_cut_sq=regge_rho_P_cut*regge_rho_P_cut;
  double regge_rho_f2_cut_sq=regge_rho_f2_cut*regge_rho_f2_cut;
  double regge_rho_rho_P_cut=regge_rho*regge_rho_P_cut
    *(cos(M_PI*(a_rho-0.5*a_rho_P))-cos(M_PI_2*a_rho_P));
  double regge_rho_rho_f2_cut=regge_rho*regge_rho_f2_cut
    *(cos(M_PI*(a_rho-0.5*a_rho_f2))-cos(M_PI_2*a_rho_f2));
  double regge_omega_rho_P_cut=regge_omega*regge_rho_P_cut
    *(cos(M_PI*(a_omega-0.5*a_rho_P))-cos(M_PI_2*a_rho_P));
  double regge_omega_rho_f2_cut=regge_omega*regge_rho_f2_cut
    *(cos(M_PI*(a_omega-0.5*a_rho_f2))-cos(M_PI_2*a_rho_f2));
  double regge_rho_P_omega_P_cuts=regge_rho_P_cut*regge_omega_P_cut
    *cos(M_PI_2*(a_rho_P-a_omega_P));
  double regge_rho_P_omega_f2_cuts=regge_rho_P_cut*regge_omega_f2_cut
    *cos(M_PI_2*(a_rho_P-a_omega_f2));
  double regge_rho_f2_omega_P_cuts=regge_rho_f2_cut*regge_omega_P_cut
    *cos(M_PI_2*(a_rho_f2-a_omega_P));
  double regge_rho_f2_omega_f2_cuts=regge_rho_f2_cut*regge_omega_f2_cut
    *cos(M_PI_2*(a_rho_f2-a_omega_f2));
  double regge_rho_f2_rho_P_cuts=regge_rho_f2_cut*regge_rho_P_cut
    *cos(M_PI_2*(a_rho_P-a_rho_f2));

  // Compute contributions to square of amplitude
  double Pi_R_sq=gRsq/(ReB*ReB+ImB*ImB);
  Pi_R_sq=1;
  double aS_aS_rho_V_and_T_sq
    =Pi_R_sq*gsq_rho_S_gamma*gsq_rho_V_and_T*regge_rho_sq;
  double aS_aS_omega_V_sq=Pi_R_sq*gsq_omega_V*gsq_omega_S_gamma
    *(regge_omega_sq 
      + C_omega_P_cut*regge_omega_omega_P_cut
      + C_omega_f2_cut*regge_omega_omega_f2_cut
      + Csq_omega_P_cut*regge_omega_P_cut_sq
      + Csq_omega_f2_cut*regge_omega_f2_cut_sq
      + 2.*C_omega_f2_C_omega_P*regge_omega_f2_omega_P_cuts);
      
  double aS_aS_omega_V_rho_V_and_T
    =Pi_R_sq*g_rho_S_gamma*g_omega_S_gamma*g_rho_V_and_T*g_omega_V
    *(0.5*regge_omega*regge_rho*cos_rho_omega_sum
      + C_omega_P_cut*regge_rho_omega_P_cut
      + C_omega_f2_cut*regge_rho_omega_f2_cut
      );
  double aS_aS_rho_V_rho_V_and_T=Pi_R_sq*gsq_rho_S_gamma*g_rho_V*g_rho_V_and_T
    *(C_rho_P_cut*regge_rho_rho_P_cut
      + C_rho_f2_cut*regge_rho_rho_f2_cut);
  double aS_aS_rho_V_omega_V
    =Pi_R_sq*g_rho_S_gamma*g_omega_S_gamma*g_rho_V*g_omega_V
    *(C_rho_P_cut*regge_omega_rho_P_cut
      + C_rho_f2_cut*regge_omega_rho_f2_cut
      + 2.*C_rho_P_cut*C_omega_P_cut*regge_rho_P_omega_P_cuts
      + 2.*C_rho_P_cut*C_omega_f2_cut*regge_rho_P_omega_f2_cuts
      + 2.*C_rho_f2_cut*C_omega_P_cut*regge_rho_f2_omega_P_cuts
      + 2.*C_rho_f2_cut*C_omega_f2_cut*regge_rho_f2_omega_f2_cuts);
  double aS_aS_rho_cut_sq=Pi_R_sq*gsq_rho_S_gamma*gsq_rho_V
    *(Csq_rho_P_cut*regge_rho_P_cut_sq
      + Csq_rho_f2_cut*regge_rho_f2_cut_sq
      + 2.*C_rho_P_cut*C_rho_f2_cut*regge_rho_f2_rho_P_cuts);
  double aS_aS=aS_aS_rho_V_and_T_sq + aS_aS_omega_V_sq + aS_aS_rho_V_rho_V_and_T
    + aS_aS_omega_V_rho_V_and_T+aS_aS_rho_V_omega_V + aS_aS_rho_cut_sq;

  double bS_bS=4*Pi_R_sq*gsq_rho_S_gamma*g_rho_T*g_rho_T*regge_rho_sq;

  double aS_bS_rho_T_rho_V_and_T
    =-4.*Pi_R_sq*gsq_rho_S_gamma*g_rho_T*g_rho_V_and_T*regge_rho_sq;
  double aS_bS_rho_cuts=-4.*Pi_R_sq*gsq_rho_S_gamma*g_rho_T*g_rho_V
    *(C_rho_P_cut*regge_rho_rho_P_cut
      + C_rho_f2_cut*regge_rho_rho_f2_cut);
  double aS_bS_omega_cuts
    =-4.*Pi_R_sq*g_rho_S_gamma*g_omega_S_gamma*g_rho_T*g_omega_V
    *(C_omega_P_cut*regge_rho_omega_P_cut
      + C_omega_f2_cut*regge_rho_omega_f2_cut);
  double aS_bS_omega_V=-Pi_R_sq*g_rho_S_gamma*g_omega_S_gamma*g_rho_T*g_omega_V
    *regge_rho*regge_omega*cos_rho_omega_sum;
  double aS_bS=aS_bS_rho_T_rho_V_and_T+aS_bS_rho_cuts+aS_bS_omega_cuts+aS_bS_omega_V;


  // Cut contribution from b1 exchange.  We ignore the pole term.
  double Tb1=-0.5*t*kin1*Pi_R_sq
    *(gsq_omega_S_gamma*gsq_omega_V*(Csq_omega_P_cut*regge_omega_P_cut_sq
				     + Csq_omega_f2_cut*regge_omega_f2_cut_sq
				     + 2.*C_omega_P_cut*C_omega_f2_cut
				     *regge_omega_f2_omega_P_cuts)
      +gsq_rho_S_gamma*gsq_rho_V*(Csq_rho_P_cut*regge_rho_P_cut_sq
				  + Csq_rho_f2_cut*regge_rho_f2_cut_sq
				  + 2.*C_rho_P_cut*C_rho_f2_cut
				  *regge_rho_f2_rho_P_cuts)
      +g_rho_V*g_omega_V*g_rho_S_gamma*g_omega_S_gamma
      *(2.*C_rho_P_cut*C_omega_P_cut*regge_rho_P_omega_P_cuts
	+ 2.*C_rho_f2_cut*C_omega_f2_cut*regge_rho_f2_omega_f2_cuts 
	+ 2.*C_rho_P_cut*C_omega_f2_cut*regge_rho_P_omega_f2_cuts
	+ 2.*C_rho_f2_cut*C_omega_P_cut*regge_rho_f2_omega_P_cuts));
  	
  double T=aS_aS*kin_aS_aS+aS_bS*kin_aS_bS+bS_bS*kin_bS_bS + Tb1;

  // Compute cross section
  double m1_plus_m2=m1+m2;
  double m1_minus_m2=m1-m2;
  double k=sqrt((ms_sq-m1_plus_m2*m1_plus_m2)*(ms_sq-m1_minus_m2*m1_minus_m2))
    /(2.*sqrt(ms_sq));
  double hbarc_sq=389.; // Convert to micro-barns
  double xsec=-hbarc_sq*k*T/(256.*M_PI*M_PI*M_PI*M_PI*mp_sq_minus_s_sq);

  return(xsec);
				 

}

// Put particle data into hddm format and output to file
void WriteEvent(unsigned int eventNumber,TLorentzVector &beam, float vert[3],
		vector<Particle_t>&particle_types,
		vector<TLorentzVector>&particle_vectors, s_iostream_t *file){  
   s_PhysicsEvents_t* pes;
   s_Reactions_t* rs;
   s_Target_t* ta;
   s_Beam_t* be;
   s_Vertices_t* vs;
   s_Origin_t* origin;
   s_Products_t* ps;
   s_HDDM_t *thisOutputEvent = make_s_HDDM();
   thisOutputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
   pes->mult = 1;
   pes->in[0].runNo = runNo;
   pes->in[0].eventNo = eventNumber;
   pes->in[0].reactions = rs = make_s_Reactions(1);
   rs->mult = 1;
   // Beam 
   rs->in[0].beam = be = make_s_Beam();
   be->type = Gamma;
   be->properties = make_s_Properties();
   be->properties->charge = ParticleCharge(be->type);
   be->properties->mass = ParticleMass(be->type);
   be->momentum = make_s_Momentum();
   be->momentum->px = 0.;
   be->momentum->py = 0.;
   be->momentum->pz = beam.Pz();
   be->momentum->E  = beam.E();
   // Target
   rs->in[0].target = ta = make_s_Target();
   ta->type = Proton;
   ta->properties = make_s_Properties();
   ta->properties->charge = ParticleCharge(ta->type);
   ta->properties->mass = ParticleMass(ta->type);
   ta->momentum = make_s_Momentum();
   ta->momentum->px = 0.;
   ta->momentum->py = 0.;
   ta->momentum->pz = 0.;
   ta->momentum->E  = ParticleMass(ta->type);
   // Primary vertex 
   rs->in[0].vertices = vs = make_s_Vertices(1);
   vs->mult = 1;
   vs->in[0].origin = origin = make_s_Origin();
   vs->in[0].products = ps = make_s_Products(particle_vectors.size());
   ps->mult = 0;
   origin->t = 0.0;
   origin->vx = vert[0];
   origin->vy = vert[1];
   origin->vz = vert[2];
   // Final state particles
   for (unsigned int i=0;i<particle_vectors.size();i++,ps->mult++){
     Particle_t my_particle=particle_types[i];
     ps->in[ps->mult].type = my_particle;
     ps->in[ps->mult].pdgtype = PDGtype(my_particle);
     ps->in[ps->mult].id = i+1; /* unique value for this particle within the event */
     ps->in[ps->mult].parentid = 0;  /* All internally generated particles have no parent */
     ps->in[ps->mult].mech = 0; // ???     
     ps->in[ps->mult].momentum = make_s_Momentum();
     ps->in[ps->mult].momentum->px = particle_vectors[i].Px();
     ps->in[ps->mult].momentum->py = particle_vectors[i].Py();
     ps->in[ps->mult].momentum->pz = particle_vectors[i].Pz();
     ps->in[ps->mult].momentum->E  = particle_vectors[i].E();
   }
   flush_s_HDDM(thisOutputEvent,file);
}

// Create some diagnostic histograms
void CreateHistograms(){

  thrown_t=new TH1D("thrown_t","Thrown t distribution",1000,-0.99,0.01);
  thrown_t->SetXTitle("t [GeV^{2}]");
  thrown_dalitzZ=new TH1D("thrown_dalitzZ","thrown dalitz Z",110,-0.05,1.05);
  thrown_Egamma=new TH1D("thrown_Egamma","Thrown E_{#gamma} distribution",
			       1000,0,12.);
  thrown_Egamma->SetTitle("E_{#gamma} [GeV]");
  thrown_dalitzXY=new TH2D("thrown_dalitzXY","Dalitz distribution Y vs X",100,-1.,1.,100,-1.,1);
  
  thrown_theta_vs_p=new TH2D("thrown_theta_vs_p","Proton #theta_{LAB} vs. p",
			       200,0,2.,180,0.,90.);
  thrown_theta_vs_p->SetXTitle("p [GeV/c]");
  thrown_theta_vs_p->SetYTitle("#theta [degrees]");
  
  cobrems_vs_E=new TH1D("cobrems_vs_E","Coherent bremsstrahlung spectrum",
			1000,Emin,Emax);

}


// Create a graph of the cross section dsigma/dt as a function of -t
void GraphCrossSection(double m1,double m2){
  // beam energy in lab
  double Egamma=Emin;
  TLorentzVector beam(0,0,Egamma,Egamma);
  TLorentzVector target(0,0,0,m_p);

  // CM energy
  double s=m_p*(m_p+2.*Egamma);
  double Ecm=sqrt(s);

  // Momenta of incoming photon and outgoing S and proton in cm frame
  double p_gamma=(s-m_p_sq)/(2.*Ecm);
  double E_S=(s+m_S_sq_R-m_p_sq)/(2.*Ecm);
  double p_S=sqrt(E_S*E_S-m_S_sq_R);
  
  // Momentum transfer t
  double p_diff=p_gamma-p_S;
  double t0=m_S_sq_R*m_S_sq_R/(4.*s)-p_diff*p_diff;
 
  // Parameters for integration over line shape	  
  double m1_plus_m2=m1+m2;
  double m_max=m_p*(sqrt(1.+2.*Egamma/m_p)-1.);
  double dmrange=m_max-m1_plus_m2;
  double dm=dmrange/1000.;
  double gR=0.,ImB=0,ReB=0; // resonance parameters
  bool got_pipi=(fabs(m1-m2)>0.01)?false:true;
  if (got_pipi){ // f0(980)
    gR=1.705/(2.*M_PI);
    m_S_sq_R=0.9783*0.9783;
    gsq_rho_S_gamma=0.239; // GeV^-2
    gsq_omega_S_gamma=0.02656;
  }
  else{ // a0(980)
    gR=2.82/(2.*M_PI);
    m_S_sq_R=0.9825*0.9825;	 
    gsq_rho_S_gamma=0.02537;
    gsq_omega_S_gamma=0.2283;
  } 
  double m1_minus_m2=m1-m2;
  double k=sqrt((m_S_sq_R-m1_plus_m2*m1_plus_m2)*(m_S_sq_R-m1_minus_m2*m1_minus_m2))
    /(2.*sqrt(m_S_sq_R));

  GetResonanceParameters(m1,m2,m_S_sq_R,ReB,ImB);
  // factor to (eventually ) get to dsigma/dt from dsigma/dt/dM/dOmega
  double scale_factor=16.*M_PI*M_PI*M_PI*(ReB*ReB+ImB*ImB)/(gR*gR*k);
  // Integral over Breit-Wigner
  double bw=0;
  for (unsigned int n=0;n<1000;n++){
    double mass=m1_plus_m2+dm*double(n);
    double m_S_sq=mass*mass;
    double rho_m1m2=sqrt((1.-m1_plus_m2*m1_plus_m2/m_S_sq)*(1-m1_minus_m2*m1_minus_m2/m_S_sq));

    GetResonanceParameters(m1,m2,m_S_sq,ReB,ImB);
    bw+=gR*gR*rho_m1m2/(ReB*ReB+ImB*ImB)*dm;
  }
  
  // differential cross section
  double sum=0.;
  double t_old=t0;
  double t_array[1000];
  double xsec_array[1000];
  for (unsigned int k=0;k<1000;k++){
    double theta_cm=M_PI*double(k)/1000.;
    double sin_theta_over_2=sin(0.5*theta_cm);
    double t=t0-4.*p_gamma*p_S*sin_theta_over_2*sin_theta_over_2;
    double xsec=scale_factor*bw*CrossSection(m1,m2,m_S_sq_R,s,t,gR,ReB,ImB);

    t_array[k]=-t;
    xsec_array[k]=xsec;

    sum-=xsec*(t-t_old);
    t_old=t;
  }
  TGraph *Gxsec=new TGraph(1000,t_array,xsec_array);
  Gxsec->Write("Cross section");
 
  cout << "Total cross section at " << Egamma << " GeV = "<< sum 
       << " micro-barns"<<endl;
}


//-----------
// main
//-----------
int main(int narg, char *argv[])
{  
  ParseCommandLineArguments(narg, argv);

  // open ROOT file
  string rootfilename="scalar_gen.root";
  TFile *rootfile=new TFile(rootfilename.c_str(),"RECREATE",
			    "Produced by genScalarRegge");

  // open HDDM file
  s_iostream_t *file = init_s_HDDM(output_file_name);
 
 
  // Initialize random number generator
  TRandom3 *myrand=new TRandom3(0);// If seed is 0, the seed is automatically computed via a TUUID object, according to TRandom3 documentation

  // Fixed target
  TLorentzVector target(0.,0.,0.,m_p);

  //----------------------------------------------------------------------------
  // Get production (Egamma range) and decay parameters from input file
  //----------------------------------------------------------------------------

  // Start reading the input file 
  ifstream infile(input_file_name);
  if (!infile.is_open()){
    cerr << "Input file missing! Exiting..." <<endl;
    exit(-1);
  } 

  // Get photon energy range
  string comment_line;
  getline(infile,comment_line);
  infile >> Emin;
  infile >> Emax;
  infile.ignore(); // ignore the '\n' at the end of this line
  // Set sensible minimum energy
  if (Emin<m_S){ 
    Emin=m_S;
    cout << "Warning:  Setting minimum beam energy to " << Emin << " [GeV]" 
	 <<endl;
  }
  cout << "Photon energy min, max [Gev] = "<< Emin <<","<<Emax <<endl;

  // Get coherent peak and collimator diameter
  getline(infile,comment_line);
  float Epeak=9.0,collDiam=0.0034;
  float Ee=12.0;
  infile >> Ee;
  infile >> Epeak;
  infile >> collDiam;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "Electron beam energy = " << Ee << " GeV, Coherent peak = " 
       << Epeak <<" GeV, collimator diameter = " 
       <<collDiam << " m" << endl;
  
  // Set number of decay particles
  int num_decay_particles=2;
  // Set up vectors of particle ids and 4-vectors
  int last_index=num_decay_particles;
  int num_final_state_particles=num_decay_particles+1;
  vector<TLorentzVector>particle_vectors(num_final_state_particles);
  vector<Particle_t>particle_types(num_final_state_particles);
  double *decay_masses =new double[num_decay_particles];
  particle_types[last_index]=Proton;

  // GEANT ids of decay particles
  getline(infile,comment_line);
  cout << "Particle id's of decay particles =";
  for (int k=0;k<num_decay_particles;k++){
    int ipart;
    infile >> ipart;
    cout << " " << ipart; 
    particle_types[k]=(Particle_t)ipart;
    decay_masses[k]=ParticleMass(particle_types[k]);
  }
  infile.ignore(); // ignore the '\n' at the end of this line
  cout << endl;

  // processes to simulate
  int num_processes=5;
  getline(infile,comment_line);
  int *generate=new int[num_processes];
  for (int k=0;k<num_processes;k++){
    infile >> generate[k];
  }
  infile.close();
 
  //----------------------------------------------------------------------------
  // Setup coherent bremsstrahlung generator
  //----------------------------------------------------------------------------
  float radColDist=76.0;// meters
  //float colDiam=0.0034; // meters
  int doPolFlux=0;  // want total flux (1 for polarized flux)
  float emitmr=10.e-9; // electron beam emittance
  float radt=20.e-6; // radiator thickness in m
  cobrems_(&Ee,&Epeak,&emitmr,&radt,&radColDist,&collDiam,&doPolFlux);
  
  // Create some diagonistic histograms
  CreateHistograms();
  // Make a TGraph of the cross section at a fixed beam energy
  GraphCrossSection(decay_masses[0],decay_masses[1]);
  
  // Fill histogram of coherent bremmsstrahlung distribution 
  for (int i=1;i<=1000;i++){
    float x=float(cobrems_vs_E->GetBinCenter(i)/Ee);
    float y=0;
    if (Epeak<Emin) y=dnidx_(&x);
    else y=dntdx_(&x);
    cobrems_vs_E->Fill(Ee*double(x),double(y));
  }

  // Set up some variables for cross section calculation
  // masses of decay particles
  double m1=decay_masses[0];
  double m2=decay_masses[1];
  bool got_pipi=(fabs(m1-m2)>0.01)?false:true;

  //----------------------------------------------------------------------------
  // Event generation loop
  //----------------------------------------------------------------------------
  for (int i=1;i<=Nevents;i++){
    // photon beam
    double Egamma=0.;
    TLorentzVector beam;

    // Maximum value for cross section 
    double xsec_max=0.3;
    double xsec=0.,xsec_test=0.;

    // Polar angle in center of mass frame
    double theta_cm=0.;

    // Eta momentum in cm
    double p_S=0.;

    // Mass squared of resonance
    double m_S_sq=m_S_sq_R;

    // Transfer 4-momentum;
    double t=0.;

    // vertex position at target
    float vert[4]={0.,0.,0.,0.};

    // masses of decay particles

    // use the rejection method to produce S's based on the cross section
    do{
      // First generate a beam photon using bremsstrahlung spectrum
      Egamma = cobrems_vs_E->GetRandom();

      // CM energy
      double s=m_p*(m_p+2.*Egamma);
      double Ecm=sqrt(s);

      // Momenta of incoming photon and outgoing S and proton in cm frame
      double p_gamma=(s-m_p_sq)/(2.*Ecm);

      // Mass of two-meson system     
      double m1_plus_m2=m1+m2;
      double m_max=m_p*(sqrt(1.+2.*Egamma/m_p)-1.);
      double m_S=m1_plus_m2+myrand->Uniform(m_max-m1_plus_m2);
      m_S_sq=m_S*m_S;

      // Momentum and energy of two-meson system
      double E_S=(s+m_S_sq-m_p_sq)/(2.*Ecm);
      p_S=sqrt(E_S*E_S-m_S_sq);
    
      // Momentum transfer t
      double p_diff=p_gamma-p_S;
      double t0=m_S_sq*m_S_sq/(4.*s)-p_diff*p_diff;
      double sin_theta_over_2=0.;
      t=t0;
      
      // Generate theta with a uniform distribution and compute the cross 
      // section at this value
      theta_cm=myrand->Uniform(M_PI);
      // compute t
      sin_theta_over_2=sin(0.5*theta_cm);
      t=t0-4.*p_gamma*p_S*sin_theta_over_2*sin_theta_over_2;

      // Generate phi using uniform distribution
      double phi_cm=myrand->Uniform(2.*M_PI);

      // beam 4-vector (ignoring px and py, which are extremely small)
      beam.SetXYZT(0.,0.,Egamma,Egamma);
   
      // Velocity of the cm frame with respect to the lab frame
      TVector3 v_cm=(1./(Egamma+m_p))*beam.Vect();
      // Four-momentum of the S in the CM frame
      double pt=p_S*sin(theta_cm);
      TLorentzVector S4(pt*cos(phi_cm),pt*sin(phi_cm),p_S*cos(theta_cm),
			sqrt(p_S*p_S+m_S_sq));
      // S4.Print();
      
      //Boost the S 4-momentum into the lab
      S4.Boost(v_cm);
      // S4.Print();
      
      
      // Compute the 4-momentum for the recoil proton
      TLorentzVector proton4=beam+target-S4;
      
      // Generate 3-body decay of S according to phase space
      TGenPhaseSpace phase_space;
      phase_space.SetDecay(S4,num_decay_particles,decay_masses);
      double weight=0.,rand_weight=1.;
      do{
	weight=phase_space.Generate();
	rand_weight=myrand->Uniform(1.);
      }
      while (rand_weight>weight);
      
      // Gather the particles in the reaction and write out event in hddm format
      particle_vectors[last_index]=proton4;
      for (int j=0;j<num_decay_particles;j++){
	particle_vectors[j]=*phase_space.GetDecay(j);
      }

      //Resonance parameters 
      double ReB=0.,ImB=0,gR=0.;

      // Cross section
      xsec=0.;
      
      // f0(600)
      if (got_pipi && generate[0]){
	/*
	double m_Sigma=0.65;
	m_S_sq_R=m_Sigma*m_Sigma; 
	gsq_rho_S_gamma=0.01*pow(0.7685*0.785-m_S_sq_R,-3); // GeV^-2
	gsq_omega_S_gamma=8.97e-4*pow(0.78265*0.78265-m_S_sq_R,-3);
	width=0.6;
	double bw=ResonanceShape(m_S_sq,width,m1,m2,0); 
	xsec+=CrossSection(s,t,m_S_sq);
	*/
      }
      // f0(980)/a0(980)
      if (generate[1]){
	GetResonanceParameters(m1,m2,m_S_sq,ReB,ImB);
	if (got_pipi){ // f0(980)
	  gR=1.705/(2.*M_PI);
	  m_S_sq_R=0.9783*0.9783;
	  gsq_rho_S_gamma=0.239; // GeV^-2
	  gsq_omega_S_gamma=0.02656;
	}
	else{ // a0(980)
	  gR=2.82/(2.*M_PI);
	  m_S_sq_R=0.9825*0.9825;	 
	  gsq_rho_S_gamma=0.02537;
	  gsq_omega_S_gamma=0.2283;
	} 
	xsec+=CrossSection(m1,m2,m_S_sq,s,t,gR,ReB,ImB);
      }
        
      if (generate[4]){ // background
	xsec+=BackGroundCrossSection(beam,particle_types,particle_vectors);

	if (generate[1]){ // interference with resonant signal

	}
      }
      


      // Generate a test value for the cross section
      xsec_test=myrand->Uniform(xsec_max);
    }
    while (xsec_test>xsec);

    // Other diagnostic histograms
    thrown_t->Fill(t);
    thrown_Egamma->Fill(Egamma);
    thrown_theta_vs_p->Fill(particle_vectors[last_index].P(),
			    180./M_PI*particle_vectors[last_index].Theta());

     // Randomly generate z position in target
    vert[2]=zmin+myrand->Uniform(zmax-zmin);

    WriteEvent(i,beam,vert,particle_types,particle_vectors,file);
    
    if ((i%100)==0) cout << 100.*double(i)/double(Nevents) << "\% done" << endl;
  }


  // Write histograms and close root file
  rootfile->Write();
  rootfile->Close();

  // Close HDDM file
  close_s_HDDM(file);
  cout<<endl<<"Closed HDDM file"<<endl;
  cout<<" "<<Nevents<<" event written to "<<output_file_name<<endl;

  // Cleanup
  delete []decay_masses;

  return 0;
}


