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

// Breit-Wigner distribution for resonance shape
double ResonanceShape(double m_S_sq,double width,double m1,double m2,
		      int shape_type){
  double bw=0.;
  double m1_plus_m2=m1+m2;
  double m1_minus_m2=m1-m2;  
  bool got_pipi=(fabs(m1_minus_m2)>0.01)?false:true;

  double M2diff=m_S_sq_R-m_S_sq;
  double rho_m1m2
    =sqrt((1.-m1_plus_m2*m1_plus_m2/m_S_sq)*(1-m1_minus_m2*m1_minus_m2/m_S_sq));

  switch(shape_type){
  case 0:  // Relativistic Breit-wigner with no channel-opening threshold effects
    {
      double fac=m_S_sq_R*sqrt(m_S_sq_R)*width
	/sqrt((1.-m1_plus_m2*m1_plus_m2/m_S_sq_R)
	      *(1-m1_minus_m2*m1_minus_m2)/m_S_sq_R);
      bw=fac*rho_m1m2/(M2diff*M2diff+m_S_sq*width*width);
    }
    break;
  case 1: // "No structure" model for a0(980)/f0(980)
    {
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

      // Real and imaginary parts of BW amplitude
      double ReDm=M2diff;
      if (M_S<2.*MKPlus){
	ReDm+=gKsq/(32.*M_PI)*sqrt(-rhoKPlussq);
      }
      if (M_S<2.*MK0){
	ReDm+=gKsq/(32.*M_PI)*sqrt(-rhoK0sq);
      }
      double ImDm=g_m1m2_sq/(16*M_PI)*rho_m1m2;
      if (M_S>2.*MKPlus){
	ImDm+=gKsq/(32.*M_PI)*sqrt(rhoKPlussq);
      }
      if (M_S>2.*MK0){
	ImDm+=gKsq/(32.*M_PI)*sqrt(rhoK0sq);
      }
      bw=g_m1m2_sq*rho_m1m2/(8*M_PI*M_PI*(ReDm*ReDm+ImDm*ImDm));
    }
    break;
  default:
    break;
  }
  return bw;
}


// Cross section dsigma/dt derived from 
double CrossSection(double s,double t,double ms_sq){
  // Coupling constants 
  double g_omega_V=15.;
  double gsq_omega_V=g_omega_V*g_omega_V;
  double g_rho_V=3.4;
  double gsq_rho_V=g_rho_V*g_rho_V;
  double g_rho_T=11.0; // GeV^-1
  double g_rho_V_and_T=g_rho_V+2.*m_p*g_rho_T;
  double g_rho_S_gamma=sqrt(gsq_rho_S_gamma);
  double g_omega_S_gamma=sqrt(gsq_omega_S_gamma);

  // s scale for regge trajectories
  double s0=1.;

  // Kinematical boundaries
  double mp_sq_minus_s=m_p_sq-s;
  double mp_sq_plus_s=m_p_sq+s;
  double mp_sq_minus_s_sq=mp_sq_minus_s*mp_sq_minus_s;
  double temp1=ms_sq*mp_sq_plus_s-mp_sq_minus_s_sq;
  double temp2=mp_sq_minus_s*sqrt(mp_sq_minus_s_sq-2.*ms_sq*mp_sq_plus_s
				  +ms_sq*ms_sq);
  double t1=(temp1+temp2)/(2.*s);
  double t2=(temp1-temp2)/(2.*s);

  // Kinematic factors for cross section
  double t_minus_ms_sq=t-ms_sq;
  double xfac1=s*(t-t1)*(t-t2);
  double xfac2=xfac1+2.*t*t_minus_ms_sq*t_minus_ms_sq;

  // Regge trajectory for omega
  double a_omega=0.44+0.9*t;
  double a_omega_prime=0.9;
  double regge_omega=pow(s/s0,a_omega-1.)*M_PI*a_omega_prime/(sin(M_PI*a_omega)*TMath::Gamma(a_omega)); // excluding phase factor

  // Regge cuts for omega
  double a_omega_P=0.52+0.196*t; // Pomeron
  double a_omega_f2=0.112+0.428*t;
  double dc=2.;
  double regge_omega_P_cut=exp(dc*t)*pow(s/s0,a_omega_P-1.);
  double regge_omega_f2_cut=exp(dc*t)*pow(s/s0,a_omega_f2-1.);
  double C_omega_P_cut=0.1;
  double C_omega_f2_cut=0.1;

  // omega amplitude squared
  double M_omega_sq=0.5*gsq_omega_S_gamma*gsq_omega_V*xfac2
    *(regge_omega*regge_omega*0.5*(1.-cos(M_PI*a_omega))
      +C_omega_P_cut*C_omega_P_cut*regge_omega_P_cut*regge_omega_P_cut
      +C_omega_f2_cut*C_omega_f2_cut*regge_omega_f2_cut*regge_omega_f2_cut
      +(2.*C_omega_f2_cut*C_omega_P_cut*regge_omega_f2_cut*regge_omega_P_cut
	*cos(M_PI_2*(a_omega_f2-a_omega_P)))
      +(C_omega_f2_cut*regge_omega*regge_omega_f2_cut
	*(cos(M_PI*(a_omega+0.5*a_omega_f2))-cos(M_PI_2*a_omega_f2)))
      +(C_omega_P_cut*regge_omega*regge_omega_P_cut
	*(cos(M_PI*(a_omega+0.5*a_omega_P))-cos(M_PI_2*a_omega_P)))
      );

  
  // Regge trajectory for rho
  double a_rho=0.55+0.8*t;
  double a_rho_prime=0.8;
  double regge_rho=pow(s/s0,a_rho-1.)*M_PI*a_rho_prime/(sin(M_PI*a_rho)*TMath::Gamma(a_rho)); // excluding phase factor

  // Regge cuts for rho
  double a_rho_P=0.64+0.16*t; // Pomeron
  double a_rho_f2=0.222+0.404*t;
  double regge_rho_P_cut=exp(dc*t)*pow(s/s0,a_rho_P-1.);
  double regge_rho_f2_cut=exp(dc*t)*pow(s/s0,a_rho_f2-1.);
  double C_rho_P_cut=0.1;
  double C_rho_f2_cut=0.1;

  // Rho amplitude squared
  double rho_cut_interference=C_rho_f2_cut*regge_rho*regge_rho_f2_cut
    *(cos(M_PI*(a_rho+0.5*a_rho_f2))-cos(M_PI_2*a_rho_f2))
    +C_rho_P_cut*regge_rho*regge_rho_P_cut
    *(cos(M_PI*(a_rho+0.5*a_rho_P))-cos(M_PI_2*a_rho_P));
  double regge_rho_amp_sq=regge_rho*regge_rho*0.5*(1.-cos(M_PI*a_rho));      
  double M_rho_sq
    =0.125*gsq_rho_S_gamma
    *(4*xfac2
      *(g_rho_V_and_T*g_rho_V_and_T*regge_rho_amp_sq
	+g_rho_V*g_rho_V_and_T*rho_cut_interference
	+gsq_rho_V*(C_rho_P_cut*C_rho_P_cut*regge_rho_P_cut*regge_rho_P_cut
			  +C_rho_f2_cut*C_rho_f2_cut*regge_rho_f2_cut*regge_rho_f2_cut
			  +(2.*C_rho_f2_cut*C_rho_P_cut*regge_rho_f2_cut*regge_rho_P_cut
			    *cos(M_PI_2*(a_rho_f2-a_rho_P)))))
      -4.*m_p*xfac1*(g_rho_T*g_rho_V_and_T*regge_rho_amp_sq
		     +0.5*g_rho_T*g_rho_T*rho_cut_interference)
      +g_rho_T*g_rho_T*xfac1*(4.*m_p_sq-t)*regge_rho_amp_sq);


  // rho-omega interference
  double regge_rho_omega_cut_P=regge_rho*regge_omega_P_cut
    *(cos(M_PI*(a_rho-0.5*a_omega_P))-cos(0.5*M_PI*a_omega_P)); 
  double regge_rho_omega_cut_f2=regge_rho*regge_omega_f2_cut
    *(cos(M_PI*(a_rho-0.5*a_omega_f2))-cos(0.5*M_PI*a_omega_f2)); 
  double regge_omega_rho_cut_P=regge_omega*regge_rho_P_cut
    *(cos(M_PI*(a_omega-0.5*a_rho_P))-cos(0.5*M_PI*a_rho_P)); 
  double regge_omega_rho_cut_f2=regge_omega*regge_rho_f2_cut
    *(cos(M_PI*(a_omega-0.5*a_rho_f2))-cos(0.5*M_PI*a_rho_f2));
  double regge_rho_cut_f2_omega_cut_f2=regge_rho_f2_cut*regge_omega_f2_cut
    *cos(M_PI_2*(a_rho_f2-a_omega_f2)); 
  double regge_rho_cut_P_omega_cut_f2=regge_rho_P_cut*regge_omega_f2_cut
    *cos(M_PI_2*(a_rho_P-a_omega_f2)); 
  double regge_rho_cut_f2_omega_cut_P=regge_rho_f2_cut*regge_omega_P_cut
    *cos(M_PI_2*(a_rho_f2-a_omega_P)); 
  double regge_rho_cut_P_omega_cut_P=regge_rho_P_cut*regge_omega_P_cut
    *cos(M_PI_2*(a_rho_P-a_omega_P));
  double regge_rho_omega=0.25*regge_rho*regge_omega*(1.+cos(M_PI*(a_omega-a_rho))
						     -cos(M_PI*a_omega)
						     -cos(M_PI*a_rho)
						     );
  double regge_rho_omega_sum=regge_rho_omega+C_omega_P_cut*regge_rho_omega_cut_P
    +C_omega_f2_cut*regge_rho_omega_cut_f2; 
  double M_rho_omega_sq=g_omega_V*g_rho_S_gamma*g_omega_S_gamma
    *(xfac2*(g_rho_V_and_T*regge_rho_omega_sum			 
	     +g_rho_V*(C_rho_P_cut*regge_omega_rho_cut_P
		       +C_rho_f2_cut*regge_omega_rho_cut_f2
		       +C_rho_P_cut*C_omega_P_cut*regge_rho_cut_P_omega_cut_P
		       +C_rho_P_cut*C_omega_f2_cut*regge_rho_cut_P_omega_cut_f2
		       +C_rho_f2_cut*C_omega_P_cut*regge_rho_cut_f2_omega_cut_P
		       +C_rho_f2_cut*C_omega_f2_cut*regge_rho_cut_f2_omega_cut_f2
		       )
	     )
      -2.*g_rho_T*xfac1*m_p*regge_rho_omega_sum);

  // Cut contribution from b1 exchange.  We ignore the pole term  
  double regge_omega_cut_P_omega_cut_f2=regge_omega_P_cut*regge_omega_f2_cut
    *cos(M_PI_2*(a_omega_P-a_omega_f2));  
  double regge_rho_cut_P_rho_cut_f2=regge_rho_P_cut*regge_rho_f2_cut
    *cos(M_PI_2*(a_rho_P-a_rho_f2)); 
  double M_b1_sq=-0.5*t*xfac1/m_S_sq_R
    *(gsq_omega_S_gamma*gsq_omega_V*(C_omega_P_cut*C_omega_P_cut
				     *regge_omega_P_cut*regge_omega_P_cut
				     +C_omega_f2_cut*C_omega_f2_cut
				     *regge_omega_f2_cut*regge_omega_f2_cut
				     +2.*C_omega_f2_cut*C_omega_f2_cut
				     *regge_omega_cut_P_omega_cut_f2)
      +gsq_rho_S_gamma*gsq_rho_V*(C_rho_P_cut*C_rho_P_cut
				     *regge_rho_P_cut*regge_rho_P_cut
				     +C_rho_f2_cut*C_rho_f2_cut
				     *regge_rho_f2_cut*regge_rho_f2_cut
				     +2.*C_rho_f2_cut*C_rho_f2_cut
				     *regge_rho_cut_P_rho_cut_f2)
      +g_rho_V*g_omega_V*g_rho_S_gamma*g_omega_S_gamma
      *(2.*C_rho_P_cut*C_omega_P_cut*regge_rho_cut_P_omega_cut_P
	+2.*C_rho_f2_cut*C_omega_f2_cut*regge_rho_cut_f2_omega_cut_f2
	+2.*C_rho_P_cut*C_omega_f2_cut*regge_rho_cut_P_omega_cut_f2
	+2.*C_rho_f2_cut*C_omega_P_cut*regge_rho_cut_f2_omega_cut_P));
 
  double M_sq =M_omega_sq+M_rho_sq+M_rho_omega_sq+M_b1_sq;

  double hbarc_sq=389.; // Convert to micro-barns
  double dsigma_dt=-hbarc_sq*M_sq/(16.*M_PI*mp_sq_minus_s_sq);

  return(dsigma_dt);
				 

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
void GraphCrossSection(double m1, double m2){
  // beam energy in lab
  double Egamma=Emin;

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
 
  // Integration over line shape	  
  double m1_plus_m2=m1+m2;
  double m_max=m_p*(sqrt(1.+2.*Egamma/m_p)-1.);
  double bw=0;
  double dmrange=m_max-m1_plus_m2;
  double dm=dmrange/1000.;
  for (int i=0;i<1000;i++){
    double m=m1_plus_m2+dm*double(i);
    bw+=ResonanceShape(m*m,width,m1,m2,1)*dm;
  }

  // differential cross section
  double sum=0.;
  double t_old=t0;
  double t_array[10000];
  double xsec_array[10000];
  for (unsigned int k=0;k<10000;k++){
    double theta_cm=M_PI*double(k)/10000.;
    double sin_theta_over_2=sin(0.5*theta_cm);
    double t=t0-4.*p_gamma*p_S*sin_theta_over_2*sin_theta_over_2;
    double xsec=2.*m_S_sq_R/M_PI*bw*CrossSection(s,t,m_S_sq_R);
    
    t_array[k]=-t;
    xsec_array[k]=xsec;

    sum-=xsec*(t-t_old);
    t_old=t;
  }
  TGraph *Gxsec=new TGraph(10000,t_array,xsec_array);
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
  
  // Create some diagonistic histographs
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
    double Egamma=0.;
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

      // Cross section
      xsec=0.;
      
      // f0(600)
      if (got_pipi && generate[0]){
	m_S_sq_R=0.6*0.6;
	gsq_rho_S_gamma=0.8155; // GeV^-2
	gsq_omega_S_gamma=0.05567;
	width=0.6;
	double bw=ResonanceShape(m_S_sq,width,m1,m2,0); 
	xsec+=2.*m_S_sq_R/M_PI*bw*CrossSection(s,t,m_S_sq);
      }
      // f0(980)/a0(980)
      if (generate[1]){
	if (got_pipi){ // f0(980)
	  m_S_sq_R=0.9783*0.9783;
	  gsq_rho_S_gamma=0.239; // GeV^-2
	  gsq_omega_S_gamma=0.02656;
	}
	else{ // a0(980)
	  m_S_sq_R=0.9825*0.9825;	 
	  gsq_rho_S_gamma=0.02537;
	  gsq_omega_S_gamma=0.2283;
	}
	double bw=ResonanceShape(m_S_sq,width,m1,m2,1); 
	xsec+=2.*m_S_sq_R/M_PI*bw*CrossSection(s,t,m_S_sq);
      }

      // Generate a test value for the cross section
      xsec_test=myrand->Uniform(xsec_max);
    }
    while (xsec_test>xsec);

    // Generate phi using uniform distribution
    double phi_cm=myrand->Uniform(2.*M_PI);

    // beam 4-vector (ignoring px and py, which are extremely small)
    TLorentzVector beam(0.,0.,Egamma,Egamma);
    thrown_Egamma->Fill(Egamma);

    // Velocity of the cm frame with respect to the lab frame
    TVector3 v_cm=(1./(Egamma+m_p))*beam.Vect();
    // Four-moementum of the S in the CM frame
    double pt=p_S*sin(theta_cm);
    TLorentzVector S4(pt*cos(phi_cm),pt*sin(phi_cm),p_S*cos(theta_cm),
			sqrt(p_S*p_S+m_S_sq));
    // S4.Print();

    //Boost the S 4-momentum into the lab
    S4.Boost(v_cm);
    // S4.Print();

  
    // Compute the 4-momentum for the recoil proton
    TLorentzVector proton4=beam+target-S4;

    //proton4.Print();
    thrown_theta_vs_p->Fill(proton4.P(),180./M_PI*proton4.Theta());

    // Generate 3-body decay of S according to phase space
    TGenPhaseSpace phase_space;
    phase_space.SetDecay(S4,num_decay_particles,decay_masses);
    double weight=0.,rand_weight=1.;
    do{
      weight=phase_space.Generate();
      rand_weight=myrand->Uniform(1.);
    }
    while (rand_weight>weight);

    // Other diagnostic histograms
    thrown_t->Fill(t);

    // Randomly generate z position in target
    vert[2]=zmin+myrand->Uniform(zmax-zmin);

    // Gather the particles in the reaction and write out event in hddm format
    particle_vectors[last_index]=proton4;
    for (int j=0;j<num_decay_particles;j++){
      particle_vectors[j]=*phase_space.GetDecay(j);
    }
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


