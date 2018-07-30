// Main program for generating eta events. 

#if HAVE_AMPTOOLS_MCGEN

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

#include <AMPTOOLS_MCGEN/CobremsGeneration.hh>

// Masses
const double m_p=0.93827; // GeV
const double m_p_sq=m_p*m_p;
double m_eta=0.54775; // GeV
double m_eta_sq=m_eta*m_eta;
// Width
double width=0.;
// Coupling constants 
double g_rho_eta_gamma=0.81;
double g_omega_eta_gamma=0.29;
double g_eta_gamma_gamma=0.0429;
double g_phi_eta_gamma=0.38;

double Emin=3.,Emax=12.0; // GeV
double zmin=50.0,zmax=80.0; // cm, target extent
int Nevents=10000;
int runNo=10000;
bool debug=false;

// Diagnostic histograms
TH1D *thrown_t;
TH1D *thrown_dalitzZ;
TH1D *thrown_Egamma;
TH2D *thrown_dalitzXY;  
TH2D *thrown_theta_vs_p;
TH1D *cobrems_vs_E;

char input_file_name[250]="eta548.in";
char output_file_name[250]="eta_gen.hddm";

void Usage(void){
  printf("genEtaRegge: generator for eta production based on Regge trajectory formalism.\n");
  printf(" Usage:  genEtaRegge <options>\n");
  printf("   Options:  -N<number of events> (number of events to generate)\n");
  printf("             -O<output.hddm>   (default: eta_gen.hddm)\n");
  printf("             -I<input.in>      (default: eta548.in)\n");
  printf("             -R<run number>    (default: 10000)\n");
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

// Cross section dsigma/dt derived from Laget(2005)
double CrossSection(double s,double t,double p_gamma,double p_eta,double theta){
  // Coupling constants 
  double c_rho_p_p=0.92/137.;
  double c_omega_p_p=6.44/137.;
  double c_gamma_p_p=1./(137.*137.);
  double c_phi_p_p=0.72/137.;

  double kappa_rho=6.1;
  double kappa_gamma=1.79;
  //double kappa_phi=0.;

  // Angular quantities
  double sintheta=sin(theta);
  double costheta=cos(theta);

  // Kinematic quantities for exchange particle
  double q0=p_gamma-sqrt(p_eta*p_eta+m_eta_sq);
  double q3=p_gamma-p_eta*costheta;
  double q0_minus_q3=q0-q3;
  double q0_minus_q3_sq=q0_minus_q3*q0_minus_q3;
  double q0_sq=q0*q0;
  double q1sq_plus_q2sq=p_eta*p_eta*sintheta*sintheta;
  // Kinematic quantities for target
  double pt3=-p_gamma;
  double pt0=sqrt(m_p_sq+pt3*pt3);
  double pt0_minus_pt3=pt0-pt3;
  double pt0_plus_pt3=pt0+pt3;
  // Kinematic quantities for beam
  double p_gamma_sq=p_gamma*p_gamma;
  // other kinematic quantities
  double pt_dot_q=pt0*q0-pt3*q3;

 // Mass scale for Regge propagators
  double s0=1.0;

  // Regge trajectory for rho
  double a_rho=0.55+0.8*t;
  double a_rho_prime=0.8;
  double regge_rho=pow(s/s0,a_rho-1.)*M_PI*a_rho_prime/(sin(M_PI*a_rho)*TMath::Gamma(a_rho));

  // Regge trajectory for omega
  double a_omega=0.44+0.9*t;
  double a_omega_prime=0.9;
  double regge_omega=pow(s/s0,a_omega-1.)*M_PI*a_omega_prime/(sin(M_PI*a_omega)*TMath::Gamma(a_omega));

  // Regge trajectory for phi(1020)
  double a_phi=0.23+0.7*t;
  double a_phi_prime=0.7;
  double regge_phi=pow(s/s0,a_phi-1.)*M_PI*a_phi_prime/(sin(M_PI*a_phi)*TMath::Gamma(a_phi));

  //  printf("> 36:  %f\n",p_gamma*p_eta*sqrt(mass_factor));

  // amplitude factors for terms involving 0, 1 and 2 powers of kappa
  double amp_factor_kappa0
    =8.*p_gamma_sq*(q1sq_plus_q2sq*(s+q0_minus_q3*pt0_minus_pt3)
		    +q0_minus_q3_sq*pt_dot_q);
  double amp_factor_kappa1=32.*p_gamma_sq*m_p*(q0_minus_q3_sq*t
					       +2.*q0_sq*q1sq_plus_q2sq);
  double amp_factor_kappa2
    =32.*p_gamma_sq*(q1sq_plus_q2sq*(2.*q0_sq*(pt_dot_q-2.*m_p_sq)
				     -t*pt0_plus_pt3*pt0_plus_pt3
				     +4.*pt_dot_q*q0*pt0_plus_pt3)
		     +q0_minus_q3_sq*(t*(pt_dot_q-2.*m_p_sq)
				      +2.*pt_dot_q*pt_dot_q)
		     );

  // rho amplitude 
  double M_rho_sq=16.*M_PI*M_PI*(g_rho_eta_gamma*g_rho_eta_gamma/m_eta_sq)*c_rho_p_p*regge_rho*regge_rho
    *(amp_factor_kappa0-kappa_rho/(4.*m_p)*amp_factor_kappa1
      +kappa_rho*kappa_rho/(16.*m_p_sq)*amp_factor_kappa2);
  double M_rho=-sqrt(M_rho_sq);

  // omega amplitude 
  double M_omega_sq=16.*M_PI*M_PI*(g_omega_eta_gamma*g_omega_eta_gamma/m_eta_sq)*c_omega_p_p*regge_omega*regge_omega
    *amp_factor_kappa0;
  double M_omega=-sqrt(M_omega_sq);

  // phi amplitude 
  double M_phi_sq=16.*M_PI*M_PI*(g_phi_eta_gamma*g_phi_eta_gamma/m_eta_sq)*c_phi_p_p*regge_phi*regge_phi
    *amp_factor_kappa0;
  double M_phi=+sqrt(M_phi_sq);

  // Primakoff amplitude 
  double M_primakoff_sq=16.*M_PI*M_PI*(g_eta_gamma_gamma*g_eta_gamma_gamma/m_eta_sq)*c_gamma_p_p/(t*t)
    *(amp_factor_kappa0-kappa_gamma/(4.*m_p)*amp_factor_kappa1+kappa_gamma*kappa_gamma/(16.*m_p_sq)*amp_factor_kappa2);
  double M_primakoff=sqrt(M_primakoff_sq);

    
  // M_primakoff=0.;
  //M_primakoff_sq=0.;
  
  //M_omega=0.;
  // M_omega_sq=0.;
  
  //M_rho=0.;
  // M_rho_sq=0.;
 
  //M_phi_sq=0.;
  //M_phi=0.;
  
  double pi_a_omega=M_PI*a_omega;
  double pi_a_rho=M_PI*a_rho;
  double pi_a_phi=M_PI*a_phi;
  double M_sq =M_omega_sq+M_rho_sq+M_primakoff_sq+M_phi_sq
    +2.*M_omega*M_phi*cos(pi_a_omega-pi_a_phi)
    +2.*M_omega*M_rho*cos(pi_a_omega-pi_a_rho)
    +2.*M_omega*M_primakoff*cos(pi_a_omega)
    +2.*M_rho*M_primakoff*cos(pi_a_rho)
    +2.*M_rho*M_phi*cos(pi_a_rho-pi_a_phi)
    +2.*M_phi*M_primakoff*cos(pi_a_phi)
    ;
  
  double hbarc_sq=389.; // Convert to micro-barns
  double dsigma_dt=hbarc_sq*M_sq/(4.*64.*M_PI*s*p_gamma_sq);
  // the extra factor for is for 2 photon spins x 2 proton spins

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

  thrown_t=new TH1D("thrown_t","Thrown -t distribution",1000,0.,2.0);
  thrown_t->SetXTitle("-t [GeV^{2}]");
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
void GraphCrossSection(){
  // beam energy in lab
  double Egamma=Emin;

  // CM energy
  double s=m_p*(m_p+2.*Egamma);
  double Ecm=sqrt(s);

  // Momenta of incoming photon and outgoing eta and proton in cm frame
  double p_gamma=(s-m_p_sq)/(2.*Ecm);
  double E_eta=(s+m_eta_sq-m_p_sq)/(2.*Ecm);
  double p_eta=sqrt(E_eta*E_eta-m_eta_sq);
  
  // Momentum transfer t
  double p_diff=p_gamma-p_eta;
  double t0=m_eta_sq*m_eta_sq/(4.*s)-p_diff*p_diff;
  
  double sum=0.;
  double t_old=t0;
  double t_array[10000];
  double xsec_array[10000];
  for (unsigned int k=0;k<10000;k++){
    double theta_cm=M_PI*double(k)/10000.;
    double sin_theta_over_2=sin(0.5*theta_cm);
    double t=t0-4.*p_gamma*p_eta*sin_theta_over_2*sin_theta_over_2;
    double xsec=CrossSection(s,t,p_gamma,p_eta,theta_cm);
    
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
  string rootfilename="eta_gen.root";
  TFile *rootfile=new TFile(rootfilename.c_str(),"RECREATE",
			    "Produced by genEta");

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
  if (Emin<m_eta){ 
    Emin=m_eta;
    cout << "Warning:  Setting minimum beam energy to " << Emin << " [GeV]" 
	 <<endl;
  }
  cout << "Photon energy min, max [Gev] = "<< Emin <<","<<Emax <<endl;

  // Get coherent peak and collimator diameter
  getline(infile,comment_line);
  float Epeak=9.0,collDiam=0.005,radThickness=50e-6;
  float Ee=12.0;
  infile >> Ee;
  infile >> Epeak;
  infile >> collDiam;
  infile >> radThickness;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "Electron beam energy = " << Ee << " GeV, Coherent peak = " 
       << Epeak <<" GeV, collimator diameter = " 
       <<collDiam << " m, radiator thickness = " 
       << radThickness << " m"
       << endl;

  // Get decaying particle mass and width
  string comment_line2;
  getline(infile,comment_line);
  infile >> m_eta;
  infile >> width;
  infile.ignore(); // ignore the '\n' at the end of this line

  m_eta_sq=m_eta*m_eta;
  cout << "Mass, width of decaying particle [GeV] = "<< m_eta <<"," << width << endl;

  // Get coupling constants for photon vertex
  getline(infile,comment_line);
  infile >> g_eta_gamma_gamma;
  infile >> g_rho_eta_gamma;
  infile >> g_omega_eta_gamma;
  infile >> g_phi_eta_gamma;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "Coupling constants:" <<endl;
  cout << " g_eta_gamma_gamma = " << g_eta_gamma_gamma <<endl;
  cout << " g_rho_eta_gamma = " << g_rho_eta_gamma <<endl;
  cout << " g_omega_eta_gamma = " << g_omega_eta_gamma << endl;
  cout << " g_phi_eta_gamma = " << g_phi_eta_gamma <<endl;
  
  // Get number of decay particles
  int num_decay_particles=0;
  getline(infile,comment_line);
  infile >> num_decay_particles;
  infile.ignore(); // ignore the '\n' at the end of this line

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
  cout << endl;

  infile.close();
 
  //----------------------------------------------------------------------------
  // Setup coherent bremsstrahlung generator
  //----------------------------------------------------------------------------
  float radColDist=76.0;// meters
  int doPolFlux=0;  // want total flux (1 for polarized flux)
  float emitmr=10.e-9; // electron beam emittance
  CobremsGeneration cobrems(Ee, Epeak);
  cobrems.setBeamEmittance(emitmr);
  cobrems.setTargetThickness(radThickness);
  cobrems.setCollimatorDistance(radColDist);
  cobrems.setCollimatorDiameter(collDiam);
  cobrems.setCollimatedFlag(true);
  cobrems.setPolarizedFlag(doPolFlux);
  cobrems.setCollimatedFlag(true);
  
  // Create some diagonistic histographs
  CreateHistograms();
  // Make a TGraph of the cross section at a fixed beam energy
  GraphCrossSection();
  
  // Fill histogram of coherent bremmsstrahlung distribution 
  for (int i=1;i<=1000;i++){
    float x=float(cobrems_vs_E->GetBinCenter(i)/Ee);
    float y=0;
    if (Epeak<Emin) y=cobrems.Rate_dNidx(x);
    else y=cobrems.Rate_dNtdx(x);
    cobrems_vs_E->Fill(Ee*double(x),double(y));
  }


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
    double p_eta=0.;

    // Transfer 4-momentum;
    double t=0.;

    // vertex position at target
    float vert[4]={0.,0.,0.,0.};

    // use the rejection method to produce eta's based on the cross section
    do{
      // First generate a beam photon using bremsstrahlung spectrum
      Egamma = cobrems_vs_E->GetRandom();

      // CM energy
      double s=m_p*(m_p+2.*Egamma);
      double Ecm=sqrt(s);

      // Momenta of incoming photon and outgoing eta and proton in cm frame
      double p_gamma=(s-m_p_sq)/(2.*Ecm);

      if (width>0){  // Take into account width of resonance
	// Use a relativistic Breit-Wigner distribution for the shape

      }
      double E_eta=(s+m_eta_sq-m_p_sq)/(2.*Ecm);
      p_eta=sqrt(E_eta*E_eta-m_eta_sq);
    
      // Momentum transfer t
      double p_diff=p_gamma-p_eta;
      double t0=m_eta_sq*m_eta_sq/(4.*s)-p_diff*p_diff;
      double sin_theta_over_2=0.;
      t=t0;
      
      // Generate cos(theta) with a uniform distribution and compute the cross 
      // section at this value
      double cos_theta_cm=-1.0+myrand->Uniform(2.);
      theta_cm=acos(cos_theta_cm);
      
      sin_theta_over_2=sin(0.5*theta_cm);
      t=t0-4.*p_gamma*p_eta*sin_theta_over_2*sin_theta_over_2;
      xsec=CrossSection(s,t,p_gamma,p_eta,theta_cm);
      
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
    // Four-moementum of the eta in the CM frame
    double pt=p_eta*sin(theta_cm);
    TLorentzVector eta4(pt*cos(phi_cm),pt*sin(phi_cm),p_eta*cos(theta_cm),
			sqrt(p_eta*p_eta+m_eta_sq));
    // eta4.Print();

    //Boost the eta 4-momentum into the lab
    eta4.Boost(v_cm);
    // eta4.Print();

  
    // Compute the 4-momentum for the recoil proton
    TLorentzVector proton4=beam+target-eta4;

    //proton4.Print();
    thrown_theta_vs_p->Fill(proton4.P(),180./M_PI*proton4.Theta());

    // Generate 3-body decay of eta according to phase space
    TGenPhaseSpace phase_space;
    phase_space.SetDecay(eta4,num_decay_particles,decay_masses);
    double weight=0.,rand_weight=1.;
    do{
      weight=phase_space.Generate();
      rand_weight=myrand->Uniform(1.);
    }
    while (rand_weight>weight);

    // Histograms of Dalitz distribution
    if (num_decay_particles==3){
      TLorentzVector one=*phase_space.GetDecay(0);  
      TLorentzVector two=*phase_space.GetDecay(1);
      TLorentzVector three=*phase_space.GetDecay(2);
      TLorentzVector one_two=one+two;
      TLorentzVector one_three=one+three;
 
      TLorentzVector eta=one_two+three;
      TVector3 boost=-eta.BoostVector();

      double eta_mass=eta.M();
      eta.Boost(boost);
      
      one.Boost(boost);
      two.Boost(boost);
      three.Boost(boost);

      double m1=one.M(),m2=two.M(),m3=three.M();
      double E1=one.E(),E2=two.E(),E3=three.E();
      double T1=E1-m1; // pi0 for charged channel
      double T2=E2-m2; // pi+ for charged channel
      double T3=E3-m3; // pi- for charged channel
      double Q_eta=eta_mass-m1-m2-m3;
      double X=sqrt(3.)*(T2-T3)/Q_eta;
      double Y=3.*T1/Q_eta-1.;
      thrown_dalitzXY->Fill(X,Y);

      double z_dalitz=X*X+Y*Y;
      //printf("z %f\n",z_dalitz);
      thrown_dalitzZ->Fill(z_dalitz);
    }
    // Other diagnostic histograms
    thrown_t->Fill(-t);

    // Randomly generate z position in target
    vert[2]=zmin+myrand->Uniform(zmax-zmin);

    // Gather the particles in the reaction and write out event in hddm format
    particle_vectors[last_index]=proton4;
    for (int j=0;j<num_decay_particles;j++){
      particle_vectors[j]=*phase_space.GetDecay(j);
    }
    WriteEvent(i,beam,vert,particle_types,particle_vectors,file);
    
    if ((i%(Nevents/10))==0) cout << 100.*double(i)/double(Nevents) << "\% done" << endl;
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

#else

// ---- The following is compiled only if AMPTOOLS was not installed ----
#include <iostream>

int main(int narg, char *argv[])
{
	std::cerr << "This executable (genEtaRegge) was built without AMPTOOLS support." << std::endl;
	std::cerr << "As such it is disabled. To use it, install AMPTOOLS and make sure" << std::endl;
	std::cerr << "your AMPTOOLS environment variable is set to point to it and then" << std::endl;
	std::cerr << "rebuild." << std::endl;

	return 0;
}

#endif   // HAVE_AMPTOOLS_MCGEN

