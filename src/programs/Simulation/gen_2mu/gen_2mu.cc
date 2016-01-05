
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

#include <HDDM/hddm_s.hpp>

#include "particleType.h"

#include <TLorentzVector.h>
#include <TRandom2.h>
#include <TF1.h>

#include "GlueXPrimaryGeneratorAction.hh"

string OUTPUT_FILENAME = "gen_2mu.hddm";
int32_t RUN_NUMBER = 2;
uint32_t MAXEVENTS = 10000;
double Z = 82.0;  // atomic number of target
double A = 208.0; // atomic weight of target
bool HDDM_USE_COMPRESSION = false;
bool HDDM_USE_INTEGRITY_CHECKS = false;
bool USE_ELECTRON_BEAM_DIRECTION = false;
double POLARIZATION_ANGLE = 0.0; // in degrees relative to x-axis

TRandom *RAND = NULL;

double Ecoherent_peak = 6.0;
double Eelectron_beam = 12.0;
double Emin = 1.0;
bool   ONLY_COHERENT = false;
bool   ONLY_INCOHERENT = false;


static ofstream *OFS = NULL;
static hddm_s::ostream *FOUT = NULL;

extern void GetMech(int &Ncoherent, int &Nincoherent);
extern int LAST_COBREMS_MECH;

void GenerateMuPair(TVector3 &pgamma, TVector3 &pol, TLorentzVector &pmuplus, TLorentzVector &pmuminus);
void AddEventToHDDM(TVector3 &pgamma, TLorentzVector &pmuplus, TLorentzVector &pmuminus);
void Usage(string message="");
void ParseCommandLineArguments(int narg, char *argv[]);

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


//-----------------------
// main
//-----------------------
int main( int narg, char* argv[] )
{
	// Parse the command line arguments
	ParseCommandLineArguments(narg, argv);

	// Random number generator
	RAND = new TRandom2(1);

	// Coherent bremstrahlung photon generator
	GlueXPrimaryGeneratorAction *photon_generator = new GlueXPrimaryGeneratorAction();
	
	// Open output HDDM file
	if( !FOUT ){

		OFS = new std::ofstream(OUTPUT_FILENAME.c_str());
		if (OFS->is_open()){
			try{
				FOUT = new hddm_s::ostream(*OFS);
			}catch(exception &e){
				cout << e.what() << endl;
			}
		}
		if( !FOUT ){
		  cout << " Error opening output file \"" << OUTPUT_FILENAME << "\"!" << endl;
		  exit(-1);
		}
	}

	cout << "Opened output file: " << OUTPUT_FILENAME << " 0x" << FOUT << ")" << endl;

	// enable on-the-fly bzip2 compression on output stream
	if (HDDM_USE_COMPRESSION == 0) {
		cout << " HDDM compression disabled" << endl;
	} else if (HDDM_USE_COMPRESSION == 1) {
		cout << " Enabling bz2 compression of output HDDM file stream" << endl;
		FOUT->setCompression(hddm_s::k_bz2_compression);
	} else {
		cout << " Enabling z compression of output HDDM file stream (default)" << endl;
		FOUT->setCompression(hddm_s::k_z_compression);
	}

	// enable a CRC data integrity check at the end of each event record
	if (HDDM_USE_INTEGRITY_CHECKS) {
		cout << " Enabling CRC data integrity check in output HDDM file stream (default)" << endl;
		FOUT->setIntegrityChecks(hddm_s::k_crc32_integrity);
	} else {
		cout << " HDDM integrity checks disabled" << endl;
	}

	// Start event processing
	cout << "Event generation starting ..." << endl;
	uint32_t Nevents_generated = 0;
	for(Nevents_generated=0; Nevents_generated<MAXEVENTS; Nevents_generated++){
	
		// Generate beam photon
		TVector3 pgamma, pol;
		photon_generator->GenerateBeamPhoton(pgamma, pol);
		
		// Generate mu+mu- pair
		TLorentzVector pmuplus, pmuminus;
		GenerateMuPair(pgamma, pol, pmuplus, pmuminus);
		
		// Write event to file
		AddEventToHDDM(pgamma, pmuplus, pmuminus);
		
		// Update ticker so user knows we're working
		if(Nevents_generated%100 == 0){
			cout << "  " << Nevents_generated << " events generated      \r";
			cout. flush();
		}
	
	}

	// Delete objects (closing output file)
	if(photon_generator) delete photon_generator;
	if(RAND) delete RAND;
	if( FOUT && Nevents_generated>0) { // (program crashes if we close without writing any events!)
		cout << "Closing file: " << OUTPUT_FILENAME << " (wrote " << Nevents_generated << " events)" << endl;
		if(FOUT) delete FOUT;
		if(OFS)  delete OFS;
		FOUT = NULL;
		OFS = NULL;
	}
	
	int Ncoherent=0, Nincoherent=0;
	GetMech(Ncoherent, Nincoherent);
	cout << "  Ncoherent = " << Ncoherent << endl;
	cout << "Nincoherent = " << Nincoherent << endl;
	cout << Nevents_generated << " events generated Total." << endl;
  
	return 0;
}

//-----------------------
// Usage
//-----------------------
void Usage(string message)
{
	cout << endl;
	cout << "Usage:" << endl;
	cout << "    gen_2mu [options]" << endl;
	cout << endl;
	cout << " -h           print this help message" << endl;
	cout << " --help       print the long form help message" << endl;
	cout << " -N events    number of events to generate" << endl;
	cout << " -p Epeak     coherent peak energy (def="<<Ecoherent_peak<<")" << endl;
	cout << " -b Ebeam     electron beam energy (def="<<Eelectron_beam<<")" << endl;
	cout << " -min Emin    minimum photon energy to generate (def="<<Emin<<")" << endl;
	cout << " -c           only generate coherent photons" << endl;
	cout << " -i           only generate incoherent photons" << endl;
	cout << " -e           let electron direction define z (def." << endl;
	cout << "              is for photon beam to define z)" << endl;
	cout << " -pol phi     set photon beam polarization direction" << endl;
	cout << "              relative to x-axis (def. is " << POLARIZATION_ANGLE << " degrees)" << endl; 
	cout << endl;
	if(message!=""){
		cout << message << endl;
		cout << endl;
		exit(-1);
	}

	exit(0);
}

//-----------------------
// ParseCommandLineArguments
//-----------------------
void ParseCommandLineArguments(int narg, char *argv[])
{
	for(int i=1; i<narg; i++){
		string arg = argv[i];
		string next = (i+1)<narg ? argv[i+1]:"";
		bool has_arg = (next!="") && (next.find("-")!=0); // true if first character is not "-"
		bool missing_arg = false; // flags if option requires second argument and it's missing
		
		if(arg=="-h"){
			Usage();
		}else if(arg=="--help"){
			cout << "Sorry, no long-form help yet! (try \"-h\")" << endl;
			exit(0);
		}else if(arg=="-N"){
			if(has_arg){
				MAXEVENTS = atoi(next.c_str());
				i++;
			}else{
				missing_arg = true;
			}
		}else if(arg=="-p"){
			if(has_arg){
				Ecoherent_peak = atof(next.c_str());
				i++;
			}else{
				missing_arg = true;
			}
		}else if(arg=="-b"){
			if(has_arg){
				Eelectron_beam = atof(next.c_str());
				i++;
			}else{
				missing_arg = true;
			}
		}else if(arg=="-min"){
			if(has_arg){
				Emin = atof(next.c_str());
				i++;
			}else{
				missing_arg = true;
			}
		}else if(arg=="-c"){
			ONLY_COHERENT = true;
		}else if(arg=="-i"){
			ONLY_INCOHERENT = true;
		}else if(arg=="-i"){
			USE_ELECTRON_BEAM_DIRECTION = true;
		}else{
			stringstream ss;
			ss << "Unknown argument \""<<arg<<"\"!" << endl;
			Usage(ss.str());
		}

		// Check if option required argument and is missing it
		if(missing_arg){
			stringstream ss;
			ss << "Option \""<<arg<<"\" requires an argument!";
			Usage(ss.str());
		}
	}
	
	if(ONLY_COHERENT && ONLY_INCOHERENT){
		Usage("Cannot specify both -c and -i !");
	}
	
	cout << endl;
	cout << "---------------------------------------------" << endl;
	cout << "     Nevents to generate: " << MAXEVENTS << endl;
	cout << "    Electron beam energy: " << Eelectron_beam << " GeV" << endl;
	cout << "           Coherent peak: " << Ecoherent_peak << " GeV" << endl;
	cout << "   Minimum photon energy: " << Emin << " GeV" << endl;
	cout << "  z-direction defined by: " << (USE_ELECTRON_BEAM_DIRECTION ? "electron":"photon") << " beam" << endl;
	cout <<"   Coherent photons only?: " << (ONLY_COHERENT ? "yes":"no") << endl;
	cout <<" Incoherent photons only?: " << (ONLY_INCOHERENT ? "yes":"no") << endl;
	cout << "---------------------------------------------" << endl;
	cout << endl;
}

//-----------------------
// GenerateMuPair
//-----------------------
void GenerateMuPair(TVector3 &pgamma, TVector3 &pol, TLorentzVector &pmuplus, TLorentzVector &pmuminus)
{
	/// This method was copied (and modified slightly) from the
	/// GEANT 4.10.02 G4GammaConversionToMuons::PostStepDoIt method.x
	/// Given an input photon momentum vector and polarization
	/// vector (pgamma, and pol), it will generate a mu+mu-
	/// pair, returning the 4-vectors in the pmuplus, pmuminus
	/// TLorentzVector objects.

//	aParticleChange.Initialize(aTrack);
//	G4Material* aMaterial = aTrack.GetMaterial();
//
//	// current Gamma energy and direction, return if energy too low
//	const G4DynamicParticle *aDynamicGamma = aTrack.GetDynamicParticle();
//	double Egam = aDynamicGamma->GetKineticEnergy();
//	if (Egam <= LowestEnergyLimit) {
//	return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
//	}
//	G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();
	
	double Mmuon = 0.1056583715;
	double electron_mass_c2 = 0.000511;
	double sqrte = 1.648721270700128;
	double pi = 3.141592653589793;

	TVector3 GammaDirection(pgamma);
	GammaDirection.SetMag(1.0);
	double Egam = pgamma.Mag();

//	// select randomly one element constituting the material
//	const G4Element* anElement = SelectRandomAtom(aDynamicGamma, aMaterial);
//	int Z = G4lrint(anElement->GetZ());
//	G4NistManager* nist = G4NistManager::Instance();

	double B,Dn;
//	double A027 = nist->GetA27(Z);
	double A027 = pow(A, 0.27);

	if(Z==1) // special case of Hydrogen
	{ B=202.4;
	  Dn=1.49;
	}
	else
	{ B=183.;
	  Dn=1.54*A027;
	}
//	double Zthird=1./nist->GetZ13(Z); // Z**(-1/3)
	double Zthird=pow(Z, -1.0/3.0); // Z**(-1/3)
	double Winfty=B*Zthird*Mmuon/(Dn*electron_mass_c2);
	double C1Num=0.35*A027;
	double C1Num2=C1Num*C1Num;
	double C2Term2=electron_mass_c2/(183.*Zthird*Mmuon);

	double GammaMuonInv=Mmuon/Egam;
	double sqrtx=sqrt(.25-GammaMuonInv);
	double xmax=.5+sqrtx;
	double xmin=.5-sqrtx;

	// generate xPlus according to the differential cross section by rejection
	double Ds2=(Dn*sqrte-2.);
	double sBZ=sqrte*B*Zthird/electron_mass_c2;
	double LogWmaxInv=1./log(Winfty*(1.+2.*Ds2*GammaMuonInv)
				   /(1.+2.*sBZ*Mmuon*GammaMuonInv));
	double xPlus,xMinus,xPM,result,W;
	int nn = 0;
	const int nmax = 1000;
	do
	{ xPlus=xmin+RAND->Rndm()*(xmax-xmin);
	xMinus=1.-xPlus;
	xPM=xPlus*xMinus;
	double del=Mmuon*Mmuon/(2.*Egam*xPM);
	W=Winfty*(1.+Ds2*del/Mmuon)/(1.+sBZ*del);
	if(W<=1. || nn > nmax) { break; } // to avoid negative cross section at xmin
	double xxp=1.-4./3.*xPM; // the main xPlus dependence
	result=xxp*log(W)*LogWmaxInv;
	if(result>1.) {
	  cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
		 << " in dSigxPlusGen, result=" << result << " > 1" << endl;
	}
	++nn;
	if(nn >= nmax) { break; }
	}
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	while (RAND->Rndm() > result);

	// now generate the angular variables via the auxilary variables t,psi,rho
	double t;
	double psi;
	double rho;

	double thetaPlus,thetaMinus,phiHalf; // final angular variables
	nn = 0;
	do      // t, psi, rho generation start  (while angle < pi)
	{
	//generate t by the rejection method
	double C1=C1Num2* GammaMuonInv/xPM;
	double f1_max=(1.-xPM) / (1.+C1);
	double f1; // the probability density
	do
	{ 
	  ++nn;
	  t=RAND->Rndm();
	  f1=(1.-2.*xPM+4.*xPM*t*(1.-t)) / (1.+C1/(t*t));
	  if(f1<0 || f1> f1_max) // should never happend
	{
	  cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
		 << "outside allowed range f1=" << f1 << " is set to zero"
		 << endl;
		  f1 = 0.0;
	}
	  if(nn > nmax) { break; }
	}
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	while ( RAND->Rndm()*f1_max > f1);
	// generate psi by the rejection method
	double f2_max=1.-2.*xPM*(1.-4.*t*(1.-t));

	// long version
	double f2;
	do
	{ 
	  ++nn;
	  psi=2.*pi*RAND->Rndm();
	  f2=1.-2.*xPM+4.*xPM*t*(1.-t)*(1.+cos(2.*psi));
	  if(f2<0 || f2> f2_max) // should never happend
	{
	  cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
		 << "outside allowed range f2=" << f2 << " is set to zero"
		 << endl;
		  f2 = 0.0;
	}
	  if(nn >= nmax) { break; }
	}
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	while ( RAND->Rndm()*f2_max > f2);

	// generate rho by direct transformation
	double C2Term1=GammaMuonInv/(2.*xPM*t);
	double C2=4./sqrt(xPM)*pow(C2Term1*C2Term1+C2Term2*C2Term2,2.);
	double rhomax=1.9/A027*(1./t-1.);
	double beta=log( (C2+rhomax*rhomax*rhomax*rhomax)/C2 );
	rho=exp(log(C2 *( exp(beta*RAND->Rndm())-1. ))*0.25);

	//now get from t and psi the kinematical variables
	double u=sqrt(1./t-1.);
	double xiHalf=0.5*rho*cos(psi);
	phiHalf=0.5*rho/u*sin(psi);

	thetaPlus =GammaMuonInv*(u+xiHalf)/xPlus;
	thetaMinus=GammaMuonInv*(u-xiHalf)/xMinus;

	// protection against infinite loop
	if(nn > nmax) {
	  if(std::abs(thetaPlus)>pi) { thetaPlus = 0.0; }
	  if(std::abs(thetaMinus)>pi) { thetaMinus = 0.0; }
	}

	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	} while ( std::abs(thetaPlus)>pi || std::abs(thetaMinus) >pi);

	// now construct the vectors
	// azimuthal symmetry, take phi0 at random between 0 and 2 pi
	double phi0=2.*pi*RAND->Rndm(); 
	double EPlus=xPlus*Egam;
	double EMinus=xMinus*Egam;
	

	// ----------------------------------------------------

	// mu+ mu- directions for gamma in z-direction
	TVector3 MuPlusDirection  ( sin(thetaPlus) *cos(phi0+phiHalf),
				   sin(thetaPlus)  *sin(phi0+phiHalf), cos(thetaPlus) );
	TVector3 MuMinusDirection (-sin(thetaMinus)*cos(phi0-phiHalf),
				  -sin(thetaMinus) *sin(phi0-phiHalf), cos(thetaMinus) );
	
	double Pplus = sqrt(EPlus*EPlus - Mmuon*Mmuon);
	double Pminus = sqrt(EMinus*EMinus - Mmuon*Mmuon);
	
	pmuplus.SetVectM(Pplus*MuPlusDirection, Mmuon);
	pmuminus.SetVectM(Pminus*MuMinusDirection, Mmuon);
	
	//-- Add phi dependence on photon polarization --
	//
	// The polarization causes an azimuthal dependence of the mu+mu-
	// system about the incident photon direction relative to the
	// polarization vector that goes like:
	//
	//  ( 1 - 2cos(2phi) )
	//
	// see https://halldweb1.jlab.org/wiki/images/a/aa/20130418_cpp_rory.pdf
	//
	// To do this, we take the current phi angle of the mu+mu- system and assume
	// it is evenly distributed over 0-2pi. We normalize this to get number from
	// 0-1 and equate it with the normalized integral of the above function of phi.
	// (The integral fraction method). This results in a transcendental equation
	// though so we use the ROOT TF1::GetX() method to find the root of a function
	// defined as the difference between the "random number" and the normalized
	// integral.
	//
	//    f(phi) = s - (phi - 0.5*sin(2*phi))/2pi
	//
	//  where:
	//      s = normalized phi from unpolarized (0 - 1)
	//    phi = azimuthal angle of mu+mu- relative to polarization direction
	//     pi = 3.14159......
	//

	if(LAST_COBREMS_MECH == 1){ // Only do this for coherently produced photons
		static TF1 *normInt = NULL;
		if(!normInt){
			// par0 is random number "s"
			// 0.1591549430919 = 1/2pi
			normInt = new TF1("normInt", "[0] - (x - 0.5*sin(2.0*x))*0.1591549430919", 0.0, TMath::TwoPi());
		}

		// direction of mu+mu- system
		TVector3 vmumu = (pmuplus+pmuminus).Vect();
		double phi_init = vmumu.Phi();
		double s = 0.5+phi_init/TMath::TwoPi(); // s is 0-1
		normInt->SetParameter(0, s);
		double deltaphi = normInt->GetX(0.0, 0.0, TMath::TwoPi()) - phi_init;
		deltaphi += POLARIZATION_ANGLE*TMath::DegToRad();
		pmuplus.RotateZ(deltaphi);
		pmuminus.RotateZ(deltaphi);
	}

	// Rotate to actual gamma direction.
	// The pmuplus and pmuminus vectors are currently relative to the
	// beam photon direction. This is almost always what you want.
	// This gives the option though of rotating to the direction where
	// z is defined by the electron beam. For coherently produced photons,
	// this option introduces an effective phi shift in the phi_mumu
	// distribution since the gamma direction is concentrated in one
	// region of phi (roughly 35 degrees)
	if(USE_ELECTRON_BEAM_DIRECTION){
		pmuplus.RotateUz(GammaDirection);
		pmuminus.RotateUz(GammaDirection);
	}
	
//	aParticleChange.SetNumberOfSecondaries(2);
//	// create G4DynamicParticle object for the particle1
//	G4DynamicParticle* aParticle1= new G4DynamicParticle(
//						   G4MuonPlus::MuonPlus(),MuPlusDirection,EPlus-Mmuon);
//	aParticleChange.AddSecondary(aParticle1);
//	// create G4DynamicParticle object for the particle2
//	G4DynamicParticle* aParticle2= new G4DynamicParticle(
//					   G4MuonMinus::MuonMinus(),MuMinusDirection,EMinus-Mmuon);
//	aParticleChange.AddSecondary(aParticle2);
//	//
//	// Kill the incident photon
//	//
//	aParticleChange.ProposeMomentumDirection( 0., 0., 0. ) ;
//	aParticleChange.ProposeEnergy( 0. ) ;
//	aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
//	//  Reset NbOfInteractionLengthLeft and return aParticleChange
//	return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

//-----------------------
// AddEventToHDDM
//-----------------------
void AddEventToHDDM(TVector3 &pgamma, TLorentzVector &pmuplus, TLorentzVector &pmuminus)
{
	using namespace hddm_s;
	static uint32_t event_number = 0;
	int mech = LAST_COBREMS_MECH; // 0=unknown, 1=coherent, 2=incoherent

	HDDM *hddmevent = new HDDM;
	hddmevent->addPhysicsEvents(1);
	PhysicsEvent &PE = hddmevent->getPhysicsEvent();
	PE.setRunNo( RUN_NUMBER );
	PE.setEventNo( ++event_number );

	ReactionList reactions = PE.addReactions();
	VertexList vertices = reactions().addVertices();
	
	// Add Beam
	BeamList beam = reactions().addBeams();
	beam().setType(Gamma);
	MomentumList momenta = beam().addMomenta();
	momenta().setE( pgamma.Mag() );
	momenta().setPx( pgamma.x() );
	momenta().setPy( pgamma.y() );
	momenta().setPz( pgamma.z() );
	PropertiesList properties = beam().addPropertiesList();
	properties().setCharge( 0 );
	properties().setMass( 0 );
		
	// Add Origin
	TVector3 pos(0.0, 0.0, 1.0);
	OriginList origins = vertices().addOrigins();
	origins().setT(0.0);
	origins().setVx(pos.x());
	origins().setVy(pos.y());
	origins().setVz(pos.z());
	
	// Add Products (particles)
	ProductList products = vertices().addProducts(2);
	ProductList::iterator it_product = products.begin();

	// Product Mu+
	Particle_t geanttype = MuonPlus;
	TVector3 mom = pmuplus.Vect(); // convert back to units of GeV
	double mass = ParticleMass(geanttype);
	it_product->setDecayVertex(0);
	it_product->setId(1);
	it_product->setMech(mech);
	it_product->setParentid(0);
	it_product->setType(geanttype);
	it_product->setPdgtype(PDGtype(geanttype));

	// Momentum Mu+
	momenta = it_product->addMomenta();
	momenta().setE( sqrt(mom.Mag2() + mass*mass) );
	momenta().setPx( mom.x() );
	momenta().setPy( mom.y() );
	momenta().setPz( mom.z() );

	// Properties Mu+
	properties = it_product->addPropertiesList();
	properties().setCharge( ParticleCharge(geanttype) );
	properties().setMass( mass );


	it_product++;
	

	// Product Mu-
	geanttype = MuonMinus;
	mom = pmuminus.Vect(); // convert back to units of GeV
	mass = ParticleMass(geanttype);
	it_product->setDecayVertex(0);
	it_product->setId(2);
	it_product->setMech(mech);
	it_product->setParentid(0);
	it_product->setType(geanttype);
	it_product->setPdgtype(PDGtype(geanttype));

	// Momentum Mu-
	momenta = it_product->addMomenta();
	momenta().setE( sqrt(mom.Mag2() + mass*mass) );
	momenta().setPx( mom.x() );
	momenta().setPy( mom.y() );
	momenta().setPz( mom.z() );

	// Properties Mu-
	properties = it_product->addPropertiesList();
	properties().setCharge( ParticleCharge(geanttype) );
	properties().setMass( mass );
	
	(*FOUT) << (*hddmevent);

	delete hddmevent;
}

