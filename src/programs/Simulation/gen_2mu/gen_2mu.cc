
#include <stdint.h>

#include <iostream>
//#include <fstream>
#include <string>
#include <vector>
using namespace std;


//#include "particleType.h"

#include <TLorentzVector.h>
#include <TRandom2.h>

#include "GlueXPrimaryGeneratorAction.hh"

TRandom *RAND = NULL;
uint32_t MAXEVENTS = 10;
double Z = 82.0;  // atomic number of target
double A = 208.0; // atomic weight of target


void GenerateMuPair(TVector3 &pgamma, TVector3 &pol, TLorentzVector &pmuplus, TLorentzVector &pmuminus);


//-----------------------
// main
//-----------------------
int main( int argc, char* argv[] )
{
	RAND = new TRandom2(1);

	GlueXPrimaryGeneratorAction *photon_generator = new GlueXPrimaryGeneratorAction();


	for(uint32_t i=0; i<MAXEVENTS; i++){
	
		// Generate beam photon
		TVector3 pgamma, pol;
		photon_generator->GenerateBeamPhoton(pgamma, pol);
		
		// Generate mu+mu- pair
		TLorentzVector pmuplus, pmuminus;
		GenerateMuPair(pgamma, pol, pmuplus, pmuminus);
		
		// Write event to file
	
	}

	if(photon_generator) delete photon_generator;
	if(RAND) delete RAND;
  
//	string  configfile("");
//	string  outname("");
//	string  hddmname("");
//	
//	bool diag = false;
//	bool genFlat = false;
//	
//	// default upper and lower bounds 
//	double lowMass = 0.2;
//	double highMass = 2.0;
//	
//	double beamMaxE   = 12.0;
//	double beamPeakE  = 9.0;
//	double beamLowE   = 0.139*2;
//	double beamHighE  = 12.0;
//	
//	int runNum = 9001;
//	int seed = 0;
//
//	int nEvents = 10000;
//	int batchSize = 10000;
//	
//	//parse command line:
//	for (int i = 1; i < argc; i++){
//		
//		string arg(argv[i]);
//		
//		if (arg == "-c"){  
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  configfile = argv[++i]; }
//		if (arg == "-o"){  
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  outname = argv[++i]; }
//		if (arg == "-hd"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  hddmname = argv[++i]; }
//		if (arg == "-l"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  lowMass = atof( argv[++i] ); }
//		if (arg == "-u"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  highMass = atof( argv[++i] ); }
//		if (arg == "-n"){  
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  nEvents = atoi( argv[++i] ); }
//		if (arg == "-m"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  beamMaxE = atof( argv[++i] ); }
//		if (arg == "-p"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//			else  beamPeakE = atof( argv[++i] ); }
//		if (arg == "-a"){
//			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//                        else  beamLowE = atof( argv[++i] ); }
//                if (arg == "-b"){
//                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//                        else  beamHighE = atof( argv[++i] ); }
//		if (arg == "-r"){
//                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//                        else  runNum = atoi( argv[++i] ); }
//		if (arg == "-s"){
//                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//                        else  seed = atoi( argv[++i] ); }
//		if (arg == "-d"){
//			diag = true; }
//		if (arg == "-f"){
//			genFlat = true; }
//		if (arg == "-h"){
//			cout << endl << " Usage for: " << argv[0] << endl << endl;
//			cout << "\t -c  <file>\t Config file" << endl;
//			cout << "\t -o  <name>\t ROOT file output name" << endl;
//			cout << "\t -hd <name>\t HDDM file output name [optional]" << endl;
//			cout << "\t -l  <value>\t Low edge of mass range (GeV) [optional]" << endl;
//			cout << "\t -u  <value>\t Upper edge of mass range (GeV) [optional]" << endl;
//			cout << "\t -n  <value>\t Minimum number of events to generate [optional]" << endl;
//			cout << "\t -m  <value>\t Electron beam energy (or photon energy endpoint) [optional]" << endl;
//                        cout << "\t -p  <value>\t Coherent peak photon energy [optional]" << endl;
//                        cout << "\t -a  <value>\t Minimum photon energy to simulate events [optional]" << endl;
//                        cout << "\t -b  <value>\t Maximum photon energy to simulate events [optional]" << endl;
//			cout << "\t -r  <value>\t Run number assigned to generated events [optional]" << endl;
//			cout << "\t -s  <value>\t Random number seed initialization [optional]" << endl;
//			cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
//			cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
//			exit(1);
//		}
//	}
//	
//	if( configfile.size() == 0 || outname.size() == 0 ){
//		cout << "No config file or output specificed:  run gen_2pi -h for help" << endl;
//		exit(1);
//	}
//	
//	// open config file and be sure only one reaction is specified
//	ConfigFileParser parser( configfile );
//	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
//	assert( cfgInfo->reactionList().size() == 1 );
//	ReactionInfo* reaction = cfgInfo->reactionList()[0];
//	
//	// random number initialization (set to 0 by default)
//	gRandom->SetSeed(seed);
//
//	// setup AmpToolsInterface
//	AmpToolsInterface::registerAmplitude( TwoPiAngles() );
//	AmpToolsInterface::registerAmplitude( BreitWigner() );
//	AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );
//	
//	ProductionMechanism::Type type =
//		( genFlat ? ProductionMechanism::kFlat : ProductionMechanism::kResonant );
//	
//	// generate over a range of mass -- the daughters are two charged pions
//	GammaPToXYP resProd( lowMass, highMass, 0.140, 0.140, beamMaxE, beamPeakE, beamLowE, beamHighE, type );
//	
//	// seed the distribution with a sum of noninterfering Breit-Wigners
//	// we can easily compute the PDF for this and divide by that when
//	// doing accept/reject -- improves efficiency if seeds are picked well
//	
//	if( !genFlat ){
//		
//		// the lines below should be tailored by the user for the particular desired
//		// set of amplitudes -- doing so will improve efficiency.  Leaving as is
//		// won't make MC incorrect, it just won't be as fast as it could be
//		
//		resProd.addResonance( 0.775, 0.146,  1.0 );
//	}
//	
//	vector< int > pTypes;
//	pTypes.push_back( Gamma );
//	pTypes.push_back( Proton );
//	pTypes.push_back( PiPlus );
//	pTypes.push_back( PiMinus );
//	
//	HDDMDataWriter* hddmOut = NULL;
//	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, runNum );
//	ROOTDataWriter rootOut( outname );
//	
//	TFile* diagOut = new TFile( "gen_2pi_diagnostic.root", "recreate" );
//	
//	TH1F* mass = new TH1F( "M", "Resonance Mass", 180, lowMass, highMass );
//	TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, lowMass, highMass );
//	massW->Sumw2();
//	TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
//	TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
//	
//	TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);
//	
//	int eventCounter = 0;
//	while( eventCounter < nEvents ){
//		
//		if( batchSize < 1E4 ){
//			
//			cout << "WARNING:  small batches could have batch-to-batch variations\n"
//			     << "          due to different maximum intensities!" << endl;
//		}
//		
//		cout << "Generating four-vectors..." << endl;
//		
//		ati.clearEvents();
//		for( int i = 0; i < batchSize; ++i ){
//			
//			Kinematics* kin = resProd.generate();
//			ati.loadEvent( kin, i, batchSize );
//			delete kin;
//		}
//		
//		cout << "Processing events..." << endl;
//		
//		// include factor of 1.5 to be safe in case we miss peak -- avoid
//		// intensity calculation of we are generating flat data
//		double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
//		
//		
//		for( int i = 0; i < batchSize; ++i ){
//			
//			Kinematics* evt = ati.kinematics( i );
//			TLorentzVector resonance( evt->particle( 2 ) + 
//                                  evt->particle( 3 ) );
//			
//			double genWeight = evt->weight();
//			
//			// cannot ask for the intensity if we haven't called process events above
//			double weightedInten = ( genFlat ? 1 : ati.intensity( i ) );
//			
//			if( !diag ){
//				
//				// obtain this by looking at the maximum value of intensity * genWeight
//				double rand = gRandom->Uniform() * maxInten;
//				
//				if( weightedInten > rand || genFlat ){
//					
//					mass->Fill( resonance.M() );
//					massW->Fill( resonance.M(), genWeight );
//					
//					intenW->Fill( weightedInten );
//					intenWVsM->Fill( resonance.M(), weightedInten );
//					
//					// calculate angular variables
//					TLorentzVector beam = evt->particle ( 0 );
//					TLorentzVector recoil = evt->particle ( 1 );
//					TLorentzVector p1 = evt->particle ( 2 );
//					
//					TLorentzRotation resonanceBoost( -resonance.BoostVector() );
//					
//					TLorentzVector beam_res = resonanceBoost * beam;
//					TLorentzVector recoil_res = resonanceBoost * recoil;
//					TLorentzVector p1_res = resonanceBoost * p1;
//					
//					TVector3 z = -recoil_res.Vect().Unit();
//					TVector3 y = beam_res.Vect().Cross(z).Unit();
//					TVector3 x = y.Cross(z).Unit();
//					
//					TVector3 angles(   (p1_res.Vect()).Dot(x),
//							   (p1_res.Vect()).Dot(y),
//							   (p1_res.Vect()).Dot(z) );
//					
//					GDouble CosTheta = angles.CosTheta();
//					
//					GDouble phi = angles.Phi();
//					GDouble Phi = recoil.Vect().Phi();
//					
//					GDouble psi = phi - Phi;
//					if(psi < -1*PI) psi += 2*PI;
//					if(psi > PI) psi -= 2*PI;
//					
//					CosTheta_psi->Fill( psi, CosTheta);
//					
//					// we want to save events with weight 1
//					evt->setWeight( 1.0 );
//					
//					if( hddmOut ) hddmOut->writeEvent( *evt, pTypes );
//					rootOut.writeEvent( *evt );
//					++eventCounter;
//				}
//			}
//			else{
//				
//				mass->Fill( resonance.M() );
//				massW->Fill( resonance.M(), genWeight );
//				
//				intenW->Fill( weightedInten );
//				intenWVsM->Fill( resonance.M(), weightedInten );
//				TLorentzVector recoil = evt->particle ( 1 );
//				
//				++eventCounter;
//			}
//			
//			delete evt;
//		}
//		
//		cout << eventCounter << " events were processed." << endl;
//	}
//	
//	mass->Write();
//	massW->Write();
//	intenW->Write();
//	intenWVsM->Write();
//	CosTheta_psi->Write();
//	diagOut->Close();
//	
//	if( hddmOut ) delete hddmOut;
	
	return 0;
}

//-----------------------
// GenerateMuPair
//-----------------------
void GenerateMuPair(TVector3 &pgamma, TVector3 &pol, TLorentzVector &pmuplus, TLorentzVector &pmuminus)
{
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

	// mu+ mu- directions for gamma in z-direction
	TVector3 MuPlusDirection  ( sin(thetaPlus) *cos(phi0+phiHalf),
				   sin(thetaPlus)  *sin(phi0+phiHalf), cos(thetaPlus) );
	TVector3 MuMinusDirection (-sin(thetaMinus)*cos(phi0-phiHalf),
				  -sin(thetaMinus) *sin(phi0-phiHalf), cos(thetaMinus) );
	// rotate to actual gamma direction
	MuPlusDirection.RotateUz(GammaDirection);
	MuMinusDirection.RotateUz(GammaDirection);
	
	double Pplus = sqrt(EPlus*EPlus - Mmuon*Mmuon);
	double Pminus = sqrt(EMinus*EMinus - Mmuon*Mmuon);
	
	pmuplus.SetVectM(Pplus*MuPlusDirection, Mmuon);
	pmuminus.SetVectM(Pminus*MuMinusDirection, Mmuon);
	
	
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


