
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "hddm_s.h"
#include "particleType.h"

char *INPUT_FILE=NULL;
char *OUTPUT_FILE = "output.hddm";

void ParseCommandLineArguments(int narg,char *argv[]);
int Str2GeantParticleID(char *str);

//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
	ParseCommandLineArguments(narg,argv);
	
	if(!INPUT_FILE){
		cerr<<"No input file!"<<endl;
	}
	
	// Open input file
	ifstream *file = new ifstream(INPUT_FILE);
	if(!file->is_open()){
		cerr<<"Unable to open file \""<<INPUT_FILE<<"\" for reading."<<endl;
		exit(-2);
	}
	
	// Open output file
	s_iostream_t* thisOutputStream = init_s_HDDM(OUTPUT_FILE);
	if(!thisOutputStream){
		cerr<<"Unable to open output file \""<<OUTPUT_FILE<<"\" for writing."<<endl;
		exit(-3);
	}
	
	// Loop over events
	while(!file->eof()){
		int runNumber=0, eventNumber=0, nParticles=0;
		(*file)>> runNumber >> eventNumber >> nParticles;
	
		// Start a new event
		s_PhysicsEvents_t* pes;
		s_Reactions_t* rs;
		s_Vertices_t* vs;
		s_Origin_t* origin;
		s_Products_t* ps;
		s_HDDM_t *thisOutputEvent = make_s_HDDM();
      thisOutputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
		pes->mult = 1;
		pes->in[0].runNo = runNumber;
		pes->in[0].eventNo = eventNumber;
		pes->in[0].reactions = rs = make_s_Reactions(1);
		rs->mult = 1;
		rs->in[0].vertices = vs = make_s_Vertices(1);
		vs->mult = 1;
		vs->in[0].origin = origin = make_s_Origin();
		vs->in[0].products = ps = make_s_Products(nParticles);
		ps->mult = 0;
		
		for(int i=0;i<nParticles;i++, ps->mult++){
			int N, charge, type;
			char typestr[256];
			float mass, px, py, pz, E;
			(*file)>> N >> typestr >> mass;
				(*file)>> charge >> px >> py >> pz >> E;

			type = Str2GeantParticleID(typestr);
			if(type<0)type = atoi(typestr);
			
			ps->in[ps->mult].type = (Particle_t)type;
			ps->in[ps->mult].momentum = make_s_Momentum();
			ps->in[ps->mult].momentum->px = px;
			ps->in[ps->mult].momentum->py = py;
			ps->in[ps->mult].momentum->pz = pz;
			ps->in[ps->mult].momentum->E  = E;

		}
		
		if(nParticles>0){
			flush_s_HDDM(thisOutputEvent, thisOutputStream);
			if(eventNumber%1000 == 0)cout<<"Wrote event "<<eventNumber<<endl;
		}
	}
	
	// Close input file
	file->close();

	// Close output file
	close_s_HDDM(thisOutputStream);
	
	return 0;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
	for(int i=1;i<narg;i++){
		if(argv[i][0]=='-'){
			char *ptr = &argv[i][1];
			switch(*ptr){
			
				default:
					cerr<<"Unknown option \""<<argv[i]<<"\""<<endl;
					exit(-1);
			}
		}else{
			INPUT_FILE = argv[i];
		}
	}
}

//-------------------------------
// Str2GeantParticleID
//-------------------------------
int Str2GeantParticleID(char *str)
{
	if(!strcmp(str, "unknown"))return Unknown;
	if(!strcmp(str, "gamma"))return Gamma;
	if(!strcmp(str, "positron"))return Positron;
	if(!strcmp(str, "electron"))return Electron;
	if(!strcmp(str, "neutrino"))return Neutrino;
	if(!strcmp(str, "mu+"))return MuonPlus;
	if(!strcmp(str, "mu-"))return MuonMinus;
	if(!strcmp(str, "pi0"))return Pi0;
	if(!strcmp(str, "pi+"))return PiPlus;
	if(!strcmp(str, "pi-"))return PiMinus;
	if(!strcmp(str, "KL"))return KLong;
	if(!strcmp(str, "K+"))return KPlus;
	if(!strcmp(str, "K-"))return KMinus;
	if(!strcmp(str, "neutron"))return Neutron;
	if(!strcmp(str, "proton"))return Proton;
	if(!strcmp(str, "pbar"))return AntiProton;
	if(!strcmp(str, "Ks"))return KShort;
	if(!strcmp(str, "eta"))return Eta;
	if(!strcmp(str, "lambda"))return Lambda;
	if(!strcmp(str, "sigma+"))return SigmaPlus;
	if(!strcmp(str, "sigma0"))return Sigma0;
	if(!strcmp(str, "sigma-"))return SigmaMinus;
	if(!strcmp(str, "Xi0"))return Xi0;
	if(!strcmp(str, "Xi-"))return XiMinus;
	if(!strcmp(str, "omega-"))return OmegaMinus;
	if(!strcmp(str, "nbar"))return AntiNeutron;
	if(!strcmp(str, "lambdabar"))return AntiLambda;
	if(!strcmp(str, "sigmabar-"))return AntiSigmaMinus;
	if(!strcmp(str, "sigmabar0"))return AntiSigma0;
	if(!strcmp(str, "sigmabar+"))return AntiSigmaPlus;
	if(!strcmp(str, "Xibar0"))return AntiXi0;
	if(!strcmp(str, "Xibar+"))return AntiXiPlus;
	if(!strcmp(str, "omegabar+"))return AntiOmegaPlus;
	if(!strcmp(str, "rho0"))return Rho0;
	if(!strcmp(str, "rho+"))return RhoPlus;
	if(!strcmp(str, "rho"))return RhoMinus;
	if(!strcmp(str, "omega"))return omega;
	if(!strcmp(str, "etaprime"))return EtaPrime;
	if(!strcmp(str, "phi"))return phiMeson;
	
	return -1;
}


