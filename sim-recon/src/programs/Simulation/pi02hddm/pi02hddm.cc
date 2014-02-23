
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "HDDM/hddm_s.h"
#include "particleType.h"

char *INPUT_FILE=NULL;
char OUTPUT_FILE[] = "pi0.hddm";

void ParseCommandLineArguments(int narg,char *argv[]);
int Str2GeantParticleID(char *str);
void Usage(void);
double randm(double, double);

float vertex[4]={0.0, 0.0, 65.0, 65.0};
Particle_t targetType = Proton;
Particle_t beamType = Gamma;

#define SQR(X) ((X)*(X))

time_t now;

//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
	ParseCommandLineArguments(narg,argv);
	
	if(!INPUT_FILE){
		cerr<<"No input file!"<<endl;
	}

	// Seed the random generator
	now=time(NULL);
	srand48(now);

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
	int Nevents = 0;
	while(!file->eof()){
		int runNumber=1, eventNumber=Nevents+1, nParticles=1;
	
		// Start a new event
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
		pes->in[0].runNo = runNumber;
		pes->in[0].eventNo = eventNumber;
		pes->in[0].reactions = rs = make_s_Reactions(1);
		rs->mult = 1;
		rs->in[0].target = ta = make_s_Target();
		ta->type = targetType;
		ta->properties = make_s_Properties();
		ta->properties->charge = ParticleCharge(targetType);
		ta->properties->mass = ParticleMass(targetType);
		ta->momentum = make_s_Momentum();
		ta->momentum->px = 0;
		ta->momentum->py = 0;
		ta->momentum->pz = 0;
		ta->momentum->E  = ParticleMass(targetType);
		rs->in[0].beam = be = make_s_Beam();
		be->type = beamType;
		be->properties = make_s_Properties();
		be->properties->charge = ParticleCharge(beamType);
		be->properties->mass = ParticleMass(beamType);
		be->momentum = make_s_Momentum();
		be->momentum->px = -ta->momentum->px;
		be->momentum->py = -ta->momentum->py;
		be->momentum->pz = -ta->momentum->pz;
		be->momentum->E  = -ta->momentum->E;
		rs->in[0].vertices = vs = make_s_Vertices(1);
		vs->mult = 1;
		vs->in[0].origin = origin = make_s_Origin();
		vs->in[0].products = ps = make_s_Products(nParticles);
		ps->mult = 0;
		
		origin->t = 0.0;
		origin->vx = vertex[0];
		origin->vy = vertex[1];

		if(vertex[2]<vertex[3]){
		  origin->vz = randm(vertex[2],vertex[3]);
		}  else {
		  origin->vz = vertex[2];
		}

		for(int i=0;i<nParticles;i++, ps->mult++){
			int N, charge=0.0, type=7;  // type=7 for pi0
			char typestr[256];
			
			float px, py, pz;
			(*file) >> px >> py >> pz;

			float mass = 0.1349766;
			float p2 = px*px + py*py + pz*pz;
			float E2 = p2 + mass*mass;
			float E = sqrt(E2);

			ps->in[ps->mult].type = (Particle_t)type;
			ps->in[ps->mult].pdgtype = 0;		/* don't bother with the PDG type here */
			ps->in[ps->mult].id = i+1;			/* unique value for this particle within the event */
			ps->in[ps->mult].parentid = 0;	/* All internally generated particles have no parent */
			ps->in[ps->mult].mech = 0;			/* maybe this should be set to something? */
			ps->in[ps->mult].momentum = make_s_Momentum();
			ps->in[ps->mult].momentum->px = px;
			ps->in[ps->mult].momentum->py = py;
			ps->in[ps->mult].momentum->pz = pz;
			ps->in[ps->mult].momentum->E  = E;
			be->momentum->px += px;
			be->momentum->py += py;
			be->momentum->pz += pz;
			be->momentum->E  += E;

		}
		be->momentum->E = sqrt(SQR(be->properties->mass)+
					SQR(be->momentum->px)+
					SQR(be->momentum->py)+
					SQR(be->momentum->pz)
                                      );
		
		if(nParticles>0){
			flush_s_HDDM(thisOutputEvent, thisOutputStream);
			if(eventNumber%1000 == 0)cout<<"Wrote event "<<eventNumber<<endl;
			Nevents++;
		}
	}
	
	// Close input file
	file->close();

	// Close output file
	close_s_HDDM(thisOutputStream);
	
	cout<<"Processed "<<Nevents<<" events"<<endl;
	
	return 0;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
	if(narg<2){
		Usage();
		exit(0);
	}

	for(int i=1;i<narg;i++){
		if(argv[i][0]=='-'){
			char *ptr = &argv[i][1];
			switch(*ptr){
				case 'V':
				  sscanf(&ptr[1], "%f %f %f %f", &vertex[0], &vertex[1], &vertex[2], &vertex[3]);
				  if(vertex[2] > vertex[3]){
				    cerr<<"Invalid parameter: z_min > z_max"<< endl;
				    exit(-1);
				  }
				  break;
				case 'b':
				  beamType = (Particle_t)Str2GeantParticleID(&ptr[1]);
				  break;
				case 't':
				  targetType = (Particle_t)Str2GeantParticleID(&ptr[1]);
				  break;
				default:
				  cerr<<"Unknown option \""<<argv[i]<<"\""<<endl;
				  Usage();
				  exit(-1);
			}
		}else{
			INPUT_FILE = argv[i];
		}
	}
}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"       pi02hddm [options] file.ascii"<<endl;
	cout<<endl;
	cout<<"Convert an ascii file of events containing momentum for single pi0"<<endl;
	cout<<" particles to HDDM for use as input to hdgeant."<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<endl;
	cout<<"  -V\"x  y  z_min  z_max\"    set the vertex for the interaction.";
	cout<<"(default: x="<<vertex[0]<<" y="<<vertex[1]<<" z_min="<<vertex[2]<<" z_max="<<vertex[3]<<")"<<endl;
	cout<<"  -b\"beam_particle_name\"    set the beam particle type [gamma]."<<endl;
	cout<<"  -t\"target_particle_name\"  set the target particle type [proton]."<<endl;
	cout<<"  -h           print this usage statement."<<endl;
	cout<<endl;
	cout<<" The input ASCII file should contain 3 numbers per line corresponding"<<endl;
	cout<<"to the momentum components of the pi0: px, py, and pz. The numbers"<<endl;
	cout<<"should be separated by spaces."<<endl;
	cout<<endl;
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

/**************************/
/*  Random generator      */
/*------------------------*/
double randm(double low, double high)
{
  return ((high - low) * drand48() + low);
}
