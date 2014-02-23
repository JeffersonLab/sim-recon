
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "HDDM/hddm_s.h"
#include "particleType.h"

char *INPUT_FILE=NULL;
string OUTPUT_FILE("output.hddm");

void ParseCommandLineArguments(int narg,char *argv[]);
int Str2GeantParticleID(char *str);
void Usage(void);
double randm(double, double);

float vertex[4]={0.0, 0.0, 65.0, 65.0};
Particle_t targetType = Proton;
Particle_t beamType = Gamma;
bool FIXED_BEAM_MOMENTUM = false;
float BEAM_MOMENTUM = 8.5;
float BEAM_MOMENTUM_SIGMA = 0.005;

#include <TRandom.h>
TRandom *rnd;

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

	// Create the random generator
	rnd = new TRandom();

	// Open input file
	ifstream *file = new ifstream(INPUT_FILE);
	if(!file->is_open()){
		cerr<<"Unable to open file \""<<INPUT_FILE<<"\" for reading."<<endl;
		exit(-2);
	}
	
	// Open output file
	s_iostream_t* thisOutputStream = init_s_HDDM((char*)OUTPUT_FILE.c_str());
	if(!thisOutputStream){
		cerr<<"Unable to open output file \""<<OUTPUT_FILE<<"\" for writing."<<endl;
		exit(-3);
	}
	
	// Loop over events
	int Nevents = 0;
	while(!file->eof()){
		int runNumber=0, eventNumber=0, nParticles=0;
		(*file)>> runNumber >> eventNumber >> nParticles;
		if(runNumber==0 && eventNumber==0 && nParticles==0)break;
	
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
			int N, charge, type;
			char typestr[256];
			float mass, px, py, pz, E;
			(*file)>> N >> typestr >> mass;
				(*file)>> charge >> px >> py >> pz >> E;

			type = Str2GeantParticleID(typestr);
			if(type<0)type = atoi(typestr);
			
			ps->in[ps->mult].type = (Particle_t)type;
			ps->in[ps->mult].pdgtype = PDGtype((Particle_t)type);
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
		
		// If a specific beam momentum was specified, overwrite
		// the calculated momentum with it.
		if(FIXED_BEAM_MOMENTUM){
			float p = BEAM_MOMENTUM;
			if(BEAM_MOMENTUM_SIGMA!=0.0){
				float delta_p = BEAM_MOMENTUM_SIGMA*rnd->Gaus();
				p += delta_p;
			}
			be->momentum->px = 0.0;
			be->momentum->py = 0.0;
			be->momentum->pz = p;
		}
		
		be->momentum->E = sqrt(SQR(be->properties->mass)+
					SQR(be->momentum->px)+
					SQR(be->momentum->py)+
					SQR(be->momentum->pz)
                                      );
		
		if(nParticles>0){
			flush_s_HDDM(thisOutputEvent, thisOutputStream);
			if(eventNumber%1000 == 0){cout<<"Wrote event "<<eventNumber<<"\r"; cout.flush();}
			Nevents++;
		}
	}
	
	// Close input file
	file->close();

	// Close output file
	close_s_HDDM(thisOutputStream);
	
	cout<<"Wrote "<<Nevents<<" events to "<<OUTPUT_FILE<<endl;
	
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
				case 'P':
					FIXED_BEAM_MOMENTUM = true;
					BEAM_MOMENTUM = atof(&ptr[1]);
					break;
				case 's':
					BEAM_MOMENTUM_SIGMA = atof(&ptr[1])/1000.0;
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
	
	// Determine output filename from input filename
	OUTPUT_FILE = INPUT_FILE;
	size_t pos = OUTPUT_FILE.find_last_of(".");
	if(pos != string::npos) OUTPUT_FILE.erase(pos);
	OUTPUT_FILE += ".hddm";
	
	if(FIXED_BEAM_MOMENTUM){
		cout<<endl;
		cout<<"Using fixed beam: "<<ParticleType(beamType)<<"  P = "<<BEAM_MOMENTUM<<" +/- "<<BEAM_MOMENTUM_SIGMA<<" GeV"<<endl;
		cout<<endl;
	}
	
	
}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"       genr8_2_hddm [options] file.ascii"<<endl;
	cout<<endl;
	cout<<"Convert an ascii file of events generated by genr8 into HDDM"<<endl;
	cout<<"for use as input to hdgeant."<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<endl;
	cout<<"  -V\"x  y  z_min  z_max\"    set the vertex for the interaction."<<endl;
	cout<<"                            (default: x="<<vertex[0]<<" y="<<vertex[1]<<" z_min="<<vertex[2]<<" z_max="<<vertex[3]<<")"<<endl;
	cout<<"  -b\"beam_particle_name\"    set the beam particle type [gamma]."<<endl;
	cout<<"  -t\"target_particle_name\"  set the target particle type [proton]."<<endl;
	cout<<"  -P#                       Set the incident particle momentum in GeV."<<endl;
	cout<<"                            (default: calculate from momentum of"<<endl;
	cout<<"                            final state particles.)"<<endl;
	cout<<"  -s#                       Set the momentum resolution of the beam"<<endl;
	cout<<"                            in MeV. [5MeV]. (Only used if -P option"<<endl;
	cout<<"                            is present.)"<<endl;
	cout<<"  -h                        print this usage statement."<<endl;
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
	if(!strcmp(str, "Pb208"))return Pb208;
	
	return -1;
}

/**************************/
/*  Random generator      */
/*------------------------*/
double randm(double low, double high)
{
  return ((high - low) * rnd->Rndm() + low);
}
