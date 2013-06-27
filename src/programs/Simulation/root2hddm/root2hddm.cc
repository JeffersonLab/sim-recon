#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "HDDM/hddm_s.h"
#include "particleType.h"

char *INPUT_FILE=NULL;
char *OUTPUT_FILE=NULL;
char *MAX_EVENTS=NULL;
string TREENAME;
vector<int> PTYPE;

void ParseCommandLineArguments(int narg,char *argv[]);
void Usage(void);
double randm(double, double);

float vertex[4]={0.0, 0.0, 65.0, 65.0};

time_t now;

//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
  //some defaults
  TREENAME="kin";
  char defOutFName[]="output.hddm",defMaxEvents[]="10000000000";
  OUTPUT_FILE=defOutFName;
  MAX_EVENTS=defMaxEvents;
  
  ParseCommandLineArguments(narg,argv);

  unsigned int NmaxEvents=atoi(MAX_EVENTS);

  // Seed the random generator
  now=time(NULL);
  srand48(now);
  
  ifstream filetest(INPUT_FILE);
  if(!filetest.good()){
    cerr<<"Unable to open input file \""<< INPUT_FILE <<"\" for reading."<<endl;
    exit(-4);
  }

  vector< string > readerArgs;
  readerArgs.push_back( INPUT_FILE );
  readerArgs.push_back( TREENAME );
  
  ROOTDataReader reader( readerArgs );
  NmaxEvents = NmaxEvents < reader.numEvents() ? NmaxEvents : reader.numEvents(); 
	
  // Open output file
  s_iostream_t* thisOutputStream = init_s_HDDM(OUTPUT_FILE);
  if(!thisOutputStream){
    cerr<<"Unable to open output file \""<<OUTPUT_FILE<<"\" for writing."<<endl;
    exit(-3);
  }
  
  // Loop over events
  const Kinematics *kin;

  int runNumber=9000, Nevents=0;
  
  for (unsigned int eventNumber=0;
       eventNumber < NmaxEvents ; ++eventNumber){
    // get event
    kin=reader.getEvent();
    int nParticles=kin->particleList().size();
    
    // Start a new event in the HDDM record
    s_PhysicsEvents_t* pes;
    s_Reactions_t* rs;
    s_Vertices_t* vs;
    s_Beam_t* bs;
    s_Target_t* ts;
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
    rs->in[0].weight = kin->weight();

    rs->in[0].beam = bs = make_s_Beam();
    bs->type = (Particle_t)1;
    bs->momentum = make_s_Momentum();
    bs->momentum->px = kin->particle(0).px();
    bs->momentum->py = kin->particle(0).py();
    bs->momentum->pz = kin->particle(0).pz();
    bs->momentum->E  = kin->particle(0).e();
    bs->properties = make_s_Properties();
    bs->properties->charge = 0;
    bs->properties->mass = 0.0;

    rs->in[0].target = ts = make_s_Target();
    ts->type = (Particle_t)14;
    ts->momentum = make_s_Momentum();
    ts->momentum->px = 0;
    ts->momentum->py = 0;
    ts->momentum->pz = 0;
    ts->momentum->E  = 0.938272;
    ts->properties = make_s_Properties();
    ts->properties->charge=+1;
    ts->properties->mass = 0.938272;

    origin->t = 0.0;
    origin->vx = vertex[0];
    origin->vy = vertex[1];
    
    if(vertex[2]<vertex[3])
      origin->vz = randm(vertex[2],vertex[3]);
    else 
      origin->vz = vertex[2];
    
    
    for(int i=1;i<nParticles;i++, ps->mult++){
      
      ps->in[ps->mult].type = (Particle_t)PTYPE[i];
      ps->in[ps->mult].pdgtype = 0;    /* don't bother with the PDG type here */
      ps->in[ps->mult].id = i+1;       /* unique value for this particle within the event */
      ps->in[ps->mult].parentid = 0;   /* All internally generated particles have no parent */
      ps->in[ps->mult].mech = 0;       /* maybe this should be set to something? */
      ps->in[ps->mult].momentum = make_s_Momentum();
      ps->in[ps->mult].momentum->px = kin->particle(i).px();
      ps->in[ps->mult].momentum->py = kin->particle(i).py();
      ps->in[ps->mult].momentum->pz = kin->particle(i).pz();
      ps->in[ps->mult].momentum->E  = kin->particle(i).e();
      
    }
    
    if(nParticles>0){
      flush_s_HDDM(thisOutputEvent, thisOutputStream);
      if(eventNumber%1000 == 0)cout<<"Wrote event "<<eventNumber<<endl;
      Nevents++;
    }
  }
  
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
  char *ptypestr;
  
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
      case 'P':
	ptypestr=strtok(&ptr[1]," ");
	while(ptypestr!=NULL){
	  PTYPE.push_back(atoi(ptypestr));
	  ptypestr=strtok(NULL," ");
	}
	break;
      case 'T':
	TREENAME=&ptr[1];
	break;
      case 'o':
	OUTPUT_FILE=&ptr[1];
	break;
      case 'n':
	MAX_EVENTS=&ptr[1];
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
  cout<<"       root2hddm [options] file.root"<<endl;
  cout<<endl;
  cout<<"Convert an ROOT file of events from an AmpTools-based generator"<<endl;
  cout<<"(4-vectors as \"Kinematics\" objects) for use as input to hdgeant."<<endl;
  cout<<endl;
  cout<<" options:"<<endl;
  cout<<endl;
  cout<<"  -P\"ptype1 ptype2 ...\"     set particle types starting from the photon..." <<endl;
  cout<<"  -o[output_file_name]      name of output file (default: output.hddm)"<<endl;
  cout<<"  -T[tree_name]             set name of tree in ROOT file default: kin)"<<endl;
  cout<<"  -V\"x  y  z_min  z_max\"    set the vertex for the interaction."<<endl;
  cout<<"(default: x="<<vertex[0]<<" y="<<vertex[1]<<" z_min="<<vertex[2]<<" z_max="<<vertex[3]<<")"<<endl;
  cout<<"  -h           print this usage statement."<<endl;
  cout<<endl;
}


/**************************/
/*  Random generator      */
/*------------------------*/
double randm(double low, double high)
{
  return ((high - low) * drand48() + low);
}
