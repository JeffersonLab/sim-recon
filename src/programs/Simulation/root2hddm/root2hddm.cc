#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

#include <stdlib.h>

#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "HDDM/hddm_s.hpp"
#include "particleType.h"

char *INPUT_FILE=NULL;
char *OUTPUT_FILE=NULL;
char *MAX_EVENTS=NULL;
string TREENAME;
vector<int> PTYPE;

void ParseCommandLineArguments(int narg,char *argv[]);
void Usage(void);
double randm(double, double);

float vertex[4] = {0.0, 0.0, 65.0, 65.0};

time_t now;

//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
  //some defaults
  TREENAME="kin";
  char defOutFName[]="output.hddm";
  char defMaxEvents[]="10000000000";
  OUTPUT_FILE=defOutFName;
  MAX_EVENTS=defMaxEvents;
  
  ParseCommandLineArguments(narg,argv);

  unsigned int NmaxEvents=atoi(MAX_EVENTS);

  // Seed the random generator
  now=time(NULL);
  srand48(now);
  
  ifstream filetest(INPUT_FILE);
  if (! filetest.good()) {
    std::cerr << "Unable to open input file \"" << INPUT_FILE
              << "\" for reading." << std::endl;
    exit(-4);
  }

  vector<string> readerArgs;
  readerArgs.push_back(INPUT_FILE);
  readerArgs.push_back(TREENAME);
  
  ROOTDataReader reader(readerArgs);
  NmaxEvents = (NmaxEvents < reader.numEvents())? 
                NmaxEvents : reader.numEvents();
   
  // Open output file
  std::ofstream *thisOutputFile = new std::ofstream(OUTPUT_FILE);
  if (! thisOutputFile->is_open()) {
    std::cerr << "Unable to open output file \"" << OUTPUT_FILE
              << "\" for writing." << std::endl;
    exit(-3);
  }
  hddm_s::ostream *thisOutputStream = new hddm_s::ostream(*thisOutputFile);
  
  // Loop over events
  const Kinematics *kin;

  int runNumber=9000, Nevents=0;
  
  for (unsigned int eventNumber=0;
       eventNumber < NmaxEvents ; ++eventNumber) {
    // get event
    kin=reader.getEvent();
    int nParticles=kin->particleList().size();
    
    // Start a new event in the HDDM record
    hddm_s::HDDM record;
    hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
    pes().setRunNo(runNumber);
    pes().setEventNo(eventNumber);
    hddm_s::ReactionList rs = pes().addReactions();
    hddm_s::VertexList vs = rs().addVertices();
    hddm_s::OriginList os = vs().addOrigins();
    hddm_s::ProductList ps = vs().addProducts(nParticles-1);
    rs().setWeight(kin->weight());
    hddm_s::BeamList bs = rs().addBeams();
    bs().setType((Particle_t)1);
    hddm_s::MomentumList bmoms = bs().addMomenta();
    bmoms().setPx(kin->particle(0).px());
    bmoms().setPy(kin->particle(0).py());
    bmoms().setPz(kin->particle(0).pz());
    bmoms().setE(kin->particle(0).e());
    hddm_s::PropertiesList bpros = bs().addPropertiesList();
    bpros().setCharge(0);
    bpros().setMass(0.0);

    hddm_s::TargetList ts = rs().addTargets();
    ts().setType((Particle_t)14);
    hddm_s::MomentumList tmoms = ts().addMomenta();
    tmoms().setPx(0);
    tmoms().setPy(0);
    tmoms().setPz(0);
    tmoms().setE(0.938272);
    hddm_s::PropertiesList tpros = ts().addPropertiesList();
    tpros().setCharge(+1);
    tpros().setMass(0.938272);

    os().setT(0.0);
    os().setVx(vertex[0]);
    os().setVy(vertex[1]);
    if (vertex[2] < vertex[3])
      os().setVz(randm(vertex[2],vertex[3]));
    else 
      os().setVz(vertex[2]);
    
    for (int i=1; i < nParticles; i++)
    {  
      ps(i-1).setType((Particle_t)PTYPE[i]);
      ps(i-1).setPdgtype(0);    /* don't bother with the PDG type here */
      ps(i-1).setId(i);         /* unique value for this particle within the event */
      ps(i-1).setParentid(0);   /* All internally generated particles have no parent */
      ps(i-1).setMech(0);       /* maybe this should be set to something? */
      hddm_s::MomentumList pmoms = ps(i-1).addMomenta();
      pmoms().setPx(kin->particle(i).px());
      pmoms().setPy(kin->particle(i).py());
      pmoms().setPz(kin->particle(i).pz());
      pmoms().setE(kin->particle(i).e());
    }
    
    if (nParticles > 0) {
      *thisOutputStream << record;
      if (eventNumber%1000 == 0)
         std::cout << "Wrote event " << eventNumber << std::endl;
      Nevents++;
    }
  }
  
  // Close output file
  delete thisOutputStream;
  delete thisOutputFile;
  
  std::cout << "Processed " << Nevents << " events" << std::endl;
  
  return 0;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
  if (narg<2) {
    Usage();
    exit(0);
  }
  char *ptypestr;
  
  for (int i=1; i<narg; i++) {
    if (argv[i][0] == '-') {
      char *ptr = &argv[i][1];
      switch(*ptr) {
      case 'V':
   sscanf(&ptr[1], "%f %f %f %f", &vertex[0], &vertex[1], &vertex[2], &vertex[3]);
   if (vertex[2] > vertex[3]) {
     std::cerr << "Invalid parameter: z_min > z_max" << std::endl;
     exit(-1);
   }
   break;
      case 'P':
   ptypestr=strtok(&ptr[1]," ");
   while (ptypestr != NULL) {
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
   std::cerr << "Unknown option \"" << argv[i] << "\"" << std::endl;
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
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "       root2hddm [options] file.root" << std::endl;
  std::cout << std::endl;
  std::cout << "Convert an ROOT file of events from an AmpTools-based generator" << std::endl;
  std::cout << "(4-vectors as \"Kinematics\" objects) for use as input to hdgeant." << std::endl;
  std::cout << std::endl;
  std::cout << " options:" << std::endl;
  std::cout << std::endl;
  std::cout << "  -P\"ptype1 ptype2 ...\"     set particle types starting from the photon..." << std::endl;
  std::cout << "  -o[output_file_name]      name of output file (default: output.hddm)" << std::endl;
  std::cout << "  -T[tree_name]             set name of tree in ROOT file default: kin)" << std::endl;
  std::cout << "  -V\"x  y  z_min  z_max\"    set the vertex for the interaction." << std::endl;
  std::cout << "(default: x=" << vertex[0] << " y=" << vertex[1] << " z_min=" << vertex[2] << " z_max=" << vertex[3] << ")" << std::endl;
  std::cout << "  -h           print this usage statement." << std::endl;
  std::cout << std::endl;
}


/**************************/
/*  Random generator      */
/*------------------------*/
double randm(double low, double high)
{
  return ((high - low) * drand48() + low);
}
