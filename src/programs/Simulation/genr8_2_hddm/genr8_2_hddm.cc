#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "HDDM/hddm_s.hpp"
#include "particleType.h"

char *INPUT_FILE=NULL;
string OUTPUT_FILE("output.hddm");

void ParseCommandLineArguments(int narg,char *argv[]);
int Str2GeantParticleID(char *str);
void Usage(void);
double randm(double, double);

float vertex[4] = {0.0, 0.0, 65.0, 65.0};
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
   
   if (!INPUT_FILE) {
      std::cerr << "No input file!" << std::endl;
   }

   // Create the random generator
   rnd = new TRandom();

   // Open input file
   ifstream *infile = new ifstream(INPUT_FILE);
   if (! infile->is_open()) {
      std::cerr << "Unable to open file \"" << INPUT_FILE << "\" for reading."
                << std::endl;
      exit(-2);
   }

   // Open output file

   std::ofstream *outfile = new ofstream(OUTPUT_FILE.c_str());
   if (! outfile->is_open()) {
      std::cerr << "Unable to open output file \"" << OUTPUT_FILE
                << "\" for writing." << std::endl;
      exit(-3);
   }
   hddm_s::ostream *outstream = new hddm_s::ostream(*outfile);

   // Loop over events
   int Nevents = 0;
   while (! infile->eof()) {
      int runNumber=0, eventNumber=0, nParticles=0;
      *infile >> runNumber >> eventNumber >> nParticles;
      if (runNumber == 0 && eventNumber == 0 && nParticles == 0)
         break;
   
      // Start a new event
      hddm_s::HDDM record;
      hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
      pes().setRunNo(runNumber);
      pes().setEventNo(eventNumber);
      hddm_s::ReactionList rs = pes().addReactions();
      hddm_s::TargetList ts = rs().addTargets();
      ts().setType(targetType);
      hddm_s::PropertiesList tpros = ts().addPropertiesList();
      tpros().setCharge(ParticleCharge(targetType));
      tpros().setMass(ParticleMass(targetType));
      hddm_s::MomentumList tmoms = ts().addMomenta();
      tmoms().setPx(0);
      tmoms().setPy(0);
      tmoms().setPz(0);
      tmoms().setE(ParticleMass(targetType));
      hddm_s::BeamList bs = rs().addBeams();
      bs().setType(beamType);
      hddm_s::PropertiesList bpros = bs().addPropertiesList();
      bpros().setCharge(ParticleCharge(beamType));
      bpros().setMass(ParticleMass(beamType));
      hddm_s::MomentumList bmoms = bs().addMomenta();
      bmoms().setPx(-tmoms().getPx());
      bmoms().setPy(-tmoms().getPy());
      bmoms().setPz(-tmoms().getPz());
      bmoms().setE(-tmoms().getE());
      hddm_s::VertexList vs = rs().addVertices();
      hddm_s::OriginList os = vs().addOrigins();
      hddm_s::ProductList ps = vs().addProducts(nParticles);
      
      os().setT(0.0);
      os().setVx(vertex[0]);
      os().setVy(vertex[1]);

      if (vertex[2] < vertex[3]) {
        os().setVz(randm(vertex[2],vertex[3]));
      }
      else {
        os().setVz(vertex[2]);
      }

      for (int i=0; i < nParticles; i++) {
         int N, charge, type;
         char typestr[256];
         float mass, px, py, pz, E;
         *infile >> N >> typestr >> mass;
         *infile >> charge >> px >> py >> pz >> E;

         type = Str2GeantParticleID(typestr);
         if (type < 0)
            type = atoi(typestr);
         
         ps(i).setType((Particle_t)type);
         ps(i).setPdgtype(PDGtype((Particle_t)type));
         ps(i).setId(i+1);         /* unique value for this particle within the event */
         ps(i).setParentid(0);     /* All internally generated particles have no parent */
         ps(i).setMech(0);         /* maybe this should be set to something? */
         hddm_s::MomentumList pmoms = ps(i).addMomenta();
         pmoms().setPx(px);
         pmoms().setPy(py);
         pmoms().setPz(pz);
         pmoms().setE(E);
         bmoms().setPx(bmoms().getPx() + px);
         bmoms().setPy(bmoms().getPy() + py);
         bmoms().setPz(bmoms().getPz() + pz);
         bmoms().setE(bmoms().getE() + E);
      }
      
      // If a specific beam momentum was specified, overwrite
      // the calculated momentum with it.
      if (FIXED_BEAM_MOMENTUM) {
         float p = BEAM_MOMENTUM;
         if (BEAM_MOMENTUM_SIGMA!=0.0) {
            float delta_p = BEAM_MOMENTUM_SIGMA*rnd->Gaus();
            p += delta_p;
         }
         bmoms().setPx(0.0);
         bmoms().setPy(0.0);
         bmoms().setPz(p);
      }
      
      bmoms().setE(sqrt(SQR(bpros().getMass()) + SQR(bmoms().getPx()) +
                        SQR(bmoms().getPy()) + SQR(bmoms().getPz())));
      
      if (nParticles > 0) {
         *outstream << record;
         if (eventNumber%1000 == 0) {
            std::cout << "Wrote event " << eventNumber << "\r";
            std::cout.flush();
         }
         Nevents++;
      }
   }
   
   // Close input file
   delete infile;

   // Close output file
   delete outstream;
   delete outfile;
   
   std::cout << "Wrote " << Nevents << " events to " << OUTPUT_FILE
             << std::endl;
   
   return 0;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
   if (narg < 2) {
      Usage();
      exit(0);
   }

   for(int i=1; i < narg; i++) {
      if (argv[i][0]=='-') {
         char *ptr = &argv[i][1];
         switch(*ptr) {
            case 'V':
              sscanf(&ptr[1], "%f %f %f %f", &vertex[0], &vertex[1],
                                             &vertex[2], &vertex[3]);
              if (vertex[2] > vertex[3]) {
                std::cerr << "Invalid parameter: z_min > z_max" << std::endl;
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
              std::cerr << "Unknown option \"" << argv[i] << "\"" << std::endl;
              Usage();
              exit(-1);
         }
      }
      else {
         INPUT_FILE = argv[i];
      }
   }
   
   // Determine output filename from input filename
   OUTPUT_FILE = INPUT_FILE;
   size_t pos = OUTPUT_FILE.find_last_of(".");
   if (pos != string::npos) OUTPUT_FILE.erase(pos);
   OUTPUT_FILE += ".hddm";
   
   if (FIXED_BEAM_MOMENTUM) {
     std::cout << std::endl;
     std::cout << "Using fixed beam: " << ParticleType(beamType)
               << "  P = " << BEAM_MOMENTUM 
               << " +/- " << BEAM_MOMENTUM_SIGMA << " GeV" << std::endl;
     std::cout << std::endl;
   }
}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "       genr8_2_hddm [options] file.ascii" << std::endl;
  std::cout << std::endl;
  std::cout << "Convert an ascii file of events generated by genr8 into HDDM"
            << std::endl;
  std::cout << "for use as input to hdgeant." << std::endl;
  std::cout << std::endl;
  std::cout << " options:" << std::endl;
  std::cout << std::endl;
  std::cout << "  -V\"x  y  z_min  z_max\"    set the vertex "
               "for the interaction." << std::endl;
  std::cout << "                            (default: x=" << vertex[0]
            << " y=" << vertex[1] << " z_min=" << vertex[2] 
            << " z_max=" << vertex[3] << ")" << std::endl;
  std::cout << "  -b\"beam_particle_name\"    "
               "set the beam particle type [gamma]." << std::endl;
  std::cout << "  -t\"target_particle_name\"  "
               "set the target particle type [proton]." << std::endl;
  std::cout << "  -P#                       "
               "Set the incident particle momentum in GeV." << std::endl;
  std::cout << "                            "
               "(default: calculate from momentum of" << std::endl;
  std::cout << "                            "
               "final state particles.)" << std::endl;
  std::cout << "  -s#                       "
               "Set the momentum resolution of the beam" << std::endl;
  std::cout << "                            in MeV. [5MeV]."
               " (Only used if -P option" << std::endl;
  std::cout << "                            is present.)"
            << std::endl;
  std::cout << "  -h                        "
               "print this usage statement." << std::endl;
  std::cout << std::endl;
}

//-------------------------------
// Str2GeantParticleID
//-------------------------------
int Str2GeantParticleID(char *str)
{
   if (strcmp(str, "unknown") == 0 || strcmp(str, "Unknown") == 0)
      return Unknown;
   if (strcmp(str, "gamma") == 0 || strcmp(str, "Gamma") == 0)
      return Gamma;
   if (strcmp(str, "positron") == 0 || strcmp(str, "Positron") == 0)
      return Positron;
   if (strcmp(str, "electron") == 0 || strcmp(str, "Electron") == 0)
      return Electron;
   if (strcmp(str, "neutrino") == 0 || strcmp(str, "Neutrino") == 0)
      return Neutrino;
   if (strcmp(str, "mu+") == 0 || strcmp(str, "Mu+") == 0)
      return MuonPlus;
   if (strcmp(str, "mu-") == 0 || strcmp(str, "Mu-") == 0)
      return MuonMinus;
   if (strcmp(str, "pi0") == 0 || strcmp(str, "Pi0") == 0)
      return Pi0;
   if (strcmp(str, "pi+") == 0 || strcmp(str, "Pi+") == 0)
      return PiPlus;
   if (strcmp(str, "pi-") == 0 || strcmp(str, "Pi-") == 0)
      return PiMinus;
   if (strcmp(str, "KL") == 0)
      return KLong;
   if (strcmp(str, "K+") == 0)
      return KPlus;
   if (strcmp(str, "K-") == 0)
      return KMinus;
   if (strcmp(str, "neutron") == 0 || strcmp(str, "Neutron") == 0)
      return Neutron;
   if (strcmp(str, "proton") == 0 || strcmp(str, "Proton") == 0)
      return Proton;
   if (strcmp(str, "pbar") == 0 || strcmp(str, "Pbar") == 0)
      return AntiProton;
   if (strcmp(str, "Ks") == 0)
      return KShort;
   if (strcmp(str, "eta") == 0 || strcmp(str, "Eta") == 0)
      return Eta;
   if (strcmp(str, "lambda") == 0 || strcmp(str, "Lambda") == 0)
      return Lambda;
   if (strcmp(str, "sigma+") == 0 || strcmp(str, "Sigma+") == 0)
      return SigmaPlus;
   if (strcmp(str, "sigma0") == 0 || strcmp(str, "Sigma0") == 0)
      return Sigma0;
   if (strcmp(str, "sigma-") == 0 || strcmp(str, "Sigma-") == 0)
      return SigmaMinus;
   if (strcmp(str, "Xi0") == 0)
      return Xi0;
   if (strcmp(str, "Xi-") == 0)
      return XiMinus;
   if (strcmp(str, "omega-") == 0 || strcmp(str, "Omega-") == 0)
      return OmegaMinus;
   if (strcmp(str, "nbar") == 0 || strcmp(str, "Nbar") == 0)
      return AntiNeutron;
   if (strcmp(str, "lambdabar") == 0 || strcmp(str, "Lambdabar") == 0)
      return AntiLambda;
   if (strcmp(str, "sigmabar-") == 0)
      return AntiSigmaMinus;
   if (strcmp(str, "sigmabar0") == 0 || strcmp(str, "Sigmabar0") == 0)
      return AntiSigma0;
   if (strcmp(str, "sigmabar+") == 0 || strcmp(str, "Sigmabar+") == 0)
      return AntiSigmaPlus;
   if (strcmp(str, "Xibar0") == 0)
      return AntiXi0;
   if (strcmp(str, "Xibar+") == 0)
      return AntiXiPlus;
   if (strcmp(str, "omegabar+") == 0 || strcmp(str, "Omegabar+") == 0)
      return AntiOmegaPlus;
   if (strcmp(str, "rho0") == 0 || strcmp(str, "Rho0") == 0)
      return Rho0;
   if (strcmp(str, "rho+") == 0 || strcmp(str, "Rho+") == 0)
      return RhoPlus;
   if (strcmp(str, "rho") == 0 || strcmp(str, "Rho") == 0)
      return RhoMinus;
   if (strcmp(str, "omega") == 0 || strcmp(str, "Omega") == 0)
      return omega;
   if (strcmp(str, "etaprime") == 0 || strcmp(str, "Etaprime") == 0)
      return EtaPrime;
   if (strcmp(str, "phi") == 0 || strcmp(str, "Phi") == 0)
      return phiMeson;
   if (!strcmp(str, "Pb208"))
      return Pb208;
   
   return -1;
}

/**************************/
/*  Random generator      */
/*------------------------*/
double randm(double low, double high)
{
  return ((high - low) * rnd->Rndm() + low);
}
