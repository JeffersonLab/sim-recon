#include <string>
#include <vector>
using std::string;
using std::vector;

#include <TLorentzVector.h>
#include <TVector3.h>

#include "HDDM/hddm_s.hpp"

class Particle{
   public:
      Particle() : type(Unknown), pdgtype(0), parentid(0), mech(0), decayVertex(0)
      {
         p.SetXYZT(0,0,0,0);
      }
   
      Particle_t type;   // Note: Particle_t is an int as opposed to Particle which is a class
      TLorentzVector p;
      
      // The following are not used for beam or target particles
      int pdgtype;
      int parentid;
      int mech;
      int decayVertex;
};

class Event{
   public:
      Event() : runNo(0), eventNo(0), reaction_type(0), reaction_weight(0)
      {
         vertex.SetXYZ(0,0,0);
      }
   
      int runNo;
      int eventNo;
      int reaction_type;       // copied to HDDM, but not used (as far as I know)
      int reaction_weight;     // copied to HDDM, but not used (as far as I know)
      TVector3 vertex;         // Set this to 0,0,0 to have GEANT distribute it in target
      Particle beam;
      Particle target;
      vector<Particle> intermediate;   // Not to be tracked by GEANT
      vector<Particle> final;          // To be tracked by GEANT
};

typedef hddm_s::ProductList::iterator PlistIter;

void open_hddm_output(string fname);
void close_hddm_output(void);
void write_hddm_event(Event &event);
void CopyParticleToProduct(int id, const Particle &part, PlistIter &prod);
Particle_t PDG_to_GEANT(int pdgtype);
int GEANT_to_PDG(Particle_t type);
void InitializePDGGEANTmap(void);
