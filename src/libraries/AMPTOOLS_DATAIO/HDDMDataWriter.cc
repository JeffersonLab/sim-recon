#include "HDDMDataWriter.h"
#include "HDDM/hddm_s.h"



HDDMDataWriter::HDDMDataWriter( const string& outFile, int runNumber)
{
  m_OutputStream = init_s_HDDM((char*)outFile.c_str());
  m_runNumber=runNumber;
  
  m_eventCounter = 0;
}

HDDMDataWriter::~HDDMDataWriter()
{
  close_s_HDDM(m_OutputStream);
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, vector<int> ptype, bool centeredVertex)
{
  if(centeredVertex) writeEvent(kin,ptype,0,0,65 /*cm*/);
  else writeEvent(kin,ptype,0,0,50 /*cm*/,80/*cm*/);
}

void HDDMDataWriter::
writeEvent( const Kinematics& kin, vector<int> ptype, 
	    float vx, float vy, float vz_min, float vz_max)
{
  if(vz_min>vz_max){
    float tmp=vz_min;
    vz_min=vz_max;
    vz_max=tmp;
  }
  writeEvent(kin, ptype,vx, vy, (vz_max - vz_min) * drand48() + vz_min);
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, vector<int> ptype, 
	    float vx, float vy, float vz)
{
  vector< HepLorentzVector > particleList = kin.particleList();
  int nParticles=kin.particleList().size();
  
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
  pes->in[0].runNo = m_runNumber;
  pes->in[0].eventNo = m_eventCounter;
  pes->in[0].reactions = rs = make_s_Reactions(1);
  rs->mult = 1;
  rs->in[0].vertices = vs = make_s_Vertices(1);
  vs->mult = 1;
  vs->in[0].origin = origin = make_s_Origin();
  vs->in[0].products = ps = make_s_Products(nParticles);
  ps->mult = 0;
  
  origin->t = 0.0;
  origin->vx = vx;
  origin->vy = vy;
  origin->vz = vz;
  

  rs->in[0].beam = bs = make_s_Beam();
  bs->type = (Particle_t)1;
  bs->momentum = make_s_Momentum();
  bs->momentum->px = kin.particle(0).px();
  bs->momentum->py = kin.particle(0).py();
  bs->momentum->pz = kin.particle(0).pz();
  bs->momentum->E  = kin.particle(0).e();
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

  
  for(int i=1;i<nParticles;i++, ps->mult++){
    
      ps->in[ps->mult].type = (Particle_t)ptype[i];
      ps->in[ps->mult].pdgtype = 0;    /* don't bother with the PDG type here */
      ps->in[ps->mult].id = i+1;       /* unique value for this particle within the event */
      ps->in[ps->mult].parentid = 0;   /* All internally generated particles have no parent */
      ps->in[ps->mult].mech = 0;       /* maybe this should be set to something? */
      ps->in[ps->mult].momentum = make_s_Momentum();
      ps->in[ps->mult].momentum->px = kin.particle(i).px();
      ps->in[ps->mult].momentum->py = kin.particle(i).py();
      ps->in[ps->mult].momentum->pz = kin.particle(i).pz();
      ps->in[ps->mult].momentum->E  = kin.particle(i).e();
      
    }
  
  if(nParticles>0) flush_s_HDDM(thisOutputEvent, m_OutputStream);
  m_eventCounter++;
}
