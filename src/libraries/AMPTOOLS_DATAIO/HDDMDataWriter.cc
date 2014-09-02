
#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"
#include "HDDM/hddm_s.hpp"

HDDMDataWriter::HDDMDataWriter(const string& outFile, int runNumber)
{
  m_OutputFile = new ofstream(outFile.c_str());
  m_OutputStream = new hddm_s::ostream(*m_OutputFile);
  m_runNumber = runNumber;
  
  m_eventCounter = 0;
}

HDDMDataWriter::~HDDMDataWriter()
{
  delete m_OutputStream;
  delete m_OutputFile;
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype, bool centeredVertex)
{
  if (centeredVertex)
    writeEvent(kin,ptype,0,0,65/*cm*/);
  else
    writeEvent(kin,ptype,0,0,50/*cm*/,80/*cm*/);
}

void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype,
	    float vx, float vy, float vz_min, float vz_max)
{
  if (vz_min > vz_max) {
    float tmp=vz_min;
    vz_min=vz_max;
    vz_max=tmp;
  }
  writeEvent(kin, ptype,vx, vy, (vz_max - vz_min) * drand48() + vz_min);
}


void HDDMDataWriter::
writeEvent( const Kinematics& kin, const vector<int>& ptype,
	    float vx, float vy, float vz)
{
  vector< HepLorentzVector > particleList = kin.particleList();
  int nParticles=kin.particleList().size();
  
  // Start a new event in the HDDM record
  hddm_s::HDDM record;
  hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
  pes().setRunNo(m_runNumber);
  pes().setEventNo(m_eventCounter);
  hddm_s::ReactionList rs = pes().addReactions();
  hddm_s::VertexList vs = rs().addVertices();
  hddm_s::OriginList os = vs().addOrigins();
  hddm_s::ProductList ps = vs().addProducts(nParticles-1);

  os().setT(0.0);
  os().setVx(vx);
  os().setVy(vy);
  os().setVz(vz);

  hddm_s::BeamList bs = rs().addBeams();
  bs().setType((Particle_t)1);
  hddm_s::MomentumList bmoms = bs().addMomenta();
  bmoms().setPx(kin.particle(0).px());
  bmoms().setPy(kin.particle(0).py());
  bmoms().setPz(kin.particle(0).pz());
  bmoms().setE(kin.particle(0).e());
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
  
  for(int i=1; i < nParticles; i++)
  {
      ps(i-1).setType((Particle_t)ptype[i]);
      ps(i-1).setPdgtype(0);    /* don't bother with the PDG type here */
      ps(i-1).setId(i);         /* unique value for this particle within the event */
      ps(i-1).setParentid(0);   /* All internally generated particles have no parent */
      ps(i-1).setMech(0);       /* maybe this should be set to something? */
      hddm_s::MomentumList pmoms = ps(i-1).addMomenta();
      pmoms().setPx(kin.particle(i).px());
      pmoms().setPy(kin.particle(i).py());
      pmoms().setPz(kin.particle(i).pz());
      pmoms().setE(kin.particle(i).e());
  }
  
  if (nParticles > 0)
    *m_OutputStream << record;
  m_eventCounter++;
}
