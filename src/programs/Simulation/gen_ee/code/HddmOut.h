/*
 * HddmOut.h
 *
 *  Created on: Nov 14, 2013
 *      Author: ben
 */

#ifndef HDDMOUT_H_
#define HDDMOUT_H_

using namespace std;

#include "HDDM/hddm_s.h"

struct tmpEvt_t {
  int nGen;
  int rxn;
  double weight;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector q1;
  TLorentzVector q2;
  TLorentzVector recoil;
};

class HddmOut {
 private:
  s_iostream_t* ostream;
  //TDatabasePDG* pdg;
  s_PhysicsEvents_t* phyEvt;
  s_Reactions_t* reactions;
  s_Reaction_t* reaction;
  s_Target_t* target;
  s_Beam_t* beam;
  s_Vertices_t* vertices;
  s_HDDM_t* hddmEvt;
  s_Origin_t* origin;
  s_Products_t* products;

  Particle_t targetType;
  Particle_t beamType;
  
 public:
  HddmOut(string filename) {
      cout << "opening HDDM file: " << filename << endl;
    ostream = init_s_HDDM((char*)filename.c_str());
    targetType = Proton;
    beamType = Gamma;
  }
  
  ~HddmOut() {
    close_s_HDDM(ostream);
  }
  
  void init(int runNo) {
    //This sets the run number and event characteristics
    //The HDDM entry has one event, which has one reaction
    hddmEvt = make_s_HDDM();
    hddmEvt->physicsEvents = phyEvt = make_s_PhysicsEvents(1);
    phyEvt->mult = 1;
    phyEvt->in[0].runNo = runNo;
    
    //We define beam and target parameters for the reaction, which
    //remain the same between events
    phyEvt->in[0].reactions = reactions = make_s_Reactions(1);
    reactions->mult = 1;
    reaction = &reactions->in[0];
    reaction->target = target = make_s_Target();
    target->type = targetType;
    target->properties = make_s_Properties();
    target->properties->charge = ParticleCharge(targetType);
    target->properties->mass = ParticleMass(targetType);
    target->momentum = make_s_Momentum();
    target->momentum->px = 0;
    target->momentum->py = 0;
    target->momentum->pz = 0;
    target->momentum->E  = ParticleMass(targetType);
    reaction->beam = beam = make_s_Beam();
    beam->type = beamType;
    beam->properties = make_s_Properties();
    beam->properties->charge = ParticleCharge(beamType);
    beam->properties->mass = ParticleMass(beamType);
    beam->momentum = make_s_Momentum();

  }
  
  void write(tmpEvt_t evt, int eventNum) {
    init(10000);
    phyEvt->in[0].eventNo = eventNum;
    reaction->vertices = vertices = make_s_Vertices(1);
    vertices->mult = 1;
    vertices->in[0].origin = origin = make_s_Origin();
    vertices->in[0].products = products = make_s_Products(evt.nGen);
    
    origin->t = 0.0;
    origin->vx = 0.0;
    origin->vy = 0.0;
    origin->vz = 0.0;

    beam->momentum->px = evt.beam.Px();
    beam->momentum->py = evt.beam.Py();
    beam->momentum->pz = evt.beam.Pz();
    beam->momentum->E  = evt.beam.E();

    products->mult = evt.nGen;
    reaction->weight = evt.weight;

    //PRODUCED ELECTRON
    products->in[0].type = Electron;
    products->in[0].pdgtype = 11;
    products->in[0].id = 1;
    products->in[0].parentid = 0;
    products->in[0].mech = 0;
    products->in[0].momentum = make_s_Momentum();
    products->in[0].momentum->px = evt.q1.Px();
    products->in[0].momentum->py = evt.q1.Py();
    products->in[0].momentum->pz = evt.q1.Pz();
    products->in[0].momentum->E = evt.q1.E();

    //PRODUCED ELECTRON
    products->in[1].type = Positron;
    products->in[1].pdgtype = -11;
    products->in[1].id = 2;
    products->in[1].parentid = 0;
    products->in[1].mech = 0;
    products->in[1].momentum = make_s_Momentum();
    products->in[1].momentum->px = evt.q2.Px();
    products->in[1].momentum->py = evt.q2.Py();
    products->in[1].momentum->pz = evt.q2.Pz();
    products->in[1].momentum->E = evt.q2.E();

    //RECOIL
    if (evt.rxn == 2) {//set type
      products->in[2].type = Proton;
      products->in[2].pdgtype = 2212;
    } else {
      products->in[2].type = Electron;
      products->in[2].pdgtype = 11;
    }
    products->in[2].id = 3;
    products->in[2].parentid = 0;
    products->in[2].mech = 0;
    products->in[2].momentum = make_s_Momentum();
    products->in[2].momentum->px = evt.recoil.Px();
    products->in[2].momentum->py = evt.recoil.Py();
    products->in[2].momentum->pz = evt.recoil.Pz();
    products->in[2].momentum->E = evt.recoil.E();

    flush_s_HDDM(hddmEvt, ostream);

  }
};

#endif /* HDDMOUT_H_ */
