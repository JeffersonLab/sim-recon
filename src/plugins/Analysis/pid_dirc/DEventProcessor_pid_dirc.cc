// -----------------------------------------
// DEventProcessor_pid_dirc.cc
// created on: 07.04.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#include "DEventProcessor_pid_dirc.h"

// Routine used to create our DEventProcessor
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_pid_dirc());
  }
}

DEventProcessor_pid_dirc::DEventProcessor_pid_dirc() {
  fTree = NULL;
  fEvent = NULL;
}

DEventProcessor_pid_dirc::~DEventProcessor_pid_dirc() {
}

jerror_t DEventProcessor_pid_dirc::init(void) {
  fRootFile = new TFile("drc.root","RECREATE");
  fTree = new TTree("dirc", "dirc tree");
  // fEvent = new TClonesArray("DrcEvent");
  fEvent = new DrcEvent();
  fTree->Branch("DrcEvent","DrcEvent",&fEvent,256000,0);
  return NOERROR;
}

jerror_t DEventProcessor_pid_dirc::evnt(JEventLoop *loop, uint64_t eventnumber) {
  vector<const DBeamPhoton*> beam_photons;
  vector<const DMCThrown*> mcthrowns;
  vector<const DMCTrackHit*> mctrackhits;
  vector<const DDIRCTruthBarHit*> dircBarHits;
  vector<const DDIRCTruthMcpHit*> dircMcpHits;

  loop->Get(beam_photons);
  loop->Get(mcthrowns);
  loop->Get(mctrackhits);
  loop->Get(dircMcpHits);
  loop->Get(dircBarHits);

  TVector3 VertexGen = TVector3(mcthrowns[0]->position().X(),
				mcthrowns[0]->position().Y(), mcthrowns[0]->position().Z());
  // Make Particle object for beam photon
  TLorentzVector beam_photon(0.0, 0.0, 9.0, 9.0);
  if (beam_photons.size() > 0) {
    const DLorentzVector &lv = beam_photons[0]->lorentzMomentum();
    beam_photon.SetPxPyPzE(lv.Px(), lv.Py(), lv.Pz(), lv.E());
  }

  // Target is proton at rest in lab frame
  TLorentzVector target(0.0, 0.0, 0.0, 0.93827);

  particle_set thr;

  // Loop over track hits
  // map first: DMCTrack index+1
  // map second: number of hits
  map<int, int> cdchits;
  map<int, int> fdchits;
  map<int, int> bcalhits;
  map<int, int> fcalhits;
  map<int, int> upvhits;
  map<int, int> tofhits;
  map<int, int> richpoints;
  map<int, int> cerepoints;

  for (unsigned int i = 0; i < mctrackhits.size(); i++) {
    const DMCTrackHit *mctrackhit = mctrackhits[i];

    switch (mctrackhit->system) {
    case SYS_CDC:
      if (mctrackhit->primary)
	cdchits[mctrackhit->track]++;
      break;
    case SYS_FDC:
      if (mctrackhit->primary)
	fdchits[mctrackhit->track]++;
      break;
    case SYS_BCAL:
      bcalhits[mctrackhit->track]++;
      break;
    case SYS_FCAL:
      fcalhits[mctrackhit->track]++;
      break;
    case SYS_UPV:
      upvhits[mctrackhit->track]++;
      break;
    case SYS_TOF:
      tofhits[mctrackhit->track]++;
      break;
      // case SYS_RICH:
      // 	richpoints[mctrackhit->track]++;
      // 	break;
    case SYS_CHERENKOV:
      cerepoints[mctrackhit->track]++;
      break;
    default:
      break;
    }
  }
  
  std::cout<<"#hits "<< dircMcpHits.size()<<std::endl;

  // loop over mc/reco tracks
  for (unsigned int j = 0; j < mcthrowns.size(); j++){

    fEvent = new DrcEvent();
    DrcHit hit;
    // loop over MCP's hits
    for (unsigned int h = 0; h < dircMcpHits.size(); h++){
      // identify bar id
      for (unsigned int j = 0; j < dircBarHits.size(); j++){
      }

      int ch=dircMcpHits[h]->ch;
      int mcp=ch/64;
      int pix=ch%64;
      
      hit.SetChannel(dircMcpHits[h]->ch);
      hit.SetMcpId(mcp);
      hit.SetPixelId(pix);
      hit.SetPosition(TVector3(dircMcpHits[h]->x,dircMcpHits[h]->y,dircMcpHits[h]->z));
      hit.SetMomentum(TVector3(0,0,0));
      hit.SetLeadTime(dircMcpHits[h]->t);
      fEvent->AddHit(hit);
      //      fHit = new DrcHit(hit);
    }
    
    //japp->RootWriteLock(); //ACQUIRE ROOT LOCK
    fTree->Fill();
    japp->RootUnLock(); //RELEASE ROOT LOCK
    //fEvent->Clear();

    // TClonesArray& cevt = *fEvent;
    // Int_t size = cevt.GetEntriesFast();
    // new (cevt[size]) DrcEvent(evt);
  }
    
  // // Although we are only filling objects local to this plugin, TTree::Fill() periodically writes to file: Global ROOT lock
  // japp->RootWriteLock(); //ACQUIRE ROOT LOCK

  // // Copy event number
  // //fTree->Fill();
  // japp->RootUnLock(); //RELEASE ROOT LOCK

  return NOERROR;
}

Particle DEventProcessor_pid_dirc::MakeParticle(const DKinematicData *kd,
						double mass, hit_set hits) {
  // Create a ROOT TLorentzVector object out of a Hall-D DKinematic Data object.
  // Here, we have the mass passed in explicitly rather than use the mass contained in
  // the DKinematicData object itself. This is because right now (Feb. 2009) the
  // PID code is not mature enough to give reasonable guesses. See above code.

  double p = kd->momentum().Mag();
  double theta = kd->momentum().Theta();
  double phi = kd->momentum().Phi();
  double px = p * sin(theta) * cos(phi);
  double py = p * sin(theta) * sin(phi);
  double pz = p * cos(theta);
  double E = sqrt(mass * mass + p * p);
  double x = kd->position().X();
  double y = kd->position().Y();
  double z = kd->position().Z();

  Particle part;
  part.p.SetPxPyPzE(px, py, pz, E);
  part.x.SetXYZ(x, y, z);
  part.P = p;
  part.E = E;
  part.Th = theta;
  part.Ph = phi;
  part.hits_cdc = hits.hits_cdc;
  part.hits_fdc = hits.hits_fdc;
  part.hits_bcal = hits.hits_bcal;
  part.hits_fcal = hits.hits_fcal;
  part.hits_upv = hits.hits_upv;
  part.hits_tof = hits.hits_tof;
  part.hits_rich = hits.hits_rich;
  part.hits_cere = hits.hits_cere;

  return part;
}


jerror_t DEventProcessor_pid_dirc::erun(void) {
  return NOERROR;
}

jerror_t DEventProcessor_pid_dirc::fini(void) {
  //japp->RootWriteLock(); //ACQUIRE ROOT LOCK
  if(fRootFile) fRootFile->Write();
  //japp->RootUnLock(); //RELEASE ROOT LOCK
  return NOERROR;
}
