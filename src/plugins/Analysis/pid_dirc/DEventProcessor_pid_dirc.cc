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
  string locOutputFileName = "hd_root.root";
  if(gPARMS->Exists("OUTPUT_FILENAME"))
    gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

  //go to file
  TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
  if(locFile != NULL)
    locFile->cd("");
  else
    gDirectory->Cd("/");

	
  
  fTree = new TTree("dirc", "dirc tree");
  fcEvent = new TClonesArray("DrcEvent");
  // fEvent = new DrcEvent();
  //fTree->Branch("DrcEvent","DrcEvent",&fEvent,256000,0);
  fTree->Branch("DrcEvent",&fcEvent,256000,0);
  return NOERROR;
}

jerror_t DEventProcessor_pid_dirc::evnt(JEventLoop *loop, uint64_t eventnumber) {
  vector<const DBeamPhoton*> beam_photons;
  vector<const DMCThrown*> mcthrowns;
  vector<const DMCTrackHit*> mctrackhits;
  vector<const DDIRCTruthBarHit*> dircBarHits;
  vector<const DDIRCTruthPmtHit*> dircPmtHits;
  
  loop->Get(beam_photons);
  loop->Get(mcthrowns);
  loop->Get(mctrackhits);
  loop->Get(dircPmtHits);
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
    case SYS_CHERENKOV:
      cerepoints[mctrackhit->track]++;
      break;
    default:
      break;
    }
  }
  
  //  std::cout<<"#hits "<< dircPmtHits.size()<<std::endl;

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK
  
  // loop over mc/reco tracks
  TClonesArray& cevt = *fcEvent;
  cevt.Clear();
  for (unsigned int m = 0; m < mcthrowns.size(); m++){
    if(dircPmtHits.size() > 0.){
      fEvent = new DrcEvent();
      DrcHit hit;
      // loop over PMT's hits
      for (unsigned int h = 0; h < dircPmtHits.size(); h++){
	int relevant(0);
	// identify bar id
	for (unsigned int j = 0; j < dircBarHits.size(); j++){
	  if(j != fabs(dircPmtHits[h]->key_bar)) continue;
	  if(mcthrowns[m]->myid == dircBarHits[j]->track){
	    // double px = mcthrowns[m]->momentum().X();
	    // double py = mcthrowns[m]->momentum().Y();
	    // double pz = mcthrowns[m]->momentum().Z();
	    
	    double px = dircBarHits[j]->px;
	    double py = dircBarHits[j]->py;
	    double pz = dircBarHits[j]->pz;

	    fEvent->SetMomentum(TVector3(px,py,pz));
	    fEvent->SetPdg(mcthrowns[m]->pdgtype);
	    fEvent->SetTime(dircBarHits[j]->t);
	    fEvent->SetParent(mcthrowns[m]->parentid);
	    fEvent->SetId(dircBarHits[j]->bar);// bar id where the particle hit the detector
	    fEvent->SetPosition(TVector3(dircBarHits[j]->x, dircBarHits[j]->y, dircBarHits[j]->z)); // position where the charged particle hit the radiator
	    relevant++;
	  }
	}
	
	if(relevant<1) continue;
	int ch=dircPmtHits[h]->ch;
	int pmt=ch/64;
	int pix=ch%64;
	
	hit.SetChannel(dircPmtHits[h]->ch);
	hit.SetPmtId(pmt);
	hit.SetPixelId(pix);
	hit.SetPosition(TVector3(dircPmtHits[h]->x,dircPmtHits[h]->y,dircPmtHits[h]->z));
	hit.SetEnergy(dircPmtHits[h]->key_bar);
	hit.SetLeadTime(dircPmtHits[h]->t);
	hit.SetPathId(dircPmtHits[h]->E);
	fEvent->AddHit(hit);
      }
      
      if(fEvent->GetHitSize()>0) new (cevt[ cevt.GetEntriesFast()]) DrcEvent(*fEvent);
      
    }
  }
  
  if(cevt.GetEntriesFast()>0) fTree->Fill();
  japp->RootUnLock(); //RELEASE ROOT LOCK
      
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
  return NOERROR;
}
