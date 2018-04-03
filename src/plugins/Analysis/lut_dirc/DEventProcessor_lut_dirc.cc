// -----------------------------------------
// DEventProcessor_lut_dirc.cc
// created on: 29.11.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#include "DEventProcessor_lut_dirc.h"

// Routine used to create our DEventProcessor
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_lut_dirc());
  }
}

DEventProcessor_lut_dirc::DEventProcessor_lut_dirc() {
  fTree = NULL;
}

DEventProcessor_lut_dirc::~DEventProcessor_lut_dirc() {
}

jerror_t DEventProcessor_lut_dirc::init(void) {
  string locOutputFileName = "hd_root.root";
  if(gPARMS->Exists("OUTPUT_FILENAME"))
    gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

  //go to file
  TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
  if(locFile != NULL)
    locFile->cd("");
  else
    gDirectory->Cd("/");

  fTree = new TTree("lut_dirc","Look-up table for the geometrical reconstruction.");
  int Nnodes = 30000;    
  for(int l=0; l<48; l++){
    fLut[l] = new TClonesArray("DrcLutNode");
    fTree->Branch(Form("LUT_%d",l),&fLut[l],256000,-1);
    TClonesArray &fLuta = *fLut[l];
    for (Long64_t n=0; n<Nnodes; n++) {
      new((fLuta)[n]) DrcLutNode(-1);
    }
  }

  
  return NOERROR;
}

jerror_t DEventProcessor_lut_dirc::evnt(JEventLoop *loop, uint64_t eventnumber) {
  vector<const DMCThrown*> mcthrowns;
  vector<const DMCTrackHit*> mctrackhits;
  vector<const DDIRCTruthBarHit*> dircBarHits;
  vector<const DDIRCTruthPmtHit*> dircPmtHits;
  
  loop->Get(mcthrowns);
  loop->Get(mctrackhits);
  loop->Get(dircPmtHits);
  loop->Get(dircBarHits);

  if(mcthrowns.size()<1) return NOERROR;
  //if(dircBarHits.size()<1) return NOERROR;
  if(dircPmtHits.size()!=1) return NOERROR;
  
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK

  TVector3 mom(0,0,0);
  if(dircBarHits.size()>0){
    mom =  TVector3(dircBarHits[0]->px,
		    dircBarHits[0]->py,
		    dircBarHits[0]->pz).Unit();
  }
  
  // loop over PMT's hits
  for (unsigned int h = 0; h < dircPmtHits.size(); h++){
    
    int ch=dircPmtHits[h]->ch;
    int pmt=ch/64;
    int pix=ch%64;
    int id = 100*pmt + pix;
    int lutId = dircPmtHits[h]->key_bar;
    TVector3 dir =  TVector3(mcthrowns[0]->momentum().X(),
			     mcthrowns[0]->momentum().Y(),
			     mcthrowns[0]->momentum().Z()).Unit();

    std::cout<<"dir.X() "<<dir.X() <<" "<<dir.Y() <<" "<<dir.Z() << " | "
	     <<mom.X() <<" "<<mom.Y() <<" "<<mom.Z() <<std::endl;
    
    
    if(lutId>=0 && lutId<48)
      ((DrcLutNode*)(fLut[lutId]->At(id)))->
	AddEntry(lutId,               // lut/bar id
		 id,                  // pixel id
		 dir,
		 dircPmtHits[h]->path,
		 dircPmtHits[h]->refl,
		 dircPmtHits[h]->t,
		 TVector3(dircPmtHits[h]->x,dircPmtHits[h]->y,dircPmtHits[h]->z),
		 TVector3(dircPmtHits[h]->x,dircPmtHits[h]->y,dircPmtHits[h]->z));
  }
  japp->RootUnLock(); //RELEASE ROOT LOCK

  return NOERROR;
}

jerror_t DEventProcessor_lut_dirc::erun(void) {
  return NOERROR;
}

jerror_t DEventProcessor_lut_dirc::fini(void) {
  japp->RootWriteLock();
  fTree->Fill();
  japp->RootUnLock();
  return NOERROR;
}
