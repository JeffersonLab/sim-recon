// -----------------------------------------
// DEventProcessor_truth_dirc.cc
// created on: 07.04.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#include "DEventProcessor_truth_dirc.h"

// Routine used to create our DEventProcessor
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_truth_dirc());
  }
}

DEventProcessor_truth_dirc::DEventProcessor_truth_dirc() {
}

DEventProcessor_truth_dirc::~DEventProcessor_truth_dirc() {
}

jerror_t DEventProcessor_truth_dirc::init(void) {
  string locOutputFileName = "hd_root.root";
  if(gPARMS->Exists("OUTPUT_FILENAME"))
    gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

  //go to file
  TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
  if(locFile != NULL)
    locFile->cd("");
  else
    gDirectory->Cd("/");

  hTruthBarHitXY = new TH2F("hTruthBarHitXY", "; Bar Hit X (cm); Bar Hit Y (cm)", 200, -100, 100, 200, -100, 100);
  hTruthBarHitBar = new TH1F("hTruthBarHitBar", "; Bar #", 48, 0.5, 47.5);
  hTruthPmtHitZY_North = new TH2F("hTruthPmtHitZY_North", "North Box; PMT Hit Z (cm); PMT Hit Y (cm)", 100, 525, 560, 110, 0., 110.);
  hTruthPmtHitZY_South = new TH2F("hTruthPmtHitZY_South", "South Box; PMT Hit Z (cm); PMT Hit Y (cm)", 100, 525, 560, 110, -110., 0.);
 
  hTruthPmtHit_North = new TH2F("hTruthPmtHit_North", "North Box; Pmt Hit Column ; Pixel Hit Row", 6, 0, 6, 18, 0, 18);
  hTruthPmtHit_South = new TH2F("hTruthPmtHit_South", "South Box; Pmt Hit Column ; Pixel Hit Row", 6, 0, 6, 18, 0, 18);
  hTruthPixelHit_North = new TH2F("hTruthPixelHit_North", "North Box; Pixel Hit X ; Pixel Hit Y", 144, 0, 144, 48, 0, 48);
  hTruthPixelHit_South = new TH2F("hTruthPixelHit_South", "South Box; Pixel Hit X ; Pixel Hit Y", 144, 0, 144, 48, 0, 48);
 
  return NOERROR;
}

jerror_t DEventProcessor_truth_dirc::brun(jana::JEventLoop *loop, int32_t runnumber)
{
   // Get the geometry
   DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
   DGeometry *geom = dapp->GetDGeometry(runnumber);

   // Outer detector geometry parameters
   vector<double>tof_face;
   geom->Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z", tof_face);
   vector<double>tof_plane;  
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", tof_plane);
   double dTOFz=tof_face[2]+tof_plane[2]; 
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", tof_plane);
   dTOFz+=tof_face[2]+tof_plane[2];
   dTOFz*=0.5;  // mid plane between tof Planes
   std::cout<<"dTOFz "<<dTOFz<<std::endl;

   double dDIRCz;
   vector<double>dirc_face;
   vector<double>dirc_plane;
   vector<double>dirc_shift;
   vector<double>bar_plane;
   geom->Get("//section/composition/posXYZ[@volume='DIRC']/@X_Y_Z", dirc_face);
   geom->Get("//composition[@name='DRCC']/mposY[@volume='DCML']/@Z_X/plane[@value='1']", dirc_plane);
   geom->Get("//composition[@name='DIRC']/posXYZ[@volume='DRCC']/@X_Y_Z", dirc_shift);
   geom->Get("//composition[@name='DCBR']/mposX[@volume='QZBL']/@Y_Z", bar_plane);
   
   dDIRCz=dirc_face[2]+dirc_plane[0]+dirc_shift[2]+bar_plane[1]; // 585.862
   std::cout<<"dDIRCz "<<dDIRCz<<std::endl;

   return NOERROR;
}

jerror_t DEventProcessor_truth_dirc::evnt(JEventLoop *loop, uint64_t eventnumber) {
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

  for (unsigned int j = 0; j < dircBarHits.size(); j++){
    double px = dircBarHits[j]->px;
    double py = dircBarHits[j]->py;
    double pz = dircBarHits[j]->pz;

    double x = dircBarHits[j]->x;
    double y = dircBarHits[j]->y;
    double z = dircBarHits[j]->z;
    int bar = dircBarHits[j]->bar;

    japp->RootWriteLock(); //ACQUIRE ROOT LOCK
    hTruthBarHitXY->Fill(x, y);
    hTruthBarHitBar->Fill(bar);
    japp->RootUnLock();
  }

  for (unsigned int h = 0; h < dircPmtHits.size(); h++){
	
     int ch=dircPmtHits[h]->ch;
     int pmt=ch/64;
     int pix=ch%64;
     double x = dircPmtHits[h]->x;
     double y = dircPmtHits[h]->y;
     double z = dircPmtHits[h]->z;
     double t = dircPmtHits[h]->t;

     // get PMT labels
     int pmt_column = pmt/18; // 0 - 5
     int pmt_row = 17 - pmt%18; // 0 - 17

     // get pixel labels
     int pixel_column = pix/8; // 0 - 7
     int pixel_row = 7 - pix%8; // 0 - 7
     int pixel_x = 8*pmt_row + pixel_row;
     int pixel_y = 47 - (8*pmt_column + pixel_column);

     // temporary geometry check
     //if(z > 550. && fabs(y) > 80) continue;
     //if(pmt!=0  && pmt!=106 && pmt!=107) continue;

     japp->RootWriteLock(); //ACQUIRE ROOT LOCK
     if(x < 0.) {
	hTruthPmtHitZY_South->Fill(z, y);
	hTruthPmtHit_South->Fill(pmt_column, pmt_row);
	hTruthPixelHit_South->Fill(pixel_x, pixel_y);
     }
     else {
	hTruthPmtHitZY_North->Fill(z, y);
	hTruthPmtHit_North->Fill(pmt_column, pmt_row);
	hTruthPixelHit_North->Fill(pixel_x, pixel_y);
     }
     japp->RootUnLock();
  }
  
  return NOERROR;
}

jerror_t DEventProcessor_truth_dirc::erun(void) {
  return NOERROR;
}

jerror_t DEventProcessor_truth_dirc::fini(void) {
  return NOERROR;
}
