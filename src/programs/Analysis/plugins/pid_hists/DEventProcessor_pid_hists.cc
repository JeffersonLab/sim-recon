#include "DEventProcessor_pid_hists.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
using namespace std;

#include <TThread.h>
#include <TDirectoryFile.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TMath.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <PID/DKinematicData.h>
#include <PID/DChargedTrack.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DTrackFitter.h>

//#define TOF_SIGMA 0.080   // TOF resolution in ns

// Routine used to create our DEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_pid_hists());
  }
} // "C"


DEventProcessor_pid_hists::DEventProcessor_pid_hists()
{
  pthread_mutex_init(&mutex, NULL);
}

DEventProcessor_pid_hists::~DEventProcessor_pid_hists()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_pid_hists::init(void)
{
  // Create PID directory
  TDirectory *dir = (TDirectory*)gROOT->FindObject("PID");
  if(!dir)dir = new TDirectoryFile("PID","PID");
  dir->cd();
  pid = (TH1F*)gROOT->FindObject("pid");
  if(!pid)pid = new TH1F("pid","Charged PID Counts", 10,0,10);
  kine = (TH2F*)gROOT->FindObject("kine");
  if(!kine)kine = new TH2F("kine","Kinematics of Misidentified Tracks",100,0,180,100,0,9);
  kineBadRec = (TH2F*)gROOT->FindObject("kineBadRec");
  if(!kineBadRec)kineBadRec = new TH2F("kineBadRec","Kinematics of Poorly Reconstructed Tracks",100,0,180,100,0,9);

  tof = (TH1F*)gROOT->FindObject("tof");
  if(!tof)tof = new TH1F("tof","TOF", 500,0,30);
  deltaBeta = (TH1F*)gROOT->FindObject("deltaBeta");
  if(!deltaBeta)deltaBeta = new TH1F("deltaBeta","#Delta #beta", 100,-1,1);
  dEdx = (TH1F*)gROOT->FindObject("dEdx");
  if(!dEdx)dEdx = new TH1F("dEdx","dE/dx", 100,0,.02);

  tdiff_Kpi = (TH1F*)gROOT->FindObject("tdiff_Kpi");
  if(!tdiff_Kpi)tdiff_Kpi = new TH1F("tdiff_Kpi","tdiff_Kpi", 1000,0,10000);
  tdiff_Kpi_mom = (TH2F*)gROOT->FindObject("tdiff_Kpi_mom");
  if(!tdiff_Kpi_mom)tdiff_Kpi_mom = new TH2F("tdiff_Kpi_mom","tdiff_Kpi_mom", 1000,0,9,1000,0,10000);
  tdiff_ppi = (TH1F*)gROOT->FindObject("tdiff_ppi");
  if(!tdiff_ppi)tdiff_ppi = new TH1F("tdiff_ppi","tdiff_ppi", 1000,0,10000);
  tdiff_ppi_mom = (TH2F*)gROOT->FindObject("tdiff_ppi_mom");
  if(!tdiff_ppi_mom)tdiff_ppi_mom = new TH2F("tdiff_ppi_mom","tdiff_ppi_mom", 1000,0,9,1000,0,10000);
  tdiff_pK = (TH1F*)gROOT->FindObject("tdiff_pK");
  if(!tdiff_pK)tdiff_pK = new TH1F("tdiff_pK","tdiff_pK", 1000,0,10000);
  tdiff_pK_mom = (TH2F*)gROOT->FindObject("tdiff_pK_mom");
  if(!tdiff_pK_mom)tdiff_pK_mom = new TH2F("tdiff_pK_mom","tdiff_pK_mom", 1000,0,9,1000,0,10000);

  dir->cd("../");

  NgoodPid=0;
  NbadRecon=0;
  NpipTruth_protonPid=0;
  NpipTruth_KpPid=0;
  NprotonTruth_pipPid=0;
  NprotonTruth_KpPid=0;
  NKpTruth_pipPid=0;
  NKpTruth_protonPid=0;
  NpimTruth_KmPid=0;
  NKmTruth_pimPid=0;

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_pid_hists::brun(JEventLoop *loop, int runnumber)
{

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_pid_hists::evnt(JEventLoop *loop, int eventnumber)
{
  // Get reconstructed and thrown objects
  vector<const DMCThrown*> mcthrowns;
  vector<const DChargedTrack*> chargedtracks;
  vector<const DTrackTimeBased*> tracks;
  loop->Get(mcthrowns);
  loop->Get(chargedtracks);
  loop->Get(tracks);
 
  double cut = 1e-3; // pid Chi^2 probability cut
  bool DEBUG = false;

  // get tracks found through truth matching
  vector<const DTrackTimeBased*> truthTracks;
  for(unsigned int j=0; j<tracks.size(); j++){
    const DTrackTimeBased *track = tracks[j];
    double FOM = GetTruthMatchingFOM(j,track,mcthrowns);
    //double FOM = DTrackTimeBased_factory::GetTruthMatchingFOM(j,track,mcthrowns);
    if (FOM>cut) truthTracks.push_back(track);
  }

  // debug
  if (DEBUG) {
    for(unsigned int j=0; j<truthTracks.size(); j++){
      const DTrackTimeBased *track = truthTracks[j];
      cout << "(pidFOM,trkChi2,trkNdof):  " << "(" << track->FOM << "," << track->chisq << "," << track->Ndof << ")" << endl;
      printInfo(track,eventnumber,j,truthTracks.size(),"truthTrack");
    }
  }

  // calculate time differences between tracks found using normal PID (best hypothesis) and the other hypotheses 
  for(unsigned int j=0; j<chargedtracks.size(); j++){
    if(chargedtracks[j]->hypotheses.size()==0)continue;
    const DTrackTimeBased *track = chargedtracks[j]->hypotheses[0];
    if (DEBUG) printInfo(track,eventnumber,j,truthTracks.size(),"track");
    double T1 = 1000.*track->pathLength()/track->lorentzMomentum().Beta()/SPEED_OF_LIGHT; // in ps
    int type1 = GetType(track);
    double mom = track->momentum().Mag(); // momentum of reconstructed best hypothesis
    for (unsigned int k=1;k<chargedtracks[j]->hypotheses.size();k++) {
      const DTrackTimeBased *track_can = chargedtracks[j]->hypotheses[k];
      if (DEBUG) printInfo(track_can,eventnumber,k,truthTracks.size(),"track_can");
      double T2 = 1000.*track_can->pathLength()/track_can->lorentzMomentum().Beta()/SPEED_OF_LIGHT; // in ps
      int type2 = GetType(track_can);
      double dT = T1 - T2;
      if (DEBUG) cout << "dT " << dT << endl;
      if ((type1==11&&type2==8)||(type1==12&&type2==9)) {tdiff_Kpi->Fill(dT);tdiff_Kpi_mom->Fill(mom,dT);}
      if ((type1==8&&type2==11)||(type1==9&&type2==12)) {tdiff_Kpi->Fill(-dT);tdiff_Kpi_mom->Fill(mom,-dT);}
      if (type1==14&&type2==8) {tdiff_ppi->Fill(dT);tdiff_ppi_mom->Fill(mom,dT);}
      if (type1==8&&type2==14) {tdiff_ppi->Fill(-dT);tdiff_ppi_mom->Fill(mom,-dT);}
      if (type1==14&&type2==11) {tdiff_pK->Fill(dT);tdiff_pK_mom->Fill(mom,dT);}
      if (type1==11&&type2==14) {tdiff_pK->Fill(-dT);tdiff_pK_mom->Fill(mom,-dT);}
    }
  }

  // loop over tracks found using normal PID and use the truth information to count the misidentified tracks
  for(unsigned int j=0; j<chargedtracks.size(); j++){
    if(chargedtracks[j]->hypotheses.size()==0)continue;
    const DTrackTimeBased *track = chargedtracks[j]->hypotheses[0];

    // skip tracks with a low FOM
    if (track->FOM<=cut)continue; 

    int type = GetType(track);
  		
    tof->Fill(track->TOF());
    deltaBeta->Fill(track->deltaBeta());
    dEdx->Fill(track->dEdx());

    bool match = false;
    // check if the track has the correct PID (according to truth matched track)
    for(unsigned int k=0; k<truthTracks.size(); k++){
      const DTrackTimeBased *truthTrack = truthTracks[k];
      if (track==truthTrack) {match=true;}
    }
    
    if (match) {
      if (DEBUG) cout << "match" << endl; 
      NgoodPid++; pid->Fill(0);
    }
    else {
      int MAX_TRACKS = (int)mcthrowns.size()+1, thrownIndex=-1; double f = 0.;
      //DTrackTimeBased_factory::GetThrownIndex(track,MAX_TRACKS,f,thrownIndex);
      GetThrownIndex(track,MAX_TRACKS,f,thrownIndex);
      //GetThrownIndex2(track,mcthrowns,thrownIndex);
      double chisq = (thrownIndex<=0 || thrownIndex>=MAX_TRACKS) ? 1e10:GetMomChisq(track,mcthrowns[thrownIndex-1]);
      double FOM = TMath::Prob(chisq,3);
      if (f>0.5&&FOM>cut) { 
	int typeTruth = mcthrowns[thrownIndex-1]->type;
	if (DEBUG) {
	  cout << chisq << endl;
	  cout << "wrong pid" << endl;
	  cout << type << "    " << typeTruth << endl;
	}
	kine->Fill(mcthrowns[thrownIndex-1]->momentum().Theta()*TMath::RadToDeg(),mcthrowns[thrownIndex-1]->momentum().Mag());
	if (typeTruth==8&&type==14) {NpipTruth_protonPid++;  pid->Fill(2);} // pip misidentified as proton
	if (typeTruth==8&&type==11) {NpipTruth_KpPid++;  pid->Fill(3);} // pip misidentified as Kp
	if (typeTruth==14&&type==8) {NprotonTruth_pipPid++;  pid->Fill(4);} // proton misidentified as pip
	if (typeTruth==14&&type==11) {NprotonTruth_KpPid++;  pid->Fill(5);} // proton misidentified as Kp
	if (typeTruth==11&&type==8) {NKpTruth_pipPid++;  pid->Fill(6);} // Kp misidentified as pip
	if (typeTruth==11&&type==14) {NKpTruth_protonPid++;  pid->Fill(7);} // Kp misidentified as proton
	if (typeTruth==9&&type==12) {NpimTruth_KmPid++;  pid->Fill(8);} // pim misidentified as Km
	if (typeTruth==12&&type==9) {NKmTruth_pimPid++;  pid->Fill(9);} // Km misidentified as pim
	//if (mcthrowns[thrownIndex-1]->charge()*track->charge()<0) { cout << "wrong charge" << endl;}
      }
      else {
	if (DEBUG) cout << "poorly reconstructed" << endl;
	if (chisq!=1e10) kineBadRec->Fill(mcthrowns[thrownIndex-1]->momentum().Theta()*TMath::RadToDeg(),mcthrowns[thrownIndex-1]->momentum().Mag());
	NbadRecon++; pid->Fill(1);
      }
    }

    // debug
    if (DEBUG) {
      cout << "(pidFOM,trkChi2,trkNdof):  " << "(" << track->FOM << "," << track->chisq << "," << track->Ndof << ")" << endl;	
      printInfo(track,eventnumber,j,chargedtracks.size(),"normalTrack");
    }

  } // particles

  // debug
  // throwns
  if (DEBUG) {
    for(unsigned int k=0; k<mcthrowns.size(); k++){
      printInfo(mcthrowns[k],eventnumber,k,mcthrowns.size(),"thrown");
    }
  }

  return NOERROR;
}

// Returns a FOM based on difference between thrown and reconstructed momentum if track matches MC truth information, 
// returns a FOM=0 otherwise;
// a match requires identical masses and charges, and that more than half of the track's hits match the truth hits 
double DEventProcessor_pid_hists::GetTruthMatchingFOM(int trackIndex,const DTrackTimeBased *track,vector<const DMCThrown*>mcthrowns)  {
  bool match=false;
  
  DLorentzVector fourMom = track->lorentzMomentum(); 
  DLorentzVector gen_fourMom[mcthrowns.size()];
  for(unsigned int i=0; i<mcthrowns.size(); i++){
    gen_fourMom[i] = mcthrowns[i]->lorentzMomentum();
  }
  
  // Get info for thrown track
  int MAX_TRACKS = (int)mcthrowns.size()+1, thrownIndex=-1; double f = 0.;
  GetThrownIndex(track,MAX_TRACKS,f,thrownIndex);
  //int thrownIndex2=-1;
  //GetThrownIndex2(track,mcthrowns,thrownIndex2);
  if(thrownIndex<=0 || thrownIndex>=MAX_TRACKS || f<=0.5) return 0.;

  double chisq = GetMomChisq(track,mcthrowns[thrownIndex-1]);

  if (fabs(track->mass()-mcthrowns[thrownIndex-1]->mass())<0.01 && track->charge()==mcthrowns[thrownIndex-1]->charge()) 
    match = true;
  
  /*if (match) {
    double trk_chi2=track->chisq;
    unsigned int ndof=track->Ndof;
    cout << "f: " << f << endl;
    cout << "trk_chi2: " << trk_chi2 << endl;
    cout << "ndof: " << ndof << endl;
    cout << "throwncharge: " << mcthrowns[thrownIndex-1]->charge() << endl;
    cout << "trackcharge: " << track->charge() << endl;
    cout << "chargediff: " << fabs(track->charge()-mcthrowns[thrownIndex-1]->charge()) << endl;
    cout << "thrownmass: " << mcthrowns[thrownIndex-1]->mass() << endl;
    cout << "trackmass: " << track->mass() << endl;
    cout << "massdiff: " << fabs(track->mass()-mcthrowns[thrownIndex-1]->mass()) << endl;
    cout << "chisq: " << chisq << endl;
    cout << "match?: " << match << endl;
    cout << "thrownIndex2: " << thrownIndex2 << "   trackIndex: " << trackIndex << endl;
    cout << "thrownIndex: " << thrownIndex << "   trackIndex: " << trackIndex << endl;
    cout<< "track   " << setprecision(4) << "Px: " << fourMom.Px() << "    Py: " << fourMom.Py() << "   Pz: " << fourMom.Pz() << "   E: " << fourMom.E() << "    M: " << fourMom.M() << "   pt: " << fourMom.Pt() << "   theta: " << fourMom.Theta() << "   phi: " << fourMom.Phi() << endl; 
    cout<< "thrown  " << setprecision(4) << "Px: " << gen_fourMom[thrownIndex-1].Px() << "    Py: " << gen_fourMom[thrownIndex-1].Py() << "   Pz: " << gen_fourMom[thrownIndex-1].Pz() << "   E: " << gen_fourMom[thrownIndex-1].E() << "    M: " << gen_fourMom[thrownIndex-1].M() << "   pt: " << gen_fourMom[thrownIndex-1].Pt() << "   theta: " << gen_fourMom[thrownIndex-1].Theta() << "   phi: " << gen_fourMom[thrownIndex-1].Phi() << endl;
    }*/

  return (match) ?  TMath::Prob(chisq,3) : 0.0; 
}

//------------------
// GetThrownIndex
//------------------
void DEventProcessor_pid_hists::GetThrownIndex(const DKinematicData *kd, int &MAX_TRACKS, double &f, int &track)
{
  vector<const DCDCTrackHit*> cdctrackhits;
  vector<const DFDCPseudo*> fdcpseudos;
	
  // The DKinematicData object should be a DTrackCandidate, DTrackWireBased, or DParticle which
  // has associated objects for the hits
  kd->Get(cdctrackhits);
  kd->Get(fdcpseudos);

  // The track number is buried in the truth hit objects of type DMCTrackHit. These should be 
  // associated objects for the individual hit objects. We need to loop through them and
  // keep track of how many hits for each track number we find

  // CDC hits
  vector<int> cdc_track_no(MAX_TRACKS, 0);
  for(unsigned int i=0; i<cdctrackhits.size(); i++){
    vector<const DMCTrackHit*> mctrackhits;
    cdctrackhits[i]->Get(mctrackhits);
    if(mctrackhits.size()==0)continue;
    if(!mctrackhits[0]->primary)continue;
    int track = mctrackhits[0]->track;
    if(track>=0 && track<MAX_TRACKS)cdc_track_no[track]++;
    //cout << "cdc:(i,trackhitssize,mchitssize,TrackNo,NhitsforTrackNo):  " << "(" << i << "," << cdctrackhits.size() << "," << mctrackhits.size() << "," << track << "," << cdc_track_no[track] << ")" << endl;
    //cout << "cdc:(system,ptype,r,phi,z):  " << "(" << mctrackhits[0]->system << "," << mctrackhits[0]->ptype << "," << mctrackhits[0]->r << "," << mctrackhits[0]->phi << "," << mctrackhits[0]->z << ")" << endl;
  }
  // FDC hits
  vector<int> fdc_track_no(MAX_TRACKS, 0);
  for(unsigned int i=0; i<fdcpseudos.size(); i++){
    vector<const DMCTrackHit*> mctrackhits;
    fdcpseudos[i]->Get(mctrackhits);
    if(mctrackhits.size()==0)continue;
    if(!mctrackhits[0]->primary)continue;
    int track = mctrackhits[0]->track;
    if(track>=0 && track<MAX_TRACKS)fdc_track_no[track]++;
    //cout << "fdc:(i,trackhitssize,mchitssize,TrackNo,NhitsforTrackNo):  " << "(" << i << "," << fdcpseudos.size() << "," << mctrackhits.size() << "," << track << "," << fdc_track_no[track] << ")" << endl;
    //cout << "fdc:(system,ptype,r,phi,z):  " << "(" << mctrackhits[0]->system << "," << mctrackhits[0]->ptype << "," << mctrackhits[0]->r << "," << mctrackhits[0]->phi << "," << mctrackhits[0]->z << ")" << endl;
  }
	
  // Find track number with most wires hit
  int track_with_max_hits = 0;
  int tot_hits_max = cdc_track_no[0] + fdc_track_no[0];
  for(int i=1; i<MAX_TRACKS; i++){
    int tot_hits = cdc_track_no[i] + fdc_track_no[i];
    if(tot_hits > tot_hits_max){
      track_with_max_hits=i;
      tot_hits_max = tot_hits;
    }
    //cout << "tot_hits_max: " << tot_hits_max << endl;
    //cout << "track_with_max_hits: " << track_with_max_hits << endl;
  }
	
  int Ncdc = cdc_track_no[track_with_max_hits];
  int Nfdc = fdc_track_no[track_with_max_hits];

  // total fraction of reconstructed hits that match truth hits
  if (cdctrackhits.size()+fdcpseudos.size()) f = 1.*(Ncdc+Nfdc)/(cdctrackhits.size()+fdcpseudos.size());
  //cout << "(Ncdc(match),Nfdc(match),Ncdc(recon),Nfdc(recon)):  " << "(" << Ncdc << "," << Nfdc << "," << cdctrackhits.size() << "," << fdcpseudos.size() << ")" << endl;
  //if(DEBUG_HISTS)hitMatchFOM->Fill(f);

  // If there are no hits on this track, then we really should report
  // a "non-track" (i.e. track=-1)
  track = tot_hits_max>0 ? track_with_max_hits:-1;
}

//---------------------------------
// GetThrownIndex2
//---------------------------------
void DEventProcessor_pid_hists::GetThrownIndex2(const DKinematicData *kd,vector<const DMCThrown*>mcthrowns,int &track)
{
  double min_chisq=1.0E8;
  for(unsigned int i=0; i<mcthrowns.size(); i++){
    const DMCThrown *thrown = mcthrowns[i];
    double chisq = GetMomChisq(kd,thrown);
    if(chisq<min_chisq){track = thrown->myid; min_chisq = chisq;}
  }
}

double DEventProcessor_pid_hists::GetMomChisq(const DKinematicData *kd,const DKinematicData *kd_thrown)
{
  DVector3 mom = kd->momentum();
  DVector3 mc_mom = kd_thrown->momentum();
  double delta_pt = (mom.Pt()-mc_mom.Pt())/mc_mom.Pt();
  double delta_theta = (mom.Theta() - mc_mom.Theta())*1000.0; // in milliradians
  double delta_phi = (mom.Phi() - mc_mom.Phi())*1000.0; // in milliradians
  double err_pt = 0.04; double err_theta = 20.0; double err_phi = 20.0;
  double pull_pt = delta_pt/err_pt; double pull_theta = delta_theta/err_theta; double pull_phi = delta_phi/err_phi;
  double chisq = pull_pt*pull_pt + pull_theta*pull_theta + pull_phi*pull_phi;
  /*double pull_px = (mom.Px()-mc_mom.Px())/sqrt(kd->errorMatrix()(0,0));
    double pull_py = (mom.Py()-mc_mom.Py())/sqrt(kd->errorMatrix()(1,1));
    double pull_pz = (mom.Pz()-mc_mom.Pz())/sqrt(kd->errorMatrix()(2,2));
    double chisq2 = pull_px*pull_px + pull_py*pull_py + pull_pz*pull_pz;*/
  return chisq;
}

void DEventProcessor_pid_hists::printInfo(const DKinematicData *kd,int event,int track,int N,TString string)
{
  if (track==0) cout << "Kinematic info of " << string << " in event " << event << ":" << endl;
  int precision = 4; 
  DLorentzVector mom = kd->lorentzMomentum();
  cout << "track " << track << endl;
  cout << "(E,Px,Py,Pz,P,M)=" << setprecision(precision) << "(" << mom.E() << "," << mom.Px() << "," << mom.Py() << "," << mom.Pz() <<
    "," << mom.P() << "," << mom.M() << ")" << " and "; 
  cout << "(Beta,Pt,Theta,Phi)=" << setprecision(precision) << "(" << mom.Beta() << "," << mom.Pt() << "," << mom.Theta()*TMath::RadToDeg() << "," << mom.Phi()*TMath::RadToDeg() << ")" << "." << endl; 
  if (string!="thrown") cout << "(tofBeta,deltaBeta,TOF,dEdx)=" << setprecision(precision) << "(" << kd->measuredBeta() << "," << kd->deltaBeta() << "," << kd->TOF() << "," << kd->dEdx() << ")" << "." << endl; 
  if (string!="thrown") cout << "(deltaInvBeta,deltaTime,pathLength,dEdx)=" << setprecision(precision) << "(" << kd->deltaInvBeta() << "," << kd->pathLength()*kd->deltaInvBeta()/SPEED_OF_LIGHT << "," << kd->pathLength() << "," << kd->dEdx() << ")" << "." << endl; 
  if (track==N-1) cout << endl;
}

int DEventProcessor_pid_hists::GetType(const DTrackTimeBased *track)
{
  // Rely on the mass of the track to decide the type. 
  int type = track->charge()<0.0 ? 9:8; // initialize to pi-(=9) or pi+(=8)
  if (track->mass() > 0.250) { 
    type = track->charge()<0.0 ? 12:11; //  K-(=12) or K+(=11)
  }
  if(fabs(track->mass() - 0.93827)<0.100)type=14;

  return type;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_pid_hists::erun(void)
{

  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_pid_hists::fini(void)
{
  cout << "NgoodPid=" << NgoodPid << endl;
  cout << "NbadRecon=" << NbadRecon<< endl;
  cout << "NpipTruth_protonPid=" << NpipTruth_protonPid<< endl;
  cout << "NpipTruth_KpPid=" << NpipTruth_KpPid<< endl;
  cout << "NprotonTruth_pipPid=" << NprotonTruth_pipPid<< endl;
  cout << "NprotonTruth_KpPid=" << NprotonTruth_KpPid<< endl;
  cout << "NKpTruth_pipPid=" << NKpTruth_pipPid<< endl;
  cout << "NKpTruth_protonPid=" << NKpTruth_protonPid<< endl;
  cout << "NpimTruth_KmPid=" << NpimTruth_KmPid<< endl;
  cout << "NKmTruth_pimPid=" << NKmTruth_pimPid<< endl;
 
  return NOERROR;
}

