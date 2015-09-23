// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TEllipse.h>
#include <TArc.h>
#include <TBox.h>
#include <TLine.h>
#include <TText.h>
#include <TVector3.h>
#include <TColor.h>

#include <particleType.h>
#include "hdview2.h"
#include "hdv_mainframe.h"
#include "hdv_debugerframe.h"
#include "MyProcessor.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DQuickFit.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "TRACKING/DTrackCandidate_factory.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DReferenceTrajectory.h"
#include "JANA/JGeometry.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "FCAL/DFCALHit.h"
#include "TOF/DTOFGeometry.h"
#include "TOF/DTOFHit.h"
#include "TOF/DTOFTDCDigiHit.h" 
#include "TOF/DTOFDigiHit.h"
#include "TOF/DTOFPaddleHit.h"
#include "TOF/DTOFPoint.h"
#include "FDC/DFDCGeometry.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "FDC/DFDCIntersection.h"
#include "HDGEOMETRY/DGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALHit.h"
#include "PID/DNeutralParticle.h"
#include "PID/DNeutralShower.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALIncidentParticle.h"
#include "TOF/DTOFPoint.h"
#include "START_COUNTER/DSCHit.h"
#include "DVector2.h"

extern hdv_mainframe *hdvmf;

// These are declared in hdv_mainframe.cc, but as static so we need to do it here as well (yechh!)
static float FCAL_Zmin = 622.8;
static float FCAL_Rmin = 6.0;
static float FCAL_Rmax = 212.0/2.0;
static float BCAL_Rmin = 65.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmin = 212.0 - BCAL_Zlen/2.0;


static vector<vector <DFDCWire *> >fdcwires;


bool DMCTrajectoryPoint_track_cmp(const DMCTrajectoryPoint *a,const DMCTrajectoryPoint *b){
	// sort by track number and then by particle type, then by z-coordinate
	// (yes, I saw the same track number have different particle types!)
	if(a->track != b->track)return a->track < b->track;
	if(a->part  != b->part )return a->part < b->part;
	return a->z < b->z;
}


MyProcessor *gMYPROC=NULL;

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	Bfield = NULL;
	loop = NULL;

	// Tell factory to keep around a few density histos
	//gPARMS->SetParameter("TRKFIND:MAX_DEBUG_BUFFERS",	16);
	
	RMAX_INTERIOR = 65.0;
	RMAX_EXTERIOR = 88.0;
	gPARMS->SetDefaultParameter("RT:RMAX_INTERIOR",	RMAX_INTERIOR, "cm track drawing Rmax inside solenoid region");
	gPARMS->SetDefaultParameter("RT:RMAX_EXTERIOR",	RMAX_EXTERIOR, "cm track drawing Rmax outside solenoid region");
	
	gMYPROC = this;
}

//------------------------------------------------------------------
// ~MyProcessor 
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{

}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// Make sure detectors have been drawn
	//if(!drew_detectors)DrawDetectors();
  
  

	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){

		vector<string> facnames;
		loops[0]->GetFactoryNames(facnames);

		hdvmf = new hdv_mainframe(gClient->GetRoot(), 1400, 700);
		hdvmf->SetCandidateFactories(facnames);
		hdvmf->SetWireBasedTrackFactories(facnames);
		hdvmf->SetTimeBasedTrackFactories(facnames);
		hdvmf->SetReconstructedFactories(facnames);
		hdvmf->SetChargedTrackFactories(facnames);
		fulllistmf = hdvmf->GetFullListFrame();
		debugermf = hdvmf->GetDebugerFrame();
		BCALHitCanvas = hdvmf->GetBcalDispFrame();

		if (BCALHitCanvas){
		  BCALHitMatrixU = new TH2F("BCALHitMatrixU","BCAL Hits Upstream",  48*4+2, -1.5, 192.5, 10, 0., 10.);
		  BCALHitMatrixD = new TH2F("BCALHitMatrixD","BCAL Hits Downstream",48*4+2, -1.5, 192.5, 10, 0., 10.);
		  BCALParticles = new TH2F("BCALParticles","BCAL Hits Downstream",(48*4+2)*4, -1.87, 361.87, 1, 0., 1.);
		  BCALHitMatrixU->SetStats(0);
		  BCALHitMatrixD->SetStats(0);
		  BCALParticles->SetStats(0);
		  BCALHitMatrixU->GetXaxis()->SetTitle("Sector number");
		  BCALHitMatrixD->GetXaxis()->SetTitle("Sector number");
		  BCALParticles->GetXaxis()->SetTitle("Phi angle [deg]");
		}
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventloop, int runnumber)
{

	// Read in Magnetic field map
	DApplication* dapp = dynamic_cast<DApplication*>(eventloop->GetJApplication());
	Bfield = dapp->GetBfield(runnumber);
	const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
	dgeom->GetFDCWires(fdcwires);

	RootGeom = dapp->GetRootGeom(runnumber);
	geom = dapp->GetDGeometry(runnumber);
	
	MATERIAL_MAP_MODEL="DGeometry";
	gPARMS->SetDefaultParameter("TRKFIT:MATERIAL_MAP_MODEL",			MATERIAL_MAP_MODEL);

	eventloop->GetCalib("PID/photon_track_matching", photon_track_matching);
	DELTA_R_FCAL = photon_track_matching["DELTA_R_FCAL"];

	return NOERROR;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
	if(!eventLoop)return NOERROR;
	loop = eventLoop;
	last_jevent.FreeEvent();
	last_jevent = loop->GetJEvent();
	
	string source = "<no source>";
	if(last_jevent.GetJEventSource())source = last_jevent.GetJEventSource()->GetSourceName();
	
	cout<<"----------- New Event "<<eventnumber<<" -------------"<<endl;
	hdvmf->SetEvent(eventnumber);
	hdvmf->SetSource(source.c_str());
	hdvmf->DoMyRedraw();	

	return NOERROR;
}

//------------------------------------------------------------------
// FillGraphics 
//------------------------------------------------------------------
void MyProcessor::FillGraphics(void)
{
	/// Create "graphics" objects for this event given the current GUI settings.
	///
	/// This method will create DGraphicSet objects that represent tracks, hits,
	/// and showers for the event. It creates objects for both hits and
	/// reconstructed entities. The "graphics" objects created here are
	/// really just collections of 3D space points with flags indicating
	/// whether they should be drawn as markers or lines and with what
	/// color and size. The actual graphics objects are created for the
	/// various views of the detector in hdv_mainframe.
		
	graphics.clear();
	graphics_xyA.clear(); // The objects placed in these will be deleted by hdv_mainframe
	graphics_xyB.clear(); // The objects placed in these will be deleted by hdv_mainframe
	graphics_xz.clear();  // The objects placed in these will be deleted by hdv_mainframe
	graphics_yz.clear();  // The objects placed in these will be deleted by hdv_mainframe
	graphics_tof_hits.clear();  // The objects placed in these will be deleted by hdv_mainframe

	if(!loop)return;

	vector<const DSCHit *>schits;
	loop->Get(schits);
	vector<const DTrackCandidate*> trCand;
	loop->Get(trCand);
	vector<const DTrackTimeBased*> trTB;
	loop->Get(trTB);
	vector<const DTrackWireBased*> trWB;
	loop->Get(trWB);
	hdv_debugerframe *p = hdvmf->GetDebugerFrame();
	p->SetNTrCand(trCand.size());
	p->SetNTrWireBased(trWB.size());
	p->SetNTrTimeBased(trTB.size());

	if (BCALHitCanvas) {
	  vector<const DBCALHit*> locBcalHits;
	  loop->Get(locBcalHits);
	  BCALHitMatrixU->Reset();
	  BCALHitMatrixD->Reset();
	  for (unsigned int k=0;k<locBcalHits.size();k++){
	    
	    const DBCALHit* hit = locBcalHits[k];
	    float idxY = (float)hit->layer-1;
	    float idxX = (float) (hit->sector-1 + (hit->module-1)*4);
	    if (hit->end == DBCALGeometry::kUpstream){
	      if (hit->layer==1){
		BCALHitMatrixU->Fill(idxX,idxY,hit->E);
	      } else if (hit->layer==2){
		BCALHitMatrixU->Fill(idxX,idxY,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+1.,hit->E);	      
	      } else if (hit->layer==3){
		BCALHitMatrixU->Fill(idxX,idxY+1,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+2.,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+3.,hit->E);	      
	      } else if (hit->layer==4){
		BCALHitMatrixU->Fill(idxX,idxY+3,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+4.,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+5.,hit->E);	      
		BCALHitMatrixU->Fill(idxX,idxY+6.,hit->E);	      
	      }
	    } else {
	      if (hit->layer==1){
		BCALHitMatrixD->Fill(idxX,idxY,hit->E);
	      } else if (hit->layer==2){
		BCALHitMatrixD->Fill(idxX,idxY,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+1.,hit->E);	      
	      } else if (hit->layer==3){
		BCALHitMatrixD->Fill(idxX,idxY+1,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+2.,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+3.,hit->E);	      
	      } else if (hit->layer==4){
		BCALHitMatrixD->Fill(idxX,idxY+3,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+4.,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+5.,hit->E);	      
		BCALHitMatrixD->Fill(idxX,idxY+6.,hit->E);	      
	      }
	    }
	  }

	  vector<const DBCALIncidentParticle*> locBcalParticles;
	  loop->Get(locBcalParticles);
	  BCALParticles->Reset();
	  BCALPLables.clear();
	  for (unsigned int k=0;k<locBcalParticles.size();k++){	    
	    const DBCALIncidentParticle* part = locBcalParticles[k];
	    
	    float p = TMath::Sqrt(part->px*part->px +  part->py*part->py +  part->pz*part->pz);
	    float phi=999;
	    if (part->x!=0){
	      phi = TMath::ATan(TMath::Abs(part->y/part->x));
	      //cout<<phi<<"   "<<part->y<<" / "<< part->x<<endl;
	      if (part->y>0){
		if (part->x<0.){
		  phi = 3.1415926 - phi;
		}
	      } else {
		if (part->x<0){
		  phi += 3.1415926;
		} else {
		  phi = 3.1415926*2. - phi;
		}
	      }
	      
	      phi = phi*180./3.1415926;
	    }
	    //cout<<phi<<"  "<<p<<endl;
	    BCALParticles->Fill(phi,0.5,p);
	    char l[20];
	    sprintf(l,"%d",part->ptype);
	    TText *t = new TText(phi,1.01,l);
	    t->SetTextSize(0.08);
	    t->SetTextFont(72);
	    t->SetTextAlign(21);
	    BCALPLables.push_back(t);
	  }
 	   
	  BCALHitCanvas->Clear();
	  BCALHitCanvas->Divide(1,3);
	  BCALHitCanvas->cd(1);
	  if (BCALHitMatrixU->GetMaximum() > 1000) {
	    BCALHitMatrixU->Scale(0.001);  // Scale histogram to GeV
	    BCALHitMatrixU->GetZaxis()->SetTitle("Energy  (GeV)");
	    printf("scaling to GeV\n");
	  } else {
	    BCALHitMatrixU->GetZaxis()->SetTitle("Energy  (MeV)");
	  }
	  BCALHitMatrixU->Draw("colz");
	  BCALHitCanvas->cd(2);
	  if (BCALHitMatrixD->GetMaximum() > 1000) {
	    BCALHitMatrixD->Scale(0.001);  // Scale histogram to GeV
	    BCALHitMatrixD->GetZaxis()->SetTitle("Energy  (GeV)");
	  } else {
	    BCALHitMatrixU->GetZaxis()->SetTitle("Energy  (MeV)");
	  }
	  BCALHitMatrixD->Draw("colz");
	  BCALHitCanvas->cd(3);
	  BCALParticles->Draw("colz");
	  for (unsigned int n=0;n<BCALPLables.size();n++){
	    BCALPLables[n]->Draw("same");
	  }
	  BCALHitCanvas->Update();
	}
	
	
	// BCAL hits
	if(hdvmf->GetCheckButton("bcal")){
	  vector<const DBCALHit*> bcalhits;
	  loop->Get(bcalhits);
	  
	  for(unsigned int i=0; i<bcalhits.size(); i++){
	    const DBCALHit *hit = bcalhits[i];
	    TPolyLine *poly = hdvmf->GetBCALPolyLine(hit->module, hit->layer, hit->sector);

	    if(!poly)continue;
	   
	    // The aim is to have a log scale in energy
	    double E = 1000.0*hit->E;      // Energy in MeV
	    double logE = log10(E);      
	    // 3 = 1 GeV, 0 = 1 MeV, use range 0 through 4
	    // 0-1 White-Yellow
	    // 1-2 Yellow-Red
	    // 2-3 Red-Cyan
	    // 3-4 Cyan-Blue

	    float r,g,b;
	    if (E<1){
	      r = 1.;
	      g = 1.;
	      b = 1.;
	    } else {
	      if (logE<1){
		r = 1.;
		g = 1.;
		b = 1.-logE;
	      } else {
		if (logE<2){
		  r = 1.;
		  g = 1.-(logE-1);
		  b = 0.;
		} else {
		  if (logE<3){
		    r = 1.;
		    g = 0.;
		    b = 1.-(logE-2);
		  } else {
		    if (logE<4){
		      r = 1.-(logE-3);
		      g = 0.;
		      b = 1.;
		    } else {
		      r = 0;
		      g = 1.;
		      b = 0.;
		      //printf("High BCAL cell reconstructed energy: E=%.1f MeV\n",E);
		    }
		  }
		}
	      }
	    }
	    if (r<0||g<0||b<0||r>1||g>1||b>1) printf("color error (r,g,b)=(%f,%f,%f)\n",r,g,b);

	    poly->SetFillColor(TColor::GetColor(r,g,b));
	    poly->SetLineColor(TColor::GetColor(r,g,b));
	    poly->SetLineWidth(1);
	    poly->SetFillStyle(1001);
	  }
	}	

	// FCAL hits
	if(hdvmf->GetCheckButton("fcal")){
	  vector<const DFCALHit*> fcalhits;
	  loop->Get(fcalhits);
	  
	  for(unsigned int i=0; i<fcalhits.size(); i++){
	    const DFCALHit *hit = fcalhits[i];
	    TPolyLine *poly = hdvmf->GetFCALPolyLine(hit->x, hit->y);
	    if(!poly)continue;
	    
#if 0			
	    double a = hit->E/0.005;
	    double f = sqrt(a>1.0 ? 1.0:a<0.0 ? 0.0:a);
	    double grey = 0.8;
	    double s = 1.0 - f;
	    
	    float r = s*grey;
	    float g = s*grey;
	    float b = f*(1.0-grey) + grey;
#endif

	    // The aim is to have a log scale in energy (see BCAL)
	    double E = 1000*hit->E;      // Change Energy to MeV
	    if(E<0.0) continue;
	    double logE = log10(E);      

	    float r,g,b;
	    if (logE<0){
	      r = 1.;
	      g = 1.;
	      b = 1.;
	    } else {
	      if (logE<1){
		r = 1.;
		g = 1.;
		b = 1.-logE;
	      } else {
		if (logE<2){
		  r = 1.;
		  g = 1.-(logE-1);
		  b = 0.;
		} else {
		  if (logE<3){
		    r = 1.;
		    g = 0.;
		    b = 1.-(logE-2);
		  } else {
		    if (logE<4){
		      r = 1.-(logE-3);
		      g = 0.;
		      b = 1.;
		    } else {
		      r = 0;
		      g = 0;
		      b = 0;
		    }
		  }
		}
	      }
	    }
	    poly->SetFillColor(TColor::GetColor(r,g,b));
	  }
	}
	// TOF hits
	if(hdvmf->GetCheckButton("tof")){

	  double hit_north[45];
	  double hit_south[45];
	  double hit_up[45];
	  double hit_down[45];

	  memset(hit_north,0,sizeof(hit_north));
	  memset(hit_south,0,sizeof(hit_south));
	  memset(hit_up,0,sizeof(hit_up));
	  memset(hit_down,0,sizeof(hit_down));

          vector<const DTOFHit*> tofhits;
          loop->Get(tofhits);

          for(unsigned int i=0; i<tofhits.size(); i++){
            const DTOFHit *tof_hit = tofhits[i];

	    int plane = tof_hit->plane;
            int bar = tof_hit->bar;
            int end = tof_hit->end;
            float t = tof_hit->t;
	    
	    int translate_side;
	    TPolyLine *pmtPline;


	    // Flash the PMTs that do have fADC hits

	    /*
	    double dE_padd = 0.2/5.2E5 * 40000;
	    double thold = 0.2/5.2E5 * 115;
	    if (tof_hit->has_fADC && (tof_hit->dE - dE_padd > thold)){
	    */

	    if (tof_hit->has_fADC){
	      switch(plane)
		{
		case 0:
		  if(end == 1){
		    //cout << "Down : " << bar << endl;
		    translate_side = 0;
		    pmtPline = hdvmf->GetTOFPolyLine(translate_side, bar);
		    pmtPline->SetFillColor(2);
		  }
		  else if(end == 0){
		    //cout << "Up : " << bar << endl;
		    translate_side = 2;
		    pmtPline = hdvmf->GetTOFPolyLine(translate_side, bar);
		    pmtPline->SetFillColor(2);
		}
		  else{
		  cerr << "Out of range TOF end" << endl;
		  }
		  break;
		case 1:
		  if(end == 0){
		    //cout << "North : " << bar << endl;
		    translate_side = 3;
		    pmtPline = hdvmf->GetTOFPolyLine(translate_side, bar);
		    pmtPline->SetFillColor(2);
		  }
		  else if(end == 1){
		    //cout << "South : " << bar << endl;
		    translate_side = 1;
		    pmtPline = hdvmf->GetTOFPolyLine(translate_side, bar);
		    pmtPline->SetFillColor(2);
		  }
		  else{
		    cerr << "Out of range TOF end" << endl;
		  }
		  break;
		default:
		  cerr << "Out of range TOF plane" << endl;
		  exit(0);
		} // close the switch plane loop
	    } // if for the fADC

	    // Draw the position from the events that do have tdc hits 
	    // with the current status of the TOFHit object those hits appear with no match from the fADC
	    //Float_t hit_dist;

	    if (tof_hit->has_TDC){
	      switch(plane)
		{
		case 0:
		  if(end == 1){
		    if (hit_down[bar]<=0 || (t < hit_down[bar]) ){
		      hit_down[bar] = t;
		    }
		  }
		  else if(end == 0){
		    if (hit_up[bar]<=0 || (t < hit_up[bar]) ){
		      hit_up[bar] = t;
		    }
		  }
		  else{
		  cerr << "Out of range TOF end" << endl;
		  }
		  break;
		case 1:
		  if(end == 0){
		    if (hit_north[bar]<=0 || (t < hit_north[bar]) ){
		      hit_north[bar] = t;
		    }
		  }
		  else if(end == 1){
		    if (hit_south[bar]<=0 || (t < hit_south[bar]) ){
		      hit_south[bar] = t;
		    }
		  }
		  else{
		    cerr << "Out of range TOF end" << endl;
		  }
		  break;
		default:
		  cerr << "Out of range TOF plane" << endl;
		  exit(0);
		}
	      
	    } // close the switch if there is a TDC Hit
	    
	  } // close the for TOFHit object

	  // Draw the TDC Points here

	  Float_t hit_dist;
	  Float_t distY_Horz = -126; // Horizontal plane start counting from the Bottom to Top
	  Float_t distX_Vert =  -126; // Vertical plane start counting from the South to North
	  int tdc_hits = 0;
	  
	  for(Int_t i_tdc = 1; i_tdc <= 44; i_tdc++){
	    if ( i_tdc == 20 || i_tdc == 21 || i_tdc == 24 || i_tdc == 25 ){
	      distY_Horz = distY_Horz + 1.5;
	      distX_Vert = distX_Vert + 1.5;
	    }
	    else{
	      distY_Horz = distY_Horz + 3.0;
	      distX_Vert = distX_Vert + 3.0;
	    }
	    if(hit_north[i_tdc] > 0 && hit_south[i_tdc] > 0){
	      hit_dist =  (15.2*(Float_t(hit_south[i_tdc] - hit_north[i_tdc])/2));
	      TArc *tdc_cir = new TArc(hit_dist,distY_Horz,2);
	      tdc_cir->SetFillColor(kGreen);

	      graphics_tof_hits.push_back(tdc_cir);
	      tdc_hits++;
	    }
	    if(hit_up[i_tdc] > 0 && hit_down[i_tdc] > 0){
	      hit_dist =  (15.2*(Float_t(hit_down[i_tdc] - hit_up[i_tdc])/2) );
	      TArc *tdc_cir = new TArc(distX_Vert,hit_dist,2);
	      tdc_cir->SetFillColor(kBlue);

	      graphics_tof_hits.push_back(tdc_cir);
	      tdc_hits++;
	    }
	    if ( i_tdc == 20 || i_tdc == 21 || i_tdc == 24 || i_tdc == 25 ){
	      distY_Horz = distY_Horz + 1.5;
	      distX_Vert = distX_Vert + 1.5;
	    }
	    else{
	      distY_Horz = distY_Horz + 3.0;
	      distX_Vert = distX_Vert + 3.0;
	    }
	  }

        } // close the if check button for the TOF

	// Start counter hits 
	for (unsigned int i=0;i<schits.size();i++){
	  DGraphicSet gset(6,kLine,2.0);
	  double r_start=7.7493;
	  double phi0=0.2094395*(schits[i]->sector-1);  // span 12 deg in phi
	  double phi1=0.2094395*(schits[i]->sector);
	  TVector3 point1(r_start*cos(phi0),r_start*sin(phi0),38.75); 
	  gset.points.push_back(point1); // upstream end of sctraight section of scint
	  TVector3 point2(r_start*cos(phi1),r_start*sin(phi1),38.75);
	  gset.points.push_back(point2); 
	  TVector3 point3(r_start*cos(phi1),r_start*sin(phi1),78.215);
	  gset.points.push_back(point3); // downstream end of sctraight section of scint
	  TVector3 point4(r_start*cos(phi0),r_start*sin(phi0),78.215);
	  gset.points.push_back(point4); 
	  TVector3 point5(r_start*cos(phi0),r_start*sin(phi0),38.75);
	  gset.points.push_back(point5);
	
	  /*
	  // ST dimensions
	  Double_t dtr = 1.74532925e-02;		// Conversion factor from degrees to radians
	  Double_t st_straight = 39.465;		// Distance of the straight section along scintillator path
	  Double_t st_arc_angle = 18.5;		// Angle of the bend arc
	  Double_t st_arc_radius = 12.0;		// Radius of the bend arc
	  Double_t st_to_beam = 7.74926;		// Distance from beam to bottom side of scintillator paddle
	  Double_t st_len_cone = 16.056;		// Length of the cone section along scintillator path
	  Double_t st_thickness = 0.3;		// Thickness of the scintillator paddles
	  Double_t st_beam_to_arc_radius = 4.25;	// Distance from the beam line to the arc radius

	  // Nose Arrays
	  Double_t tp_nose_z[5];
	  Double_t tp_nose_y[5];
	  Double_t bm_nose_z[5];
	  Double_t bm_nose_y[5];
	
	  // Offsets for Hall coordinates
	  Double_t z_center = 65.0;	// Target Center (0,0,65)
	  Double_t us_end_pt = -26.25;    // Distance of upstream end relative to target
	  Double_t ds_end_pt = 97.4;	// Distance of downstream end relative

	  // Top start counter paddle straight section
	  TPolyLine *top_paddle_straight;
	  Double_t tp_straight_z[5] = {us_end_pt + z_center, us_end_pt + z_center, us_end_pt + st_straight + z_center, us_end_pt + st_straight + z_center, us_end_pt + z_center};
	  Double_t tp_straight_y[5] = {st_to_beam, st_to_beam + st_thickness, st_to_beam + st_thickness, st_to_beam, st_to_beam};
	  */
	  graphics.push_back(gset);
	}
	  
	// CDC hits
	if(hdvmf->GetCheckButton("cdc")){
		vector<const DCDCTrackHit*> cdctrackhits;
		loop->Get(cdctrackhits);
		
		for(unsigned int i=0; i<cdctrackhits.size(); i++){
			const DCDCWire *wire = cdctrackhits[i]->wire;
			
			int color = (cdctrackhits[i]->tdrift>-20 && cdctrackhits[i]->tdrift<800) ? kCyan:kYellow;
			DGraphicSet gset(color, kLine, 1.0);
			DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
			TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			dpoint=wire->origin+(wire->L/2.0)*wire->udir;
			tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			graphics.push_back(gset);
			
			// Rings for drift times.
			// NOTE: These are not perfect since they still have TOF in them
			if(hdvmf->GetCheckButton("cdcdrift") && fabs(wire->stereo)<0.05){
				double x = wire->origin.X();
				double y = wire->origin.Y();
				double dist1 = cdctrackhits[i]->dist;
				TEllipse *e = new TEllipse(x, y, dist1, dist1);
				e->SetLineColor(38);
				e->SetFillStyle(0);
				graphics_xyA.push_back(e);

				double dist2 = dist1 - 4.0*55.0E-4; // what if our TOF was 4ns?
				e = new TEllipse(x, y, dist2, dist2);
				e->SetLineColor(38);
				e->SetFillStyle(0);
				graphics_xyA.push_back(e);
			}
		}
	}	

	// FDC wire
	if(hdvmf->GetCheckButton("fdcwire")){
		vector<const DFDCHit*> fdchits;
		loop->Get(fdchits);

		for(unsigned int i=0; i<fdchits.size(); i++){
			const DFDCHit *fdchit = fdchits[i];
			if(fdchit->type!=0)continue;
			const DFDCWire *wire =fdcwires[fdchit->gLayer-1][fdchit->element-1];
			if(!wire){
				_DBG_<<"Couldn't find wire for gLayer="<<fdchit->gLayer<<" and element="<<fdchit->element<<endl;
				continue;
			}
			
			// Wire
			int color = (fdchit->t>-50 && fdchit->t<2000) ? kCyan:kYellow;
			DGraphicSet gset(color, kLine, 1.0);
			DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
			TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			dpoint=wire->origin+(wire->L/2.0)*wire->udir;
			tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			graphics.push_back(gset);
		}		
	}

	// FDC intersection hits
	if(hdvmf->GetCheckButton("fdcintersection")){
		vector<const DFDCIntersection*> fdcints;
		loop->Get(fdcints);
		DGraphicSet gsetp(46, kMarker, 0.5);
		
		for(unsigned int i=0; i<fdcints.size(); i++){
		  TVector3 tpos(fdcints[i]->pos.X(),fdcints[i]->pos.Y(),
				fdcints[i]->pos.Z());
		  gsetp.points.push_back(tpos);
		}
		graphics.push_back(gsetp);
	}
	
	// FDC psuedo hits
	if(hdvmf->GetCheckButton("fdcpseudo")){
		vector<const DFDCPseudo*> fdcpseudos;
		loop->Get(fdcpseudos);
		DGraphicSet gsetp(38, kMarker, 0.5);
		
		for(unsigned int i=0; i<fdcpseudos.size(); i++){
			const DFDCWire *wire = fdcpseudos[i]->wire;
			
			// Pseudo point
			TVector3 pos(fdcpseudos[i]->xy.X(), 
				     fdcpseudos[i]->xy.Y(), wire->origin.Z());
			gsetp.points.push_back(pos);
		}
		graphics.push_back(gsetp);
	}

	// DMCThrown
	if(hdvmf->GetCheckButton("thrown")){
		vector<const DMCThrown*> mcthrown;
		loop->Get(mcthrown);
		for(unsigned int i=0; i<mcthrown.size(); i++){
			int color=14;
			double size=1.5;
			if(mcthrown[i]->charge()==0.0) color = kGreen;
			if(mcthrown[i]->charge() >0.0) color = kBlue;
			if(mcthrown[i]->charge() <0.0) color = kRed;
			switch(mcthrown[i]->type){
				case Gamma:
				case Positron:
				case Electron:
					size = 1.0;
					break;
				case Pi0:
				case PiPlus:
				case PiMinus:
					size = 2.0;
					break;
				case Neutron:
				case Proton:
				case AntiProton:
					size = 3.0;
					break;
			}
			AddKinematicDataTrack(mcthrown[i], color, size);
		}
	}
	
	// CDC Truth points
	if(hdvmf->GetCheckButton("cdctruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(12, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_CDC)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}
	
	// FDC Truth points
	if(hdvmf->GetCheckButton("fdctruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(12, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_FDC)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// Track Hits for Track Candidates and Candidate trajectory in Debuger Window
	for(unsigned int n=0; n<trCand.size(); n++){
	  if (n>9)
	    break;
	  char str1[128];
	  sprintf(str1,"Candidate%d",n+1);
	  
	  if(hdvmf->GetCheckButton(str1)){	

	    int color = n+1;
	    if (color > 4)
	      color++;
	    if (color > 6)
	      color++;	      
	    
	    AddKinematicDataTrack(trCand[n], color, 1.5);

	    vector<const DCDCTrackHit*> cdctrackhits;
	    trCand[n]->Get(cdctrackhits);
	    for(unsigned int i=0; i<cdctrackhits.size(); i++){
	      const DCDCWire *wire = cdctrackhits[i]->wire;
	      DGraphicSet gset(color, kLine, 1.0);
	      DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
	      TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      dpoint=wire->origin+(wire->L/2.0)*wire->udir;
	      tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      graphics.push_back(gset);
	      
	    } // end loop of cdc hits of track candidate
	    vector<const DFDCPseudo*> fdcpseudos;
	    trCand[n]->Get(fdcpseudos);
	    DGraphicSet gsetp(color, kMarker, 0.5);
	    
	    for(unsigned int i=0; i<fdcpseudos.size(); i++){
	      const DFDCWire *wire = fdcpseudos[i]->wire;
	      
	      // Pseudo point
	      TVector3 pos(fdcpseudos[i]->xy.X(), 
			   fdcpseudos[i]->xy.Y(), wire->origin.Z());
	      gsetp.points.push_back(pos);
	    }	
	    graphics.push_back(gsetp);
	    
	  }
	}

	// Wire Based Track Hits and trajectory for Debuger Window
	for(unsigned int n=0; n<trWB.size(); n++){
	  if (n>=MaxWireTracks)
	    break;
	  char str1[128];
	  sprintf(str1,"WireBased%d",n+1);
	  
	  if(hdvmf->GetCheckButton(str1)){	

	    int color = trWB[n]->candidateid;
	    if (color > 4)
	      color++;
	    if (color > 6)
	      color++;	      

	    AddKinematicDataTrack(trWB[n], color, 1.5);

	    vector<const DCDCTrackHit*> cdctrackhits;
	    trWB[n]->Get(cdctrackhits);
	    for(unsigned int i=0; i<cdctrackhits.size(); i++){
	      const DCDCWire *wire = cdctrackhits[i]->wire;
	      DGraphicSet gset(color, kLine, 1.0);
	      DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
	      TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      dpoint=wire->origin+(wire->L/2.0)*wire->udir;
	      tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      graphics.push_back(gset);
	      
	    } // end loop of cdc hits of track candidate
	    vector<const DFDCPseudo*> fdcpseudos;
	    trWB[n]->Get(fdcpseudos);
	    DGraphicSet gsetp(color, kMarker, 0.5);
	    
	    for(unsigned int i=0; i<fdcpseudos.size(); i++){
	      const DFDCWire *wire = fdcpseudos[i]->wire;
	      
	      // Pseudo point
	      TVector3 pos(fdcpseudos[i]->xy.X(), 
			   fdcpseudos[i]->xy.Y(), wire->origin.Z());
	      gsetp.points.push_back(pos);
	    }
	    graphics.push_back(gsetp);
	    
	  }
	}

	// Time Based Track Hits and trajectory for Debuger Window
	for(unsigned int n=0; n<trTB.size(); n++){
	  if (n>=MaxTimeTracks)
	    break;
	  char str1[128];
	  sprintf(str1,"TimeBased%d",n+1);
	  
	  if(hdvmf->GetCheckButton(str1)){	

	    int color = trTB[n]->candidateid;
	    if (color > 4)
	      color++;
	    if (color > 6)
	      color++;	      

	    AddKinematicDataTrack(trTB[n], color, 1.5);

	    vector<const DCDCTrackHit*> cdctrackhits;
	    trTB[n]->Get(cdctrackhits);
	    for(unsigned int i=0; i<cdctrackhits.size(); i++){
	      const DCDCWire *wire = cdctrackhits[i]->wire;
	      DGraphicSet gset(color, kLine, 1.0);
	      DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
	      TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      dpoint=wire->origin+(wire->L/2.0)*wire->udir;
	      tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
	      gset.points.push_back(tpoint);
	      graphics.push_back(gset);
	      
	    } // end loop of cdc hits of track candidate
	    vector<const DFDCPseudo*> fdcpseudos;
	    trTB[n]->Get(fdcpseudos);
	    DGraphicSet gsetp(color, kMarker, 0.5);
	    
	    for(unsigned int i=0; i<fdcpseudos.size(); i++){
	      const DFDCWire *wire = fdcpseudos[i]->wire;
	      
	      // Pseudo point
	      TVector3 pos(fdcpseudos[i]->xy.X(), 
			   fdcpseudos[i]->xy.Y(), wire->origin.Z());
	      gsetp.points.push_back(pos);
	    }
	    graphics.push_back(gsetp);
	    
	  }
	}
	
	// TOF reconstructed points 
	if (hdvmf->GetCheckButton("tof")){
	  vector<const DTOFPoint *>tofpoints;
	  loop->Get(tofpoints);
	  DGraphicSet gset(kRed, kMarker, 0.5);
	  for(unsigned int i=0; i<tofpoints.size(); i++){
	    const DTOFPoint *hit = tofpoints[i];
	    TVector3 pos(hit->pos.x(),hit->pos.y(),hit->pos.z());
	    gset.points.push_back(pos);
	  }
	  graphics.push_back(gset);
	}
	
	// TOF Truth points
	if(hdvmf->GetCheckButton("toftruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(kBlack, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_TOF)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// BCAL Truth points
	if(hdvmf->GetCheckButton("bcaltruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(kBlack, kMarker, 1.0);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_BCAL)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);

			TMarker *m = new TMarker(pos.X(), pos.Y(), 2);	
			graphics_xyA.push_back(m);
		}
		//graphics.push_back(gset);
	}

	// FCAL Truth points
	if(hdvmf->GetCheckButton("fcaltruth")){
		vector<const DFCALGeometry*> fcalgeometries;
		vector<const DFCALHit*> mcfcalhits;
		loop->Get(fcalgeometries);
		loop->Get(mcfcalhits, "TRUTH");
		if(fcalgeometries.size()>0){
			const DFCALGeometry *fgeom = fcalgeometries[0];

			DGraphicSet gset(kBlack, kMarker, 0.25);
			for(unsigned int i=0; i<mcfcalhits.size(); i++){
				const DFCALHit *hit = mcfcalhits[i];

				DVector2 pos_face = fgeom->positionOnFace(hit->row, hit->column);
				TVector3 pos(pos_face.X(), pos_face.Y(), FCAL_Zmin);
				gset.points.push_back(pos);
				
				TMarker *m = new TMarker(pos.X(), pos.Y(), 2);	
				//m->SetColor(kGreen);
				//m->SetLineWidth(1);
				graphics_xyB.push_back(m);

				TMarker *m1 = new TMarker(pos.Z(), pos.X(), 2);	
				graphics_xz.push_back(m1);
				TMarker *m2 = new TMarker(pos.Z(), pos.Y(), 2);	
				graphics_yz.push_back(m2);
			}
			//graphics.push_back(gset);
		}
	}

	// BCAL reconstructed photons
	if(hdvmf->GetCheckButton("recon_photons_bcal")){
		vector<const DNeutralParticle*> neutrals;
		loop->Get(neutrals);
	
		DGraphicSet gset(kYellow+2, kMarker, 1.25);
		gset.marker_style=21;
		for(unsigned int i=0; i<neutrals.size(); i++){
		  vector<const DNeutralShower*> locNeutralShowers;
		  neutrals[i]->GetT(locNeutralShowers);
		  DetectorSystem_t locDetectorSystem = locNeutralShowers[0]->dDetectorSystem;
		  if(locDetectorSystem == SYS_BCAL){
		    TVector3 pos( locNeutralShowers[0]->dSpacetimeVertex.X(), 
				  locNeutralShowers[0]->dSpacetimeVertex.Y(), 
				  locNeutralShowers[0]->dSpacetimeVertex.Z());
		    gset.points.push_back(pos);
		    
		    double dist2 = 2.0 + 5.0*locNeutralShowers[0]->dEnergy;
		    TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
		    e->SetLineColor(kGreen);
		    e->SetFillStyle(0);
		    e->SetLineWidth(2);
		    graphics_xyA.push_back(e);
		  }
		}
		//graphics.push_back(gset);
	}

	// FCAL reconstructed photons
	if(hdvmf->GetCheckButton("recon_photons_fcal")){
		vector<const DNeutralParticle*> neutrals;
		loop->Get(neutrals);
		DGraphicSet gset(kOrange, kMarker, 1.25);
		gset.marker_style=2;
		for(unsigned int i=0; i<neutrals.size(); i++){
		  vector<const DNeutralShower*> locNeutralShowers;
		  neutrals[i]->GetT(locNeutralShowers);
		  DetectorSystem_t locDetectorSystem = locNeutralShowers[0]->dDetectorSystem;
		  if(locDetectorSystem == SYS_FCAL){
			
			TVector3 pos( locNeutralShowers[0]->dSpacetimeVertex.X(), 
				      locNeutralShowers[0]->dSpacetimeVertex.Y(), 
				      locNeutralShowers[0]->dSpacetimeVertex.Z());
			gset.points.push_back(pos);
			
			double dist2 = 2.0 + 10.0*locNeutralShowers[0]->dEnergy;
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(kGreen);
			e->SetFillStyle(0);
			e->SetLineWidth(2);
			graphics_xyB.push_back(e);

			TEllipse *e1 = new TEllipse(pos.Z(), pos.X(), dist2, dist2);
			e1->SetLineColor(kGreen);
			e1->SetFillStyle(0);
			e1->SetLineWidth(2);
			graphics_xz.push_back(e1);
			TEllipse *e2 = new TEllipse(pos.Z(), pos.Y(), dist2, dist2);
			e2->SetLineColor(kGreen);
			e2->SetFillStyle(0);
			e2->SetLineWidth(2);
			graphics_yz.push_back(e2);

		  }
		}
		//graphics.push_back(gset);
	}

	// Reconstructed photons matched with tracks
	if(hdvmf->GetCheckButton("recon_photons_track_match")){
		vector<const DChargedTrack*> ctracks;
		loop->Get(ctracks);
		for(unsigned int i=0; i<ctracks.size(); i++){
		  const DChargedTrack *locCTrack = ctracks[i];
		  vector<const DNeutralShower*> locNeutralShowers;
		  locCTrack->GetT(locNeutralShowers);

		  if (!locNeutralShowers.size()) continue;
		  
		  
		  // Decide if this hit BCAL of FCAL based on z of position on calorimeter
		  bool is_bcal = (locNeutralShowers[0]->dDetectorSystem == SYS_BCAL );
		  
		  // Draw on all frames except FCAL frame
		  DGraphicSet gset(kRed, kMarker, 1.25);
		  gset.marker_style = is_bcal ? 22:3;
		  TVector3 tpos( locNeutralShowers[0]->dSpacetimeVertex.X(), 
				 locNeutralShowers[0]->dSpacetimeVertex.Y(), 
				 locNeutralShowers[0]->dSpacetimeVertex.Z());
		  gset.points.push_back(tpos);
		  graphics.push_back(gset);
		  
		  // For BCAL hits, don't draw them on FCAL pane
		  if(is_bcal)continue;
		  
		  // Draw on FCAL pane
		  double dist2 = 2.0 + 2.0*locNeutralShowers[0]->dEnergy;
		  TEllipse *e = new TEllipse(tpos.X(), tpos.Y(), dist2, dist2);
		  e->SetLineColor(gset.color);
		  e->SetFillStyle(0);
		  e->SetLineWidth(1);
		  TMarker *m = new TMarker(tpos.X(), tpos.Y(), gset.marker_style);
		  m->SetMarkerColor(gset.color);
		  m->SetMarkerSize(1.75);
		  graphics_xyB.push_back(e);
		  graphics_xyB.push_back(m);
		}
	}
	
	// FCAL and BCAL thrown photon projections
	if(hdvmf->GetCheckButton("thrown_photons_fcal") || hdvmf->GetCheckButton("thrown_photons_bcal")){
		vector<const DMCThrown*> throwns;
		loop->Get(throwns);
		DGraphicSet gset(kSpring, kMarker, 1.25);
		for(unsigned int i=0; i<throwns.size(); i++){
			const DMCThrown *thrown = throwns[i];
			if(thrown->charge()!=0.0)continue;
			
			// This may seem a little funny, but it saves swimming the reference trajectory
			// multiple times. The GetIntersectionWithCalorimeter() method will find the
			// intersection point of the photon with whichever calorimeter it actually hits
			// or if it doesn't hit either of them. Then, we decide to draw the marker based
			// on whether the flag is set to draw the calorimeter it hit or not.
			DVector3 pos;
			DetectorSystem_t who;
			GetIntersectionWithCalorimeter(thrown, pos, who);
			
			if(who!=SYS_FCAL && who!=SYS_BCAL)continue;
			if(who==SYS_FCAL && !hdvmf->GetCheckButton("thrown_photons_fcal"))continue;
			if(who==SYS_BCAL && !hdvmf->GetCheckButton("thrown_photons_bcal"))continue;
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			
			// Only draw on FCAL pane if photon hits FCAL
			if(who==SYS_BCAL)continue;
			
			double dist2 = 2.0 + 2.0*thrown->energy();
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(kSpring);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
		graphics.push_back(gset);
	}
	
	// FCAL and BCAL thrown charged particle projections
	if(hdvmf->GetCheckButton("thrown_charged_fcal") || hdvmf->GetCheckButton("thrown_charged_bcal")){
		vector<const DMCThrown*> throwns;
		loop->Get(throwns);
		
		for(unsigned int i=0; i<throwns.size(); i++){
			const DMCThrown *thrown = throwns[i];
			if(thrown->charge()==0.0)continue;
			
			// This may seem a little funny, but it saves swimming the reference trajectory
			// multiple times. The GetIntersectionWithCalorimeter() method will find the
			// intersection point of the photon with whichever calorimeter it actually hits
			// or if it doesn't hit either of them. Then, we decide to draw the marker based
			// on whether the flag is set to draw the calorimeter it hit or not.
			DVector3 pos;
			DetectorSystem_t who;
			GetIntersectionWithCalorimeter(thrown, pos, who);
			
			if(who!=SYS_FCAL && who!=SYS_BCAL)continue;
			if(who==SYS_FCAL && !hdvmf->GetCheckButton("thrown_charged_fcal"))continue;
			if(who==SYS_BCAL && !hdvmf->GetCheckButton("thrown_charged_bcal"))continue;
			
			DGraphicSet gset(thrown->charge()>0.0 ? kBlue:kRed, kMarker, 1.25);
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			graphics.push_back(gset);
			
			//double dist2 = 6.0 + 2.0*thrown->momentum().Mag();
			double dist2 = DELTA_R_FCAL;
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(thrown->charge()>0.0 ? kBlue:kRed);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
	}
	
	// FCAL and BCAL reconstructed charged particle projections
	if(hdvmf->GetCheckButton("recon_charged_bcal") || hdvmf->GetCheckButton("recon_charged_fcal")){
	
		// Here we used the full time-based track list, even though it includes multiple
		// hypotheses for each candidate. This is because currently, the photon/track
		// matching code in PID/DPhoton_factory.cc uses the DTrackTimeBased objects and
		// the current purpose of drawing these is to see matching of reconstructed
		// charged tracks with calorimeter clusters.
		vector<const DTrackTimeBased*> tracks;
		loop->Get(tracks, hdvmf->GetFactoryTag("DTrackTimeBased"));
		
		for(unsigned int i=0; i<tracks.size(); i++){
			const DTrackTimeBased *track = tracks[i];
			
			// See notes for above sections
			DVector3 pos;
			DetectorSystem_t who;
			GetIntersectionWithCalorimeter(track, pos, who);
			
			if(who!=SYS_FCAL && who!=SYS_BCAL)continue;
			if(who==SYS_FCAL && !hdvmf->GetCheckButton("recon_charged_fcal"))continue;
			if(who==SYS_BCAL && !hdvmf->GetCheckButton("recon_charged_bcal"))continue;
			
			DGraphicSet gset(track->charge()>0.0 ? kBlue:kRed, kMarker, 1.25);
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			graphics.push_back(gset);
			
			if(who==SYS_BCAL)continue; // Don't draw tracks hitting BCAL on FCAL pane
			
			//double dist2 = 6.0 + 2.0*track->momentum().Mag();
			double dist2 = DELTA_R_FCAL;
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(track->charge()>0.0 ? kBlue:kRed);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
	}

	// DMCTrajectoryPoints
	if(hdvmf->GetCheckButton("trajectories")){
		vector<const DMCTrajectoryPoint*> mctrajectorypoints;
		loop->Get(mctrajectorypoints);
		//sort(mctrajectorypoints.begin(), mctrajectorypoints.end(), DMCTrajectoryPoint_track_cmp);
		
		poly_type drawtype = hdvmf->GetCheckButton("trajectories_lines") ? kLine:kMarker;
		double drawsize  = hdvmf->GetCheckButton("trajectories_lines") ? 1.0:0.3;
		DGraphicSet gset(kBlack, drawtype, drawsize);
		//gset.marker_style = 7;
		TVector3 last_point;
		int last_track=-1;
		int last_part=-1;
		for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
			const DMCTrajectoryPoint *pt = mctrajectorypoints[i];
			
			switch(pt->part){
				case	Gamma:
					if(!hdvmf->GetCheckButton("trajectories_photon"))continue;
					break;
				case	Electron:
					if(!hdvmf->GetCheckButton("trajectories_electron"))continue;
					break;
				case	Positron:
					if(!hdvmf->GetCheckButton("trajectories_positron"))continue;
					break;
				case	Proton:
					if(!hdvmf->GetCheckButton("trajectories_proton"))continue;
					break;
				case	Neutron:
					if(!hdvmf->GetCheckButton("trajectories_neutron"))continue;
					break;
				case	PiPlus:
					if(!hdvmf->GetCheckButton("trajectories_piplus"))continue;
					break;
				case	PiMinus:
					if(!hdvmf->GetCheckButton("trajectories_piminus"))continue;
					break;
				default:
					if(!hdvmf->GetCheckButton("trajectories_other"))continue;
					break;
			}
						
			TVector3 v(pt->x, pt->y, pt->z);

			if(i>0){
				//if((v-last_point).Mag() > 10.0){
				if(pt->track!=last_track || pt->part!=last_part){
					if(hdvmf->GetCheckButton("trajectories_colors")){
			 			switch(last_part){
							case	Gamma:
				 				gset.color = kOrange;
					 			break;
				 			case	Electron:
				 			case	PiMinus:
				 				gset.color = kRed+2;
					 			break;
				 			case	Positron:
				 			case	Proton:
				 			case	PiPlus:
				 				gset.color = kBlue+1;
					 			break;
				 			case	Neutron:
				 				gset.color = kGreen+2;
					 			break;
							default:
				 				gset.color = kBlack;
					 			break;
			 			}
					}else{
						gset.color = kBlack;
					}
					graphics.push_back(gset);
					gset.points.clear();
				}
			}
			
			gset.points.push_back(v);
			last_point = v;
			last_track = pt->track;
			last_part = pt->part;
		}
		
		if(hdvmf->GetCheckButton("trajectories_colors")){
			switch(last_part){
				case	Gamma:
					gset.color = kOrange;
					 break;
				case	Electron:
				case	PiMinus:
					gset.color = kRed+2;
					break;
				case	Positron:
				case	Proton:
				case	PiPlus:
				 	gset.color = kBlue+1;
					break;
				case	Neutron:
				 	gset.color = kGreen+2;
					break;
				default:
				 	gset.color = kBlack;
					break;
			 }
		}else{
			gset.color = kBlack;
		}
		graphics.push_back(gset);
	}

	// DTrackCandidate
	if(hdvmf->GetCheckButton("candidates")){
		vector<const DTrackCandidate*> trackcandidates;
		loop->Get(trackcandidates, hdvmf->GetFactoryTag("DTrackCandidate"));
		for(unsigned int i=0; i<trackcandidates.size(); i++){
			int color=i+1;
			double size=2.0;
			//if(trackcandidates[i]->charge() >0.0) color += 100; // lighter shade
			//if(trackcandidates[i]->charge() <0.0) color += 150; // darker shade
			AddKinematicDataTrack(trackcandidates[i], color, size);
		}
	}

	// DTrackWireBased
	if(hdvmf->GetCheckButton("wiretracks")){
		vector<const DTrackWireBased*> wiretracks;
		loop->Get(wiretracks, hdvmf->GetFactoryTag("DTrackWireBased"));
		for(unsigned int i=0; i<wiretracks.size(); i++){
			AddKinematicDataTrack(wiretracks[i], (wiretracks[i]->charge()>0.0 ? kBlue:kRed)+2, 1.25);
		}
	}

	// DTrackTimeBased
	if(hdvmf->GetCheckButton("timetracks")){
		vector<const DTrackTimeBased*> timetracks;
		loop->Get(timetracks, hdvmf->GetFactoryTag("DTrackTimeBased"));
		for(unsigned int i=0; i<timetracks.size(); i++){
			AddKinematicDataTrack(timetracks[i], (timetracks[i]->charge()>0.0 ? kBlue:kRed)+0, 1.00);
		}
	}

	// DChargedTrack
	if(hdvmf->GetCheckButton("chargedtracks")){
	  vector<const DChargedTrack*> chargedtracks;
		loop->Get(chargedtracks, hdvmf->GetFactoryTag("DChargedTrack"));
		for(unsigned int i=0; i<chargedtracks.size(); i++){
		  int color=kViolet-3;
		  double size=1.25;
		  if (chargedtracks[i]->Get_Charge() > 0) color=kMagenta;
		
		  if (chargedtracks[i]->Get_BestFOM()->mass() > 0.9) size=2.5;
		  AddKinematicDataTrack(chargedtracks[i]->Get_BestFOM(),color,size);
		}
	}

	// DNeutralParticles
	if(hdvmf->GetCheckButton("neutrals")){
		vector<const DNeutralParticle*> photons;
		loop->Get(photons, hdvmf->GetFactoryTag("DNeutralParticle"));
    
		for(unsigned int i=0; i<photons.size(); i++){
		  int color = kBlack;
		  vector<const DNeutralShower*> locNeutralShowers;
		  photons[i]->GetT(locNeutralShowers);
		  DetectorSystem_t locDetSys = locNeutralShowers[0]->dDetectorSystem;
		  if(locDetSys==SYS_FCAL)color = kOrange;
		  if(locDetSys==SYS_BCAL)color = kYellow+2;
		  //if(locDetSys==DPhoton::kCharge)color = kRed;
		  AddKinematicDataTrack(photons[i]->Get_BestFOM(), color, 1.00);
		}
	}
}
void MyProcessor::UpdateBcalDisp(void)
{
  BCALHitCanvas = hdvmf->GetBcalDispFrame();
  BCALHitMatrixU = new TH2F("BCALHitMatrixU","BCAL Hits Upstream Energy;Sector number;Layer;Energy  (MeV)",  48*4+2, -1.5, 192.5, 10, 0., 10.);
  BCALHitMatrixD = new TH2F("BCALHitMatrixD","BCAL Hits Downstream Energy;Sector number;Layer;Energy  (MeV)",48*4+2, -1.5, 192.5, 10, 0., 10.);
  BCALParticles = new TH2F("BCALParticles","BCAL Particle Hits Type;Phi angle [deg];;Particle Momentum",(48*4+2)*4, -1.87, 361.87, 1, 0., 1.);
  BCALHitMatrixU->SetStats(0);
  BCALHitMatrixD->SetStats(0);
  BCALParticles->SetStats(0);
  Float_t size = 0.06;
  BCALHitMatrixU->GetXaxis()->SetTitleSize(size);
  BCALHitMatrixU->GetXaxis()->SetTitleOffset(0.8);
  BCALHitMatrixU->GetXaxis()->SetLabelSize(size);
  BCALHitMatrixU->GetYaxis()->SetTitleSize(size);
  BCALHitMatrixU->GetYaxis()->SetTitleOffset(0.4);
  BCALHitMatrixU->GetYaxis()->SetLabelSize(size);
  BCALHitMatrixU->GetZaxis()->SetTitleSize(size);
  BCALHitMatrixU->GetZaxis()->SetTitleOffset(0.4);
  BCALHitMatrixU->GetZaxis()->SetLabelSize(size);
  BCALHitMatrixD->GetXaxis()->SetTitleSize(size);
  BCALHitMatrixD->GetXaxis()->SetTitleOffset(0.8);
  BCALHitMatrixD->GetXaxis()->SetLabelSize(size);
  BCALHitMatrixD->GetYaxis()->SetTitleSize(size);
  BCALHitMatrixD->GetYaxis()->SetTitleOffset(0.4);
  BCALHitMatrixD->GetYaxis()->SetLabelSize(size);
  BCALHitMatrixD->GetZaxis()->SetTitleSize(size);
  BCALHitMatrixD->GetZaxis()->SetTitleOffset(0.4);
  BCALHitMatrixD->GetZaxis()->SetLabelSize(size);
  BCALParticles->GetXaxis()->SetTitleSize(size);
  BCALParticles->GetXaxis()->SetTitleOffset(0.8);
  BCALParticles->GetXaxis()->SetLabelSize(size);
  BCALParticles->GetYaxis()->SetTitleSize(size);
  BCALParticles->GetYaxis()->SetTitleOffset(0.4);
  BCALParticles->GetYaxis()->SetLabelSize(size);
  BCALParticles->GetZaxis()->SetTitleSize(size);
  BCALParticles->GetZaxis()->SetTitleOffset(0.4);
  BCALParticles->GetZaxis()->SetLabelSize(size);

  if (BCALHitCanvas) {
    vector<const DBCALHit*> locBcalHits;
    loop->Get(locBcalHits);
    BCALHitMatrixU->Reset();
    BCALHitMatrixD->Reset();
    for (unsigned int k=0;k<locBcalHits.size();k++){
      
      const DBCALHit* hit = locBcalHits[k];
      float idxY = (float)hit->layer-1;
      float idxX = (float) (hit->sector-1 + (hit->module-1)*4);
      if (hit->end == DBCALGeometry::kUpstream){
	if (hit->layer==1){
	  BCALHitMatrixU->Fill(idxX,idxY,hit->E);
	} else if (hit->layer==2){
	  BCALHitMatrixU->Fill(idxX,idxY,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+1.,hit->E);	      
	} else if (hit->layer==3){
	  BCALHitMatrixU->Fill(idxX,idxY+1,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+2.,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+3.,hit->E);	      
	} else if (hit->layer==4){
	  BCALHitMatrixU->Fill(idxX,idxY+3,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+4.,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+5.,hit->E);	      
	  BCALHitMatrixU->Fill(idxX,idxY+6.,hit->E);	      
	}
      } else {
	if (hit->layer==1){
	  BCALHitMatrixD->Fill(idxX,idxY,hit->E);
	} else if (hit->layer==2){
	  BCALHitMatrixD->Fill(idxX,idxY,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+1.,hit->E);	      
	} else if (hit->layer==3){
	  BCALHitMatrixD->Fill(idxX,idxY+1,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+2.,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+3.,hit->E);	      
	} else if (hit->layer==4){
	  BCALHitMatrixD->Fill(idxX,idxY+3,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+4.,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+5.,hit->E);	      
	  BCALHitMatrixD->Fill(idxX,idxY+6.,hit->E);	      
	}
      }
    }

    vector<const DBCALIncidentParticle*> locBcalParticles;
    loop->Get(locBcalParticles);
    BCALParticles->Reset();
    BCALPLables.clear();
    for (unsigned int k=0;k<locBcalParticles.size();k++){
      
      const DBCALIncidentParticle* part = locBcalParticles[k];

      float p = TMath::Sqrt(part->px*part->px +  part->py*part->py +  part->pz*part->pz);

	    float phi=999;
	    if (part->x!=0){
	      phi = TMath::ATan(TMath::Abs(part->y/part->x));
	      //cout<<phi<<"   "<<part->y<<" / "<< part->x<<endl;
	      if (part->y>0){
		if (part->x<0.){
		  phi = 3.1415926 - phi;
		}
	      } else {
		if (part->x<0){
		  phi += 3.1415926;
		} else {
		  phi = 3.1415926*2. - phi;
		}
	      }
	      
	      phi = phi*180./3.1415926;
	    }
      BCALParticles->Fill(phi,0.5,p);
      char l[20];
      sprintf(l,"%d",part->ptype);
      TText *t = new TText(phi,1.01,l);
      t->SetTextSize(0.08);
      t->SetTextFont(72);
      t->SetTextAlign(21);
      BCALPLables.push_back(t);
    }
    
    
    BCALHitCanvas->Clear();
    BCALHitCanvas->Divide(1,3);
    BCALHitCanvas->cd(1);
    BCALHitMatrixU->Draw("colz");
    BCALHitCanvas->cd(2);
    BCALHitMatrixD->Draw("colz");
    BCALHitCanvas->cd(3);
    BCALParticles->Draw("colz");
    for (unsigned int n=0;n<BCALPLables.size();n++){
      BCALPLables[n]->Draw("same");
    }
    BCALHitCanvas->Update();
  }
  
}
//------------------------------------------------------------------
// UpdateTrackLabels 
//------------------------------------------------------------------
void MyProcessor::UpdateTrackLabels(void)
{
	// Get the label pointers
	string name, tag;
	map<string, vector<TGLabel*> > &thrownlabs = hdvmf->GetThrownLabels();
	map<string, vector<TGLabel*> > &reconlabs = hdvmf->GetReconstructedLabels();
	hdvmf->GetReconFactory(name, tag);
	
	// Get Thrown particles
	vector<const DMCThrown*> throwns;
	if(loop)loop->Get(throwns);

	// Get the track info as DKinematicData objects
	vector<const DKinematicData*> trks;
	vector<const DKinematicData*> TrksCand;
	vector<const DTrackWireBased*> TrksWireBased;
	vector<const DTrackTimeBased*> TrksTimeBased;
	vector<const DTrackCandidate*> cand;
	if(loop)loop->Get(cand);
	for(unsigned int i=0; i<cand.size(); i++)TrksCand.push_back(cand[i]);

	if(loop)loop->Get(TrksWireBased);
	if(loop)loop->Get(TrksTimeBased);
	
	if(name=="DChargedTrack"){
		vector<const DChargedTrack*> chargedtracks;
		if(loop)loop->Get(chargedtracks, tag.c_str());
		for(unsigned int i=0; i<chargedtracks.size(); i++){
		  trks.push_back(chargedtracks[i]->Get_BestFOM());
		}
	}	
	if(name=="DTrackTimeBased"){
		vector<const DTrackTimeBased*> timetracks;
		if(loop)loop->Get(timetracks, tag.c_str());
		for(unsigned int i=0; i<timetracks.size(); i++)trks.push_back(timetracks[i]);
	}
	if(name=="DTrackWireBased"){
		vector<const DTrackWireBased*> wiretracks;
		if(loop)loop->Get(wiretracks, tag.c_str());
		for(unsigned int i=0; i<wiretracks.size(); i++)trks.push_back(wiretracks[i]);
	}
	if(name=="DTrackCandidate"){
		vector<const DTrackCandidate*> candidates;
		if(loop)loop->Get(candidates, tag.c_str());
		for(unsigned int i=0; i<candidates.size(); i++)trks.push_back(candidates[i]);
	}
	if(name=="DNeutralParticle"){
		vector<const DNeutralParticle*> photons;
		if(loop)loop->Get(photons, tag.c_str());
		for(unsigned int i=0; i<photons.size(); i++) {
		  trks.push_back(photons[i]->Get_BestFOM());
		}
	}
	
	// Clear all labels (i.e. draw ---- in them)
	map<string, vector<TGLabel*> >::iterator iter;
	for(iter=reconlabs.begin(); iter!=reconlabs.end(); iter++){
		vector<TGLabel*> &labs = iter->second;
		for(unsigned int i=1; i<labs.size(); i++){
			labs[i]->SetText("--------");
		}
	}
	for(iter=thrownlabs.begin(); iter!=thrownlabs.end(); iter++){
		vector<TGLabel*> &labs = iter->second;
		for(unsigned int i=1; i<labs.size(); i++){
			labs[i]->SetText("--------");
		}
	}

	// Loop over thrown particles and fill in labels
	int ii=0;
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *trk = throwns[i];
		if(trk->type==0)continue;
		int row = thrownlabs["trk"].size()-(ii++)-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z;
		trkno<<setprecision(4)<<i+1;
		thrownlabs["trk"][row]->SetText(trkno.str().c_str());

		thrownlabs["type"][row]->SetText(ParticleType((Particle_t)trk->type));

		p<<setprecision(4)<<trk->momentum().Mag();
		thrownlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		thrownlabs["theta"][row]->SetText(theta.str().c_str());

		double myphi = trk->momentum().Phi();
		if(myphi<0.0)myphi+=2.0*M_PI;
		phi<<setprecision(4)<<myphi;
		thrownlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		thrownlabs["z"][row]->SetText(z.str().c_str());
	}

	// Loop over tracks and fill in labels
	for(unsigned int i=0; i<trks.size(); i++){
		const DKinematicData *trk = trks[i];
		int row = reconlabs["trk"].size()-i-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z, chisq_per_dof, Ndof,cand;
		stringstream fom;
		trkno<<setprecision(4)<<i+1;
		reconlabs["trk"][row]->SetText(trkno.str().c_str());
		
		double mass = trk->mass();
		if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
		else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
		else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
		else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
		else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
		else type<<"q=";
		if (fabs(trk->charge())>1.e-4){
		  type<<(trk->charge()>0 ? "+":"-");
		}
		reconlabs["type"][row]->SetText(type.str().c_str());

		p<<setprecision(4)<<trk->momentum().Mag();
		reconlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		reconlabs["theta"][row]->SetText(theta.str().c_str());

		phi<<setprecision(4)<<trk->momentum().Phi()*TMath::RadToDeg();
		reconlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		reconlabs["z"][row]->SetText(z.str().c_str());

		// Get chisq and Ndof for DTrackTimeBased or DTrackWireBased objects
		const DTrackTimeBased *timetrack=dynamic_cast<const DTrackTimeBased*>(trk);
		const DTrackWireBased *track=dynamic_cast<const DTrackWireBased*>(trk);	
	
		const DTrackCandidate *candidate=dynamic_cast<const DTrackCandidate*>(trk);
		const DChargedTrackHypothesis *chargedtrack=dynamic_cast<const DChargedTrackHypothesis *>(trk);

		if(timetrack){
			chisq_per_dof<<setprecision(4)<<timetrack->chisq/timetrack->Ndof;
			Ndof<<timetrack->Ndof;
			fom << timetrack->FOM;
		}else if(track){
			chisq_per_dof<<setprecision(4)<<track->chisq/track->Ndof;
			Ndof<<track->Ndof;
			fom << "N/A";
		}else if (candidate){
			chisq_per_dof<<setprecision(4)<<candidate->chisq/candidate->Ndof;
			Ndof<<candidate->Ndof;
			fom << "N/A";
		}
		else if (chargedtrack){
		  chisq_per_dof<<setprecision(4)<<chargedtrack->dChiSq_Track/chargedtrack->dNDF_Track;
			Ndof<<chargedtrack->dNDF_Track;
			fom << chargedtrack->dFOM;
		}
		else{
		  chisq_per_dof << "--------";
		  Ndof << "--------";
		  fom << "--------";
		}
		
		reconlabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
		reconlabs["Ndof"][row]->SetText(Ndof.str().c_str());
		reconlabs["FOM"][row]->SetText(fom.str().c_str());
		
		if (timetrack){
		  cand << timetrack->candidateid;
		}
		else if (track){
		  cand << track->candidateid;
		}
		else {
		  cand << "--------";
		}
		reconlabs["cand"][row]->SetText(cand.str().c_str());
	}
	
	// Have the pop-up window with the full particle list update it's labels
	fulllistmf->UpdateTrackLabels(throwns, trks);
	debugermf->SetTrackCandidates(TrksCand);
	debugermf->SetTrackWireBased(TrksWireBased);
	debugermf->SetTrackTimeBased(TrksTimeBased);
	debugermf->UpdateTrackLabels();
}                 

//------------------------------------------------------------------
// AddKinematicDataTrack 
//------------------------------------------------------------------
void MyProcessor::AddKinematicDataTrack(const DKinematicData* kd, int color, double size)
{
	// Create a reference trajectory with the given kinematic data and swim
	// it through the detector.
	DReferenceTrajectory rt(Bfield);
	rt.Rmax_interior = RMAX_INTERIOR;
	rt.Rmax_exterior = RMAX_EXTERIOR;

	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt.SetDRootGeom(RootGeom);
		rt.SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt.SetDRootGeom(NULL);
		rt.SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}

	rt.SetMass(kd->mass());
	rt.Swim(kd->position(), kd->momentum(), kd->charge());

	// Create a new graphics set and fill it with all of the trajectory points
	DGraphicSet gset(color, kLine, size);
	DReferenceTrajectory::swim_step_t *step = rt.swim_steps;
	for(int i=0; i<rt.Nswim_steps; i++, step++){
	  TVector3 tpoint(step->origin.X(),step->origin.Y(),step->origin.Z());
		gset.points.push_back(tpoint);
	}
	
	// Push the graphics set onto the stack
	graphics.push_back(gset);
}

//------------------------------------------------------------------
// GetIntersectionWithCalorimeter 
//------------------------------------------------------------------
void MyProcessor::GetIntersectionWithCalorimeter(const DKinematicData* kd, DVector3 &pos, DetectorSystem_t &who)
{
	// Create a reference trajectory with the given kinematic data and swim
	// it through the detector.
	DReferenceTrajectory rt(Bfield);
	rt.Rmax_interior = RMAX_INTERIOR;
	rt.Rmax_exterior = RMAX_EXTERIOR;

	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt.SetDRootGeom(RootGeom);
		rt.SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt.SetDRootGeom(NULL);
		rt.SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}

	rt.SetMass(kd->mass());
	rt.Swim(kd->position(), kd->momentum(), kd->charge());

	// Find intersection with FCAL
	DVector3 pos_fcal;
	double s_fcal = 1.0E6;
	DVector3 origin(0.0, 0.0, FCAL_Zmin);
	DVector3 norm(0.0, 0.0, -1.0);
	//rt.GetIntersectionWithPlane(origin, norm, pos_fcal, &s_fcal); // This gives different answer than below!
	DVector3 p_at_intersection;
	rt.GetIntersectionWithPlane(origin, norm, pos_fcal, p_at_intersection, &s_fcal, NULL, NULL, SYS_FCAL);
	if(pos_fcal.Perp()<FCAL_Rmin || pos_fcal.Perp()>FCAL_Rmax || !isfinite(pos_fcal.Z()))s_fcal = 1.0E6;
	
	// Find intersection with BCAL
	DVector3 pos_bcal;
	double s_bcal = 1.0E6;
	rt.GetIntersectionWithRadius(BCAL_Rmin, pos_bcal, &s_bcal);
	if(pos_bcal.Z()<BCAL_Zmin || pos_bcal.Z()>(BCAL_Zmin+BCAL_Zlen) || !isfinite(pos_bcal.Z()))s_bcal = 1.0E6;

	if(s_fcal>10000.0 && s_bcal>10000.0){
		// neither calorimeter hit
		who = SYS_NULL;
		pos.SetXYZ(0.0,0.0,0.0);
	}else if(s_fcal<s_bcal){
		// FCAL hit
		who = SYS_FCAL;
		pos = pos_fcal;
	}else{
		// BCAL hit
		who = SYS_BCAL;
		pos = pos_bcal;
	}
}

//------------------------------------------------------------------
// GetFactoryNames 
//------------------------------------------------------------------
void MyProcessor::GetFactoryNames(vector<string> &facnames)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		vector<string> facnames;
		loops[0]->GetFactoryNames(facnames);
	}
}

//------------------------------------------------------------------
// GetFactories 
//------------------------------------------------------------------
void MyProcessor::GetFactories(vector<JFactory_base*> &factories)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		factories = loops[0]->GetFactories();
	}
}

//------------------------------------------------------------------
// GetNrows 
//------------------------------------------------------------------
unsigned int MyProcessor::GetNrows(const string &factory, string tag)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		// Here is something a little tricky. The GetFactory() method of JEventLoop
		// gets the factory of the specified data name and tag, but without trying
		// to substitute a user-specified tag (a'la -PDEFTAG:XXX=YYY) as is done
		// on normal Get() method calls. Therefore, we have to check for the default
		// tags ourselves and substitute it "by hand".
		if(tag==""){
			map<string,string> default_tags = loops[0]->GetDefaultTags();
			map<string, string>::const_iterator iter = default_tags.find(factory);
			if(iter!=default_tags.end())tag = iter->second.c_str();
		}
		JFactory_base *fac = loops[0]->GetFactory(factory, tag.c_str());

		// Since calling GetNrows will cause the program to quit if there is
		// not a valid event, then first check that there is one before calling it
		if(loops[0]->GetJEvent().GetJEventSource() == NULL)return 0;
		
		return fac==NULL ? 0:(unsigned int)fac->GetNrows();
	}
	
	return 0;
}

//------------------------------------------------------------------
// GetDReferenceTrajectory 
//------------------------------------------------------------------
void MyProcessor::GetDReferenceTrajectory(string dataname, string tag, unsigned int index, DReferenceTrajectory* &rt, vector<const DCDCTrackHit*> &cdchits)
{
_DBG__;
	// initialize rt to NULL in case we don't find the one requested
	rt = NULL;
	cdchits.clear();

	// Get pointer to the JEventLoop so we can get at the data
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()==0)return;
	JEventLoop* &loop = loops[0];
	
	// Variables to hold track parameters
	DVector3 pos, mom(0,0,0);
	double q=0.0;
	double mass;

	// Find the specified track	
	if(dataname=="DChargedTrack"){
		vector<const DChargedTrack*> chargedtracks;
		vector<const DTrackTimeBased*> timebasedtracks;
		loop->Get(chargedtracks, tag.c_str());
		if(index>=chargedtracks.size())return;
		q = chargedtracks[index]->Get_Charge();
		pos = chargedtracks[index]->Get_BestFOM()->position();
		mom = chargedtracks[index]->Get_BestFOM()->momentum();
		chargedtracks[index]->Get_BestFOM()->GetT(timebasedtracks);
		timebasedtracks[0]->Get(cdchits);
		mass = chargedtracks[index]->Get_BestFOM()->mass();
	}

	if(dataname=="DTrackTimeBased"){
		vector<const DTrackTimeBased*> timetracks;
		loop->Get(timetracks, tag.c_str());
		if(index>=timetracks.size())return;
		q = timetracks[index]->charge();
		pos = timetracks[index]->position();
		mom = timetracks[index]->momentum();
		timetracks[index]->Get(cdchits);
		mass = timetracks[index]->mass();
	}

	if(dataname=="DTrackWireBased"){
		vector<const DTrackWireBased*> wiretracks;
		loop->Get(wiretracks, tag.c_str());
		if(index>=wiretracks.size())return;
		q = wiretracks[index]->charge();
		pos = wiretracks[index]->position();
		mom = wiretracks[index]->momentum();
		wiretracks[index]->Get(cdchits);
		mass = wiretracks[index]->mass();
	}

	if(dataname=="DTrackCandidate"){
		vector<const DTrackCandidate*> tracks;
		loop->Get(tracks, tag.c_str());
		if(index>=tracks.size())return;
		q = tracks[index]->charge();
		pos = tracks[index]->position();
		mom = tracks[index]->momentum();
		tracks[index]->Get(cdchits);
		mass = tracks[index]->mass();
	}

	if(dataname=="DMCThrown"){
		vector<const DMCThrown*> tracks;
		loop->Get(tracks, tag.c_str());
		if(index>=tracks.size())return;
		const DMCThrown *t = tracks[index];
		q = t->charge();
		pos = t->position();
		mom = t->momentum();
		tracks[index]->Get(cdchits);
		mass = tracks[index]->mass();
_DBG_<<"mass="<<mass<<endl;
	}

	// Make sure we found a charged particle we can track
	if(q==0.0 || mom.Mag()<0.01)return;
	
	// Create a new DReference trajectory object. The caller takes
	// ownership of this and so they are responsible for deleting it.
	rt = new DReferenceTrajectory(Bfield);
	rt->Rmax_interior = RMAX_INTERIOR;
	rt->Rmax_exterior = RMAX_EXTERIOR;
	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt->SetDRootGeom(RootGeom);
		rt->SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt->SetDRootGeom(NULL);
		rt->SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}
	rt->Swim(pos, mom, q);
}

//------------------------------------------------------------------
// GetAllWireHits
//------------------------------------------------------------------
void MyProcessor::GetAllWireHits(vector<pair<const DCoordinateSystem*,double> > &allhits)
{
	/// Argument is vector of pairs that contain a pointer to the
	/// DCoordinateSystem representing a wire and a double that
	/// represents the drift distance. To get info on the specific
	/// wire, one needs to attempt a dynamic_cast to both a DCDCWire
	/// and a DFDCWire and access the parameters of whichever one succeeds.
	
	// Get pointer to the JEventLoop so we can get at the data
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()==0)return;
	JEventLoop* &loop = loops[0];

	// Get CDC wire hits
	vector<const DCDCTrackHit*> cdchits;
	loop->Get(cdchits);
	for(unsigned int i=0; i<cdchits.size(); i++){
		pair<const DCoordinateSystem*,double> hit;
		hit.first = cdchits[i]->wire;
		hit.second = cdchits[i]->dist;
		allhits.push_back(hit);
	}
	
	// Get FDC wire hits
	vector<const DFDCPseudo*> fdchits;
	loop->Get(fdchits);
	for(unsigned int i=0; i<fdchits.size(); i++){
		pair<const DCoordinateSystem*,double> hit;
		hit.first = fdchits[i]->wire;
		hit.second = 0.0055*fdchits[i]->time;
		allhits.push_back(hit);
	}
}



