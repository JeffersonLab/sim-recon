// $Id$
//
//    File: JEventProcessor_BCAL_TimeCalibration.cc
// Created: Mon Apr 18 15:28:52 CST 2016
// Creator: semenov (on Linux selene.phys.uregina.ca 2.6.32-573.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_BCAL_TimeCalibration.h"
using namespace jana;

#include <TLorentzVector.h>
#include "TMath.h"

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "BCAL/DBCALGeometry.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DEventRFBunch.h"
#include "PID/DDetectorMatches.h"
#include "DAQ/Df250PulsePedestal.h" // Needed for pulse peak information
#include "DAQ/Df250PulseIntegral.h" // Needed for pulse integral !!! added
#include "BCAL/DBCALDigiHit.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "PID/DVertex.h"
#include "TAGGER/DTAGHHit.h"
#include "TAGGER/DTAGMHit.h"
//#include "TRACKING/DTrackFinder.h"
#include "TRIGGER/DL1Trigger.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TChain.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPaveLabel.h>
#include <TGraphErrors.h>
#include <TLine.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TROOT.h"

#include <mutex>
#include <vector>
#include <deque>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

static mutex mtx;

static int ievent;

static int iRunNumber;
 
static double R0[5]={0., 65.089, 67.1464, 71.2612, 77.4334}; //layer-dependent radius of the beginning cell position
static double DD[5]={0., 2.0574, 4.1148, 6.1722, 8.23}; //layer-dependent cell thickness 

static int busyShower[1000];

static double cvel = 29.9792458 ; //Speed-of-Light (cm/ns)

static unsigned short int aU[48][4][3][5000],aD[48][4][3][5000],tU[48][4][3][5000],tD[48][4][3][5000];
static unsigned short int Z0[48][4][3][5000], zH[48][4][3][5000], RUD[48][4][3][5000];
static float Ene[48][4][3][5000];
		
static double r2U,r2D,rUD,zUD, cosPHI,sinPHI,PHIp,PHIs,PHIdiff;
static unsigned short int ipo[48][4][3];
static double itU, itD;
static double Ze[2]={17.,407.}; //CORRECT Z-position of BCAL ends (cm) May 19, 2016
static int imo,ise,ila,ip;
static int pdfps = 0; //PDF(0) or PostScript(1) picture
	
static TGraphErrors * gr01;
static TGraphErrors * gr11;
static TGraphErrors * gr02;
static TGraphErrors * gr12;

static TGraph * gr1;
static TGraph * gr2;
static TGraph * gr3;

static TH1F * h1;
static TH1F * h2;
static TH1F * h3;
static TH1F * h4;
static TH1F * his4;
static TH1F * his44;

static TF1 * fgg10;
static TF1 * fgg20;
static TF1 * fgau;
static TF1 * fgau2;
static TF1 * fz1;
static TF1 * fz2;

static char * ftxt;
static char * ftxt2;
static char * ftxt3;
static char * ftxt4;
static char * hisname;
static char * hisname2;

static TLatex * txt1;
static TLatex * txt2;
static TLatex * txt3;
static TLatex * txt4;

static TLine * lin1;

static         double tADC_U, tADC_D;
static 	double E1[100], E3[100], E13[100];
	
static 	double dZs1[3]={-11., -9.5, -23.};
static 	double dZs2[3]={ 9.0, 10.5,  37.};

//*********************************************************************************
static Double_t tzfunc (Double_t *x, Double_t *par) {
   Double_t x1=x[0]-212.;
   Double_t func = par[0]+par[1]*x1;
   return func;
}
//*********************************************************************************
static Double_t twfunc2 (Double_t *x, Double_t *par) {
   Double_t t0=par[0];
   Double_t x1=x[0]-t0;

   Double_t func = par[1] + par[2]*pow(x1,par[3]) + par[4]*pow(x1,par[5]);
   return func;
}
//*********************************************************************************
static Double_t gfunc (Double_t *x, Double_t *par) {
   Double_t x1=x[0]-par[1];

   Double_t func;
 
   func = par[0]*exp(-x1*x1/(2.*par[2]*par[2]));
   return func;
}
//*********************************************************************************

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_TimeCalibration());
}
} // "C"

//------------------
// JEventProcessor_BCAL_TimeCalibration (Constructor)
//------------------
JEventProcessor_BCAL_TimeCalibration::JEventProcessor_BCAL_TimeCalibration()
{

}

//------------------
// ~JEventProcessor_BCAL_TimeCalibration (Destructor)
//------------------
JEventProcessor_BCAL_TimeCalibration::~JEventProcessor_BCAL_TimeCalibration()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_TimeCalibration::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

	ievent = 0;
	
	for (imo=0; imo<48; imo++) {
	  for (ise=0; ise<4; ise++) {
	    for (ila=0; ila<3; ila++) {
	      ipo[imo][ise][ila] = 0;  // Null the counter
	    }
	  }
	}
		
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_TimeCalibration::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	iRunNumber = runnumber;
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_TimeCalibration::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	
static 	int goodFlag = 0;
static 	double ttU, ttD;
static 	int aaU,aaD, hend,hmodule,hsector,hlayer;
static 	double t0,x0,y0,z0, tt0;
static 	double Xs,Ys,Zs, Xh,Yh,Zh, Rh,rU,rD,zU,zD ;
static 	double x,y,z,R, dPhi,dZ, E;

	goodFlag = 0;
	
//****************	
//CHECK TRIGGER TYPE
//       const DL1Trigger* locTrigger = NULL;
//       locEventLoop->GetSingle(locTrigger);

         uint32_t trig_mask=0;
         const DL1Trigger *trig_words = NULL;
	  try {
	    locEventLoop->GetSingle(trig_words);
	  } catch(...) {};
	  if (trig_words) {
	    trig_mask = trig_words->trig_mask;
	  } else {
	    trig_mask=0;
	  }    
	    
	    
	    
//****************
	
	vector<const DTrackTimeBased*> locTrackTimeBased;
	  locEventLoop->Get(locTrackTimeBased);
	vector<const DChargedTrackHypothesis*> locTrackHypothesis;
	  locEventLoop->Get(locTrackHypothesis);
	
	vector<const DBCALShower*> locBCALShowers;
	  locEventLoop->Get(locBCALShowers);
	  	
	vector<const DBCALShower*> matchedShowers;
	
	vector<const DBCALCluster*> locBCALClusters;
//	  locEventLoop->Get(locBCALClusters);
	  
	vector<const DBCALPoint*> matchedBCALPoints;
	  
	vector<const DBCALPoint*> locBCALPoints;
	  locEventLoop->Get(locBCALPoints);
	  
	vector<const DVertex*> kinfitVertex;
	  locEventLoop->Get(kinfitVertex);
	  
	vector<const DEventRFBunch*> RFbunch;
	  locEventLoop->Get(RFbunch);
	int Nvotes;  

	for (unsigned int ij=0; ij<RFbunch.size(); ij++) {
	  Nvotes = RFbunch[ij]->dNumParticleVotes;
	}    
	
	vector<const DBCALUnifiedHit *> locBCALUnifiedHit;
	  locEventLoop->Get(locBCALUnifiedHit);
	  
	vector<const DBCALUnifiedHit *> matchedBCALUnifiedHit;

        const DBCALHit * thisADCHit;
	
        const DBCALTDCHit * thisTDCHit;
	
// Here we want to match tracks to the BCAL so that we can throw away charged particles in the BCAL to decrease combinatorics. To do this we look over time based tracks and swim each track in an event out to each BCAL shower in each event. Then check to see if the BCAL shower and charged track are within some common z and phi region of the calorimeter. If they are then we fill a vector with those showers that are geometrically matched to a charged track. This vector will be called later to skip such showers.

	DVector3 mypos(0.0,0.0,0.0);
	
	// Only one thread at a time from here to end of method
	lock_guard<mutex> lck(mtx);
	
	for (int ij=0; ij<1000; ij++) busyShower[ij]=0;
	
	for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
	  tt0 = locTrackHypothesis[i]->t0(); if (tt0<-10000.||tt0>10000.) goodFlag = 1;
	  for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	
	  	x = locBCALShowers[j]->x;
		y = locBCALShowers[j]->y;
		z = locBCALShowers[j]->z;
		DVector3 pos_bcal(x,y,z);  // Vector that has the BCAL shower centroid position so we know where to try and match charged track
		R = pos_bcal.Perp();
		//locTrackTimeBased[i]->rt->GetIntersectionWithRadius(R, mypos);  // This will swim a track out of the drift chambers to some radius that is defined by pos_bcal and returns a vector mypos that will give the location of the charged track at the redius we input (R)
		locTrackTimeBased[i]->GetProjection(SYS_BCAL, mypos);  // This will swim a track out of the drift chambers to some radius that is defined by pos_bcal and returns a vector mypos that will give the location of the charged track at the redius we input (R)
                                
		dPhi = TMath::Abs(mypos.Phi()-pos_bcal.Phi()); // find how far away in phi and z that the shower and track are
		dZ = TMath::Abs(mypos.Z() - z);
		if(dZ<60.0&&dPhi<0.50&&busyShower[j]==0) { // if the track and shower are close enough in phi and z then fill the matchedshowers vector
		  busyShower[j]=1;         // Mark the shower to be NOT used 
		}
	      	
          }
	}
	
	for (unsigned int j=0; j< locBCALShowers.size(); ++j){
	  if (busyShower[j]==0) matchedShowers.push_back(locBCALShowers[j]); //Showers that are NOT matched with tracks
	}
	
	t0 = 99999999.;						  
	for (unsigned int ij=0; ij<kinfitVertex.size(); ij++) {
	  t0= kinfitVertex[ij]->dSpacetimeVertex.T();
	  x0= kinfitVertex[ij]->dSpacetimeVertex.X();
	  y0= kinfitVertex[ij]->dSpacetimeVertex.Y();
	  z0= kinfitVertex[ij]->dSpacetimeVertex.Z();
	} 
	
//Check on the signal in the first layer and E1/E3 ratio
	for (unsigned int ij=0; ij< matchedShowers.size(); ++ij){ //<------------------------- ij loop on Neutral Showers
	  E1[ij]=0;
	  E3[ij]=0;
	  E13[ij] = -999999.;
	  
	  if (goodFlag==0) {
	    matchedShowers[ij]->Get(matchedBCALPoints);
	    for (unsigned int il=0; il<matchedBCALPoints.size(); il++) { //<------------------------- il loop on Points
	      matchedBCALPoints[il]->Get(matchedBCALUnifiedHit);
	      for (unsigned int ik=0; ik<matchedBCALUnifiedHit.size(); ik++) { //<---------------------- ik loop on Hits
	        matchedBCALUnifiedHit[ik]->GetSingle(thisADCHit);
		hlayer = thisADCHit->layer;
		aaU = thisADCHit->pulse_peak;
	      }//<------------------------------------------------------------------------ End of ik loop on Hits
	      
	      if (hlayer==1) E1[ij] += matchedBCALPoints[il]->E();
	      if (hlayer==3) E3[ij] += matchedBCALPoints[il]->E();
	      
	    }//<------------------------------------------------------------------------ End of il loop on Points
          }
	  if (E3[ij]>0) {
	    E13[ij] = E1[ij]/E3[ij];
	  } else {
	    E13[ij] = 999999.;
	  }     
	} //<----------------------------------------------------------------------------- End of ij loop on Neutral Showers
//End of check on the signal in the first layer	 
	  	  	  
	for (unsigned int ij=0; ij< matchedShowers.size(); ++ij){ //<------------------------- ij loop on Neutral Showers
	  if (goodFlag==0) {
	    Xs  = matchedShowers[ij]->x;
	    Ys  = matchedShowers[ij]->y;
	    Zs  = matchedShowers[ij]->z;
            PHIs = atan2(Ys,Xs);
	  
	    matchedShowers[ij]->Get(matchedBCALPoints);
	    
	    for (unsigned int il=0; il<matchedBCALPoints.size(); il++) { //<------------------------- il loop on Points
	    
	      Rh = matchedBCALPoints[il]->r();
	      PHIp = matchedBCALPoints[il]->phi();
	      cosPHI = cos(PHIp);
	      sinPHI = sin(PHIp);
	    
	      Zh = matchedBCALPoints[il]->z();
	      
	      E = matchedBCALPoints[il]->E();
	      
	      matchedBCALPoints[il]->Get(matchedBCALUnifiedHit);
	  
	      ttU=-9999.;
	      ttD=-9999.;
	      aaU=-9999.;
	      aaD=-9999.;
	      tADC_U = -9999.;
	      tADC_D = -9999.;
	      
	    for (unsigned int ik=0; ik<matchedBCALUnifiedHit.size(); ik++) { //<---------------------- ik loop on Hits
	      matchedBCALUnifiedHit[ik]->GetSingle(thisADCHit);
	      matchedBCALUnifiedHit[ik]->GetSingle(thisTDCHit);
	      
	      hend = thisADCHit->end;
	      hmodule = thisADCHit->module;
	      hsector = thisADCHit->sector;
	      hlayer = thisADCHit->layer;	  
	      
	      if (hend==0) tADC_U = matchedBCALUnifiedHit[ik]->t_ADC;
	      if (hend==1) tADC_D = matchedBCALUnifiedHit[ik]->t_ADC;
	      
	      if (hend==0) {
	        aaU = thisADCHit->pulse_peak;
		rU=0.;
		if (aaU>0) rU = R0[hlayer]+DD[hlayer]*2./aaU; //<---- Threshold of 2 counts
		zU = z0 + (Zh+65.-z0) * rU/Rh;
	        Xh = rU*cosPHI;
	        Yh = rU*sinPHI;
		r2U = (Xh-x0)*(Xh-x0)+(Yh-y0)*(Yh-y0); 
	      }	
		
	      if (hend==1) {
	        aaD = thisADCHit->pulse_peak;
		rD=0.;
		if (aaD>0) rD = R0[hlayer]+DD[hlayer]*2./aaD; //<---- Threshold of 2 counts
		zD = z0 + (Zh+65.-z0) * rD/Rh;
	        Xh = rD*cosPHI;
	        Yh = rD*sinPHI;
		r2D = (Xh-x0)*(Xh-x0)+(Yh-y0)*(Yh-y0); 
	      }	
	      
	      if (matchedBCALUnifiedHit[ik]->has_TDC_hit && Nvotes>=2) {
	        if (hend==0) ttU = thisTDCHit->t;
	        if (hend==1) ttD = thisTDCHit->t;
	      }	
	    } //end of ik loop
	    
	      imo = hmodule-1;
	      ise = hsector-1;
	      ila = hlayer -1;
	      ip  = ipo[imo][ise][ila];
	      
	      PHIdiff = abs(PHIp-PHIs);
	      if (PHIdiff>3.1415927) PHIdiff = 6.2831854 - PHIdiff;
	      itU = 100.*(ttU - tt0);
	      itD = 100.*(ttD - tt0);
	      zUD = 0.5*(zU+zD);
	      rUD = 0.5*(sqrt(r2U)+sqrt(r2D));

//japp->RootFillLock(this);         
	      if (ip<5000                                                      // To not overflow array
	          && itU>0. && itD>0. && itU<60000. && ttD<60000. 
		  && aaU>5. && aaD>5. && aaU<60000. && aaD<60000.              // Both ends made signals
	          && matchedBCALPoints.size()>2 && matchedBCALPoints.size()<30 // Reasonable multiplicity
		  && hlayer>0 && hlayer<4                                      // Layers 1-3
		  && trig_mask==4                                              // BCAL trigger
		  && zUD>(Ze[0]+10.) && zUD<(Ze[1]-10.)                        // Z inside BCAL
		  && rUD>60. && rUD<100.
		  && E13[ij]>0.1 && E13[ij]<10000.                             // Suppressed neutrons as on 170706 
		  && abs(z0-65.)<15.                                           // Vertex is inside the target
                  && PHIdiff<0.035
		  && (zUD-Zs)>dZs1[ila] && (zUD-Zs)<dZs2[ila]) {               // Hit is not far away from the Shower
	      
		aU[imo][ise][ila][ip] = 1*aaU;
		aD[imo][ise][ila][ip] = 1*aaD;
		tU[imo][ise][ila][ip] = 1*itU;
		tD[imo][ise][ila][ip] = 1*itD;
		Z0[imo][ise][ila][ip] = 100*z0;
		zH[imo][ise][ila][ip] = 100*zUD;
		RUD[imo][ise][ila][ip] = 100*rUD;
		Ene[imo][ise][ila][ip] = E;
		
		ipo[imo][ise][ila]++;
	        ievent++;
	      }
//japp->RootFillUnLock(this);
	    
	    } //end of il loop
	    
	  } //if(GoodFlag end 
	} //end of ij loop
        
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_TimeCalibration::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
//*********************************************************************************
jerror_t JEventProcessor_BCAL_TimeCalibration::fini(void)
{
	// Called before program exit after event processing is finished.
	
//*********************************************************************************
static double f1p0,f1p1,f1p2,f1p3,f1p4,f1p5;
static double ff1p0,ff1p1,ff1p2,ff1p3,ff1p4,ff1p5;

static double f2p0,f2p1,f2p2,f2p3,f2p4,f2p5;
static double ff2p0,ff2p1,ff2p2,ff2p3,ff2p4,ff2p5;

static double f3p0,f3p1,f3p2,f3p3,f3p4,f3p5;
static double ff3p0,ff3p1,ff3p2,ff3p3,ff3p4,ff3p5;

static double f4p0,f4p1,f4p2,f4p3,f4p4,f4p5;
static double ff4p0,ff4p1,ff4p2,ff4p3,ff4p4,ff4p5;

static double xp10,xp11,xp12, xp20,xp21,xp22;
static double xp1,xp2,xp3,xp4,xp5, ytmp, ytmp2, Zh, pL;
static double veff2, zp0,zp1, veff22, zzp0,zzp1, vef3,vef33;

int ipo3;
//************************************************************************************
//double cvel = 29.9792458 ; //Speed-of-Light (cm/ns)
//double veff = 16.5;       //Init veff (cm/ns)
static double vef1 = 15.9;       //Init veff (cm/ns)
static double vef2 = 16.9;       //Init veff (cm/ns)
//double Ze[2]={17.,407.}; //CORRECT Z-position of BCAL ends (cm) May 19, 2016

static double y1[5000],x1[5000],ey1[5000],ex1[5000],z1[5000],y10[5000],y20[5000];
static double yy1[5000],yy2[5000],xx1[5000],eyy1[5000],exx1[5000];
static double y2[5000],x2[5000],ey2[5000],ex2[5000];
static int ipo2 = 0;

static double xaU,xaD,xtU,xtD,xZ0,xZH,xRUD;

static double gme,gsi,egme,egsi, gme2,gsi2,egme2,egsi2,  cont1,cont4;

static float itt1, itt2, itt3, itt4, xkoe;
//************************************************************************************

static   char *outcalname = new char[90];
    sprintf(outcalname, "TDCcalib_run%d.vec",iRunNumber);
    ofstream dat100(outcalname);

static   char *outcalname2 = new char[90];
    sprintf(outcalname2, "ccdb-BCAL_effective_velocities_regina-run%d.dat",iRunNumber);
    ofstream datveff(outcalname2);

static   char *outcalname3 = new char[90];
    sprintf(outcalname3, "ccdb-BCAL_timewalk_constants_regina-run%d.dat",iRunNumber);
    ofstream datwalk(outcalname3);

static   char *outcalname4 = new char[90];
    sprintf(outcalname4, "WARNINGS-run%d.vec",iRunNumber);
    ofstream datwarn(outcalname4);

//************************************************************************************
//************************************************************************************

//Canvas Definition  
static   TCanvas *c1 = new TCanvas("c1", " c1",0,0,699,499);

//************************************************************************************

static   char *pdfname = new char[90];
static   char *pdfname1 = new char[90];
static   char *pdfname2 = new char[90];
  if (pdfps==0) {
    sprintf(pdfname, "TDCcalib_run%d.pdf",iRunNumber);
    sprintf(pdfname1, "TDCcalib_run%d.pdf[",iRunNumber);
    sprintf(pdfname2, "TDCcalib_run%d.pdf]",iRunNumber);
  } else {  
    sprintf(pdfname, "TDCcalib_run%d.ps",iRunNumber);
    sprintf(pdfname1, "TDCcalib_run%d.ps[",iRunNumber);
    sprintf(pdfname2, "TDCcalib_run%d.ps]",iRunNumber);
  }
  c1->Update(); c1->Print(pdfname1);  
//***************************************************    
 
for (int imo=0; imo<48; imo++) {  //Start of imo loop   
for (int ila=0; ila<3; ila++) {  //Start of ila loop   
for (int ise=0; ise<4; ise++) {  //Start of ise loop
//************************************************************************************

h1 = new TH1F("h1","Middle U",150,-50.,150.);
h2 = new TH1F("h2","High U",  150,-50.,150.);

h3 = new TH1F("h3","Middle D",150,-50.,150.);
h4 = new TH1F("h4","High D",  150,-50.,150.);

  hisname = new char[90];
    sprintf(hisname, "#gamma Hits: Run %d, Module=%02i, Sector=%i, Layer=%i, End=0",iRunNumber,imo+1,ise+1,ila+1);
  hisname2 = new char[90];
    sprintf(hisname2, "#gamma Hits: Run %d, Module=%02i, Sector=%i, Layer=%i, End=1",iRunNumber,imo+1,ise+1,ila+1);
    
//*******
  his4 = new TH1F("his4"," ",200,-8.,12.);
   his4->GetXaxis()->SetTitle("Corrected #DeltaTime (ns)");
   his4->GetYaxis()->SetTitle("Counts");
   his4->SetTitle(hisname);
//   his4->SetMaximum(20.);
   his4->SetLineColor(1);
   his4->SetLineWidth(1.9);
   his4->SetFillStyle(1001);
   his4->SetFillColor(kGreen-7);
   
//*******
  his44 = new TH1F("his44"," ",200,-8.,12.);
   his44->GetXaxis()->SetTitle("Corrected #DeltaTime (ns)");
   his44->GetYaxis()->SetTitle("Counts");
   his44->SetTitle(hisname2);
//   his44->SetMaximum(20.);
   his44->SetLineColor(1);
   his44->SetLineWidth(1.9);
   his44->SetFillStyle(1001);
   his44->SetFillColor(kGreen-7);
   
//************************************************************************************
//************************************************************************************
//1st Iteration
//************************************************************************************
//************************************************************************************
{ cout<<endl<<"************ Imo="<<imo<<"  Ise="<<ise<<"  Ila="<<ila<<" *******************************"<<endl<<endl;

  ipo2 = 0;
  for (unsigned short int ip=0; ip<ipo[imo][ise][ila]; ip++) {
  
		xaU = 1. *aU[imo][ise][ila][ip];
		xaD = 1. *aD[imo][ise][ila][ip];
		xtU = 0.01 *(tU[imo][ise][ila][ip]);
		xtD = 0.01 *(tD[imo][ise][ila][ip]);
		xZ0 = 0.01 *Z0[imo][ise][ila][ip];
		xZH = 0.01 *zH[imo][ise][ila][ip];
		xRUD= 0.01 *RUD[imo][ise][ila][ip];

	pL =(xZH-xZ0)*(xZH-xZ0);
	pL =sqrt(pL+xRUD*xRUD);

	y1[ipo2]=xtU-pL/cvel-abs(xZH-Ze[0])/vef1;
	y2[ipo2]=xtD-pL/cvel-abs(xZH-Ze[1])/vef2;
	if (y1[ipo2]>-20.&&y1[ipo2]<30.&&y2[ipo2]>-20.&&y2[ipo2]<30.&&Ene[imo][ise][ila][ip]>0.
	    &&xZH>Ze[0]&&xZH<Ze[1]) {
	  x1[ipo2]=xaU;
//	  ey1[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey1[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex1[ipo2]= 0.1;
	  z1[ipo2] = xZH;
	  
	  x2[ipo2]=xaD;
//	  ey2[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey2[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex2[ipo2]= 0.1;
	  
	  ipo2++;
	}
  }
  
  itt1 = 1.*ipo2;
   
  for (int ip=0; ip<ipo2; ip++) {
	  if (x1[ip]>10.&&x1[ip]<20.)    h1->Fill(y1[ip]);
	  if (x1[ip]>50.&&x1[ip]<100.)   h2->Fill(y1[ip]);
  
	  if (x2[ip]>10.&&x2[ip]<20.)    h3->Fill(y2[ip]);
	  if (x2[ip]>50.&&x2[ip]<100.) h4->Fill(y2[ip]);
  }
  
  cout<<endl<<"1st Iteration: "<<ipo2<<" events in the plot."<<endl<<endl;
  cont1 = ipo2;
  
  xp10 = (h2->GetMean()) - 7.;
  xp12 = log((h2->GetMean()-xp10)/(h1->GetMean()-xp10))/log(60./15.) -0.1;
  xp11 = (h2->GetMean()-xp10)/pow(60.,xp12);
  
  cout<<endl<<"xp10="<<xp10<<"  xp11="<<xp11<<"   xp12="<<xp12<<endl;
  
  xp20 = (h4->GetMean()) - 7.;
  xp22 = log((h4->GetMean()-xp20)/(h3->GetMean()-xp20))/log(60./15.);
  xp21 = (h4->GetMean()-xp20)/pow(60.,xp22);
   
  cout<<"xp20="<<xp20<<"  xp21="<<xp21<<"   xp22="<<xp22<<endl;
  
 
//*******
//********
   fgg10 = new TF1("fgg10",twfunc2,5.,2500.,6); 
       
   fgg10->SetParameter(0,2.5);  //  Pole parameter
     fgg10->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=xp10 ;
   fgg10->SetParameter(1,xp1); //  Level
       fgg10->SetParLimits(1,xp1-4.,xp1+4.);
    
   xp2=30.;    
   fgg10->SetParameter(2,xp2); //  Amp1
       fgg10->SetParLimits(2,10.,50.);
       
   xp3=-0.4;
   fgg10->SetParameter(3,xp3); //  Pow1
       fgg10->SetParLimits(3,-0.55,-0.2);
       
   xp4 = 10000.;
   fgg10->SetParameter(4,xp4); //  Amp2
       fgg10->SetParLimits(4,100.,1000000.);
     
   xp5 = -3.5;  
   fgg10->SetParameter(5,xp5); //  Pow2
       fgg10->SetParLimits(5,-5.5,-2.5);
     
   fgg10->SetLineColor(4);
   fgg10->SetLineWidth(1.);
   
   gr01 = new TGraphErrors(ipo2,x1,y1,ex1,ey1);
   gr01->Fit("fgg10","MER");
     f1p0 = fgg10->GetParameter(0); 
     f1p1 = fgg10->GetParameter(1); 
     f1p2 = fgg10->GetParameter(2); 
     f1p3 = fgg10->GetParameter(3); 
     f1p4 = fgg10->GetParameter(4); 
     f1p5 = fgg10->GetParameter(5); 
     
     cout<<"-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-"<<endl;
//********   
   fgg20 = new TF1("fgg20",twfunc2,5.,2500.,6); 
       
   fgg20->SetParameter(0,2.5);  //  Pole parameter
     fgg20->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=xp20 ;
   fgg20->SetParameter(1,xp1); //  Level
       fgg20->SetParLimits(1,xp1-4.,xp1+4.);
    
   xp2=30.;    
   fgg20->SetParameter(2,xp2); //  Amp1
       fgg20->SetParLimits(2,10.,50.);
       
   xp3=-0.4;
   fgg20->SetParameter(3,xp3); //  Pow1
       fgg20->SetParLimits(3,-0.55,-0.2);
       
   xp4 = 10000.;
   fgg20->SetParameter(4,xp4); //  Amp2
       fgg20->SetParLimits(4,100.,1000000.);
     
   xp5 = -3.5;  
   fgg20->SetParameter(5,xp5); //  Pow2
       fgg20->SetParLimits(5,-5.5,-2.5);
     
   fgg20->SetLineColor(4);
   fgg20->SetLineWidth(1.);
   
   gr11 = new TGraphErrors(ipo2,x2,y2,ex2,ey2);
   gr11->Fit("fgg20","MER");
     ff1p0 = fgg20->GetParameter(0); 
     ff1p1 = fgg20->GetParameter(1); 
     ff1p2 = fgg20->GetParameter(2); 
     ff1p3 = fgg20->GetParameter(3); 
     ff1p4 = fgg20->GetParameter(4); 
     ff1p5 = fgg20->GetParameter(5); 
     
//********   
//********

   int ipo3 = 0;
   for (int ii=0; ii<ipo2; ii++) {
     ytmp=y1[ii]-f1p1-f1p2*pow((x1[ii]-f1p0),f1p3)-f1p4*pow((x1[ii]-f1p0),f1p5);
     ytmp2=y2[ii]-ff1p1-ff1p2*pow((x2[ii]-ff1p0),ff1p3)-ff1p4*pow((x2[ii]-ff1p0),ff1p5);
     
     if (abs(ytmp)<3.&&abs(ytmp2)<3.) {
       y1[ipo3]=ytmp;
       y2[ipo3]=ytmp2;
       ey1[ipo3]=ey1[ii];
       x1[ipo3]=z1[ii];
       ex1[ipo3]=0.01;
       ipo3++; 
     }
   }
   
   fz1 = new TF1("fz1",tzfunc, 10.,420.,2); 
   fz1->SetLineColor(2);
   fz1->SetLineWidth(1.);
   
   fz1->SetParameter(0,0.);  //  Extra-Shift
     fz1->SetParLimits(0,-0.9,0.9);
   fz1->SetParameter(1,0.0); //  Slope
     fz1->SetParLimits(1,-0.02,0.02);
        
   gr02 = new TGraphErrors(ipo3,x1,y1,ex1,ey1);
   gr02->Fit("fz1","MER");
     zp0 = fz1->GetParameter(0); 
     zp1 = fz1->GetParameter(1);  
     veff2 = vef1/(1+vef1*zp1); 
     
//********

   fz2 = new TF1("fz2",tzfunc, 10.,420.,2); 
   fz2->SetLineColor(2);
   fz2->SetLineWidth(1.);
   
   fz2->SetParameter(0,0.);  //  Extra-Shift
     fz2->SetParLimits(0,-0.9,0.9);
   fz2->SetParameter(1,0.0); //  Slope
     fz2->SetParLimits(1,-0.02,0.02);
        
   gr12 = new TGraphErrors(ipo3,x1,y2,ex1,ey1);
   gr12->Fit("fz2","MER");
     zzp0 = fz2->GetParameter(0); 
     zzp1 = fz2->GetParameter(1);  
     veff22 = vef2/(1-vef2*zzp1);
     
     
delete gr01;  
delete gr02;  
delete gr11;  
delete gr12;
delete fgg10;
delete fgg20;
delete fz1;
delete fz2;

}
     if (veff2<15.7) veff2=16.2;
     if (veff2>16.7) veff2=16.2;
     if (veff22<16.2) veff22=16.8;
     if (veff22>17.2) veff22=16.8;
        
     cout<<endl<<"1st veff2="<<veff2<<"  veff22="<<veff22<<endl<<endl;

//*********************************************************************************
//*********************************************************************************
// 2nd Iteration
//*********************************************************************************
//*********************************************************************************
{   
	xkoe = 1.;
metka1:  ipo2 = 0;

    for (unsigned short int ip=0; ip<ipo[imo][ise][ila]; ip++) {
    
		xaU = 1. *aU[imo][ise][ila][ip];
		xaD = 1. *aD[imo][ise][ila][ip];
		xtU = 0.01 *(tU[imo][ise][ila][ip]);
		xtD = 0.01 *(tD[imo][ise][ila][ip]);
		xZ0 = 0.01 *Z0[imo][ise][ila][ip];
		xZH = 0.01 *zH[imo][ise][ila][ip];
		xRUD= 0.01 *RUD[imo][ise][ila][ip];

	ytmp = xtU-f1p1-f1p2*pow((xaU-f1p0),f1p3)-f1p4*pow((xaU-f1p0),f1p5);
	ytmp2= xtD-ff1p1-ff1p2*pow((xaD-ff1p0),ff1p3)-ff1p4*pow((xaD-ff1p0),ff1p5);
	
	Zh = ((ytmp-ytmp2)*veff2*veff22 + Ze[0]*veff22 + Ze[1]*veff2)/(veff2 + veff22);
	pL = (Zh-xZ0)*(Zh-xZ0);
	pL =sqrt(pL + xRUD*xRUD);
	
	ytmp = ytmp-pL/cvel-abs(Zh-Ze[0])/veff2;
	ytmp2 = ytmp2-pL/cvel-abs(Zh-Ze[1])/veff22;
	
	y1[ipo2]=xtU-pL/cvel-abs(Zh-Ze[0])/veff2;
	y2[ipo2]=xtD-pL/cvel-abs(Zh-Ze[1])/veff22;
	
			
	if (y1[ipo2]>-20.&&y1[ipo2]<30.&&y2[ipo2]>-20.&&y2[ipo2]<30.&&Ene[imo][ise][ila][ip]>0.&&abs(ytmp)<8.*xkoe&&abs(ytmp2)<8.*xkoe
	     &&Zh>Ze[0]&&Zh<Ze[1]) {
	  x1[ipo2]=xaU;
//	  ey1[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey1[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex1[ipo2]= 0.1;
	  z1[ipo2] = Zh;
	  
	  x2[ipo2]=xaD;
//	  ey2[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey2[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex2[ipo2]= 0.1;

	  ipo2++;
	}
  }
  
  itt2 = 1.*ipo2;
  if ((itt2/itt1)<0.9) {
    f2p0 = 2.5;
    f2p1 = xp10;
    f2p2 = 30.;
    f2p3 = -0.4;
    f2p4 = 100000.;
    f2p5 = -4.0;
    
    ff2p0 = 2.5;
    ff2p1 = xp20;
    ff2p2 = 30.;
    ff2p3 = -0.4;
    ff2p4 = 100000.;
    ff2p5 = -4.0;
    xkoe = xkoe * 1.5;
    cout<<"Itt2  reset: xkoe = "<<xkoe<<endl;
    goto metka1;
  }
  cout<<endl<<"2nd Iteration: "<<ipo2<<" events in the plot."<<endl<<endl;
  
//*******
   fgg10 = new TF1("fgg10",twfunc2,5.,2500.,6); 
       
   fgg10->SetParameter(0,2.5);  //  Pole parameter
     fgg10->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=f1p1 ;
   fgg10->SetParameter(1,xp1); //  Level
       fgg10->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=f1p2;    
   fgg10->SetParameter(2,xp2); //  Amp1
       fgg10->SetParLimits(2,0.8*xp2,1.2*xp2);
       
   xp3=f1p3;
   fgg10->SetParameter(3,xp3); //  Pow1
       fgg10->SetParLimits(3,xp3-0.1,xp3+0.1);
       
   xp4 = f1p4;
   fgg10->SetParameter(4,xp4); //  Amp2
       fgg10->SetParLimits(4,0.8*xp4,1.2*xp4);
     
   xp5 = f1p5;  
   fgg10->SetParameter(5,xp5); //  Pow2
       fgg10->SetParLimits(5,xp5-0.1,xp5+0.1);
     
   fgg10->SetLineColor(4);
   fgg10->SetLineWidth(1.);
   
   gr01 = new TGraphErrors(ipo2,x1,y1,ex1,ey1);
   gr01->Fit("fgg10","MER");
     f2p0 = fgg10->GetParameter(0); 
     f2p1 = fgg10->GetParameter(1); 
     f2p2 = fgg10->GetParameter(2); 
     f2p3 = fgg10->GetParameter(3); 
     f2p4 = fgg10->GetParameter(4); 
     f2p5 = fgg10->GetParameter(5);
      
     cout<<"-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-"<<endl;
//********   
   fgg20 = new TF1("fgg20",twfunc2,5.,2500.,6); 
       
   fgg20->SetParameter(0,2.5);  //  Pole parameter
     fgg20->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=ff1p1 ;
   fgg20->SetParameter(1,xp1); //  Level
       fgg20->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=ff1p2;    
   fgg20->SetParameter(2,xp2); //  Amp1
       fgg20->SetParLimits(2,0.8*xp2,1.2*xp2);
       
   xp3=ff1p3;
   fgg20->SetParameter(3,xp3); //  Pow1
       fgg20->SetParLimits(3,xp3-0.1,xp3+0.1);
       
   xp4 = ff1p4;
   fgg20->SetParameter(4,xp4); //  Amp2
       fgg20->SetParLimits(4,0.8*xp4,1.2*xp4);
     
   xp5 = ff1p5;  
   fgg20->SetParameter(5,xp5); //  Pow2
       fgg20->SetParLimits(5,xp5-0.1,xp5+0.1);
     
   fgg20->SetLineColor(4);
   fgg20->SetLineWidth(1.);
   
   gr11 = new TGraphErrors(ipo2,x2,y2,ex2,ey2);
   gr11->Fit("fgg20","MER");
     ff2p0 = fgg20->GetParameter(0); 
     ff2p1 = fgg20->GetParameter(1); 
     ff2p2 = fgg20->GetParameter(2); 
     ff2p3 = fgg20->GetParameter(3); 
     ff2p4 = fgg20->GetParameter(4); 
     ff2p5 = fgg20->GetParameter(5); 
     
//********   
//********   

   ipo3 = 0;
   for (int ii=0; ii<ipo2; ii++) {
     ytmp=y1[ii]-f2p1-f2p2*pow((x1[ii]-f2p0),f2p3)-f2p4*pow((x1[ii]-f2p0),f2p5);
     ytmp2=y2[ii]-ff2p1-ff2p2*pow((x2[ii]-ff2p0),ff2p3)-ff2p4*pow((x2[ii]-ff2p0),ff2p5);
     
     if (abs(ytmp)<3. && abs(ytmp2)<3.) {
       yy1[ipo3]=ytmp;
       yy2[ipo3]=ytmp2;
       eyy1[ipo3]=ey1[ii];
       xx1[ipo3]=z1[ii];
       exx1[ipo3]=0.1;
       ipo3++; 
     }
   }
   
   fz1 = new TF1("fz1",tzfunc, 10.,420.,2); 
   fz1->SetLineColor(2);
   fz1->SetLineWidth(1.);
   
   fz1->SetParameter(0,0.);  //  Extra-Shift
     fz1->SetParLimits(0,-0.9,0.9);
   fz1->SetParameter(1,0.0); //  Slope
     fz1->SetParLimits(1,-0.1,0.1);
        
   gr02 = new TGraphErrors(ipo3,xx1,yy1,exx1,eyy1);
   gr02->Fit("fz1","MER");
     zp0 = fz1->GetParameter(0); 
     zp1 = fz1->GetParameter(1);  
     veff2 = veff2/(1+veff2*zp1); 
     
//********

   fz2 = new TF1("fz2",tzfunc, 10.,420.,2); 
   fz2->SetLineColor(2);
   fz2->SetLineWidth(1.);
   
   fz2->SetParameter(0,0.);  //  Extra-Shift
     fz2->SetParLimits(0,-0.9,0.9);
   fz2->SetParameter(1,0.0); //  Slope
     fz2->SetParLimits(1,-0.02,0.02);
        
   gr12 = new TGraphErrors(ipo3,xx1,yy2,exx1,eyy1);
   gr12->Fit("fz2","MER");
     zzp0 = fz2->GetParameter(0); 
     zzp1 = fz2->GetParameter(1);  
     veff22 = veff22/(1-veff22*zzp1); 
     
delete gr01;  
delete gr02;  
delete gr11;  
delete gr12;
delete fgg10;
delete fgg20;
delete fz1;
delete fz2;

}
     
     if (veff2<15.7) veff2=16.2;
     if (veff2>16.7) veff2=16.2;
     if (veff22<16.2) veff22=16.8;
     if (veff22>17.2) veff22=16.8;
        
     cout<<endl<<"2nd veff2="<<veff2<<"  veff22="<<veff22<<endl<<endl;
     

//*********************************************************************************
//*********************************************************************************
// 3rd Iteration
//*********************************************************************************
//*********************************************************************************
{ 
    xkoe = 1.;
metka2: ipo2 = 0;

    for (unsigned short int ip=0; ip<ipo[imo][ise][ila]; ip++) {
    
		xaU = 1. *aU[imo][ise][ila][ip];
		xaD = 1. *aD[imo][ise][ila][ip];
		xtU = 0.01 *(tU[imo][ise][ila][ip]);
		xtD = 0.01 *(tD[imo][ise][ila][ip]);
		xZ0 = 0.01 *Z0[imo][ise][ila][ip];
		xZH = 0.01 *zH[imo][ise][ila][ip];
		xRUD= 0.01 *RUD[imo][ise][ila][ip];

	ytmp = xtU-f2p1-f2p2*pow((xaU-f2p0),f2p3)-f2p4*pow((xaU-f2p0),f2p5);
	ytmp2= xtD-ff2p1-ff2p2*pow((xaD-ff2p0),ff2p3)-ff2p4*pow((xaD-ff2p0),ff2p5);
	
	Zh = ((ytmp-ytmp2)*veff2*veff22 + Ze[0]*veff22 + Ze[1]*veff2)/(veff2 + veff22);
	pL = (Zh-xZ0)*(Zh-xZ0);
	pL =sqrt(pL + xRUD*xRUD);
	
	ytmp = ytmp-pL/cvel-abs(Zh-Ze[0])/veff2;
	ytmp2 = ytmp2-pL/cvel-abs(Zh-Ze[1])/veff22;
	
	y1[ipo2]=xtU-pL/cvel-abs(Zh-Ze[0])/veff2;
	y2[ipo2]=xtD-pL/cvel-abs(Zh-Ze[1])/veff22;
	
		
	if (y1[ipo2]>-20.&&y1[ipo2]<30.&&y2[ipo2]>-20.&&y2[ipo2]<30.&&Ene[imo][ise][ila][ip]>0.&&abs(ytmp)<5.*xkoe&&abs(ytmp2)<5.*xkoe
	    &&Zh>Ze[0]&&Zh<Ze[1]) {
	  x1[ipo2]=xaU;
//	  ey1[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey1[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex1[ipo2]= 0.1;
	  z1[ipo2] = Zh;
	  
	  x2[ipo2]=xaD;
//	  ey2[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey2[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex2[ipo2]= 0.1;

	  ipo2++;
	}
  }
  
  itt3 = 1.*ipo2;
  if ((itt3/itt2)<0.9) {
    xkoe = xkoe * 1.5;
    f2p0 = f1p0;
    f2p1 = f1p1;
    f2p2 = f1p2;
    f2p3 = f1p3;
    f2p4 = f1p4;
    f2p5 = f1p5;
    ff2p0 = ff1p0;
    ff2p1 = ff1p1;
    ff2p2 = ff1p2;
    ff2p3 = ff1p3;
    ff2p4 = ff1p4;
    ff2p5 = ff1p5;
    cout<<"Itt3  reset: xkoe = "<<xkoe<<endl;
    goto metka2;
  }
  cout<<endl<<"3rd Iteration: "<<ipo2<<" events in the plot."<<endl<<endl;
  
//*******
   fgg10 = new TF1("fgg10",twfunc2,5.,2500.,6); 
       
   fgg10->SetParameter(0,f2p0);  //  Pole parameter
     fgg10->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=f2p1 ;
   fgg10->SetParameter(1,xp1); //  Level
       fgg10->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=f2p2;    
   fgg10->SetParameter(2,xp2); //  Amp1
       fgg10->SetParLimits(2,0.9*xp2,1.1*xp2);
       
   xp3=f2p3;
   fgg10->SetParameter(3,xp3); //  Pow1
       fgg10->SetParLimits(3,xp3-0.1,xp3+0.1);
       
   xp4 = f2p4;
   fgg10->SetParameter(4,xp4); //  Amp2
       fgg10->SetParLimits(4,0.9*xp4,1.1*xp4);
     
   xp5 = f2p5;  
   fgg10->SetParameter(5,xp5); //  Pow2
       fgg10->SetParLimits(5,xp5-0.1,xp5+0.1);
     
   fgg10->SetLineColor(4);
   fgg10->SetLineWidth(1.);
   
   gr01 = new TGraphErrors(ipo2,x1,y1,ex1,ey1);
   gr01->Fit("fgg10","MER");
     f3p0 = fgg10->GetParameter(0); 
     f3p1 = fgg10->GetParameter(1); 
     f3p2 = fgg10->GetParameter(2); 
     f3p3 = fgg10->GetParameter(3); 
     f3p4 = fgg10->GetParameter(4); 
     f3p5 = fgg10->GetParameter(5);
      
     cout<<"-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-"<<endl;
//********   
   fgg20 = new TF1("fgg20",twfunc2,5.,2500.,6); 
       
   fgg20->SetParameter(0,ff2p0);  //  Pole parameter
     fgg20->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=ff2p1 ;
   fgg20->SetParameter(1,xp1); //  Level
       fgg20->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=ff2p2;    
   fgg20->SetParameter(2,xp2); //  Amp1
       fgg20->SetParLimits(2,0.9*xp2,1.1*xp2);
       
   xp3=ff2p3;
   fgg20->SetParameter(3,xp3); //  Pow1
       fgg20->SetParLimits(3,xp3-0.1,xp3+0.1);
       
   xp4 = ff2p4;
   fgg20->SetParameter(4,xp4); //  Amp2
       fgg20->SetParLimits(4,0.9*xp4,1.1*xp4);
     
   xp5 = ff2p5;  
   fgg20->SetParameter(5,xp5); //  Pow2
       fgg20->SetParLimits(5,xp5-0.1,xp5+0.1);
     
   fgg20->SetLineColor(4);
   fgg20->SetLineWidth(1.);
   
   gr11 = new TGraphErrors(ipo2,x2,y2,ex2,ey2);
   gr11->Fit("fgg20","MER");
     ff3p0 = fgg20->GetParameter(0); 
     ff3p1 = fgg20->GetParameter(1); 
     ff3p2 = fgg20->GetParameter(2); 
     ff3p3 = fgg20->GetParameter(3); 
     ff3p4 = fgg20->GetParameter(4); 
     ff3p5 = fgg20->GetParameter(5); 
     
//********   
//********   

   ipo3 = 0;
   for (int ii=0; ii<ipo2; ii++) {
     ytmp=y1[ii]-f3p1-f3p2*pow((x1[ii]-f3p0),f3p3)-f3p4*pow((x1[ii]-f3p0),f3p5);
     ytmp2=y2[ii]-ff3p1-ff3p2*pow((x2[ii]-ff3p0),ff3p3)-ff3p4*pow((x2[ii]-ff3p0),ff3p5);
     
     if (abs(ytmp)<2.5 && abs(ytmp2)<2.5) {
       yy1[ipo3]=ytmp;
       yy2[ipo3]=ytmp2;
       eyy1[ipo3]=ey1[ii];
       xx1[ipo3]=z1[ii];
       exx1[ipo3]=0.1;
       ipo3++; 
     }
   }
   
   fz1 = new TF1("fz1",tzfunc, 10.,420.,2); 
   fz1->SetLineColor(2);
   fz1->SetLineWidth(1.);
   
   fz1->SetParameter(0,0.);  //  Extra-Shift
     fz1->SetParLimits(0,-0.9,0.9);
   fz1->SetParameter(1,0.0); //  Slope
     fz1->SetParLimits(1,-0.1,0.1);
        
   gr02 = new TGraphErrors(ipo3,xx1,yy1,exx1,eyy1);
   gr02->Fit("fz1","MER");
     zp0 = fz1->GetParameter(0); 
     zp1 = fz1->GetParameter(1);  
     veff2 = veff2/(1+veff2*zp1); 
     
//********

   fz2 = new TF1("fz2",tzfunc, 10.,420.,2); 
   fz2->SetLineColor(2);
   fz2->SetLineWidth(1.);
   
   fz2->SetParameter(0,0.);  //  Extra-Shift
     fz2->SetParLimits(0,-0.9,0.9);
   fz2->SetParameter(1,0.0); //  Slope
     fz2->SetParLimits(1,-0.02,0.02);
        
   gr12 = new TGraphErrors(ipo3,xx1,yy2,exx1,eyy1);
   gr12->Fit("fz2","MER");
     zzp0 = fz2->GetParameter(0); 
     zzp1 = fz2->GetParameter(1);  
     veff22 = veff22/(1-veff22*zzp1);
     
delete gr01;  
delete gr02;  
delete gr11;  
delete gr12;
delete fgg10;
delete fgg20;
delete fz1;
delete fz2;

} 
     
     if (veff2<15.7) veff2=16.2;
     if (veff2>16.7) veff2=16.2;
     if (veff22<16.2) veff22=16.8;
     if (veff22>17.2) veff22=16.8;
        
     cout<<endl<<"3rd veff2="<<veff2<<"  veff22="<<veff22<<endl<<endl;
     
     vef3 = vef2;
     vef33= veff22;

//*********************************************************************************
//*********************************************************************************
// 4th Iteration
//*********************************************************************************
//*********************************************************************************
{ 

    xkoe = 1.;
metka3:  ipo2 = 0;

    for (unsigned short int ip=0; ip<ipo[imo][ise][ila]; ip++) {
    
		xaU = 1. *aU[imo][ise][ila][ip];
		xaD = 1. *aD[imo][ise][ila][ip];
		xtU = 0.01 *(tU[imo][ise][ila][ip]);
		xtD = 0.01 *(tD[imo][ise][ila][ip]);
		xZ0 = 0.01 *Z0[imo][ise][ila][ip];
		xZH = 0.01 *zH[imo][ise][ila][ip];
		xRUD= 0.01 *RUD[imo][ise][ila][ip];

	ytmp = xtU-f3p1-f3p2*pow((xaU-f3p0),f3p3)-f3p4*pow((xaU-f3p0),f3p5);
	ytmp2= xtD-ff3p1-ff3p2*pow((xaD-ff3p0),ff3p3)-ff3p4*pow((xaD-ff3p0),ff3p5);
	
	Zh = ((ytmp-ytmp2)*veff2*veff22 + Ze[0]*veff22 + Ze[1]*veff2)/(veff2 + veff22);
	pL = (Zh-xZ0)*(Zh-xZ0);
	pL =sqrt(pL + xRUD*xRUD);
	
	ytmp = ytmp-pL/cvel-abs(Zh-Ze[0])/veff2;
	ytmp2 = ytmp2-pL/cvel-abs(Zh-Ze[1])/veff22;
	
	y1[ipo2]=xtU-pL/cvel-abs(Zh-Ze[0])/veff2;
	y2[ipo2]=xtD-pL/cvel-abs(Zh-Ze[1])/veff22;
	
		
	if (y1[ipo2]>-20.&&y1[ipo2]<30.&&y2[ipo2]>-20.&&y2[ipo2]<30.&&Ene[imo][ise][ila][ip]>0.&&abs(ytmp)<5.*xkoe&&abs(ytmp2)<5.*xkoe
	     &&Zh>Ze[0]&&Zh<Ze[1]) {
	  x1[ipo2]=xaU;
//	  ey1[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey1[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex1[ipo2]= 0.1;
	  z1[ipo2] = Zh;
	  
	  x2[ipo2]=xaD;
//	  ey2[ipo2]=sqrt(0.2*0.2 + 0.09*0.09/Ene[ip]);
	  ey2[ipo2]=0.1+0.2/sqrt(Ene[imo][ise][ila][ip]);
	  ex2[ipo2]= 0.1;

	  ipo2++;
	}
  }
  
  itt4 = 1.*ipo2;
  if ((itt4/itt3)<0.9) {
    xkoe = xkoe * 1.5;
    f3p0 = f2p0;
    f3p1 = f2p1;
    f3p2 = f2p2;
    f3p3 = f2p3;
    f3p4 = f2p4;
    f3p5 = f2p5;
    ff3p0 = ff2p0;
    ff3p1 = ff2p1;
    ff3p2 = ff2p2;
    ff3p3 = ff2p3;
    ff3p4 = ff2p4;
    ff3p5 = ff2p5;
    cout<<"Itt4  reset: xkoe = "<<xkoe<<endl;
    goto metka3;
  }
  cout<<endl<<"4th Iteration: "<<ipo2<<" events in the plot."<<endl<<endl;
  cont4 = ipo2;  

//*******
   fgg10 = new TF1("fgg10",twfunc2,5.,2500.,6); 
       
   fgg10->SetParameter(0,f3p0);  //  Pole parameter
     fgg10->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=f3p1 ;
   fgg10->SetParameter(1,xp1); //  Level
       fgg10->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=f3p2;    
   fgg10->SetParameter(2,xp2); //  Amp1
       fgg10->SetParLimits(2,0.9*xp2,1.1*xp2);
       
   xp3=f3p3;
   fgg10->SetParameter(3,xp3); //  Pow1
       fgg10->SetParLimits(3,xp3-0.05,xp3+0.05);
       
   xp4 = f3p4;
   fgg10->SetParameter(4,xp4); //  Amp2
       fgg10->SetParLimits(4,0.9*xp4,1.1*xp4);
     
   xp5 = f3p5;  
   fgg10->SetParameter(5,xp5); //  Pow2
       fgg10->SetParLimits(5,xp5-0.05,xp5+0.05);
     
   fgg10->SetLineColor(4);
   fgg10->SetLineWidth(1.);
   
   gr01 = new TGraphErrors(ipo2,x1,y1,ex1,ey1);
   gr01->Fit("fgg10","MER");
     f4p0 = fgg10->GetParameter(0); 
     f4p1 = fgg10->GetParameter(1); 
     f4p2 = fgg10->GetParameter(2); 
     f4p3 = fgg10->GetParameter(3); 
     f4p4 = fgg10->GetParameter(4); 
     f4p5 = fgg10->GetParameter(5);
      
     cout<<"-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-"<<endl;
//********   
   fgg20 = new TF1("fgg20",twfunc2,5.,2500.,6); 
       
   fgg20->SetParameter(0,ff3p0);  //  Pole parameter
     fgg20->SetParLimits(0,0.,4.9);  //  Pole parameter

   xp1=ff3p1 ;
   fgg20->SetParameter(1,xp1); //  Level
       fgg20->SetParLimits(1,xp1-1.,xp1+1.);
    
   xp2=ff3p2;    
   fgg20->SetParameter(2,xp2); //  Amp1
       fgg20->SetParLimits(2,0.9*xp2,1.1*xp2);
       
   xp3=ff3p3;
   fgg20->SetParameter(3,xp3); //  Pow1
       fgg20->SetParLimits(3,xp3-0.05,xp3+0.05);
       
   xp4 = ff3p4;
   fgg20->SetParameter(4,xp4); //  Amp2
       fgg20->SetParLimits(4,0.9*xp4,1.1*xp4);
     
   xp5 = ff3p5;  
   fgg20->SetParameter(5,xp5); //  Pow2
       fgg20->SetParLimits(5,xp5-0.05,xp5+0.05);
     
   fgg20->SetLineColor(4);
   fgg20->SetLineWidth(1.);
   
   gr11 = new TGraphErrors(ipo2,x2,y2,ex2,ey2);
   gr11->Fit("fgg20","MER");
     ff4p0 = fgg20->GetParameter(0); 
     ff4p1 = fgg20->GetParameter(1); 
     ff4p2 = fgg20->GetParameter(2); 
     ff4p3 = fgg20->GetParameter(3); 
     ff4p4 = fgg20->GetParameter(4); 
     ff4p5 = fgg20->GetParameter(5); 
     
//********   
//********   

   ipo3 = 0;
   for (int ii=0; ii<ipo2; ii++) {
     ytmp=y1[ii]-f4p1-f4p2*pow((x1[ii]-f4p0),f4p3)-f4p4*pow((x1[ii]-f4p0),f4p5);
     ytmp2=y2[ii]-ff4p1-ff4p2*pow((x2[ii]-ff4p0),ff4p3)-ff4p4*pow((x2[ii]-ff4p0),ff4p5);
     
     if (abs(ytmp)<2.5 && abs(ytmp2)<2.5) {
       yy1[ipo3]=ytmp;
       yy2[ipo3]=ytmp2;
       eyy1[ipo3]=ey1[ii];
       xx1[ipo3]=z1[ii];
       exx1[ipo3]=0.1;
       ipo3++; 
     }
   }
   
   fz1 = new TF1("fz1",tzfunc, 10.,420.,2); 
   fz1->SetLineColor(2);
   fz1->SetLineWidth(1.);
   
   fz1->SetParameter(0,0.);  //  Extra-Shift
     fz1->SetParLimits(0,-0.9,0.9);
   fz1->SetParameter(1,0.0); //  Slope
     fz1->SetParLimits(1,-0.1,0.1);
        
   gr02 = new TGraphErrors(ipo3,xx1,yy1,exx1,eyy1);
   gr02->Fit("fz1","MER");
     zp0 = fz1->GetParameter(0); 
     zp1 = fz1->GetParameter(1);  
     veff2 = veff2/(1+veff2*zp1); 
     
//********

   fz2 = new TF1("fz2",tzfunc, 10.,420.,2); 
   fz2->SetLineColor(2);
   fz2->SetLineWidth(1.);
   
   fz2->SetParameter(0,0.);  //  Extra-Shift
     fz2->SetParLimits(0,-0.9,0.9);
   fz2->SetParameter(1,0.0); //  Slope
     fz2->SetParLimits(1,-0.02,0.02);
        
   gr12 = new TGraphErrors(ipo3,xx1,yy2,exx1,eyy1);
   gr12->Fit("fz2","MER");
     zzp0 = fz2->GetParameter(0); 
     zzp1 = fz2->GetParameter(1);  
     veff22 = veff22/(1-veff22*zzp1); 
     
     cout<<endl<<"4th veff2="<<veff2<<"  veff22="<<veff22<<endl<<endl;
     
//************************************************************************************
     
   for (int ii=0; ii<ipo2; ii++) {
//     y7[ii]=y1[ii]-zp0-zp1*(z1[ii]-212.)-fp1-fp2*pow((x1[ii]-fp0),fp3)-fp4*pow((x1[ii]-fp0),fp5);
     y10[ii]=y1[ii]-f4p1-f4p2*pow((x1[ii]-f4p0),f4p3)-f4p4*pow((x1[ii]-f4p0),f4p5);
     his4->Fill(y10[ii]);
     y20[ii]=y2[ii]-ff4p1-ff4p2*pow((x2[ii]-ff4p0),ff4p3)-ff4p4*pow((x2[ii]-ff4p0),ff4p5);
     his44->Fill(y20[ii]);
   }
//*******************************   
      
   fgau = new TF1("fgau",gfunc,-0.55,0.85,3); 
   fgau->SetNpx(100);   
   fgau->SetParameter(1,0.);  //  Mean parameter
   fgau->SetParameter(2,0.4); //  Sigma parameter
   his4->Fit("fgau","","MERQ",-0.75,0.75);
     gme = fgau->GetParameter(1);
     gsi = fgau->GetParameter(2);
     egme = fgau->GetParError(1);
     egsi = fgau->GetParError(2);
      
   fgau2 = new TF1("fgau2",gfunc,-0.55,0.85,3); 
   fgau2->SetNpx(100);   
   fgau2->SetParameter(1,0.);  //  Mean parameter
   fgau2->SetParameter(2,0.4); //  Sigma parameter
   his44->Fit("fgau2","","MERQ",-0.75,0.75);
     gme2 = fgau2->GetParameter(1);
     gsi2 = fgau2->GetParameter(2);
     egme2 = fgau2->GetParError(1);
     egsi2 = fgau2->GetParError(2);
//************************************************************************************
  dat100<<imo+1<<" "<<ise+1<<" "<<ila+1<<" 0 "
        <<f4p0<<" "<<f4p1<<" "<<f4p2<<" "<<f4p3<<" "<<f4p4<<" "<<f4p5<<" "
	<<zp0<<" "<<zp1<<" "
	<<gme<<" "<<egme<<" "<<gsi<<" "<<egsi<<" "
	<<vef3<<" "<<veff2<<endl;
  dat100<<imo+1<<" "<<ise+1<<" "<<ila+1<<" 1 "
        <<ff4p0<<" "<<ff4p1<<" "<<ff4p2<<" "<<ff4p3<<" "<<ff4p4<<" "<<ff4p5<<" "
	<<zzp0<<" "<<zzp1<<" "
	<<gme2<<" "<<egme2<<" "<<gsi2<<"  "<<egsi2<<" "
	<<vef33<<" "<<veff22<<endl;

 datveff<<imo+1<<" "<<ila+1<<" "<<ise+1<<" "<<veff2<<" "<<veff22<<endl;

 datwalk<<imo+1<<" "<<ise+1<<" "<<ila+1<<" 0 "
        <<f4p0<<" "<<f4p1+gme<<" "<<f4p2<<" "<<f4p3<<" "<<f4p4<<" "<<f4p5<<endl;
 datwalk<<imo+1<<" "<<ise+1<<" "<<ila+1<<" 1 "
        <<ff4p0<<" "<<ff4p1+gme2<<" "<<ff4p2<<" "<<ff4p3<<" "<<ff4p4<<" "<<ff4p5<<endl;

 if (abs(gme)>0.250) datwarn<<"Shifted Mean: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=0 "<<endl;
 if (abs(gme2)>0.250) datwarn<<"Shifted Mean: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=1 "<<endl;

 if (abs(gsi)>0.60) datwarn<<"Big Sigma: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=0 "<<endl;
 if (abs(gsi2)>0.60) datwarn<<"Big Sigma: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=1 "<<endl;

 if (abs(gsi)<0.10) datwarn<<"Small Sigma: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=0 "<<endl;
 if (abs(gsi2)<0.10) datwarn<<"Small Sigma: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<" E=1 "<<endl;

 if (cont4<0.5*cont1) datwarn<<"Lost Statistics: M="<<imo+1<<" S="<<ise+1<<" L="<<ila+1<<endl;

//************************************************************************************ 
  c1->Clear();
  c1->Divide(2,2);
  gStyle->SetOptFit(1110);
  c1->Range(0.1875,-81.7688,3.3125,735.919);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTitle("Datafile: ");
//  gStyle->SetPadTopMargin(.5);
//  gStyle->SetPadLeftMargin(2.5);
//  gStyle->SetPadRightMargin(2.5);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetOptFit(220);
  gStyle->SetOptStat(2220);
  gStyle->SetPalette(1);
//  gStyle->SetStatFontSize(0.5);
  
//  gStyle->SetStatW(0.35);
//  gStyle->SetStatH(0.10);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.12);
  gStyle->SetStatX(0.93);
  gStyle->SetStatY(0.92);     
  
  gStyle->SetTitleW(0.9); 
  gStyle->SetTitleH(0.08); 

// Set a bunch of parameters to make the plot look nice

/*  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);       // still leaves red frame bottom and right
  canvas->SetFrameBorderMode(0);   // need this to turn off red hist frame!
  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle(); */

//  gStyle->SetOptStat(0);
//  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.6,"t");
//  gStyle->SetTitleFontSize(0.3);
//  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
//  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)




  c1->Update();  
  c1->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(1);
  
  gr1 = new TGraph(ipo2,x1,y1);
  gr1->SetTitle(hisname);
  gr1->GetXaxis()->SetTitle("Pulse Peak");
  gr1->GetXaxis()->SetLimits(4.5,3000.0);
  gr1->SetMinimum(-15.0);
  gr1->SetMaximum(25.0);
  gr1->GetYaxis()->SetNdivisions(810);
  gr1->GetYaxis()->SetTitle("Uncorrected #DeltaTime (ns)");
  gr1->SetMarkerColor(2);
//  gr1->SetLineColor(2);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.3);
  gr1->Draw("AP");
//  fgg5->Draw("Same");
  fgg10->Draw("Same");
//hl1->Draw();
  
  
  ftxt = new char[90];
    sprintf(ftxt, "f(X) = %5.2f + %4.2f#upoint(X-%4.2f)^{%5.2f}",f4p1,f4p2,f4p0,f4p3);
  ftxt2 = new char[90];
    sprintf(ftxt2, " + %4.2f#upoint(X-%4.2f)^{%5.2f}",f4p4,f4p0,f4p5);
  txt1 = new TLatex(40.,21.,ftxt);
  txt1->SetTextFont(52);
  txt1->SetTextColor(4);
  txt1->Draw("Same");
  txt2 = new TLatex(70.,17.0,ftxt2);
  txt2->SetTextFont(52);
  txt2->SetTextColor(4);
  txt2->Draw("Same");
    
  c1->cd(2);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(1);
  
  gr3 = new TGraph(ipo2,x1,y10);
  gr3->SetTitle(hisname);
  gr3->GetXaxis()->SetTitle("Pulse Peak");
  gr3->GetXaxis()->SetLimits(4.5,3000.);
  gr3->SetMinimum(-10.0);
  gr3->SetMaximum(14.0);
  gr3->GetYaxis()->SetNdivisions(810);
  gr3->GetYaxis()->SetTitle("Corrected #DeltaTime (ns)");
  gr3->SetMarkerColor(4);
//  gr3->SetLineColor(2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.3);
  gr3->Draw("AP");
  
  lin1 = new TLine(5.,0.,2000.,0.);
  lin1->SetLineColor(2);
  lin1->Draw();
//  his1->Draw();
//hl2->Draw();
    
  c1->cd(3);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(0);
  
  gr2 = new TGraph(ipo3,xx1,yy1);
  gr2->SetTitle(hisname);
  gr2->GetXaxis()->SetTitle("Z-Coordinate (cm)");
//  gr2->GetXaxis()->SetLimits(-5.0,5.0);
  gr2->SetMinimum(-5.0);
  gr2->SetMaximum(5.0);
  gr2->GetYaxis()->SetNdivisions(810);
  gr2->GetYaxis()->SetTitle("Corrected #DeltaTime (ns)");
  gr2->SetMarkerColor(4);
//  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.3);
  gr2->Draw("AP");
  fz1->Draw("Same");
//hl3->Draw();

  ftxt3 = new char[90];
    sprintf(ftxt3, "Optimal v_{eff} = %5.2f cm/ns",veff2);
  txt3 = new TLatex(30.,3.0,ftxt3);
  txt3->SetTextFont(52);
  txt3->SetTextColor(2);
  txt3->Draw("Same");

  ftxt4 = new char[90];
    sprintf(ftxt4, "Initial v_{eff} = %5.2f cm/ns",vef1);
  txt4 = new TLatex(30.,4.0,ftxt4);
  txt4->SetTextFont(52);
  txt4->SetTextColor(4);
  txt4->Draw("Same");
  
  c1->cd(4);
//  his4->Fit("gaus","","MER",-1.,0.8);
  his4->Draw();
//  gPad->SetLogy(1);
//hl0->Draw();
    
  c1->Print(pdfname);
  
//***************************************************    
  c1->Update();  
    
  c1->cd(1);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(1);
  
  gr1 = new TGraph(ipo2,x2,y2);
  gr1->SetTitle(hisname2);
  gr1->GetXaxis()->SetTitle("Pulse Peak");
    gr1->GetXaxis()->SetLimits(4.5,3000.0);
  gr1->SetMinimum(-15.0);
  gr1->SetMaximum(25.0);
  gr1->GetYaxis()->SetNdivisions(810);
  gr1->GetYaxis()->SetTitle("Uncorrected #DeltaTime (ns)");
  gr1->SetMarkerColor(2);
//  gr1->SetLineColor(2);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.3);
  gr1->Draw("AP");
//  fgg5->Draw("Same");
  fgg20->Draw("Same");
//hl1->Draw();
  
  
  ftxt = new char[90];
    sprintf(ftxt, "f(X) = %5.2f + %4.2f#upoint(X-%4.2f)^{%5.2f}",ff4p1,ff4p2,ff4p0,ff4p3);
  ftxt2 = new char[90];
    sprintf(ftxt2, " + %4.2f#upoint(X-%4.2f)^{%5.2f}",ff4p4,ff4p0,ff4p5);
  txt1 = new TLatex(40.,21.,ftxt);
  txt1->SetTextFont(52);
  txt1->SetTextColor(4);
  txt1->Draw("Same");
  txt2 = new TLatex(70.,17.,ftxt2);
  txt2->SetTextFont(52);
  txt2->SetTextColor(4);
  txt2->Draw("Same");
    
  c1->cd(2);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(1);
  
  TGraph *gr3 = new TGraph(ipo2,x2,y20);
  gr3->SetTitle(hisname2);
  gr3->GetXaxis()->SetTitle("Pulse Peak");
  gr3->GetXaxis()->SetLimits(4.5,3000.0);
  gr3->SetMinimum(-10.0);
  gr3->SetMaximum(14.0);
  gr3->GetYaxis()->SetNdivisions(810);
  gr3->GetYaxis()->SetTitle("Corrected #DeltaTime (ns)");
  gr3->SetMarkerColor(4);
//  gr3->SetLineColor(2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.3);
  gr3->Draw("AP");
  
  lin1 = new TLine(5.,0.,2000.,0.);
  lin1->SetLineColor(2);
  lin1->Draw();
    
  c1->cd(3);
  gPad->SetGridy(0);
  gPad->SetGridx(0);
  gPad->SetFillColor(0);
  gPad->SetTicks();
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.17);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogx(0);
  
  gr2 = new TGraph(ipo3,xx1,yy2);
  gr2->SetTitle(hisname2);
  gr2->GetXaxis()->SetTitle("Z-Coordinate (cm)");
//  gr2->GetXaxis()->SetLimits(-5.0,5.0);
  gr2->SetMinimum(-5.0);
  gr2->SetMaximum(5.0);
  gr2->GetYaxis()->SetNdivisions(810);
  gr2->GetYaxis()->SetTitle("Corrected #DeltaTime (ns)");
  gr2->SetMarkerColor(4);
//  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.3);
  gr2->Draw("AP");
  fz1->Draw("Same");
//hl3->Draw();

  ftxt3 = new char[90];
    sprintf(ftxt3, "Optimal v_{eff} = %5.2f cm/ns",veff22);
  txt3 = new TLatex(30.,3.0,ftxt3);
  txt3->SetTextFont(52);
  txt3->SetTextColor(2);
  txt3->Draw("Same");

  ftxt4 = new char[90];
    sprintf(ftxt4, "Initial v_{eff} = %5.2f cm/ns",vef2);
  txt4 = new TLatex(30.,4.0,ftxt4);
  txt4->SetTextFont(52);
  txt4->SetTextColor(4);
  txt4->Draw("Same");
  
  c1->cd(4);
//  his4->Fit("gaus","","MER",-1.,0.8);
  his44->Draw();
//  gPad->SetLogy(1);
//hl0->Draw();
    
  c1->Print(pdfname);
  
delete gr01;  
delete gr02;  
delete gr11;  
delete gr12;
delete fgg10;
delete fgg20;
delete fz1;
delete fz2;

 
delete fgau;
delete fgau2;
delete h1; 
delete h2; 
delete h3; 
delete h4; 
delete his4; 
delete his44; 
delete gr1;  
delete gr2;  
delete gr3; 
delete ftxt; 
delete ftxt2; 
delete ftxt3; 
delete ftxt4;
delete txt1; 
delete txt2; 
delete txt3; 
delete txt4;
delete hisname; 
delete hisname2;
delete lin1; 
  
} //End of 4th Iteration

} //End of ise loop
} //End of ila loop
} //End of imo loop
  
//*************************************************************  
  c1->Update();    
  c1->Print(pdfname2);
  dat100.close(); 
  datveff.close();
  datwalk.close();
  datwarn.close();  
//************************************************************************************
	
	return NOERROR;
}
