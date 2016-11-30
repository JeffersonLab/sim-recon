// $Id$
//
//    File: JEventProcessor_FDC_Residuals.cc
// Created: Wed Nov 30 10:08:00 EDT 2016
// Creator: aaustreg
//

#include "JEventProcessor_FDC_Residuals.h"
using namespace jana;
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "HistogramTools.h"

// extracted from run11553
const double deflect[24]={0.256449, 0.257153, 0.259054, 0.257982, 0.257663, 0.252944,
			  0.260737, 0.265413, 0.262766, 0.263953, 0.265173, 0.264403,
			  0.260384, 0.268576, 0.2578, 0.262166, 0.264653, 0.256705,
			  0.226083, 0.233847, 0.222489, 0.220759, 0.222032, 0.212103};

const double correct[24]={0.00555047, 0.00471467, 0.00451803, 0.00695173, 0.00838218, 0.00916391,
			  0.015701, 0.0179723, 0.0186109, 0.0179447, 0.0186627, 0.0151424,
			  0.0224091, 0.0278104, 0.0297215, 0.0257817, 0.0248605, 0.0266139, 
			  0.0244911, 0.0248676, 0.0294895,  0.0273482, 0.019603, 0.0214564};

//For extraction of magnetic field slope, split the detectors into bins in radius
const unsigned int rad = 1; // 1, 5 or 9

static TH1I *hChi2OverNDF;
static TH1I *hChi2OverNDF_accepted;
static TH1I *hPseudoRes;
// static TH1I *hPseudoResX[25];
// static TH1I *hPseudoResY[25];
static TH1I *hPseudoResU[25];
static TH1I *hPseudoResV[25];
static TH2I *hPseudoResUvsV[25][rad];
static TH2I *hPseudoResVsT[25];
static TH2I *hPseudoResSvsQ;
static TH2I *hResVsT[25];
static TH1I *hMom;
static TH1I *hMom_accepted;
static TH1I *hTheta;
static TH1I *hTheta_accepted;
static TH1I *hPhi;
static TH1I *hPhi_accepted;
static TH1I *hCellsHit;
static TH1I *hCellsHit_accepted;
static TH1I *hRingsHit;
static TH1I *hRingsHit_accepted;
static TH2I *hRingsHit_vs_P;

static TH1I *hWireTime[25];
static TH1I *hWireTime_accepted[25];
static TH1I *hCathodeTime[25];
static TH1I *hCathodeTime_accepted[25];
static TH1I *hPseudoTime[25];
static TH1I *hPseudoTime_accepted[25];
static TH1I *hPullTime[25];
static TH1I *hDeltaTime[25];

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FDC_Residuals());
}
} // "C"


//------------------
// JEventProcessor_FDC_Residuals (Constructor)
//------------------
JEventProcessor_FDC_Residuals::JEventProcessor_FDC_Residuals()
{
    ;
}

//------------------
// ~JEventProcessor_FDC_Residuals (Destructor)
//------------------
JEventProcessor_FDC_Residuals::~JEventProcessor_FDC_Residuals()
{
    ;
}

//------------------
// init
//------------------
jerror_t JEventProcessor_FDC_Residuals::init(void)
{
  // create root folder for fdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("FDC_Residuals")->cd();

  gDirectory->mkdir("Track_Quality")->cd();

  hChi2OverNDF = new TH1I("hChi2OverNDF","hChi2OverNDF", 500, 0.0, 1.0);
  hChi2OverNDF_accepted = new TH1I("hChi2OverNDF_accepted","hChi2OverNDF_accepted", 500, 0.0, 1.0);
  hMom = new TH1I("hMom","Momentum", 500, 0.0, 10);
  hMom_accepted = new TH1I("hMom_accepted","Momentum (accepted)", 500, 0.0, 10);
  hTheta = new TH1I("hTheta","Theta of Track", 500, 0.0, 50);
  hTheta_accepted = new TH1I("hTheta_accepted","Theta of Track (accepted)", 500, 0.0, 50);
  hPhi = new TH1I("hPhi","Phi of Track", 360, -180, 180);
  hPhi_accepted = new TH1I("hPhi_accepted","Phi of Track (accepted)", 360, -180, 180);
  hCellsHit = new TH1I("hCellsHit", "Cells of FDC Contributing to Track", 25, -0.5, 24.5);
  hCellsHit_accepted = new TH1I("hCellsHit_accepted", "Cells of FDC Contributing to Track (accepted)", 25, -0.5, 24.5);
  hRingsHit = new TH1I("hRingsHit", "Rings of CDC Contributing to Track", 29, -0.5, 28.5);
  hRingsHit_accepted = new TH1I("hRingsHit_accepted", "Rings of CDC Contributing to Track (accepted)", 25, -0.5, 24.5);

  hRingsHit_vs_P = new TH2I("hRingsHit_vs_P", "Rings of CDC Contributing to Track vs P (accepted)", 25, -0.5, 24.5, 34, 0.6, 4);

  gDirectory->cd("/FDC_Residuals");
  gDirectory->mkdir("Residuals")->cd();
  hPseudoRes = new TH1I("hPseudoRes","Pseudo Residual in R", 500, 0, 5);

  hPseudoResSvsQ = new TH2I("hPseudoResSvsQ","Pseudo Residual along Wire vs. Charge", 500, 0, 50, 200, -0.1, 0.1);

  for(int icell=0; icell<24; icell++){

    char hname_X[256];
    char hname_Y[256];
    char hname_XY[256];
   
    // sprintf(hname_X, "hPseudoResX_cell[%d]", icell+1);
    // sprintf(hname_Y, "hPseudoResY_cell[%d]", icell+1);
    // hPseudoResX[icell+1] = new TH1I(hname_X,"Pseudo Residual in X", 600, -3, 3);
    // hPseudoResY[icell+1] = new TH1I(hname_Y,"Pseudo Residual in Y", 600, -3, 3);

    sprintf(hname_X, "hPseudoResU_cell[%d]", icell+1);
    sprintf(hname_Y, "hPseudoResV_cell[%d]", icell+1);
    hPseudoResU[icell+1] = new TH1I(hname_X,"Pseudo Residual along Wire", 600, -3, 3);
    hPseudoResV[icell+1] = new TH1I(hname_Y,"Pseudo Residual perp. to Wire", 600, -3, 3);
    for (unsigned int r=0; r<rad; r++){
      sprintf(hname_XY, "hPseudoResUvsV_cell[%d]_radius[%d]", icell+1, (r+1)*45/rad);
      hPseudoResUvsV[icell+1][r] = new TH2I(hname_XY,"Pseudo Residual 2D", 200, -1, 1, 200, -1, 1);
    }

    sprintf(hname_X, "hResVsT_cell[%d]", icell+1);
    sprintf(hname_Y, "hPseudoResVsT[%d]", icell+1);
    hResVsT[icell+1] = new TH2I(hname_X,"Tracking Residual (Biased) Vs Wire Drift Time; DOCA (cm); Drift Time (ns)", 100, 0, 0.5, 400, -100, 300);
    hPseudoResVsT[icell+1] = new TH2I(hname_Y,"Tracking Residual (Biased) Vs Pseudo Time; Residual (cm); Drift Time (ns)", 200, -0.5, 0.5, 400, -100, 300);

    sprintf(hname_X, "hWireTime_cell[%d]", icell+1);
    sprintf(hname_Y, "hWireTime_acc_cell[%d]", icell+1);
    hWireTime[icell+1] = new TH1I(hname_X, "Timing of Hits", 600, -200, 400);
    hWireTime_accepted[icell+1] = new TH1I(hname_Y, "Timing of Hits (accepted)", 600, -200, 400);
    
    sprintf(hname_X, "hCathodeTime_cell[%d]", icell+1);
    sprintf(hname_Y, "hCathodeTime_acc_cell[%d]", icell+1);
    hCathodeTime[icell+1] = new TH1I(hname_X, "Timing of Cathode Clusters", 600, -200, 400);
    hCathodeTime_accepted[icell+1] = new TH1I(hname_Y, "Timing of Cathode Clusters (accepted)", 600, -200, 400);
    
    sprintf(hname_X, "hPseudoTime_cell[%d]", icell+1);
    sprintf(hname_Y, "hPseudoTime_acc_cell[%d]", icell+1);
    hPseudoTime[icell+1] = new TH1I(hname_X, "Timing of Pseudos", 600, -200, 400);
    hPseudoTime_accepted[icell+1] = new TH1I(hname_Y, "Timing of Pseudos (accepted)", 600, -200, 400);

    sprintf(hname_X, "hPullTime_cell[%d]", icell+1);
    hPullTime[icell+1] = new TH1I(hname_X, "Timing of Wire Hits (from pulls)", 600, -200, 400);
    
    sprintf(hname_X, "hDeltaTime_cell[%d]", icell+1);
    hDeltaTime[icell+1] = new TH1I(hname_X, "Time Difference between Wires and Cathode Clusters", 200, -50, 50);

  }

  main->cd();

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FDC_Residuals::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);
  //JCalibration *jcalib = dapp->GetJCalibration(runnumber);
  dgeom  = dapp->GetDGeometry(runnumber);
  //bfield = dapp->GetBfield();

  // Get the position of the FDC wires from DGeometry
  dgeom->GetFDCWires(fdcwires);

  // Get the z position of the FDC planes from DGeometry
  dgeom->GetFDCZ(fdcz);

  // Get inner radii for FDC packages
  dgeom->GetFDCRmin(fdcrmin);

  // Get outer radius of FDC 
  dgeom->GetFDCRmax(fdcrmax);
  fdcrmax = 48; // fix to 48cm from DFDCPseudo_factory
		  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FDC_Residuals::evnt(JEventLoop *loop, uint64_t eventnumber){

  vector< const DFDCHit *> locFDCHitVector;
  loop->Get(locFDCHitVector);
  vector< const DFDCPseudo *> locFDCPseudoVector;
  loop->Get(locFDCPseudoVector);

  //Pre-sort Hits to save time
  //only need to search within the given cell, wire
  map<int, map<int, set<const DFDCHit*> > > locSortedFDCHits; //first int: cell //second int: wire
  for( unsigned int hitNum = 0; hitNum < locFDCHitVector.size(); hitNum++){
    const DFDCHit * locHit = locFDCHitVector[hitNum];
    if (locHit->plane != 2) continue; // only wires!
    
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    hWireTime[locHit->gLayer]->Fill(locHit->t);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    // cut on timing of hits
    //if (-100 > locHit->t || locHit->t > 300) continue;
    
    locSortedFDCHits[locHit->gLayer][locHit->element].insert(locHit);

  }

  //Pre-sort PseudoHits to save time
  //only need to search within the given cell
  map<int, set<const DFDCPseudo*> > locSortedFDCPseudos; // int: cell
  for( unsigned int hitNum = 0; hitNum < locFDCPseudoVector.size(); hitNum++){
    const DFDCPseudo * locPseudo = locFDCPseudoVector[hitNum];
    locSortedFDCPseudos[locPseudo->wire->layer].insert(locPseudo);

    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    hPseudoTime[locPseudo->wire->layer]->Fill(locPseudo->time);
    hCathodeTime[locPseudo->wire->layer]->Fill(locPseudo->t_u);
    hCathodeTime[locPseudo->wire->layer]->Fill(locPseudo->t_v);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  }
  
  const DDetectorMatches *detMatches;
  vector<DSCHitMatchParams> SCMatches;
  vector<DTOFHitMatchParams> TOFMatches;
  vector<DBCALShowerMatchParams> BCALMatches;

  if(!dIsNoFieldFlag){
    loop->GetSingle(detMatches);
  }

  vector <const DChargedTrack *> chargedTrackVector;
  loop->Get(chargedTrackVector);

  for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

    const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

    // Cut very loosely on the track quality
    const DTrackTimeBased *thisTimeBasedTrack;

    bestHypothesis->GetSingle(thisTimeBasedTrack);

    // First loop over the FDC/CDC hits of the track to find out how many cells/rings were hit
    // For the first track selection, use already pseudo hits from fdc
    // Determine how many cells (FDC) and rings (CDC) contribute to track in total, and in which package
    vector< int > cellsHit;
    vector< int > ringsHit;
    vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
    for (unsigned int i = 0; i < pulls.size(); i++){
      const DFDCPseudo * thisTrackFDCHit = pulls[i].fdc_hit;
      if (thisTrackFDCHit != NULL){      
	hPullTime[thisTrackFDCHit->wire->layer]->Fill(pulls[i].tdrift);
	if ( find(cellsHit.begin(), cellsHit.end(), thisTrackFDCHit->wire->layer) == cellsHit.end())
	  cellsHit.push_back(thisTrackFDCHit->wire->layer);
      }
      
      const DCDCTrackHit * thisTrackCDCHit = pulls[i].cdc_hit;
      if (thisTrackCDCHit != NULL){
	if ( find(ringsHit.begin(), ringsHit.end(), thisTrackCDCHit->wire->ring) == ringsHit.end())
	  ringsHit.push_back(thisTrackCDCHit->wire->ring);
      }
    }
    unsigned int cells = cellsHit.size();
    unsigned int rings = ringsHit.size();
    
    if (cells == 0) continue;

    double theta_deg = thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg();
    double phi_deg = thisTimeBasedTrack->momentum().Phi()*TMath::RadToDeg();
    double tmom = thisTimeBasedTrack->pmag();
    
    // Fill Histograms for all Tracks    
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    hCellsHit->Fill(cells);
    hRingsHit->Fill(rings);
    hChi2OverNDF->Fill(thisTimeBasedTrack->FOM);
    hMom->Fill(tmom);
    hTheta->Fill(theta_deg);
    hPhi->Fill(phi_deg);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    // All Cuts on Track Quality:
    
    // if (thisTimeBasedTrack->FOM < 5.73303E-7) 
    //   continue; // 5 sigma cut from Paul
    if (thisTimeBasedTrack->FOM < 1E-20)
      continue;  // from CDC analysis

    if(!dIsNoFieldFlag){ // Quality cuts for Field on runs.
      if(tmom < 0.6)
       	continue; // Cut on the reconstructed momentum below 600MeV, no curlers
      
      // TRY TO IMPROVE RECONSTUCTION FOR LOW MOMENTA WITH CDC INFO, DISABLED
      // Tuned on 1345A data, not tuned for 1200A
      // if (tmom < 1.2 && rings < 12)
      //  	continue; // for Low-Momentum Tracks < 2 GeV, check number of hits in CDC
      // if (tmom < 1.4 && rings < 10)
      //  	continue; // for Low-Momentum Tracks < 2 GeV, check number of hits in CDC
      // if (tmom < 1.6 && rings < 8)
      //  	continue; // for Low-Momentum Tracks < 2 GeV, check number of hits in CDC
      // if (tmom < 1.8 && rings < 6)
      //  	continue; // for Low-Momentum Tracks < 2 GeV, check number of hits in CDC
      // if (tmom < 2 && rings < 4)
      //  	continue; // for Low-Momentum Tracks < 2 GeV, check number of hits in CDC
      
      
      if(!detMatches->Get_SCMatchParams(thisTimeBasedTrack, SCMatches)) {
	continue; // At least one match to the Start Counter
      }
      
      if(!detMatches->Get_TOFMatchParams(thisTimeBasedTrack, TOFMatches) && !detMatches->Get_BCALMatchParams(thisTimeBasedTrack,BCALMatches)) {
      	continue; // At least one match either to the Time-Of-Flight (forward) OR the BCAL (large angles)
      }
    }
    else return NOERROR; // Field off not supported for now !!
    
    // count how many hits in each package (= 6 cells)
    unsigned int packageHit[4] = {0,0,0,0};
    for (unsigned int i = 0; i < cells; i++)
      packageHit[(cellsHit[i] - 1) / 6]++;
    
    unsigned int minCells = 4; //At least 4 cells hit in any package for relatively "unbiased" efficiencies
    if (packageHit[0] < minCells && packageHit[1] < minCells && packageHit[2] < minCells && packageHit[3] < minCells) continue;
    
    // Fill Histograms for accepted Tracks
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    hCellsHit_accepted->Fill(cells);
    //if (tmom < 2)
    hRingsHit_accepted->Fill(rings);
    hChi2OverNDF_accepted->Fill(thisTimeBasedTrack->FOM);
    hMom_accepted->Fill(tmom);
    hTheta_accepted->Fill(theta_deg);
    hPhi_accepted->Fill(phi_deg);
    hRingsHit_vs_P->Fill(rings, tmom);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
    
    // Start efficiency computation with these tracks
    
    // Loop over cells (24)
    for (unsigned int cellIndex = 0; cellIndex < fdcwires.size(); cellIndex ++){
      unsigned int cellNum = cellIndex +1;
      vector< DFDCWire * > wireByNumber = fdcwires[cellIndex];

      // Use only tracks that have at least 4 (or other # of) hits in the package this cell is in
      // OR 12 hits in total
      if (packageHit[cellIndex / 6] < minCells /*&& cells < 12*/) continue;

      // Interpolate track to layer
      DVector3 plane_origin(0.0, 0.0, fdcz[cellIndex]);
      DVector3 plane_normal(0.0, 0.0, 1.0);
      DVector3 interPosition;
      thisTimeBasedTrack->rt->GetIntersectionWithPlane(plane_origin, plane_normal, interPosition);
      
      // cut out central hole and intersections at too large radii
      int packageIndex = (cellNum - 1) / 6;
      if (interPosition.Perp() < fdcrmin[packageIndex] || interPosition.Perp() > fdcrmax) continue;
      

      /////////////////// WIRE ANALYSIS /////////////////////////////////
	
      for (unsigned int wireIndex = 0; wireIndex < wireByNumber.size(); wireIndex++){
	unsigned int wireNum = wireIndex+1;
	DFDCWire * wire = wireByNumber[wireIndex]; 
	double wireLength = wire->L;
	double distanceToWire = thisTimeBasedTrack->rt->DistToRT(wire, &wireLength);
	bool expectHit = false;

	// starting from here, only histograms with distance to wire < 0.5, maybe change later
	if (distanceToWire < 0.5 )
	  expectHit = true;

	//SKIP IF NOT CLOSE (from Paul)
	if(distanceToWire > 50.0)
	  {
	    wireIndex += 30;
	    continue;
	  }
	if(distanceToWire > 20.0)
	  {
	    wireIndex += 10;
	    continue;
	  }
	if(distanceToWire > 10.0)
	  {
	    wireIndex += 5;
	    continue;
	  }

	  
	if (expectHit && cellNum < 25){
	  
	  // look in the presorted FDC Hits for a match
	  if(!locSortedFDCHits[cellNum][wireNum].empty() || !locSortedFDCHits[cellNum][wireNum-1].empty() || !locSortedFDCHits[cellNum][wireNum+1].empty()){
	    // Look not only in one wire, but also in adjacent ones (?)
	    // This can remove the dependence on the track error

	    // Fill the histograms only for the cell in the center
	    if(!locSortedFDCHits[cellNum][wireNum].empty()){
	      // Loop over multiple hits in the wire
	      for(set<const DFDCHit*>::iterator locIterator = locSortedFDCHits[cellNum][wireNum].begin();  locIterator !=  locSortedFDCHits[cellNum][wireNum].end(); ++locIterator){
	      	const DFDCHit* locHit = * locIterator;
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		hWireTime_accepted[cellNum]->Fill(locHit->t);
		hResVsT[cellNum]->Fill(distanceToWire, locHit->t);
		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	      }
	    }
	    
	  }
	  
	  break; // break if 1 expected hit was found, 2 are geometrically not possible (speedup!)
	  
	} // End expected
      } // End wire loop
      
	// END WIRE ANALYSIS
      

      
	///////////// START PSEUDO ANALYSIS /////////////////////////////////////////
      
      
      bool foundPseudo = false;
      if(!locSortedFDCPseudos[cellNum].empty()){

	set<const DFDCPseudo*>& locPlanePseudos = locSortedFDCPseudos[cellNum];
      	// loop over the FDC Pseudo Hits in that cell to look for a match
	for(set<const DFDCPseudo*>::iterator locIterator = locPlanePseudos.begin();  locIterator !=  locPlanePseudos.end(); ++locIterator)
	  {
	    const DFDCPseudo * locPseudo = * locIterator;
	    if (locPseudo == NULL) continue;

	    DVector2 pseudoPosition = locPseudo->xy;
	    DVector2 interPosition2D(interPosition.X(), interPosition.Y());

	    DVector2 residual2D = (pseudoPosition - interPosition2D);
	    double residualR = residual2D.Mod();
	    // double residualX = pseudoPosition.X() - interPosition2D.X();
	    // double residualY = pseudoPosition.Y() - interPosition2D.Y();

	    // rotate to wire coordinate system
	    const DFDCWire * wire = locPseudo->wire;
	    double residualU = -1*(residual2D.Rotate(wire->angle - deflect[cellIndex] - correct[cellIndex])).X();
	    double residualV = -1*(residual2D.Rotate(wire->angle - deflect[cellIndex] - correct[cellIndex])).Y();

	    // these can be used for background studies/correction
	    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	    hPseudoRes->Fill(residualR);
	    // hPseudoResX[cellNum]->Fill(residualX);
	    // hPseudoResY[cellNum]->Fill(residualY);
	    hPseudoResU[cellNum]->Fill(residualU);
	    hPseudoResV[cellNum]->Fill(residualV);
	    
	    if (locPseudo->status == 6)
	      hPseudoResSvsQ->Fill(locPseudo->q,residualV);
	    
	    unsigned int radius = interPosition2D.Mod()/(45/rad);
	    if (radius<rad)
	      hPseudoResUvsV[cellNum][radius]->Fill(residualU, residualV);
	    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	    
	    if (foundPseudo) continue; 
	    // to avoid double conting
	    	    
	    //if (residualR < 1.6){ // = 5 sigma in x and in y
	    if (residualR < 2){ // = to account for non-gaussian tails and tracking/extrapolation errors
	      foundPseudo = true;
	      
	      if(cellNum < 25){
		// fill histogramms with the predicted, not with the measured position
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		hPseudoTime_accepted[cellNum]->Fill(locPseudo->time);
		hCathodeTime_accepted[cellNum]->Fill(locPseudo->t_u);
		hCathodeTime_accepted[cellNum]->Fill(locPseudo->t_v);
		hDeltaTime[cellNum]->Fill(locPseudo->time - locPseudo->t_u);
		hDeltaTime[cellNum]->Fill(locPseudo->time - locPseudo->t_v);
		hPseudoResVsT[cellNum]->Fill(residualU, locPseudo->time);
		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	      }
	    }
	    
	    
	  } // End Pseudo Hit loop
      }
      
      // END PSEUDO ANALYSIS
      
    } // End cell loop
  } // End track loop
   
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FDC_Residuals::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FDC_Residuals::fini(void)
{
  return NOERROR;
}

