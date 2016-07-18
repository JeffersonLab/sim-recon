// $Id$
//
//    File: JEventProcessor_FDC_Efficiency.cc
// Created: Thu May 26 15:57:38 EDT 2016
// Creator: aaustreg
//

#include "JEventProcessor_FDC_Efficiency.h"
using namespace jana;
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "HistogramTools.h"

static TH1D *fdc_wire_measured_cell[25]; //Filled with total actually detected before division at end
static TH1D *fdc_wire_expected_cell[25]; // Contains total number of expected hits by DOCA

static TH2D *fdc_pseudo_measured_cell[25]; //Filled with total actually detected before division at end
static TH2D *fdc_pseudo_expected_cell[25]; // Contains total number of expected hits by DOCA

static TH2D *fdc_pseudo_residualU_cell[25]; //Filled with residual mean perp to Wire
static TH2D *fdc_pseudo_residualV_cell[25]; //Filled with residual mean along Wire

const unsigned int rad=1; // 1, 5 or 9

static TH1I *hChi2OverNDF;
static TH1I *hChi2OverNDF_accepted;
static TH1I *hPseudoRes;
static TH1I *hPseudoResX[25];
static TH1I *hPseudoResY[25];
static TH1I *hPseudoResU[25];
static TH1I *hPseudoResV[25];
static TH2I *hPseudoResUvsV[25][rad];
static TH2I *hPseudoResVsT;
static TH2I *hResVsT;
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

static TH1I *hWireTime;
static TH1I *hWireTime_accepted;
static TH1I *hPseudoTime;
static TH1I *hPseudoTime_accepted;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FDC_Efficiency());
}
} // "C"


//------------------
// JEventProcessor_FDC_Efficiency (Constructor)
//------------------
JEventProcessor_FDC_Efficiency::JEventProcessor_FDC_Efficiency()
{
    ;
}

//------------------
// ~JEventProcessor_FDC_Efficiency (Destructor)
//------------------
JEventProcessor_FDC_Efficiency::~JEventProcessor_FDC_Efficiency()
{
    ;
}

//------------------
// init
//------------------
jerror_t JEventProcessor_FDC_Efficiency::init(void)
{
  // create root folder for fdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("FDC_Efficiency")->cd();
  gDirectory->mkdir("FDC_View")->cd();

  for(int icell=0; icell<24; icell++){
      
    char hname_measured[256];
    char hname_expected[256];
    sprintf(hname_measured, "fdc_wire_measured_cell[%d]", icell+1);
    sprintf(hname_expected, "fdc_wire_expected_cell[%d]", icell+1);
    fdc_wire_measured_cell[icell+1] = new TH1D(hname_measured, "", 96, 0, 95);
    fdc_wire_expected_cell[icell+1] = new TH1D(hname_expected, "", 96, 0, 95);

    sprintf(hname_measured, "fdc_pseudo_measured_cell[%d]", icell+1);
    sprintf(hname_expected, "fdc_pseudo_expected_cell[%d]", icell+1);
    fdc_pseudo_measured_cell[icell+1] = new TH2D(hname_measured, "", 100, -50, 50, 100, -50, 50);
    fdc_pseudo_expected_cell[icell+1] = new TH2D(hname_expected, "", 100, -50, 50, 100, -50, 50);

  }	

  gDirectory->cd("/FDC_Efficiency");
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

  hWireTime = new TH1I("hWireTime", "Timing of Hits", 500, -500, 1000);
  hWireTime_accepted = new TH1I("hWireTime_accepted", "Timing of Hits (accepted)", 500, -500, 1000);

  hPseudoTime = new TH1I("hPseudoTime", "Timing of Pseudos", 500, -500, 1000);
  hPseudoTime_accepted = new TH1I("hPseudoTime_accepted", "Timing of Pseudos (accepted)", 500, -500, 1000);


  gDirectory->cd("/FDC_Efficiency");
  gDirectory->mkdir("Residuals")->cd();
  hResVsT = new TH2I("hResVsT","Tracking Residual (Biased) Vs Wire Drift Time; DOCA (cm); Drift Time (ns)", 100, 0, 0.5, 600, -100, 500);
  hPseudoResVsT = new TH2I("hPseudoResVsT","Tracking Residual (Biased) Vs Pseudo Time; Residual (cm); Drift Time (ns)", 100, 0, 0.5, 400, -100, 300);
  hPseudoRes = new TH1I("hPseudoRes","Pseudo Residual in R", 500, 0, 5);

  for(int icell=0; icell<24; icell++){

    char hname_X[256];
    char hname_Y[256];
    char hname_XY[256];
   
    sprintf(hname_X, "hPseudoResX_cell[%d]", icell+1);
    sprintf(hname_Y, "hPseudoResY_cell[%d]", icell+1);
    hPseudoResX[icell+1] = new TH1I(hname_X,"Pseudo Residual in X", 600, -3, 3);
    hPseudoResY[icell+1] = new TH1I(hname_Y,"Pseudo Residual in Y", 600, -3, 3);

    sprintf(hname_X, "hPseudoResU_cell[%d]", icell+1);
    sprintf(hname_Y, "hPseudoResV_cell[%d]", icell+1);
    hPseudoResU[icell+1] = new TH1I(hname_X,"Pseudo Residual along Wire", 600, -3, 3);
    hPseudoResV[icell+1] = new TH1I(hname_Y,"Pseudo Residual perp. to Wire", 600, -3, 3);
    for (int r=0; r<rad; r++){
      sprintf(hname_XY, "hPseudoResUvsV_cell[%d]_radius[%d]", icell+1, (r+1)*45/rad);
      hPseudoResUvsV[icell+1][r] = new TH2I(hname_XY,"Pseudo Residual 2D", 200, -1, 1, 200, -1, 1);
    }

    sprintf(hname_X, "fdc_pseudo_residualU_cell[%d]", icell+1);
    sprintf(hname_Y, "fdc_pseudo_residualV_cell[%d]", icell+1);
    fdc_pseudo_residualU_cell[icell+1] = new TH2D(hname_X, "", 100, -50, 50, 100, -50, 50);
    fdc_pseudo_residualV_cell[icell+1] = new TH2D(hname_Y, "", 100, -50, 50, 100, -50, 50);

  }

  main->cd();

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FDC_Efficiency::brun(JEventLoop *eventLoop, int32_t runnumber)
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
		  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FDC_Efficiency::evnt(JEventLoop *loop, uint64_t eventnumber){

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
    hWireTime->Fill(locHit->t);
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
    hPseudoTime->Fill(locPseudo->time);
    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  }
  
  const DDetectorMatches *detMatches;
  //vector<DSCHitMatchParams> SCMatches;
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
      if(tmom < .6)
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
      
      
      // if(!detMatches->Get_SCMatchParams(thisTimeBasedTrack, SCMatches)) {
      //  	continue; // At least one match to the Start Counter DISABLED
      // }
      
      if(!detMatches->Get_TOFMatchParams(thisTimeBasedTrack, TOFMatches) && !detMatches->Get_BCALMatchParams(thisTimeBasedTrack,BCALMatches)) {
	continue; // At least one match either to the Time-Of-Flight (forward) OR the BCAL (large angles)
      }
    }
    else return NOERROR; // Field off not supported for now !!
    
    // count how many hits in each package (= 6 cells)
    unsigned int packageHit[4] = {0,0,0,0};
    for (unsigned int i = 0; i < cells; i++)
      packageHit[(cellsHit[i] - 1) / 6]++;

    // for alignment studies (Lubomir)
    // if (cells < 17) continue;
    
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
      if (packageHit[cellIndex / 6] < minCells && cells < 12) continue;

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

	  
	if (expectHit && fdc_wire_expected_cell[cellNum] != NULL && cellNum < 25){
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs DOCA", distanceToWire, "Expected Hits", 100, 0 , 0.5);
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs Tracking FOM", thisTimeBasedTrack->FOM, "Expected Hits", 100, 0 , 1.0);
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs theta", theta_deg, "Expected Hits", 100, 0, 180);
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs phi", phi_deg, "Measured Hits", 100, -180, 180);
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs p", tmom, "Expected Hits", 100, 0 , 10.0);
	  Fill1DHistogram("FDC_Efficiency", "Offline", "Expected Hits Vs Hit Cells", cells, "Expected Hits", 25, -0.5 , 24.5);
	  
	  Double_t w, v;
	  if(fdc_wire_expected_cell[cellNum] != NULL && cellNum < 25){
	    // FILL HISTOGRAMS
	    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	    w = fdc_wire_expected_cell[cellNum]->GetBinContent(wireNum, 1) + 1.0;
	    fdc_wire_expected_cell[cellNum]->SetBinContent(wireNum, 1, w);
	    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	  }
	  
	  
	  // look in the presorted FDC Hits for a match
	  if(!locSortedFDCHits[cellNum][wireNum].empty() || !locSortedFDCHits[cellNum][wireNum-1].empty() || !locSortedFDCHits[cellNum][wireNum+1].empty()){
	    // Look not only in one wire, but also in adjacent ones (?)
	    // This can remove the dependence on the track error

	    if(!locSortedFDCHits[cellNum][wireNum].empty()){
	      const DFDCHit* locHit = *(locSortedFDCHits[cellNum][wireNum].begin());
	      // Fill the histograms only for the cell in the center
	      japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	      hWireTime_accepted->Fill(locHit->t);
	      hResVsT->Fill(distanceToWire, locHit->t);
	      japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	    }
	    	    
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs DOCA", distanceToWire, "Measured Hits", 100, 0 , 0.5);
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs Tracking FOM", thisTimeBasedTrack->FOM, "Measured Hits", 100, 0 , 1.0);
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs theta", theta_deg, "Measured Hits", 100, 0, 180);
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs phi", phi_deg, "Measured Hits", 100, -180, 180);
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs p", tmom, "Measured Hits", 100, 0 , 10.0);
	    Fill1DHistogram("FDC_Efficiency", "Offline", "Measured Hits Vs Hit Cells", cells, "Measured Hits", 25, -0.5 , 24.5);
	    
	    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	    v = fdc_wire_measured_cell[cellNum]->GetBinContent(wireNum, 1) + 1.0;
	    fdc_wire_measured_cell[cellNum]->SetBinContent(wireNum, 1, v);
	    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	  }

	  break; // break if 1 expected hit was found, 2 are geometrically not possible (speedup!)
	  
	} // End expected
      } // End wire loop
      
	// END WIRE ANALYSIS
      

      
	///////////// START PSEUDO ANALYSIS /////////////////////////////////////////
      
      
      if(fdc_pseudo_expected_cell[cellNum] != NULL && cellNum < 25){
	// FILL HISTOGRAMS
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	fdc_pseudo_expected_cell[cellNum]->Fill(interPosition.X(), interPosition.Y());
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
      }

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
	    double residualX = pseudoPosition.X() - interPosition2D.X();
	    double residualY = pseudoPosition.Y() - interPosition2D.Y();

	    // rotate to wire coordinate system
	    const DFDCWire * wire = locPseudo->wire;
	    double residualU = -1*(residual2D.Rotate(wire->angle - 0.0)).X();
	    double residualV = -1*(residual2D.Rotate(wire->angle - 0.0)).Y();

	    // these can be used for background studies/correction
	    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	    hPseudoRes->Fill(residualR);
	    hPseudoResX[cellNum]->Fill(residualX);
	    hPseudoResY[cellNum]->Fill(residualY);
	    hPseudoResU[cellNum]->Fill(residualU);
	    hPseudoResV[cellNum]->Fill(residualV);
	    unsigned int radius = interPosition2D.Mod()/(45/rad);
	    if (radius<rad)
	      hPseudoResUvsV[cellNum][radius]->Fill(residualU, residualV);
	    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	    
	    if (foundPseudo) continue; 
	    // to avoid double conting
	    	    
	    //if (residualR < 1.6){ // = 5 sigma in x and in y
	    if (residualR < 2){ // = to account for non-gaussian tails and tracking/extrapolation errors
	      foundPseudo = true;
	      
	      if(fdc_pseudo_measured_cell[cellNum] != NULL && cellNum < 25){
		// fill histogramms with the predicted, not with the measured position
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		fdc_pseudo_measured_cell[cellNum]->Fill(interPosition.X(), interPosition.Y());
		hPseudoTime_accepted->Fill(locPseudo->time);
		hPseudoResVsT->Fill(residualR, locPseudo->time);

		fdc_pseudo_residualU_cell[cellNum]->Fill(interPosition.X(), interPosition.Y(), residualU);
		fdc_pseudo_residualV_cell[cellNum]->Fill(interPosition.X(), interPosition.Y(), residualV);
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
jerror_t JEventProcessor_FDC_Efficiency::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FDC_Efficiency::fini(void)
{
  return NOERROR;
}

