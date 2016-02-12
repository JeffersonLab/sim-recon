// $Id$
//
//    File: JEventProcessor_CDC_Efficiency.cc
// Created: Tue Sep  9 15:41:38 EDT 2014
// Creator: hdcdcops (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_CDC_Efficiency.h"
using namespace jana;
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "HistogramTools.h"

static TH2D *cdc_measured_ring[29]; //Filled with total actually detected before division at end
static TH2D *cdc_measured_ring_DOCA0[29];
static TH2D *cdc_measured_ring_DOCA1[29];
static TH2D *cdc_measured_ring_DOCA2[29];
static TH2D *cdc_measured_ring_DOCA3[29];
static TH2D *cdc_measured_ring_DOCA4[29];
static TH2D *cdc_measured_ring_DOCA5[29];
static TH2D *cdc_measured_ring_DOCA6[29];
static TH2D *cdc_measured_ring_DOCA7[29];
static TH2D *cdc_expected_ring[29]; // Contains total number of expected hits by DOCA
static TH2D *cdc_expected_ring_DOCA0[29];
static TH2D *cdc_expected_ring_DOCA1[29];
static TH2D *cdc_expected_ring_DOCA2[29];
static TH2D *cdc_expected_ring_DOCA3[29];
static TH2D *cdc_expected_ring_DOCA4[29];
static TH2D *cdc_expected_ring_DOCA5[29];
static TH2D *cdc_expected_ring_DOCA6[29];
static TH2D *cdc_expected_ring_DOCA7[29];
static TH2I *ChargeVsTrackLength;
static TH1I * hChi2OverNDF;
static TH2I *hResVsT;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_Efficiency());
}
} // "C"


//------------------
// JEventProcessor_CDC_Efficiency (Constructor)
//------------------
JEventProcessor_CDC_Efficiency::JEventProcessor_CDC_Efficiency()
{
    ;
}

//------------------
// ~JEventProcessor_CDC_Efficiency (Destructor)
//------------------
JEventProcessor_CDC_Efficiency::~JEventProcessor_CDC_Efficiency()
{
    ;
}

//------------------
// init
//------------------
jerror_t JEventProcessor_CDC_Efficiency::init(void)
{
    // For the overall 2D plots, thus DOCA cut is used
    DOCACUT = 0.35;
    if(gPARMS){
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:DOCACUT", DOCACUT, "DOCA Cut on Efficiency Measurement");
    }  

    // Some information

    int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
    double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};

    // Define a different 2D histogram for each ring. X-axis is phi, Y-axis is radius (to plot correctly with "pol" option)
    japp->RootWriteLock();
    // create root folder for cdc and cd to it, store main dir
    TDirectory *main = gDirectory;
    gDirectory->mkdir("CDC_Efficiency")->cd();
    gDirectory->mkdir("CDC_View")->cd();

    for(int iring=0; iring<28; iring++){
        double r_start = radius[iring] - 0.8;
        double r_end = radius[iring] + 0.8;
        double phi_start = phi[iring]; // this is for center of straw. Need additional calculation for phi at end plate
        double phi_end = phi_start + TMath::TwoPi();

        char hname_measured[256];
        char hname_expected[256];
        sprintf(hname_measured, "cdc_measured_ring[%d]", iring+1);
        sprintf(hname_expected, "cdc_expected_ring[%d]", iring+1);
        if(gDirectory->Get(hname_measured) == NULL)
            cdc_measured_ring[iring+1] = new TH2D(hname_measured, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA0[iring+1] = new TH2D((TString)hname_measured +"DOCA0", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA1[iring+1] = new TH2D((TString)hname_measured +"DOCA1", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA2[iring+1] = new TH2D((TString)hname_measured +"DOCA2", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA3[iring+1] = new TH2D((TString)hname_measured +"DOCA3", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA4[iring+1] = new TH2D((TString)hname_measured +"DOCA4", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA5[iring+1] = new TH2D((TString)hname_measured +"DOCA5", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA6[iring+1] = new TH2D((TString)hname_measured +"DOCA6", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_measured_ring_DOCA7[iring+1] = new TH2D((TString)hname_measured +"DOCA7", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
        if(gDirectory->Get(hname_expected) == NULL)
            cdc_expected_ring[iring+1] = new TH2D(hname_expected, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA0[iring+1] = new TH2D((TString)hname_expected + "DOCA0", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA1[iring+1] = new TH2D((TString)hname_expected + "DOCA1", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA2[iring+1] = new TH2D((TString)hname_expected + "DOCA2", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA3[iring+1] = new TH2D((TString)hname_expected + "DOCA3", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA4[iring+1] = new TH2D((TString)hname_expected + "DOCA4", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA5[iring+1] = new TH2D((TString)hname_expected + "DOCA5", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA6[iring+1] = new TH2D((TString)hname_expected + "DOCA6", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
            cdc_expected_ring_DOCA7[iring+1] = new TH2D((TString)hname_expected + "DOCA7", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
    }	

    gDirectory->cd("/CDC_Efficiency");
    gDirectory->mkdir("Track_Quality")->cd();
    if(gDirectory->Get("hChi2OverNDF") == NULL)
        hChi2OverNDF = new TH1I("hChi2OverNDF","hChi2OverNDF", 500, 0.0, 1.0);
    if(gDirectory->Get("ChargeVsTrackLength") == NULL)
        ChargeVsTrackLength = new TH2I("ChargeVsTrackLength", "ChargeVsTrackLength", 1000, 0, 8.0, 2000, 0, 5000000);
    if(gDirectory->Get("hResVsT") == NULL)
        hResVsT = new TH2I("hResVsT","Tracking Residual (Biased) Vs Drift Time; Drift Time [ns]; Residual [cm]", 500, 0.0, 700.0, 1000, -0.5, 0.5);
    main->cd();

    japp->RootUnLock();


    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CDC_Efficiency::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);
    JCalibration *jcalib = dapp->GetJCalibration(runnumber);
    dgeom  = dapp->GetDGeometry(runnumber);
    //bfield = dapp->GetBfield();

    // Get the position of the CDC downstream endplate from DGeometry
    //double endplate_z,endplate_dz,endplate_rmin,endplate_rmax;
    //dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
    dgeom->GetCDCWires(cdcwires);
    unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
        135,135,146,146,158,158,170,170,182,182,197,197,
        209,209};

    // Get the straw sag parameters from the database
    vector< map<string, double> > tvals;
    max_sag.clear();
    sag_phi_offset.clear();
    unsigned int straw_count=0,ring_count=0;
    if (jcalib->Get("CDC/sag_parameters", tvals)==false){
        vector<double>temp,temp2;
        for(unsigned int i=0; i<tvals.size(); i++){
            map<string, double> &row = tvals[i];

            temp.push_back(row["offset"]);
            temp2.push_back(row["phi"]);

            straw_count++;
            if (straw_count==numstraws[ring_count]){
                max_sag.push_back(temp);
                sag_phi_offset.push_back(temp2);
                temp.clear();
                temp2.clear();
                straw_count=0;
                ring_count++;
            }
        }
    }

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_CDC_Efficiency::evnt(JEventLoop *loop, uint64_t eventnumber){

    vector< const DCDCHit *> locCDCHitVector;
    loop->Get(locCDCHitVector);

    const DDetectorMatches *detMatches;
    vector<DSCHitMatchParams> SCMatches;
    if(!dIsNoFieldFlag){
        loop->GetSingle(detMatches);
    }

    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Require Single track events
        //if (trackCandidateVector.size() != 1) return NOERROR;
        //const DTrackCandidate* thisTrackCandidate = trackCandidateVector[0];
        // Cut very loosely on the track quality
        const DTrackTimeBased *thisTimeBasedTrack;
        bestHypothesis->GetSingle(thisTimeBasedTrack);
        hChi2OverNDF->Fill(thisTimeBasedTrack->FOM);
        if (thisTimeBasedTrack->FOM < 1E-20) return NOERROR;
        //The cuts used for track quality
        if(!dIsNoFieldFlag){ // Quality cuts for Field on runs.
            if(thisTimeBasedTrack->ddEdx_CDC > 1E-3) {
                //cout << "Cut on dEdX" << endl;
                continue; // Trying to cut out "proton" candidates
            }
            if(thisTimeBasedTrack->pmag() < 0.5 || thisTimeBasedTrack->pmag() > 6.0) {
                //cout << "Cut on momentum" << endl;
                continue; // Cut on the reconstructed momentum to make sure we have the right
            }
            if(!detMatches->Get_SCMatchParams(thisTimeBasedTrack, SCMatches)) {
                //cout << "Cut on detector matches" << endl;
                continue; // Require there to be at least one match to the Start Counter
            }
            if(fabs(thisTimeBasedTrack->position().Z() - 65.0) > 15.0 || thisTimeBasedTrack->position().Perp() > 1.0) {
                //cout << " Cut on vertex " << endl;
                continue; // Cut on reconstructed vertex location
            }
            //Let's try a cut on the angle of the track here
            //if (TMath::Abs(((*trackIter)->momentum().Theta() * TMath::RadToDeg() ) - 90 ) > 5) continue;
        }

        // There are two loops over the CDC hits, first to impose some additional cuts, then to calculate the efficiency
        vector< int > ringsHit;
        bool hasRing1 = false, hasRing2 = false, hasRing27 = false, hasRing28 = false;

        // Loop over the pulls to get the appropriate information for our ring
        vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
        for (unsigned int i = 0; i < pulls.size(); i++){
            const DCDCTrackHit * thisTrackHit = pulls[i].cdc_hit;
            if (thisTrackHit == NULL) continue;
            //cout << "Checking ring number" << thisTrackHit->wire->ring << endl;
            if ( find(ringsHit.begin(), ringsHit.end(), thisTrackHit->wire->ring) == ringsHit.end()) ringsHit.push_back(thisTrackHit->wire->ring);
            if (thisTrackHit->wire->ring == 1) hasRing1 = true;
            else if (thisTrackHit->wire->ring == 2) hasRing2 = true;
            else if (thisTrackHit->wire->ring == 27) hasRing27 = true;
            else if (thisTrackHit->wire->ring == 28) hasRing28 = true;
        }
        if(!dIsNoFieldFlag){
            if ( !(hasRing1 || hasRing2 ) || !(hasRing27 || hasRing28)) continue; // Has a hit in one of the two inner layers and one of the two outer
        }
        if (ringsHit.size() < 15) continue; //At least half of the rings hit
        //cout << "Passed Number of hits cuts" << endl;

        // Alright now we truly have the tracks we are interested in for calculating the efficiency
        // git-r-dun

        for (unsigned int ringIndex = 0; ringIndex < cdcwires.size(); ringIndex ++){
            int ringNum = ringIndex +1;
            vector< DCDCWire * > wireByNumber = cdcwires[ringIndex];
            for (unsigned int wireIndex = 0; wireIndex < wireByNumber.size(); wireIndex++){
                int wireNum = wireIndex+1;
                DCDCWire * wire = wireByNumber[wireIndex]; 
                double wireLength = wire->L;
                double distanceToWire = thisTimeBasedTrack->rt->DistToRT(wire, &wireLength);
                bool expectHit = false;
                double delta = 0.0;
                double dz = 0.0;
                if (distanceToWire < 1.2 ) {
                    // Loose cut before delta information
                    // Need to get phi_doca for each of the wires that pass this cut
                    DVector3 pos, mom;
                    thisTimeBasedTrack->rt->GetLastDOCAPoint(pos, mom);
                    // Form the vector between the wire and the DOCA point
                    DVector3 DOCA = (-1) * ((wire->origin - pos) - (wire->origin - pos).Dot(wire->udir) * wire->udir);

                    double docaphi = DOCA.Phi();
                    dz = (pos - wire->origin).Z();
                    //cout << "distanceToWire = " << distanceToWire << " DOCA = " << DOCA.Mag() << endl;
                    // Get delta at this location for this straw
                    int ring_index = wire->ring - 1;
                    int straw_index = wire->straw - 1;
                    delta = max_sag[ring_index][straw_index] * ( 1. - (dz*dz/5625.)) * TMath::Cos(docaphi + sag_phi_offset[ring_index][straw_index]);
                    if (distanceToWire < (0.78 + delta) && fabs(dz) < 65.0) expectHit = true;
                }

                if (expectHit && cdc_expected_ring[ringNum] != NULL && ringNum < 29){
                    double dx = thisTimeBasedTrack->rt->Straw_dx(wire, 0.78);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs Path Length",
                            dx,
                            "Expected Hits", 
                            100, 0 , 4.0);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs DOCA",
                            distanceToWire,
                            "Expected Hits",
                            100, 0 , 0.78);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs Tracking FOM",
                            thisTimeBasedTrack->FOM,
                            "Expected Hits",
                            100, 0 , 1.0);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs theta",
                            thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(),
                            "Expected Hits",
                            100, 0, 180);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs p",
                            thisTimeBasedTrack->pmag(),
                            "Expected Hits",
                            100, 0 , 4.0);
                    Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs delta",
                            delta,
                            "Expected Hits",
                            100, -0.3 , 0.3); 
                    Fill2DHistogram("CDC_Efficiency", "Offline", "Expected hits p Vs Theta",
                            thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), thisTimeBasedTrack->pmag(),
                            "Expected Hits",
                            100, 0, 180, 100, 0 , 4.0);

                    bool foundHit = false;
                    // loop over the CDC Hits to look for a match
                    for( unsigned int hitNum = 0; hitNum < locCDCHitVector.size(); hitNum++){
                        const DCDCHit * locHit = locCDCHitVector[hitNum];
                        if(locHit->ring == ringNum && locHit->straw == wireNum){
                            foundHit = true;
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs Path Length",
                                    dx,
                                    "Measured Hits",
                                    100, 0 , 4.0);
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs DOCA",
                                    distanceToWire,
                                    "Measured Hits",
                                    100, 0 , 0.78);
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs Tracking FOM",
                                    thisTimeBasedTrack->FOM,
                                    "Measured Hits",
                                    100, 0 , 1.0);
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs theta",
                                    thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(),
                                    "Measured Hits",
                                    100, 0, 180);
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs p",
                                    thisTimeBasedTrack->pmag(),
                                    "Measured Hits",
                                    100, 0 , 4.0);
                            Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs delta",
                                    delta,
                                    "Measured Hits",
                                    100, -0.3 , 0.3);
                            Fill2DHistogram("CDC_Efficiency", "Offline", "Measured hits p Vs Theta",
                                    thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), thisTimeBasedTrack->pmag(),
                                    "Measured Hits",
                                    100, 0, 180, 100, 0 , 4.0);
                            break;
                        }

                    }
                    if (foundHit){
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Path Length",
                                dx,1.0,
                                "Efficiency; dx [cm]; Efficiency",
                                100, 0 , 4.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs DOCA",
                                distanceToWire,1.0,
                                "Efficiency; DOCA [cm]; Efficiency",
                                100, 0 , 0.78);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Tracking FOM",
                                thisTimeBasedTrack->FOM,1.0,
                                "Efficiency; Tracking FOM; Efficiency",
                                100, 0 , 1.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs theta",
                                thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(),1.0,
                                "Efficiency; Track #Theta [deg]; Efficiency",
                                100, 0, 180);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs p",
                                thisTimeBasedTrack->pmag(),1.0,
                                "Efficiency; Momentum [GeV]; Efficiency",
                                100, 0 , 4.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs delta",
                                delta,1.0,
                                "Efficiency; #delta [cm]; Efficiency",
                                100, -0.3 , 0.3);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency p Vs Theta",
                                thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), thisTimeBasedTrack->pmag(),1.0,
                                "Efficiency; Track #Theta [deg]; Momentum [GeV]",
                                100, 0, 180, 100, 0 , 4.0);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency distance Vs delta",
                                delta,distanceToWire,1.0,
                                "Efficiency;#delta [cm]; DOCA [cm]",
                                100, -0.3, 0.3, 100, 0 , 1.2);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency z Vs delta",
                                delta,dz,1.0,
                                "Efficiency;#delta [cm]; z [cm] (CDC local coordinates)",
                                100, -0.3, 0.3, 150, -75 , 75);
                    }
                    else{
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Path Length",
                                dx,0.0,
                                "Efficiency; dx [cm]; Efficiency",
                                100, 0 , 4.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs DOCA",
                                distanceToWire,0.0,
                                "Efficiency; DOCA [cm]; Efficiency",
                                100, 0 , 0.78);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Tracking FOM",
                                thisTimeBasedTrack->FOM,0.0,
                                "Efficiency; Tracking FOM; Efficiency",
                                100, 0 , 1.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs theta",
                                thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(),0.0,
                                "Efficiency; Track #Theta [deg]; Efficiency",
                                100, 0, 180);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs p",
                                thisTimeBasedTrack->pmag(),0.0,
                                "Efficiency; Momentum [GeV]; Efficiency",
                                100, 0 , 4.0);
                        Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs delta",
                                delta,0.0,
                                "Efficiency; #delta [cm]; Efficiency",
                                100, -0.3 , 0.3);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency p Vs Theta",
                                thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), thisTimeBasedTrack->pmag(),0.0,
                                "Efficiency; Track #Theta [deg]; Momentum [GeV]",
                                100, 0, 180, 100, 0 , 4.0);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency distance Vs delta",
                                delta,distanceToWire,0.0, 
                                "Efficiency;#delta [cm]; DOCA [cm]",
                                100, -0.3, 0.3, 100, 0 , 1.2);
                        Fill2DProfile("CDC_Efficiency", "Online", "Efficiency z Vs delta",
                                delta,dz,0.0,
                                "Efficiency;#delta [cm]; z [cm] (CDC local coordinates)",
                                100, -0.3, 0.3, 150, -75 , 75);
                    }
                }

                if (distanceToWire < 0.78 && expectHit){
                    //Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
                    //cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
                    //Fill the expected number of hits histogram
                    Double_t w, v;
                    if(cdc_expected_ring[ringNum] != NULL && ringNum < 29){
                        if (distanceToWire < DOCACUT){
                            japp->RootWriteLock();
                            w = cdc_expected_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                            cdc_expected_ring[ringNum]->SetBinContent(wireNum, 1, w);
                            japp->RootUnLock();
                        }
                        switch ( (int) (distanceToWire * 10) % 8){
                            case 0:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA0[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA0[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 1:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA1[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA1[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 2:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA2[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA2[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 3:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA3[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA3[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 4:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA4[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA4[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 5:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA5[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA5[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 6:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA6[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA6[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            case 7:
                                japp->RootWriteLock();
                                w = cdc_expected_ring_DOCA7[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_expected_ring_DOCA7[ringNum]->SetBinContent(wireNum, 1, w);
                                japp->RootUnLock();
                                break;
                            default:
                                cout << "Unknown DOCA Lookup?" << endl;

                        }
                    }

                    // loop over the CDC Hits to look for a match
                    for( unsigned int hitNum = 0; hitNum < locCDCHitVector.size(); hitNum++){
                        const DCDCHit * locHit = locCDCHitVector[hitNum];
                        if(locHit->ring == ringNum && locHit->straw == wireNum){
                            if (distanceToWire < DOCACUT){
                                japp->RootWriteLock();
                                //printf("Matching Hit!!!!!\n");
                                v = cdc_measured_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_measured_ring[ringNum]->SetBinContent(wireNum, 1, v);
                                double dx = thisTimeBasedTrack->rt->Straw_dx(wire, 0.78);
                                ChargeVsTrackLength->Fill(dx,locHit->q);
                                japp->RootUnLock();
                            }
                            switch ( (int) (distanceToWire * 10) % 8){
                                case 0:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA0[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA0[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 1:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA1[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA1[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 2:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA2[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA2[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 3:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA3[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA3[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 4:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA4[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA4[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 5:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA5[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA5[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 6:
                                    japp->RootWriteLock();
                                    v = cdc_measured_ring_DOCA6[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA6[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                case 7:
                                    japp->RootWriteLock();
                                    v  = cdc_measured_ring_DOCA7[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                    cdc_measured_ring_DOCA7[ringNum]->SetBinContent(wireNum, 1, v);
                                    japp->RootUnLock();
                                    break;
                                default:
                                    cout << "Unknown DOCA Lookup?" << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_CDC_Efficiency::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_CDC_Efficiency::fini(void)
{
    return NOERROR;
}

