// $Id$
//
//    File: JEventProcessor_cdc_efficiency.cc
// Created: Tue Sep  9 15:41:38 EDT 2014
// Creator: hdcdcops (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_cdc_efficiency.h"
using namespace jana;
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
    app->AddProcessor(new JEventProcessor_cdc_efficiency());
}
} // "C"


//------------------
// JEventProcessor_cdc_efficiency (Constructor)
//------------------
JEventProcessor_cdc_efficiency::JEventProcessor_cdc_efficiency()
{
    ;
}

//------------------
// ~JEventProcessor_cdc_efficiency (Destructor)
//------------------
JEventProcessor_cdc_efficiency::~JEventProcessor_cdc_efficiency()
{
    ;
}

//------------------
// init
//------------------
jerror_t JEventProcessor_cdc_efficiency::init(void)
{
    // This is called once at program startup. If you are creating
    // and filling historgrams in this plugin, you should lock the
    // ROOT mutex like this:
    BFIELD = 0;
    COSMICS = 0;
    DOCACUT = 0.35;
    USE_TIMEBASEDTRACKS = 1;

    REQUIRE_BEAM = 1;
    BEAM_CURRENT = 50;
    BEAM_EVENTS_TO_KEEP = 99999999;

    if(gPARMS){
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:BFIELD", BFIELD, "Set to > 0 to use DTrackTimeBased instead of DTrackCandidate_StraightLine");
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:COSMICS", COSMICS, "Set to 1 to assume cosmic tracks");
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:DOCACUT", DOCACUT, "DOCA Cut on Efficiency Measurement");
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:USE_TIMEBASEDTRACKS", USE_TIMEBASEDTRACKS, "Set to 0 to use the other track types");
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:REQUIRE_BEAM", REQUIRE_BEAM);
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:BEAM_EVENTS_TO_KEEP", BEAM_EVENTS_TO_KEEP);
    }  

    if (COSMICS > 0) REQUIRE_BEAM = 0;

    // Some information

    int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
    double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};

    // Define a different 2D histogram for each ring. X-axis is phi, Y-axis is radius (to plot correctly with "pol" option)
    japp->RootWriteLock();
    // create root folder for cdc and cd to it, store main dir
    TDirectory *main = gDirectory;
    gDirectory->mkdir("cdc_efficiency")->cd();

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
            cdc_measured_ring_DOCA0[iring+1] = new TH2D((TString) hname_measured +"DOCA0", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
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
jerror_t JEventProcessor_cdc_efficiency::brun(JEventLoop *eventLoop, int runnumber)
{
    // This is called whenever the run number changes
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    dgeom  = dapp->GetDGeometry(runnumber);
    //bfield = dapp->GetBfield();

    // Get the position of the CDC downstream endplate from DGeometry
    //double endplate_z,endplate_dz,endplate_rmin,endplate_rmax;
    //dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
    dgeom->GetCDCWires (cdcwires);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_cdc_efficiency::evnt(JEventLoop *loop, int eventnumber)
{
    if (REQUIRE_BEAM){
        vector<const DEPICSvalue *> epicsValues;
        loop->Get(epicsValues);
        for(unsigned int j = 0; j < epicsValues.size(); j++){
            const DEPICSvalue *thisValue = epicsValues[j];
            if (strcmp((thisValue->name).c_str(), "IBCAD00CRCUR6") == 0){
                BEAM_CURRENT = thisValue->fval;
                Fill1DHistogram("HLDetectorTiming", "", "Beam Current",
                        BEAM_CURRENT,
                        "Beam Current; Beam Current [nA]; Entries",
                        100, 0, 200);
            }
            //cout << "EPICS Name " <<  (thisValue->name).c_str() << " Value " << thisValue->fval << endl;
        }
        if (BEAM_CURRENT < 10.0) {
            return NOERROR; // Skip events where we can't verify the beam current
        }
    }

    if (USE_TIMEBASEDTRACKS != 0 ) {
        TimeBasedAnalysis(loop);
        return NOERROR;
    }
    //The idea is to treat the cases where there is BField and no BField seperately. 
    // This idea is no longer necessary...
    vector<const DTrackCandidate*> locTrackCandidates;

    if (BFIELD == 0 && COSMICS ==0){
        loop->Get(locTrackCandidates, "StraightLine");
    }
    else loop->Get(locTrackCandidates, "CDC");

    vector<const DCDCHit*> locCDCHits;
    loop->Get(locCDCHits);

    vector<const DCDCTrackHit *> locCDCTrackHitVector;
    loop->Get(locCDCTrackHitVector);

    vector<const DTrackWireBased *> wireBasedTracks;
    loop->Get(wireBasedTracks);

    for( unsigned int i = 0; i<locTrackCandidates.size(); i++){
        //if (locTrackCandidates.size() != 1) continue;
        // We will only do this calculation if there is a hit in the innermost and outermost rings 1 and 28. 
        // In order to include these layers we will also use tracks that have a hit in both rings 2 and 27.
        const DTrackCandidate * locTrack = locTrackCandidates[i];

        bool hasRing1 = false, hasRing2 = false, hasRing27 = false, hasRing28 = false;
        vector< int > usedCDCTrackHitVector = locTrack->used_cdc_indexes;
        vector< int > ringsHit;
        //cout << " StartTrack ============ " << (locTrack->used_cdc_indexes).size() << " CDC hits" << endl;
        for ( vector< int >::const_iterator index = usedCDCTrackHitVector.begin(); index != usedCDCTrackHitVector.end(); index++){
            const DCDCTrackHit * thisTrackHit = locCDCTrackHitVector[(*index)];
            //cout << "Checking ring number" << thisTrackHit->wire->ring << endl;
            if ( find(ringsHit.begin(), ringsHit.end(), thisTrackHit->wire->ring) == ringsHit.end()) ringsHit.push_back(thisTrackHit->wire->ring);  
            if (thisTrackHit->wire->ring == 1) hasRing1 = true;
            else if (thisTrackHit->wire->ring == 2) hasRing2 = true;
            else if (thisTrackHit->wire->ring == 27) hasRing27 = true;
            else if (thisTrackHit->wire->ring == 28) hasRing28 = true; 
        }
        if ( !(hasRing1 || hasRing2 ) || !(hasRing27 || hasRing28)) continue;
        if (ringsHit.size() < 15) continue; //At least half of the rings hit

        const DTrackWireBased *thisWireBasedTrack = NULL;
        for (unsigned int iTrack=0 ; iTrack < wireBasedTracks.size(); iTrack++){
            if( wireBasedTracks[iTrack]->candidateid == i) {
                thisWireBasedTrack = wireBasedTracks[iTrack];
                //cout << "Found the time based track " << endl;
                break;
            }
        }

        float chisq = locTrack->chisq;
        int Ndof = locTrack->Ndof;
        //printf("chisq = %f, Ndof = %i \n", chisq, Ndof);
        hChi2OverNDF->Fill(chisq / Ndof);
        //const DReferenceTrajectory *rt = locTrack->rt;
        //		if (rt == NULL) printf("rt null \n");
        if(chisq / Ndof > 5) continue;	
        DVector3 trackPosition = locTrack->position();
        DVector3 trackMomentum = locTrack->momentum();
        //printf("New Good Track\n");
        //vector< vector< DCDCWire * > > cdcwires
        for (unsigned int ringIndex = 0; ringIndex < cdcwires.size(); ringIndex ++){
            int ringNum = ringIndex +1;
            vector< DCDCWire * > wireByNumber = cdcwires[ringIndex];
            for (unsigned int wireIndex = 0; wireIndex < wireByNumber.size(); wireIndex++){
                int wireNum = wireIndex+1;
                DCDCWire * wire = wireByNumber[wireIndex];
                //float DOCACut = 0.6; // cm
                if (BFIELD == 0 && COSMICS ==0){					
                    DVector3 wirePosition = wire->origin;
                    DVector3 wireDirection = wire->udir;
                    //wirePosition.Print(); wireDirection.Print();
                    Float_t a = trackMomentum.Dot(trackMomentum);
                    Float_t b = trackMomentum.Dot(wireDirection);
                    Float_t c = wireDirection.Dot(wireDirection);
                    DVector3 w0 = trackPosition - wirePosition;
                    Float_t d = trackMomentum.Dot(w0);
                    Float_t e = wireDirection.Dot(w0);
                    Float_t sc = ((b*e - c*d)/(a*c-b*b));
                    if (sc < 0) continue; // Track must come from location away from origin
                    DVector3 POCAOnTrack = trackPosition + sc * trackMomentum;
                    DVector3 LOCA = w0 + ((b*e - c*d)/(a*c-b*b))*trackMomentum - ((a*e - b*d)/(a*c-b*b))*wireDirection;
                    if (LOCA.Z() > 50.0 || LOCA.Z() < -50.0) continue; // Skip if outside the cdc
                    Float_t DOCA = LOCA.Mag();
                    if (DOCA < DOCACUT) {
                        japp->RootWriteLock();
                        //Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
                        //cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
                        if(cdc_expected_ring[ringNum] != NULL){
                            Double_t w = cdc_expected_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                            cdc_expected_ring[ringNum]->SetBinContent(wireNum, 1, w);
                        }

                        //printf("DOCA = %f\n", DOCA);
                        for( unsigned int hitNum = 0; hitNum < locCDCHits.size(); hitNum++){
                            const DCDCHit * locHit = locCDCHits[hitNum];
                            if(locHit->ring == ringNum && locHit->straw == wireNum){
                                //printf("Matching Hit!!!!!\n");
                                Double_t v = cdc_measured_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_measured_ring[ringNum]->SetBinContent(wireNum, 1, v);
                                // Now we want to calculate the track length through the CDC straw
                                // From this we will plot the track length vs charge (Integral)
                                DVector3 trackMomentumInWireFrame = trackMomentum;
                                DVector3 POCAOnTrackInWireFrame = POCAOnTrack;
                                //wire->FromLab(trackMomentumInWireFrame);
                                trackMomentumInWireFrame.SetXYZ(trackMomentumInWireFrame.Dot(wire->sdir), trackMomentumInWireFrame.Dot(wire->tdir), trackMomentumInWireFrame.Dot(wire->udir));
                                wire->FromLab(POCAOnTrackInWireFrame);
                                // Now to get our two points where the track intersects the edge of the straw
                                DVector3 intersect1;
                                DVector3 intersect2;
                                Float_t strawRadius = 0.8;

                                Float_t a = trackMomentumInWireFrame.X() * trackMomentumInWireFrame.X() + trackMomentumInWireFrame.Y() * trackMomentumInWireFrame.Y();
                                Float_t b = 2.0 * (POCAOnTrackInWireFrame.X() * trackMomentumInWireFrame.X() + POCAOnTrackInWireFrame.Y() * trackMomentumInWireFrame.Y());
                                Float_t c = POCAOnTrackInWireFrame.X() * POCAOnTrackInWireFrame.X() + POCAOnTrackInWireFrame.Y() * POCAOnTrackInWireFrame.Y() - strawRadius;

                                //Just use quadratic equation to solve for the two points
                                if ((b * b - 4.0 * a * c) < 0) continue; //Result imaginary... It shouldnt be if we made it here
                                Float_t t1 = ((-1.0) * b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
                                Float_t t2 = ((-1.0) * b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

                                DVector3 p1 = POCAOnTrackInWireFrame + t1 * trackMomentumInWireFrame;
                                DVector3 p2 = POCAOnTrackInWireFrame + t2 * trackMomentumInWireFrame;

                                Float_t trackLength = (p1 - p2).Mag();
                                ChargeVsTrackLength->Fill(trackLength,locHit->q);							
                                //printf("DOCA = %f , Track Length = %f , Theta(Lab) = %f , Theta(Straw) = %f \n", DOCA, trackLength, trackMomentum.Theta(), trackMomentumInWireFrame.Theta());

                            }
                        } 
                        japp->RootUnLock();
                    }
                }
                else{
                    if (thisWireBasedTrack == NULL) {
                        continue;
                    }
                    if (TMath::Prob(thisWireBasedTrack->chisq, thisWireBasedTrack->Ndof) < 0.01){ //1% cut on tracking FOM
                        continue;
                    }
                    if(fabs(thisWireBasedTrack->position().Z() - 65.0) > 3.0 || thisWireBasedTrack->position().Perp() > 1.0){ // Cut on vertex used for track
                        continue;
                    }
                    //rt->FastSwim(thisTimeBasedTrack->position(), thisTimeBasedTrack->momentum(), thisTimeBasedTrack->charge());

                    if (thisWireBasedTrack->rt->Nswim_steps == 0) {
                        continue;
                    }
                    double wireLength = 50;
                    double distanceToWire = thisWireBasedTrack->rt->DistToRT(wire, &wireLength);
                    if (distanceToWire < DOCACUT){
                        //Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
                        //cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
                        if(cdc_expected_ring[ringNum] != NULL){
                            japp->RootWriteLock();
                            Double_t w = cdc_expected_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                            cdc_expected_ring[ringNum]->SetBinContent(wireNum, 1, w);
                            japp->RootUnLock();
                        }
                        else {
                            continue;
                        }

                        for( unsigned int hitNum = 0; hitNum < locCDCHits.size(); hitNum++){
                            const DCDCHit * locHit = locCDCHits[hitNum];
                            if(locHit->ring == ringNum && locHit->straw == wireNum){
                                japp->RootWriteLock();
                                //printf("Matching Hit!!!!!\n");
                                Double_t v = cdc_measured_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
                                cdc_measured_ring[ringNum]->SetBinContent(wireNum, 1, v);
                                double dx = thisWireBasedTrack->rt->Straw_dx(wire, 0.8);
                                ChargeVsTrackLength->Fill(dx,locHit->q);
                                japp->RootUnLock();
                            }
                        }

                    }
                }
            }
        }
        //printf("==================");
        //trackPosition.Print(); trackMomentum.Print(); 
        //printf("==================");		
        //double z = position.Z();
        //printf("z = %f\n", z);
    }

    //
    // japp->RootWriteLock();
    //  ... fill historgrams or trees ...
    // japp->RootUnLock();


    return NOERROR;
}

void JEventProcessor_cdc_efficiency::TimeBasedAnalysis(JEventLoop *loop){

    vector< const DCDCHit *> locCDCHitVector;
    loop->Get(locCDCHitVector);

    vector< const DCDCTrackHit *> locCDCTrackHitVector; // Get from the track

    const DDetectorMatches *detMatches;
    loop->GetSingle(detMatches);
    vector<DSCHitMatchParams> SCMatches;

    vector< const DTrackTimeBased *> locTrackTimeBasedVector;
    loop->Get(locTrackTimeBasedVector);

    //Loop over the DTrackCandidates and get the CDCTrackHits associated witht the tracks
    vector< const DTrackTimeBased *>::const_iterator trackIter;
    for (trackIter = locTrackTimeBasedVector.begin(); trackIter != locTrackTimeBasedVector.end(); trackIter++){
        hChi2OverNDF->Fill(TMath::Prob((*trackIter)->chisq, (*trackIter)->Ndof));
        //The cuts used for track quality
        if(TMath::Prob((*trackIter)->chisq, (*trackIter)->Ndof) < 0.01) continue; //Traking FOM cut
        if((*trackIter)->ddEdx_CDC > 1E-3) {
            //cout << "Cut on dEdX" << endl;
            continue; // Trying to cut out "proton" candidates
        }
        if((*trackIter)->pmag() < 0.3 || (*trackIter)->pmag() > 6.0) {
            //cout << "Cut on momentum" << endl;
            continue; // Cut on the reconstructed momentum to make sure we have the right
        }
        if(!detMatches->Get_SCMatchParams((*trackIter), SCMatches)) {
            //cout << "Cut on detector matches" << endl;
            continue; // Require there to be at least one match to the Start Counter
        }
        if(fabs((*trackIter)->position().Z() - 65.0) > 3.0 || (*trackIter)->position().Perp() > 1.0) {
            //cout << " Cut on vertex " << endl;
            continue; // Cut on reconstructed vertex location
        }

        //cout << "Passed Track Quality Cuts" << endl;
        // Let's try a cut on the angle of the track here
        //if (TMath::Abs(((*trackIter)->momentum().Theta() * TMath::RadToDeg() ) - 90 ) > 5) continue;

        // Get the CDCTrackHits Associated with the track
        (*trackIter)->Get(locCDCTrackHitVector);
        //cout << "There are " << locCDCTrackHitVector.size() << " CDCTrackHits associated with this track" << endl;
        // There are two loops over the CDC hits, first to impose some additional cuts, then to calculate the efficiency
        vector< int > ringsHit;
        bool hasRing1 = false, hasRing2 = false, hasRing27 = false, hasRing28 = false;

        vector<const DCDCTrackHit *>::const_iterator iHit;
        for (iHit = locCDCTrackHitVector.begin(); iHit != locCDCTrackHitVector.end(); iHit++){
            const DCDCTrackHit * thisTrackHit = (*iHit);
            //cout << "Checking ring number" << thisTrackHit->wire->ring << endl;
            if ( find(ringsHit.begin(), ringsHit.end(), thisTrackHit->wire->ring) == ringsHit.end()) ringsHit.push_back(thisTrackHit->wire->ring);
            if (thisTrackHit->wire->ring == 1) hasRing1 = true;
            else if (thisTrackHit->wire->ring == 2) hasRing2 = true;
            else if (thisTrackHit->wire->ring == 27) hasRing27 = true;
            else if (thisTrackHit->wire->ring == 28) hasRing28 = true;
        }

        if ( !(hasRing1 || hasRing2 ) || !(hasRing27 || hasRing28)) continue; // Has a hit in one of the two inner layers and one of the two outer
        if (ringsHit.size() < 15) continue; //At least half of the rings hit
        //cout << "Passed Number of hits cuts" << endl;

        // Alright now we truly have the tracks we are interested in for calculating the efficiency
        // git-r-dun
        //rt->FastSwim(thisTimeBasedTrack->position(), thisTimeBasedTrack->momentum(), thisTimeBasedTrack->charge());

        for (unsigned int ringIndex = 0; ringIndex < cdcwires.size(); ringIndex ++){
            int ringNum = ringIndex +1;
            vector< DCDCWire * > wireByNumber = cdcwires[ringIndex];
            for (unsigned int wireIndex = 0; wireIndex < wireByNumber.size(); wireIndex++){
                int wireNum = wireIndex+1;
                DCDCWire * wire = wireByNumber[wireIndex]; 

                double wireLength = 50;
                double distanceToWire = (*trackIter)->rt->DistToRT(wire, &wireLength);
                //Loop over the track hits. If there is a match for ring and straw, plot residual vs drift time.
                for(unsigned int iTrackHit = 0; iTrackHit < locCDCTrackHitVector.size(); iTrackHit++){
                    double residual;
                    const DCDCTrackHit * thisTrackHit = locCDCTrackHitVector[iTrackHit];
                    if (thisTrackHit->wire->ring == ringNum && thisTrackHit->wire->straw == wireNum){
                        residual = distanceToWire - thisTrackHit->dist;
                        hResVsT->Fill(thisTrackHit->tdrift, residual);
                    }
                }
                if (distanceToWire < 0.78 && cdc_expected_ring[ringNum] != NULL && ringNum < 29){
                    double dx = (*trackIter)->rt->Straw_dx(wire, 0.78);
                    Fill1DHistogram("cdc_efficiency", "", "Expected Hits Vs Path Length",
                            dx,
                            "Expected Hits", 
                            100, 0 , 4.0);
                    Fill1DHistogram("cdc_efficiency", "", "Expected Hits Vs DOCA",
                            distanceToWire,
                            "Expected Hits",
                            100, 0 , 0.78);
                    Fill1DHistogram("cdc_efficiency", "", "Expected Hits Vs Tracking FOM",
                            TMath::Prob((*trackIter)->chisq, (*trackIter)->Ndof),
                            "Expected Hits",
                            100, 0 , 1.0);
                    Fill1DHistogram("cdc_efficiency", "", "Expected Hits Vs theta",
                            (*trackIter)->momentum().Theta()*TMath::RadToDeg(),
                            "Expected Hits",
                            100, 0, 180);
                    Fill1DHistogram("cdc_efficiency", "", "Expected Hits Vs p",
                            (*trackIter)->pmag(),
                            "Expected Hits",
                            100, 0 , 4.0); 
                    Fill2DHistogram("cdc_efficiency", "", "Expected hits p Vs Theta",
                            (*trackIter)->momentum().Theta()*TMath::RadToDeg(), (*trackIter)->pmag(),
                            "Expected Hits",
                            100, 0, 180, 100, 0 , 4.0);

                    // loop over the CDC Hits to look for a match
                    for( unsigned int hitNum = 0; hitNum < locCDCHitVector.size(); hitNum++){
                        const DCDCHit * locHit = locCDCHitVector[hitNum];
                        if(locHit->ring == ringNum && locHit->straw == wireNum){
                            Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs Path Length",
                                    dx,
                                    "Measured Hits",
                                    100, 0 , 4.0);
                            Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs DOCA",
                                    distanceToWire,
                                    "Measured Hits",
                                    100, 0 , 0.78);
                            Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs Tracking FOM",
                                    TMath::Prob((*trackIter)->chisq, (*trackIter)->Ndof),
                                    "Measured Hits",
                                    100, 0 , 1.0);
                            Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs theta",
                                    (*trackIter)->momentum().Theta()*TMath::RadToDeg(),
                                    "Measured Hits",
                                    100, 0, 180);
                            Fill1DHistogram("cdc_efficiency", "", "Measured Hits Vs p",
                                    (*trackIter)->pmag(),
                                    "Measured Hits",
                                    100, 0 , 4.0);
                            Fill2DHistogram("cdc_efficiency", "", "Measured hits p Vs Theta",
                                    (*trackIter)->momentum().Theta()*TMath::RadToDeg(), (*trackIter)->pmag(),
                                    "Measured Hits",
                                    100, 0, 180, 100, 0 , 4.0);

                        }
                    }
                }
                if (distanceToWire < 0.78){
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
                                double dx = (*trackIter)->rt->Straw_dx(wire, 0.78);
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
    return;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_cdc_efficiency::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_cdc_efficiency::fini(void)
{
    return NOERROR;
}

