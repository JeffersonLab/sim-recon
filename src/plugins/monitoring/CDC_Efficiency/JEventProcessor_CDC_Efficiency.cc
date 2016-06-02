// $Id$
//
//    File: JEventProcessor_CDC_Efficiency.cc
// Created: Tue Sep  9 15:41:38 EDT 2014
// Creator: hdcdcops (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_CDC_Efficiency.h"
using namespace jana;
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125CDCPulse.h"
#include "HistogramTools.h"


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
	dMinTrackingFOM = 5.73303E-7; // +/- 5 sigma
	dMinNumRingsToEvalSuperlayer = 3; //technically, reconstruction minimum is 2 (TRKFIND:MIN_SEED_HITS) //but trust more with 3

    // For the overall 2D plots, thus DOCA cut is used
    DOCACUT = 0.35;
    if(gPARMS){
        gPARMS->SetDefaultParameter("CDC_EFFICIENCY:DOCACUT", DOCACUT, "DOCA Cut on Efficiency Measurement");
    }  
    
    for (unsigned int i=0; i < 28; i++){
        for (unsigned int j=0; j < 209; j++){
            ChannelFromRingStraw[i][j] = -1;
            SlotFromRingStraw[i][j] = -1;
            ChannelFromRingStraw[i][j] = -1;
        }
    } 
    // Some information

    int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
    double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};

    // Define a different 2D histogram for each ring. X-axis is phi, Y-axis is radius (to plot correctly with "pol" option)

    // create root folder for cdc and cd to it, store main dir
    TDirectory *main = gDirectory;
    gDirectory->mkdir("CDC_Efficiency")->cd();
    gDirectory->mkdir("CDC_View")->cd();

    cdc_measured_ring.resize(29);
    cdc_expected_ring.resize(29);
    for(int locDOCABin = 0; locDOCABin < 8; ++locDOCABin)
    {
    	cdc_measured_ringmap[locDOCABin].resize(29);
    	cdc_expected_ringmap[locDOCABin].resize(29);
    }

    for(int iring=0; iring<28; iring++){
        double r_start = radius[iring] - 0.8;
        double r_end = radius[iring] + 0.8;
        double phi_start = phi[iring]; // this is for center of straw. Need additional calculation for phi at end plate
        double phi_end = phi_start + TMath::TwoPi();

        char hname_measured[256];
        char hname_expected[256];
        sprintf(hname_measured, "cdc_measured_ring[%d]", iring+1);
        sprintf(hname_expected, "cdc_expected_ring[%d]", iring+1);

		cdc_measured_ring[iring+1] = new TH2D(hname_measured, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);

		cdc_measured_ringmap[0][iring+1] = new TH2D((TString)hname_measured +"DOCA0", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[1][iring+1] = new TH2D((TString)hname_measured +"DOCA1", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[2][iring+1] = new TH2D((TString)hname_measured +"DOCA2", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[3][iring+1] = new TH2D((TString)hname_measured +"DOCA3", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[4][iring+1] = new TH2D((TString)hname_measured +"DOCA4", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[5][iring+1] = new TH2D((TString)hname_measured +"DOCA5", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[6][iring+1] = new TH2D((TString)hname_measured +"DOCA6", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_measured_ringmap[7][iring+1] = new TH2D((TString)hname_measured +"DOCA7", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);

		cdc_expected_ring[iring+1] = new TH2D(hname_expected, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[0][iring+1] = new TH2D((TString)hname_expected + "DOCA0", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[1][iring+1] = new TH2D((TString)hname_expected + "DOCA1", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[2][iring+1] = new TH2D((TString)hname_expected + "DOCA2", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[3][iring+1] = new TH2D((TString)hname_expected + "DOCA3", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[4][iring+1] = new TH2D((TString)hname_expected + "DOCA4", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[5][iring+1] = new TH2D((TString)hname_expected + "DOCA5", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[6][iring+1] = new TH2D((TString)hname_expected + "DOCA6", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
		cdc_expected_ringmap[7][iring+1] = new TH2D((TString)hname_expected + "DOCA7", "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
    }	

    gDirectory->cd("/CDC_Efficiency");
    gDirectory->mkdir("Track_Quality")->cd();
    hChi2OverNDF = new TH1I("hChi2OverNDF","hChi2OverNDF", 500, 0.0, 1.0);
    ChargeVsTrackLength = new TH2I("ChargeVsTrackLength", "ChargeVsTrackLength", 1000, 0, 8.0, 2000, 0, 5000000);
    hResVsT = new TH2I("hResVsT","Tracking Residual (Biased) Vs Drift Time; Drift Time [ns]; Residual [cm]", 500, 0.0, 700.0, 1000, -0.5, 0.5);
    main->cd();

	dTargetCenterZ = 65.0;
	dTargetLength = 30.0;

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

	//Get Target Center Z, length
	dgeom->GetTargetZ(dTargetCenterZ);
	dgeom->GetTargetLength(dTargetLength);

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

    //Pre-sort hits by ring to save time //only need to search within the given ring, straw
    map<int, map<int, set<const DCDCHit*> > > locSortedCDCHits; //first int: ring //second int: straw
    for(auto& locHit : locCDCHitVector)
    	locSortedCDCHits[locHit->ring][locHit->straw].insert(locHit);

    const DDetectorMatches *detMatches = nullptr;
    if(!dIsNoFieldFlag)
        loop->GetSingle(detMatches);

    const DParticleID *locParticleID = nullptr;
    loop->GetSingle(locParticleID);

    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Cut very loosely on the track quality
        const DTrackTimeBased *thisTimeBasedTrack = nullptr;
        bestHypothesis->GetSingle(thisTimeBasedTrack);

        japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
        hChi2OverNDF->Fill(thisTimeBasedTrack->FOM);
        japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

        if (thisTimeBasedTrack->FOM < dMinTrackingFOM)
        	continue;

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
            if(!detMatches->Get_IsMatchedToDetector(thisTimeBasedTrack, SYS_TOF) && !detMatches->Get_IsMatchedToDetector(thisTimeBasedTrack, SYS_BCAL))
            {
                //cout << "Cut on detector matches" << endl;
                continue; // Require there to be at least one match to BCAL or TOF //not SC: lights up like xmas tree
            }
            if(fabs(thisTimeBasedTrack->position().Z() - dTargetCenterZ) > dTargetLength/2.0 || thisTimeBasedTrack->position().Perp() > 1.0) {
                //cout << " Cut on vertex " << endl;
                continue; // Cut on reconstructed vertex location
            }
        }

        // Require hits on at least 2 axial layers and at least 2 stereo layers:
        // necessary to trust reconstructed phi & theta: respectable projection
        set<int> locCDCRings;
    	locParticleID->Get_CDCRings(thisTimeBasedTrack->dCDCRings, locCDCRings);

        map<int, int> locNumHitRingsPerSuperlayer; //key: superlayer (1 -> 7) //axial: 1, 4, 7
        locParticleID->Get_CDCNumHitRingsPerSuperlayer(locCDCRings, locNumHitRingsPerSuperlayer);

        int locNumSuperLayersWith2Hits_Axial = 0;
        int locNumSuperLayersWith2Hits_Stereo = 0;
        for(auto& locSuperlayerPair : locNumHitRingsPerSuperlayer)
        {
        	if(locSuperlayerPair.second < 2)
        		continue;
        	if((locSuperlayerPair.first == 1) || (locSuperlayerPair.first == 4) || (locSuperlayerPair.first == 7))
        		++locNumSuperLayersWith2Hits_Axial;
        	else
        		++locNumSuperLayersWith2Hits_Stereo;
        }

        if((locNumSuperLayersWith2Hits_Axial < 2) || (locNumSuperLayersWith2Hits_Stereo < 2))
        	continue; //don't trust the track projections

        // Alright now we truly have the tracks we are interested in for calculating the efficiency
        //BUT, we need to make sure that we aren't biased by the fact that the track was reconstructed in the first place
        	//AND by our requirement above that there be at least 2 hits in a few superlayers
        for(int locCDCSuperlayer = 1; locCDCSuperlayer <= 7; ++locCDCSuperlayer)
        {
			if(locNumHitRingsPerSuperlayer[locCDCSuperlayer] < dMinNumRingsToEvalSuperlayer)
				continue;

        	int locFirstRing = 4*(locCDCSuperlayer - 1) + 1;
        	if(locNumHitRingsPerSuperlayer[locCDCSuperlayer] == dMinNumRingsToEvalSuperlayer)
			{
				//All hits required: Can only evaluate the rings that do NOT have hits
	        	for (int locRing = locFirstRing; locRing < locFirstRing + 4; ++locRing)
	        	{
	        		if(locCDCRings.find(locRing) == locCDCRings.end())
	        			GitRDone(locRing, thisTimeBasedTrack, locSortedCDCHits); // git-r-dun
	        	}
	        	continue;
			}
        	//so many hits that no individual ring was required: evaluate for all
        	for (int locRing = locFirstRing; locRing < locFirstRing + 4; ++locRing)
       			GitRDone(locRing, thisTimeBasedTrack, locSortedCDCHits); // git-r-dun
    	}
    }
    return NOERROR;
}

void JEventProcessor_CDC_Efficiency::GitRDone(unsigned int ringNum, const DTrackTimeBased *thisTimeBasedTrack, map<int, map<int, set<const DCDCHit*> > >& locSortedCDCHits)
{
	vector< DCDCWire * > wireByNumber = cdcwires[ringNum - 1];
	for (unsigned int wireIndex = 0; wireIndex < wireByNumber.size(); wireIndex++)
	{
		int wireNum = wireIndex+1;
		DCDCWire * wire = wireByNumber[wireIndex];
		double wireLength = wire->L;
		double distanceToWire = thisTimeBasedTrack->rt->DistToRT(wire, &wireLength);

		//SKIP IF NOT CLOSE
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

		double delta = 0.0, dz = 0.0;
		if(!Expect_Hit(thisTimeBasedTrack, wire, distanceToWire, delta, dz))
			continue;

		//FILL EXPECTED HISTOGRAMS
		double dx = thisTimeBasedTrack->rt->Straw_dx(wire, 0.78);
		double locTheta = thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg();
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs Path Length", dx, "Expected Hits", 100, 0 , 4.0);
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs DOCA", distanceToWire, "Expected Hits", 100, 0 , 0.78);
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs Tracking FOM", thisTimeBasedTrack->FOM, "Expected Hits", 100, 0 , 1.0);
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs theta", locTheta, "Expected Hits", 100, 0, 180);
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs p", thisTimeBasedTrack->pmag(), "Expected Hits", 100, 0 , 4.0);
		Fill1DHistogram("CDC_Efficiency", "Offline", "Expected Hits Vs delta", delta, "Expected Hits", 100, -0.3 , 0.3);
		Fill2DHistogram("CDC_Efficiency", "Offline", "Expected hits p Vs Theta", locTheta, thisTimeBasedTrack->pmag(), "Expected Hits", 100, 0, 180, 100, 0 , 4.0);

		// look for a CDC hit match
		// We need a backwards map from ring/straw to flash channel. Unfortunately there is no easy way
		// Will construct the map manually
		const DCDCHit* locHit = Find_Hit(ringNum, wireNum, locSortedCDCHits[ringNum]);
		bool foundHit = (locHit != nullptr);
		if(foundHit)
		{
			const DCDCDigiHit *thisDigiHit = NULL;
			const Df125CDCPulse *thisPulse = NULL;
			locHit->GetSingle(thisDigiHit);
			if (thisDigiHit != NULL)
				thisDigiHit->GetSingle(thisPulse);
			if (thisPulse != NULL)
			{
				ROCIDFromRingStraw[ringNum - 1][wireNum - 1] = thisPulse->rocid;
				SlotFromRingStraw[ringNum - 1][wireNum - 1] = thisPulse->slot;
				ChannelFromRingStraw[ringNum - 1][wireNum - 1] = thisPulse->channel;
			}
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs Path Length", dx, "Measured Hits", 100, 0 , 4.0);
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs DOCA", distanceToWire, "Measured Hits", 100, 0 , 0.78);
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs Tracking FOM", thisTimeBasedTrack->FOM, "Measured Hits", 100, 0 , 1.0);
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs theta", locTheta, "Measured Hits", 100, 0, 180);
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs p", thisTimeBasedTrack->pmag(), "Measured Hits", 100, 0 , 4.0);
			Fill1DHistogram("CDC_Efficiency", "Offline", "Measured Hits Vs delta", delta, "Measured Hits", 100, -0.3 , 0.3);
			Fill2DHistogram("CDC_Efficiency", "Offline", "Measured hits p Vs Theta", locTheta, thisTimeBasedTrack->pmag(), "Measured Hits", 100, 0, 180, 100, 0 , 4.0);
		}

		//FILL PROFILES: BASED ON FOUND OR NOT
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Path Length", dx,foundHit, "Efficiency; dx [cm]; Efficiency", 100, 0 , 4.0);
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs DOCA", distanceToWire,foundHit, "Efficiency; DOCA [cm]; Efficiency", 100, 0 , 0.78);
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Tracking FOM", thisTimeBasedTrack->FOM,foundHit, "Efficiency; Tracking FOM; Efficiency", 100, 0 , 1.0);
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs theta", locTheta,foundHit, "Efficiency; Track #Theta [deg]; Efficiency", 100, 0, 180);
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs p", thisTimeBasedTrack->pmag(),foundHit, "Efficiency; Momentum [GeV]; Efficiency", 100, 0 , 4.0);
		Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs delta", delta,foundHit, "Efficiency; #delta [cm]; Efficiency", 100, -0.3 , 0.3);
		if( ChannelFromRingStraw[ringNum - 1][wireNum - 1] != -1)
		{
			Fill1DProfile("CDC_Efficiency", "Online", "Efficiency Vs Channel Number", ChannelFromRingStraw[ringNum - 1][wireNum - 1],foundHit, "Efficiency; Channel Number; Efficiency", 73, -0.5 , 72.5);
			char name [200];
			sprintf(name, "Slot Efficiency ROCID %.2i", ROCIDFromRingStraw[ringNum - 1][wireNum - 1]);
			Fill1DProfile("CDC_Efficiency", "Online", name, SlotFromRingStraw[ringNum - 1][wireNum - 1],foundHit, "Efficiency; Slot Number; Efficiency", 21, -0.5 , 20.5);
			sprintf(name, "Channel Efficiency ROCID %.2i", ROCIDFromRingStraw[ringNum - 1][wireNum - 1]);
			Fill1DProfile("CDC_Efficiency", "Online", name, SlotFromRingStraw[ringNum - 1][wireNum - 1] * 100 + ChannelFromRingStraw[ringNum - 1][wireNum - 1],foundHit, "Efficiency; Channel; Efficiency", 1501, 299.5 , 1800.5);
		}

		Fill2DProfile("CDC_Efficiency", "Online", "Efficiency p Vs Theta", locTheta, thisTimeBasedTrack->pmag(),foundHit, "Efficiency; Track #Theta [deg]; Momentum [GeV]", 100, 0, 180, 100, 0 , 4.0);
		Fill2DProfile("CDC_Efficiency", "Online", "Efficiency distance Vs delta", delta,distanceToWire,foundHit, "Efficiency;#delta [cm]; DOCA [cm]", 100, -0.3, 0.3, 100, 0 , 1.2);
		Fill2DProfile("CDC_Efficiency", "Online", "Efficiency z Vs delta", delta,dz,foundHit, "Efficiency;#delta [cm]; z [cm] (CDC local coordinates)", 100, -0.3, 0.3, 150, -75 , 75);

		//FILL AS FUNCTION OF DOCA
		if (distanceToWire < 0.78)
		{
			Fill_ExpectedHit(ringNum, wireNum, distanceToWire);
			if(foundHit)
				Fill_MeasuredHit(ringNum, wireNum, distanceToWire, thisTimeBasedTrack, wire, locHit);
		}
	}
}

bool JEventProcessor_CDC_Efficiency::Expect_Hit(const DTrackTimeBased* thisTimeBasedTrack, DCDCWire* wire, double distanceToWire, double& delta, double& dz)
{
	delta = 0.0;
	dz = 0.0;
	if (distanceToWire >= 1.2 )
		return false;

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

	return (distanceToWire < (0.78 + delta) && fabs(dz) < 65.0);
}

void JEventProcessor_CDC_Efficiency::Fill_MeasuredHit(int ringNum, int wireNum, double distanceToWire, const DTrackTimeBased* thisTimeBasedTrack, DCDCWire* wire, const DCDCHit* locHit)
{
	//Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
	//cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
	//Fill the expected number of hits histogram
	if (distanceToWire < DOCACUT)
	{
		//printf("Matching Hit!!!!!\n");
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		{
			Double_t v = cdc_measured_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
			cdc_measured_ring[ringNum]->SetBinContent(wireNum, 1, v);
			double dx = thisTimeBasedTrack->rt->Straw_dx(wire, 0.78);
			ChargeVsTrackLength->Fill(dx, locHit->q);

			Double_t w = cdc_expected_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
			cdc_measured_ring[ringNum]->SetBinContent(wireNum, 1, w);
		}
		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	}

	int locDOCABin = (int) (distanceToWire * 10) % 8;
	TH2D* locHistToFill = cdc_measured_ringmap[locDOCABin][ringNum];
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	{
		Double_t w = locHistToFill->GetBinContent(wireNum, 1) + 1.0;
		locHistToFill->SetBinContent(wireNum, 1, w);
	}
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
}

void JEventProcessor_CDC_Efficiency::Fill_ExpectedHit(int ringNum, int wireNum, double distanceToWire)
{
	//Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
	//cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
	//Fill the expected number of hits histogram
	if (distanceToWire < DOCACUT)
	{
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		{
			Double_t w = cdc_expected_ring[ringNum]->GetBinContent(wireNum, 1) + 1.0;
			cdc_expected_ring[ringNum]->SetBinContent(wireNum, 1, w);
		}
		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	}

	int locDOCABin = (int) (distanceToWire * 10) % 8;
	TH2D* locHistToFill = cdc_expected_ringmap[locDOCABin][ringNum];
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	{
		Double_t w = locHistToFill->GetBinContent(wireNum, 1) + 1.0;
		locHistToFill->SetBinContent(wireNum, 1, w);
	}
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
}

const DCDCHit* JEventProcessor_CDC_Efficiency::Find_Hit(int locRing, int locProjectedStraw, map<int, set<const DCDCHit*> >& locSortedCDCHits)
{
	if(!locSortedCDCHits[locProjectedStraw].empty())
		return *(locSortedCDCHits[locProjectedStraw].begin());

	int locNumStraws = cdcwires[locRing - 1].size();

	//previous straw
	int locSearchStraw = locProjectedStraw - 1;
	if(locSearchStraw <= 0)
		locSearchStraw += locNumStraws;
	if(!locSortedCDCHits[locSearchStraw].empty())
		return *(locSortedCDCHits[locProjectedStraw].begin());

	//next straw
	locSearchStraw = locProjectedStraw + 1;
	if(locSearchStraw > locNumStraws)
		locSearchStraw -= locNumStraws;
	if(!locSortedCDCHits[locSearchStraw].empty())
		return *(locSortedCDCHits[locProjectedStraw].begin());

	return nullptr;
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

