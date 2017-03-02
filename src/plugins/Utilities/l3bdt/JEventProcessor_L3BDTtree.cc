// $Id$
//
//    File: JEventProcessor_L3BDTtree.cc
// Created: Wed May 11 22:26:46 EDT 2016
// Creator: davidl (on Linux gluon49.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//


#include "JEventProcessor_L3BDTtree.h"
using namespace jana;

#include <TMath.h>



// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_L3BDTtree());
}
} // "C"




//------------------
// JEventProcessor_L3BDTtree (Constructor)
//------------------
JEventProcessor_L3BDTtree::JEventProcessor_L3BDTtree()
{

}

//------------------
// ~JEventProcessor_L3BDTtree (Destructor)
//------------------
JEventProcessor_L3BDTtree::~JEventProcessor_L3BDTtree()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_L3BDTtree::init(void)
{

	// Create Tree
	l3tree = new TTree("l3tree", "L3 tree for BDT");

	// Add branches for all JANA object counts
	#define AddNBranch(A) l3tree->Branch("N" #A, &bdt.N##A, #A "/F");
	MyTypes(AddNBranch)

	// Add branches for all sorted float types
	#define AddBranch(A) l3tree->Branch(#A, &bdt.A, #A "/F");
	MyDerivedTypes(AddBranch)

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_L3BDTtree::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_L3BDTtree::evnt(JEventLoop *loop, uint64_t eventnumber)
{

	// Get all JANA objects of interest
	#define GetObjs(A) vector<const A*> v##A; loop->Get(v##A);
	MyTypes(GetObjs)
	
	// Trigger
	Int_t   trig_mask = 0;
	Int_t   fp_trig_mask = 0;
	Float_t trig1 = 0.0;
	Float_t trig3 = 0.0;
	Float_t trig4 = 0.0;
	if(!vDL1Trigger.empty()){
		auto trig  = vDL1Trigger[0];
		trig_mask    = trig->trig_mask;
		fp_trig_mask = trig->fp_trig_mask;
		if(trig_mask&0x01) trig1 = 1.0;
		if(trig_mask&0x04) trig3 = 1.0;
		if(trig_mask&0x08) trig4 = 1.0;
	}

	// CDC hits by superlayer
	Int_t NCDC_superlayer[5] = {0,0,0,0,0};
	vector<int> superlayer_maxring = { 4, 12, 16, 24, 28 };
	for(auto hit : vDCDCDigiHit){
		for(int superlayer=0; superlayer<5; superlayer++){
			if(hit->ring <= superlayer_maxring[superlayer]){
				NCDC_superlayer[superlayer]++;
				break;
			}
		}
	}
	
	// FDC wire hits by package
	Int_t NFDCWire_package[4] = {0,0,0,0};
	for(auto hit : vDFDCWireDigiHit){
		if(hit->package>=1 && hit->package<=4) NFDCWire_package[hit->package-1]++;
	}
	
	// FDC cathode hits by package
	Int_t NFDCCathodes_package[4] = {0,0,0,0};
	for(auto hit : vDFDCCathodeDigiHit){
		if(hit->package>=1 && hit->package<=4) NFDCCathodes_package[hit->package-1]++;
	}
	
	// TOF half-length, half-width bars
	Int_t NTOF_half_length = 0;
	Int_t NTOF_half_width  = 0;
	for(auto hit : vDTOFDigiHit){
		auto &bar = hit->bar;
		if(bar==22 || bar==23 || bar==45 || bar==46) NTOF_half_length++;
		if(bar==20 || bar==21 || bar==24 || bar==25) NTOF_half_width++;
	}

	// Various total energies
	Float_t Esc_tot        = 0.0;
	Float_t Etof_tot       = 0.0;
	Float_t Ebcal_points   = 0.0;
	Float_t Ebcal_clusters = 0.0;
	Float_t Efcal_clusters = 0.0;
	for(auto sc : vDSCHit      ) Esc_tot += sc->dE;
	for(auto th : vDTOFPoint   ) Etof_tot += th->dE;
	for(auto bp : vDBCALPoint  ) Ebcal_points   += bp->E();
	for(auto bc : vDBCALCluster) Ebcal_clusters += bc->E();
	for(auto fc : vDFCALCluster) Efcal_clusters += fc->getEnergy();

	// FCAL Rmin and Rmax (for Eugene)
	Float_t Rfcal_min = 10000.0;
	Float_t Rfcal_max = 0.0;
	for(auto fc : vDFCALCluster){
		Float_t r = fc->getCentroid().Perp();
		if( r < Rfcal_min ) Rfcal_min = r;
		if( r > Rfcal_max ) Rfcal_max = r;
	}

	// Beam photons
	Int_t Nbeam_photons_coherent = 0.0;
	vector<Int_t> Nbeam_photons(13, 0.0);
	double Eelectron = 11.668;
	double Ecoherent_min = 8.4*Eelectron/12.0;
	double Ecoherent_max = 9.0*Eelectron/12.0;
	for(auto p : vDBeamPhoton){
		double Ephoton = p->energy();
		int idx = floor(Ephoton);
		if(idx>=0 && idx<=12) Nbeam_photons[idx]++;
		if( (Ephoton>=Ecoherent_min) && (Ephoton<=Ecoherent_max) ) Nbeam_photons_coherent++;
	}

	// Ptot for candidates
	double Ptot_candidates = 0.0;
	for(auto tc : vDTrackCandidate) Ptot_candidates += tc->momentum().Mag();

	// Visible energy
	Float_t Evisible_neutrals = 0.0;
	Float_t Evisible_tracks = 0.0;
	Float_t Evisible_charged_Kaons = 0.0;
	Float_t Evisible_charged_pions = 0.0;
	Float_t Evisible_protons = 0.0;
	Float_t Evisible = 0.0;
	for(auto ct : vDChargedTrack){
		auto cth = ct->Get_BestTrackingFOM();
		double confidence_level = TMath::Prob((double)cth->dChiSq, (double)cth->dNDF);
		if(confidence_level > 0.0001){
			Float_t E = cth->energy();
			switch(cth->PID()){
				case PiPlus:
				case PiMinus:
					Evisible_charged_pions += E;
					break;
				case KPlus:
				case KMinus:
					Evisible_charged_Kaons += E;
					break;
				case Proton:
					Evisible -= cth->mass(); // Do not include proton mass in visible energy
				case AntiProton:
					Evisible_protons += E;
					break;
				default:
					break;			
			}
			
			Evisible_tracks += E;
			Evisible += E;
		}
	}
	
	// Add neutral energy
	for(auto ns : vDNeutralShower){
		Evisible_neutrals +=  ns->dEnergy;
		Evisible += ns->dEnergy;
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Get ROOT lock
	japp->RootWriteLock();
	
	// Copy all object counts
	#define CopyNobjs(A) bdt.N##A = (Float_t)v##A.size();
	MyTypes(CopyNobjs)
	
	// L1 trigger
	bdt.trig_mask    = (Float_t)trig_mask;
	bdt.fp_trig_mask = (Float_t)fp_trig_mask;
	bdt.trig1        = trig1;
	bdt.trig3        = trig3;
	bdt.trig4        = trig4;

	// CDC hits by superlayer
	bdt.NCDC_superlayer1 = (Float_t)NCDC_superlayer[0];
	bdt.NCDC_superlayer2 = (Float_t)NCDC_superlayer[1];
	bdt.NCDC_superlayer3 = (Float_t)NCDC_superlayer[2];
	bdt.NCDC_superlayer4 = (Float_t)NCDC_superlayer[3];
	bdt.NCDC_superlayer5 = (Float_t)NCDC_superlayer[4];

	// FDC wire hits by package
	bdt.NFDCwires_package1 = (Float_t)NFDCWire_package[0];
	bdt.NFDCwires_package2 = (Float_t)NFDCWire_package[1];
	bdt.NFDCwires_package3 = (Float_t)NFDCWire_package[2];
	bdt.NFDCwires_package4 = (Float_t)NFDCWire_package[3];

	// FDC cathode hits by package
	bdt.NFDCCathodes_package1 = (Float_t)NFDCCathodes_package[0];
	bdt.NFDCCathodes_package2 = (Float_t)NFDCCathodes_package[1];
	bdt.NFDCCathodes_package3 = (Float_t)NFDCCathodes_package[2];
	bdt.NFDCCathodes_package4 = (Float_t)NFDCCathodes_package[3];
	
	// TOF half-length, half-width bars
	bdt.NTOF_half_length = (Float_t)NTOF_half_length;
	bdt.NTOF_half_width  = (Float_t)NTOF_half_width;

	// Various total energies
	bdt.Esc_tot        = Esc_tot;
	bdt.Etof_tot       = Etof_tot;
	bdt.Ebcal_points   = Ebcal_points;
	bdt.Ebcal_clusters = Ebcal_clusters;
	bdt.Efcal_clusters = Efcal_clusters;

	// FCAL Rmin and Rmax (for Eugene)
	bdt.Rfcal_min = Rfcal_min;
	bdt.Rfcal_max = Rfcal_max;

	// Beam photons
	bdt.Nbeam_photons_coherent = Nbeam_photons_coherent;
	bdt.Nbeam_photons_3_4   = (Float_t)Nbeam_photons[ 3];
	bdt.Nbeam_photons_4_5   = (Float_t)Nbeam_photons[ 4];
	bdt.Nbeam_photons_5_6   = (Float_t)Nbeam_photons[ 5];
	bdt.Nbeam_photons_6_7   = (Float_t)Nbeam_photons[ 6];
	bdt.Nbeam_photons_7_8   = (Float_t)Nbeam_photons[ 7];
	bdt.Nbeam_photons_8_9   = (Float_t)Nbeam_photons[ 8];
	bdt.Nbeam_photons_9_10  = (Float_t)Nbeam_photons[ 9];
	bdt.Nbeam_photons_10_11 = (Float_t)Nbeam_photons[10];
	bdt.Nbeam_photons_11_12 = (Float_t)Nbeam_photons[11];

	// Ptot for candidates
	bdt.Ptot_candidates = Ptot_candidates;

	// Visible energy
	bdt.Evisible_neutrals      = Evisible_neutrals;
	bdt.Evisible_tracks        = Evisible_tracks;
	bdt.Evisible_charged_Kaons = Evisible_charged_Kaons;
	bdt.Evisible_charged_pions = Evisible_charged_pions;
	bdt.Evisible_protons       = Evisible_protons;
	bdt.Evisible               = Evisible;

	// Fill tree
	l3tree->Fill();

	// Release ROOT lock
	japp->RootUnLock();

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_L3BDTtree::erun(void)
{
	japp->RootWriteLock();
	l3tree->FlushBaskets();
	japp->RootUnLock();

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_L3BDTtree::fini(void)
{
	japp->RootWriteLock();
	l3tree->FlushBaskets();
	japp->RootUnLock();
		
	return NOERROR;
}

