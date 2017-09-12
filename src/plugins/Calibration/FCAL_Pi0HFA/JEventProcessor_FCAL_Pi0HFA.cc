// $Id$
//
//    File: JEventProcessor_FCAL_Pi0HFA.cc
// Created: Wed Aug 30 16:23:49 EDT 2017
// Creator: mstaib (on Linux egbert 2.6.32-696.10.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_FCAL_Pi0HFA.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FCAL_Pi0HFA());
}
} // "C"


//------------------
// JEventProcessor_FCAL_Pi0HFA (Constructor)
//------------------
JEventProcessor_FCAL_Pi0HFA::JEventProcessor_FCAL_Pi0HFA()
{

}

//------------------
// ~JEventProcessor_FCAL_Pi0HFA (Destructor)
//------------------
JEventProcessor_FCAL_Pi0HFA::~JEventProcessor_FCAL_Pi0HFA()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCAL_Pi0HFA::init(void)
{
	// This is called once at program startup. 

   gDirectory->mkdir("FCAL_Pi0HFA");
   gDirectory->cd("FCAL_Pi0HFA");
   hCurrentGainConstants = new TProfile("CurrentGainConstants", "Current Gain Constants", 2800, -0.5, 2799.5);
   gDirectory->cd("..");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCAL_Pi0HFA::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes

   // Put the current gain constants into the output file
   vector< double > raw_gains;
   // This is called whenever the run number changes
   eventLoop->GetCalib("/FCAL/gains", raw_gains);
   for (unsigned int i=0; i<raw_gains.size(); i++){
      hCurrentGainConstants->Fill(i,raw_gains[i]);
   }


   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCAL_Pi0HFA::evnt(JEventLoop *loop, uint64_t eventnumber)
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
   // japp->RootFillLock(this);
   //  ... fill historgrams or trees ...
   // japp->RootFillUnLock(this);

   vector<const DNeutralParticle *> neutralParticleVector;
   loop->Get(neutralParticleVector);

   // Cut at most 6 neutral particles
   if (neutralParticleVector.size() > 6 || neutralParticleVector.size() < 2) return NOERROR;

   for (unsigned int i=0; i< neutralParticleVector.size() - 1; i++){
      const DNeutralParticleHypothesis *photon1 = neutralParticleVector[i]->Get_Hypothesis(Gamma);
      // Go into the FCAL shower and find the largest energy deposition
      const DNeutralShower *shower1;
      photon1->GetSingle(shower1);
      if(shower1->dDetectorSystem != SYS_FCAL) continue;
      DFCALShower *fcalShower1 = (DFCALShower *) shower1->dBCALFCALShower;
      const DFCALCluster *fcalCluster1;
      fcalShower1->GetSingle(fcalCluster1);
      int ch1 = fcalCluster1->getChannelEmax();
      double frac1 = fcalCluster1->getEmax()/fcalCluster1->getEnergy();
      if(fcalCluster1->getEnergy() < 0.8) continue;
      for (unsigned int j=i+1; j< neutralParticleVector.size(); j++){
         const DNeutralParticleHypothesis *photon2 = neutralParticleVector[j]->Get_Hypothesis(Gamma);
         const DNeutralShower *shower2;
         photon2->GetSingle(shower2);
         if(shower2->dDetectorSystem != SYS_FCAL) continue;
         DFCALShower *fcalShower2 = (DFCALShower *) shower2->dBCALFCALShower;
         const DFCALCluster *fcalCluster2;
         fcalShower2->GetSingle(fcalCluster2);
         int ch2 = fcalCluster2->getChannelEmax();
         double frac2 = fcalCluster2->getEmax()/fcalCluster2->getEnergy();
         if(fcalCluster2->getEnergy() < 0.8) continue;

         double pi0Mass = (photon1->lorentzMomentum() + photon2->lorentzMomentum()).M();
         double avgE = 0.5*fcalCluster1->getEnergy() + 0.5*fcalCluster2->getEnergy();


         Fill1DHistogram("FCAL_Pi0HFA","","Pi0Mass",
               pi0Mass,
               "#pi^{0} Mass; #pi^{0} Mass;",
               100, 0.05, 0.25);

         if(frac1 > 0.5){
            Fill2DHistogram("FCAL_Pi0HFA","","Pi0MassVsChNum",
                  ch1, pi0Mass,
                  "#pi^{0} Mass Vs. Channel Number; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);
         }
         if(frac2 > 0.5){
            Fill2DHistogram("FCAL_Pi0HFA","","Pi0MassVsChNum",
                  ch2, pi0Mass,
                  "#pi^{0} Mass Vs. Channel Number; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);
         }
         // Energy Weighted Pi0 Mass
         for(auto hit : fcalCluster1->GetHits()){
            Fill2DWeightedHistogram("FCAL_Pi0HFA","","Pi0MassVsChNumWeighted",
                  hit.ch, pi0Mass, hit.E / fcalCluster1->getEnergy(),
                  "#pi^{0} Mass Vs. Channel Number Weighted; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);
            Fill2DWeightedHistogram("FCAL_Pi0HFA","","Pi0MassVsChNumWeightedSquared",
                  hit.ch, pi0Mass, (hit.E / fcalCluster1->getEnergy())*(hit.E / fcalCluster1->getEnergy()),
                  "#pi^{0} Mass Vs. Channel Number Weighted; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);

         }

         for(auto hit : fcalCluster2->GetHits()){
            Fill2DWeightedHistogram("FCAL_Pi0HFA","","Pi0MassVsChNumWeighted",
                  hit.ch, pi0Mass, hit.E / fcalCluster2->getEnergy(),
                  "#pi^{0} Mass Vs. Channel Number Weighted; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);
            Fill2DWeightedHistogram("FCAL_Pi0HFA","","Pi0MassVsChNumWeightedSquared",
                  hit.ch, pi0Mass, (hit.E / fcalCluster2->getEnergy())*(hit.E / fcalCluster2->getEnergy()),
                  "#pi^{0} Mass Vs. Channel Number Weighted; CCDB Index; #pi^{0} Mass",
                  2800, -0.5, 2799.5, 200, 0.00, 0.5);

         }

         if (fabs(fcalCluster1->getEnergy() - fcalCluster1->getEnergy()) < 0.25){
            Fill2DHistogram("FCAL_Pi0HFA","","Pi0MassVsE",
                  avgE, pi0Mass,
                  "#pi^{0} Mass Vs. Average Shower Energy; CCDB Index; #pi^{0} Mass",
                  100, 0.0, 10.0, 100, 0.05, 0.25);
         }
      }
   }

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCAL_Pi0HFA::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCAL_Pi0HFA::fini(void)
{
   // Called before program exit after event processing is finished.
   return NOERROR;
}

