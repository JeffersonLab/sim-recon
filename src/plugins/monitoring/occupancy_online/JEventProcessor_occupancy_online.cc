// $Id$
//
//    File: JEventProcessor_occupancy_online.cc
// Created: Tue Apr 12 09:43:54 EDT 2016
// Creator: zihlmann (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_occupancy_online.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_occupancy_online());
  }
} // "C"


//------------------
// JEventProcessor_occupancy_online (Constructor)
//------------------
JEventProcessor_occupancy_online::JEventProcessor_occupancy_online()
{
  
}

//------------------
// ~JEventProcessor_occupancy_online (Destructor)
//------------------
JEventProcessor_occupancy_online::~JEventProcessor_occupancy_online()
{
  
}

//------------------
// init
//------------------
jerror_t JEventProcessor_occupancy_online::init(void)
{
  // This is called once at program startup. If you are creating
  // and filling historgrams in this plugin, you should lock the
  // ROOT mutex like this:
  //
  // japp->RootWriteLock();
  //  ... fill historgrams or trees ...
  // japp->RootUnLock();
  //

  // Hit occupancy
  int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 
		     135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
  double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 
		       23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 
		       34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 
		       45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
  double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 
		    0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 
		    0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 
		    0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 
		    0.005348855, 0.005997358, 0.021018761}; 
  
  // Define a different 2D histogram for each ring. 
  // X-axis is phi, Y-axis is radius (to plot correctly with "pol" option)
  for(int iring=0; iring<28; iring++){
    double r_start = radius[iring] - 0.8;
    double r_end = radius[iring] + 0.8;
    double phi_start = phi[iring]; // this is for center of straw. Need additional calculation for phi at end plate
    double phi_end = phi_start + TMath::TwoPi();
    
    char hname[256];
    sprintf(hname, "cdc_occ_ring[%d]", iring);
    cdc_occ_ring[iring] = new TH2F(hname, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
  }
  
  

  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_occupancy_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes


  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_occupancy_online::evnt(JEventLoop *loop, uint64_t eventnumber)
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
 
  vector<const DCDCDigiHit*> digihits;
  loop->Get(digihits);
  
  for(uint32_t i=0; i<digihits.size(); i++) {
    
    const DCDCDigiHit *digihit = digihits[i];  
    int ring     = digihit->ring-1; // start counting from zero! 
    int straw    = digihit->straw;  // first bin is one
    int n = straw_offset[ring] + straw;

    Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
    cdc_occ_ring[ring]->SetBinContent(straw, 1, w);

  }
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_occupancy_online::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_occupancy_online::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

