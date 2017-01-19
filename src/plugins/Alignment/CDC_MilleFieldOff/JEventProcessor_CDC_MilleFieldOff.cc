// $Id$
//
//    File: JEventProcessor_CDC_MilleFieldOff.cc
// Created: Tue Jan 17 19:32:32 Local time zone must be set--see zic manual page 2017
// Creator: mstaib (on Linux egbert 2.6.32-642.6.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_CDC_MilleFieldOff.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "TRACKING/DTrackCandidate.h"
#include "CDC/DCDCTrackHit.h"
#include "CDC/DCDCWire.h"
#include "TDirectory.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
   InitJANAPlugin(app);
   app->AddProcessor(new JEventProcessor_CDC_MilleFieldOff());
}
} // "C"


//------------------
// JEventProcessor_CDC_MilleFieldOff (Constructor)
//------------------
JEventProcessor_CDC_MilleFieldOff::JEventProcessor_CDC_MilleFieldOff()
{

}

//------------------
// ~JEventProcessor_CDC_MilleFieldOff (Destructor)
//------------------
JEventProcessor_CDC_MilleFieldOff::~JEventProcessor_CDC_MilleFieldOff()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_CDC_MilleFieldOff::init(void)
{
   // This is called once at program startup.
   milleWriter = new Mille("cdc_mille_out.mil");

   gDirectory->mkdir("AlignmentConstants");
   gDirectory->cd("AlignmentConstants");
   // We need the constants used for this iteration
   // Use a TProfile to avoid problems adding together multiple root files...
   HistCurrentConstants = new TProfile("CDCAlignmentConstants", "Constants Used for CDC Alignment (In MILLEPEDE Order)", 16000 ,0.5, 16000.5);

   gDirectory->cd("..");

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CDC_MilleFieldOff::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   // Get the current set of constants and sve them in the histogram
   // This is called whenever the run number changes
   // Check for magnetic field
   DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
   bool dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);

   // This plugin is designed for field off data. If this is used for field on data, Abort...
   if (!dIsNoFieldFlag){
      jerr << " Plugin FDC_MilleFieldOff Must be run with zero magnetic field!!! Aborting " << endl;
      jerr << " Use -PBFIELD_TYPE=NoField " << endl;
      japp->Quit();
   }

   // Store the current values of the alignment constants
   JCalibration * jcalib = eventLoop->GetJCalibration();
   vector<map<string,double> >vals;
   if (jcalib->Get("CDC/global_alignment",vals)==false){
      map<string,double> &row = vals[0];
      // Get the offsets from the calibration database
      HistCurrentConstants->Fill(1,row["dX"]);
      HistCurrentConstants->Fill(2,row["dY"]);
      HistCurrentConstants->Fill(3,row["dZ"]);
      HistCurrentConstants->Fill(4,row["dPhiX"]);
      HistCurrentConstants->Fill(5,row["dPhiY"]);
      HistCurrentConstants->Fill(6,row["dPhiZ"]);
   }

   if (jcalib->Get("CDC/wire_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstants->Fill(1000+(i*4+1),row["dxu"]);
         HistCurrentConstants->Fill(1000+(i*4+2),row["dyu"]);
         HistCurrentConstants->Fill(1000+(i*4+3),row["dxd"]);
         HistCurrentConstants->Fill(1000+(i*4+4),row["dyd"]);
      }
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_CDC_MilleFieldOff::evnt(JEventLoop *loop, uint64_t eventnumber)
{
   int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};
   // Loop over the tracks, get the tracking pulls, and fill some histograms. Easy peasy
   vector<const DTrackCandidate *> trackVector;
   loop->Get(trackVector,"StraightLine");

   for (size_t i = 0; i < trackVector.size(); i++){
      const DTrackCandidate *track = trackVector[i];

      double trackingFOM = TMath::Prob(track->chisq, track->Ndof);

      // Some quality cuts for the tracks we will use
      // Keep this minimal for now and investigate later
      double trackingFOMCut = 0.001;
      int trackingNDFCut = 5;

      if(trackingFOM < trackingFOMCut) continue;
      if(track->Ndof < trackingNDFCut) continue;

      //double phi = bestHypothesis->momentum().Phi()*TMath::RadToDeg();
      //double theta = bestHypothesis->momentum().Theta()*TMath::RadToDeg();

      if (!track->IsSmoothed) continue;

      vector<DTrackFitter::pull_t> pulls = track->pulls;
      japp->RootWriteLock(); // Just use the root lock as a temporary
      for (size_t iPull = 0; iPull < pulls.size(); iPull++){
         // Here is all of the information currently stored in the pulls from the fit
         // From TRACKING/DTrackFitter.h
         double resi                 = pulls[iPull].resi;   // residual of measurement
         double err                  = pulls[iPull].err;      // estimated error of measurement
         //double s                    = pulls[iPull].s;
         //double tdrift               = pulls[iPull].tdrift;      // drift time of this measurement
         //double d                    = pulls[iPull].d;  // doca to wire
         const DCDCTrackHit *cdc_hit = pulls[iPull].cdc_hit;
         if (cdc_hit == NULL) continue;

         vector<double> trackDerivatives = pulls[iPull].trackDerivatives;
         vector<double> stateVector = pulls[iPull].stateVector;

         const DCDCWire *constWire = cdc_hit->wire;
         DCDCWire *thisWire = const_cast<DCDCWire *>(constWire);

         vector<double> wireDerivatives = thisWire->GetCDCWireDerivatives();

         //Now we need to fill the Mille structure once for our wire measurment and once for our cathode measurement
         const int NLC = 4; // Number of local parameters
         const int NGL = 10; // Number of global parameters for wire alignment
         float localDer[NLC];
         float globalDer[NGL];
         int label[NGL];

         localDer[0]=trackDerivatives[CDCTrackD::dDOCAdx]; localDer[1]=trackDerivatives[CDCTrackD::dDOCAdy];
         localDer[2]=trackDerivatives[CDCTrackD::dDOCAdtx]; localDer[3]=trackDerivatives[CDCTrackD::dDOCAdty];

         globalDer[0]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaX]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaX]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaX]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaX]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaX]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaX];
         label[0]=1;

         globalDer[1]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaY]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaY]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaY]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaY]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaY]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaY];
         label[1]=2;

         globalDer[2]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaZ]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaZ]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaZ];
         label[2]=3;

         globalDer[3]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaPhiX]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaPhiX]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaPhiX]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaPhiX]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaPhiX]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaPhiX];
         label[3]=4;

         globalDer[4]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaPhiY]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaPhiY]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaPhiY]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaPhiY]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaPhiY]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaPhiY];
         label[4]=5;

         globalDer[5]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaPhiZ]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaPhiZ]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaPhiZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaPhiZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaPhiZ]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaPhiZ];
         label[5]=6;

         globalDer[6]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaXu]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaXu]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaXu]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaXu]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaXu]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaXu];
         label[6]=1000 + (straw_offset[thisWire->ring]+(thisWire->straw-1))*4 + 1;

         globalDer[7]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaYu]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaYu]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaYu]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaYu]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaYu]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaYu];
         label[7]=1000 + (straw_offset[thisWire->ring]+(thisWire->straw-1))*4 + 2;

         globalDer[8]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaXd]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaXd]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaXd]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaXd]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaXd]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaXd];
         label[8]=1000 + (straw_offset[thisWire->ring]+(thisWire->straw-1))*4 + 3;

         globalDer[9]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaYd]
            +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaYd]
            +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaYd]
            +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaYd]
            +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaYd]
            +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaYd];
         label[9]=1000 + (straw_offset[thisWire->ring]+(thisWire->straw-1))*4 + 4;

         milleWriter->mille(NLC, localDer, NGL, globalDer, label, resi, err);

      }
      milleWriter->end();

      japp->RootUnLock();

   }// End loop over tracks

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_CDC_MilleFieldOff::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_CDC_MilleFieldOff::fini(void)
{
   // Called before program exit after event processing is finished.
   delete milleWriter;
   return NOERROR;
}

