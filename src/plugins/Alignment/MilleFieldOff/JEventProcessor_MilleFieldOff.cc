// $Id$
//
//    File: JEventProcessor_MilleFieldOff.cc
// Created: Tue Jan 17 19:32:32 Local time zone must be set--see zic manual page 2017
// Creator: mstaib (on Linux egbert 2.6.32-642.6.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_MilleFieldOff.h"
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
   app->AddProcessor(new JEventProcessor_MilleFieldOff());
}
} // "C"


//------------------
// JEventProcessor_MilleFieldOff (Constructor)
//------------------
JEventProcessor_MilleFieldOff::JEventProcessor_MilleFieldOff()
{

}

//------------------
// ~JEventProcessor_MilleFieldOff (Destructor)
//------------------
JEventProcessor_MilleFieldOff::~JEventProcessor_MilleFieldOff()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_MilleFieldOff::init(void)
{
   // This is called once at program startup.
   milleWriter = new Mille("nofield_mille_out.mil");

   gDirectory->mkdir("AlignmentConstants");
   gDirectory->cd("AlignmentConstants");
   // We need the constants used for this iteration
   // Use a TProfile to avoid problems adding together multiple root files...
   HistCurrentConstantsCDC = new TProfile("CDCAlignmentConstants", "Constants Used for CDC Alignment (In MILLEPEDE Order)", 16000 ,0.5, 16000.5);
   HistCurrentConstantsFDC = new TProfile("FDCAlignmentConstants", "Constants Used for FDC Alignment (In MILLEPEDE Order)", 26000 ,0.5, 26000.5);

   gDirectory->cd("..");

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_MilleFieldOff::brun(JEventLoop *eventLoop, int32_t runnumber)
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
      HistCurrentConstantsCDC->Fill(1,row["dX"]);
      HistCurrentConstantsCDC->Fill(2,row["dY"]);
      HistCurrentConstantsCDC->Fill(3,row["dZ"]);
      HistCurrentConstantsCDC->Fill(4,row["dPhiX"]);
      HistCurrentConstantsCDC->Fill(5,row["dPhiY"]);
      HistCurrentConstantsCDC->Fill(6,row["dPhiZ"]);
   }

   if (jcalib->Get("CDC/wire_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstantsCDC->Fill(1000+(i*4+1),row["dxu"]);
         HistCurrentConstantsCDC->Fill(1000+(i*4+2),row["dyu"]);
         HistCurrentConstantsCDC->Fill(1000+(i*4+3),row["dxd"]);
         HistCurrentConstantsCDC->Fill(1000+(i*4+4),row["dyd"]);
      }
   }

   if (jcalib->Get("FDC/cathode_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstantsFDC->Fill((i+1)*1000+101,row["dU"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000+102,row["dV"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000+103,row["dPhiU"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000+104,row["dPhiV"]);
      }
   }

   if (jcalib->Get("FDC/cell_offsets",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstantsFDC->Fill((i+1)*1000+1,row["xshift"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000+100,row["yshift"]);
      }
   }

   if (jcalib->Get("FDC/wire_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstantsFDC->Fill((i+1)*1000+2,row["dPhi"]);
      }
   }

   if (jcalib->Get("FDC/package1/strip_gains_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         char name[32];
         for (unsigned int j=1; j<= 216; j++)
         {
            sprintf(name,"strip%i", j);
            if (i%2){ // V
               HistCurrentConstantsFDC->Fill((i/2+1)*1000+600+j,row[name]);
            }
            else {// U
               HistCurrentConstantsFDC->Fill((i/2+1)*1000+300+j,row[name]);
            }
         }
      }
   }

   if (jcalib->Get("FDC/package2/strip_gains_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         char name[32];
         for (unsigned int j=1; j<= 216; j++)
         {
            sprintf(name,"strip%i", j);
            if (i%2){ // V
               HistCurrentConstantsFDC->Fill((i/2+7)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstantsFDC->Fill((i/2+7)*1000+300+j,row[name]);
         }
      }
   }

   if (jcalib->Get("FDC/package3/strip_gains_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         char name[32];
         for (unsigned int j=1; j<= 216; j++)
         {
            sprintf(name,"strip%i", j);
            if (i%2){ // V
               HistCurrentConstantsFDC->Fill((i/2+13)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstantsFDC->Fill((i/2+13)*1000+300+j,row[name]);
         }
      }
   }

   if (jcalib->Get("FDC/package4/strip_gains_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         char name[32];
         for (unsigned int j=1; j<= 216; j++)
         {
            sprintf(name,"strip%i", j);
            if (i%2){ // V
               HistCurrentConstantsFDC->Fill((i/2+19)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstantsFDC->Fill((i/2+19)*1000+300+j,row[name]);
         }
      }
   }

   if (jcalib->Get("FDC/strip_pitches_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 200, row["U_SP_1"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 201, row["U_G_1"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 202, row["U_SP_2"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 203, row["U_G_2"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 204, row["U_SP_3"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 205, row["V_SP_1"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 206, row["V_G_1"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 207, row["V_SP_2"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 208, row["V_G_2"]);
         HistCurrentConstantsFDC->Fill((i+1)*1000 + 209, row["V_SP_3"]);
      }
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_MilleFieldOff::evnt(JEventLoop *loop, uint64_t eventnumber)
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
      int trackingNDFCut = 8;

      if(trackingFOM < trackingFOMCut) continue;
      if(track->Ndof < trackingNDFCut) continue;

      //double phi = bestHypothesis->momentum().Phi()*TMath::RadToDeg();
      //double theta = bestHypothesis->momentum().Theta()*TMath::RadToDeg();

      if (!track->IsSmoothed) continue;

      vector<DTrackFitter::pull_t> pulls = track->pulls;
      // Determine TrackType
      bool isCDCOnly=true; //bool isFDCOnly=true;
      for (size_t iPull = 0; iPull < pulls.size(); iPull++){
         const DCDCTrackHit *cdc_hit = pulls[iPull].cdc_hit;
         //const DFDCPseudo *fdc_hit   = pulls[iPull].fdc_hit;
         if (cdc_hit == NULL) isCDCOnly=false;
         //if (fdc_hit == NULL) isFDCOnly=false;
      }

      if (isCDCOnly && track->Ndof < 16) continue;

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
         const DFDCPseudo *fdc_hit   = pulls[iPull].fdc_hit;
         //double docaphi              = pulls[iPull].docaphi; // phi of doca in CDC straws
         //double z                    = pulls[iPull].z;// z position at doca
         //double tcorr                = pulls[iPull].tcorr; // drift time with correction for B
         double resic                = pulls[iPull].resic; // residual for FDC cathode measuremtns
         double errc                 = pulls[iPull].errc;

         vector<double> trackDerivatives = pulls[iPull].trackDerivatives;

         if (fdc_hit != NULL && fdc_hit->status == 6) {
            // Add fdc hit
            DFDCPseudo *thisHit = const_cast<DFDCPseudo *>(fdc_hit);

            vector<double> pseudoAlignmentDerivatives = thisHit->GetFDCPseudoAlignmentDerivatives();
            vector<double> fdcStripGainDerivatives    = thisHit->GetFDCStripGainDerivatives();
            vector<double> fdcStripPitchDerivatives   = thisHit->GetFDCStripPitchDerivatives();

            unsigned int layerOffset = 100000 + thisHit->wire->layer * 1000;

            //Now we need to fill the Mille structure once for our wire measurment and once for our cathode measurement
            const int NLC = 4; // Number of local parameters
            const int NGL_W = 2; // Number of global parameters for wire alignment
            float localDerW[NLC];
            float globalDerW[NGL_W];
            int labelW[NGL_W];

            localDerW[0]=trackDerivatives[FDCTrackD::dDOCAW_dx]; localDerW[1]=trackDerivatives[FDCTrackD::dDOCAW_dy];
            localDerW[2]=trackDerivatives[FDCTrackD::dDOCAW_dtx]; localDerW[3]=trackDerivatives[FDCTrackD::dDOCAW_dty];

            globalDerW[0] = trackDerivatives[FDCTrackD::dDOCAW_dDeltaX]; labelW[0] = layerOffset + 1;
            globalDerW[1] = trackDerivatives[FDCTrackD::dDOCAW_dDeltaPhiX]; labelW[1] = layerOffset + 2;

            milleWriter->mille(NLC, localDerW, NGL_W, globalDerW, labelW, resi, err);

            // Now for the cathode measurement, there are more alignment parameters.
            const int NGL_C = 23; // Number of global parameters for wire alignment
            float localDerC[NLC];
            float globalDerC[NGL_C];
            int labelC[NGL_C];

            localDerC[0]=trackDerivatives[FDCTrackD::dDOCAC_dx]; localDerC[1]=trackDerivatives[FDCTrackD::dDOCAC_dy];
            localDerC[2]=trackDerivatives[FDCTrackD::dDOCAC_dtx]; localDerC[3]=trackDerivatives[FDCTrackD::dDOCAC_dty];

            globalDerC[0] = -1.0;                                                labelC[0] = layerOffset + 100;
            globalDerC[1] = trackDerivatives[FDCTrackD::dDOCAC_dDeltaX];         labelC[1] = layerOffset + 1;
            globalDerC[2] = trackDerivatives[FDCTrackD::dDOCAC_dDeltaPhiX];      labelC[2] = layerOffset + 2;

            // Cathode U and V offsets
            globalDerC[3] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU];    labelC[3] = layerOffset + 101;
            globalDerC[4] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV];    labelC[4] = layerOffset + 102;
            globalDerC[5] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaPhiU]; labelC[5] = layerOffset + 103;
            globalDerC[6] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaPhiV]; labelC[6] = layerOffset + 104;

            // Strip Pitch Calibration
            globalDerC[7]  = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[0]; labelC[7] = layerOffset + 200;
            globalDerC[8]  = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[1]; labelC[8] = layerOffset + 201;
            globalDerC[9]  = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[2]; labelC[9] = layerOffset + 202;
            globalDerC[10] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[3]; labelC[10] = layerOffset + 203;
            globalDerC[11] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[4]; labelC[11] = layerOffset + 204;
            globalDerC[12] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[5]; labelC[12] = layerOffset + 205;
            globalDerC[13] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[6]; labelC[13] = layerOffset + 206;
            globalDerC[14] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[7]; labelC[14] = layerOffset + 207;
            globalDerC[15] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[8]; labelC[15] = layerOffset + 208;
            globalDerC[16] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[9]; labelC[16] = layerOffset + 209;

            // Strip Gain Calibration
            vector<const DFDCCathodeCluster *> cathodeClusters;
            thisHit->Get(cathodeClusters);
            unsigned int gainLabels[6]={};
            for(unsigned int j = 0; j< cathodeClusters.size(); j++){
               if(cathodeClusters[j]->plane == 3) { // U strips
                  unsigned int k=0;
                  for (; k < cathodeClusters[j]->members.size(); k++){
                     if (thisHit->cluster_u.index(0) == cathodeClusters[j]->members[k]->element) break;
                  }
                  gainLabels[0] = cathodeClusters[j]->members[k]->type != 3   ? thisHit->cluster_u.index(0) : thisHit->cluster_u.index(0)+108;
                  gainLabels[1] = cathodeClusters[j]->members[k+1]->type != 3 ? thisHit->cluster_u.index(1) : thisHit->cluster_u.index(1)+108;
                  gainLabels[2] = cathodeClusters[j]->members[k+2]->type != 3 ? thisHit->cluster_u.index(2) : thisHit->cluster_u.index(2)+108;
               }
               else if (cathodeClusters[j]->plane == 1) { // V strips
                  unsigned int k=0;
                  for (; k < cathodeClusters[j]->members.size(); k++){
                     if (thisHit->cluster_v.index(0) == cathodeClusters[j]->members[k]->element) break;
                  }
                  gainLabels[3] = cathodeClusters[j]->members[k]->type != 3   ? thisHit->cluster_v.index(0) : thisHit->cluster_v.index(0)+108;
                  gainLabels[4] = cathodeClusters[j]->members[k+1]->type != 3 ? thisHit->cluster_v.index(1) : thisHit->cluster_v.index(1)+108;
                  gainLabels[5] = cathodeClusters[j]->members[k+2]->type != 3 ? thisHit->cluster_v.index(2) : thisHit->cluster_v.index(2)+108;
               }
            }

            globalDerC[17] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[0]; labelC[17] = layerOffset + 300 + gainLabels[0];
            globalDerC[18] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[1]; labelC[18] = layerOffset + 300 + gainLabels[1];
            globalDerC[19] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[2]; labelC[19] = layerOffset + 300 + gainLabels[2];
            globalDerC[20] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[3]; labelC[20] = layerOffset + 600 + gainLabels[3];
            globalDerC[21] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[4]; labelC[21] = layerOffset + 600 + gainLabels[4];
            globalDerC[22] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[5]; labelC[22] = layerOffset + 600 + gainLabels[5];

            milleWriter->mille(NLC, localDerC, NGL_C, globalDerC, labelC, resic, errc);
         }

         if (cdc_hit != NULL){

            const DCDCWire *constWire = cdc_hit->wire;
            DCDCWire *thisWire = const_cast<DCDCWire *>(constWire);

            vector<double> wireDerivatives = thisWire->GetCDCWireDerivatives();

            //Now we need to fill the Mille structure once for our wire measurment and once for our cathode measurement
            const int NLC = 4; // Number of local parameters
            const int NGL = 10; // Number of global parameters for wire alignment
            float localDer[NLC];
            float globalDer[NGL];
            int label[NGL];

            localDer[0]=trackDerivatives[CDCTrackD::dDOCAdS0]; localDer[1]=trackDerivatives[CDCTrackD::dDOCAdS1];
            localDer[2]=trackDerivatives[CDCTrackD::dDOCAdS2]; localDer[3]=trackDerivatives[CDCTrackD::dDOCAdS3];

            if (isCDCOnly){ // Global shifts will not affect residuals
               globalDer[0]=0.0; label[0]=1;
               globalDer[1]=0.0; label[0]=2;
               globalDer[2]=0.0; label[0]=3;
               globalDer[3]=0.0; label[0]=4;
               globalDer[4]=0.0; label[0]=5;
               globalDer[5]=0.0; label[0]=6;
            }
            else{
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
            }
            globalDer[6]=trackDerivatives[CDCTrackD::dDOCAdOriginX]*wireDerivatives[CDCWireD::dOriginXddeltaXu]
               +trackDerivatives[CDCTrackD::dDOCAdOriginY]*wireDerivatives[CDCWireD::dOriginYddeltaXu]
               +trackDerivatives[CDCTrackD::dDOCAdOriginZ]*wireDerivatives[CDCWireD::dOriginZddeltaXu]
               +trackDerivatives[CDCTrackD::dDOCAdDirX]*wireDerivatives[CDCWireD::dDirXddeltaXu]
               +trackDerivatives[CDCTrackD::dDOCAdDirY]*wireDerivatives[CDCWireD::dDirYddeltaXu]
               +trackDerivatives[CDCTrackD::dDOCAdDirZ]*wireDerivatives[CDCWireD::dDirZddeltaXu];
            label[6]=1000 + (straw_offset[thisWire->ring]+(thisWire->straw-1))*4 + 1;

            if (false){
               jout << " Dumping deltaXu derivatives============ Wire " << thisWire->ring << " Straw " << thisWire->straw << endl;
               jout << " Total = " << globalDer[6] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdOriginX] " << trackDerivatives[CDCTrackD::dDOCAdOriginX] << " wireDerivatives[CDCWireD::dOriginXddeltaXu]" << wireDerivatives[CDCWireD::dOriginXddeltaXu] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdOriginY] " << trackDerivatives[CDCTrackD::dDOCAdOriginY] << " wireDerivatives[CDCWireD::dOriginYddeltaXu]" << wireDerivatives[CDCWireD::dOriginYddeltaXu] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdOriginZ] " << trackDerivatives[CDCTrackD::dDOCAdOriginZ] << " wireDerivatives[CDCWireD::dOriginZddeltaXu]" << wireDerivatives[CDCWireD::dOriginZddeltaXu] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdDirX] " << trackDerivatives[CDCTrackD::dDOCAdDirX] << " wireDerivatives[CDCWireD::dDirXddeltaXu]" << wireDerivatives[CDCWireD::dDirXddeltaXu] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdDirY] " << trackDerivatives[CDCTrackD::dDOCAdDirY] << " wireDerivatives[CDCWireD::dDirYddeltaXu]" << wireDerivatives[CDCWireD::dDirYddeltaXu] << endl;
               jout << "trackDerivatives[CDCTrackD::dDOCAdDirZ] " << trackDerivatives[CDCTrackD::dDOCAdDirZ] << " wireDerivatives[CDCWireD::dDirZddeltaXu]" << wireDerivatives[CDCWireD::dDirZddeltaXu] << endl;
            }

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

      }
      milleWriter->end();

      japp->RootUnLock();

   }// End loop over tracks

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_MilleFieldOff::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_MilleFieldOff::fini(void)
{
   // Called before program exit after event processing is finished.
   delete milleWriter;
   return NOERROR;
}

