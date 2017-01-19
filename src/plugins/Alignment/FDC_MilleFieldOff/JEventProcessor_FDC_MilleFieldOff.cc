// $Id$
//
//    File: JEventProcessor_FDC_MilleFieldOff.cc
// Created: Tue Dec 20 17:43:35 Local time zone must be set--see zic manual page 2016
// Creator: mstaib (on Linux egbert 2.6.32-642.6.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_FDC_MilleFieldOff.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "TRACKING/DTrackCandidate.h"
#include "FDC/DFDCPseudo.h"
#include "TDirectory.h"

using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FDC_MilleFieldOff());
}
} // "C"


//------------------
// JEventProcessor_FDC_MilleFieldOff (Constructor)
//------------------
JEventProcessor_FDC_MilleFieldOff::JEventProcessor_FDC_MilleFieldOff()
{

}

//------------------
// ~JEventProcessor_FDC_MilleFieldOff (Destructor)
//------------------
JEventProcessor_FDC_MilleFieldOff::~JEventProcessor_FDC_MilleFieldOff()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FDC_MilleFieldOff::init(void)
{
	// This is called once at program startup. 
   milleWriter = new Mille("fdc_mille_out.mil");

   gDirectory->mkdir("AlignmentConstants");
   gDirectory->cd("AlignmentConstants");
   // We need the constants used for this iteration
   // Use a TProfile to avoid problems adding together multiple root files...
   HistCurrentConstants = new TProfile("FDCAlignmentConstants", "Constants Used for FDC Alignment (In MILLEPEDE Order)", 26000 ,0.5, 26000.5);

   gDirectory->cd("..");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FDC_MilleFieldOff::brun(JEventLoop *eventLoop, int32_t runnumber)
{
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
   if (jcalib->Get("FDC/cathode_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstants->Fill((i+1)*1000+101,row["dU"]);
         HistCurrentConstants->Fill((i+1)*1000+102,row["dV"]);
         HistCurrentConstants->Fill((i+1)*1000+103,row["dPhiU"]);
         HistCurrentConstants->Fill((i+1)*1000+104,row["dPhiV"]);
      }
   }

   if (jcalib->Get("FDC/cell_offsets",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstants->Fill((i+1)*1000+1,row["xshift"]);
         HistCurrentConstants->Fill((i+1)*1000+100,row["yshift"]);
      }
   }

   if (jcalib->Get("FDC/wire_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database
         HistCurrentConstants->Fill((i+1)*1000+2,row["dPhi"]);
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
               HistCurrentConstants->Fill((i/2+1)*1000+600+j,row[name]);
            }
            else {// U
               HistCurrentConstants->Fill((i/2+1)*1000+300+j,row[name]);
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
               HistCurrentConstants->Fill((i/2+7)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstants->Fill((i/2+7)*1000+300+j,row[name]);
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
               HistCurrentConstants->Fill((i/2+13)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstants->Fill((i/2+13)*1000+300+j,row[name]);
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
               HistCurrentConstants->Fill((i/2+19)*1000+600+j,row[name]);
            }
            else // U
               HistCurrentConstants->Fill((i/2+19)*1000+300+j,row[name]);
         }
      }
   }

   if (jcalib->Get("FDC/strip_pitches_v2",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         HistCurrentConstants->Fill((i+1)*1000 + 200, row["U_SP_1"]);
         HistCurrentConstants->Fill((i+1)*1000 + 201, row["U_G_1"]);
         HistCurrentConstants->Fill((i+1)*1000 + 202, row["U_SP_2"]);
         HistCurrentConstants->Fill((i+1)*1000 + 203, row["U_G_2"]);
         HistCurrentConstants->Fill((i+1)*1000 + 204, row["U_SP_3"]);
         HistCurrentConstants->Fill((i+1)*1000 + 205, row["V_SP_1"]);
         HistCurrentConstants->Fill((i+1)*1000 + 206, row["V_G_1"]);
         HistCurrentConstants->Fill((i+1)*1000 + 207, row["V_SP_2"]);
         HistCurrentConstants->Fill((i+1)*1000 + 208, row["V_G_2"]);
         HistCurrentConstants->Fill((i+1)*1000 + 209, row["V_SP_3"]);
      }
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FDC_MilleFieldOff::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   // Loop over the tracks, get the tracking pulls
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
         //const DCDCTrackHit *cdc_hit = pulls[iPull].cdc_hit;
         const DFDCPseudo *fdc_hit   = pulls[iPull].fdc_hit;
         if (fdc_hit == NULL) continue;
         //double docaphi              = pulls[iPull].docaphi; // phi of doca in CDC straws
         //double z                    = pulls[iPull].z;// z position at doca
         //double tcorr                = pulls[iPull].tcorr; // drift time with correction for B
         double resic                = pulls[iPull].resic; // residual for FDC cathode measuremtns
         double errc                 = pulls[iPull].errc;

         if (fdc_hit->status != 6) continue;

         vector<double> trackDerivatives = pulls[iPull].trackDerivatives;
         vector<double> stateVector = pulls[iPull].stateVector;

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

         localDerW[0]=trackDerivatives[0]; localDerW[1]=trackDerivatives[1];
         localDerW[2]=trackDerivatives[2]; localDerW[3]=trackDerivatives[3];

         // need some things from the state vector
         double x = stateVector[0]; double y = stateVector[1];
         double tx = stateVector[2]; double ty = stateVector[3];
         double cosa=thisHit->wire->udir.y();
         double sina=thisHit->wire->udir.x();
         double w=thisHit->w;
         double wpred=x*cosa-y*sina;
         double tu=tx*cosa-ty*sina;
         double alpha=atan(tu);
         double cosalpha=cos(alpha);

         globalDerW[0] = -cosalpha; labelW[0] = layerOffset + 1;
         globalDerW[1] = -cosalpha*(y*cosa+x*sina)+(ty*cosa+tx*sina)*(w-wpred)/(1+tu*tu); labelW[1] = layerOffset + 2;

         milleWriter->mille(NLC, localDerW, NGL_W, globalDerW, labelW, resi, err);
         // Now for the cathode measurement, there are more alignment parameters.
         const int NGL_C = 21; // Number of global parameters for wire alignment
         float localDerC[NLC];
         float globalDerC[NGL_C];
         int labelC[NGL_C];

         localDerC[0]=trackDerivatives[4]; localDerC[1]=trackDerivatives[5];
         localDerC[2]=trackDerivatives[6]; localDerC[3]=trackDerivatives[7];

         globalDerC[0] = -1.0;
         labelC[0] = layerOffset + 100; // delta_y

         // Cathode U and V offsets
         globalDerC[1] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU];    labelC[1] = layerOffset + 101;
         globalDerC[2] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV];    labelC[2] = layerOffset + 102;
         globalDerC[3] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaPhiU]; labelC[3] = layerOffset + 103;
         globalDerC[4] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaPhiV]; labelC[4] = layerOffset + 104;

         // Strip Pitch Calibration
         globalDerC[5] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[0]; labelC[5] = layerOffset + 200;
         globalDerC[6] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[1]; labelC[6] = layerOffset + 201;
         globalDerC[7] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[2]; labelC[7] = layerOffset + 202;
         globalDerC[8] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[3]; labelC[8] = layerOffset + 203;
         globalDerC[9] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripPitchDerivatives[4]; labelC[9] = layerOffset + 204;
         globalDerC[10] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[5]; labelC[10] = layerOffset + 205;
         globalDerC[11] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[6]; labelC[11] = layerOffset + 206;
         globalDerC[12] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[7]; labelC[12] = layerOffset + 207;
         globalDerC[13] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[8]; labelC[13] = layerOffset + 208;
         globalDerC[14] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripPitchDerivatives[9]; labelC[14] = layerOffset + 209;

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

         globalDerC[15] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[0]; labelC[15] = layerOffset + 300 + gainLabels[0];
         globalDerC[16] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[1]; labelC[16] = layerOffset + 300 + gainLabels[1];
         globalDerC[17] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaU]*fdcStripGainDerivatives[2]; labelC[17] = layerOffset + 300 + gainLabels[2];
         globalDerC[18] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[3]; labelC[18] = layerOffset + 600 + gainLabels[3];
         globalDerC[19] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[4]; labelC[19] = layerOffset + 600 + gainLabels[4];
         globalDerC[20] = -pseudoAlignmentDerivatives[FDCPseudoD::dSddeltaV]*fdcStripGainDerivatives[5]; labelC[20] = layerOffset + 600 + gainLabels[5];

         milleWriter->mille(NLC, localDerC, NGL_C, globalDerC, labelC, resic, errc);

      }
      milleWriter->end();

      japp->RootUnLock();
   } 
   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FDC_MilleFieldOff::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FDC_MilleFieldOff::fini(void)
{
   // Called before program exit after event processing is finished.
   delete milleWriter; // Closes Mille output file
   return NOERROR;
}

