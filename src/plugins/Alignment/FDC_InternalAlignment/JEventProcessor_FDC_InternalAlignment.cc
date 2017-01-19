// $Id$
//
//    File: JEventProcessor_FDC_InternalAlignment.cc
// Created: Sun Nov 27 16:10:26 EST 2016
// Creator: mstaib (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_FDC_InternalAlignment.h"
using namespace jana;
#include "HistogramTools.h"
#include "FDC/DFDCPseudo.h"

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FDC_InternalAlignment());
}
} // "C"


//------------------
// JEventProcessor_FDC_InternalAlignment (Constructor)
//------------------
JEventProcessor_FDC_InternalAlignment::JEventProcessor_FDC_InternalAlignment()
{

}

//------------------
// ~JEventProcessor_FDC_InternalAlignment (Destructor)
//------------------
JEventProcessor_FDC_InternalAlignment::~JEventProcessor_FDC_InternalAlignment()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FDC_InternalAlignment::init(void)
{
	// This is called once at program startup. 

   gDirectory->mkdir("FDC_InternalAlignment");
   gDirectory->cd("FDC_InternalAlignment");
   for (int i=1; i<=24;i++){

      char thisName[256];
      sprintf(thisName, "Plane %.2i Wire Position Vs XY", i);

      Hist3D[i-1]=new TH3I(thisName, thisName, 100, -50., 50., 100, -50., 50., 100, -0.5, 0.5);
   }
   // We need the constants used for this iteration
   // Use a TProfile to avoid problems adding together multiple root files...
   HistCurrentConstants = new TProfile("CathodeAlignmentConstants", "Constants Used for Cathode Alignment (In CCDB Order)", 96,0.5,96.5);

   gDirectory->cd("..");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FDC_InternalAlignment::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
   // Get curent cathode alignment constants

   JCalibration * jcalib = eventLoop->GetJCalibration();
   vector<map<string,double> >vals;
   if (jcalib->Get("FDC/cathode_alignment",vals)==false){
      for(unsigned int i=0; i<vals.size(); i++){
         map<string,double> &row = vals[i];
         // Get the offsets from the calibration database 
         HistCurrentConstants->Fill(i*4+1,row["dPhiU"]);
         HistCurrentConstants->Fill(i*4+2,row["dU"]);
         HistCurrentConstants->Fill(i*4+3,row["dPhiV"]);
         HistCurrentConstants->Fill(i*4+4,row["dV"]);
      }
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FDC_InternalAlignment::evnt(JEventLoop *loop, uint64_t eventnumber)
{
   vector <const DFDCPseudo*> fdcPseudoVector;
   loop->Get(fdcPseudoVector);

   for (unsigned int i = 0; i<fdcPseudoVector.size(); i++){
      const DFDCPseudo* thisPseudo = fdcPseudoVector[i];
      if (thisPseudo->status != 6) continue;

      japp->RootWriteLock();
      Hist3D[thisPseudo->wire->layer - 1]->Fill(thisPseudo->w, thisPseudo->s, thisPseudo->w_c - thisPseudo->w);
      japp->RootUnLock();

      char thisName[256];
      sprintf(thisName, "Plane %.2i Wire Position Vs XY", thisPseudo->wire->layer);
      Fill2DProfile("FDC_InternalAlignment", "Profile2D", thisName,
            thisPseudo->w, thisPseudo->s, thisPseudo->w_c - thisPseudo->w,
            thisName,
            100, -50.,50., 100,-50., 50.);
   }


   // Test For Lubomir...
   /*
      double constants[192] =  { -0.0043154,   -0.1691660,   0.0007645 ,  -0.0002993 ,  5.0059700 ,  5.0065600 ,  -0.2965900 ,  0.5452100   ,  
      0.3594690,   0.2248910 ,  0.0004884  ,  0.0001668  ,  5.0040600 ,  5.0063800 ,  0.3391100  , 0.5998900    ,
      -0.2362050,   -0.1891460,   0.0005648 ,  -0.0005000 ,  5.0055400 ,  5.0069200 ,  0.5505000  , 0.1661200    ,
      -0.1562530,   -0.1499240,   0.0000257 ,  0.0000318  ,  5.0071500 ,  5.0048300 ,  0.1676100  , -0.2970800   , 
      -0.2040510,   -0.5034870,   -0.0002976,  -0.0009806 ,  5.0056300 ,  5.0060900 ,  -0.4721900 ,  -0.252970   , 
      -0.3032560,   -0.2500840,   0.0002768 ,  -0.0003581 ,  5.0128400 ,  5.0098100 ,  -0.6313500 ,  0.4893600   , 
      0.1563260,   -0.0514443,   0.0001716 ,  0.0002051  ,  5.0078500 ,  5.0084900 ,  -0.3773300 ,  0.2163800   , 
      0.7422930,   0.3868080 ,  -0.0004678 ,  -0.0006898 ,  5.0062200 ,  5.0058600 ,  0.2586500  , 0.0381900    ,
      -0.1862160,   -0.3921980,   -0.0003152,   0.0005399 ,  5.0058100 ,  5.0070700 ,  0.4594900  , 0.3248800    ,
      -0.8467560,   0.0483033 ,  -0.0003287 ,  0.0005756  ,  5.0100700 ,  5.0068000 ,  0.0434000  , -0.2726000   , 
      -0.1340560,   -0.0933582,   0.0012310 ,  0.0000364  ,  5.0065800 ,  5.0076400 ,  -0.3061300 ,  -0.277280   , 
      0.0450700,   -0.0604768,   -0.0010056,   0.0001712 ,  5.0093500 ,  5.0060800 ,  -0.2725000 ,  0.1725300   , 
      -0.2025840,   -0.1141940,   0.0007098 ,  -0.0001474 ,  5.0055900 ,  5.0052500 ,  -0.0582500 ,  0.3434700   , 
      -0.0267204,   -0.0445084,   0.0002086 ,  0.0011202  ,  5.0062000 ,  5.0040900 ,  0.2537100  , 0.2302600    ,
      0.1929130,   -0.2376780,   0.0002160 ,  0.0005313  ,  5.0057000 ,  5.0057500 ,  0.2489000  , 0.0218399    ,
      -0.2891600,   0.1794350 ,   0.0011607 ,  0.0007243  ,  5.0054800 ,  5.0059400 ,  0.1409400  , -0.0889300   , 
      0.3780580,   -0.0605780,   0.0003887 ,  0.0006382  ,  5.0094100 ,  5.0071100 ,  -0.1349100 , 0.0094002    ,
      -0.1119370,   -0.0205372,   0.0003893 ,  0.0000034  ,  5.0099600 ,  5.0073600 ,  0.0673701  , 0.0904398    ,
      -0.1003550,   0.1026090 ,  -0.0003057 ,  0.0000636  ,  5.0070700 ,  5.0071000 ,  -0.1667400 , -0.3741500   , 
      0.4315830,   0.2629620 ,   0.0005938 ,  0.0007167  ,  5.0070100 ,  5.0063000 ,  0.1928000  , -0.1254600   , 
      0.0109732,   0.0740713 ,   0.0006309 ,  0.0000190  ,  5.0058600 ,  5.0063100 ,  -0.2017900 ,  0.1079100   , 
      -0.7136310,   -0.1934810,   0.0006850 ,  0.0000467  ,  5.0024400 ,  5.0048500 ,  -0.0082700 ,  0.4051300   , 
      0.2190650,   -0.0636425,   0.0006076 ,  0.0011785  ,  5.0050100 ,  5.0049200 ,  0.0951700  ,  0.1409100   , 
      -0.0848179,   0.1256750 ,  -0.0000582 ,  -0.0000461 ,  5.0078600 ,  5.0059700 ,  0.2675400  ,  -0.1691500};

      int npseudo=fdcPseudoVector.size();
      for (unsigned int i = 0; i<npseudo; i++){
      const DFDCPseudo* thisPseudo = fdcPseudoVector[i];

      int lpseudo = thisPseudo->wire->layer;
      double zpseudo = thisPseudo->wire->origin.Z();
      double upseudo = (thisPseudo->u)*10.;
      double vpseudo = (thisPseudo->v)*10.;
      double spseudo = (thisPseudo->s)*10.;
      double wpseudo = (thisPseudo->w)*10.;

      int layer=lpseudo;

   // get constants for this layer

   float pxwl=wpseudo;
   float pycl=-spseudo;
   double u0 = constants[(layer-1)*8+0];
   double d0 = constants[(layer-1)*8+1];
   double udel = constants[(layer-1)*8+2];
   double ddel = constants[(layer-1)*8+3];
   double strip_pitch_u = constants[(layer-1)*8+4];
   double strip_pitch_d = constants[(layer-1)*8+5];
   double xshift = constants[(layer-1)*8+6];
   double yshift = constants[(layer-1)*8+7];

   float v=vpseudo*strip_pitch_u/5.005-u0;
   float u=upseudo*strip_pitch_d/5.005-d0;
   float myangle=15.*M_PI/180.;
   float strip_angle=75.*M_PI/180.;
   float phi_u=strip_angle+ddel;
   float phi_v=M_PI-strip_angle+udel;
   float cosPhiU=cos(phi_u);
   float cosPhiV=cos(phi_v);
   float sinPhiU=sin(phi_u);
   float sinPhiV=sin(phi_v);
   float pxcl=(u*sinPhiV-v*sinPhiU)/(cosPhiV*sinPhiU-cosPhiU*sinPhiV);//-xshift;
   pycl=(u*cosPhiV-v*cosPhiU)/(cosPhiU*sinPhiV-cosPhiV*sinPhiU);//+yshift;

   // The wire position projected from the wire
   japp->RootWriteLock();
   Hist3D[layer-1]->Fill(pxwl/10.,pycl/10.,(pxwl-pxcl)/10.); 
   japp->RootUnLock();

   char thisName[256];
   sprintf(thisName, "Plane %.2i Wire Position Vs XY", layer);
   Fill2DProfile("FDC_InternalAlignment", "Profile2D", thisName,
         pxwl/10., pycl/10., (pxwl-pxcl)/10.,
         thisName,
         100, -50.,50., 100,-50., 50.);
}
*/
return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FDC_InternalAlignment::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FDC_InternalAlignment::fini(void)
{
   // Called before program exit after event processing is finished.
   return NOERROR;
}

