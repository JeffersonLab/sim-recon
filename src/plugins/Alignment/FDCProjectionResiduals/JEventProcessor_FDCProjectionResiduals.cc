// $Id$
//
//    File: JEventProcessor_FDCProjectionResiduals.cc
// Created: Wed Oct 26 14:07:16 EDT 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_FDCProjectionResiduals.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DReferenceTrajectory.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "HistogramTools.h"

using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FDCProjectionResiduals());
}
} // "C"


//------------------
// JEventProcessor_FDCProjectionResiduals (Constructor)
//------------------
JEventProcessor_FDCProjectionResiduals::JEventProcessor_FDCProjectionResiduals()
{

}

//------------------
// ~JEventProcessor_FDCProjectionResiduals (Destructor)
//------------------
JEventProcessor_FDCProjectionResiduals::~JEventProcessor_FDCProjectionResiduals()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FDCProjectionResiduals::init(void)
{
   // This is called once at program startup. 

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FDCProjectionResiduals::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   // This is called whenever the run number changes
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
   dgeom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
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

   typedef map<string,double>::iterator iter_double;
   if (jcalib->Get("CDC/cdc_drift_table::NoBField", tvals)==false){    
      for(unsigned int i=0; i<tvals.size(); i++){
         map<string, double> &row = tvals[i];
         iter_double iter = row.find("t");
         cdc_drift_table.push_back(1000.*iter->second);
      }
   }

   if (jcalib->Get("CDC/drift_parameters::NoBField", tvals)==false){
      map<string, double> &row = tvals[0]; //long drift side
      long_drift_func[0][0]=row["a1"]; 
      long_drift_func[0][1]=row["a2"];
      long_drift_func[0][2]=row["a3"];  
      long_drift_func[1][0]=row["b1"];
      long_drift_func[1][1]=row["b2"];
      long_drift_func[1][2]=row["b3"];
      long_drift_func[2][0]=row["c1"];
      long_drift_func[2][1]=row["c2"];
      long_drift_func[2][2]=row["c3"];

      row = tvals[1]; // short drift side
      short_drift_func[0][0]=row["a1"];
      short_drift_func[0][1]=row["a2"];
      short_drift_func[0][2]=row["a3"];  
      short_drift_func[1][0]=row["b1"];
      short_drift_func[1][1]=row["b2"];
      short_drift_func[1][2]=row["b3"];
      short_drift_func[2][0]=row["c1"];
      short_drift_func[2][1]=row["c2"];
      short_drift_func[2][2]=row["c3"];
   }


   MAX_DRIFT_TIME = 1000.0; //ns: from TRKFIND:MAX_DRIFT_TIME in DTrackCandidate_factory_CDC
   PLANE_TO_SKIP = 0;
   //Make sure it gets initialize first, in case we want to change it:
   vector<const DTrackCandidate*> locTrackCandidates;
   eventLoop->Get(locTrackCandidates,"StraightLine");
   gPARMS->GetParameter("TRKFIT:PLANE_TO_SKIP", PLANE_TO_SKIP);

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FDCProjectionResiduals::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   vector <const DCDCHit *> cdcHitVector;
   loop->Get(cdcHitVector);

   vector <const DChargedTrack *> chargedTrackVector;
   loop->Get(chargedTrackVector);

   unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
      135,135,146,146,158,158,170,170,182,182,197,197,
      209,209};

   for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

      const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

      
      // Cut very loosely on the track quality
      const DTrackTimeBased *thisTimeBasedTrack = nullptr;
      bestHypothesis->GetSingle(thisTimeBasedTrack);

      double t0 = thisTimeBasedTrack->t0();

      Fill1DHistogram("FDCProjectionResiduals", "", "Tracking FOM", thisTimeBasedTrack->FOM, "TrackingFOM", 200, 0.0, 1.0);
      if (thisTimeBasedTrack->FOM < 0.0027 || thisTimeBasedTrack->Ndof < 4) continue;

      // Loop through the pulls to try to check on the FDC calibtrations
      vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
      for (size_t iPull = 0; iPull < pulls.size(); iPull++){
         if ( pulls[iPull].fdc_hit != nullptr){
            if (PLANE_TO_SKIP == 0 || PLANE_TO_SKIP == pulls[iPull].fdc_hit->wire->layer){
               double DOCA = pulls[iPull].d;
               double residual = pulls[iPull].resi;
               double cathode_resi = pulls[iPull].resic;
               double tdrift = pulls[iPull].tdrift;
               Fill2DHistogram("FDCProjectionResiduals", "FDCReco","Distance Vs Time",
                     tdrift, DOCA,
                     "Distance Vs. Time; Time [ns]; Distance [cm]",
                     300, 0.0, 300., 100, 0.0, 0.5);
               Fill2DHistogram("FDCProjectionResiduals", "FDCReco","Residual Vs Time",
                     tdrift, residual,
                     "Residual Vs. Time; Time [ns]; Residual [cm]",
                     200, 0.0, 200., 160, -0.2, 0.2); 
               Fill1DHistogram("FDCProjectionResiduals", "FDCReco","Cathode Residuals",
                     cathode_resi, "Cathode Residual; Residual [cm];", 160,-0.2, 0.2);

            }
         }
      }

      for (auto ringPtr=cdcwires.begin(); ringPtr < cdcwires.end(); ringPtr++){
         vector< DCDCWire * > wireByNumber = (*ringPtr);
         for (auto wirePtr = wireByNumber.begin(); wirePtr < wireByNumber.end(); wirePtr++)
         {
            DCDCWire * wire = *wirePtr;
            //double wireLength = wire->L;
            //double distanceToWire = thisTimeBasedTrack->rt->DistToRT(wire, &wireLength);
            DVector3 POCAOnTrack(0.0, 0.0, 0.0), POCAOnWire(0.0, 0.0, 0.0);
            double zVertex = thisTimeBasedTrack->position().Z();
            double distanceToBeamline = thisTimeBasedTrack->position().Perp();
            double distanceToWire = GetDOCA(wire->origin, wire->udir, thisTimeBasedTrack->position(), thisTimeBasedTrack->momentum(), POCAOnTrack, POCAOnWire);
            double zPOCA = POCAOnTrack.Z();
            DVector3 LOCA = POCAOnTrack - POCAOnWire;

            if(distanceToWire > 1.2 || distanceToBeamline > 1.0 || zPOCA < zVertex || POCAOnWire.Z() > endplate_z) continue;

            double delta = 0.0, dz = 0.0;
            if(!Expect_Hit(thisTimeBasedTrack, wire, distanceToWire, delta, dz))
               continue;
            // Check for a CDC Hit on this wire
            for (auto cdcHit = cdcHitVector.begin(); cdcHit != cdcHitVector.end(); cdcHit++){
               const DCDCHit *thisHit = (*cdcHit);
               if (thisHit->ring == wire->ring && thisHit->straw == wire->straw){
                  // We found the hit on the wire and have all of the information we need.
                  // Get the corrected drift time
                  //DReferenceTrajectory::swim_step_t* swimstep = thisTimeBasedTrack->rt->FindClosestSwimStep(wire);
                  double tdrift = thisHit->t - t0;
                  double measurement = CDCDriftDistance(delta, tdrift);
                  double signedResidual = measurement - (POCAOnTrack.Phi() - POCAOnWire.Phi())*POCAOnTrack.Perp();
                  double residual = measurement - distanceToWire;
                  //jout << "evnt" << eventnumber << " ring " << wire->ring << " straw " << wire->straw << " t " << thisHit->t << " t0 " << t0 << " measurement " << measurement << " residual " << residual << endl;
                  char name[200];
                  char title[200];
                  sprintf(name,"Ring %i Residual Vs. Straw Number", thisHit->ring);
                  sprintf(title,"Ring %i Residual Vs. Straw Number; Straw Number; Residual [cm]", thisHit->ring);
                  Fill2DHistogram("FDCProjectionResiduals","ResidualVsStrawNumber",name,
                        thisHit->straw, residual,
                        title,
                        numstraws[thisHit->ring-1], 0.5, numstraws[thisHit->ring-1] + 0.5, 1000, -0.5, 0.5);
                  sprintf(name,"Ring %i rPhi Residual Vs. phi", thisHit->ring);
                  sprintf(title,"Ring %i #Deltar#phi Vs. #phi; Straw Number; Residual [cm]", thisHit->ring);
                  Fill2DHistogram("FDCProjectionResiduals","ResidualVsPhi",name,
                        thisTimeBasedTrack->momentum().Phi(), signedResidual,
                        title,
                        numstraws[thisHit->ring-1], -3.14, 3.14, 1000, -0.5, 0.5);
                  if (thisHit->ring == 1){
                     sprintf(name,"Ring %i Straw %i Distance Vs. Time", thisHit->ring, thisHit->straw);
                     sprintf(title,"Ring %i Straw %i Distance Vs. Time ; Time [ns]; Distance [cm]", thisHit->ring, thisHit->straw);
                     Fill2DHistogram("FDCProjectionResiduals","DistanceVsTimeRing1",name,
                           tdrift, distanceToWire,
                           title,
                           500, -50.0, 1000, 120, 0.0, 1.2);
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
jerror_t JEventProcessor_FDCProjectionResiduals::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FDCProjectionResiduals::fini(void)
{
   // Called before program exit after event processing is finished.
   return NOERROR;
}

// Convert time to distance for the cdc
double JEventProcessor_FDCProjectionResiduals::CDCDriftDistance(double delta, double t){
   double d=0.;

   if (t>0){
      double f_0=0.;
      double f_delta=0.;

      if (delta > 0){
         double a1=long_drift_func[0][0];
         double a2=long_drift_func[0][1];
         double b1=long_drift_func[1][0];
         double b2=long_drift_func[1][1];
         double c1=long_drift_func[2][0];
         double c2=long_drift_func[2][1];
         double c3=long_drift_func[2][2];

         // use "long side" functional form
         double my_t=0.001*t;
         double sqrt_t=sqrt(my_t);
         double t3=my_t*my_t*my_t;
         double delta_mag=fabs(delta);
         f_delta=(a1+a2*delta_mag)*sqrt_t+(b1+b2*delta_mag)*my_t
            +(c1+c2*delta_mag+c3*delta*delta)*t3;
         f_0=a1*sqrt_t+b1*my_t+c1*t3;
      }
      else{
         double my_t=0.001*t;
         double sqrt_t=sqrt(my_t);
         double delta_mag=fabs(delta);

         // use "short side" functional form
         double a1=short_drift_func[0][0];
         double a2=short_drift_func[0][1];
         double a3=short_drift_func[0][2];
         double b1=short_drift_func[1][0];
         double b2=short_drift_func[1][1];
         double b3=short_drift_func[1][2];

         double delta_sq=delta*delta;
         f_delta= (a1+a2*delta_mag+a3*delta_sq)*sqrt_t
            +(b1+b2*delta_mag+b3*delta_sq)*my_t;
         f_0=a1*sqrt_t+b1*my_t;
      }

      unsigned int max_index=cdc_drift_table.size()-1;
      if (t>cdc_drift_table[max_index]){
         //_DBG_ << "t: " << t <<" d " << f_delta <<endl;
         d=f_delta;
         //cout << "t = " << t << " > " << cdc_drift_table[max_index] << " d = " << f_delta << endl;
         return d;
      }

      // Drift time is within range of table -- interpolate...
      unsigned int index=0;
      index=Locate(cdc_drift_table,t);
      double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
      double frac=(t-cdc_drift_table[index])/dt;
      double d_0=0.01*(double(index)+frac); 

      double P=0.;
      double tcut=250.0; // ns
      if (t<tcut) {
         P=(tcut-t)/tcut;
      }
      d=f_delta*(d_0/f_0*P+1.-P);
   }
   return d;
}

bool JEventProcessor_FDCProjectionResiduals::Expect_Hit(const DTrackTimeBased* thisTimeBasedTrack, DCDCWire* wire, double distanceToWire, double& delta, double& dz)
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

// Locate a position in vector xx given x
unsigned int JEventProcessor_FDCProjectionResiduals::Locate(vector<double>&xx,
      double x){
   int n=xx.size();
   if (x==xx[0]) return 0;
   else if (x==xx[n-1]) return n-2;

   int jl=-1;
   int ju=n;
   int ascnd=(xx[n-1]>=xx[0]);
   while(ju-jl>1){
      int jm=(ju+jl)>>1;
      if ( (x>=xx[jm])==ascnd)
         jl=jm;
      else
         ju=jm;
   } 
   return jl;
}

double JEventProcessor_FDCProjectionResiduals::GetDOCA(DVector3 wirePosition, DVector3 wireDirection, DVector3 trackPosition, DVector3 trackMomentum, DVector3 &POCAOnTrack, DVector3 &POCAOnWire){
   // Get the vector pointing from the wire to the doca point
   Float_t a = trackMomentum.Dot(trackMomentum);
   Float_t b = trackMomentum.Dot(wireDirection);
   Float_t c = wireDirection.Dot(wireDirection);
   DVector3 w0 = trackPosition - wirePosition;
   Float_t d = trackMomentum.Dot(w0);
   Float_t e = wireDirection.Dot(w0);
   Float_t sc = ((b*e - c*d)/(a*c-b*b));
   Float_t tc = ((a*e - b*d)/(a*c-b*b));
   //if (sc < 0) continue; // Track must come from location away from origin
   POCAOnTrack = trackPosition + sc * trackMomentum;
   POCAOnWire  = wirePosition + tc * wireDirection;
   DVector3 LOCA = w0 + sc*trackMomentum - tc*wireDirection;
   return LOCA.Mag();
}
