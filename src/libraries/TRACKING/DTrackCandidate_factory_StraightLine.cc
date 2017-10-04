// $Id$
//
//    File: DTrackCandidate_factory_StraightLine.cc
// Created: Fri Aug 15 09:14:04 EDT 2014
// Creator: staylor (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <map>

#include "DTrackCandidate_factory_StraightLine.h"
using namespace jana;
#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
#include <DANA/DApplication.h>
#include <BCAL/DBCALShower.h>
#include "HDGEOMETRY/DGeometry.h"

bool DTrackCandidate_StraightLine_cdc_hit_cmp(const DCDCTrackHit *a,
      const DCDCTrackHit *b){

   return(a->wire->origin.Y()>b->wire->origin.Y());
}

bool DTrackCandidate_StraightLine_cdc_hit_reverse_cmp(const DCDCTrackHit *a,
      const DCDCTrackHit *b){

   return(a->wire->origin.Y()<b->wire->origin.Y());
}

bool DTrackCandidate_StraightLine_cdc_hit_radius_cmp(const DCDCTrackHit *a,
      const DCDCTrackHit *b){
   if (a==NULL || b==NULL){
      cout << "Null pointer in CDC hit list??" << endl;
      return false;
   }
   const DCDCWire *wire_a= a->wire;
   const DCDCWire *wire_b= b->wire;
   if(wire_a->ring == wire_b->ring){
      return wire_a->straw < wire_b->straw;
   }

   return (wire_a->ring<wire_b->ring);
}


bool DTrackCandidate_StraightLine_fdc_hit_cmp(const DFDCPseudo *a,
      const DFDCPseudo *b){

   return(a->wire->origin.z()<b->wire->origin.z());
}

// parametrization of time-to-distance for FDC
double DTrackCandidate_factory_StraightLine::fdc_drift_distance(double time){
   if (time<0.) return 0.;
   double tsq=time*time;
   double d=DRIFT_FUNC_PARMS[0]*sqrt(time)+DRIFT_FUNC_PARMS[1]*time
      +DRIFT_FUNC_PARMS[2]*tsq+DRIFT_FUNC_PARMS[3]*tsq*time;

   return d;
}

// Crude approximation for the variance in drift distance due to smearing
double DTrackCandidate_factory_StraightLine::fdc_drift_variance(double t){
   //return FDC_ANODE_VARIANCE;
   if (t<5.) t=5.;
   double sigma=DRIFT_RES_PARMS[0]/(t+1.)+DRIFT_RES_PARMS[1]+DRIFT_RES_PARMS[2]*t*t;

   return sigma*sigma;
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_StraightLine::init(void)
{
   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_StraightLine::brun(jana::JEventLoop *loop, int runnumber)
{
   DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());

   double dz;
   const DGeometry *geom = dapp->GetDGeometry(runnumber);
   geom->GetCDCEndplate(cdc_endplate_z,dz,cdc_endplate_rmin,cdc_endplate_rmax);
   geom->GetCDCAxialLength(cdc_length);

   JCalibration *jcalib = dapp->GetJCalibration(runnumber);
   size_t found=jcalib->GetContext().find("variation=mc");
   isMC = found != string::npos;

   const char *ccdbRequest;
   vector< map<string, double> > tvals;
   cdc_drift_table.clear();

   ccdbRequest="CDC/cdc_drift_table::NoBField";  // The "NoBField" table is used in simulation

   if (jcalib->Get(ccdbRequest, tvals)==false){    
      for(unsigned int i=0; i<tvals.size(); i++){
         map<string, double> &row = tvals[i];
         cdc_drift_table.push_back(1000.*row["t"]);
      }
   }
   else{
      jerr << " CDC time-to-distance table not available... bailing..." << endl;
      exit(0);
   }

   map<string, double> cdc_res_parms;
   ccdbRequest="CDC/cdc_resolution_parms_v2::NoBField"; // There is a difference between errors in Field on and off
   jcalib->Get(ccdbRequest, cdc_res_parms);
   CDC_RES_PAR1 = cdc_res_parms["res_par1"];
   CDC_RES_PAR2 = cdc_res_parms["res_par2"];
   CDC_RES_PAR3 = cdc_res_parms["res_par3"];

   // Get the straw sag parameters from the database
   max_sag.clear();
   sag_phi_offset.clear();
   unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
      135,135,146,146,158,158,170,170,182,182,197,197,
      209,209};
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

   if (isMC) ccdbRequest="CDC/drift_parameters";
   else ccdbRequest="CDC/drift_parameters::NoBField";

   if (jcalib->Get(ccdbRequest, tvals)==false){
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

   // Time-to-distance function parameters for FDC
   map<string,double>drift_func_parms;
   if (isMC) ccdbRequest="FDC/drift_function_parms";
   else ccdbRequest="FDC/drift_function_parms::NoBField";
   jcalib->Get(ccdbRequest,drift_func_parms); 
   DRIFT_FUNC_PARMS[0]=drift_func_parms["p0"];   
   DRIFT_FUNC_PARMS[1]=drift_func_parms["p1"];
   DRIFT_FUNC_PARMS[2]=drift_func_parms["p2"]; 
   DRIFT_FUNC_PARMS[3]=drift_func_parms["p3"];

   // Parameters for accounting for variation in drift distance from FDC
   map<string,double>drift_res_parms;
   if (isMC) ccdbRequest="FDC/drift_resolution_parms";
   else ccdbRequest="FDC/drift_resolution_parms::NoBField";
   jcalib->Get(ccdbRequest,drift_res_parms); 
   DRIFT_RES_PARMS[0]=drift_res_parms["p0"];   
   DRIFT_RES_PARMS[1]=drift_res_parms["p1"];
   DRIFT_RES_PARMS[2]=drift_res_parms["p2"]; 

   COSMICS=false;
   gPARMS->SetDefaultParameter("TRKFIND:COSMICS",COSMICS);

   CHI2CUT = 15.0; 
   gPARMS->SetDefaultParameter("TRKFIT:CHI2CUT",CHI2CUT);    

   DO_PRUNING = 1;
   gPARMS->SetDefaultParameter("TRKFIT:DO_PRUNING",DO_PRUNING);

   DEBUG_HISTS=false;
   gPARMS->SetDefaultParameter("TRKFIND:DEBUG_HISTS",DEBUG_HISTS);


   USE_FDC_DRIFT_TIMES=true;
   gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_DRIFT_TIMES",
         USE_FDC_DRIFT_TIMES);

   PLANE_TO_SKIP=0;
   gPARMS->SetDefaultParameter("TRKFIT:PLANE_TO_SKIP",PLANE_TO_SKIP);

   SKIP_CDC=false;
   gPARMS->SetDefaultParameter("TRKFIT:SKIP_CDC",SKIP_CDC);

   SKIP_FDC=false;
   gPARMS->SetDefaultParameter("TRKFIT:SKIP_FDC",SKIP_FDC);

   VERBOSE=0;
   gPARMS->SetDefaultParameter("TRKFIT:VERBOSE",VERBOSE);

   CDC_MATCH_DOCA=0.78;
   gPARMS->SetDefaultParameter("TRKFIT:CDC_MATCH_DOCA",CDC_MATCH_DOCA);

   // Get pointer to TrackFinder object 
   vector<const DTrackFinder *> finders;
   eventLoop->Get(finders);

   if(finders.size()<1){
      _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
      return RESOURCE_UNAVAILABLE;
   }

   // Drop the const qualifier from the DTrackFinder pointer
   finder = const_cast<DTrackFinder*>(finders[0]);

   if (DEBUG_HISTS){
      Hvres=(TH2F *)gROOT->FindObject("Hvres");
      if (!Hvres) Hvres=new TH2F("Hvres","Residual along wire",100,-0.25,0.25,24,0.5,24.5);
      hFDCOccTrkFind=new TH1I("Occ form track finding", "Occ per plane", 24,0.5,24.5);
      hFDCOccTrkFit=new TH1I("Occ form track fitting", "Occ per plane", 24,0.5,24.5);
      hFDCOccTrkSmooth=new TH1I("Occ form track smoothing", "Occ per plane", 24,0.5,24.5);
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_StraightLine::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   // Reset the track finder
   finder->Reset();

   vector<const DCDCTrackHit*>cdcs;
   vector<const DFDCPseudo*>pseudos;
   loop->Get(cdcs);
   loop->Get(pseudos);

   set<unsigned int> used_cdc;

   // Look for tracks in the FDC.
   if(!SKIP_FDC){
      if (pseudos.size()>4){
         for (size_t i=0;i<pseudos.size();i++) finder->AddHit(pseudos[i]);
         finder->FindFDCSegments();
         finder->LinkFDCSegments();

         // Get the list of linked segments and fit the hits to lines
         const vector<DTrackFinder::fdc_segment_t>tracks=finder->GetFDCTracks();
         for (size_t i=0;i<tracks.size();i++){
            // list of FDC hits
            vector<const DFDCPseudo *>hits=tracks[i].hits;
            if(DEBUG_HISTS){
               for (size_t j=0; j< hits.size(); j++){
                  hFDCOccTrkFind->Fill(hits[j]->wire->layer);
               }
            }
            sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_fdc_hit_cmp);

            // Initial guess for state vector
            DMatrix4x1 S(tracks[i].S);

            // Move x and y to just before the first hit
            double my_z=hits[0]->wire->origin.z()-1.;
            S(state_x)+=my_z*S(state_tx);
            S(state_y)+=my_z*S(state_ty);

            // Use earliest fdc time to estimate t0
            double t0=1e6;
            double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
            for (unsigned int m=0;m<hits.size();m++){
               if (hits[m]->time<t0){
                  double L=(hits[m]->wire->origin.z()-my_z)*dsdz;
                  t0=hits[m]->time-L/29.98; // assume moving at speed of light
               }
            }


            //Run the Kalman Filter algorithm
            DoFilter(t0,my_z,S,hits,cdcs,used_cdc);
         }
      }
   }

   if(!SKIP_CDC){
      if (cdcs.size()>4){
         for (size_t i=0;i<cdcs.size();i++) {
            // If the CDC hit had not been grabbed by the FDC fit, try to find CDC only tracks.
            if(used_cdc.find(i) == used_cdc.end()) finder->AddHit(cdcs[i]);
         }
         finder->FindAxialSegments();
         finder->LinkCDCSegments();

         // Get the list of linked segments and fit the hits to lines
         const vector<DTrackFinder::cdc_track_t>tracks=finder->GetCDCTracks();
         if (VERBOSE > 0) jout << "Looping over " << tracks.size() << " CDC tracks..." << endl;
         for (size_t i=0;i<tracks.size();i++){

            // start z position and direction of propagation (default = +z direction)
            double z0=tracks[i].z,dzsign=1.;
            // Shift z0 towars target to limit LR ambiguity issues
            if(!COSMICS) z0 = z0 - (z0-65.0)/2;

            // Initial guess for state vector
            DMatrix4x1 S(tracks[i].S);

            // list of axial and stereo hits for this track
            vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
            hits.insert(hits.end(),tracks[i].stereo_hits.begin(),
                  tracks[i].stereo_hits.end());

            if (COSMICS){
               // Step track to front plate.
               if (VERBOSE) {
                  jout<< " S before step. z = " << z0 << endl;
                  S.Print();
               }
               S(state_x)-=(z0-(cdc_endplate_z-cdc_length))*S(state_tx);
               S(state_y)-=(z0-(cdc_endplate_z-cdc_length))*S(state_ty);
               z0=cdc_endplate_z-cdc_length;
               if (VERBOSE){
                  jout << " S after step z = " << z0 << endl;
               }

               if (S(state_ty) > 0.) sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_cdc_hit_reverse_cmp);
               else sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_cdc_hit_cmp);
            }
            else{	
               DVector3 pos,origin,dir(0,0,1.);
               finder->FindDoca(z0,S,dir,origin,&pos);
               S(state_x)=pos.x();
               S(state_y)=pos.y();
               if (z0<pos.z()) dzsign=-1.;
               z0=pos.z();

               sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_cdc_hit_radius_cmp);

            }

            // Use earliest cdc time to estimate t0
            double t0=1e6;
            for (unsigned int j=0;j<hits.size();j++){
               double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
               double t_test=hits[j]->tdrift-L/29.98;
               if (t_test<t0) t0=t_test;
            }

            // Run the Kalman Filter algorithm
            DoFilter(t0,z0,S,hits,dzsign);	 
         }
      }
   }

   // Set CDC ring & FDC plane hit patterns
   for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
   {
      vector<const DCDCTrackHit*> locCDCTrackHits;
      _data[loc_i]->Get(locCDCTrackHits);

      vector<const DFDCPseudo*> locFDCPseudos;
      _data[loc_i]->Get(locFDCPseudos);

      _data[loc_i]->dCDCRings = dParticleID->Get_CDCRingBitPattern(locCDCTrackHits);
      _data[loc_i]->dFDCPlanes = dParticleID->Get_FDCPlaneBitPattern(locFDCPseudos);
   }

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_StraightLine::erun(void)
{
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_StraightLine::fini(void)
{
   return NOERROR;
}


// Steering routine for the kalman filter
jerror_t 
DTrackCandidate_factory_StraightLine::DoFilter(double t0,double OuterZ,
      DMatrix4x1 &S,
      vector<const DCDCTrackHit *>&hits,
      double dzsign){
   unsigned int numhits=hits.size();
   unsigned int maxindex=numhits-1;

   // vectors of indexes to cdc hits used in the fit
   vector<int> used_cdc_hits(numhits);
   vector<int> used_cdc_hits_best_fit(numhits);

   // vectors of residual information 
   vector<cdc_update_t>updates(numhits);
   vector<cdc_update_t>best_updates(numhits);

   // deque to store reference trajectory
   deque<trajectory_t>trajectory;
   deque<trajectory_t>best_trajectory;

   // State vector to store "best" values
   DMatrix4x1 Sbest;

   // Covariance matrix
   DMatrix4x4 C0,C,Cbest;
   C0(state_x,state_x)=C0(state_y,state_y)=10.0;     
   C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.1;

   double chi2=1e16,chi2_old=1e16;
   unsigned int ndof=0,ndof_old=0;
   unsigned int iter=0;
   double z0=OuterZ;

   // Perform a wire-based pass
   for(iter=0;iter<20;iter++){
      if (VERBOSE) jout << " Performing Wire Based Pass iter " << iter << endl;
      chi2_old=chi2; 
      ndof_old=ndof;

      if(!COSMICS){
         DVector3 pos,origin,dir(0,0,1.);
         finder->FindDoca(z0,S,dir,origin,&pos);
         S(state_x)=pos.x();
         S(state_y)=pos.y();
         z0=pos.z();
      }

      trajectory.clear();
      if (SetReferenceTrajectory(t0,z0,S,trajectory,hits[maxindex],dzsign)
            !=NOERROR) break;

      if (VERBOSE) jout << " Reference Trajectory Set " << endl;
      C=C0;
      if (KalmanFilter(S,C,hits,used_cdc_hits,trajectory,updates,chi2,ndof,false,iter)!=NOERROR) break;
      if (VERBOSE) jout << " Wire Based Filter returns NOERROR" << endl;
      if (fabs(chi2_old-chi2)<0.1 || chi2>chi2_old) break;  

      // Save the current state and covariance matrixes
      Cbest=C;
      Sbest=S; 
   }
   if (iter>0){
      // Perform a time-based pass
      S=Sbest;
      chi2=1e16;
      ndof=0;

      for (iter=0;iter<20;iter++){
         if (VERBOSE) jout << " Performing Time Based Pass iter " << iter << endl;
         chi2_old=chi2; 
         ndof_old=ndof;
         if (!COSMICS){
            DVector3 pos,origin,dir(0,0,1.);
            finder->FindDoca(z0,S,dir,origin,&pos);
            S(state_x)=pos.x();
            S(state_y)=pos.y();
            z0=pos.z();
         }
         trajectory.clear();
         if (SetReferenceTrajectory(t0,z0,S,trajectory,hits[maxindex],dzsign)
               ==NOERROR){
            if (VERBOSE) jout << " Set Reference Trajectory" << endl;
            C=C0;
            if (KalmanFilter(S,C,hits,used_cdc_hits,trajectory,updates,chi2,ndof,true,iter)!=NOERROR) break;
            if (VERBOSE) jout << " Fit Succeeded chi2 " << chi2 << " prob " << TMath::Prob(chi2,ndof) << " chi2_old " << chi2_old << " prob " << TMath::Prob(chi2_old,ndof_old) << endl;

            //printf("chi2 %f %f\n",chi2_old,chi2);

            if (fabs(chi2-chi2_old)<0.1  
                  || TMath::Prob(chi2,ndof)<TMath::Prob(chi2_old,ndof_old)) break;

            Sbest=S;
            Cbest=C;

            used_cdc_hits_best_fit=used_cdc_hits;
            best_updates=updates;
            best_trajectory=trajectory;
         }
         else {
            if (VERBOSE) jout << " Set Reference Trajectory Failed" << endl;
            break;
         }
      }
      if (iter>0 && trajectory.size()>1){
         // Create a new track candidate
         if (VERBOSE) jout << " Method converged on time based iter " << iter << endl;
         DTrackCandidate *cand = new DTrackCandidate;

         double sign=1.;
         unsigned int last_index=trajectory.size()-1;
         if (COSMICS==false){
            DVector3 pos,origin,dir(0,0,1.); 
            finder->FindDoca(trajectory[last_index].z,Sbest,dir,origin,&pos);
            cand->setPosition(pos);

            if (trajectory[0].z<pos.z()) sign=-1.;
         }
         else{ 
            cand->setPosition(DVector3(Sbest(state_x),Sbest(state_y),
                     trajectory[last_index].z));
         }

         double tx=Sbest(state_tx),ty=Sbest(state_ty);
         double phi=atan2(ty,tx);
         if (sign<0) phi+=M_PI;
         if (VERBOSE) jout << " phi " << phi << endl;
         double tanl=sign/sqrt(tx*tx+ty*ty);
         double pt=10.*cos(atan(tanl));      
         cand->setMomentum(DVector3(pt*cos(phi),pt*sin(phi),pt*tanl));

         cand->Ndof=ndof_old;
         cand->chisq=chi2_old;
         cand->setPID(PiPlus);

         // Smooth the result
         if (Smooth(best_trajectory,best_updates,hits,cand) == NOERROR) cand->IsSmoothed=true;
         if(VERBOSE){
            if (cand->IsSmoothed) jout << "Smooth Success" << endl;
            else jout << " Smooth Failed!!!" << endl;
         }

         // Add hits used in the fit as associated objects and add best pull 
         // vector to the candidate
         for (unsigned int k=0;k<used_cdc_hits_best_fit.size();k++){
            if (used_cdc_hits_best_fit[k]==1){
               cand->AddAssociatedObject(hits[k]);
            }
         }


         _data.push_back(cand);
      }
   } // Check that the wire-based fit did not fail


   return NOERROR;
}


//Reference trajectory for the track for cdc tracks
jerror_t DTrackCandidate_factory_StraightLine
::SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
      deque<trajectory_t>&trajectory,
      const DCDCTrackHit *last_cdc,double &dzsign){ 
   DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

   double ds=1.0;
   double t=t0;

   // last y position of hit (approximate, using center of wire)
   double last_y=last_cdc->wire->origin.y();
   double last_r2=last_cdc->wire->origin.Perp2();
   // Check that track points towards last wire, otherwise swap deltaz
   DVector3 dir(S(state_tx),S(state_ty),dzsign);
   if (!COSMICS){
      double dphi=dir.Phi()-last_cdc->wire->origin.Phi(); 
      while (dphi>M_PI) dphi-=2*M_PI;
      while (dphi<-M_PI) dphi+=2*M_PI;
      if (fabs(dphi) > M_PI/2.) dzsign*=-1.;
   }
   if (fabs(dir.Theta() - M_PI/2.) < 0.2) ds = 0.1;

   //jout << "dPhi " << dphi << " theta " << dir.Theta() << endl;
   double dz=dzsign*ds/sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

   if (VERBOSE) {
      if (COSMICS) jout << " Swimming Reference Trajectory last CDC y " << last_y << " dz "<< dz << endl;
      else jout << " Swimming Reference Trajectory last CDC r2 " << last_r2 << " dz "<< dz << endl;
      jout << " S" << endl; S.Print(); jout << "z= "<< z << endl;
      jout << "  Last CDC ring " << last_cdc->wire->ring << " straw " << last_cdc->wire->straw << endl;
   }
   unsigned int numsteps=0;
   const unsigned int MAX_STEPS=5000;
   double upstreamEndplate = cdc_endplate_z - cdc_length;
   bool done=false;
   do{
      numsteps++;
      z+=dz;
      J(state_x,state_tx)=-dz;
      J(state_y,state_ty)=-dz;
      // Flight time: assume particle is moving at the speed of light
      t+=ds/29.98;
      //propagate the state to the next z position
      S(state_x)+=S(state_tx)*dz;
      S(state_y)+=S(state_ty)*dz;
      if (z > cdc_endplate_z && dz < 0) continue;
      if (z < upstreamEndplate && dz > 0) continue;
      trajectory.push_front(trajectory_t(z,t,S,J,DMatrix4x1(),DMatrix4x4()));

      if (COSMICS) done=((z>cdc_endplate_z) | (z<cdc_endplate_z-cdc_length));
      else{
         double r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
         if (VERBOSE >= 5) jout << " r2=" << r2 << endl;
         done=(r2>last_r2) | (z>cdc_endplate_z) | (z<cdc_endplate_z-cdc_length); 
      }
   }while (!done && numsteps<MAX_STEPS);

   if (VERBOSE)
   {
      if (VERBOSE > 10){
         printf("Trajectory:\n");
         for (unsigned int i=0;i<trajectory.size();i++){
            printf(" x %f y %f z %f\n",trajectory[i].S(state_x),
                  trajectory[i].S(state_y),trajectory[i].z); 
         }
      }
      else{
         printf("%i step trajectory Begin/End:\n", numsteps);
         printf(" x %f y %f z %f\n",trajectory[0].S(state_x), trajectory[0].S(state_y), trajectory[0].z);
         if (trajectory.size() > 1) printf(" x %f y %f z %f\n",trajectory[trajectory.size()-1].S(state_x),
               trajectory[trajectory.size()-1].S(state_y),trajectory[trajectory.size()-1].z);
      }
   }
   if (trajectory.size()<2) return UNRECOVERABLE_ERROR;
   return NOERROR;
}
// Perform the Kalman Filter for the current set of cdc hits
jerror_t 
DTrackCandidate_factory_StraightLine::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
      vector<const DCDCTrackHit *>&hits,
      vector<int>&used_hits,
      deque<trajectory_t>&trajectory,
      vector<cdc_update_t>&updates,
      double &chi2,unsigned int &ndof,
      bool timebased,
      unsigned int iter){
   DMatrix1x4 H;  // Track projection matrix
   DMatrix4x1 H_T; // Transpose of track projection matrix 
   DMatrix4x1 K;  // Kalman gain matrix
   DMatrix4x4 I; // identity matrix
   DMatrix4x4 J; // Jacobian matrix
   DMatrix4x1 S0; // State vector from reference trajectory
   double V=1.15*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor

   const double d_EPS=1e-8;

   // Zero out the vector of used hit flags
   for (unsigned int i=0;i<used_hits.size();i++) used_hits[i]=0;

   //Initialize chi2 and ndof
   chi2=0.;
   ndof=0;

   //Save the starting values for C and S
   trajectory[0].Skk=S;
   trajectory[0].Ckk=C;

   double doca2=0.;

   // CDC index and wire position variables
   unsigned int cdc_index=hits.size()-1;
   bool more_hits=true;
   const DCDCWire *wire=hits[cdc_index]->wire;
   DVector3 origin=wire->origin;
   double z0=origin.z();
   double vz=wire->udir.z();
   if (VERBOSE) jout << " Starting in Ring " << wire->ring << endl;
   DVector3 wdir=(1./vz)*wire->udir;
   DVector3 wirepos=origin+(trajectory[0].z-z0)*wdir;

   /// compute initial doca^2 to first wire
   double dx=S(state_x)-wirepos.X();
   double dy=S(state_y)-wirepos.Y();
   double old_doca2=dx*dx+dy*dy;

   // Loop over all steps in the trajectory
   S0=trajectory[0].S;
   J=trajectory[0].J;
   for (unsigned int k=1;k<trajectory.size();k++){
      if (!C.IsPosDef()) return UNRECOVERABLE_ERROR;

      // Propagate the state and covariance matrix forward in z
      S=trajectory[k].S+J*(S-S0);
      C=J*C*J.Transpose();

      // Save the current state and covariance matrix in the deque
      trajectory[k].Skk=S;
      trajectory[k].Ckk=C;

      // Save S and J for the next step
      S0=trajectory[k].S;
      J=trajectory[k].J;

      // Position along wire
      wirepos=origin+(trajectory[k].z-z0)*wdir;

      // New doca^2
      dx=S(state_x)-wirepos.X();
      dy=S(state_y)-wirepos.Y();
      doca2=dx*dx+dy*dy;
      if (VERBOSE > 10) jout<< "At Position " << S(state_x) << " " << S(state_y) << " " << trajectory[k].z << " doca2 " << doca2 << endl;

      if (doca2>old_doca2 && more_hits){

         // zero-position and direction of line describing particle trajectory
         double tx=S(state_tx),ty=S(state_ty);
         DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
         DVector3 tdir(tx,ty,1.);

         // Find the true doca to the wire
         DVector3 diff=pos0-origin;
         double dx0=diff.x(),dy0=diff.y();
         double wdir_dot_diff=diff.Dot(wdir);
         double tdir_dot_diff=diff.Dot(tdir);
         double tdir_dot_wdir=tdir.Dot(wdir);
         double tdir2=tdir.Mag2();
         double wdir2=wdir.Mag2();
         double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
         double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
         double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
         double scale=1./D;
         double s=scale*N;
         double t=scale*N1;
         diff+=s*tdir-t*wdir;
         double d=diff.Mag()+d_EPS; // prevent division by zero

         // The next measurement and its variance
         double tdrift=hits[cdc_index]->tdrift-trajectory[k].t;
         double dmeas=0.39; 
         double delta=0.0;
         if (timebased){
            V=CDCDriftVariance(tdrift);	

            double phi_d=diff.Phi();
            double dphi=phi_d-origin.Phi();  
            while (dphi>M_PI) dphi-=2*M_PI;
            while (dphi<-M_PI) dphi+=2*M_PI;

            int ring_index=hits[cdc_index]->wire->ring-1;
            int straw_index=hits[cdc_index]->wire->straw-1;
            double dz=t*wdir.z();
            delta=max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
               *cos(phi_d+sag_phi_offset[ring_index][straw_index]);
            dmeas=CDCDriftDistance(dphi,delta,tdrift);
         }

         // residual
         double res=dmeas-d;
         if (VERBOSE>5) jout << " Residual " << res << endl;

         // Track projection

         double one_over_d=1./d;
         double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
         double wx=wdir.x(),wy=wdir.y();

         double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
         double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
         double dtdtx=scale*(dN1dtx-t*dDdtx);

         double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
         double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
         double dtdty=scale*(dN1dty-t*dDdty);

         double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
         double dsdtx=scale*(dNdtx-s*dDdtx);

         double dNdty=wy*wdir_dot_diff-wdir2*dy0;
         double dsdty=scale*(dNdty-s*dDdty);

         H(state_tx)=H_T(state_tx)
            =one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
                  +diffz*(dsdtx-dtdtx));
         H(state_ty)=H_T(state_ty)
            =one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
                  +diffz*(dsdty-dtdty));

         double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
         double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
         double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
         double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);

         H(state_x)=H_T(state_x)
            =one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
                  +diffz*(dsdx-dtdx));
         H(state_y)=H_T(state_y)
            =one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
                  +diffz*(dsdy-dtdy));

         /*
            DVector3 WirePosNew = wirepos+s*wdir;
            double dy = S(state_y)+S(state_ty)*s-WirePosNew.Y();
            double dx = S(state_x)+S(state_tx)*s-WirePosNew.X();
            double cosstereo=cos(hits[cdc_index]->wire->stereo);
            double thisd = sqrt(dx*dx+dy*dy)*cosstereo+d_EPS;
            double cosstereo2_over_d=cosstereo*cosstereo/thisd;
            H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
            H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;     
            H(state_tx)=H_T(state_tx)=s*H(state_x);
            H(state_ty)=H_T(state_ty)=s*H(state_y);
            */ 
         double InvV=1./(V+H*C*H_T);

         // Check how far this hit is from the projection
         double chi2check=res*res*InvV;
         if (chi2check < CHI2CUT || DO_PRUNING == 0 || (COSMICS && iter == 0)){
            if (VERBOSE>5) jout << "Hit Added to track " << endl;
            // Compute Kalman gain matrix
            K=InvV*(C*H_T);

            // Update state vector covariance matrix
            DMatrix4x4 Ctest=C-K*(H*C);

            //C.Print();
            //K.Print();
            //Ctest.Print();

            // Check that Ctest is positive definite
            if (!Ctest.IsPosDef()) return VALUE_OUT_OF_RANGE;
            C=Ctest;
            if(VERBOSE>10) C.Print();
            // Update the state vector 
            //S=S+res*K;
            S+=res*K;
            if(VERBOSE>10) S.Print();

            // Compute new residual 
            //d=finder->FindDoca(trajectory[k].z,S,wdir,origin);
            res=res-H*K*res;

            // Update chi2 
            double fit_V=V-H*C*H_T;
            chi2+=res*res/fit_V;
            ndof++;

            // Flag that we used this hit
            used_hits[cdc_index]=1;

            // fill updates
            updates[cdc_index].resi=res;
            updates[cdc_index].d=d;
            updates[cdc_index].delta=delta;
            updates[cdc_index].S=S;
            updates[cdc_index].C=C;
            updates[cdc_index].V=V;
            updates[cdc_index].tdrift=tdrift;
            updates[cdc_index].ddrift=dmeas;
            updates[cdc_index].s=29.98*trajectory[k].t; // assume beta=1
            trajectory[k].id=cdc_index+1;
         }
         // move to next cdc hit
         if (cdc_index>0){
            cdc_index--;

            //New wire position
            wire=hits[cdc_index]->wire;
            if (VERBOSE>5) {
               jout << " Next Wire ring " << wire->ring << " straw " << wire->straw << " origin udir" << endl;
               wire->origin.Print(); wire->udir.Print();
            }
            origin=wire->origin;
            z0=origin.z();
            vz=wire->udir.z();
            wdir=(1./vz)*wire->udir;
            wirepos=origin+((trajectory[k].z-z0))*wdir;

            // New doca^2
            dx=S(state_x)-wirepos.x();
            dy=S(state_y)-wirepos.y();
            doca2=dx*dx+dy*dy;	

         }
         else more_hits=false;
      }

      old_doca2=doca2;
   }
   if (ndof<=4) return VALUE_OUT_OF_RANGE;

   ndof-=4;

   return NOERROR;
}

// Smooth the CDC only tracks
jerror_t
DTrackCandidate_factory_StraightLine::Smooth(deque<trajectory_t>&trajectory,
      vector<cdc_update_t>&cdc_updates,
      vector<const DCDCTrackHit *>&hits,
      DTrackCandidate *cand){
   unsigned int max=trajectory.size()-1;
   DMatrix4x1 S=(trajectory[max].Skk);
   DMatrix4x4 C=(trajectory[max].Ckk);
   DMatrix4x4 JT=trajectory[max].J.Transpose();
   DMatrix4x1 Ss=S;
   DMatrix4x4 Cs=C;
   DMatrix4x4 A,dC;
   DMatrix1x4 H;  // Track projection matrix
   DMatrix4x1 H_T; // Transpose of track projection matrix

   const double d_EPS=1e-8;

   for (unsigned int m=max-1;m>0;m--){
      if (trajectory[m].id>0){
         unsigned int id=trajectory[m].id-1;
         A=cdc_updates[id].C*JT*C.Invert();
         Ss=cdc_updates[id].S+A*(Ss-S);

         dC=A*(Cs-C)*A.Transpose();
         Cs=cdc_updates[id].C+dC;
         if (VERBOSE > 10) {
            jout << " In Smoothing Step Using ID " << id << "/" << cdc_updates.size() << " for ring " << hits[id]->wire->ring << endl;
            jout << " A cdc_updates[id].C Ss Cs " << endl;
            A.Print(); cdc_updates[id].C.Print(); Ss.Print(); Cs.Print();
         }
         if(!Cs.IsPosDef()) {
            if (VERBOSE) jout << "Cs is not PosDef!" << endl;
            return VALUE_OUT_OF_RANGE;
         }

         const DCDCWire *wire=hits[id]->wire;
         DVector3 origin=wire->origin;
         double z0=origin.z();
         double vz=wire->udir.z();
         DVector3 wdir=(1./vz)*wire->udir;
         DVector3 wirepos=origin+(trajectory[m].z-z0)*wdir;
         // Position and direction from state vector
         double x=Ss(state_x);
         double y=Ss(state_y);
         double tx=Ss(state_tx);
         double ty=Ss(state_ty);

         DVector3 pos0(x,y,trajectory[m].z);
         DVector3 tdir(tx,ty,1.);

         // Find the true doca to the wire
         DVector3 diff=pos0-origin;
         double dx0=diff.x(),dy0=diff.y();
         double wdir_dot_diff=diff.Dot(wdir);
         double tdir_dot_diff=diff.Dot(tdir);
         double tdir_dot_wdir=tdir.Dot(wdir);
         double tdir2=tdir.Mag2();
         double wdir2=wdir.Mag2();
         double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
         double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
         double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
         double scale=1./D;
         double s=scale*N;
         double t=scale*N1;
         diff+=s*tdir-t*wdir;
         double d=diff.Mag()+d_EPS; // prevent division by zero
         double ddrift = cdc_updates[id].ddrift;

         double resi = ddrift - d;
         /* 
            DVector3 WirePosNew = wirepos+s*wdir;
            double dy = Ss(state_y)+Ss(state_ty)*s-WirePosNew.Y();
            double dx = Ss(state_x)+Ss(state_tx)*s-WirePosNew.X();
            double cosstereo=cos(hits[id]->wire->stereo);
            double thisd = sqrt(dx*dx+dy*dy)*cosstereo+d_EPS;
            double cosstereo2_over_d=cosstereo*cosstereo/thisd;
            H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
            H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;
            H(state_tx)=H_T(state_tx)=s*H(state_x);
            H(state_ty)=H_T(state_ty)=s*H(state_y);
            */  

         // Track projection

         {
            double one_over_d=1./d;
            double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
            double wx=wdir.x(),wy=wdir.y();

            double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
            double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
            double dtdtx=scale*(dN1dtx-t*dDdtx);

            double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
            double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
            double dtdty=scale*(dN1dty-t*dDdty);

            double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
            double dsdtx=scale*(dNdtx-s*dDdtx);

            double dNdty=wy*wdir_dot_diff-wdir2*dy0;
            double dsdty=scale*(dNdty-s*dDdty);

            H(state_tx)=H_T(state_tx)
               =one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
                     +diffz*(dsdtx-dtdtx));
            H(state_ty)=H_T(state_ty)
               =one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
                     +diffz*(dsdty-dtdty));

            double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
            double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
            double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
            double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);

            H(state_x)=H_T(state_x)
               =one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
                     +diffz*(dsdx-dtdx));
            H(state_y)=H_T(state_y)
               =one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
                     +diffz*(dsdy-dtdy));
         }
         double V=cdc_updates[id].V;


         if (VERBOSE > 10) jout << " d " << d << " H*S " << H*S << endl;
         V=V-H*Cs*H_T;

         if (V<0) return VALUE_OUT_OF_RANGE;

         // Add the pull
         DTrackFitter::pull_t thisPull(resi,sqrt(V),
               trajectory[m].t*SPEED_OF_LIGHT,
               cdc_updates[id].tdrift,
               d,
               hits[id], NULL,
               diff.Phi(), //docaphi
               trajectory[m].z,
               cdc_updates[id].tdrift);

         // Derivatives for alignment
         double wtx=wire->udir.X(), wty=wire->udir.Y(), wtz=wire->udir.Z();
         double wx=wire->origin.X(), wy=wire->origin.Y(), wz=wire->origin.Z();

         double z=trajectory[m].z;
         double tx2=tx*tx, ty2=ty*ty;
         double wtx2=wtx*wtx, wty2=wty*wty, wtz2=wtz*wtz;
         double denom=(1 + ty2)*wtx2 + (1 + tx2)*wty2 - 2*ty*wty*wtz + (tx2 + ty2)*wtz2 - 2*tx*wtx*(ty*wty + wtz) +d_EPS;
         double denom2=denom*denom;
         double c1=-(wtx - tx*wtz)*(wy - y);
         double c2=wty*(wx - tx*wz - x + tx*z);
         double c3=ty*(-(wtz*wx) + wtx*wz + wtz*x - wtx*z);
         double dscale=0.5*(1./d);

         vector<double> derivatives(11);

         derivatives[CDCTrackD::dDOCAdOriginX]=dscale*(2*(wty - ty*wtz)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdOriginY]=dscale*(2*(-wtx + tx*wtz)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdOriginZ]=dscale*(2*(ty*wtx - tx*wty)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdDirX]=dscale*(2*(wty - ty*wtz)*(c1 + c2 + c3)*
               (tx*(ty*wty + wtz)*(wx - x) + (wty - ty*wtz)*(-wy + y + ty*(wz - z)) +
                wtx*(-((1 + ty2)*wx) + (1 + ty2)*x + tx*(ty*wy + wz - ty*y - z)) + tx2*(wty*(-wy + y) + wtz*(-wz + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdDirY]=dscale*(-2*(wtx - tx*wtz)*(c1 + c2 + c3)*
               (tx*(ty*wty + wtz)*(wx - x) + (wty - ty*wtz)*(-wy + y + ty*(wz - z)) +
                wtx*(-((1 + ty2)*wx) + (1 + ty2)*x + tx*(ty*wy + wz - ty*y - z)) + tx2*(wty*(-wy + y) + wtz*(-wz + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdDirZ]=dscale*(-2*(ty*wtx - tx*wty)*(c1 + c2 + c3)*
               (-(tx*(ty*wty + wtz)*(wx - x)) + tx2*(wty*(wy - y) + wtz*(wz - z)) + (wty - ty*wtz)*(wy - y + ty*(-wz + z)) +
                wtx*((1 + ty2)*wx - (1 + ty2)*x + tx*(-(ty*wy) - wz + ty*y + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdS0]=-derivatives[CDCTrackD::dDOCAdOriginX];

         derivatives[CDCTrackD::dDOCAdS1]=-derivatives[CDCTrackD::dDOCAdOriginY];

         derivatives[CDCTrackD::dDOCAdS2]=dscale*(2*(wty - ty*wtz)*(-c1 - c2 - c3)*
               (-(wtx*wtz*wx) - wty*wtz*wy + wtx2*wz + wty2*wz + wtx*wtz*x + wty*wtz*y - wtx2*z - wty2*z +
                tx*(wty2*(wx - x) + wtx*wty*(-wy + y) + wtz*(wtz*wx - wtx*wz - wtz*x + wtx*z)) +
                ty*(wtx*wty*(-wx + x) + wtx2*(wy - y) + wtz*(wtz*wy - wty*wz - wtz*y + wty*z))))/denom2;

         derivatives[CDCTrackD::dDOCAdS3]=dscale*(2*(wtx - tx*wtz)*(c1 + c2 + c3)*
               (-(wtx*wtz*wx) - wty*wtz*wy + wtx2*wz + wty2*wz + wtx*wtz*x + wty*wtz*y - wtx2*z - wty2*z +
                tx*(wty2*(wx - x) + wtx*wty*(-wy + y) + wtz*(wtz*wx - wtx*wz - wtz*x + wtx*z)) +
                ty*(wtx*wty*(-wx + x) + wtx2*(wy - y) + wtz*(wtz*wy - wty*wz - wtz*y + wty*z))))/denom2;

         thisPull.AddTrackDerivatives(derivatives);

         cand->pulls.push_back(thisPull);

      }
      else{
         A=trajectory[m].Ckk*JT*C.Invert();
         Ss=trajectory[m].Skk+A*(Ss-S);
         Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
      }

      S=trajectory[m].Skk;
      C=trajectory[m].Ckk;
      JT=trajectory[m].J.Transpose();
   }

   return NOERROR;
}

// Locate a position in vector xx given x
unsigned int DTrackCandidate_factory_StraightLine::Locate(vector<double>&xx,
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


// Convert time to distance for the cdc
double DTrackCandidate_factory_StraightLine::CDCDriftDistance(double t){
   double d=0.;
   if (t>cdc_drift_table[cdc_drift_table.size()-1]) return 0.78;
   if (t>0){
      unsigned int index=0;
      index=Locate(cdc_drift_table,t);
      double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
      double frac=(t-cdc_drift_table[index])/dt;
      d=0.01*(double(index)+frac); 
   }
   return d;
}

// Convert time to distance for the cdc
double DTrackCandidate_factory_StraightLine::CDCDriftDistance(double dphi, 
      double delta,double t){
   double d=0.;
   if (t>0){
      double f_0=0.;
      double f_delta=0.;

      if (delta>0){
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

         return d;
      }

      // Drift time is within range of table -- interpolate...
      unsigned int index=0;
      index=Locate(cdc_drift_table,t);
      double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
      double frac=(t-cdc_drift_table[index])/dt;
      double d_0=0.01*(double(index)+frac); 

      if (isMC) return d_0;

      double P=0.;
      double tcut=250.0; // ns
      if (t<tcut) {
         P=(tcut-t)/tcut;
      }
      d=f_delta*(d_0/f_0*P+1.-P);
   }
   return d;
}




// Smearing function derived from fitting residuals
inline double DTrackCandidate_factory_StraightLine::CDCDriftVariance(double t){ 
   if (t<0.) t=0.;
   double sigma=CDC_RES_PAR1/(t+1.)+CDC_RES_PAR2 + CDC_RES_PAR3*t;
   return sigma*sigma;
}


// Steering routine for the kalman filter
jerror_t 
DTrackCandidate_factory_StraightLine::DoFilter(double t0,double start_z,
      DMatrix4x1 &S,
      vector<const DFDCPseudo *>&hits,
      vector<const DCDCTrackHit *>&cdc_hits,
      set<unsigned int> &used_cdc_hits){
   // vectors of indexes to fdc hits used in the fit
   unsigned int numhits=hits.size();
   vector<int> used_fdc_hits(numhits);
   vector<int> used_fdc_hits_best_fit(numhits);

   // Best guess for state vector at the beginning of the trajectory
   DMatrix4x1 Sbest;

   // Use the result from the initial line fit to form a reference trajectory 
   // for the track. 
   deque<trajectory_t>trajectory;
   deque<trajectory_t>best_trajectory;

   // vectors of residual information 
   vector<fdc_update_t>updates(numhits);
   vector<fdc_update_t>best_updates(numhits);
   vector<cdc_update_t>cdc_updates;
   vector<cdc_update_t>best_cdc_updates;

   vector<const DCDCTrackHit *> matchedCDCHits;

   // Intial guess for covariance matrix
   DMatrix4x4 C,C0,Cbest;
   C0(state_x,state_x)=C0(state_y,state_y)=1.;
   C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;

   // Chi-squared and degrees of freedom
   double chi2=1e16,chi2_old=1e16;
   unsigned int ndof=0,ndof_old=0;
   unsigned iter=0;
   // First pass
   for(iter=0;iter<20;iter++){
      chi2_old=chi2; 
      ndof_old=ndof;

      trajectory.clear();
      if (SetReferenceTrajectory(t0,start_z,S,trajectory,hits)!=NOERROR) break;

      C=C0;
      if (KalmanFilter(S,C,hits,used_fdc_hits,matchedCDCHits,trajectory,updates,cdc_updates,chi2,ndof
               )!=NOERROR) break;

      // printf(" == iter %d =====chi2 %f ndof %d \n",iter,chi2,ndof);
      if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1) break;  

      // Save the current state and covariance matrixes
      Cbest=C;
      Sbest=S;

      used_fdc_hits_best_fit=used_fdc_hits;
      best_trajectory=trajectory;
      best_updates=updates;
   }

   // Take these best fit values and try to grab CDC hits that may be associated with this track.
   if (iter>0 && trajectory.size()>1 && !SKIP_CDC){

      // Get intersection of track with CDC endplate.
      double tx=Sbest(state_tx),ty=Sbest(state_ty);
      double phi=atan2(ty,tx);
      double tanl=1./sqrt(tx*tx+ty*ty);
      double pt=10.*cos(atan(tanl));
      DVector3 mom(pt*cos(phi),pt*sin(phi),pt*tanl);

      unsigned int last_index=trajectory.size()-1;
      DVector3 pos,origin,dir(0,0,1.);
      double z=trajectory[last_index].z;
      finder->FindDoca(z,Sbest,dir,origin,&pos);

      DVector3 norm(0,0,1);
      DVector3 pointInPlane(0.,0.,cdc_endplate_z);
      DVector3 intersection;
      finder->FindIntersectionWithPlane(pointInPlane,norm,
            pos,mom,intersection);

      double intersectionR=intersection.Perp();
      if(intersectionR > cdc_endplate_rmin && intersectionR < cdc_endplate_rmax){ // We might have some CDC hits
         for (size_t i=0; i<cdc_hits.size(); i++){
            origin = cdc_hits[i]->wire->origin;
            if (origin.Perp() > intersectionR) continue; // Assume the track is coming from the target
            dir = cdc_hits[i]->wire->udir;
            DVector3 pos;
            double doca = finder->FindDoca(z,Sbest,dir,origin,&pos);     
            if (doca < CDC_MATCH_DOCA){
               if (VERBOSE) jout << " Matched CDC hit R" << cdc_hits[i]->wire->ring << " S" << cdc_hits[i]->wire->straw << "to FDC track " << endl;
               matchedCDCHits.push_back(cdc_hits[i]);
               used_cdc_hits.insert(i);
            }
         }
      }
      if (matchedCDCHits.size()  > 0) { 
         if (VERBOSE) jout << matchedCDCHits.size() << " CDC hits have been found for this FDC track" << endl;
         cdc_updates.resize(matchedCDCHits.size());
         best_cdc_updates.resize(matchedCDCHits.size());
         sort(matchedCDCHits.begin(),matchedCDCHits.end(),DTrackCandidate_StraightLine_cdc_hit_radius_cmp);

         // Perform fit including the information from the CDC hits
         // Chi-squared and degrees of freedom
         chi2=1e16;chi2_old=1e16;
         ndof=0;ndof_old=0;
         // Pass including CDC hit information
         for(iter=0;iter<20;iter++){
            last_index=trajectory.size()-1;
            z=trajectory[last_index].z;
            DVector3 pos,origin,dir(0,0,1.);
            finder->FindDoca(z,S,dir,origin,&pos);
            if (VERBOSE) {jout << " Finding DOCA z=" << z << " S DOCA " << endl; S.Print(); pos.Print();}
            S(state_x)=pos.x();
            S(state_y)=pos.y();
            double z0=pos.z();

            // Use earliest fdc time to estimate t0
            t0=1e6;
            double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
            for (unsigned int m=0;m<hits.size();m++){
               if (hits[m]->time<t0){
                  double L=(hits[m]->wire->origin.z()-z0)*dsdz;
                  t0=hits[m]->time-L/29.98; // assume moving at speed of light
               }
            }

            if (VERBOSE) jout << " FDC/CDC RT starting at z=" << z0 << " t0=" << t0 << endl;;
            chi2_old=chi2;
            ndof_old=ndof;

            trajectory.clear();
            if (SetReferenceTrajectory(t0,z0,S,trajectory,hits)!=NOERROR) break;

            C=C0;
            if (KalmanFilter(S,C,hits,used_fdc_hits,matchedCDCHits,trajectory,updates,cdc_updates,chi2,ndof
                     )!=NOERROR) break;

            // printf(" == iter %d =====chi2 %f ndof %d \n",iter,chi2,ndof);
            if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1) break;

            // Save the current state and covariance matrixes
            Cbest=C;
            Sbest=S;

            used_fdc_hits_best_fit=used_fdc_hits;
            best_trajectory=trajectory;
            best_updates=updates;
            best_cdc_updates=cdc_updates;
         }

      }
   }


   if (iter>0 && trajectory.size()>1){  
      // Create new track candidate
      DTrackCandidate *cand = new DTrackCandidate;

      double tx=Sbest(state_tx),ty=Sbest(state_ty);
      double phi=atan2(ty,tx);
      double tanl=1./sqrt(tx*tx+ty*ty);
      double pt=10.*cos(atan(tanl));    
      cand->setMomentum(DVector3(pt*cos(phi),pt*sin(phi),pt*tanl));

      unsigned int last_index=trajectory.size()-1;
      DVector3 pos,origin,dir(0,0,1.);
      double z=trajectory[last_index].z;
      finder->FindDoca(z,Sbest,dir,origin,&pos);
      cand->setPosition(pos);

      // Run the smoother 
      if (Smooth(best_trajectory,best_updates,hits,best_cdc_updates,matchedCDCHits,cand) == NOERROR) cand->IsSmoothed=true;

      for (unsigned int k=0;k<used_fdc_hits_best_fit.size();k++){
         if (used_fdc_hits_best_fit[k]==1){
            cand->AddAssociatedObject(hits[k]);
         }
      }

      if (DEBUG_HISTS){
         for (unsigned int id=0;id<hits.size();id++){	  
            double cospsi=hits[id]->wire->udir.y();
            double sinpsi=hits[id]->wire->udir.x();

            DVector3 norm(0,0,1);
            DVector3 intersection;
            finder->FindIntersectionWithPlane(hits[id]->wire->origin,norm,
                  pos,cand->momentum(),intersection);
            // To transform from (x,y) to (u,v), need to do a rotation:
            double v = intersection.y()*cospsi+intersection.x()*sinpsi;

            Hvres->Fill(v-hits[id]->s,hits[id]->wire->layer);

         }
      }


      cand->Ndof=ndof_old;
      cand->chisq=chi2_old;
      cand->setPID(PiMinus);

      _data.push_back(cand);

   }


   return NOERROR;
}


// Reference trajectory for the track
jerror_t 
DTrackCandidate_factory_StraightLine::SetReferenceTrajectory(double t0,double z,
      DMatrix4x1 &S,
      deque<trajectory_t>&trajectory,
      vector<const DFDCPseudo *>&pseudos){
   const double EPS=1e-3;

   // Jacobian matrix 
   DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

   double dz=1.1;
   double t=t0;
   trajectory.push_front(trajectory_t(z,t,S,J,DMatrix4x1(),DMatrix4x4()));

   double zhit=z;
   double old_zhit=z;
   for (unsigned int i=0;i<pseudos.size();i++){  
      zhit=pseudos[i]->wire->origin.z();
      dz=1.1;

      if (fabs(zhit-old_zhit)<EPS && i>0){
         trajectory[0].numhits++;
         continue;
      }
      // propagate until we would step beyond the FDC hit plane
      bool done=false;
      while (!done){	    
         double new_z=z+dz;	      

         if (new_z>zhit){
            dz=zhit-z;
            new_z=zhit;
            done=true;
         }
         J(state_x,state_tx)=-dz;
         J(state_y,state_ty)=-dz;
         // Flight time: assume particle is moving at the speed of light
         t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;
         //propagate the state to the next z position
         S(state_x)+=S(state_tx)*dz;
         S(state_y)+=S(state_ty)*dz;


         trajectory.push_front(trajectory_t(new_z,t,S,J,DMatrix4x1(),
                  DMatrix4x4())); 
         if (done){
            trajectory[0].id=i+1;
            trajectory[0].numhits=1;
         }

         z=new_z;
      }	   
      old_zhit=zhit;
   }
   // One last step
   dz=1.1;
   J(state_x,state_tx)=-dz;
   J(state_y,state_ty)=-dz;

   // Flight time: assume particle is moving at the speed of light
   t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;

   //propagate the state to the next z position
   S(state_x)+=S(state_tx)*dz;
   S(state_y)+=S(state_ty)*dz;
   trajectory.push_front(trajectory_t(z+dz,t,S,J,DMatrix4x1(),DMatrix4x4()));

   if (false)
   {
      printf("Trajectory:\n");
      for (unsigned int i=0;i<trajectory.size();i++){
         printf(" x %f y %f z %f first hit %d num in layer %d\n",trajectory[i].S(state_x),
               trajectory[i].S(state_y),trajectory[i].z,trajectory[i].id,
               trajectory[i].numhits); 
      }
   }

   return NOERROR;
}

// Perform Kalman Filter for the current trajectory
jerror_t 
DTrackCandidate_factory_StraightLine::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
      vector<const DFDCPseudo *>&hits,
      vector<int>&used_hits,
      vector<const DCDCTrackHit *>&cdc_hits,
      deque<trajectory_t>&trajectory,
      vector<fdc_update_t>&updates,
      vector<cdc_update_t>&cdc_updates,
      double &chi2,unsigned int &ndof){
   DMatrix2x4 H;  // Track projection matrix
   DMatrix4x2 H_T; // Transpose of track projection matrix 
   DMatrix4x2 K;  // Kalman gain matrix
   DMatrix2x2 V(0.0833,0.,0.,0.000256);  // Measurement variance 
   DMatrix2x2 Vtemp,InvV;
   DMatrix2x1 Mdiff;
   DMatrix4x4 I; // identity matrix
   DMatrix4x4 J; // Jacobian matrix
   DMatrix4x1 S0; // State vector from reference trajectory

   DMatrix1x4 H_CDC;  // Track projection matrix
   DMatrix4x1 H_T_CDC; // Transpose of track projection matrix
   DMatrix4x1 K_CDC;  // Kalman gain matrix
   double V_CDC;

   const double d_EPS=1e-8;

   // Zero out the vector of used hit flags
   for (unsigned int i=0;i<used_hits.size();i++) used_hits[i]=0;

   //Initialize chi2 and ndof
   chi2=0.;
   ndof=0;

   // Save the starting values for C and S in the deque
   trajectory[0].Skk=S;
   trajectory[0].Ckk=C;

   // Loop over all steps in the trajectory
   S0=trajectory[0].S;
   J=trajectory[0].J;

   // CDC index and wire position variables
   bool more_hits = cdc_hits.size() == 0 ? false: true;
   bool firstCDCStep=true;
   unsigned int cdc_index=0;
   const DCDCWire *wire;
   DVector3 origin,wdir,wirepos;
   double doca2=0.0, old_doca2=0.0;
   if(more_hits){
      cdc_index=cdc_hits.size()-1;
      wire=cdc_hits[cdc_index]->wire;
      origin=wire->origin;
      double vz=wire->udir.z();
      if (VERBOSE) jout << " Additional CDC Hits in FDC track Starting in Ring " << wire->ring << endl;
      wdir=(1./vz)*wire->udir;
   }

   for (unsigned int k=1;k<trajectory.size();k++){
      if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
         return UNRECOVERABLE_ERROR;

      // Propagate the state and covariance matrix forward in z
      S=trajectory[k].S+J*(S-S0);
      C=J*C*J.Transpose();

      // Save the current state and covariance matrix in the deque
      trajectory[k].Skk=S;
      trajectory[k].Ckk=C;

      // Save S and J for the next step
      S0=trajectory[k].S;
      J=trajectory[k].J;

      // Correct S and C for the hit 
      if (trajectory[k].id>0){
         unsigned int id=trajectory[k].id-1;

         double cospsi=cos(hits[id]->wire->angle);
         double sinpsi=sin(hits[id]->wire->angle);

         // State vector
         double x=S(state_x);
         double y=S(state_y);
         double tx=S(state_tx);
         double ty=S(state_ty);

         // Small angle alignment correction
         x = x + hits[id]->wire->angles.Z()*y;
         y = y - hits[id]->wire->angles.Z()*x;
         //tz = 1. + my_fdchits[id]->phiY*tx - my_fdchits[id]->phiX*ty;
         tx = (tx + hits[id]->wire->angles.Z()*ty - hits[id]->wire->angles.Y());
         ty = (ty - hits[id]->wire->angles.Z()*tx + hits[id]->wire->angles.X());

         if (std::isnan(x) || std::isnan(y)) return UNRECOVERABLE_ERROR;

         // x,y and tx,ty in local coordinate system	
         // To transform from (x,y) to (u,v), need to do a rotation:
         //   u = x*cos(psi)-y*sin(psi)
         //   v = y*cos(psi)+x*sin(psi)
         // (without alignment offsets)
         double vpred_wire_plane=y*cospsi+x*sinpsi;
         double upred_wire_plane=x*cospsi-y*sinpsi;
         double tu=tx*cospsi-ty*sinpsi;
         double tv=tx*sinpsi+ty*cospsi;

         // Variables for angle of incidence with respect to the z-direction in
         // the u-z plane
         double alpha=atan(tu);
         double cosalpha=cos(alpha);
         double cos2_alpha=cosalpha*cosalpha;
         double sinalpha=sin(alpha);
         double sin2_alpha=sinalpha*sinalpha;
         double cos2_alpha_minus_sin2_alpha=cos2_alpha-sin2_alpha;

         // Difference between measurement and projection
         for (int m=trajectory[k].numhits-1;m>=0;m--){
            unsigned int my_id=id+m;
            double uwire=hits[my_id]->w;
            // (signed) distance of closest approach to wire
            double du=upred_wire_plane-uwire;
            double doca=du*cosalpha;

            // Predicted avalanche position along the wire
            double vpred=vpred_wire_plane;

            // Measured position of hit along wire
            double v=hits[my_id]->s; 

            // Difference between measurements and predictions
            double drift=0.; // assume hit at wire position
            if (USE_FDC_DRIFT_TIMES){
               double drift_time=hits[my_id]->time-trajectory[k].t; 
               drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time);

               V(0,0)=fdc_drift_variance(drift_time);
            }
            Mdiff(0)=drift-doca;
            Mdiff(1)=v-vpred;

            // Matrix for transforming from state-vector space to measurement space
            H_T(state_x,0)=cospsi*cosalpha;
            H_T(state_y,0)=-sinpsi*cosalpha;
            double temp=-du*sinalpha*cos2_alpha;
            H_T(state_tx,0)=cospsi*temp;
            H_T(state_ty,0)=-sinpsi*temp;
            double temp2=cosalpha*sinalpha*tv;
            H_T(state_x,1)=sinpsi-temp2*cospsi;
            H_T(state_y,1)=cospsi+temp2*sinpsi;
            double temp4=sinalpha*doca;
            double temp5=tv*cos2_alpha*du*cos2_alpha_minus_sin2_alpha;
            H_T(state_tx,1)=-sinpsi*temp4-cospsi*temp5;
            H_T(state_ty,1)=-cospsi*temp4+sinpsi*temp5;

            // Matrix transpose H_T -> H
            H(0,state_x)=H_T(state_x,0);
            H(0,state_y)=H_T(state_y,0);
            H(0,state_tx)=H_T(state_tx,0);
            H(0,state_ty)=H_T(state_ty,0);
            H(1,state_x)=H_T(state_x,1);
            H(1,state_y)=H_T(state_y,1);
            H(1,state_tx)=H_T(state_tx,1);
            H(1,state_ty)=H_T(state_ty,1);

            // Variance for this hit
            InvV=(V+H*C*H_T).Invert();

            // Compute Kalman gain matrix
            K=(C*H_T)*InvV;

            if (hits[my_id]->wire->layer!=PLANE_TO_SKIP){        	
               if(DEBUG_HISTS){
                  hFDCOccTrkFit->Fill(hits[my_id]->wire->layer);
               }
               // Update the state vector 
               S+=K*Mdiff;
               if(VERBOSE) S.Print();
               // Update state vector covariance matrix
               C=C-K*(H*C);    

               // Update the filtered measurement covariane matrix and put results in 
               // update vector
               DMatrix2x2 RC=V-H*C*H_T;
               DMatrix2x1 res=Mdiff-H*K*Mdiff;

               chi2+=RC.Chi2(res);
               ndof+=2;

               // fill pull vector entries
               updates[my_id].V=RC;
            }
            else{
               updates[my_id].V=V;
            }

            used_hits[my_id]=1;

            // fill pull vector
            updates[my_id].d=doca;
            updates[my_id].S=S;
            updates[my_id].C=C;
            updates[my_id].tdrift=hits[my_id]->time-trajectory[k].t;
            updates[my_id].s=29.98*trajectory[k].t; // assume beta=1
         } 
      }

      if (more_hits && trajectory[k].z < cdc_endplate_z){
         // Position along wire
         double z0=origin.Z();
         wirepos=origin+(trajectory[k].z-z0)*wdir;

         // New doca^2
         double dx=S(state_x)-wirepos.X();
         double dy=S(state_y)-wirepos.Y();
         doca2=dx*dx+dy*dy;
         if (VERBOSE > 10) jout<< "At Position " << S(state_x) << " " << S(state_y) << " " << trajectory[k].z << " doca2 " << doca2 << endl;

         if (doca2>old_doca2 && more_hits && !firstCDCStep){

            // zero-position and direction of line describing particle trajectory
            double tx=S(state_tx),ty=S(state_ty);
            DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
            DVector3 tdir(tx,ty,1.);

            // Find the true doca to the wire
            DVector3 diff=pos0-origin;
            double dx0=diff.x(),dy0=diff.y();
            double wdir_dot_diff=diff.Dot(wdir);
            double tdir_dot_diff=diff.Dot(tdir);
            double tdir_dot_wdir=tdir.Dot(wdir);
            double tdir2=tdir.Mag2();
            double wdir2=wdir.Mag2();
            double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
            double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
            double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
            double scale=1./D;
            double s=scale*N;
            double t=scale*N1;
            diff+=s*tdir-t*wdir;
            double d=diff.Mag()+d_EPS; // prevent division by zero

            // The next measurement and its variance
            double tdrift=cdc_hits[cdc_index]->tdrift-trajectory[k].t;
            V_CDC=CDCDriftVariance(tdrift);

            double phi_d=diff.Phi();
            double dphi=phi_d-origin.Phi();
            while (dphi>M_PI) dphi-=2*M_PI;
            while (dphi<-M_PI) dphi+=2*M_PI;

            int ring_index=cdc_hits[cdc_index]->wire->ring-1;
            int straw_index=cdc_hits[cdc_index]->wire->straw-1;
            double dz=t*wdir.z();
            double delta=max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
               *cos(phi_d+sag_phi_offset[ring_index][straw_index]);
            double dmeas=CDCDriftDistance(dphi,delta,tdrift);

            // residual
            double res=dmeas-d;
            if (VERBOSE>5) jout << " Residual " << res << endl;
            /*
               DVector3 WirePosNew = wirepos+s*wdir;
               double dy = S(state_y)+S(state_ty)*s-WirePosNew.Y();
               double dx = S(state_x)+S(state_tx)*s-WirePosNew.X();
               double cosstereo=cos(cdc_hits[cdc_index]->wire->stereo);
               double thisd = sqrt(dx*dx+dy*dy)*cosstereo+d_EPS;
               double cosstereo2_over_d=cosstereo*cosstereo/thisd;
               H_CDC(state_x)=H_T_CDC(state_x)=dx*cosstereo2_over_d;
               H_CDC(state_y)=H_T_CDC(state_y)=dy*cosstereo2_over_d;
               H_CDC(state_tx)=H_T_CDC(state_tx)=s*H_CDC(state_x);
               H_CDC(state_ty)=H_T_CDC(state_ty)=s*H_CDC(state_y);
               */
            // Track projection
            double one_over_d=1./d;
            double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
            double wx=wdir.x(),wy=wdir.y();

            double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
            double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
            double dtdtx=scale*(dN1dtx-t*dDdtx);

            double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
            double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
            double dtdty=scale*(dN1dty-t*dDdty);

            double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
            double dsdtx=scale*(dNdtx-s*dDdtx);

            double dNdty=wy*wdir_dot_diff-wdir2*dy0;
            double dsdty=scale*(dNdty-s*dDdty);

            H_CDC(state_tx)=H_T_CDC(state_tx)
               =one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
                     +diffz*(dsdtx-dtdtx));
            H_CDC(state_ty)=H_T_CDC(state_ty)
               =one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
                     +diffz*(dsdty-dtdty));

            double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
            double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
            double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
            double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);

            H_CDC(state_x)=H_T_CDC(state_x)
               =one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
                     +diffz*(dsdx-dtdx));
            H_CDC(state_y)=H_T_CDC(state_y)
               =one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
                     +diffz*(dsdy-dtdy));

            double InvV=1./(V_CDC+H_CDC*C*H_T_CDC);

            // Check how far this hit is from the projection
            double chi2check=res*res*InvV;
            if (chi2check < CHI2CUT || DO_PRUNING == 0){
               if (VERBOSE) jout << "CDC Hit Added to FDC track " << endl;
               // Compute Kalman gain matrix
               K_CDC=InvV*(C*H_T_CDC);
               // Update state vector covariance matrix
               DMatrix4x4 Ctest=C-K_CDC*(H_CDC*C);

               //C.Print();
               //K.Print();
               //Ctest.Print();

               // Check that Ctest is positive definite
               if (!Ctest.IsPosDef()) return VALUE_OUT_OF_RANGE;
               C=Ctest;
               if(VERBOSE>10) C.Print();
               // Update the state vector
               //S=S+res*K;
               S+=res*K_CDC;
               if(VERBOSE) {jout << "traj[z]=" << trajectory[k].z<< endl; S.Print();} 

               // Compute new residual
               //d=finder->FindDoca(trajectory[k].z,S,wdir,origin);
               res=res-H_CDC*K_CDC*res;

               // Update chi2
               double fit_V=V_CDC-H_CDC*C*H_T_CDC;
               chi2+=res*res/fit_V;
               ndof++;

               // fill updates
               cdc_updates[cdc_index].resi=res;
               cdc_updates[cdc_index].d=d;
               cdc_updates[cdc_index].delta=delta;
               cdc_updates[cdc_index].S=S;
               cdc_updates[cdc_index].C=C;
               cdc_updates[cdc_index].V=V_CDC;
               cdc_updates[cdc_index].tdrift=tdrift;
               cdc_updates[cdc_index].ddrift=dmeas;
               cdc_updates[cdc_index].s=29.98*trajectory[k].t; // assume beta=1
               trajectory[k].id=cdc_index+1000;

            }
            // move to next cdc hit
            if (cdc_index>0){
               cdc_index--;

               //New wire position
               wire=cdc_hits[cdc_index]->wire;
               if (VERBOSE>5) jout << " Next Wire ring " << wire->ring << endl;
               origin=wire->origin;
               double vz=wire->udir.z();
               wdir=(1./vz)*wire->udir;
               wirepos=origin+((trajectory[k].z-z0))*wdir;

               // New doca^2
               dx=S(state_x)-wirepos.x();
               dy=S(state_y)-wirepos.y();
               doca2=dx*dx+dy*dy;

            }
            else more_hits=false;
         }
         firstCDCStep=false;
         old_doca2=doca2;
      }
   }

   ndof-=4;

   return NOERROR;
}


// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.

jerror_t 
DTrackCandidate_factory_StraightLine::Smooth(deque<trajectory_t>&trajectory,
      vector<fdc_update_t>&fdc_updates,
      vector<const DFDCPseudo *>&hits,
      vector<cdc_update_t>&cdc_updates,
      vector<const DCDCTrackHit *>&cdc_hits,
      DTrackCandidate *cand){ 
   unsigned int max=trajectory.size()-1;
   DMatrix4x1 S=(trajectory[max].Skk);
   DMatrix4x4 C=(trajectory[max].Ckk);
   DMatrix4x4 JT=trajectory[max].J.Transpose();
   DMatrix4x1 Ss=S;
   DMatrix4x4 Cs=C;
   DMatrix4x4 A,dC;

   const double d_EPS=1e-8;

   for (unsigned int m=max-1;m>0;m--){
      if (trajectory[m].id>0 && trajectory[m].id<1000){ // FDC Hit
         unsigned int id=trajectory[m].id-1;
         A=fdc_updates[id].C*JT*C.Invert();
         Ss=fdc_updates[id].S+A*(Ss-S);

         dC=A*(Cs-C)*A.Transpose();
         Cs=fdc_updates[id].C+dC;

         double cosa=cos(hits[id]->wire->angle);
         double cos2a=cos(2*hits[id]->wire->angle);
         double sina=sin(hits[id]->wire->angle);
         double u=hits[id]->w;
         double v=hits[id]->s;

         // Position and direction from state vector
         double x=Ss(state_x);
         double y=Ss(state_y);
         double tx=Ss(state_tx);
         double ty=Ss(state_ty);

         // Small angle alignment correction
         x = x + hits[id]->wire->angles.Z()*y;
         y = y - hits[id]->wire->angles.Z()*x;
         //tz = 1. + my_fdchits[id]->phiY*tx - my_fdchits[id]->phiX*ty;
         tx = (tx + hits[id]->wire->angles.Z()*ty - hits[id]->wire->angles.Y());
         ty = (ty - hits[id]->wire->angles.Z()*tx + hits[id]->wire->angles.X());

         // Projected position along the wire 
         double vpred=x*sina+y*cosa;

         // Projected position in the plane of the wires transverse to the wires
         double upred=x*cosa-y*sina;

         // Direction tangent in the u-z plane
         double tu=tx*cosa-ty*sina;
         double alpha=atan(tu);
         double cosalpha=cos(alpha);
         //double cosalpha2=cosalpha*cosalpha;
         double sinalpha=sin(alpha);

         // (signed) distance of closest approach to wire
         double du=upred-u;
         double doca=du*cosalpha;
         // Difference between measurement and projection for the cathodes
         double tv=tx*sina+ty*cosa;
         double resi_c=v-vpred;

         // Difference between measurement and projection perpendicular to the wire
         double drift=0.; // assume hit at wire position
         if (USE_FDC_DRIFT_TIMES){
            double drift_time=fdc_updates[id].tdrift;
            drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time);
         }
         double resi_a=drift-doca;

         // Variance from filter step
         DMatrix2x2 V=fdc_updates[id].V;
         // Compute projection matrix and find the variance for the residual
         DMatrix4x2 H_T;
         double temp2=-tv*sinalpha;
         H_T(state_x,1)=sina+cosa*cosalpha*temp2;	
         H_T(state_y,1)=cosa-sina*cosalpha*temp2;	

         double cos2_minus_sin2=cosalpha*cosalpha-sinalpha*sinalpha;
         double doca_cosalpha=doca*cosalpha;
         H_T(state_tx,1)=-doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2);
         H_T(state_ty,1)=-doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2);

         H_T(state_x,0)=cosa*cosalpha;
         H_T(state_y,0)=-sina*cosalpha;
         double one_plus_tu2=1.+tu*tu;
         double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
         H_T(state_ty,0)=sina*factor;
         H_T(state_tx,0)=-cosa*factor;

         // Matrix transpose H_T -> H
         DMatrix2x4 H;
         H(0,state_x)=H_T(state_x,0);
         H(0,state_y)=H_T(state_y,0);
         H(0,state_tx)=H_T(state_tx,0);
         H(0,state_ty)=H_T(state_ty,0);
         H(1,state_x)=H_T(state_x,1);
         H(1,state_y)=H_T(state_y,1);
         H(1,state_tx)=H_T(state_tx,1);
         H(1,state_ty)=H_T(state_ty,1);

         if (hits[id]->wire->layer==PLANE_TO_SKIP){
            //V+=Cs.SandwichMultiply(H_T);
            V=V+H*Cs*H_T;
         }
         else{
            //V-=dC.SandwichMultiply(H_T);
            V=V-H*dC*H_T;
         }

         if(DEBUG_HISTS){
            hFDCOccTrkSmooth->Fill(hits[id]->wire->layer);
         }

         // Implement derivatives wrt track parameters needed for millepede alignment

         DTrackFitter::pull_t thisPull(resi_a,sqrt(V(0,0)),
               trajectory[m].t*SPEED_OF_LIGHT,
               fdc_updates[id].tdrift,
               fdc_updates[id].d,
               NULL,hits[id],
               0.0, //docaphi
               trajectory[m].z, 
               0.0, //tcorr
               resi_c, sqrt(V(1,1))
               );

         if (hits[id]->wire->layer!=PLANE_TO_SKIP){
            vector<double> derivatives;
            derivatives.resize(FDCTrackD::size);

            //dDOCAW/dDeltaX
            derivatives[FDCTrackD::dDOCAW_dDeltaX] = -(1/sqrt(1 + pow(tx*cosa - ty*sina,2)));

            //dDOCAW/dDeltaZ
            derivatives[FDCTrackD::dDOCAW_dDeltaZ] = (tx*cosa - ty*sina)/sqrt(1 + pow(tx*cosa - ty*sina,2));

            //dDOCAW/ddeltaPhiX
            derivatives[FDCTrackD::dDOCAW_dDeltaPhiX] = (sina*(-(tx*cosa) + ty*sina)*(u - x*cosa + y*sina))/pow(1 + pow(tx*cosa - ty*sina,2),1.5);

            //dDOCAW/ddeltaphiY
            derivatives[FDCTrackD::dDOCAW_dDeltaPhiY] = (cosa*(tx*cosa - ty*sina)*(-u + x*cosa - y*sina))/pow(1 + pow(tx*cosa - ty*sina,2),1.5);

            //dDOCAW/ddeltaphiZ
            derivatives[FDCTrackD::dDOCAW_dDeltaPhiZ] = (tx*ty*u*cos2a + (x + pow(ty,2)*x - tx*ty*y)*sina + 
                  cosa*(-(tx*ty*x) + y + pow(tx,2)*y + (pow(tx,2) - pow(ty,2))*u*sina))/
               pow(1 + pow(tx*cosa - ty*sina,2),1.5);

            // dDOCAW/dx
            derivatives[FDCTrackD::dDOCAW_dx] = cosa/sqrt(1 + pow(tx*cosa - ty*sina,2));

            // dDOCAW/dy
            derivatives[FDCTrackD::dDOCAW_dy] = -(sina/sqrt(1 + pow(tx*cosa - ty*sina,2)));

            // dDOCAW/dtx
            derivatives[FDCTrackD::dDOCAW_dtx] = -((cosa*(tx*cosa - ty*sina)*(-u + x*cosa - y*sina))/pow(1 + pow(tx*cosa - ty*sina,2),1.5));

            // dDOCAW/dty
            derivatives[FDCTrackD::dDOCAW_dty] = (sina*(-(tx*cosa) + ty*sina)*(u - x*cosa + y*sina))/pow(1 + pow(tx*cosa - ty*sina,2),1.5); 

            // And the cathodes
            //dDOCAW/ddeltax
            derivatives[FDCTrackD::dDOCAC_dDeltaX] = 0.;

            //dDOCAW/ddeltax
            derivatives[FDCTrackD::dDOCAC_dDeltaZ] = ty*cosa + tx*sina;

            //dDOCAW/ddeltaPhiX
            derivatives[FDCTrackD::dDOCAC_dDeltaPhiX] = 0.;

            //dDOCAW/ddeltaPhiX
            derivatives[FDCTrackD::dDOCAC_dDeltaPhiY] = 0.;

            //dDOCAW/ddeltaPhiX
            derivatives[FDCTrackD::dDOCAC_dDeltaPhiZ] = -(x*cosa) + y*sina;

            // dDOCAW/dx
            derivatives[FDCTrackD::dDOCAC_dx] = sina;

            // dDOCAW/dy
            derivatives[FDCTrackD::dDOCAW_dy] = cosa;

            // dDOCAW/dtx
            derivatives[FDCTrackD::dDOCAW_dtx] = 0.;

            // dDOCAW/dty
            derivatives[FDCTrackD::dDOCAW_dty] = 0.;

            thisPull.AddTrackDerivatives(derivatives);
         }

         cand->pulls.push_back(thisPull);

      }
      else if (trajectory[m].id>=1000){ // CDC Hit
         unsigned int id=trajectory[m].id-1000;
         A=cdc_updates[id].C*JT*C.Invert();
         Ss=cdc_updates[id].S+A*(Ss-S);

         dC=A*(Cs-C)*A.Transpose();
         Cs=cdc_updates[id].C+dC;
         if (VERBOSE > 10) {
            jout << " In Smoothing Step Using ID " << id << "/" << cdc_updates.size() << " for ring " << cdc_hits[id]->wire->ring << endl;
            jout << " A cdc_updates[id].C Ss Cs " << endl;
            A.Print(); cdc_updates[id].C.Print(); Ss.Print(); Cs.Print();
         }
         if(!Cs.IsPosDef()) {
            if (VERBOSE) jout << "Cs is not PosDef!" << endl;
            return VALUE_OUT_OF_RANGE;
         }

         const DCDCWire *wire=cdc_hits[id]->wire;
         DVector3 origin=wire->origin;
         double z0=origin.z();
         double vz=wire->udir.z();
         DVector3 wdir=(1./vz)*wire->udir;
         DVector3 wirepos=origin+(trajectory[m].z-z0)*wdir;
         // Position and direction from state vector
         double x=Ss(state_x);
         double y=Ss(state_y);
         double tx=Ss(state_tx);
         double ty=Ss(state_ty);

         DVector3 pos0(x,y,trajectory[m].z);
         DVector3 tdir(tx,ty,1.);

         // Find the true doca to the wire
         DVector3 diff=pos0-origin;
         double dx0=diff.x(),dy0=diff.y();
         double wdir_dot_diff=diff.Dot(wdir);
         double tdir_dot_diff=diff.Dot(tdir);
         double tdir_dot_wdir=tdir.Dot(wdir);
         double tdir2=tdir.Mag2();
         double wdir2=wdir.Mag2();
         double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
         double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
         double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
         double scale=1./D;
         double s=scale*N;
         double t=scale*N1;
         diff+=s*tdir-t*wdir;
         double d=diff.Mag()+d_EPS; // prevent division by zero
         double ddrift = cdc_updates[id].ddrift;

         double resi = ddrift - d;


         // Track projection
         DMatrix1x4 H; DMatrix4x1 H_T;
         {
            double one_over_d=1./d;
            double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
            double wx=wdir.x(),wy=wdir.y();

            double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
            double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
            double dtdtx=scale*(dN1dtx-t*dDdtx);

            double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
            double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
            double dtdty=scale*(dN1dty-t*dDdty);

            double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
            double dsdtx=scale*(dNdtx-s*dDdtx);

            double dNdty=wy*wdir_dot_diff-wdir2*dy0;
            double dsdty=scale*(dNdty-s*dDdty);

            H(state_tx)=H_T(state_tx)
               =one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
                     +diffz*(dsdtx-dtdtx));
            H(state_ty)=H_T(state_ty)
               =one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
                     +diffz*(dsdty-dtdty));

            double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
            double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
            double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
            double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);

            H(state_x)=H_T(state_x)
               =one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
                     +diffz*(dsdx-dtdx));
            H(state_y)=H_T(state_y)
               =one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
                     +diffz*(dsdy-dtdy));
         }
         /*
            DVector3 WirePosNew = wirepos+s*wdir;
            double dy = Ss(state_y)+Ss(state_ty)*s-WirePosNew.Y();
            double dx = Ss(state_x)+Ss(state_tx)*s-WirePosNew.X();
            double cosstereo=cos(cdc_hits[id]->wire->stereo);
            double thisd = sqrt(dx*dx+dy*dy)*cosstereo+d_EPS;
            double cosstereo2_over_d=cosstereo*cosstereo/thisd;
            DMatrix1x4 H; DMatrix4x1 H_T;
            H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
            H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;
            H(state_tx)=H_T(state_tx)=s*H(state_x);
            H(state_ty)=H_T(state_ty)=s*H(state_y);
            */
         double V=cdc_updates[id].V;

         if (VERBOSE > 10) jout << " d " << d << " H*S " << H*S << endl;
         V=V-H*Cs*H_T;
         if (V<0) return VALUE_OUT_OF_RANGE;

         // Add the pull
         DTrackFitter::pull_t thisPull(resi,sqrt(V),
               trajectory[m].t*SPEED_OF_LIGHT,
               cdc_updates[id].tdrift,
               d,
               cdc_hits[id], NULL,
               diff.Phi(), //docaphi
               trajectory[m].z,
               cdc_updates[id].tdrift);

         // Derivatives for alignment
         double wtx=wire->udir.X(), wty=wire->udir.Y(), wtz=wire->udir.Z();
         double wx=wire->origin.X(), wy=wire->origin.Y(), wz=wire->origin.Z();

         double z=trajectory[m].z;
         double tx2=tx*tx, ty2=ty*ty;
         double wtx2=wtx*wtx, wty2=wty*wty, wtz2=wtz*wtz;
         double denom=(1 + ty2)*wtx2 + (1 + tx2)*wty2 - 2*ty*wty*wtz + (tx2 + ty2)*wtz2 - 2*tx*wtx*(ty*wty + wtz)+d_EPS;
         double denom2=denom*denom;
         double c1=-(wtx - tx*wtz)*(wy - y);
         double c2=wty*(wx - tx*wz - x + tx*z);
         double c3=ty*(-(wtz*wx) + wtx*wz + wtz*x - wtx*z);
         double dscale=0.5*(1/d);

         vector<double> derivatives(11);

         derivatives[CDCTrackD::dDOCAdOriginX]=dscale*(2*(wty - ty*wtz)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdOriginY]=dscale*(2*(-wtx + tx*wtz)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdOriginZ]=dscale*(2*(ty*wtx - tx*wty)*(c1 + c2 + c3))/denom;

         derivatives[CDCTrackD::dDOCAdDirX]=dscale*(2*(wty - ty*wtz)*(c1 + c2 + c3)*
               (tx*(ty*wty + wtz)*(wx - x) + (wty - ty*wtz)*(-wy + y + ty*(wz - z)) +
                wtx*(-((1 + ty2)*wx) + (1 + ty2)*x + tx*(ty*wy + wz - ty*y - z)) + tx2*(wty*(-wy + y) + wtz*(-wz + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdDirY]=dscale*(-2*(wtx - tx*wtz)*(c1 + c2 + c3)*
               (tx*(ty*wty + wtz)*(wx - x) + (wty - ty*wtz)*(-wy + y + ty*(wz - z)) +
                wtx*(-((1 + ty2)*wx) + (1 + ty2)*x + tx*(ty*wy + wz - ty*y - z)) + tx2*(wty*(-wy + y) + wtz*(-wz + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdDirZ]=dscale*(-2*(ty*wtx - tx*wty)*(c1 + c2 + c3)*
               (-(tx*(ty*wty + wtz)*(wx - x)) + tx2*(wty*(wy - y) + wtz*(wz - z)) + (wty - ty*wtz)*(wy - y + ty*(-wz + z)) +
                wtx*((1 + ty2)*wx - (1 + ty2)*x + tx*(-(ty*wy) - wz + ty*y + z))))/denom2;

         derivatives[CDCTrackD::dDOCAdS0]=-derivatives[CDCTrackD::dDOCAdOriginX];

         derivatives[CDCTrackD::dDOCAdS1]=-derivatives[CDCTrackD::dDOCAdOriginY];

         derivatives[CDCTrackD::dDOCAdS2]=dscale*(2*(wty - ty*wtz)*(-c1 - c2 - c3)*
               (-(wtx*wtz*wx) - wty*wtz*wy + wtx2*wz + wty2*wz + wtx*wtz*x + wty*wtz*y - wtx2*z - wty2*z +
                tx*(wty2*(wx - x) + wtx*wty*(-wy + y) + wtz*(wtz*wx - wtx*wz - wtz*x + wtx*z)) +
                ty*(wtx*wty*(-wx + x) + wtx2*(wy - y) + wtz*(wtz*wy - wty*wz - wtz*y + wty*z))))/denom2;

         derivatives[CDCTrackD::dDOCAdS3]=dscale*(2*(wtx - tx*wtz)*(c1 + c2 + c3)*
               (-(wtx*wtz*wx) - wty*wtz*wy + wtx2*wz + wty2*wz + wtx*wtz*x + wty*wtz*y - wtx2*z - wty2*z +
                tx*(wty2*(wx - x) + wtx*wty*(-wy + y) + wtz*(wtz*wx - wtx*wz - wtz*x + wtx*z)) +
                ty*(wtx*wty*(-wx + x) + wtx2*(wy - y) + wtz*(wtz*wy - wty*wz - wtz*y + wty*z))))/denom2;

         thisPull.AddTrackDerivatives(derivatives);

         cand->pulls.push_back(thisPull);

      }
      else{
         A=trajectory[m].Ckk*JT*C.Invert();
         Ss=trajectory[m].Skk+A*(Ss-S);
         Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
      }

      S=trajectory[m].Skk;
      C=trajectory[m].Ckk;
      JT=trajectory[m].J.Transpose();
   }

   return NOERROR;
}

