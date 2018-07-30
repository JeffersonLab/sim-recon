// $Id: DTrackWireBased_factory.cc 5612 2009-10-15 20:51:25Z staylor $
//
//    File: DTrackWireBased_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
#include <cmath>
#include <mutex>
using namespace std;

#include "DTrackWireBased_factory.h"
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <SplitString.h>

#include <TROOT.h>

#define C_EFFECTIVE     15.

using namespace jana;

//------------------
// CDCSortByRincreasing
//------------------
bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) {
   // use the ring number to sort by R(decreasing) and then straw(increasing)
   if(hit1->wire->ring == hit2->wire->ring){
      return hit1->wire->straw < hit2->wire->straw;
   }
   return hit1->wire->ring < hit2->wire->ring;
}

//------------------
// FDCSortByZincreasing
//------------------
bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2) {
   // use the layer number to sort by Z(decreasing) and then wire(increasing)
   if(hit1->wire->layer == hit2->wire->layer){
      return hit1->wire->wire < hit2->wire->wire;
   }
   return hit1->wire->layer < hit2->wire->layer;
}


//------------------
// count_common_members
//------------------
   template<typename T>
static unsigned int count_common_members(vector<T> &a, vector<T> &b)
{
   unsigned int n=0;
   for(unsigned int i=0; i<a.size(); i++){
      for(unsigned int j=0; j<b.size(); j++){
         if(a[i]==b[j])n++;
      }
   }

   return n;
}

bool DTrackWireBased_cmp(DTrackWireBased *a,DTrackWireBased *b){
  if (a->candidateid==b->candidateid) return a->mass()<b->mass();
  return a->candidateid<b->candidateid;
}

//------------------
// init
//------------------
jerror_t DTrackWireBased_factory::init(void)
{
   fitter = NULL;
   MAX_DReferenceTrajectoryPoolSize = 50;

   //DEBUG_HISTS = true;	
   DEBUG_HISTS = false;
   DEBUG_LEVEL = 0;

   gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",DEBUG_LEVEL);

   SKIP_MASS_HYPOTHESES_WIRE_BASED=false; 
   gPARMS->SetDefaultParameter("TRKFIT:SKIP_MASS_HYPOTHESES_WIRE_BASED",
         SKIP_MASS_HYPOTHESES_WIRE_BASED);

   mNumHypPlus=1;
   mNumHypMinus=1;
   if(!SKIP_MASS_HYPOTHESES_WIRE_BASED)
   {
      vector<int> hypotheses;
      hypotheses.push_back(Positron);
      hypotheses.push_back(PiPlus);
      hypotheses.push_back(KPlus);
      hypotheses.push_back(Proton);
      hypotheses.push_back(Electron);
      hypotheses.push_back(PiMinus);
      hypotheses.push_back(KMinus);
      hypotheses.push_back(AntiProton);

      ostringstream locMassStream;
      for(size_t loc_i = 0; loc_i < hypotheses.size(); ++loc_i)
      {
         locMassStream << hypotheses[loc_i];
         if(loc_i != (hypotheses.size() - 1))
            locMassStream << ",";
      }

      string HYPOTHESES = locMassStream.str();
      gPARMS->SetDefaultParameter("TRKFIT:HYPOTHESES", HYPOTHESES);

      // Parse MASS_HYPOTHESES strings to make list of masses to try
      hypotheses.clear();
      SplitString(HYPOTHESES, hypotheses, ",");
      for(size_t loc_i = 0; loc_i < hypotheses.size(); ++loc_i)
      {
         if(ParticleCharge(Particle_t(hypotheses[loc_i])) > 0)
            mass_hypotheses_positive.push_back(hypotheses[loc_i]);
         else if(ParticleCharge(Particle_t(hypotheses[loc_i])) < 0)
            mass_hypotheses_negative.push_back(hypotheses[loc_i]);
      }

      if(mass_hypotheses_positive.empty()){
         static once_flag pwarn_flag;
         call_once(pwarn_flag, [](){
               jout << endl;
               jout << "############# WARNING !! ################ " <<endl;
               jout << "There are no mass hypotheses for positive tracks!" << endl;
               jout << "Be SURE this is what you really want!" << endl;
               jout << "######################################### " <<endl;
               jout << endl;
               });
      }
      if(mass_hypotheses_negative.empty()){
         static once_flag nwarn_flag;
         call_once(nwarn_flag, [](){
               jout << endl;
               jout << "############# WARNING !! ################ " <<endl;
               jout << "There are no mass hypotheses for negative tracks!" << endl;
               jout << "Be SURE this is what you really want!" << endl;
               jout << "######################################### " <<endl;
               jout << endl;
               });
      }
      mNumHypPlus=mass_hypotheses_positive.size();
      mNumHypMinus=mass_hypotheses_negative.size();
   }

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackWireBased_factory::brun(jana::JEventLoop *loop, int32_t runnumber)
{
   // Get the geometry
   DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
   geom = dapp->GetDGeometry(runnumber);
   // Check for magnetic field
   dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);

   // Get pointer to DTrackFitter object that actually fits a track
   vector<const DTrackFitter *> fitters;
   loop->Get(fitters);
   if(fitters.size()<1){
      _DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
      return RESOURCE_UNAVAILABLE;
   }

   // Drop the const qualifier from the DTrackFitter pointer (I'm surely going to hell for this!)
   fitter = const_cast<DTrackFitter*>(fitters[0]);

   // Warn user if something happened that caused us NOT to get a fitter object pointer
   if(!fitter){
      _DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
      return RESOURCE_UNAVAILABLE;
   }

   USE_HITS_FROM_CANDIDATE=false;
   gPARMS->SetDefaultParameter("TRKFIT:USE_HITS_FROM_CANDIDATE",
         USE_HITS_FROM_CANDIDATE);

   MIN_FIT_P = 0.050; // GeV
   gPARMS->SetDefaultParameter("TRKFIT:MIN_FIT_P", MIN_FIT_P, "Minimum fit momentum in GeV/c for fit to be considered successful");

   if(DEBUG_HISTS){
      dapp->Lock();

      // Histograms may already exist. (Another thread may have created them)
      // Try and get pointers to the existing ones.


      dapp->Unlock();
   }

   // Get the particle ID algorithms
   loop->GetSingle(dPIDAlgorithm);

   //Pre-allocate memory for DReferenceTrajectory objects early
   //The swim-step objects of these arrays take up a significant amount of memory, and it can be difficult to find enough free contiguous space for them.
   //Therefore, allocate them at the beginning before the available memory becomes randomly populated
   while(rtv.size() < MAX_DReferenceTrajectoryPoolSize)
      rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));

   // Outer detector geometry parameters
   geom->GetFCALZ(dFCALz); 
   vector<double>tof_face;
   geom->Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z",
	      tof_face);
   vector<double>tof_plane;  
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", tof_plane);
   dTOFz=tof_face[2]+tof_plane[2]; 
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", tof_plane);
   dTOFz+=tof_face[2]+tof_plane[2];
   dTOFz*=0.5;  // mid plane between tof planes
   
    // Get start counter geometry;
   if (geom->GetStartCounterGeom(sc_pos,sc_norm)){
     // Create vector of direction vectors in scintillator planes
     for (int i=0;i<30;i++){
       vector<DVector3>temp;
       for (unsigned int j=0;j<sc_pos[i].size()-1;j++){
	 double dx=sc_pos[i][j+1].x()-sc_pos[i][j].x();
	 double dy=sc_pos[i][j+1].y()-sc_pos[i][j].y();
	 double dz=sc_pos[i][j+1].z()-sc_pos[i][j].z();
	 temp.push_back(DVector3(dx/dz,dy/dz,1.));
       }
     sc_dir.push_back(temp);
     }
     SC_END_NOSE_Z=sc_pos[0][12].z();
     SC_BARREL_R=sc_pos[0][0].Perp();
     SC_PHI_SECTOR1=sc_pos[0][0].Phi();
   }

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackWireBased_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   if(!fitter)return NOERROR; 
   
   if(rtv.size() > MAX_DReferenceTrajectoryPoolSize){
     for(size_t loc_i = MAX_DReferenceTrajectoryPoolSize; loc_i < rtv.size(); ++loc_i)
       delete rtv[loc_i];
     rtv.resize(MAX_DReferenceTrajectoryPoolSize);
   }


   // Get candidates and hits
   vector<const DTrackCandidate*> candidates;
   loop->Get(candidates);

   if (candidates.size()==0) return NOERROR;

   // Reset the number of used reference trajectories from the pool
   num_used_rts=0;

   if (dIsNoFieldFlag){
      // Copy results over from the StraightLine or CDCCOSMIC candidate and
     // add extrapolations
      for (unsigned int i=0;i<candidates.size();i++){
         const DTrackCandidate *cand=candidates[i];

         // Make a new wire-based track
         DTrackWireBased *track = new DTrackWireBased(); //share the memory: isn't changed below
         *static_cast<DTrackingData*>(track) = *static_cast<const DTrackingData*>(cand);
         track->IsSmoothed = cand->IsSmoothed;

         // candidate id
         track->candidateid=i+1;

         // Track quality properties
         track->Ndof=cand->Ndof;
         track->chisq=cand->chisq;	
         track->FOM = TMath::Prob(track->chisq, track->Ndof);

         // pull vector
         track->pulls = cand->pulls;

         // Lists of hits used in the previous pass
         vector<const DCDCTrackHit *>cdchits;
         cand->GetT(cdchits);
         vector<const DFDCPseudo *>fdchits;
         cand->GetT(fdchits);

         for (unsigned int k=0;k<cdchits.size();k++){
            track->AddAssociatedObject(cdchits[k]);
         }
         for (unsigned int k=0;k<fdchits.size();k++){
            track->AddAssociatedObject(fdchits[k]);
         }
         track->dCDCRings = dPIDAlgorithm->Get_CDCRingBitPattern(cdchits);
         track->dFDCPlanes = dPIDAlgorithm->Get_FDCPlaneBitPattern(fdchits);

	 // Create the extrapolation vectors
	 vector<DTrackFitter::Extrapolation_t>myvector;
	 track->extrapolations.emplace(SYS_BCAL,myvector);
	 track->extrapolations.emplace(SYS_TOF,myvector);
	 track->extrapolations.emplace(SYS_FCAL,myvector);
	 track->extrapolations.emplace(SYS_FDC,myvector);
	 track->extrapolations.emplace(SYS_CDC,myvector);
	 track->extrapolations.emplace(SYS_START,myvector);	

	 // Extrapolate to TOF
	 DVector3 dir=track->momentum();
	 dir.SetMag(1.);
	 DVector3 pos0=track->position();
	 double z0=track->position().z();
	 double z=z0;
	 double uz=dir.z();
	 DVector3 diff=((dTOFz-z0)/uz)*dir;
	 DVector3 pos=pos0+diff;
	 double s=diff.Mag();
	 double t=s/29.98;
	 track->extrapolations[SYS_TOF].push_back(DTrackFitter::Extrapolation_t(pos,dir,t,s));	 
	 // Extrapolate to FCAL
	 diff=((dFCALz-z0)/uz)*dir;
	 pos=pos0+diff;
	 s=diff.Mag();
	 t=s/29.98;
	 track->extrapolations[SYS_FCAL].push_back(DTrackFitter::Extrapolation_t(pos,dir,t,s));  
	 // extrapolate to exit of FCAL
	 diff=((dFCALz+45.-z0)/uz)*dir;
	 pos=pos0+diff;
	 s=diff.Mag();
	 t=s/29.98;
	 track->extrapolations[SYS_FCAL].push_back(DTrackFitter::Extrapolation_t(pos,dir,t,s));
	  
	 // Extrapolate to Start Counter and BCAL
	 double R=pos0.Perp();
	 diff.SetMag(0.);
	 while (R<89.0 && z>17. && z<410.){
	   diff+=(1./dir.z())*dir;
	   pos=pos0+diff;
	   R=pos.Perp();
	   z=pos.z();
	   s=diff.Mag();
	   t=s/29.98;
	   //	   printf("R %f z %f\n",R,z);
	   // start counter
	   if (sc_pos.empty()==false && R<SC_BARREL_R && z<SC_END_NOSE_Z){
	     double d_old=1000.,d=1000.;
	     unsigned int index=0;
	     for (unsigned int m=0;m<12;m++){
	       double dphi=pos.Phi()-SC_PHI_SECTOR1;
	       if (dphi<0) dphi+=2.*M_PI;
	       index=int(floor(dphi/(2.*M_PI/30.)));
	       if (index>29) index=0;
	       d=sc_norm[index][m].Dot(pos-sc_pos[index][m]);
	       if (d*d_old<0){ // break if we cross the current plane  
		 // Find the new distance to the start counter (which could 
		 // now be to a plane in the one adjacent to the one before the
		 // step...)
		 int count=0;
		 while (fabs(d)>0.05 && count<20){ 
		   // Find the index for the nearest start counter paddle
		   double dphi=pos.Phi()-SC_PHI_SECTOR1;
		   if (dphi<0) dphi+=2.*M_PI;
		   index=int(floor(dphi/(2.*M_PI/30.)));
		   d=sc_norm[index][m].Dot(pos-sc_pos[index][m]);
		   pos+=d*dir;
		   count++;
		 }
		 track->extrapolations[SYS_START].push_back(DTrackFitter::Extrapolation_t(pos,dir,t,s));
		 break;
	       }
	       d_old=d;
	     }
	   }
	   if (R>64.){	 
	     track->extrapolations[SYS_BCAL].push_back(DTrackFitter::Extrapolation_t(pos,dir,t,s));
	   }
	 }

         _data.push_back(track);
      }
      return NOERROR;
   }




   // Loop over candidates
   for(unsigned int i=0; i<candidates.size(); i++){
      const DTrackCandidate *candidate = candidates[i];

      // Skip candidates with momentum below some cutoff
      if (candidate->momentum().Mag()<MIN_FIT_P){
         continue;
      }

      if (SKIP_MASS_HYPOTHESES_WIRE_BASED){
         // Make sure there are enough DReferenceTrajectory objects
         unsigned int locNumInitialReferenceTrajectories = rtv.size();
         while(rtv.size()<=num_used_rts){
            //printf("Adding %d %d\n",rtv.size(),_data.size());
            rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
         }
         DReferenceTrajectory *rt = rtv[num_used_rts];
         if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
            rt->Reset();
         rt->SetDGeometry(geom);
         rt->q = candidate->charge();

         // Increment the number of used reference trajectories
         num_used_rts++;

         DoFit(i,candidate,rt,loop,0.13957);
      }
      else{
         // Choose list of mass hypotheses based on charge of candidate
         vector<int> mass_hypotheses;
         if(candidate->charge()<0.0){
            mass_hypotheses = mass_hypotheses_negative;
         }else{
            mass_hypotheses = mass_hypotheses_positive;
         }

         if ((!isfinite(candidate->momentum().Mag())) || (!isfinite(candidate->position().Mag())))
            _DBG_ << "Invalid seed data for event "<< eventnumber <<"..."<<endl;

         // Loop over potential particle masses
         for(unsigned int j=0; j<mass_hypotheses.size(); j++){
            if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting wire based fit with id: "<<mass_hypotheses[j]<<endl;}
            // Make sure there are enough DReferenceTrajectory objects
            unsigned int locNumInitialReferenceTrajectories = rtv.size();
            while(rtv.size()<=num_used_rts){
               //printf("Adding %d\n",rtv.size());
               rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
            }
            DReferenceTrajectory *rt = rtv[num_used_rts];
            if(locNumInitialReferenceTrajectories == rtv.size()){ //didn't create a new one
               rt->Reset();
	    }
	 
            rt->SetDGeometry(geom);
            rt->q = candidate->charge();

            // Increment the number of used reference trajectories
            num_used_rts++;

            DoFit(i,candidate,rt,loop,ParticleMass(Particle_t(mass_hypotheses[j])));
         }

      }
   }

   // Filter out duplicate tracks
   FilterDuplicates();

   // Add any missing hypotheses
   InsertMissingHypotheses();

   return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DTrackWireBased_factory::erun(void)
{
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackWireBased_factory::fini(void)
{
   for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
   rtv.clear();

   return NOERROR;
}

//------------------
// FilterDuplicates
//------------------
void DTrackWireBased_factory::FilterDuplicates(void)
{
   /// Look through all current DTrackWireBased objects and remove any
   /// that have all of their hits in common with another track

   if(_data.size()==0)return;

   if(DEBUG_LEVEL>2)_DBG_<<"Looking for clones of wire-based tracks ..."<<endl;

   set<unsigned int> indexes_to_delete;
   for(unsigned int i=0; i<_data.size()-1; i++){
      DTrackWireBased *dtrack1 = _data[i];

      vector<const DCDCTrackHit*> cdchits1;
      vector<const DFDCPseudo*> fdchits1;
      dtrack1->Get(cdchits1);
      dtrack1->Get(fdchits1);

      JObject::oid_t cand1=dtrack1->candidateid;
      for(unsigned int j=i+1; j<_data.size(); j++){
         DTrackWireBased *dtrack2 = _data[j];
         if (dtrack2->candidateid==cand1) continue;
         if (dtrack2->mass() != dtrack1->mass())continue;

         vector<const DCDCTrackHit*> cdchits2;
         vector<const DFDCPseudo*> fdchits2;
         dtrack2->Get(cdchits2);
         dtrack2->Get(fdchits2);

         // Count number of cdc and fdc hits in common
         unsigned int Ncdc = count_common_members(cdchits1, cdchits2);
         unsigned int Nfdc = count_common_members(fdchits1, fdchits2);
         unsigned int total = Ncdc + Nfdc;

         if (total==0) continue;
         if(Ncdc!=cdchits1.size() && Ncdc!=cdchits2.size())continue;
         if(Nfdc!=fdchits1.size() && Nfdc!=fdchits2.size())continue;

         unsigned int total1 = cdchits1.size()+fdchits1.size();
         unsigned int total2 = cdchits2.size()+fdchits2.size();
         if(total!=total1 && total!=total2)continue;

         if(total1<total2){
            // The two track candidates actually correspond to 
            // a single track.  Set the candidate id for this 
            // track to the candidate id from the clone match to 
            // prevent multiple clone tracks confusing matters 
            // at a later stage of the reconstruction...
            _data[j]->candidateid=cand1;
            indexes_to_delete.insert(i);
         }else{
            indexes_to_delete.insert(j);
         }

      }
   }

   if(DEBUG_LEVEL>2)_DBG_<<"Found "<<indexes_to_delete.size()<<" wire-based clones"<<endl;

   // Return now if we're keeping everyone
   if(indexes_to_delete.size()==0)return;

   // Copy pointers that we want to keep to a new container and delete
   // the clone objects
   vector<DTrackWireBased*> new_data;
   for(unsigned int i=0; i<_data.size(); i++){
      if(indexes_to_delete.find(i)==indexes_to_delete.end()){
         new_data.push_back(_data[i]);
      }else{
         delete _data[i];
         if(DEBUG_LEVEL>1)_DBG_<<"Deleting clone wire-based track "<<i<<endl;
      }
   }	
   _data = new_data;
}

// Routine to find the hits, do the fit, and fill the list of wire-based tracks
void DTrackWireBased_factory::DoFit(unsigned int c_id,
      const DTrackCandidate *candidate,
      DReferenceTrajectory *rt,
      JEventLoop *loop, double mass){ 
   // Do the fit
   DTrackFitter::fit_status_t status = DTrackFitter::kFitNotDone;
   if (USE_HITS_FROM_CANDIDATE) {
      fitter->Reset();
      fitter->SetFitType(DTrackFitter::kWireBased);	

      // Get the hits from the track candidate
      vector<const DFDCPseudo*>myfdchits;
      candidate->GetT(myfdchits);
      fitter->AddHits(myfdchits);
      vector<const DCDCTrackHit *>mycdchits;
      candidate->GetT(mycdchits);
      fitter->AddHits(mycdchits);

      status=fitter->FitTrack(candidate->position(),candidate->momentum(),
            candidate->charge(),mass,0.);
   }
   else{
     fitter->Reset();
      fitter->SetFitType(DTrackFitter::kWireBased);
      // Swim a reference trajectory using the candidate starting momentum
      // and position
      rt->SetMass(mass);
      //rt->Swim(candidate->position(),candidate->momentum(),candidate->charge());
      rt->FastSwimForHitSelection(candidate->position(),candidate->momentum(),candidate->charge());

      status=fitter->FindHitsAndFitTrack(*candidate,rt,loop,mass,candidate->Ndof+3);
      if (/*false && */status==DTrackFitter::kFitNotDone){
         if (DEBUG_LEVEL>1)_DBG_ << "Using hits from candidate..." << endl;
         fitter->Reset();

         // Get the hits from the candidate
         vector<const DFDCPseudo*>myfdchits;
         candidate->GetT(myfdchits);
         fitter->AddHits(myfdchits);
         vector<const DCDCTrackHit *>mycdchits;
         candidate->GetT(mycdchits);
         fitter->AddHits(mycdchits);

         status=fitter->FitTrack(candidate->position(),candidate->momentum(),
               candidate->charge(),mass,0.);
      }
   }

   // if the fit returns chisq=-1, something went terribly wrong... 
   if (fitter->GetChisq()<0){
     status=DTrackFitter::kFitFailed;
   }

   // Check the status of the fit
   switch(status){
      case DTrackFitter::kFitNotDone:
         //_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
      case DTrackFitter::kFitFailed:
         break;
      case DTrackFitter::kFitNoImprovement:	
      case DTrackFitter::kFitSuccess:
         if(!isfinite(fitter->GetFitParameters().position().X())) break;
         {    
            // Make a new wire-based track
             DTrackWireBased *track = new DTrackWireBased();
             *static_cast<DTrackingData*>(track) = fitter->GetFitParameters();

            track->chisq = fitter->GetChisq();
            track->Ndof = fitter->GetNdof();
            track->FOM = TMath::Prob(track->chisq, track->Ndof);
            track->pulls =std::move(fitter->GetPulls()); 
	    track->extrapolations=std::move(fitter->GetExtrapolations());
            track->candidateid = c_id+1;

            // Add hits used as associated objects
            vector<const DCDCTrackHit*> cdchits = fitter->GetCDCFitHits();
            vector<const DFDCPseudo*> fdchits = fitter->GetFDCFitHits();
            sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
            sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
            for(unsigned int m=0; m<cdchits.size(); m++)track->AddAssociatedObject(cdchits[m]);
            for(unsigned int m=0; m<fdchits.size(); m++)track->AddAssociatedObject(fdchits[m]);

	    // Set CDC ring & FDC plane hit patterns before candidate tracks are associated
	    vector<const DCDCTrackHit*> tempCDCTrackHits;
	    vector<const DFDCPseudo*> tempFDCPseudos;
	    track->Get(tempCDCTrackHits);
	    track->Get(tempFDCPseudos);
	    track->dCDCRings = dPIDAlgorithm->Get_CDCRingBitPattern(tempCDCTrackHits);
	    track->dFDCPlanes = dPIDAlgorithm->Get_FDCPlaneBitPattern(tempFDCPseudos);

            // Add DTrackCandidate as associated object
            track->AddAssociatedObject(candidate);

            _data.push_back(track);
            break;
         }
      default:
         break;
   }
}

// If the fit failed for certain hypotheses, fill in the gaps using data from
// successful fits for each candidate.
bool DTrackWireBased_factory::InsertMissingHypotheses(void){
  if (_data.size()==0) return false;
  
  // Make sure the tracks are ordered by candidate id
  sort(_data.begin(),_data.end(),DTrackWireBased_cmp);
  
  JObject::oid_t old_id=_data[0]->candidateid;
  unsigned int mass_bits=0;
  double q=_data[0]->charge();
  vector<DTrackWireBased*>myhypotheses;
  vector<DTrackWireBased*>tracks_to_add;
  for (size_t i=0;i<_data.size();i++){
    if (_data[i]->candidateid!=old_id){
      int num_hyp=myhypotheses.size();
      if ((q<0 && num_hyp!=mNumHypMinus)||(q>0 && num_hyp!=mNumHypPlus)){ 
	AddMissingTrackHypotheses(mass_bits,tracks_to_add,myhypotheses,q);
      }
      
      // Clear the myhypotheses vector for the next track
      myhypotheses.clear();
      // Reset charge 
      q=_data[i]->charge();	
   
      // Set the bit for this mass hypothesis
      mass_bits = 1<<_data[i]->PID();

      // Add the data to the myhypotheses vector
      myhypotheses.push_back(_data[i]);
    }
    else{
      myhypotheses.push_back(_data[i]);
      
      // Set the bit for this mass hypothesis
      mass_bits |= 1<< _data[i]->PID();
    }
    
    old_id=_data[i]->candidateid;
  }
  // Deal with last track candidate	
  int num_hyp=myhypotheses.size();
  if ((q<0 && num_hyp!=mNumHypMinus)||(q>0 && num_hyp!=mNumHypPlus)){
    AddMissingTrackHypotheses(mass_bits,tracks_to_add,myhypotheses,q);
  }
    
  // Add the new list of tracks to the output list
  if (tracks_to_add.size()>0){
    for (size_t i=0;i<tracks_to_add.size();i++){ 
      _data.push_back(tracks_to_add[i]);
    }
    // Make sure the tracks are ordered by candidate id
    sort(_data.begin(),_data.end(),DTrackWireBased_cmp);
  }

  return true;
}

// Create a track with a mass hypothesis that was not present in the list of 
// fitted tracks from an existing fitted track.
void DTrackWireBased_factory::AddMissingTrackHypothesis(vector<DTrackWireBased*>&tracks_to_add,
				      const DTrackWireBased *src_track,
							double my_mass,
							double q){
  // Create a new wire-based track object
  DTrackWireBased *wirebased_track = new DTrackWireBased();
  *static_cast<DTrackingData*>(wirebased_track) = *static_cast<const DTrackingData*>(src_track);

  // Copy over DKinematicData part from the result of a successful fit
  wirebased_track->setPID(IDTrack(q, my_mass));
  wirebased_track->chisq = src_track->chisq;
  wirebased_track->Ndof = src_track->Ndof;
  wirebased_track->pulls = src_track->pulls;
  wirebased_track->extrapolations = src_track->extrapolations;
  wirebased_track->candidateid=src_track->candidateid;
  wirebased_track->FOM=src_track->FOM;
  wirebased_track->IsSmoothed=src_track->IsSmoothed;
  wirebased_track->dCDCRings=src_track->dCDCRings;
  wirebased_track->dFDCPlanes=src_track->dFDCPlanes;

  // (Partially) compensate for the difference in energy loss between the 
  // source track and a particle of mass my_mass 
  DVector3 position,momentum;
  if (wirebased_track->extrapolations.at(SYS_CDC).size()>0){
    unsigned int index=wirebased_track->extrapolations.at(SYS_CDC).size()-1;
    position=wirebased_track->extrapolations[SYS_CDC][index].position;
    momentum=wirebased_track->extrapolations[SYS_CDC][index].momentum;
  }
  else if (wirebased_track->extrapolations.at(SYS_FDC).size()>0){
    unsigned int index=wirebased_track->extrapolations.at(SYS_FDC).size()-1;
    position=wirebased_track->extrapolations[SYS_FDC][index].position;
    momentum=wirebased_track->extrapolations[SYS_FDC][index].momentum;
  }
  if (momentum.Mag()>0.){
    CorrectForELoss(position,momentum,q,my_mass);
    
    wirebased_track->setMomentum(momentum);
    wirebased_track->setPosition(position);
  }

  // Get the hits used in the fit and add them as associated objects 
  vector<const DCDCTrackHit *>cdchits;
  src_track->GetT(cdchits);
  vector<const DFDCPseudo *>fdchits;
  src_track->GetT(fdchits);
  for(unsigned int m=0; m<fdchits.size(); m++)
    wirebased_track->AddAssociatedObject(fdchits[m]); 
  for(unsigned int m=0; m<cdchits.size(); m++)
    wirebased_track->AddAssociatedObject(cdchits[m]);
   
  tracks_to_add.push_back(wirebased_track);
}

// Use the FastSwim method in DReferenceTrajectory to propagate back to the 
// POCA to the beam line, adding a bit of energy at each step that would have 
// been lost had the particle emerged from the target.
void DTrackWireBased_factory::CorrectForELoss(DVector3 &position,DVector3 &momentum,double q,double my_mass){  
  // Make sure there are enough DReferenceTrajectory objects
  unsigned int locNumInitialReferenceTrajectories = rtv.size();
  while(rtv.size()<=num_used_rts){
    //printf("Adding %d\n",rtv.size());
    rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
  }
  DReferenceTrajectory *rt = rtv[num_used_rts];
  if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
    rt->Reset();
  rt->SetDGeometry(geom);
  rt->q = q;
  rt->SetMass(my_mass);
  rt->SetPLossDirection(DReferenceTrajectory::kBackward);
  DVector3 last_pos,last_mom;
  DVector3 origin(0.,0.,65.);
  DVector3 dir(0.,0.,1.);
  rt->FastSwim(position,momentum,last_pos,last_mom,rt->q,origin,dir,300.);   
  position=last_pos;
  momentum=last_mom;   
    
  // Increment the number of used reference trajectories
  num_used_rts++;
}


// Fill in all missing hypotheses for a given track candidate
void DTrackWireBased_factory::AddMissingTrackHypotheses(unsigned int mass_bits,
							vector<DTrackWireBased*>&tracks_to_add,
							vector<DTrackWireBased *>&myhypotheses,
							double q){ 
  Particle_t negative_particles[3]={KMinus,PiMinus,Electron};
  Particle_t positive_particles[3]={KPlus,PiPlus,Positron};

  unsigned int last_index=myhypotheses.size()-1;
  if (q>0){
    if ((mass_bits & (1<<Proton))==0){
      AddMissingTrackHypothesis(tracks_to_add,myhypotheses[last_index],
				ParticleMass(Proton),+1.);  
    } 
    for (int i=0;i<3;i++){
      if ((mass_bits & (1<<positive_particles[i]))==0){
	AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				  ParticleMass(positive_particles[i]),+1.); 
      } 
    }    
  }
  else{
    if ((mass_bits & (1<<AntiProton))==0){
      AddMissingTrackHypothesis(tracks_to_add,myhypotheses[last_index],
				ParticleMass(Proton),-1.);  
    } 	
    for (int i=0;i<3;i++){
      if ((mass_bits & (1<<negative_particles[i]))==0){
	AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				  ParticleMass(negative_particles[i]),-1.);  
      } 
    }
  }
} 
