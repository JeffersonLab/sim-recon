// $Id$
//
//    File: DTrackHitSelectorALT1.cc
// Created: Fri Feb  6 08:22:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include <TRACKING/DReferenceTrajectory.h>

#include "DTrackHitSelectorALT1.h"

#ifndef ansi_escape
#define ansi_escape			((char)0x1b)
#define ansi_bold 			ansi_escape<<"[1m"
#define ansi_normal			ansi_escape<<"[0m"
#define ansi_red				ansi_escape<<"[31m"
#define ansi_green			ansi_escape<<"[32m"
#define ansi_blue				ansi_escape<<"[34m"
#endif // ansi_escape

#define ONE_OVER_SQRT12  0.288675

bool static DTrackHitSelector_cdchit_cmp(pair<double,const DCDCTrackHit *>a,
				      pair<double,const DCDCTrackHit *>b){
  if (a.second->wire->ring!=b.second->wire->ring) 
    return (a.second->wire->ring>b.second->wire->ring);
  return (a.first>b.first);
}
bool static DTrackHitSelector_fdchit_cmp(pair<double,const DFDCPseudo *>a,
				      pair<double,const DFDCPseudo *>b){
  if (a.second->wire->layer!=b.second->wire->layer) 
    return (a.second->wire->layer>b.second->wire->layer);
  return (a.first>b.first);
}


//---------------------------------
// DTrackHitSelectorALT1    (Constructor)
//---------------------------------
DTrackHitSelectorALT1::DTrackHitSelectorALT1(jana::JEventLoop *loop):DTrackHitSelector(loop)
{
	HS_DEBUG_LEVEL = 0;
	gPARMS->SetDefaultParameter("TRKFIT:HS_DEBUG_LEVEL", HS_DEBUG_LEVEL);

}

//---------------------------------
// ~DTrackHitSelectorALT1    (Destructor)
//---------------------------------
DTrackHitSelectorALT1::~DTrackHitSelectorALT1()
{

}

//---------------------------------
// GetCDCHits
//---------------------------------
void DTrackHitSelectorALT1::GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const
{
  // Vector of pairs storing the hit with the probability it is on the track
  vector<pair<double,const DCDCTrackHit*> >cdchits_tmp;

  /// Determine the probability that for each CDC hit that it came from the 
  /// track with the given trajectory.
  ///
  /// This will calculate a probability for each CDC hit that
  /// it came from the track represented by the given
  /// DReference trajectory. The probability is based on
  /// the residual between the distance of closest approach
  /// of the trajectory to the wire and the drift time for
  /// time-based tracks and the distance to the wire for
  /// wire-based tracks.
  
  // Calculate beta of particle.
  double my_mass=rt->GetMass();
  double one_over_beta =sqrt(1.0+my_mass*my_mass/rt->swim_steps[0].mom.Mag2());
  
  // The error on the residual. This will be different based on the
  // quality of the track and whether MULS is on or not etc.
  // In principle, this could also depend on the momentum parameters
  // of the track.
  double sigma;
  switch(fit_type){
  case kTimeBased:
    sigma = 0.8*ONE_OVER_SQRT12;
    break;
  case kWireBased:
    sigma = 1.6*ONE_OVER_SQRT12;
    break;
  case kHelical:
  default:
    sigma = 8.0*ONE_OVER_SQRT12;
  }
  
  // Low-momentum tracks are more poorly defined than high-momentum tracks.
  // We account for that here by increasing the error as a function of momentum
  double g = 0.350/sqrt(log(2.0)); // total guess
  double p_over_g=rt->swim_steps[0].mom.Mag()/g;
  //sigma *= 1.0 + exp(-pow(rt->swim_steps[0].mom.Mag()/g,2.0));
  sigma*=1.0+exp(-p_over_g*p_over_g);
  
  // Minimum probability of hit belonging to wire and still be accepted
  double MIN_HIT_PROB = 0.05;
  vector<const DCDCTrackHit*>::const_iterator iter;
  for(iter=cdchits_in.begin(); iter!=cdchits_in.end(); iter++){
    const DCDCTrackHit *hit = *iter;
    
    // Find the DOCA to this wire
    double s;
    double doca = rt->DistToRT(hit->wire, &s);
    
    // Get "measured" distance to wire. For time-based tracks
    // this is calculated from the drift time. For all other
    // tracks, this is assumed to be half a cell size
    double dist;
    if(fit_type == kTimeBased){
      // Distance using drift time
      // NOTE: Right now we assume pions for the TOF
      // and a constant drift velocity of 55um/ns
      double tof = s*one_over_beta/29.98;
      dist = (hit->tdrift - tof)*55E-4;
    }else{
      dist = 0.4; // =0.8/2.0; half cell-size
    }
    
    // For time-based and wire-based tracks, the fit was
    // weighted for multiple scattering by material times 
    // angle giving preference to the begining of the 
    // track. Take this into account here by enhancing the
    // error for hits further from the vertex
    double sigma_total = sigma;
    if(fit_type == kTimeBased || fit_type == kWireBased){
      sigma_total *= 1.0 + s/50.0; // double error at 50cm out (guess for now)
    }
    
    // Residual
    double resi = dist - doca;
    //double chisq = pow(resi/sigma_total, 2.0);
    double chisq=resi*resi/(sigma_total*sigma_total);
    
    // Use chi-sq probability function with Ndof=1 to calculate probability
    double probability = TMath::Prob(chisq, 1);
    if(probability>=MIN_HIT_PROB){
      pair<double,const DCDCTrackHit*>myhit;
      myhit.first=probability;
      myhit.second=hit;
      cdchits_tmp.push_back(myhit);
    }
    
    if(HS_DEBUG_LEVEL>10){
      _DBG_;
      if(probability>=MIN_HIT_PROB)jerr<<ansi_bold<<ansi_green;
      jerr<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" sigma="<<sigma_total<<" prob="<<probability<<endl;
      jerr<<ansi_normal;
    }
  }

  // Order according to ring number and probability, then put the hits in the 
  // output list with the following algorithm:  hits with the highest 
  // probability in a given ring are automatically put in the output list, 
  // but if there is more than one hit in a given ring, only those hits 
  // that are within +/-1 of the straw # of the most probable hit are added 
  // to the list.
  sort(cdchits_tmp.begin(),cdchits_tmp.end(),DTrackHitSelector_cdchit_cmp);
  int old_straw=1000,old_ring=1000;
  for (unsigned int i=0;i<cdchits_tmp.size();i++){
    if (cdchits_tmp[i].second->wire->ring!=old_ring || 
	abs(cdchits_tmp[i].second->wire->straw-old_straw)==1){
      cdchits_out.push_back(cdchits_tmp[i].second);   
    }
    old_straw=cdchits_tmp[i].second->wire->straw;
    old_ring=cdchits_tmp[i].second->wire->ring;
  }
}

//---------------------------------
// GetFDCHits
//---------------------------------
void DTrackHitSelectorALT1::GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const
{
  // Vector of pairs storing the hit with the probability it is on the track
  vector<pair<double,const DFDCPseudo*> >fdchits_tmp;

  /// Determine the probability that for each FDC hit that it came from the 
  /// track with the given trajectory.
  ///
  /// This will calculate a probability for each FDC hit that
  /// it came from the track represented by the given
  /// DReference trajectory. The probability is based on
  /// the residual between the distance of closest approach
  /// of the trajectory to the wire and the drift time
  /// and the distance along the wire.
  
  // Calculate beta of particle assuming its a pion for now. If the
  // particles is really a proton or an electron, the residual
  // calculated below will only be off by a little.
  double my_mass=rt->GetMass();
  double one_over_beta =sqrt(1.0+my_mass*my_mass/rt->swim_steps[0].mom.Mag2());
  
  // The error on the residual. This will be different based on the
  // quality of the track and whether MULS is on or not etc.
  // In principle, this could also depend on the momentum parameters
  // of the track.
  double sigma_anode = 0.5*ONE_OVER_SQRT12;
  double sigma_cathode = 0.5*ONE_OVER_SQRT12;
  switch(fit_type){
  case kTimeBased:
    sigma_anode = 0.5*ONE_OVER_SQRT12;
    sigma_cathode = 0.5*ONE_OVER_SQRT12;
    break;
  case kWireBased:
    sigma_anode = ONE_OVER_SQRT12;
    sigma_cathode = ONE_OVER_SQRT12;
    break;
  case kHelical:
  default:
    sigma_anode = 5.0*ONE_OVER_SQRT12;
    sigma_cathode = 5.0*ONE_OVER_SQRT12;
  }
  
  // Scale the errors up by 1/p for p<1GeV/c. This accounts for the poorer knowledge
  // of the track parameters contained in the candidate for low momentum tracks
  double ptot = rt->swim_steps[0].mom.Mag();
  if(ptot<1.0){
    sigma_anode /= ptot;
    sigma_cathode /= ptot;
  }
  
  // Higher mass particles (protons) lose energy quicker and therefore have even
  // less helical shapes for candidates and larger errors due to eloss straggling.
  // These next 2 lines are here mainly to improve the situation for protons so
  // we take a stab at scaling up the errors by the number of pion masses since
  // this worked reasonably well for pions before.
  double mass_scale=my_mass/0.13957;
  sigma_anode *= mass_scale;
  sigma_cathode *= mass_scale;
  
  // Minimum probability of hit belonging to wire and still be accepted
  double MIN_HIT_PROB = 0.01;
  
  vector<const DFDCPseudo*>::const_iterator iter;
  for(iter=fdchits_in.begin(); iter!=fdchits_in.end(); iter++){
    const DFDCPseudo *hit = *iter;
    
    // Find the DOCA to this wire
    double s;
    double doca = rt->DistToRT(hit->wire, &s);  
       
    // Get "measured" distance to wire. For time-based tracks
    // this is calculated from the drift time. For all other
    // tracks, this is assumed to be half a cell size
    double dist;
    if(kTimeBased){
      // Distance using drift time
      // NOTE: Right now we assume pions for the TOF
      // and a constant drift velocity of 55um/ns
      double tof = s*one_over_beta/29.98;
      dist = (hit->time - tof)*55E-4;
    }else{
      dist = 0.25; //= 0.5/2.0; half cell-size
    }
    
    // Anode Residual
    double resi = dist - doca;		
    
    // Cathode Residual
    double u=rt->GetLastDistAlongWire();
    double resic = u - hit->s;
    
    // Probability of this hit being on the track
    //double chisq = pow(resi/sigma_anode, 2.0) + pow(resic/sigma_cathode, 2.0);
    double chisq
      =resi*resi/(sigma_anode*sigma_anode)
      +resic*resic/(sigma_cathode*sigma_cathode);
      
    double probability = TMath::Prob(chisq, 2);
    if(probability>=MIN_HIT_PROB){
      pair<double,const DFDCPseudo*>myhit;
      myhit.first=probability;
      myhit.second=hit;
      fdchits_tmp.push_back(myhit);
    }
    if(HS_DEBUG_LEVEL>10){
      _DBG_;
      if(probability>=MIN_HIT_PROB)jerr<<ansi_bold<<ansi_blue;
      jerr<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" resic="<<resic<<" chisq="<<chisq<<" prob="<<probability<<endl;
      jerr<<ansi_normal;
    }
  }
  // Order according to layer number and probability,then put the hits in the 
  // output list with the following algorithm:  hits with the highest 
  // probability in a given layer are automatically put in the output list, 
  // but if there is more than one hit in a given layer, only those hits 
  // that are within +/-1 of the wire # of the most probable hit are added 
  // to the list.
  sort(fdchits_tmp.begin(),fdchits_tmp.end(),DTrackHitSelector_fdchit_cmp);
  int old_layer=1000,old_wire=1000;
  for (unsigned int i=0;i<fdchits_tmp.size();i++){
    if (fdchits_tmp[i].second->wire->layer!=old_layer || 
	abs(fdchits_tmp[i].second->wire->wire-old_wire)==1){
      fdchits_out.push_back(fdchits_tmp[i].second);   
    }
    old_wire=fdchits_tmp[i].second->wire->wire;
    old_layer=fdchits_tmp[i].second->wire->layer;
  }
}
