// $Id$
//
//    File: DTrackHitSelectorALT1.cc
// Created: Fri Feb  6 08:22:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include <TROOT.h>

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
	MAKE_DEBUG_TREES = false;

	gPARMS->SetDefaultParameter("TRKFIT:HS_DEBUG_LEVEL", HS_DEBUG_LEVEL, "Debug verbosity level for hit selector used in track fitting (0=no debug messages)");
	gPARMS->SetDefaultParameter("TRKFIT:MAKE_DEBUG_TREES", MAKE_DEBUG_TREES, "Create a TTree with debugging info on hit selection for the FDC and CDC");
	
	cdchitsel = NULL;
	fdchitsel = NULL;
	if(MAKE_DEBUG_TREES){
		loop->GetJApplication()->Lock();
		
		 cdchitsel= (TTree*)gROOT->FindObject("cdchitsel");
		if(!cdchitsel){
			cdchitsel = new TTree("cdchitsel", "CDC Hit Selector");
			cdchitsel->Branch("H", &cdchitdbg, "fit_type/I:p/F:theta:mass:sigma:mom_factor:x:y:z:s:s_factor:itheta02:itheta02s:itheta02s2:dist:doca:resi:sigma_total:chisq:prob");
		}else{
			_DBG__;
			jerr<<" !!! WARNING !!!"<<endl;
			jerr<<"It appears that the cdchitsel TTree is already defined."<<endl;
			jerr<<"This is probably means you are running with multiple threads."<<endl;
			jerr<<"To avoid complication, filling of the hit selector debug"<<endl;
			jerr<<"trees will be disabled for this thread."<<endl;
			_DBG__;
			cdchitsel = NULL;
		}

		fdchitsel= (TTree*)gROOT->FindObject("fdchitsel");
		if(!fdchitsel){
			fdchitsel = new TTree("fdchitsel", "FDC Hit Selector");
			fdchitsel->Branch("H", &fdchitdbg, "fit_type/I:p/F:theta:mass:sigma_anode:sigma_cathode:mom_factor_anode:mom_factor_cathode:x:y:z:s:s_factor_anode:s_factor_cathode:itheta02:itheta02s:itheta02s2:dist:doca:resi:u:u_cathodes:resic:sigma_anode_total:sigma_cathode_total:chisq:prob:prob_anode:prob_cathode:pull_anode:pull_cathode");
		}else{
			_DBG__;
			jerr<<" !!! WARNING !!!"<<endl;
			jerr<<"It appears that the fdchitsel TTree is already defined."<<endl;
			jerr<<"This is probably means you are running with multiple threads."<<endl;
			jerr<<"To avoid complication, filling of the hit selector debug"<<endl;
			jerr<<"trees will be disabled for this thread."<<endl;
			_DBG__;
			fdchitsel = NULL;
		}

		loop->GetJApplication()->Unlock();
	}
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
  double p = rt->swim_steps[0].mom.Mag();
  double my_mass=rt->GetMass();
  double one_over_beta =sqrt(1.0+my_mass*my_mass/(p*p));
  
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
  double p_over_g=p/g;
  //sigma *= 1.0 + exp(-pow(rt->swim_steps[0].mom.Mag()/g,2.0));
  double mom_factor = 1.0+exp(-p_over_g*p_over_g);
  sigma *= mom_factor;
  
  int old_straw=1000,old_ring=1000;

  // Minimum probability of hit belonging to wire and still be accepted
  double MIN_HIT_PROB = 0.05;
  vector<const DCDCTrackHit*>::const_iterator iter;
  for(iter=cdchits_in.begin(); iter!=cdchits_in.end(); iter++){
    const DCDCTrackHit *hit = *iter;

    // Skip hit if it is on the same wire as the previous hit
    if (fit_type!=kTimeBased){
      if (hit->wire->ring == old_ring && hit->wire->straw==old_straw) continue;
      old_ring=hit->wire->ring;
      old_straw=hit->wire->straw;
    }
    // Find the DOCA to this wire
    double s;
    double doca = rt->DistToRT(hit->wire, &s);
    if(!finite(doca))
      continue;

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
	double s_factor = 1.0;
    if(fit_type == kTimeBased || fit_type == kWireBased){
      s_factor = 1.0 + s/50.0; // double error at 50cm out (guess for now)
      sigma_total *= s_factor;
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

	// Optionally fill debug tree
    if(cdchitsel){
		DVector3 pos = rt->GetLastDOCAPoint();
		const DReferenceTrajectory::swim_step_t *last_step = rt->GetLastSwimStep();

		cdchitdbg.fit_type = fit_type;
		cdchitdbg.p = p;
		cdchitdbg.theta = rt->swim_steps[0].mom.Theta();
		cdchitdbg.mass = my_mass;
		cdchitdbg.sigma = sigma;
		cdchitdbg.mom_factor = mom_factor;
        cdchitdbg.x = pos.X();
        cdchitdbg.y = pos.Y();
        cdchitdbg.z = pos.Z();
        cdchitdbg.s = s;
        cdchitdbg.s_factor = s_factor;
        cdchitdbg.itheta02 = last_step->itheta02;
        cdchitdbg.itheta02s = last_step->itheta02s;
        cdchitdbg.itheta02s2 = last_step->itheta02s2;
        cdchitdbg.dist = dist;
        cdchitdbg.doca = doca;
        cdchitdbg.resi = resi;
        cdchitdbg.sigma_total = sigma_total;
        cdchitdbg.chisq = chisq;
        cdchitdbg.prob = probability;
		
		cdchitsel->Fill();
		
		static bool printed_first = false;
		if(!printed_first){
			_DBG_<<"=== Printing first entry for CDC hit selector debug tree ==="<<endl;
			_DBG_<<"   fit_type = "<<cdchitdbg.fit_type<<endl;
			_DBG_<<"          p = "<<cdchitdbg.p<<endl;
			_DBG_<<"      theta = "<<cdchitdbg.theta<<endl;
			_DBG_<<"       mass = "<<cdchitdbg.mass<<endl;
			_DBG_<<"      sigma = "<<cdchitdbg.sigma<<endl;
			_DBG_<<" mom_factor = "<<cdchitdbg.mom_factor<<endl;
			_DBG_<<"          x = "<<cdchitdbg.x<<endl;
			_DBG_<<"          y = "<<cdchitdbg.y<<endl;
			_DBG_<<"          z = "<<cdchitdbg.z<<endl;
			_DBG_<<"          s = "<<cdchitdbg.s<<endl;
			_DBG_<<"   s_factor = "<<cdchitdbg.s_factor<<endl;
			_DBG_<<"   itheta02 = "<<cdchitdbg.itheta02<<endl;
			_DBG_<<"  itheta02s = "<<cdchitdbg.itheta02s<<endl;
			_DBG_<<" itheta02s2 = "<<cdchitdbg.itheta02s2<<endl;
			_DBG_<<"       dist = "<<cdchitdbg.dist<<endl;
			_DBG_<<"       doca = "<<cdchitdbg.doca<<endl;
			_DBG_<<"       resi = "<<cdchitdbg.resi<<endl;
			_DBG_<<"sigma_total = "<<cdchitdbg.sigma_total<<endl;
			_DBG_<<"      chisq = "<<cdchitdbg.chisq<<endl;
			_DBG_<<"       prob = "<<cdchitdbg.prob<<endl;
			
			printed_first = true;
		}
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
  old_straw=1000,old_ring=1000;
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
  double p = rt->swim_steps[0].mom.Mag();
  //double one_over_beta =sqrt(1.0+my_mass*my_mass/(p*p));
  
  // Scale the errors up by 1/p for p<1GeV/c. This accounts for the poorer knowledge
  // of the track parameters contained in the candidate for low momentum tracks
  double mom_factor = 1.0;
  if(p<1.0)mom_factor = 1.0/p;
  
  // Higher mass particles (protons) lose energy quicker and therefore have even
  // less helical shapes for candidates and larger errors due to eloss straggling.
  // These next 2 lines are here mainly to improve the situation for protons so
  // we take a stab at scaling up the errors by the number of pion masses since
  // this worked reasonably well for pions before.
  double mass_factor=my_mass/0.13957;
  
  // Minimum probability of hit belonging to wire and still be accepted
  double MIN_HIT_PROB = 0.01;
  
  vector<const DFDCPseudo*>::const_iterator iter;
  for(iter=fdchits_in.begin(); iter!=fdchits_in.end(); iter++){
    const DFDCPseudo *hit = *iter;

    // Find the DOCA to this wire
    double s;
    DVector3 fdc_pos(hit->xy.X(),hit->xy.Y(),hit->wire->origin.z());
    double doca=rt->DistToRT(fdc_pos,&s);
    if(!finite(doca))
      continue;
    double fdc_var=mom_factor*mom_factor*mass_factor*mass_factor*(1.0*1.0+0.3*0.3)/12.;

    if (fit_type==kHelical) fdc_var*=25.;
    double chisq=doca*doca/fdc_var;
    double MAX_HIT_DOCA = 30;
    double probability = (doca<MAX_HIT_DOCA)?TMath::Prob(chisq, 1):0.;
    if(probability>=MIN_HIT_PROB){
      pair<double,const DFDCPseudo*>myhit;
      myhit.first=probability;
      myhit.second=hit;
      fdchits_tmp.push_back(myhit);
    }

    // Optionally fill debug tree
    /*
    if(fdchitsel){
      DVector3 pos = rt->GetLastDOCAPoint();
      const DReferenceTrajectory::swim_step_t *last_step = rt->GetLastSwimStep();
      

      fdchitdbg.fit_type = fit_type;
      fdchitdbg.p = p;
      fdchitdbg.theta = rt->swim_steps[0].mom.Theta();
      fdchitdbg.mass = my_mass;
      fdchitdbg.sigma_anode = sigma_anode;
      fdchitdbg.sigma_cathode = sigma_cathode;
      fdchitdbg.mom_factor_anode = mom_factor;
      fdchitdbg.mom_factor_cathode = mom_factor;
        fdchitdbg.x = pos.X();
        fdchitdbg.y = pos.Y();
        fdchitdbg.z = pos.Z();
        fdchitdbg.s = s;
        fdchitdbg.s_factor_anode = 1.0;
        fdchitdbg.s_factor_cathode = 1.0;
        fdchitdbg.itheta02 = last_step->itheta02;
        fdchitdbg.itheta02s = last_step->itheta02s;
        fdchitdbg.itheta02s2 = last_step->itheta02s2;
        fdchitdbg.dist = dist;
        fdchitdbg.doca = doca;
        fdchitdbg.resi = resi;
        fdchitdbg.u = u;
        fdchitdbg.u_cathodes = u_cathodes;
        fdchitdbg.resic = resic;
        fdchitdbg.sigma_anode_total = sigma_anode_total;
        fdchitdbg.sigma_cathode_total = sigma_cathode_total;
        fdchitdbg.chisq = chisq;
        fdchitdbg.prob = probability;
		fdchitdbg.prob_anode = TMath::Prob(pull_anode*pull_anode, 1);
        fdchitdbg.prob_cathode = TMath::Prob(pull_cathode*pull_cathode, 1);
		fdchitdbg.pull_anode = pull_anode;
		fdchitdbg.pull_cathode = pull_cathode;
//_DBG_<<"fdchitdbg.pull_anode="<<fdchitdbg.pull_anode<<" pull_anode="<<pull_anode<<" resi/sigma_anode_total="<<resi/sigma_anode_total<<endl;
		
		fdchitsel->Fill();
	}
    */

    if(HS_DEBUG_LEVEL>10){
      _DBG_;
      if(probability>=MIN_HIT_PROB)jerr<<ansi_bold<<ansi_blue;
      jerr<<"s="<<s<<" doca="<<doca<<" chisq="<<chisq<<" prob="<<probability<<endl;
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
