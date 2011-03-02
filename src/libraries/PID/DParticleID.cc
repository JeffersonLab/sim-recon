// $Id$
//
//    File: DParticleID.cc
// Created: Mon Feb 28 14:48:56 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID.h"
#include "TRACKING/DReferenceTrajectory.h"

// Routine for sorting dEdx data
bool static DParticleID_dedx_cmp(DParticleID::dedx_t a,DParticleID::dedx_t b){
  double dEdx1=a.dE/a.dx;
  double dEdx2=b.dE/b.dx;
  return dEdx1<dEdx2;
}

// Routine for sorting hypotheses according to FOM
bool static DParticleID_hypothesis_cmp(const DTrackTimeBased *a,
				       const DTrackTimeBased *b){
  return (a->FOM>b->FOM);
}


//---------------------------------
// DParticleID    (Constructor)
//---------------------------------
DParticleID::DParticleID(JEventLoop *loop)
{
  
  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  if(!dapp){
    _DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
  }
  RootGeom = dapp->GetRootGeom();

  // Get material properties for chamber gas
  double rho_Z_over_A_LnI=0,radlen=0;
  RootGeom->FindMat("CDchamberGas",mRhoZoverAGas,rho_Z_over_A_LnI,
		    radlen);
  mLnIGas=rho_Z_over_A_LnI/mRhoZoverAGas;
  mKRhoZoverAGas=0.1535*mRhoZoverAGas;

}

//---------------------------------
// ~DParticleID    (Destructor)
//---------------------------------
DParticleID::~DParticleID()
{

}

// Group fitted tracks according to candidate id
jerror_t DParticleID::GroupTracks(vector<const DTrackTimeBased *> tracks,
			      vector<vector<const DTrackTimeBased*> >&grouped_tracks){ 
  if (tracks.size()==0) return RESOURCE_UNAVAILABLE;
 
  JObject::oid_t old_id=tracks[0]->candidateid;
  vector<const DTrackTimeBased *>hypotheses;
  for (unsigned int i=0;i<tracks.size();i++){
    const DTrackTimeBased *track=tracks[i];
    
    if (old_id != track->candidateid){
      // Sort according to FOM
      sort(hypotheses.begin(),hypotheses.end(),DParticleID_hypothesis_cmp);

      // Add to the grouped_tracks vector
      grouped_tracks.push_back(hypotheses);

      // Clear the hypothesis list for the next track
      hypotheses.clear();
    }
    hypotheses.push_back(track);
    old_id=track->candidateid;
  }
  // Final set
  sort(hypotheses.begin(),hypotheses.end(),DParticleID_hypothesis_cmp);
  grouped_tracks.push_back(hypotheses);

  return NOERROR;
}



// Calculate the most probable energy loss per unit length in units of 
// MeV/cm in the FDC or CDC gas for a particle of momentum p and mass "mass"
double DParticleID::GetMostProbabledEdx(double p,double mass,double dx){
  double betagamma=p/mass;
  double beta2=1./(1.+1./betagamma/betagamma);
  if (beta2<1e-6) beta2=1e-6;
  
  // Electron mass 
  double Me=0.000511; //GeV

  // First (non-logarithmic) term in Bethe-Bloch formula
  double mean_dedx=mKRhoZoverAGas/beta2;
 
  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
  return mean_dedx*(log(mean_dedx*dx/1000.)
		    -log((1.-beta2)/2./Me/beta2)-2.*mLnIGas-beta2+0.198);
}

// Empirical form for sigma/mean for gaseous detectors with num_dedx 
// samples and sampling thickness path_length.  Assumes that the number of 
// hits has already been converted from an (unsigned) int to a double.
#define TCUT 100.e-6 // energy cut for bethe-bloch 
double DParticleID::GetdEdxSigma(double num_hits,double p,double mass,
				  double mean_path_length){
  // kinematic quantities
  double betagamma=p/mass;
  double betagamma2=betagamma*betagamma;
  double beta2=1./(1.+1./betagamma2);
  if (beta2<1e-6) beta2=1e-6;

  double Me=0.000511; //GeV
  double m_ratio=Me/mass;
  double two_Me_betagamma_sq=2.*Me*betagamma2;
  double Tmax
    =two_Me_betagamma_sq/(1.+2.*sqrt(1.+betagamma2)*m_ratio+m_ratio*m_ratio);
  // Energy truncation for knock-on electrons
  double T0=(Tmax>TCUT)?TCUT:Tmax;
  
  // Bethe-Bloch
  double mean_dedx=mKRhoZoverAGas/beta2
    *(log(two_Me_betagamma_sq*T0)-2.*mLnIGas-beta2*(1.+T0/Tmax));
  
  return 0.41*mean_dedx*pow(num_hits,-0.43)*pow(mean_path_length,-0.32);
  //return 0.41*mean_dedx*pow(double(num_hits),-0.5)*pow(mean_path_length,-0.32);
}

// Compute the energy losses and the path lengths in the chambers for each hit 
// on the track. Returns a list of dE and dx pairs with the momentum at the 
// hit.
jerror_t DParticleID::GetdEdx(const DTrackTimeBased *track,
			       vector<dedx_t>&dEdx_list){
  // Position and momentum
  DVector3 pos,mom;
  
  //dE and dx pairs
  pair<double,double>de_and_dx;

  // We cast away the const-ness of the reference trajectory so that we can use the DisToRT method
  DReferenceTrajectory *my_rt=const_cast<DReferenceTrajectory*>(track->rt);

  //Get the list of cdc hits used in the fit
  vector<const DCDCTrackHit*>cdchits;
  track->GetT(cdchits);

  // Loop over cdc hits
  for (unsigned int i=0;i<cdchits.size();i++){
    my_rt->DistToRT(cdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);

    // Create the dE,dx pair from the position and momentum using a helical approximation for the path 
    // in the straw and keep track of the momentum in the active region of the detector
    if (CalcdEdxHit(mom,pos,cdchits[i],de_and_dx)==NOERROR){
      dEdx_list.push_back(dedx_t(de_and_dx.first,de_and_dx.second,
				 mom.Mag()));
    }
  }
  
  //Get the list of fdc hits used in the fit
  vector<const DFDCPseudo*>fdchits;
  track->GetT(fdchits);

  // loop over fdc hits
  for (unsigned int i=0;i<fdchits.size();i++){
    my_rt->DistToRT(fdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);
   
    double gas_thickness=1.0; // cm
    dEdx_list.push_back(dedx_t(1000.*fdchits[i]->dE,
			       gas_thickness/cos(mom.Theta()),
			       mom.Mag()));
  }
    
  // Sort the dEdx entries from smallest to largest
  sort(dEdx_list.begin(),dEdx_list.end(),DParticleID_dedx_cmp);  
 
  return NOERROR;
}





// Calculate the path length for a single hit in a straw and return ds and the 
// energy deposition in the straw.  It returns dE as the first element in the 
// dedx pair and ds as the second element in the dedx pair.
jerror_t DParticleID::CalcdEdxHit(const DVector3 &mom,
				  const DVector3 &pos,
				  const DCDCTrackHit *hit,
				  pair <double,double> &dedx){
  if (hit==NULL || hit->wire==NULL) return RESOURCE_UNAVAILABLE;
  
  // Track direction parameters
  double phi=mom.Phi();
  double lambda=M_PI_2-mom.Theta();
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double tanl=tan(lambda);
  
  //Position relative to wire origin
  double dz=pos.z()-hit->wire->origin.z();
  double dx=pos.x()-hit->wire->origin.x();
  double dy=pos.y()-hit->wire->origin.y();
  
  // square of straw radius
  double rs2=0.776*0.776;
  
  // Useful temporary variables related to the direction of the wire
  double ux=hit->wire->udir.x();
  double uy=hit->wire->udir.y();
  double uz=hit->wire->udir.z();
  double A=1.-ux*ux;
  double B=-2.*ux*uy;
  double C=-2.*ux*uz;
  double D=-2.*uy*uz;
  double E=1.-uy*uy;
  double F=1.-uz*uz;  

  // The path length in the straw is given by  s=sqrt(b*b-4*a*c)/a/cosl.
  // a, b, and c follow.
  double a=A*cosphi*cosphi+B*cosphi*sinphi+C*cosphi*tanl+D*sinphi*tanl
    +E*sinphi*sinphi+F*tanl*tanl;
  double b=2.*A*dx*cosphi+B*dx*sinphi+B*dy*cosphi+C*dx*tanl+C*cosphi*dz
    +D*dy*tanl+D*sinphi*dz+2.*E*dy*sinphi+2.*F*dz*tanl;
  double c=A*dx*dx+B*dx*dy+C*dx*dz+D*dy*dz+E*dy*dy+F*dz*dz-rs2;
  
  // Check for valid arc length and compute dEdx
  double temp=b*b-4.*a*c;
  if (temp>0){
    double cosl=fabs(cos(lambda));
    //    double gas_density=0.0018;

    // arc length and energy deposition
    //dedx.second=gas_density*sqrt(temp)/a/cosl; // g/cm^2
    dedx.second=sqrt(temp)/a/cosl;
    dedx.first=1000.*hit->dE; //MeV

    return NOERROR;
  }
  
  return VALUE_OUT_OF_RANGE;
}
  
