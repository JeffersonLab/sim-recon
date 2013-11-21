// $Id$
//
//    File: DParticleID.cc
// Created: Mon Feb 28 14:48:56 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID.h"
#include <TRACKING/DTrackFitter.h>
#include <TMath.h>
#include "FCAL/DFCALGeometry.h"

#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

#define C_EFFECTIVE 15. // start counter light propagation speed
#define ATTEN_LENGTH 150.   // Start counter attenuation length
#define OUT_OF_TIME_CUT 200.

// Routine for sorting dEdx data
bool static DParticleID_dedx_cmp(DParticleID::dedx_t a,DParticleID::dedx_t b){
  return a.dEdx < b.dEdx;
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
  dRFBunchFrequency = 2.004;

  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  if(!dapp){
    _DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
  }
  const DRootGeom *RootGeom = dapp->GetRootGeom(loop->GetJEvent().GetRunNumber());
  bfield = dapp->GetBfield(); 
  stepper= new DMagneticFieldStepper(bfield);

  // Get material properties for chamber gas
  double rho_Z_over_A_LnI=0,radlen=0;

  RootGeom->FindMat("CDchamberGas", dRhoZoverA_CDC, rho_Z_over_A_LnI, radlen);
  dLnI_CDC = rho_Z_over_A_LnI/dRhoZoverA_CDC;
  dKRhoZoverA_CDC = 0.1535E-3*dRhoZoverA_CDC;

  RootGeom->FindMat("FDchamberGas", dRhoZoverA_FDC, rho_Z_over_A_LnI, radlen);
  dLnI_FDC = rho_Z_over_A_LnI/dRhoZoverA_FDC;
  dKRhoZoverA_FDC = 0.1535E-3*dRhoZoverA_FDC;

  RootGeom->FindMat("Scintillator", dRhoZoverA_Scint, rho_Z_over_A_LnI, radlen);
  dLnI_Scint = rho_Z_over_A_LnI/dRhoZoverA_Scint;
  dKRhoZoverA_Scint = 0.1535E-3*dRhoZoverA_Scint;


  // Get the geometry
  geom = dapp->GetDGeometry(loop->GetJEvent().GetRunNumber());

  vector<double>sc_origin;
  bool got_sc=geom->Get("//posXYZ[@volume='StartCntr']/@X_Y_Z",sc_origin);
  if (got_sc){  // Check if Start Counter geometry is present
    double num_paddles;
    geom->Get("//mposPhi[@volume='STRC']/@ncopy",num_paddles); 
    dSCdphi=M_TWO_PI/num_paddles;
    
    double Phi0;
    geom->Get("///mposPhi[@volume='STRC']/@Phi0",Phi0);
    dSCphi0=Phi0*M_PI/180.;

    vector<vector<double> > sc_rioz;
    geom->GetMultiple("//pgon[@name='STRC']/polyplane/@Rio_Z", sc_rioz);

    for (unsigned int k=0;k<sc_rioz.size()-1;k++){
      if(sc_origin.size() < 3)continue; // in case start counter is comment out in XML
      if(sc_rioz[k].size() < 3)continue; // in case start counter is comment out in XML
      DVector3 pos((sc_rioz[k][0]+sc_rioz[k][1])/2.,0.,sc_rioz[k][2]+sc_origin[2]);
      DVector3 dir(sc_rioz[k+1][2]-sc_rioz[k][2],0,
		   -sc_rioz[k+1][0]+sc_rioz[k][0]);
      dir.SetMag(1.);
      
      sc_pos.push_back(pos);
      sc_norm.push_back(dir);    
    }
  
    if(sc_pos.size() >= 1){// in case start counter is comment out in XML
      // sc_leg_tcor=(sc_light_guide[2]-sc_pos[0].z())/C_EFFECTIVE;
      sc_leg_tcor=-sc_pos[0].z()/C_EFFECTIVE;
      double theta=sc_norm[sc_norm.size()-1].Theta();
      sc_angle_cor=1./cos(M_PI-theta);
    }
  }

  //Get calibration constants
  map<string, double> locPIDParams;
  if ( !loop->GetCalib("PID/photon_track_matching", locPIDParams)){
    cout<<"DParticleID: loading values from PID data base"<<endl;
    DELTA_R_BCAL = locPIDParams["DELTA_R_BCAL"];
    DELTA_R_FCAL = locPIDParams["DELTA_R_FCAL"];
  } else {
    cout << "DParticleID: Error loading values from PID data base" <<endl;
    DELTA_R_BCAL = 15.0;
    DELTA_R_FCAL = 15.0;
  }

	dTargetZCenter = 0.0;
	geom->GetTargetZ(dTargetZCenter);
}

//---------------------------------
// ~DParticleID    (Destructor)
//---------------------------------
DParticleID::~DParticleID()
{
  if (stepper) delete stepper;
}

// Group fitted tracks according to candidate id
jerror_t DParticleID::GroupTracks(vector<const DTrackTimeBased *> &tracks,
			      vector<vector<const DTrackTimeBased*> >&grouped_tracks) const{ 
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

// Compute the energy losses and the path lengths in the chambers for each hit 
// on the track. Returns a list of dE and dx pairs with the momentum at the 
// hit.
jerror_t DParticleID::GetDCdEdxHits(const DTrackTimeBased *track, vector<dedx_t>& dEdxHits_CDC, vector<dedx_t>& dEdxHits_FDC) const{
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
    double locReturnValue = my_rt->DistToRT(cdchits[i]->wire);
    if(!((locReturnValue >= 0.0) || (locReturnValue <= 0.0)))
      continue; //NaN

    my_rt->GetLastDOCAPoint(pos, mom);

    // Create the dE,dx pair from the position and momentum using a helical approximation for the path 
    // in the straw and keep track of the momentum in the active region of the detector
    if (CalcdEdxHit(mom,pos,cdchits[i],de_and_dx)==NOERROR)
      dEdxHits_CDC.push_back(dedx_t(de_and_dx.first, de_and_dx.second, mom.Mag()));
  }
  
  //Get the list of fdc hits used in the fit
  vector<const DFDCPseudo*>fdchits;
  track->GetT(fdchits);

  // loop over fdc hits
  for (unsigned int i=0;i<fdchits.size();i++){
    double locReturnValue = my_rt->DistToRT(fdchits[i]->wire);
    if(!((locReturnValue >= 0.0) || (locReturnValue <= 0.0)))
      continue; //NaN
    my_rt->GetLastDOCAPoint(pos, mom);
   
    double gas_thickness = 1.0; // cm
    dEdxHits_FDC.push_back(dedx_t(fdchits[i]->dE, gas_thickness/cos(mom.Theta()), mom.Mag()));
  }
    
  // Sort the dEdx entries from smallest to largest
  sort(dEdxHits_FDC.begin(),dEdxHits_FDC.end(),DParticleID_dedx_cmp);  
  sort(dEdxHits_CDC.begin(),dEdxHits_CDC.end(),DParticleID_dedx_cmp);  
 
  return NOERROR;
}

jerror_t DParticleID::CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const{
	vector<dedx_t> locdEdxHits_CDC, locdEdxHits_FDC;
	jerror_t locReturnStatus = GetDCdEdxHits(locTrackTimeBased, locdEdxHits_CDC, locdEdxHits_FDC);
	if(locReturnStatus != NOERROR){
		locdEdx_FDC = numeric_limits<double>::quiet_NaN();
		locdx_FDC = numeric_limits<double>::quiet_NaN();
		locNumHitsUsedFordEdx_FDC = 0;
		locdEdx_CDC = numeric_limits<double>::quiet_NaN();
		locdx_CDC = numeric_limits<double>::quiet_NaN();
		locNumHitsUsedFordEdx_CDC = 0;
		return locReturnStatus;
	}
	return CalcDCdEdx(locTrackTimeBased, locdEdxHits_CDC, locdEdxHits_FDC, locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC, locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC);
}

jerror_t DParticleID::CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, const vector<dedx_t>& locdEdxHits_CDC, const vector<dedx_t>& locdEdxHits_FDC, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const{
	locdx_CDC = 0.0;
	locdEdx_CDC = 0.0;
	locNumHitsUsedFordEdx_CDC = locdEdxHits_CDC.size()/2;
	if(locNumHitsUsedFordEdx_CDC > 0){
		for(unsigned int loc_i = 0; loc_i < locNumHitsUsedFordEdx_CDC; ++loc_i){
			locdEdx_CDC += locdEdxHits_CDC[loc_i].dE; //weight is ~ #e- (scattering sites): dx!
			locdx_CDC += locdEdxHits_CDC[loc_i].dx;
		}
		locdEdx_CDC /= locdx_CDC;
	}

	locdx_FDC = 0.0;
	locdEdx_FDC = 0.0;
	locNumHitsUsedFordEdx_FDC = locdEdxHits_FDC.size()/2;
	if(locNumHitsUsedFordEdx_FDC > 0){
		for(unsigned int loc_i = 0; loc_i < locNumHitsUsedFordEdx_FDC; ++loc_i){
			locdEdx_FDC += locdEdxHits_FDC[loc_i].dE; //weight is ~ #e- (scattering sites): dx!
			locdx_FDC += locdEdxHits_FDC[loc_i].dx;
		}
		locdEdx_FDC /= locdx_FDC; //weight is dx/dx_total
	}
	return NOERROR;
}

// Calculate the path length for a single hit in a straw and return ds and the 
// energy deposition in the straw.  It returns dE as the first element in the 
// dedx pair and ds as the second element in the dedx pair.
jerror_t DParticleID::CalcdEdxHit(const DVector3 &mom,
				  const DVector3 &pos,
				  const DCDCTrackHit *hit,
				  pair <double,double> &dedx) const{
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
    dedx.first=hit->dE; //GeV

    return NOERROR;
  }
  
  return VALUE_OUT_OF_RANGE;
}
  
//------------------
// MatchToTOF
//------------------
// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match. 
//
// NOTE:  an initial guess for tproj is expected as input so that out-of-time 
// hits can be skipped
jerror_t DParticleID::MatchToTOF(const DReferenceTrajectory *rt, DTrackFitter::fit_type_t fit_type, vector<const DTOFPoint*>&tof_points, double &tproj, unsigned int &tof_match_id, double &locPathLength, double &locFlightTime,pair<double,double>*dEdx) const{
  //tproj=NaN;
  tof_match_id=0;
  if (tof_points.size()==0){
    tproj=numeric_limits<double>::quiet_NaN();
    return RESOURCE_UNAVAILABLE;
  }

  double dmin=10000.;
  // loop over tof points
  double locTempPathLength=0.;
  for (unsigned int k=0;k<tof_points.size();k++){

    // Get the TOF cluster position and normal vector for the TOF plane
    DVector3 tof_pos=tof_points[k]->pos;
    DVector3 norm(0,0,1);
    DVector3 proj_pos,proj_mom;
    
    // Find the distance of closest approach between the track trajectory
    // and the tof cluster position, looking for the minimum
    double locTempFlightTime=0.;
    if (rt->GetIntersectionWithPlane(tof_pos,norm,proj_pos,proj_mom,&locTempPathLength,&locTempFlightTime)==NOERROR){
      double d=(tof_pos-proj_pos).Mag();

      // Check that the hit is not out of time with respect to the track
      if (fabs(tof_points[k]->t - locTempFlightTime - tproj) > OUT_OF_TIME_CUT) continue;
      
      if (d<dmin){
	dmin=d;
	if (dEdx!=NULL){
	  double p=proj_mom.Mag();
	  double dx=2.54*p/proj_mom.Dot(norm);
	  double dE=tof_points[k]->dE;
	  double dE_sigma,dE_mp;
	  GetScintMPdEandSigma(p,rt->GetMass(),dx,dE_mp,dE_sigma);
	  
	  dEdx->first=dE/dx;
	  dEdx->second=(dE-dE_mp)/dE_sigma;

	}
	locPathLength = locTempPathLength;
	locFlightTime = locTempFlightTime;
	tof_match_id=k;
      }
    }
  }
    
  // Check for a match 
  //  double p=rt->swim_steps[0].mom.Mag();
  double match_cut=0.;
  match_cut = 6.15; //current dPositionMatchCut_DoubleEnded variable in DTOFPoint_factory.cc
  //  if (fit_type==DTrackFitter::kTimeBased) match_cut=3.624+0.488/p;
  //  else match_cut=4.0+0.488/p;
  if (dmin<match_cut){
    // Projected time at the target
    tproj=tof_points[tof_match_id]->t - locFlightTime;
      
    return NOERROR;
  }
  tproj=numeric_limits<double>::quiet_NaN();
  return VALUE_OUT_OF_RANGE;
}

//------------------
// MatchToBCAL
//------------------
// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match. 
//
// NOTE:  an initial guess for tproj is expected as input so that out-of-time 
// hits can be skipped
jerror_t DParticleID::MatchToBCAL(const DReferenceTrajectory *rt,const vector<const DBCALShower*>& locInputBCALShowers, deque<const DBCALShower*>& locMatchedBCALShowers, double& locProjectedTime, double& locPathLength, double& locFlightTime) const{
  if (locInputBCALShowers.size() == 0){
    locProjectedTime = numeric_limits<double>::quiet_NaN();
    return RESOURCE_UNAVAILABLE;
  }
  // Clear the list of matched showers
  locMatchedBCALShowers.clear();

  // Minimum matching distance
  double d_min=1000.;
  
  //Loop over bcal showers
  double locTempPathLength, dphi, dz;
  for (unsigned int k = 0; k < locInputBCALShowers.size(); ++k){
    // Get the BCAL cluster position and normal
    const DBCALShower *shower = locInputBCALShowers[k];
    DVector3 bcal_pos(shower->x, shower->y, shower->z); 
    DVector3 proj_pos;
    // and the bcal cluster position, looking for the minimum
    double locTempFlightTime = 0.0;
    double d = rt->DistToRTwithTime(bcal_pos, &locTempPathLength, &locTempFlightTime);

    if(!isfinite(d))
      continue;
    // Check that the hit is not out of time with respect to the track
    if (fabs(locInputBCALShowers[k]->t - locTempFlightTime - locProjectedTime) > OUT_OF_TIME_CUT) continue;
    
    proj_pos = rt->GetLastDOCAPoint();
    
    dz = proj_pos.z() - bcal_pos.z();
    dphi = proj_pos.Phi() - bcal_pos.Phi();
    double p = rt->swim_steps[0].mom.Mag();
    dphi += 0.002 + 8.314e-3/(p + 0.3788)/(p + 0.3788);
    while(dphi >  M_PI) dphi -= M_TWO_PI;
    while(dphi < -M_PI) dphi += M_TWO_PI;
    double phi_sigma = 0.025 + 5.8e-4/(p*p*p);
    if (fabs(dz) < 10. && fabs(dphi) < 3.*phi_sigma){
      // Keep all possible matches, but keep the one with the smallest separation 
      // at the beginning of the deque.
      if (d<d_min){
	d_min=d;

	locPathLength = locTempPathLength;
	locFlightTime = locTempFlightTime;
	
	locMatchedBCALShowers.push_front(shower);
      }
      else {
	locMatchedBCALShowers.push_back(shower);
      }
    }
  }
  
  if(locMatchedBCALShowers.size() > 0){
    // Projected time at the target
    locProjectedTime = locMatchedBCALShowers[0]->t - locFlightTime;
    // locProjectedTime += 0.03; // correction determined from observation of simulated data //Paul Mattione
    
    return NOERROR;
  }
  
  locProjectedTime = numeric_limits<double>::quiet_NaN();
  return VALUE_OUT_OF_RANGE;
}

//------------------
// MatchToFCAL
//------------------
// Loop over fcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match. 
//
// NOTE:  an initial guess for locProjectedTime is expected as input so that out-of-time 
// hits can be skipped
jerror_t DParticleID::MatchToFCAL(const DReferenceTrajectory *rt, const vector<const DFCALShower*>& locInputFCALShowers, deque<const DFCALShower*>& locMatchedFCALShowers, double& locProjectedTime, double& locPathLength, double& locFlightTime,double *dEdx) const{
  if (locInputFCALShowers.size() == 0){
    locProjectedTime = numeric_limits<double>::quiet_NaN();
    return RESOURCE_UNAVAILABLE;
  }
  // Clear the list of matched showers
  locMatchedFCALShowers.clear();
  
  // Minimum matching distance
  double d_min=1000.;
  
  // Local flight path and flight time variables
  double locTempPathLength=0.,locTempFlightTime=0.;

  // loop over clusters
  for (unsigned int k = 0;k < locInputFCALShowers.size(); k++){
    
    // Get the FCAL cluster position and normal vector for the FCAL plane
    DVector3 fcal_pos = locInputFCALShowers[k]->getPosition();
    DVector3 norm(0,0,1);
    DVector3 proj_pos, proj_mom;
    
    // Find the distance of closest approach between the track trajectory
    // and the cluster position, looking for the minimum
    if (rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locTempPathLength, &locTempFlightTime)==NOERROR){
      double d = (fcal_pos - proj_pos).Mag();
      
      // Check that the hit is not out of time with respect to the track
      if (fabs(locInputFCALShowers[k]->getTime() - locTempFlightTime - locProjectedTime) > OUT_OF_TIME_CUT) continue;
      if(d < DELTA_R_FCAL){
	// Keep all possible matches, but keep the one with the smallest separation 
	// at the beginning of the deque.
	if (d<d_min){
	  d_min=d;
	  
	  locPathLength = locTempPathLength;
	  locFlightTime = locTempFlightTime;
  
	  if (dEdx!=NULL){
	    double p=proj_mom.Mag();
	    double dx=45.0*p/proj_mom.Dot(norm);
	    *dEdx=locInputFCALShowers[k]->getEnergy()/dx;
	  }
	  locMatchedFCALShowers.push_front(locInputFCALShowers[k]);
	}
	else{
	  locMatchedFCALShowers.push_back(locInputFCALShowers[k]);
	}
      }
    }
  }
  
  // Check for a match 
  if(locMatchedFCALShowers.size() > 0){
    locProjectedTime = locMatchedFCALShowers[0]->getTime() - locFlightTime;
    return NOERROR;
  }
  
  locProjectedTime = numeric_limits<double>::quiet_NaN();
  return VALUE_OUT_OF_RANGE;
}


//------------------
// MatchToSC
//------------------
// Match track to the start counter paddles with hits.  If a match
// is found, use the z-position of the track projection to the start counter 
// planes to correct for the light propagation within the scintillator and 
// estimate the "vertex" time.
//
// Unlike the next method, this method does not use the reference trajectory.
//
// NOTE:  an initial guess for tproj is expected as input so that out-of-time 
// hits can be skipped
jerror_t DParticleID::MatchToSC(const DKinematicData &parms, 
				vector<const DSCHit*>&sc_hits, 
				double &tproj,unsigned int &sc_match_id) const{ 
  sc_match_id=0;
  if (sc_hits.size()==0){
    tproj=numeric_limits<double>::quiet_NaN();
    return RESOURCE_UNAVAILABLE;
  }
  double myz=0.;
  double dphi_min=10000.,myphi=0.;
  unsigned int num=sc_norm.size()-1;

  DVector3 pos(parms.position());
  DVector3 mom(parms.momentum());
  stepper->SetCharge(parms.charge());

  // Swim to barrel representing the start counter straight portion
  double ds=0.,myds=0;
  if (stepper->SwimToRadius(pos,mom,sc_pos[1].x(),&ds)){
    tproj=numeric_limits<double>::quiet_NaN();
    return VALUE_OUT_OF_RANGE;
  }
  // Change the sign of ds -- most of the time we will need to add a small 
  // amount of time to the time calculated at the start counter
  ds*=-1.;

  // Position along z and phi at intersection
  double proj_z=pos.z();
  double proj_phi=pos.Phi();
  if (proj_phi<0) proj_phi+=M_TWO_PI;

  // Position of transition between start counter nose and leg
  double sc_pos1=sc_pos[1].z();
  // position in z at the end of the nose
  double sc_posn=sc_pos[num].z();

  // loop over sc hits, looking for the one with closest phi value to the 
  // projected phi value
  for (unsigned int i=0;i<sc_hits.size();i++){
    // Check that the hit is not out of time with respect to the track
    if (fabs(tproj-sc_hits[i]->t)>OUT_OF_TIME_CUT) continue;

    double phi=dSCphi0+dSCdphi*(sc_hits[i]->sector-1);
    double dphi=phi-proj_phi;
 
    // If the z position is in the nose region, match to the appropriate start
    // counter plane
    if (proj_z>sc_pos1){
      double cosphi=cos(phi);
      double sinphi=sin(phi);    
      for (unsigned int i=1;i<num;i++){
	DVector3 mymom=(-1.)*mom;
	DVector3 mypos=pos;
	double xhat=sc_norm[i].x(); 
	double r=sc_pos[i].X();
	double x=r*cosphi;
	double y=r*sinphi;
	double z=sc_pos[i].z();
	DVector3 norm(x*xhat,y*xhat,z);  
	DVector3 plane(x,y,z);
	double ds2=0.;
	if (stepper->SwimToPlane(mypos,mymom,plane,norm,&ds2)) continue;

	proj_z=mypos.z();
	if (proj_z<sc_pos[i+1].z()){
	  proj_phi=mypos.Phi();
	  if (proj_phi<0) proj_phi+=M_TWO_PI;
	  dphi=proj_phi-phi;
	  ds+=ds2;

	  break;
	}
      }
    }

    // Look for smallest difference in phi
    if (fabs(dphi)<dphi_min){
      dphi_min=fabs(dphi);
      myphi=phi;
      myz=proj_z;
      myds=ds;
      sc_match_id=i;
    }
  }
  
  // Look for a match in phi
  if (dphi_min<0.16){
    // Find the time at the start counter, before correcting for position
    // along the start counter
    tproj=sc_hits[sc_match_id]->t-sc_leg_tcor;
    
    // Check position along z
    if (myz<sc_pos[0].z()) myz=sc_pos[0].z();
    if (myz<sc_pos1){
      // Leg region
      tproj-=myz/C_EFFECTIVE;
    }
    else if (myz<sc_posn){
      // Nose region
      tproj-=((myz-sc_pos1)*sc_angle_cor+sc_pos1)/C_EFFECTIVE;
    }
    else{
      // Assume that the particle hit the most downstream z position of the
      // start counter
      tproj-=((sc_posn-sc_pos1)*sc_angle_cor+sc_pos1)/C_EFFECTIVE;
    }
	   
    // Adjust for flight time to the position pos
    double mass=parms.mass();
    double p2=mom.Mag2();
    double one_over_beta=sqrt(1.+mass*mass/p2);   
    tproj += myds*one_over_beta/SPEED_OF_LIGHT;

    return NOERROR;
  }

  tproj=numeric_limits<double>::quiet_NaN();
  return VALUE_OUT_OF_RANGE;

}


//------------------
// MatchToSC
//------------------
// Match track to the start counter paddles with hits.  If a match
// is found, use the z-position of the track projection to the start counter 
// planes to correct for the light propagation within the scintillator and 
// estimate the "vertex" time.
//
// NOTE:  an initial guess for tproj is expected as input so that out-of-time 
// hits can be skipped
jerror_t DParticleID::MatchToSC(const DReferenceTrajectory *rt, DTrackFitter::fit_type_t fit_type, vector<const DSCHit*>&sc_hits, double &tproj,unsigned int &sc_match_id, double &locPathLength, double &locFlightTime,pair<double,double>*dEdx) const{

  //tproj=numeric_limits<double>::quiet_NaN();
  sc_match_id=0;
  if (sc_hits.size()==0 || sc_pos.size()==0 || sc_norm.size()==0){
    tproj=numeric_limits<double>::quiet_NaN();
    return RESOURCE_UNAVAILABLE;
  }
  double myz=0.;
  double dphi_min=10000.,myphi=0.;
  double locTempPathLength, locTempFlightTime = 0.0;
  DVector3 proj_pos,proj_mom;
  
  // Find intersection with a "barrel" approximation for the start counter
  DVector3 dir; // direction of the momentum at the intersection point
  rt->GetIntersectionWithRadius(sc_pos[1].x(),proj_pos,&locTempPathLength,&locTempFlightTime,&proj_mom);
  double proj_phi=proj_pos.Phi();
  if (proj_phi<0) proj_phi+=M_TWO_PI;

  // loop over sc hits, looking for the one with closest phi value to the 
  // projected phi value
  for (unsigned int i=0;i<sc_hits.size();i++){
    // Check that the hit is not out of time with respect to the track
    if (fabs(tproj-sc_hits[i]->t)>OUT_OF_TIME_CUT) continue;

    double phi=dSCphi0+dSCdphi*(sc_hits[i]->sector-1);
    double dphi=phi-proj_phi;
    if (fabs(dphi)<dphi_min){
      dphi_min=fabs(dphi);
      myphi=phi;
      myz=proj_pos.z();
      sc_match_id=i;
    }
  }
  // Look for a match in phi
  if (dphi_min<0.21){
    // Length along scintillator
    double L=0.;

    // Initialize the normal vector for the SC paddle to the long, unbent region
    DVector3 norm(cos(myphi),sin(myphi),0.);

    // Now check to see if the intersection is in the nose region and find the
    // start time
    tproj=sc_hits[sc_match_id]->t-sc_leg_tcor;
    double sc_pos0=sc_pos[0].z();
    if (myz<sc_pos0) myz=sc_pos0;
    if (myz>sc_pos[1].z()){
      unsigned int num=sc_norm.size()-1;
      for (unsigned int i=1;i<num;i++){
	double xhat=sc_norm[i].x();
	norm.SetXYZ(cos(myphi)*xhat,sin(myphi)*xhat,sc_norm[i].z());
	double r=sc_pos[i].X();
	DVector3 pos(r*cos(myphi),r*sin(myphi),sc_pos[i].z());
	rt->GetIntersectionWithPlane(pos,norm,proj_pos,proj_mom,&locTempPathLength,&locTempFlightTime);
	myz=proj_pos.z();
	if (myz<sc_pos[i+1].z())
	  break;
      }
      double sc_pos1=sc_pos[1].z();
      // Note: in the following code, L does not include a correction for where the start counter starts
      // in z...  This is absorbed into tproj, above.
      if (myz<sc_pos1){
	L=sc_pos1;
	tproj-=locTempFlightTime+L/C_EFFECTIVE;
	locFlightTime = locTempFlightTime;
	locPathLength = locTempPathLength;
      }else if (myz>sc_pos[num].z()){
	// Assume that the particle hit the most downstream z position of the
	// start counter
	double costheta=rt->swim_steps[0].mom.CosTheta();
	double s=(sc_pos[num].z()-rt->swim_steps[0].origin.z())/costheta;
	double mass=rt->GetMass();
	double p2=rt->swim_steps[0].mom.Mag2();
	double one_over_beta=sqrt(1.+mass*mass/p2);
	L= (sc_pos[num].z()-sc_pos1)*sc_angle_cor+sc_pos1;
	
	tproj -= s*one_over_beta/SPEED_OF_LIGHT + L/C_EFFECTIVE;
	locFlightTime = s*one_over_beta/SPEED_OF_LIGHT;
	locPathLength = s;
      }else{
	L=(myz-sc_pos1)*sc_angle_cor+sc_pos1;
	tproj-=locTempFlightTime+L/C_EFFECTIVE;
	locFlightTime = locTempFlightTime;
	locPathLength = locTempPathLength;
      }
    }else{
      L=myz;
      tproj-=locTempFlightTime+L/C_EFFECTIVE;
      locFlightTime = locTempFlightTime;
      locPathLength = locTempPathLength;
    }
    if (dEdx!=NULL){ 
      double p=proj_mom.Mag();
      double dx=0.3*p/proj_mom.Dot(norm);
      // For the dEdx measurement we now need to take into account that L does not 
      // compensate for the position in z at which the start counter paddle starts
      double dE=(sc_hits[sc_match_id]->dE)*exp((L-sc_pos0)/ATTEN_LENGTH);	 
      double dE_sigma,dE_mp;
      GetScintMPdEandSigma(p,rt->GetMass(),dx,dE_mp,dE_sigma);
	  
      dEdx->first=dE/dx;
      dEdx->second=(dE-dE_mp)/dE_sigma;
    }
 
    return NOERROR;
  }
  
  tproj=numeric_limits<double>::quiet_NaN();
  return VALUE_OUT_OF_RANGE;
}

void DParticleID::GetScintMPdEandSigma(double p,double M,double x,
					 double &most_probable_dE,
					 double &sigma_dE) const{
  double one_over_beta_sq=1.+M*M/(p*p);
  double betagamma=p/M;
  double Xi=dKRhoZoverA_Scint*one_over_beta_sq*x;
  double Me=0.000511;

  // Density effect
  double X=log10(betagamma);
  double X0,X1;
  double Cbar=2.*(dLnI_Scint-log(28.816e-9*sqrt(dRhoZoverA_Scint)))+1.;
  if (dLnI_Scint<-1.6118){ // I<100
    if (Cbar<=3.681) X0=0.2;
    else X0=0.326*Cbar-1.;
    X1=2.;
  }
  else{
    if (Cbar<=5.215) X0=0.2;
    else X0=0.326*Cbar-1.5;
    X1=3.;
  }
  double delta=0;
  if (X>=X0 && X<X1)
    delta=4.606*X-Cbar+(Cbar-4.606*X0)*pow((X1-X)/(X1-X0),3.);
  else if (X>=X1)
    delta= 4.606*X-Cbar;  
  
  most_probable_dE=Xi*(log(Xi)-2.*dLnI_Scint-1./one_over_beta_sq
		       -log((one_over_beta_sq-1)/(2.*Me))+0.2-delta);
  sigma_dE=4.*Xi/2.354;
}

bool DParticleID::Calc_PropagatedRFTime(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, double& locPropagatedRFTime, bool locRFTimeFixedFlag) const
{
	locPropagatedRFTime = locEventRFBunch->dTime + (locChargedTrackHypothesis->z() - dTargetZCenter)/SPEED_OF_LIGHT;
	return true; //disabled until timing/etc. issues sorted out. (locEventRFBunch->dTime always = 0)

	if(locRFTimeFixedFlag)
		return locEventRFBunch->dMatchedToTracksFlag;

	//Propagates RF time to the track vertex-z, and then selects the closest RF bunch
	//Method: match track to RF bunch.  If cannot reliably match (e.g. no TOF or start counter hit) then use the best guess for this event (from locEventRFBunch) if available
		//First use TOF hit if any, then use a BCAL hit if track is fast enough (resolution low enough), else use start counter hit
	double locProjectedHitTime;

	//Use TOF hit if any
	bool locMatchFlag = false;
	double locP = locChargedTrackHypothesis->momentum().Mag();
	if(locChargedTrackHypothesis->t1_detector() == SYS_TOF)
	{
		locProjectedHitTime = locChargedTrackHypothesis->time();
		locMatchFlag = true;
	}
	else if((locChargedTrackHypothesis->t1_detector() == SYS_BCAL) && (locP > 0.25))
	{
		//if p < 250 MeV/c: too slow for the BCAL to distinguish RF bunch
			//at 225 MeV/c, BCAL t-resolution is ~333ps (3sigma = 999ps), BCAL delta-t error is ~40ps: ~1040ps: bad
			//at 250 MeV/c, BCAL t-resolution is ~310ps (3sigma = 920ps), BCAL delta-t error is ~40ps: ~960 ps < 1 ns: OK!!
		locProjectedHitTime = locChargedTrackHypothesis->time();
		locMatchFlag = true;
	}

	if(!locMatchFlag) // this track can't distinguish which RF bunch: use the propagated RF time from locEventRFBunch if it was succesfully matched to other tracks, else abort
		return locEventRFBunch->dMatchedToTracksFlag;
			
	//have a matching hit in either TOF/BCAL/ST: match this track to the closest RF bunch
	//If using an ST hit, this assumes that the start counter resolution will be good enough to unambiguously determine the RF bunch for each track individually.
		//Otherwise just delete the ST matching section above
	while((locPropagatedRFTime - locProjectedHitTime) > (0.5*dRFBunchFrequency))
		locPropagatedRFTime -= dRFBunchFrequency;
	while((locPropagatedRFTime - locProjectedHitTime) < (-0.5*dRFBunchFrequency))
		locPropagatedRFTime += dRFBunchFrequency;
	return true;
}

bool DParticleID::Calc_TrackStartTime(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, double& locStartTime, double& locStartTimeVariance, bool& locUsedRFTimeFlag, bool locRFTimeFixedFlag) const
{
	//use RF bunch if available, else use t0 detector if available (e.g. start counter, CDC)
	double locPropagatedRFTime;
	locUsedRFTimeFlag = Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch, locPropagatedRFTime, locRFTimeFixedFlag);
	if(locUsedRFTimeFlag)
	{
		locStartTime = locPropagatedRFTime;
		locStartTimeVariance = locEventRFBunch->dTimeVariance;
	}
	else //no confidence in selecting the RF bunch
	{
		locStartTime = locChargedTrackHypothesis->t0(); //BEWARE, THIS MAY CHANGE AS SC SYSTEM IS UPDATED!! (may need to propagate time to beamline here)
		locStartTimeVariance = locChargedTrackHypothesis->t0_err()*locChargedTrackHypothesis->t0_err();
		//try to use the start counter or CDC time instead of the RF time
		if(locChargedTrackHypothesis->t0_detector() == SYS_NULL)
			return false;
	}
	return true;
}

void DParticleID::Calc_TimingChiSq(DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag) const
{
	double locStartTime, locStartTimeVariance;
	bool locUsedRFTimeFlag;
	if(!Calc_TrackStartTime(locChargedTrackHypothesis, locEventRFBunch, locStartTime, locStartTimeVariance, locUsedRFTimeFlag, locRFTimeFixedFlag))
	{
		locChargedTrackHypothesis->dChiSq_Timing = 0.0;
		locChargedTrackHypothesis->dNDF_Timing = 0;
		return;
	}

	if((!locUsedRFTimeFlag) && (locChargedTrackHypothesis->t0_detector() == locChargedTrackHypothesis->t1_detector())) //e.g. both CDC (not matched to any hits)
	{
		locChargedTrackHypothesis->dChiSq_Timing = 0.0;
		locChargedTrackHypothesis->dNDF_Timing = 0;
		return;
	}

	double locTimeDifference = locStartTime - locChargedTrackHypothesis->time();
	double locTimeDifferenceVariance = (locChargedTrackHypothesis->errorMatrix())(6, 6) + locStartTimeVariance;

	double locTimingChiSq = locTimeDifference*locTimeDifference/locTimeDifferenceVariance;
	locChargedTrackHypothesis->dChiSq_Timing = locTimingChiSq;
	locChargedTrackHypothesis->dNDF_Timing = 1;
}

Particle_t DParticleID::IDTrack(float locCharge, float locMass) const
{
	float locMassTolerance = 0.010;
	if (locCharge > 0.1) // Positive particles
	{
		if (fabs(locMass - ParticleMass(Proton)) < locMassTolerance) return Proton;
		if (fabs(locMass - ParticleMass(PiPlus)) < locMassTolerance) return PiPlus;
		if (fabs(locMass - ParticleMass(KPlus)) < locMassTolerance) return KPlus;
		if (fabs(locMass - ParticleMass(Positron)) < locMassTolerance) return Positron;
		if (fabs(locMass - ParticleMass(MuonPlus)) < locMassTolerance) return MuonPlus;
	}
	else if(locCharge < -0.1) // Negative particles
	{
		if (fabs(locMass - ParticleMass(PiMinus)) < locMassTolerance) return PiMinus;
		if (fabs(locMass - ParticleMass(KMinus)) < locMassTolerance) return KMinus;
		if (fabs(locMass - ParticleMass(MuonMinus)) < locMassTolerance) return MuonMinus;
		if (fabs(locMass - ParticleMass(Electron)) < locMassTolerance) return Electron;
		if (fabs(locMass - ParticleMass(AntiProton)) < locMassTolerance) return AntiProton;
	}
	else //Neutral Track
	{
		if (fabs(locMass - ParticleMass(Gamma)) < locMassTolerance) return Gamma;
		if (fabs(locMass - ParticleMass(Neutron)) < locMassTolerance) return Neutron;
	}
	return Unknown;
}

void DParticleID::Calc_ChargedPIDFOM(DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag) const
{
	CalcDCdEdxChiSq(locChargedTrackHypothesis);
	Calc_TimingChiSq(locChargedTrackHypothesis, locEventRFBunch, locRFTimeFixedFlag);

	unsigned int locNDF_Total = locChargedTrackHypothesis->dNDF_Timing + locChargedTrackHypothesis->dNDF_DCdEdx;
	double locChiSq_Total = locChargedTrackHypothesis->dChiSq_Timing + locChargedTrackHypothesis->dChiSq_DCdEdx;

	/* Disable inclusion of TOF and SC dEdx.  dE/dx is Landau-distributed, not-Gaussian.
		Also, In the region (below 1 GeV/c) where some proton-pion separation is possible, the flight time is a much better tool.

	   if (locChargedTrackHypothesis->dTOFdEdx>0.){
	   locNDF_Total++;
	   locChiSq_Total+=locChargedTrackHypothesis->dTOFdEdx_norm_residual*locChargedTrackHypothesis->dTOFdEdx_norm_residual;
	   }

	if (locChargedTrackHypothesis->dStartCounterdEdx>0.){
	  locNDF_Total++;
	  locChiSq_Total+=locChargedTrackHypothesis->dStartCounterdEdx_norm_residual*locChargedTrackHypothesis->dStartCounterdEdx_norm_residual;
	}
	*/

	locChargedTrackHypothesis->dChiSq = locChiSq_Total;
	locChargedTrackHypothesis->dNDF = locNDF_Total;
	locChargedTrackHypothesis->dFOM = (locNDF_Total > 0) ? TMath::Prob(locChiSq_Total, locNDF_Total) : numeric_limits<double>::quiet_NaN();

}

