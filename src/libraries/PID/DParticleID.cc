// $Id$
//
//    File: DParticleID.cc
// Created: Mon Feb 28 14:48:56 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID.h"

#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

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

	C_EFFECTIVE = 15.0;
	ATTEN_LENGTH = 150.0;
	OUT_OF_TIME_CUT = 200.0;

  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  if(!dapp){
    _DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
  }
  const DRootGeom *RootGeom = dapp->GetRootGeom(loop->GetJEvent().GetRunNumber());

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
	DGeometry* locGeometry = dapp->GetDGeometry(loop->GetJEvent().GetRunNumber());

	// Check if Start Counter geometry is present
	vector<double> sc_origin;
	bool got_sc = locGeometry->Get("//posXYZ[@volume='StartCntr']/@X_Y_Z", sc_origin);
	if(got_sc)
	{
		double num_paddles;
		locGeometry->Get("//mposPhi[@volume='STRC']/@ncopy",num_paddles); 
		dSCdphi = M_TWO_PI/num_paddles;
		
		double Phi0;
		locGeometry->Get("///mposPhi[@volume='STRC']/@Phi0",Phi0);
		dSCphi0 = Phi0*M_PI/180.;

		vector<vector<double> > sc_rioz;
		locGeometry->GetMultiple("//pgon[@name='STRC']/polyplane/@Rio_Z", sc_rioz);

		for(unsigned int k = 0; k < sc_rioz.size() - 1; ++k)
		{
			if(sc_origin.size() < 3)
				continue; // in case start counter is comment out in XML
			if(sc_rioz[k].size() < 3)
				continue; // in case start counter is comment out in XML
			DVector3 pos((sc_rioz[k][0]+sc_rioz[k][1])/2.,0.,sc_rioz[k][2]+sc_origin[2]);
			DVector3 dir(sc_rioz[k+1][2]-sc_rioz[k][2],0,-sc_rioz[k+1][0]+sc_rioz[k][0]);
			dir.SetMag(1.);

			sc_pos.push_back(pos);
			sc_norm.push_back(dir);		
		}
	
		// in case start counter is comment out in XML
		if(!sc_pos.empty())
		{
			// sc_leg_tcor=(sc_light_guide[2]-sc_pos[0].z())/C_EFFECTIVE;
			sc_leg_tcor = -sc_pos[0].z()/C_EFFECTIVE;
			double theta = sc_norm[sc_norm.size() - 1].Theta();
			sc_angle_cor = 1./cos(M_PI - theta);
		}
	}

	//Get calibration constants
	map<string, double> locPIDParams;
	if(!loop->GetCalib("PID/photon_track_matching", locPIDParams))
	{
		static bool printed_message = false;
		if(!printed_message){
			printed_message = true;
			cout<<"DParticleID: loading values from PID data base"<<endl;
		}
		DELTA_R_FCAL = locPIDParams["DELTA_R_FCAL"];
	}
	else
	{
		cout << "DParticleID: Error loading values from PID data base" <<endl;
		DELTA_R_FCAL = 15.0;
	}
	FCAL_CUT_PAR1=3.3;
	gPARMS->SetDefaultParameter("FCAL:CUT_PAR1",FCAL_CUT_PAR1);

	FCAL_CUT_PAR2=0.88;
	gPARMS->SetDefaultParameter("FCAL:CUT_PAR2",FCAL_CUT_PAR2);

	BCAL_Z_CUT=10.;
	gPARMS->SetDefaultParameter("BCAL:Z_CUT",BCAL_Z_CUT);

	BCAL_PHI_CUT_PAR1=0.0175;
	gPARMS->SetDefaultParameter("BCAL:PHI_CUT_PAR1",BCAL_PHI_CUT_PAR1);

	BCAL_PHI_CUT_PAR2=1.76e-3;
	gPARMS->SetDefaultParameter("BCAL:PHI_CUT_PAR2",BCAL_PHI_CUT_PAR2);

	dTargetZCenter = 0.0;
	locGeometry->GetTargetZ(dTargetZCenter);
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

//------------------
// MatchToBCAL
//------------------
// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match. 
//
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// hits can be skipped

//to be called by track reconstruction
bool DParticleID::MatchToBCAL(const DReferenceTrajectory* rt, const vector<const DBCALShower*>& locBCALShowers, double& locStartTime, double& locTimeVariance) const
{
	double d_min = 9.9E9;

	//Loop over bcal showers
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		DShowerMatchParams locShowerMatchParams;
		if(!MatchToBCAL(NULL, rt, locBCALShowers[loc_i], locStartTime, locShowerMatchParams))
			continue;
		if(locShowerMatchParams.dDOCAToShower >= d_min)
			continue;
		d_min = locShowerMatchParams.dDOCAToShower;
		locStartTime = locBCALShowers[loc_i]->t - locShowerMatchParams.dFlightTime;
//		locTimeVariance = sqrt(locShowerMatchParams.dFlightTimeVariance) - locBCALShowers[loc_i]->dCovarianceMatrix(4, 4); //uncomment when ready!!
//		locTimeVariance *= locTimeVariance; //uncomment when ready!!
		locTimeVariance = 0.5*0.5;
		locMatchFlag = true;
	}

	return locMatchFlag;
}

bool DParticleID::MatchToBCAL(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DBCALShower* locBCALShower, double locInputStartTime, DShowerMatchParams& locShowerMatchParams) const
{
	// NOTE: locTrackTimeBased is NULL if calling from track reconstruction!!!
	// Get the BCAL cluster position and normal
	DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z); 

	double locFlightTime = 9.9E9, locPathLength = 9.9E9;
	double d = rt->DistToRTwithTime(bcal_pos, &locPathLength, 
					&locFlightTime,SYS_BCAL);

	if(!isfinite(d))
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locBCALShower->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	DVector3 proj_pos = rt->GetLastDOCAPoint();
	if (proj_pos.Perp()<65.) return false;  // not inside BCAL!
	
	double dz = proj_pos.z() - bcal_pos.z();
	double dphi = proj_pos.Phi() - bcal_pos.Phi();
	double p = rt->swim_steps[0].mom.Mag();
	//dphi += 0.002 + 8.314e-3/(p + 0.3788)/(p + 0.3788);
	while(dphi >	M_PI)
		dphi -= M_TWO_PI;
	while(dphi < -M_PI)
		dphi += M_TWO_PI;
	double phi_cut = BCAL_PHI_CUT_PAR1 + BCAL_PHI_CUT_PAR2/(p*p);

	if((fabs(dz) >= BCAL_Z_CUT) || (fabs(dphi) >= phi_cut))
		return false; //not close enough

	//successful match
	locShowerMatchParams.dTrackTimeBased = locTrackTimeBased;
	locShowerMatchParams.dShowerObject = locBCALShower;
	locShowerMatchParams.dx = 0.0; //SET ME!!!!
	locShowerMatchParams.dFlightTime = locFlightTime;
	locShowerMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!!
	locShowerMatchParams.dPathLength = locPathLength;
	locShowerMatchParams.dDOCAToShower = d; //DOCA of track to shower

	return true;
}

//------------------
// MatchToTOF
//------------------
// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match. 
//
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// hits can be skipped

//to be called by track reconstruction
bool DParticleID::MatchToTOF(const DReferenceTrajectory* rt, const vector<const DTOFPoint*>& locTOFPoints, double& locStartTime, double& locTimeVariance) const
{
	double d_min = 9.9E9;

	//Loop over tof points
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
	{
		DTOFHitMatchParams locTOFHitMatchParams;
		if(!MatchToTOF(NULL, rt, locTOFPoints[loc_i], locStartTime, locTOFHitMatchParams))
			continue;
		if(locTOFHitMatchParams.dDOCAToHit >= d_min)
			continue;

		d_min = locTOFHitMatchParams.dDOCAToHit;
		locStartTime = locTOFPoints[loc_i]->t - locTOFHitMatchParams.dFlightTime;
//		locTimeVariance = sqrt(locTOFHitMatchParams.dFlightTimeVariance) - locTOFPoints[loc_i]->tErr; //uncomment when ready!
//		locTimeVariance *= locTimeVariance; //uncomment when ready!
		locTimeVariance = 0.1*0.1;
		locMatchFlag = true;
	}

	return locMatchFlag;
}

bool DParticleID::MatchToTOF(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DTOFPoint* locTOFPoint, double locInputStartTime, DTOFHitMatchParams& locTOFHitMatchParams) const
{
	// NOTE: locTrackTimeBased is NULL if calling from track reconstruction!!!

	// Find the distance of closest approach between the track trajectory
	// and the tof cluster position, looking for the minimum
	DVector3 tof_pos = locTOFPoint->pos;
	DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,SYS_TOF) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locTOFPoint->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	double d = (tof_pos - proj_pos).Mag();
	double match_cut = 6.15; //current dPositionMatchCut_DoubleEnded variable in DTOFPoint_factory.cc
	if(d >= match_cut)
		return false;

	//successful match
	double dx = 2.54*proj_mom.Mag()/proj_mom.Dot(norm);
	locTOFHitMatchParams.dTrackTimeBased = locTrackTimeBased;
	locTOFHitMatchParams.dTOFPoint = locTOFPoint;
	locTOFHitMatchParams.dEdx = (locTOFPoint->dE)/dx;
	locTOFHitMatchParams.dFlightTime = locFlightTime;
	locTOFHitMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!
	locTOFHitMatchParams.dPathLength = locPathLength;
	locTOFHitMatchParams.dDOCAToHit = d; //DOCA of track to hit

	return true;
}

//------------------
// MatchToFCAL
//------------------
// Loop over fcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match. 
//
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// hits can be skipped

//to be called by track reconstruction
bool DParticleID::MatchToFCAL(const DReferenceTrajectory* rt, const vector<const DFCALShower*>& locFCALShowers, double& locStartTime, double& locTimeVariance) const
{
	double d_min = 9.9E9;

	//Loop over FCAL showers
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		DShowerMatchParams locShowerMatchParams;
		if(!MatchToFCAL(NULL, rt, locFCALShowers[loc_i], locStartTime, locShowerMatchParams))
			continue;
		if(locShowerMatchParams.dDOCAToShower >= d_min)
			continue;

		d_min = locShowerMatchParams.dDOCAToShower;
		locStartTime = locFCALShowers[loc_i]->getTime() - locShowerMatchParams.dFlightTime;
//		locTimeVariance = sqrt(locShowerMatchParams.dFlightTimeVariance) - sqrt(locFCALShowers[loc_i]->dCovarianceMatrix(4, 4)); //uncomment when ready!
//		locTimeVariance *= locTimeVariance; //uncomment when ready!
		locTimeVariance = 0.5*0.5;
		locMatchFlag = true;
	}

	return locMatchFlag;
}

bool DParticleID::MatchToFCAL(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DFCALShower* locFCALShower, double locInputStartTime, DShowerMatchParams& locShowerMatchParams) const
{
	// NOTE: locTrackTimeBased is NULL if calling from track reconstruction!!!

	// Find the distance of closest approach between the track trajectory
	// and the cluster position, looking for the minimum
	DVector3 fcal_pos = locFCALShower->getPosition();
	DVector3 norm(0.0, 0.0, 1.0); //normal vector for the FCAL plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,SYS_FCAL) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locFCALShower->getTime() - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	double d = (fcal_pos - proj_pos).Mag();
	double p=proj_mom.Mag();
	double cut=FCAL_CUT_PAR1+FCAL_CUT_PAR2/p;
	if(d >= cut)
		return false;

	locShowerMatchParams.dTrackTimeBased = locTrackTimeBased;
	locShowerMatchParams.dShowerObject = locFCALShower;
	locShowerMatchParams.dx = 45.0*p/(proj_mom.Dot(norm));
	locShowerMatchParams.dFlightTime = locFlightTime;
	locShowerMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!
	locShowerMatchParams.dPathLength = locPathLength;
	locShowerMatchParams.dDOCAToShower = d;

	return true;
}

//------------------
// MatchToSC
//------------------
// Match track to the start counter paddles with hits.	If a match
// is found, use the z-position of the track projection to the start counter 
// planes to correct for the light propagation within the scintillator and 
// estimate the "vertex" time.
//
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// hits can be skipped
//to be called by track reconstruction
bool DParticleID::MatchToSC(const DReferenceTrajectory* rt, const vector<const DSCHit*>& locSCHits, double& locStartTime, double& locTimeVariance) const
{
	double dphi_min = 10000.0;

	//Loop over SC points
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	{
		DSCHitMatchParams locSCHitMatchParams;
		if(!MatchToSC(NULL, rt, locSCHits[loc_i], locStartTime, locSCHitMatchParams))
			continue;
		if(locSCHitMatchParams.dDeltaPhiToHit >= dphi_min)
			continue;

		dphi_min = locSCHitMatchParams.dDeltaPhiToHit;
		locStartTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime;
//		locTimeVariance = sqrt(locSCHitMatchParams.dFlightTimeVariance) - sqrt(locSCHitMatchParams.dHitTimeVariance); //uncomment when ready!
//		locTimeVariance *= locTimeVariance; //uncomment when ready!
		locTimeVariance = 0.3*0.3;
		locMatchFlag = true;
	}

	return locMatchFlag;
}

bool DParticleID::MatchToSC(const DTrackTimeBased* locTrackTimeBased, const DReferenceTrajectory* rt, const DSCHit* locSCHit, double locInputStartTime, DSCHitMatchParams& locSCHitMatchParams) const
{
	// NOTE: locTrackTimeBased is NULL if calling from track reconstruction!!!
	if(sc_pos.empty() || sc_norm.empty())
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locInputStartTime - locSCHit->t) > OUT_OF_TIME_CUT)
		return false;

	// Find intersection with a "barrel" approximation for the start counter
	DVector3 proj_pos(NaN,NaN,NaN), proj_mom(NaN,NaN,NaN);
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithRadius(sc_pos[1].x(), proj_pos, &locPathLength, &locFlightTime, &proj_mom) != NOERROR)
		return false;
	double proj_phi = proj_pos.Phi();
	if(proj_phi < 0.0)
		proj_phi += M_TWO_PI;

	// Look for a match in phi
	double phi = dSCphi0 + dSCdphi*(locSCHit->sector - 1);
	double dphi = phi - proj_phi;
	if(fabs(dphi) >= 0.21)
		return false; //no match

	//match successful
	double myphi = phi;
	double myz = proj_pos.z();

	// Length along scintillator
	double L = 0.;

	// Initialize the normal vector for the SC paddle to the long, unbent region
	DVector3 norm(cos(myphi), sin(myphi), 0.);

	// Now check to see if the intersection is in the nose region and find the
	// start time
	double locCorrectedHitTime = locSCHit->t - sc_leg_tcor;
	double sc_pos0 = sc_pos[0].z();
	if(myz < sc_pos0)
		myz = sc_pos0;

	if(myz <= sc_pos[1].z())
	{
		L=myz;
		locCorrectedHitTime -= L/C_EFFECTIVE;
	}
	else
	{
		unsigned int num = sc_norm.size() - 1;
		for (unsigned int loc_i = 1; loc_i < num; ++loc_i)
		{
			double xhat = sc_norm[loc_i].x();
			norm.SetXYZ(cos(myphi)*xhat, sin(myphi)*xhat, sc_norm[loc_i].z());
			double r = sc_pos[loc_i].X();
			DVector3 pos(r*cos(myphi), r*sin(myphi), sc_pos[loc_i].z());
			locPathLength = 9.9E9;
			locFlightTime = 9.9E9;
			if(rt->GetIntersectionWithPlane(pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime) != NOERROR)
           continue;
			myz = proj_pos.z();
			if(myz < sc_pos[loc_i + 1].z())
				break;
		}
		double sc_pos1 = sc_pos[1].z();

		// Note: in the following code, L does not include a correction for where the start counter starts
		// in z...	This is absorbed into locCorrectedHitTime, above.
		if(myz < sc_pos1)
		{
			L = sc_pos1;
			locCorrectedHitTime -= L/C_EFFECTIVE;
		}
		else if (myz > sc_pos[num].z())
		{
			// Assume that the particle hit the most downstream z position of the
			// start counter
			double costheta = rt->swim_steps[0].mom.CosTheta();
			double s = (sc_pos[num].z() - rt->swim_steps[0].origin.z())/costheta;
			double mass = rt->GetMass();
			double p2 = rt->swim_steps[0].mom.Mag2();
			double one_over_beta = sqrt(1. + mass*mass/p2);
			L = (sc_pos[num].z() - sc_pos1)*sc_angle_cor + sc_pos1;

			locCorrectedHitTime -= s*one_over_beta/SPEED_OF_LIGHT + L/C_EFFECTIVE;
			locFlightTime = s*one_over_beta/SPEED_OF_LIGHT;
			locPathLength = s;
		}
		else
		{
			L = (myz - sc_pos1)*sc_angle_cor + sc_pos1;
			locCorrectedHitTime -= L/C_EFFECTIVE;
		}
	}

	double dx = 0.3*proj_mom.Mag()/(proj_mom.Dot(norm));

	// For the dEdx measurement we now need to take into account that L does not 
	// compensate for the position in z at which the start counter paddle starts
	locSCHitMatchParams.dTrackTimeBased = locTrackTimeBased;
	locSCHitMatchParams.dSCHit = locSCHit;
	locSCHitMatchParams.dHitEnergy = (locSCHit->dE)*exp((L - sc_pos0)/ATTEN_LENGTH);
	locSCHitMatchParams.dEdx = locSCHitMatchParams.dHitEnergy/dx;
	locSCHitMatchParams.dHitTime = locCorrectedHitTime;
	locSCHitMatchParams.dHitTimeVariance = 0.0; //SET ME!!!
	locSCHitMatchParams.dFlightTime = locFlightTime;
	locSCHitMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!
	locSCHitMatchParams.dPathLength = locPathLength;
	locSCHitMatchParams.dDeltaPhiToHit = dphi;

	return true;
}

//------------------
// MatchToTrack
//------------------
// Loop over time-based tracks, looking for the closest one to the given shower
//
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// tracks can be skipped
bool DParticleID::MatchToTrack(const DBCALShower* locBCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const
{
	// Get the BCAL cluster position and normal
	DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z); 

	double locFlightTime = 9.9E9, locPathLength = 9.9E9;
	locDistance = rt->DistToRTwithTime(bcal_pos, &locPathLength, 
					   &locFlightTime,SYS_BCAL);
	if(!isfinite(locDistance))
		return false;
	// Check that the hit is not out of time with respect to the track
	if(fabs(locBCALShower->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;
	return true;
}

bool DParticleID::MatchToTrack(const DFCALShower* locFCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const
{
	DVector3 fcal_pos = locFCALShower->getPosition();
	DVector3 norm(0.0, 0.0, 1.0); //normal vector for the FCAL plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,SYS_FCAL) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locFCALShower->getTime() - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	locDistance = (fcal_pos - proj_pos).Mag();
	return true;
}

bool DParticleID::Get_BestSCMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DSCHitMatchParams& locBestMatchParams) const
{
	//choose the "best" detector hit to use for computing quantities
	vector<DSCHitMatchParams> locSCHitMatchParams;
	if(!locDetectorMatches->Get_SCMatchParams(locTrackTimeBased, locSCHitMatchParams))
		return false;

	double locMinDeltaPhi = 9.9E9;
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locSCHitMatchParams.size(); ++loc_i)
	{
		if(locSCHitMatchParams[loc_i].dDeltaPhiToHit >= locMinDeltaPhi)
			continue;
		locMinDeltaPhi = locSCHitMatchParams[loc_i].dDeltaPhiToHit;
		locBestMatchParams = locSCHitMatchParams[loc_i];
		locMatchFlag = true;
	}
	return locMatchFlag;
}

bool DParticleID::Get_BestBCALMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DShowerMatchParams& locBestMatchParams) const
{
	//choose the "best" shower to use for computing quantities
	vector<DShowerMatchParams> locShowerMatchParams;
	if(!locDetectorMatches->Get_BCALMatchParams(locTrackTimeBased, locShowerMatchParams))
		return false;

	double locMinDistance = 9.9E9;
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
	{
		if(locShowerMatchParams[loc_i].dDOCAToShower >= locMinDistance)
			continue;
		locMinDistance = locShowerMatchParams[loc_i].dDOCAToShower;
		locBestMatchParams = locShowerMatchParams[loc_i];
		locMatchFlag = true;
	}
	return locMatchFlag;
}

bool DParticleID::Get_BestTOFMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DTOFHitMatchParams& locBestMatchParams) const
{
	//choose the "best" hit to use for computing quantities
	vector<DTOFHitMatchParams> locTOFHitMatchParams;
	if(!locDetectorMatches->Get_TOFMatchParams(locTrackTimeBased, locTOFHitMatchParams))
		return false;

	double locMinDistance = 9.9E9;
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locTOFHitMatchParams.size(); ++loc_i)
	{
		if(locTOFHitMatchParams[loc_i].dDOCAToHit >= locMinDistance)
			continue;
		locMinDistance = locTOFHitMatchParams[loc_i].dDOCAToHit;
		locBestMatchParams = locTOFHitMatchParams[loc_i];
		locMatchFlag = true;
	}
	return locMatchFlag;
}

bool DParticleID::Get_BestFCALMatchParams(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches, DShowerMatchParams& locBestMatchParams) const
{
	//choose the "best" shower to use for computing quantities
	vector<DShowerMatchParams> locShowerMatchParams;
	if(!locDetectorMatches->Get_FCALMatchParams(locTrackTimeBased, locShowerMatchParams))
		return false;

	double locMinDistance = 9.9E9;
	bool locMatchFlag = false;
	for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
	{
		if(locShowerMatchParams[loc_i].dDOCAToShower >= locMinDistance)
			continue;
		locMinDistance = locShowerMatchParams[loc_i].dDOCAToShower;
		locBestMatchParams = locShowerMatchParams[loc_i];
		locMatchFlag = true;
	}
	return locMatchFlag;
}

double DParticleID::Calc_BCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const
{
	DShowerMatchParams locBCALShowerMatchParams;
	if(!Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
		return numeric_limits<double>::quiet_NaN();
//	const DBCALShower* locBCALShower = static_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_FCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const
{
	DShowerMatchParams locFCALShowerMatchParams;
	if(!Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
		return numeric_limits<double>::quiet_NaN();
//	const DFCALShower* locFCALShower = static_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_TOFFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const
{
	DTOFHitMatchParams locTOFHitMatchParams;
	if(!Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
		return numeric_limits<double>::quiet_NaN();
	//locTOFHitMatchParams->dTOFPoint;
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_SCFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, const DDetectorMatches* locDetectorMatches) const
{
	DSCHitMatchParams locSCHitMatchParams;
	if(!Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
		return numeric_limits<double>::quiet_NaN();
	//locSCHitMatchParams->dSCHit;
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
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
	//use RF bunch if available
	double locPropagatedRFTime = 0.0;
	locUsedRFTimeFlag = Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch, locPropagatedRFTime, locRFTimeFixedFlag);
	if(locUsedRFTimeFlag)
	{
		locStartTime = locPropagatedRFTime;
		locStartTimeVariance = locEventRFBunch->dTimeVariance;
	}
	else //no confidence in selecting the RF bunch: use t0 detector if available (e.g. start counter, CDC)
	{
		locStartTime = locChargedTrackHypothesis->t0();
		locStartTimeVariance = locChargedTrackHypothesis->t0_err()*locChargedTrackHypothesis->t0_err();
	}
	return true;
}

double DParticleID::Calc_TimingChiSq(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch, bool locRFTimeFixedFlag, unsigned int &locNDF, double& locPull) const
{
	double locStartTime = 0.0, locStartTimeVariance = 0.0;
	bool locUsedRFTimeFlag = false;
	if(!Calc_TrackStartTime(locChargedTrackHypothesis, locEventRFBunch, locStartTime, locStartTimeVariance, locUsedRFTimeFlag, locRFTimeFixedFlag))
	{
		locNDF = 0;
		locPull = 0.0;
		return 0.0;
	}

	if((!locUsedRFTimeFlag) && (locChargedTrackHypothesis->t0_detector() == locChargedTrackHypothesis->t1_detector())) //e.g. both CDC (not matched to any hits)
	{
		locNDF = 0;
		locPull = 0.0;
		return 0.0;
	}

	double locTimeDifference = locStartTime - locChargedTrackHypothesis->time();
	double locTimeDifferenceVariance = (locChargedTrackHypothesis->errorMatrix())(6, 6) + locStartTimeVariance;

	locNDF = 1;
	locPull = locTimeDifference/sqrt(locTimeDifferenceVariance);
	return locPull*locPull;
}

double DParticleID::Calc_TimingChiSq(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DEventRFBunch* locEventRFBunch, unsigned int &locNDF, double& locPull) const
{
	double locStartTimeVariance = locEventRFBunch->dTimeVariance;
	double locStartTime = locEventRFBunch->dMatchedToTracksFlag ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();

	double locTimeDifferenceVariance = (locNeutralParticleHypothesis->errorMatrix())(6, 6) + locStartTimeVariance;
	locPull = (locStartTime - locNeutralParticleHypothesis->time())/sqrt(locTimeDifferenceVariance);
	locNDF = 1;
	return locPull*locPull;
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

	unsigned int locTimingNDF = 0;
	double locTimingPull = 0.0;
	double locTimingChiSq = Calc_TimingChiSq(locChargedTrackHypothesis, locEventRFBunch, locRFTimeFixedFlag, locTimingNDF, locTimingPull);
	locChargedTrackHypothesis->dChiSq_Timing = locTimingChiSq;
	locChargedTrackHypothesis->dNDF_Timing = locTimingNDF;

	unsigned int locNDF_Total = locChargedTrackHypothesis->dNDF_Timing + locChargedTrackHypothesis->dNDF_DCdEdx;
	double locChiSq_Total = locChargedTrackHypothesis->dChiSq_Timing + locChargedTrackHypothesis->dChiSq_DCdEdx;

	locChargedTrackHypothesis->dChiSq = locChiSq_Total;
	locChargedTrackHypothesis->dNDF = locNDF_Total;
	locChargedTrackHypothesis->dFOM = (locNDF_Total > 0) ? TMath::Prob(locChiSq_Total, locNDF_Total) : numeric_limits<double>::quiet_NaN();
}

