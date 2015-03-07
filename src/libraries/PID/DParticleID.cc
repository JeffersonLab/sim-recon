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

// Routine for sorting hypotheses accorpding to FOM
bool static DParticleID_hypothesis_cmp(const DTrackTimeBased *a,
				       const DTrackTimeBased *b){
  return (a->FOM>b->FOM);
}


//---------------------------------
// DParticleID    (Constructor)
//---------------------------------
DParticleID::DParticleID(JEventLoop *loop)
{
	vector<double> locRFFrequencyVector;
	loop->GetCalib("PHOTON_BEAM/rf_frequency", locRFFrequencyVector);
	dRFBunchFrequency = locRFFrequencyVector[0];

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
  if(got_sc){
    // z-position at upstream face of scintillators.
    double z0=sc_origin[2];

    // Get rotation angles
    vector<double>sc_rot_angles;
    locGeometry->Get("//posXYZ[@volume='StartCntr']/@rot", sc_rot_angles);
    double ThetaX=sc_rot_angles[0]*M_PI/180.;
    double ThetaY=sc_rot_angles[1]*M_PI/180.;  
    double ThetaZ=sc_rot_angles[2]*M_PI/180.;
    //ThetaX=0.;
    //ThetaY=0.;

    double num_paddles;
    locGeometry->Get("//mposPhi[@volume='STRC']/@ncopy",num_paddles); 
    dSCdphi = M_TWO_PI/num_paddles;
		
    double Phi0;
    locGeometry->Get("///mposPhi[@volume='STRC']/@Phi0",Phi0);
    dSCphi0 =Phi0*M_PI/180.;
    
    vector<vector<double> > sc_rioz;
    locGeometry->GetMultiple("//pgon[@name='STRC']/polyplane/@Rio_Z", sc_rioz);

    // Create vectors of positions and normal vectors for each paddle
    for (unsigned int i=0;i<30;i++){
      double phi=dSCphi0+dSCdphi*double(i);
      double sinphi=sin(phi);
      double cosphi=cos(phi);
      double r=0.5*(sc_rioz[0][0]+sc_rioz[0][1]);
      DVector3 oldray;
      // Rotate by phi and take into account the tilt
      DVector3 ray(r*cosphi,r*sinphi,sc_rioz[0][2]);
      ray.RotateX(ThetaX);
      ray.RotateY(ThetaY);
      ray.RotateZ(ThetaZ);

      // Create stl-vectors to store positions and norm vectors
      vector<DVector3>posvec;
      vector<DVector3>dirvec;
      // Loop over radial/z positions describing start counter geometry from xml
      for(unsigned int k = 1; k < sc_rioz.size(); ++k){
	oldray=ray;
	r=0.5*(sc_rioz[k][0]+sc_rioz[k][1]);
	// Point in midplane of scintillator
	ray.SetXYZ(r*cosphi,r*sinphi,sc_rioz[k][2]);
	// Second point in the plane of the scintillator
	DVector3 ray2(r*cosphi-10.0*sinphi,r*sinphi+10.0*cosphi,sc_rioz[k][2]);
	// Take into account tilt
	ray.RotateX(ThetaX);
	ray.RotateY(ThetaY);
	ray.RotateZ(ThetaZ);
	ray2.RotateX(ThetaX);
	ray2.RotateY(ThetaY);
	ray2.RotateZ(ThetaZ);
	// Store one position on current plane
	posvec.push_back(DVector3(oldray.X(),oldray.Y(),oldray.Z()+z0));
	// Compute normal vector to plane
	DVector3 dir=(ray-oldray).Cross(ray2-oldray);
	dir.SetMag(1.);
	dirvec.push_back(dir);		
      }
      sc_pos.push_back(posvec);
      sc_norm.push_back(dirvec);
		  
      posvec.clear();
      dirvec.clear();
    }

	
		// in case start counter is comment out in XML
		if(!sc_pos.empty())
		{
			// sc_leg_tcor=(sc_light_guide[2]-sc_pos[0].z())/C_EFFECTIVE;
			sc_leg_tcor = -sc_pos[0][0].z()/C_EFFECTIVE;
			double theta = sc_norm[0][sc_norm[0].size()-2].Theta();
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

	SC_DPHI_CUT=0.105;
	gPARMS->SetDefaultParameter("SC:DPHI_CUT",SC_DPHI_CUT);
	
	SC_DPHI_CUT_WB=0.21;
	gPARMS->SetDefaultParameter("SC:DPHI_CUT_WB",SC_DPHI_CUT_WB);

	dTargetZCenter = 0.0;
	locGeometry->GetTargetZ(dTargetZCenter);

	
  // Track finder helper class
  vector<const DTrackFinder *> finders;
  loop->Get(finders);

  if(finders.size()<1){
    _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
    return;
  }

  finder = finders[0];

	//TOF calibration constants & geometry
	if(loop->GetCalib("TOF/propagation_speed", propagation_speed))
		jout << "Error loading /TOF/propagation_speed !" << endl;

	map<string, double> tofparms;
 	loop->GetCalib("TOF/tof_parms", tofparms);
	TOF_ATTEN_LENGTH = tofparms["TOF_ATTEN_LENGTH"];

	loop->GetSingle(dTOFGeometry);
	double locHalfPaddle_OneSided = dTOFGeometry->SHORTBARLENGTH/2.0; //GET FROM GEOMETRY??
	double locBeamHoleWidth = dTOFGeometry->LONGBARLENGTH - 2.0*dTOFGeometry->SHORTBARLENGTH;
	ONESIDED_PADDLE_MIDPOINT_MAG = locHalfPaddle_OneSided + locBeamHoleWidth/2.0;

	// Start counter calibration constants
	vector<map<string,double> >tvals;

	if(loop->GetCalib("/START_COUNTER/propagation_speed",tvals))
	  jout << "Error loading /START_COUNTER/propagation_speed !" << endl;
	else{
	  for(unsigned int i=0; i<tvals.size(); i++){
            map<string, double> &row = tvals[i];
	    sc_veff[SC_STRAIGHT].push_back(row["SC_STRAIGHT_PROPAGATION_B"]);
	    sc_veff[SC_BEND].push_back(row["SC_BEND_PROPAGATION_B"]);
	    sc_veff[SC_NOSE].push_back(row["SC_NOSE_PROPAGATION_B"]);
	  }
	}

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

jerror_t DParticleID::CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const
{
	vector<dedx_t> locdEdxHits_CDC, locdEdxHits_FDC;
	jerror_t locReturnStatus = GetDCdEdxHits(locTrackTimeBased, locdEdxHits_CDC, locdEdxHits_FDC);
	if(locReturnStatus != NOERROR)
	{
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

jerror_t DParticleID::CalcDCdEdx(const DTrackTimeBased *locTrackTimeBased, const vector<dedx_t>& locdEdxHits_CDC, const vector<dedx_t>& locdEdxHits_FDC, double& locdEdx_FDC, double& locdx_FDC, double& locdEdx_CDC, double& locdx_CDC, unsigned int& locNumHitsUsedFordEdx_FDC, unsigned int& locNumHitsUsedFordEdx_CDC) const
	{
	locdx_CDC = 0.0;
	locdEdx_CDC = 0.0;
	locNumHitsUsedFordEdx_CDC = locdEdxHits_CDC.size()/2;
	if(locNumHitsUsedFordEdx_CDC > 0)
	{
		for(unsigned int loc_i = 0; loc_i < locNumHitsUsedFordEdx_CDC; ++loc_i)
		{
			locdEdx_CDC += locdEdxHits_CDC[loc_i].dE; //weight is ~ #e- (scattering sites): dx!
			locdx_CDC += locdEdxHits_CDC[loc_i].dx;
		}
		locdEdx_CDC /= locdx_CDC;
	}

	locdx_FDC = 0.0;
	locdEdx_FDC = 0.0;
	locNumHitsUsedFordEdx_FDC = locdEdxHits_FDC.size()/2;
	if(locNumHitsUsedFordEdx_FDC > 0)
	{
		for(unsigned int loc_i = 0; loc_i < locNumHitsUsedFordEdx_FDC; ++loc_i)
		{
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
	//Loop over bcal showers
	vector<DBCALShowerMatchParams> locShowerMatchParamsVector;
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		DBCALShowerMatchParams locShowerMatchParams;
		if(MatchToBCAL(NULL, rt, locBCALShowers[loc_i], locStartTime, locShowerMatchParams))
			locShowerMatchParamsVector.push_back(locShowerMatchParams);
	}
	if(locShowerMatchParamsVector.empty())
		return false;

	DBCALShowerMatchParams locBestMatchParams;
	Get_BestBCALMatchParams(rt->swim_steps[0].mom, locShowerMatchParamsVector, locBestMatchParams);
	locStartTime = locBestMatchParams.dBCALShower->t - locBestMatchParams.dFlightTime;
//	locTimeVariance = locBestMatchParams.dFlightTimeVariance + locBestMatchParams.dBCALShower->dCovarianceMatrix(4, 4); //uncomment when ready!!
	locTimeVariance = 0.5*0.5;

	return true;
}

bool DParticleID::MatchToBCAL(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DBCALShower* locBCALShower, double locInputStartTime, DBCALShowerMatchParams& locShowerMatchParams) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!
	// Get the BCAL cluster position and normal
	DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z); 

	double locFlightTime = 9.9E9, locPathLength = 9.9E9;
	double d = rt->DistToRTwithTime(bcal_pos, &locPathLength, &locFlightTime, SYS_BCAL);

	if(!isfinite(d))
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locBCALShower->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	DVector3 proj_pos = rt->GetLastDOCAPoint();
	if(proj_pos.Perp() < 65.0)
		return false;  // not inside BCAL!

	double dz = bcal_pos.z() - proj_pos.z();
	double dphi = bcal_pos.Phi() - proj_pos.Phi();
	double p = rt->swim_steps[0].mom.Mag();
	//dphi += 0.002 + 8.314e-3/(p + 0.3788)/(p + 0.3788);
	while(dphi >	M_PI)
		dphi -= M_TWO_PI;
	while(dphi < -M_PI)
		dphi += M_TWO_PI;
	double phi_cut = BCAL_PHI_CUT_PAR1 + BCAL_PHI_CUT_PAR2/(p*p);

	if((fabs(dz) >= BCAL_Z_CUT) || (fabs(dphi) >= phi_cut))
		return false; //not close enough

//	if (locPathLength<0.) _DBG_ << " s " << locPathLength << " t " << locFlightTime <<endl;

	//successful match
	locShowerMatchParams.dTrack = locTrack;
	locShowerMatchParams.dBCALShower = locBCALShower;
	locShowerMatchParams.dx = 0.0; //SET ME!!!!
	locShowerMatchParams.dFlightTime = locFlightTime;
	locShowerMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!!
	locShowerMatchParams.dPathLength = locPathLength;
	locShowerMatchParams.dDeltaPhiToShower = dphi;
	locShowerMatchParams.dDeltaZToShower = dz;

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
	//Loop over tof points
	vector<DTOFHitMatchParams> locTOFHitMatchParamsVector;
	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
	{
		DTOFHitMatchParams locTOFHitMatchParams;
		if(MatchToTOF(NULL, rt, locTOFPoints[loc_i], locStartTime, locTOFHitMatchParams))
			locTOFHitMatchParamsVector.push_back(locTOFHitMatchParams);
	}
	if(locTOFHitMatchParamsVector.empty())
		return false;

	DTOFHitMatchParams locBestMatchParams;
	Get_BestTOFMatchParams(locTOFHitMatchParamsVector, locBestMatchParams);
	locStartTime = locBestMatchParams.dHitTime - locBestMatchParams.dFlightTime;
//	locTimeVariance = locBestMatchParams.dFlightTimeVariance + locBestMatchParams.dHitTimeVariance; //uncomment when ready!
	locTimeVariance = 0.1*0.1;

	return true;
}

bool DParticleID::MatchToTOF(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DTOFPoint* locTOFPoint, double locInputStartTime, DTOFHitMatchParams& locTOFHitMatchParams) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!

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

	//If the position in one dimension is not well-defined, compare distance only in the other direction
	//Otherwise, cut in R
	double locMatchCut_1D = 6.15;
	double locMatchCut_2D = 6.15;

	double locDeltaX = locTOFPoint->Is_XPositionWellDefined() ? tof_pos.X() - proj_pos.X() : 999.0;
	double locDeltaY = locTOFPoint->Is_YPositionWellDefined() ? tof_pos.Y() - proj_pos.Y() : 999.0;
	if(!locTOFPoint->Is_XPositionWellDefined())
	{
		//Is unmatched horizontal paddle with only one hit above threshold: Only compare y-distance
		if(fabs(locDeltaY) > locMatchCut_1D)
			return false;
	}
	else if(!locTOFPoint->Is_YPositionWellDefined())
	{
		//Is unmatched vertical paddle with only one hit above threshold: Only compare x-distance
		if(fabs(locDeltaX) > locMatchCut_1D)
			return false;
	}
	else
	{
		//Both are good, cut on R
		double locDistance = (tof_pos - proj_pos).Perp();
		if(locDistance > locMatchCut_2D)
			return false;
	}

	//SUCCESSFUL MATCH

	//If position was not well-defined, correct deposited energy due to attenuation, and time due to propagation along paddle
		//These values were reported at the midpoint of the paddle
	float locHitEnergy = locTOFPoint->dE;
	double locHitTime = locTOFPoint->t;
	double locHitTimeVariance = locTOFPoint->tErr*locTOFPoint->tErr;
	if(!locTOFPoint->Is_XPositionWellDefined())
	{
		//Is unmatched horizontal paddle with only one hit above threshold
		bool locNorthIsGoodHit = (locTOFPoint->dHorizontalBarStatus == 1); //+x
		int locBar = locTOFPoint->dHorizontalBar;
		bool locIsDoubleEndedBar = ((locBar < dTOFGeometry->FirstShortBar) || (locBar > dTOFGeometry->LastShortBar));

		//Paddle midpoint
		double locPaddleMidPoint = 0.0; //is 0 except when is single-ended bar (22 & 23)
		if(!locIsDoubleEndedBar)
			locPaddleMidPoint = locNorthIsGoodHit ? ONESIDED_PADDLE_MIDPOINT_MAG : -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;

		//delta_x = delta_x_actual - delta_x_mid
			//if end.x > 0: delta_x = (end.x - track.x) - (end.x - mid.x) = mid.x - track.x //if track.x > mid.x, delta_x < 0: decrease energy & increase time
			//if end.x < 0: delta_x = (track.x - end.x) - (mid.x - end.x) = track.x - mid.x //if track.x > mid.x, delta_x > 0: increase energy & decrease time
		double locDistanceToMidPoint = locNorthIsGoodHit ? locPaddleMidPoint - proj_pos.X() : proj_pos.X() - locPaddleMidPoint;

		//Energy
		locHitEnergy *= exp(locDistanceToMidPoint/TOF_ATTEN_LENGTH);

		//Time
		int id = 44 + locBar - 1;
		locHitTime -= locDistanceToMidPoint/propagation_speed[id];
		//locHitTimeVariance = //UPDATE ME!!!
	}
	else if(!locTOFPoint->Is_YPositionWellDefined())
	{
		//Is unmatched vertical paddle with only one hit above threshold
		bool locNorthIsGoodHit = (locTOFPoint->dVerticalBarStatus == 1); //+y
		int locBar = locTOFPoint->dVerticalBar;
		bool locIsDoubleEndedBar = ((locBar < dTOFGeometry->FirstShortBar) || (locBar > dTOFGeometry->LastShortBar));

		//Paddle midpoint
		double locPaddleMidPoint = 0.0; //is 0 except when is single-ended bar (22 & 23)
		if(!locIsDoubleEndedBar)
			locPaddleMidPoint = locNorthIsGoodHit ? ONESIDED_PADDLE_MIDPOINT_MAG : -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;

		//delta_x = delta_x_actual - delta_x_mid
			//if end.x > 0: delta_x = (end.x - track.x) - (end.x - mid.x) = mid.x - track.x //if track.x > mid.x, delta_x < 0: decrease energy & increase time
			//if end.x < 0: delta_x = (track.x - end.x) - (mid.x - end.x) = track.x - mid.x //if track.x > mid.x, delta_x > 0: increase energy & decrease time
		double locDistanceToMidPoint = locNorthIsGoodHit ? locPaddleMidPoint - proj_pos.Y() : proj_pos.Y() - locPaddleMidPoint;

		//Energy
		locHitEnergy *= exp(locDistanceToMidPoint/TOF_ATTEN_LENGTH);

		//Time
		int id = locBar - 1;
		locHitTime -= locDistanceToMidPoint/propagation_speed[id];
		//locHitTimeVariance = //UPDATE ME!!!
	}

	//Fill out match info
	double dx = 2.54*proj_mom.Mag()/proj_mom.Dot(norm);
	locTOFHitMatchParams.dTrack = locTrack;
	locTOFHitMatchParams.dTOFPoint = locTOFPoint;

	locTOFHitMatchParams.dHitTime = locHitTime;
	locTOFHitMatchParams.dHitTimeVariance = locHitTimeVariance;
	locTOFHitMatchParams.dHitEnergy = locHitEnergy;

	locTOFHitMatchParams.dEdx = locHitEnergy/dx;
	locTOFHitMatchParams.dFlightTime = locFlightTime;
	locTOFHitMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!
	locTOFHitMatchParams.dPathLength = locPathLength;
	locTOFHitMatchParams.dDeltaXToHit = locDeltaX;
	locTOFHitMatchParams.dDeltaYToHit = locDeltaY;

	return true;
}

// Given a reference trajectory from a track, predict which TOF paddles should
// fire due to the charged particle passing through the TOF planes.
bool DParticleID::PredictTOFPaddles(const DReferenceTrajectory *rt,
				    unsigned int &hbar,unsigned int &vbar,
				    DVector3 *intersection) const{
  // Initialize output variables
  vbar=0;
  hbar=0;
  // Find intersection with TOF plane given by tof_pos
  DVector3 tof_pos(0,0,dTOFGeometry->CenterMPlane);
  DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
  DVector3 proj_mom,proj_pos;
  double locPathLength = 9.9E9, locFlightTime = 9.9E9;
  if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,SYS_TOF) != NOERROR)
    return false;

  double x=proj_pos.x();
  double y=proj_pos.y();

  vbar=dTOFGeometry->y2bar(x);
  hbar=dTOFGeometry->y2bar(y);

  if (intersection) *intersection=proj_pos;

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
	//Loop over FCAL showers
	vector<DFCALShowerMatchParams> locShowerMatchParamsVector;
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		DFCALShowerMatchParams locShowerMatchParams;
		if(MatchToFCAL(NULL, rt, locFCALShowers[loc_i], locStartTime, locShowerMatchParams))
			locShowerMatchParamsVector.push_back(locShowerMatchParams);
	}
	if(locShowerMatchParamsVector.empty())
		return false;

	DFCALShowerMatchParams locBestMatchParams;
	Get_BestFCALMatchParams(locShowerMatchParamsVector, locBestMatchParams);

	locStartTime = locBestMatchParams.dFCALShower->getTime() - locBestMatchParams.dFlightTime;
//	locTimeVariance = locBestMatchParams.dFlightTimeVariance + locBestMatchParams.dFCALShower->dCovarianceMatrix(4, 4); //uncomment when ready!
	locTimeVariance = 0.5*0.5;

	return true;
}

bool DParticleID::MatchToFCAL(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DFCALShower* locFCALShower, double locInputStartTime, DFCALShowerMatchParams& locShowerMatchParams) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!

	// Find the distance of closest approach between the track trajectory
	// and the cluster position, looking for the minimum
	DVector3 fcal_pos = locFCALShower->getPosition();
	DVector3 norm(0.0, 0.0, 1.0); //normal vector for the FCAL plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, SYS_FCAL) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locFCALShower->getTime() - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	double d = (fcal_pos - proj_pos).Mag();
	double p=proj_mom.Mag();
	double cut=FCAL_CUT_PAR1+FCAL_CUT_PAR2/p;
	if(d >= cut)
		return false;

	locShowerMatchParams.dTrack = locTrack;
	locShowerMatchParams.dFCALShower = locFCALShower;
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
bool DParticleID::MatchToSC(const DReferenceTrajectory* rt, const vector<const DSCHit*>& locSCHits, double& locStartTime, double& locTimeVariance, bool locIsTimeBased) const
{
	//Loop over SC points
	vector<DSCHitMatchParams> locSCHitMatchParamsVector;
	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	{
		DSCHitMatchParams locSCHitMatchParams;
		if(MatchToSC(NULL, rt, locSCHits[loc_i], locStartTime, locSCHitMatchParams,locIsTimeBased))
			locSCHitMatchParamsVector.push_back(locSCHitMatchParams);
	}
	if(locSCHitMatchParamsVector.empty())
		return false;

	DSCHitMatchParams locBestMatchParams;
	Get_BestSCMatchParams(locSCHitMatchParamsVector, locBestMatchParams);

	locStartTime = locBestMatchParams.dHitTime - locBestMatchParams.dFlightTime;
//	locTimeVariance = locBestMatchParams.dFlightTimeVariance - locBestMatchParams.dHitTimeVariance; //uncomment when ready!
	locTimeVariance = 0.3*0.3;

	return true;
}

bool DParticleID::MatchToSC(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DSCHit* locSCHit, double locInputStartTime, DSCHitMatchParams& locSCHitMatchParams, bool locIsTimeBased, DVector3 *IntersectionPoint, DVector3 *IntersectionDir) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!
	if(sc_pos.empty() || sc_norm.empty())
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locInputStartTime - locSCHit->t) > OUT_OF_TIME_CUT)
		return false;

	// Find intersection with a "barrel" approximation for the start counter
	DVector3 proj_pos(NaN,NaN,NaN), proj_mom(NaN,NaN,NaN);
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	unsigned int sc_index=locSCHit->sector-1;
	
	if(rt->GetIntersectionWithPlane(sc_pos[sc_index][0],
					sc_norm[sc_index][0], 
					proj_pos, proj_mom, 
					&locPathLength, 
					&locFlightTime) != NOERROR)
		return false;
	// Check that the intersection isn't upstream of the paddle
	double myz = proj_pos.z();
	if (myz<sc_pos[sc_index][0].z()) return false;

	double proj_phi = proj_pos.Phi();
	//if(proj_phi < 0.0)
	//proj_phi += M_TWO_PI;

	// Look for a match in phi
	//double phi = dSCphi0 + dSCdphi*(locSCHit->sector - 1);
	double sc_dphi_cut=(locIsTimeBased)?SC_DPHI_CUT:SC_DPHI_CUT_WB;
	DVector3 average_pos=0.5*(sc_pos[sc_index][0]+sc_pos[sc_index][1]);
	double phi=average_pos.Phi();
	double dphi = phi - proj_phi; //phi could be 0 degrees & proj_phi could be 359 degrees
	while(dphi > TMath::Pi())
		dphi -= M_TWO_PI;
	while(dphi < -1.0*TMath::Pi())
		dphi += M_TWO_PI;
	if(fabs(dphi) >= sc_dphi_cut)
		return false; //no match
	
	// Match in phi successful, refine match in nose region were applicable

	// Length along scintillator
	double L = 0.;

	// Initialize the normal vector for the SC paddle to the long, unbent region
	DVector3 norm=sc_norm[sc_index][0];

	// Now check to see if the intersection is in the nose region and find the
	// start time
	double locCorrectedHitTime = locSCHit->t - sc_leg_tcor;
	double sc_pos0 = sc_pos[sc_index][0].z();	
	double sc_pos1 = sc_pos[sc_index][1].z();
	if(myz <= sc_pos1)
	{
	  L=myz;
	  locCorrectedHitTime -= L/C_EFFECTIVE;
	}
	else
	{
	  unsigned int num = sc_norm[sc_index].size() - 1;
	  for (unsigned int loc_i = 1; loc_i < num; ++loc_i)
	    {
	      locPathLength = 9.9E9;
	      locFlightTime = 9.9E9;
	      if(rt->GetIntersectionWithPlane(sc_pos[sc_index][loc_i],
					      sc_norm[sc_index][loc_i], 
					      proj_pos, proj_mom, 
					      &locPathLength, 
					      &locFlightTime) != NOERROR)
		continue;
	      myz = proj_pos.z();
	      norm=sc_norm[sc_index][loc_i];
	      if(myz < sc_pos[sc_index][loc_i + 1].z()){
		average_pos
		  =0.5*(sc_pos[sc_index][loc_i]+sc_pos[sc_index][loc_i+1]);
		dphi=average_pos.Phi()-proj_pos.Phi();
		while(dphi > TMath::Pi())
		  dphi -= M_TWO_PI;
		while(dphi < -1.0*TMath::Pi())
		  dphi += M_TWO_PI;
		if (fabs(dphi)>sc_dphi_cut) return false;
		break;
	      }
	    }	
	  // Check for intersection point beyond nose
	  if (myz> sc_pos[sc_index][num].z()) return false;

	  
		// Note: in the following code, L does not include a correction for where the start counter starts
		// in z...	This is absorbed into locCorrectedHitTime, above.
		if(myz < sc_pos1)
		{
		        L = sc_pos1;
			locCorrectedHitTime -= L/C_EFFECTIVE;	
		}
		else
		{
			L = (myz - sc_pos1)*sc_angle_cor + sc_pos1;
			locCorrectedHitTime -= L/C_EFFECTIVE;
		}
	}

	double dx = 0.3*proj_mom.Mag()/fabs(proj_mom.Dot(norm));

	// For the dEdx measurement we now need to take into account that L does not 
	// compensate for the position in z at which the start counter paddle starts
	locSCHitMatchParams.dTrack = locTrack;
	locSCHitMatchParams.dSCHit = locSCHit;
	locSCHitMatchParams.dHitEnergy = (locSCHit->dE)*exp((L - sc_pos0)/ATTEN_LENGTH);
	locSCHitMatchParams.dEdx = locSCHitMatchParams.dHitEnergy/dx;
	locSCHitMatchParams.dHitTime = locCorrectedHitTime;
	locSCHitMatchParams.dHitTimeVariance = 0.0; //SET ME!!!
	locSCHitMatchParams.dFlightTime = locFlightTime;
	locSCHitMatchParams.dFlightTimeVariance = 0.0; //SET ME!!!
	locSCHitMatchParams.dPathLength = locPathLength;
	locSCHitMatchParams.dDeltaPhiToHit = dphi;

	// Optionally output intersection position	
	if (IntersectionPoint!=NULL){
	  *IntersectionPoint=proj_pos;
	}
	if (IntersectionDir!=NULL){
	  *IntersectionDir=proj_mom;
	  IntersectionDir->SetMag(1.);
	}

	return true;
}

// Predict the start counter paddle that would match a track whose reference 
// trajectory is given by rt.
unsigned int DParticleID::PredictSCSector(const DReferenceTrajectory* rt, 
					  const double dphi_cut) const
{
  if(sc_pos.empty() || sc_norm.empty())
    return 0;

  unsigned int best_sc_index=0;
  double min_dphi=1e6;
  // loop over geometry for all SC paddles looking for track intersections
  for (unsigned int sc_index=0;sc_index<30;sc_index++){
    // Find intersection with the leg region of the start counter
    DVector3 proj_pos, proj_mom(NaN,NaN,NaN);
    if(rt->GetIntersectionWithPlane(sc_pos[sc_index][0],sc_norm[sc_index][0], 
				    proj_pos, proj_mom) != NOERROR)
      continue;

    // Check that the intersection isn't upstream of the paddle
    double myz = proj_pos.z();
    if (myz<sc_pos[sc_index][0].z()) continue;

    // Compute phi difference
    double proj_phi = proj_pos.Phi();
    DVector3 average_pos=0.5*(sc_pos[sc_index][0]+sc_pos[sc_index][1]);
    double phi=average_pos.Phi();
    double dphi = phi - proj_phi; //phi could be 0 degrees & proj_phi could be 359 degrees

    while(dphi > TMath::Pi())
      dphi -= M_TWO_PI;
    while(dphi < -1.0*TMath::Pi())
      dphi += M_TWO_PI;
    if(fabs(dphi) >= SC_DPHI_CUT_WB)
      continue;
    
    // Match in phi successful, refine match in nose region were applicable
    // Initialize the normal vector for the SC paddle to the long, unbent region
    DVector3 norm=sc_norm[sc_index][0];
    double sc_pos1 = sc_pos[sc_index][1].z();
    if(myz <= sc_pos1){
      if (fabs(dphi)<min_dphi){
	best_sc_index=sc_index;
	min_dphi=fabs(dphi);
      }
    }
    else{
      bool got_match=true;
      unsigned int num = sc_norm[sc_index].size() - 1;
      for (unsigned int loc_i = 1; loc_i < num; ++loc_i){
	if(rt->GetIntersectionWithPlane(sc_pos[sc_index][loc_i],
					sc_norm[sc_index][loc_i], 
					proj_pos, proj_mom) != NOERROR)
	  continue;
	myz = proj_pos.z();
	norm=sc_norm[sc_index][loc_i];
	if(myz < sc_pos[sc_index][loc_i + 1].z()){
	  average_pos
	    =0.5*(sc_pos[sc_index][loc_i]+sc_pos[sc_index][loc_i+1]);
	  dphi=average_pos.Phi()-proj_pos.Phi();
	  while(dphi > TMath::Pi())
	    dphi -= M_TWO_PI;
	  while(dphi < -1.0*TMath::Pi())
	    dphi += M_TWO_PI;
	  if (fabs(dphi)>SC_DPHI_CUT_WB){
	    got_match=false;
	    break;
	  }
	}
      }	
      if (got_match==false) continue; //doesn't match well with nose region

      // Check for intersection point beyond nose
      if (myz> sc_pos[sc_index][num].z()) continue;

      if (fabs(dphi)<min_dphi){
	best_sc_index=sc_index;
	min_dphi=fabs(dphi);
      }
    }
  }
  if (min_dphi<dphi_cut) return best_sc_index+1;

  return 0;
}




// Match to Start Counter using straight-line projection to start counter planes
bool 
DParticleID::MatchToSC(const DKinematicData *kd,
		       const vector<const DSCHit*>& locSCHits,
		       vector<DSCHitMatchParams>& locSCHitMatchParams) const{
  
  if(sc_pos.empty() || sc_norm.empty())
    return false;

  // Ends in z of the straight portion	
  double sc_pos0 = sc_pos[0][0].z();
  double sc_pos1 = sc_pos[0][1].z();

  DVector3 mom=kd->momentum();
  DVector3 pos=kd->position();
  DVector3 cylpos[2];

  if (finder->FindIntersectionsWithCylinder(sc_pos[0][1].x(),mom,pos,cylpos[0],
					    cylpos[1])){
    for (unsigned int j=0;j<2;j++){
      double cyl_phi=cylpos[j].Phi();
      for (size_t i=0;i<locSCHits.size();i++){
	const DSCHit *locSCHit=locSCHits[i];
	// Look for a match in phi
	double phi = dSCphi0 + dSCdphi*(locSCHit->sector - 1);
	
	// First intersection point of line with cylinder
	double dphi = cyl_phi - phi; //phi could be 0 degrees & cyl_phi could be 359 degrees
	while(dphi > TMath::Pi())
	  dphi -= M_TWO_PI;
	while(dphi < -1.0*TMath::Pi())
	  dphi += M_TWO_PI;
	if(fabs(dphi) < SC_DPHI_CUT_WB){
	  double myz=cylpos[j].z();
	  double sc_time=locSCHit->t - sc_leg_tcor;
	  if (myz>=sc_pos0){
	    DSCHitMatchParams scmatch;
	    scmatch.dSCHit=locSCHit;
	    scmatch.dDeltaPhiToHit=dphi; 
	    if (myz<sc_pos1){ // intersection in leg 
	      scmatch.dHitTime=sc_time-myz/C_EFFECTIVE;
	    }
	    else{ // intersection in nose
	      unsigned int num = sc_norm[0].size() - 1;
	      for (unsigned int k = 1; k < num; ++k){
		double xhat = sc_norm[0][k].x();
		double cosphi=cos(phi);
		double sinphi=sin(phi);
		DVector3 norm(cosphi*xhat, sinphi*xhat, sc_norm[0][k].z());
		double r = sc_pos[0][k].X();
		DVector3 origin(r*cosphi, r*sinphi, sc_pos[0][k].z());
		if (finder->FindIntersectionWithPlane(origin,norm,pos,mom,
						      cylpos[j])){
		  myz = cylpos[j].z();
		  if(myz < sc_pos[0][k + 1].z()){
		    break;
		  }
		}
	      }
	      if (myz<sc_pos0) continue;
	      if (myz>sc_pos[0][num].z()) continue;

	      // Note: in the following code, L does not include a correction
	      // for where the start counter starts in z...	
	      // This is absorbed into sc_time, above.
	      if(myz < sc_pos1){
		scmatch.dHitTime=sc_time-sc_pos1/C_EFFECTIVE;
	      }
	      else{
		double L = (myz - sc_pos1)*sc_angle_cor + sc_pos1;
		scmatch.dHitTime=sc_time-L/C_EFFECTIVE;
	      }
	    }
	    locSCHitMatchParams.push_back(scmatch);
	  } // back position check	
	} // dphi check
      } // loop over start counter hits
    } // two intersection points    
  } // intersection with cylinder

  return true;
}
  
//------------------
// Distance_ToTrack
//------------------
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// tracks can be skipped
bool DParticleID::Distance_ToTrack(const DBCALShower* locBCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance, double& locDeltaPhi, double& locDeltaZ) const
{
	// Get the BCAL cluster position and normal
	DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z); 

	double locFlightTime = 9.9E9, locPathLength = 9.9E9;
	locDistance = rt->DistToRTwithTime(bcal_pos, &locPathLength, &locFlightTime,SYS_BCAL);
	if(!isfinite(locDistance))
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locBCALShower->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	DVector3 proj_pos = rt->GetLastDOCAPoint();
	if(proj_pos.Perp() < 65.0)
		return false;  // not inside BCAL!

	locDeltaZ = proj_pos.z() - bcal_pos.z();
	locDeltaPhi = proj_pos.Phi() - bcal_pos.Phi();
	while(locDeltaPhi >	M_PI)
		locDeltaPhi -= M_TWO_PI;
	while(locDeltaPhi < -M_PI)
		locDeltaPhi += M_TWO_PI;

	return true;
}

//------------------
// Distance_ToTrack
//------------------
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// tracks can be skipped
bool DParticleID::Distance_ToTrack(const DFCALShower* locFCALShower, const DReferenceTrajectory* rt, double locInputStartTime, double& locDistance) const
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

//------------------
// Distance_ToTrack
//------------------
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// tracks can be skipped
bool DParticleID::Distance_ToTrack(const DTOFPoint* locTOFPoint, const DReferenceTrajectory* rt, double locInputStartTime, double& locDeltaX, double& locDeltaY) const
{
	DVector3 tof_pos = locTOFPoint->pos;
	DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, SYS_TOF) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locTOFPoint->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	locDeltaX = locTOFPoint->Is_XPositionWellDefined() ? tof_pos.X() - proj_pos.X() : 999.0;
	locDeltaY = locTOFPoint->Is_YPositionWellDefined() ? tof_pos.Y() - proj_pos.Y() : 999.0;
	return true;
}

//------------------
// Distance_ToTrack
//------------------
// NOTE: an initial guess for start time is expected as input so that out-of-time 
// tracks can be skipped
bool DParticleID::Distance_ToTrack(const DSCHit* locSCHit, const DReferenceTrajectory* rt, double locInputStartTime, double& locDeltaPhi) const
{
	if(sc_pos.empty() || sc_norm.empty())
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locInputStartTime - locSCHit->t) > OUT_OF_TIME_CUT)
		return false;

	// Find intersection with a "barrel" approximation for the start counter
	DVector3 proj_pos(NaN,NaN,NaN), proj_mom(NaN,NaN,NaN);
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithRadius(sc_pos[0][1].x(), proj_pos, &locPathLength, &locFlightTime, &proj_mom) != NOERROR)
		return false;

	double proj_phi = proj_pos.Phi();
	if(proj_phi < 0.0)
		proj_phi += M_TWO_PI;

	double phi = dSCphi0 + dSCdphi*(locSCHit->sector - 1);
	locDeltaPhi = phi - proj_phi; //phi could be 0 degrees & proj_phi could be 359 degrees
	while(locDeltaPhi > TMath::Pi())
		locDeltaPhi -= M_TWO_PI;
	while(locDeltaPhi < -1.0*TMath::Pi())
		locDeltaPhi += M_TWO_PI;

	return true;
}

bool DParticleID::Get_BestSCMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DSCHitMatchParams& locBestMatchParams) const
{
	//choose the "best" detector hit to use for computing quantities
	vector<DSCHitMatchParams> locSCHitMatchParams;
	if(!locDetectorMatches->Get_SCMatchParams(locTrack, locSCHitMatchParams))
		return false;

	Get_BestSCMatchParams(locSCHitMatchParams, locBestMatchParams);
	return true;
}

void DParticleID::Get_BestSCMatchParams(vector<DSCHitMatchParams>& locSCHitMatchParams, DSCHitMatchParams& locBestMatchParams) const
{
	double locMinDeltaPhi = 9.9E9;
	for(size_t loc_i = 0; loc_i < locSCHitMatchParams.size(); ++loc_i)
	{
		if(locSCHitMatchParams[loc_i].dDeltaPhiToHit >= locMinDeltaPhi)
			continue;
		locMinDeltaPhi = locSCHitMatchParams[loc_i].dDeltaPhiToHit;
		locBestMatchParams = locSCHitMatchParams[loc_i];
	}
}

bool DParticleID::Get_BestBCALMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DBCALShowerMatchParams& locBestMatchParams) const
{
	//choose the "best" shower to use for computing quantities
	vector<DBCALShowerMatchParams> locShowerMatchParams;
	if(!locDetectorMatches->Get_BCALMatchParams(locTrack, locShowerMatchParams))
		return false;

	Get_BestBCALMatchParams(locTrack->momentum(), locShowerMatchParams, locBestMatchParams);
	return true;
}

void DParticleID::Get_BestBCALMatchParams(DVector3 locMomentum, vector<DBCALShowerMatchParams>& locShowerMatchParams, DBCALShowerMatchParams& locBestMatchParams) const
{
	double locMinChiSq = 9.9E9;
	double locP = locMomentum.Mag();
	for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
	{
		double locDeltaPhiCut = BCAL_PHI_CUT_PAR1 + BCAL_PHI_CUT_PAR2/(locP*locP);
		double locDeltaPhiError = locDeltaPhiCut/3.0; //Cut is "3 sigma"
		double locMatchChiSq = locShowerMatchParams[loc_i].dDeltaPhiToShower*locShowerMatchParams[loc_i].dDeltaPhiToShower/(locDeltaPhiError*locDeltaPhiError);

		double locDeltaZError = BCAL_Z_CUT/3.0; //Cut is "3 sigma"
		locMatchChiSq += locShowerMatchParams[loc_i].dDeltaZToShower*locShowerMatchParams[loc_i].dDeltaZToShower/(locDeltaZError*locDeltaZError);

		if(locMatchChiSq >= locMinChiSq)
			continue;

		locMinChiSq = locMatchChiSq;
		locBestMatchParams = locShowerMatchParams[loc_i];
	}
}

bool DParticleID::Get_BestTOFMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DTOFHitMatchParams& locBestMatchParams) const
{
	//choose the "best" hit to use for computing quantities
	vector<DTOFHitMatchParams> locTOFHitMatchParams;
	if(!locDetectorMatches->Get_TOFMatchParams(locTrack, locTOFHitMatchParams))
		return false;

	Get_BestTOFMatchParams(locTOFHitMatchParams, locBestMatchParams);
	return true;
}

void DParticleID::Get_BestTOFMatchParams(vector<DTOFHitMatchParams>& locTOFHitMatchParams, DTOFHitMatchParams& locBestMatchParams) const
{
	double locMinDistance = 9.9E9;
	for(size_t loc_i = 0; loc_i < locTOFHitMatchParams.size(); ++loc_i)
	{
		double locDeltaR = sqrt(locTOFHitMatchParams[loc_i].dDeltaXToHit*locTOFHitMatchParams[loc_i].dDeltaXToHit + locTOFHitMatchParams[loc_i].dDeltaYToHit*locTOFHitMatchParams[loc_i].dDeltaYToHit);
		if(locDeltaR >= locMinDistance)
			continue;
		locMinDistance = locDeltaR;
		locBestMatchParams = locTOFHitMatchParams[loc_i];
	}
}

bool DParticleID::Get_BestFCALMatchParams(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches, DFCALShowerMatchParams& locBestMatchParams) const
{
	//choose the "best" shower to use for computing quantities
	vector<DFCALShowerMatchParams> locShowerMatchParams;
	if(!locDetectorMatches->Get_FCALMatchParams(locTrack, locShowerMatchParams))
		return false;

	Get_BestFCALMatchParams(locShowerMatchParams, locBestMatchParams);
	return true;
}

void DParticleID::Get_BestFCALMatchParams(vector<DFCALShowerMatchParams>& locShowerMatchParams, DFCALShowerMatchParams& locBestMatchParams) const
{
	double locMinDistance = 9.9E9;
	for(size_t loc_i = 0; loc_i < locShowerMatchParams.size(); ++loc_i)
	{
		if(locShowerMatchParams[loc_i].dDOCAToShower >= locMinDistance)
			continue;
		locMinDistance = locShowerMatchParams[loc_i].dDOCAToShower;
		locBestMatchParams = locShowerMatchParams[loc_i];
	}
}

const DBCALShower* DParticleID::Get_ClosestToTrack_BCAL(const DKinematicData* locTrack, vector<const DBCALShower*>& locBCALShowers, double& locBestMatchDeltaPhi, double& locBestMatchDeltaZ) const
{
	const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locTrack);
	const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locTrack);
	if((locTrackTimeBased == NULL) && (locTrackWireBased == NULL))
		return NULL;
	const DReferenceTrajectory* locReferenceTrajectory = (locTrackTimeBased != NULL) ? locTrackTimeBased->rt : locTrackWireBased->rt;
	if(locReferenceTrajectory == NULL)
		return NULL;

	double locInputStartTime = locTrack->t0();
	double locMinDistance = 999.0;
	const DBCALShower* locBestBCALShower = NULL;
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		double locDistance = 0.0, locDeltaPhi = 0.0, locDeltaZ = 0.0;
		if(!Distance_ToTrack(locBCALShowers[loc_i], locReferenceTrajectory, locInputStartTime, locDistance, locDeltaPhi, locDeltaZ))
			continue;
		if(locDistance > locMinDistance)
			continue;
		locBestBCALShower = locBCALShowers[loc_i];
		locMinDistance = locDistance;
		locBestMatchDeltaPhi = locDeltaPhi;
		locBestMatchDeltaZ = locDeltaZ;
	}
	return locBestBCALShower;
}

const DFCALShower* DParticleID::Get_ClosestToTrack_FCAL(const DKinematicData* locTrack, vector<const DFCALShower*>& locFCALShowers, double& locBestDistance) const
{
	const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locTrack);
	const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locTrack);
	if((locTrackTimeBased == NULL) && (locTrackWireBased == NULL))
		return NULL;
	const DReferenceTrajectory* locReferenceTrajectory = (locTrackTimeBased != NULL) ? locTrackTimeBased->rt : locTrackWireBased->rt;
	if(locReferenceTrajectory == NULL)
		return NULL;

	double locInputStartTime = locTrack->t0();
	locBestDistance = 999.0;
	const DFCALShower* locBestFCALShower = NULL;
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		double locDistance = 0.0;
		if(!Distance_ToTrack(locFCALShowers[loc_i], locReferenceTrajectory, locInputStartTime, locDistance))
			continue;
		if(locDistance > locBestDistance)
			continue;
		locBestDistance = locDistance;
		locBestFCALShower = locFCALShowers[loc_i];
	}
	return locBestFCALShower;
}

const DTOFPoint* DParticleID::Get_ClosestToTrack_TOF(const DKinematicData* locTrack, vector<const DTOFPoint*>& locTOFPoints, double& locBestDeltaX, double& locBestDeltaY) const
{
	const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locTrack);
	const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locTrack);
	if((locTrackTimeBased == NULL) && (locTrackWireBased == NULL))
		return NULL;
	const DReferenceTrajectory* locReferenceTrajectory = (locTrackTimeBased != NULL) ? locTrackTimeBased->rt : locTrackWireBased->rt;
	if(locReferenceTrajectory == NULL)
		return NULL;

	double locMinDistance = 9.9E9;
	const DTOFPoint* locClosestTOFPoint = NULL;
	double locInputStartTime = locTrack->t0();
	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
	{
		double locDistance = 0.0, locDeltaX = 0.0, locDeltaY = 0.0;
		if(!Distance_ToTrack(locTOFPoints[loc_i], locReferenceTrajectory, locInputStartTime, locDeltaX, locDeltaY))
			continue;
		locDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
		if(locDistance > locMinDistance)
			continue;
		locMinDistance = locDistance;
		locBestDeltaX = locDeltaX;
		locBestDeltaY = locDeltaY;
		locClosestTOFPoint = locTOFPoints[loc_i];
	}

	return locClosestTOFPoint;
}

const DSCHit* DParticleID::Get_ClosestToTrack_SC(const DKinematicData* locTrack, vector<const DSCHit*>& locSCHits, double& locBestDeltaPhi) const
{
	const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locTrack);
	const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locTrack);
	if((locTrackTimeBased == NULL) && (locTrackWireBased == NULL))
		return NULL;
	const DReferenceTrajectory* locReferenceTrajectory = (locTrackTimeBased != NULL) ? locTrackTimeBased->rt : locTrackWireBased->rt;
	if(locReferenceTrajectory == NULL)
		return NULL;

	double locInputStartTime = locTrack->t0();
	const DSCHit* locBestSCHit = NULL;

	locBestDeltaPhi = TMath::Pi();
	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	{
		double locDeltaPhi = 0.0;
		if(!Distance_ToTrack(locSCHits[loc_i], locReferenceTrajectory, locInputStartTime, locDeltaPhi))
			continue;
		if(locDeltaPhi > locBestDeltaPhi)
			continue;
		locBestDeltaPhi = locDeltaPhi;
		locBestSCHit = locSCHits[loc_i];
	}

	return locBestSCHit;
}

double DParticleID::Calc_BCALFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const
{
	DBCALShowerMatchParams locBCALShowerMatchParams;
	if(!Get_BestBCALMatchParams(locTrack, locDetectorMatches, locBCALShowerMatchParams))
		return numeric_limits<double>::quiet_NaN();
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_FCALFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const
{
	DFCALShowerMatchParams locFCALShowerMatchParams;
	if(!Get_BestFCALMatchParams(locTrack, locDetectorMatches, locFCALShowerMatchParams))
		return numeric_limits<double>::quiet_NaN();
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_TOFFlightTimePCorrelation(const DKinematicData* locTrack, DDetectorMatches* locDetectorMatches) const
{
	DTOFHitMatchParams locTOFHitMatchParams;
	if(!Get_BestTOFMatchParams(locTrack, locDetectorMatches, locTOFHitMatchParams))
		return numeric_limits<double>::quiet_NaN();
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_SCFlightTimePCorrelation(const DKinematicData* locTrack, const DDetectorMatches* locDetectorMatches) const
{
	DSCHitMatchParams locSCHitMatchParams;
	if(!Get_BestSCMatchParams(locTrack, locDetectorMatches, locSCHitMatchParams))
		return numeric_limits<double>::quiet_NaN();
	double locFlightTimePCorrelation = 0.0; //SET ME!!!
	return locFlightTimePCorrelation;
}

double DParticleID::Calc_PropagatedRFTime(const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch) const
{
	//Propagate RF time to the track vertex-z
	return locEventRFBunch->dTime + (locKinematicData->z() - dTargetZCenter)/SPEED_OF_LIGHT;
}

double DParticleID::Calc_TimingChiSq(const DKinematicData* locKinematicData, unsigned int &locNDF, double& locPull) const
{
	if((locKinematicData->t0_detector() == SYS_NULL) || (locKinematicData->t1_detector() == SYS_NULL))
	{
		// not matched to any hits
		locNDF = 0;
		locPull = 0.0;
		return 0.0;
	}

	double locStartTimeError = locKinematicData->t0_err();
	double locTimeDifferenceVariance = (locKinematicData->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
	locPull = (locKinematicData->t0() - locKinematicData->time())/sqrt(locTimeDifferenceVariance);
	locNDF = 1;
	return locPull*locPull;
}

void DParticleID::Calc_ChargedPIDFOM(DChargedTrackHypothesis* locChargedTrackHypothesis, const DEventRFBunch* locEventRFBunch) const
{
	CalcDCdEdxChiSq(locChargedTrackHypothesis);

	unsigned int locTimingNDF = 0;
	double locTimingPull = 0.0;
	double locTimingChiSq = Calc_TimingChiSq(locChargedTrackHypothesis, locTimingNDF, locTimingPull);
	locChargedTrackHypothesis->dChiSq_Timing = locTimingChiSq;
	locChargedTrackHypothesis->dNDF_Timing = locTimingNDF;

	unsigned int locNDF_Total = locChargedTrackHypothesis->dNDF_Timing + locChargedTrackHypothesis->dNDF_DCdEdx;
	double locChiSq_Total = locChargedTrackHypothesis->dChiSq_Timing + locChargedTrackHypothesis->dChiSq_DCdEdx;

	locChargedTrackHypothesis->dChiSq = locChiSq_Total;
	locChargedTrackHypothesis->dNDF = locNDF_Total;
	locChargedTrackHypothesis->dFOM = (locNDF_Total > 0) ? TMath::Prob(locChiSq_Total, locNDF_Total) : numeric_limits<double>::quiet_NaN();
}

unsigned int DParticleID::Get_CDCRingBitPattern(vector<const DCDCTrackHit*>& locCDCTrackHits) const
{
	unsigned int locBitPattern = 0;
	for(size_t loc_i = 0; loc_i < locCDCTrackHits.size(); ++loc_i)
	{
		//bit-shift to get bit for this ring, then bitwise-or
		unsigned int locBitShift = locCDCTrackHits[loc_i]->wire->ring - 1;
		unsigned int locBit = (1 << locBitShift); //if ring # is 1, shift 1 by 0
		locBitPattern |= locBit;
	}
	return locBitPattern;
}

unsigned int DParticleID::Get_FDCPlaneBitPattern(vector<const DFDCPseudo*>& locFDCPseudos) const
{
	unsigned int locBitPattern = 0;
	for(size_t loc_i = 0; loc_i < locFDCPseudos.size(); ++loc_i)
	{
		//bit-shift to get bit for this ring, then bitwise-or
		unsigned int locBitShift = locFDCPseudos[loc_i]->wire->layer - 1;
		unsigned int locBit = (1 << locBitShift); //if ring # is 1, shift 1 by 0
		locBitPattern |= locBit;
	}
	return locBitPattern;
}

void DParticleID::Get_CDCRings(int locBitPattern, set<int>& locCDCRings) const
{
	locCDCRings.clear();
	for(unsigned int locRing = 1; locRing <= 28; ++locRing)
	{
		unsigned int locBitShift = locRing - 1;
		unsigned int locBit = (1 << locBitShift); //if ring # is 1, shift 1 by 0
		if((locBitPattern & locBit) != 0)
			locCDCRings.insert(locRing);
	}
}

void DParticleID::Get_FDCPlanes(int locBitPattern, set<int>& locFDCPlanes) const
{
	locFDCPlanes.clear();
	for(unsigned int locPlane = 1; locPlane <= 24; ++locPlane)
	{
		unsigned int locBitShift = locPlane - 1;
		unsigned int locBit = (1 << locBitShift); //if ring # is 1, shift 1 by 0
		if((locBitPattern & locBit) != 0)
			locFDCPlanes.insert(locPlane);
	}
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

