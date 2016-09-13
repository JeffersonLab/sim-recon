// $Id$
//
//    File: DParticleID.cc
// Created: Mon Feb 28 14:48:56 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID.h"
#include "START_COUNTER/DSCHit_factory.h"

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
  dSCdphi=12.0*M_PI/180.;  // 12 degrees

	C_EFFECTIVE = 15.0;
	ATTEN_LENGTH = 150.0;
	OUT_OF_TIME_CUT = 35.0; // Changed 200 -> 35 ns, March 2016
    gPARMS->SetDefaultParameter("PID:OUT_OF_TIME_CUT",OUT_OF_TIME_CUT);

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

  // Get z position of face of FCAL
  locGeometry->GetFCALZ(dFCALz);

  // Get start counter geometry;
  if (locGeometry->GetStartCounterGeom(sc_pos,sc_norm)){
    dSCphi0=sc_pos[0][0].Phi();
    //sc_leg_tcor = -sc_pos[0][0].z()/C_EFFECTIVE;
    double theta = sc_norm[0][sc_norm[0].size()-2].Theta();
    sc_angle_cor = 1./cos(M_PI_2 - theta);

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
	FCAL_CUT_PAR1=4.5;
	gPARMS->SetDefaultParameter("FCAL:CUT_PAR1",FCAL_CUT_PAR1);

	FCAL_CUT_PAR2=0.0;
	gPARMS->SetDefaultParameter("FCAL:CUT_PAR2",FCAL_CUT_PAR2);

	BCAL_Z_CUT=50.;
	gPARMS->SetDefaultParameter("BCAL:Z_CUT",BCAL_Z_CUT);

	BCAL_PHI_CUT_PAR1=0.021;
	gPARMS->SetDefaultParameter("BCAL:PHI_CUT_PAR1",BCAL_PHI_CUT_PAR1);

	BCAL_PHI_CUT_PAR2=0.01;
	gPARMS->SetDefaultParameter("BCAL:PHI_CUT_PAR2",BCAL_PHI_CUT_PAR2);

	BCAL_PHI_CUT_PAR3=1280.;
	gPARMS->SetDefaultParameter("BCAL:PHI_CUT_PAR3",BCAL_PHI_CUT_PAR3);

	SC_DPHI_CUT=0.125;
	gPARMS->SetDefaultParameter("SC:DPHI_CUT",SC_DPHI_CUT);
	
	SC_DPHI_CUT_SLOPE=0.004;
	gPARMS->SetDefaultParameter("SC:DPHI_CUT_SLOPE",SC_DPHI_CUT_SLOPE);

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
  
  // FCAL geometry
  loop->GetSingle(dFCALGeometry);

	//TOF calibration constants & geometry
	if(loop->GetCalib("TOF/propagation_speed", propagation_speed))
		jout << "Error loading /TOF/propagation_speed !" << endl;

	map<string, double> tofparms;
 	loop->GetCalib("TOF/tof_parms", tofparms);
	TOF_ATTEN_LENGTH = tofparms["TOF_ATTEN_LENGTH"];
	TOF_E_THRESHOLD = tofparms["TOF_E_THRESHOLD"];
	TOF_HALFPADDLE = tofparms["TOF_HALFPADDLE"];

	loop->GetSingle(dTOFGeometry);
	dHalfPaddle_OneSided = dTOFGeometry->SHORTBARLENGTH/2.0; //GET FROM GEOMETRY??
	double locBeamHoleWidth = dTOFGeometry->LONGBARLENGTH - 2.0*dTOFGeometry->SHORTBARLENGTH;
	ONESIDED_PADDLE_MIDPOINT_MAG = dHalfPaddle_OneSided + locBeamHoleWidth/2.0;

	// Start counter calibration constants
	// vector<map<string,double> > tvals;
	vector<map<string,double> > pt_vals;
	vector<map<string,double> > attn_vals;

	// if(loop->GetCalib("/START_COUNTER/propagation_speed",tvals))
	//   jout << "Error loading /START_COUNTER/propagation_speed !" << endl;
	// else{
	//   for(unsigned int i=0; i<tvals.size(); i++){
        //     map<string, double> &row = tvals[i];
	//     sc_veff[SC_STRAIGHT].push_back(row["SC_STRAIGHT_PROPAGATION_B"]);
	//     sc_veff[SC_BEND].push_back(row["SC_BEND_PROPAGATION_B"]);
	//     sc_veff[SC_NOSE].push_back(row["SC_NOSE_PROPAGATION_B"]);
	//   }
	// }

	// Individual propagation speed calibrations (beam data)
	if(loop->GetCalib("/START_COUNTER/propagation_speed", pt_vals))
	  jout << "Error loading /START_COUNTER/propagation_speed !" << endl;
	else
	  {
	    for(unsigned int i = 0; i < pt_vals.size(); i++)
	      {
		// Functional form is: A + B*x
		map<string, double> &row = pt_vals[i];
		sc_pt_yint[SC_STRAIGHT].push_back(row["SC_STRAIGHT_PROPAGATION_A"]);
		sc_pt_yint[SC_BEND].push_back(row["SC_BEND_PROPAGATION_A"]);
		sc_pt_yint[SC_NOSE].push_back(row["SC_NOSE_PROPAGATION_A"]);
		    
		sc_pt_slope[SC_STRAIGHT].push_back(row["SC_STRAIGHT_PROPAGATION_B"]);
		sc_pt_slope[SC_BEND].push_back(row["SC_BEND_PROPAGATION_B"]);
		sc_pt_slope[SC_NOSE].push_back(row["SC_NOSE_PROPAGATION_B"]);
	      }
	  }

	// Individual attneuation calibrations (FIU bench mark data) 
	if(loop->GetCalib("START_COUNTER/attenuation_factor", attn_vals))
	  jout << "Error in loading START_COUNTER/attenuation_factor !" << endl;
	else
	  {
	    for(unsigned int i = 0; i < attn_vals.size(); i++)
	      {
		// Functional form is: A*exp(B*x) + C
		map<string, double> &row = attn_vals[i];
		sc_attn_A[SC_STRAIGHT_ATTN].push_back(row["SC_STRAIGHT_ATTENUATION_A"]);
		sc_attn_A[SC_BENDNOSE_ATTN].push_back(row["SC_BENDNOSE_ATTENUATION_A"]); 

		sc_attn_B[SC_STRAIGHT_ATTN].push_back(row["SC_STRAIGHT_ATTENUATION_B"]);
		sc_attn_B[SC_BENDNOSE_ATTN].push_back(row["SC_BENDNOSE_ATTENUATION_B"]); 

		sc_attn_C[SC_STRAIGHT_ATTN].push_back(row["SC_STRAIGHT_ATTENUATION_C"]);
		sc_attn_C[SC_BENDNOSE_ATTN].push_back(row["SC_BENDNOSE_ATTENUATION_C"]); 
	      }
	  }

    // Start counter individual paddle resolutions
    if(loop->GetCalib("START_COUNTER/time_resol_paddle", sc_paddle_resols))
        jout << "Error in loading START_COUNTER/time_resol_paddle !" << endl;
	else {
        if(sc_paddle_resols.size() != (unsigned int)DSCHit_factory::MAX_SECTORS)
            jerr << "Start counter paddle resolutions table has wrong number of entries:" << endl
                 << "  loaded = " << sc_paddle_resols.size() 
                 << "  expexted = " << DSCHit_factory::MAX_SECTORS << endl;
    }

	//be sure that DRFTime_factory::init() and brun() are called
	vector<const DTOFPoint*> locTOFPoints;
	loop->Get(locTOFPoints);

	dTOFPointFactory = static_cast<DTOFPoint_factory*>(loop->GetFactory("DTOFPoint"));
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

// Calculate the most probable energy loss per unit length in units of
// GeV/cm in the FDC or CDC gas for a particle of momentum p and mass "mass"
double DParticleID::GetMostProbabledEdx_DC(double p,double mass,double dx, bool locIsCDCFlag) const{
  double betagamma=p/mass;
  double beta2=1./(1.+1./betagamma/betagamma);
  if (beta2<1e-6) beta2=1e-6;

  // Electron mass
  double Me=0.000511; //GeV

  double locKRhoZoverAGas = locIsCDCFlag ? dKRhoZoverA_CDC : dKRhoZoverA_FDC;
  double locLnIGas = locIsCDCFlag ? dLnI_CDC : dLnI_FDC;

  // First (non-logarithmic) term in Bethe-Bloch formula
  double mean_dedx=locKRhoZoverAGas/beta2;

  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
  return mean_dedx*(log(mean_dedx*dx)
		    -log((1.-beta2)/2./Me/beta2)-2.*locLnIGas-beta2+0.200);
}

// Empirical form for sigma/mean for gaseous detectors with num_dedx
// samples and sampling thickness path_length.  Assumes that the number of
// hits has already been converted from an (unsigned) int to a double.
double DParticleID::GetdEdxSigma_DC(double num_hits,double p,double mass,
				  double mean_path_length, bool locIsCDCFlag) const{
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
  double T0=(Tmax>100.e-6)?100.e-6:Tmax;  //100.e-6: energy cut for bethe-bloch

  // Bethe-Bloch
  double locKRhoZoverAGas = locIsCDCFlag ? dKRhoZoverA_CDC : dKRhoZoverA_FDC;
  double locLnIGas = locIsCDCFlag ? dLnI_CDC : dLnI_FDC;
  double mean_dedx=locKRhoZoverAGas/beta2
    *(log(two_Me_betagamma_sq*T0)-2.*locLnIGas-beta2*(1.+T0/Tmax));

  return 0.41*mean_dedx*pow(num_hits,-0.43)*pow(mean_path_length,-0.32);
  //return 0.41*mean_dedx*pow(double(num_hits),-0.5)*pow(mean_path_length,-0.32);
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
	locTimeVariance = 0.3*0.3+locBestMatchParams.dFlightTimeVariance;

	return true;
}


bool DParticleID::MatchToBCAL(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DBCALShower* locBCALShower, double locInputStartTime, DBCALShowerMatchParams& locShowerMatchParams) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!
	// Get the BCAL cluster position and normal
	DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z); 

	double locFlightTime = 9.9E9, locPathLength = 9.9E9;
	double locFlightTimeVariance=9.9E9;
	double d = rt->DistToRTwithTime(bcal_pos, &locPathLength,
					&locFlightTime, &locFlightTimeVariance,
					SYS_BCAL);

	if(!isfinite(d))
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locBCALShower->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;
	
	DVector3 proj_pos = rt->GetLastDOCAPoint();
	//if(proj_pos.Perp() < 64.0)
	//  return false;  // not inside BCAL!

	// Rough phi cut
	double dphi=bcal_pos.Phi()-proj_pos.Phi();
	while(dphi > M_PI)
	  dphi -= M_TWO_PI;
	while(dphi < -M_PI)
	  dphi += M_TWO_PI;
	
	if (fabs(dphi)>0.26) //0.13 radians = 7.5 degrees = one bcal module
	  return false;
	
	// rough z cut
	double dz=proj_pos.z()-locBCALShower->z;
	if (fabs(dz)>BCAL_Z_CUT)
	  return false; 

	// Get clusters associated with this shower
	vector<const DBCALCluster*>clusters;
	locBCALShower->Get(clusters);

    // make list of points associated with the shower
    vector<const DBCALPoint*> points;
    if(clusters.size() > 0) {
        // classic BCAL shower objects are built from the output of the clusterizer
        // so the points need to be accessed as shower -> cluster -> points
        for (unsigned int k=0;k<clusters.size();k++){
            vector<const DBCALPoint*> cluster_points=clusters[k]->points();
            points.insert(points.end(), cluster_points.begin(), cluster_points.end());
        }
    } else {
        // other BCAL shower objects directly keep a list of the points associated with the shower
        // (e.g. "CURVATURE" showers)
        locBCALShower->Get(points);
    }

    // loop over points associated with this shower, finding 
    // the closest match between a point and the track
    double dphi_min=1000.;
    double z_for_dphi_min=1000.;
    for (unsigned int m=0;m<points.size();m++){
      double rpoint=points[m]->r();
      if (rt->GetIntersectionWithRadius(rpoint,proj_pos)==NOERROR){
        dphi=points[m]->phi()-proj_pos.Phi();
        while(dphi > M_PI)
	  dphi -= M_TWO_PI;
	    while(dphi < -M_PI)
		  dphi += M_TWO_PI;
	    if (fabs(dphi)<fabs(dphi_min)){
	      dphi_min=dphi;
		  z_for_dphi_min=proj_pos.z();
	    }
      }
    }

	double p = rt->swim_steps[0].mom.Mag();
	double phi_cut = (BCAL_PHI_CUT_PAR1 + BCAL_PHI_CUT_PAR2/(p*p))
	*(1.+BCAL_PHI_CUT_PAR3*pow(450.-z_for_dphi_min,-2));
	// Look for a match in phi 
	if (fabs(dphi_min) >= phi_cut)
	  return false; //not close enough

//	if (locPathLength<0.) _DBG_ << " s " << locPathLength << " t " << locFlightTime <<endl;

	//successful match
	locShowerMatchParams.dTrack = locTrack;
	locShowerMatchParams.dBCALShower = locBCALShower;
	locShowerMatchParams.dx = 0.0; //SET ME!!!!
	locShowerMatchParams.dFlightTime = locFlightTime;
	locShowerMatchParams.dFlightTimeVariance = locFlightTimeVariance;
	locShowerMatchParams.dPathLength = locPathLength;
	locShowerMatchParams.dDeltaPhiToShower = dphi_min;
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
	locTimeVariance = 0.1*0.1+locBestMatchParams.dFlightTimeVariance;

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
	double locFlightTimeVariance=9.9E9;
	if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,&locFlightTimeVariance,SYS_TOF) != NOERROR)
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

	//GEOMETRIC MATCH

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

	// Check that the hit is not out of time with respect to the track
	double locDeltaT = locHitTime - locFlightTime - locInputStartTime;
	if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
		return false;

	//SUCCESSFUL MATCH

	//Fill out match info
	double dx = 2.54*proj_mom.Mag()/proj_mom.Dot(norm);
	locTOFHitMatchParams.dTrack = locTrack;
	locTOFHitMatchParams.dTOFPoint = locTOFPoint;

	locTOFHitMatchParams.dHitTime = locHitTime;
	locTOFHitMatchParams.dHitTimeVariance = locHitTimeVariance;
	locTOFHitMatchParams.dHitEnergy = locHitEnergy;

	locTOFHitMatchParams.dEdx = locHitEnergy/dx;
	locTOFHitMatchParams.dFlightTime = locFlightTime;
	locTOFHitMatchParams.dFlightTimeVariance = locFlightTimeVariance;
	locTOFHitMatchParams.dPathLength = locPathLength;
	locTOFHitMatchParams.dDeltaXToHit = locDeltaX;
	locTOFHitMatchParams.dDeltaYToHit = locDeltaY;

	return true;
}

// Given a reference trajectory from a track, predict which FCAL block (if any)
// should fire.
bool DParticleID::PredictFCALHit(const DReferenceTrajectory *rt,
				 unsigned int &row, unsigned int &col,
				 DVector3 *intersection) const{
   // Initialize output variables
  row=0;
  col=0;
  if(rt == NULL)
    return false;
  // Find intersection with FCAL plane given by fcal_pos
  DVector3 fcal_pos(0,0,dFCALz);
  DVector3 norm(0.0, 0.0, 1.0); //normal vector to FCAL plane
  DVector3 proj_mom,proj_pos;
  if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom,
				  NULL,NULL,NULL,SYS_FCAL) != NOERROR)
    return false;

  if (intersection) *intersection=proj_pos;

  double x=proj_pos.x();
  double y=proj_pos.y();
  row=dFCALGeometry->row(float(y));
  col=dFCALGeometry->column(float(x));
  return (dFCALGeometry->isBlockActive(row,col));
}

// Given a reference trajectory from a track, predict which BCAL wedge should
// have a hit
bool DParticleID::PredictBCALWedge(const DReferenceTrajectory *rt,
				   unsigned int &module,unsigned int &sector,
				   DVector3 *intersection) const{
  //initialize output variables 
  sector=0;
  module=0;
  if(rt == NULL)
    return false;
  // Find intersection of track with inner radius of BCAL
  DVector3 proj_pos;
  if (rt->GetIntersectionWithRadius(65.0,proj_pos)==NOERROR){
    double phi=180./M_PI*proj_pos.Phi();
    if (phi<0)phi+=360.;
    double slice=phi/7.5;
    double mid_slice=round(slice);
    module=int(mid_slice)+1;
    sector=int(floor((phi-7.5*mid_slice+3.75)/1.875))+1;
    
    if (intersection) *intersection=proj_pos;
    
    return true;
  }
  
  return false;
}
				   


// Given a reference trajectory from a track, predict which TOF paddles should
// fire due to the charged particle passing through the TOF planes.
bool DParticleID::PredictTOFPaddles(const DReferenceTrajectory *rt,
				    unsigned int &hbar,unsigned int &vbar,
				    DVector3 *intersection) const{
  // Initialize output variables
  vbar=0;
  hbar=0;
  if(rt == NULL)
    return false;
  // Find intersection with TOF plane given by tof_pos
  DVector3 tof_pos(0,0,dTOFGeometry->CenterMPlane);
  DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
  DVector3 proj_mom,proj_pos;
  if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, 
				  NULL,NULL,NULL,SYS_TOF) != NOERROR)
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
	locTimeVariance = 0.5*0.5+locBestMatchParams.dFlightTimeVariance;

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
	double locFlightTimeVariance=9.9E9;
	if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, &locFlightTimeVariance,SYS_FCAL) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locFCALShower->getTime() - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	// Find minimum distance between track projection and each of the hits
	// associated with the shower.
	double d2min=100000.;
	double xproj=proj_pos.x();
	double yproj=proj_pos.y();
	vector<const DFCALCluster*>clusters;
	locFCALShower->Get(clusters);
	for (unsigned int k=0;k<clusters.size();k++){
	  vector<DFCALCluster::DFCALClusterHit_t>hits=clusters[k]->GetHits();
	  for (unsigned int m=0;m<hits.size();m++){
	    double dx=hits[m].x-xproj;
	    double dy=hits[m].y-yproj;
	    double d2=dx*dx+dy*dy;
	    if (d2<d2min){
	      d2min=d2;
	    }
	  }
	}

	double d = sqrt(d2min);
	//double d = (fcal_pos - proj_pos).Mag();
	double p=proj_mom.Mag();
	double cut=FCAL_CUT_PAR1+FCAL_CUT_PAR2/p;
	if(d >= cut)
		return false;

	locShowerMatchParams.dTrack = locTrack;
	locShowerMatchParams.dFCALShower = locFCALShower;
	locShowerMatchParams.dx = 45.0*p/(proj_mom.Dot(norm));
	locShowerMatchParams.dFlightTime = locFlightTime;
	locShowerMatchParams.dFlightTimeVariance = locFlightTimeVariance;
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
		if(MatchToSC(NULL, rt, locSCHits[loc_i], locStartTime, locSCHitMatchParams, locIsTimeBased))
			locSCHitMatchParamsVector.push_back(locSCHitMatchParams);
	}
	if(locSCHitMatchParamsVector.empty())
		return false;

	DSCHitMatchParams locBestMatchParams;
	Get_BestSCMatchParams(locSCHitMatchParamsVector, locBestMatchParams);

	locStartTime = locBestMatchParams.dHitTime - locBestMatchParams.dFlightTime;
	locTimeVariance = locBestMatchParams.dFlightTimeVariance + locBestMatchParams.dHitTimeVariance; 
	//locTimeVariance = 0.3*0.3+locBestMatchParams.dFlightTimeVariance;

	return true;
}

bool DParticleID::MatchToSC(const DKinematicData* locTrack, const DReferenceTrajectory* rt, const DSCHit* locSCHit, double locInputStartTime, DSCHitMatchParams& locSCHitMatchParams, bool locIsTimeBased, const float* locInputDeltaPhiCut, DVector3 *IntersectionPoint, DVector3 *IntersectionDir) const
{
	// NOTE: locTrack is NULL if calling from track reconstruction!!!
	if(sc_pos.empty() || sc_norm.empty())
		return false;

	// Find intersection with a "barrel" approximation for the start counter
	DVector3 proj_pos(NaN,NaN,NaN), proj_mom(NaN,NaN,NaN);
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	double locFlightTimeVariance = 9.9E9;
	unsigned int sc_index=locSCHit->sector - 1;
	if(rt->GetIntersectionWithPlane(sc_pos[sc_index][0], sc_norm[sc_index][0], proj_pos, proj_mom, &locPathLength, &locFlightTime, &locFlightTimeVariance) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	if(fabs(locSCHit->t - locFlightTime - locInputStartTime) > OUT_OF_TIME_CUT)
		return false;

	// Check that the intersection isn't upstream of the paddle
	double myz = proj_pos.z();
	if (myz < sc_pos[sc_index][0].z() + 1e-4)
		return false;
	
	// Look for a match in phi
	//double phi = dSCphi0 + dSCdphi*(locSCHit->sector - 1);
	double sc_dphi_cut = 0.0;
	if(locInputDeltaPhiCut != NULL)
		sc_dphi_cut = *locInputDeltaPhiCut;
	else
		sc_dphi_cut = (locIsTimeBased) ? SC_DPHI_CUT : SC_DPHI_CUT_WB;

	double L = 0.; // Length along scintillator
	double dphi = 0.;

	// Initialize the normal vector for the SC paddle to the long, unbent region
	DVector3 norm = sc_norm[sc_index][0];

	// Now check to see if the intersection is in the nose region and find the
	// start time
	//double locCorrectedHitTime = locSCHit->t - sc_leg_tcor;
	// double sc_pos0 = sc_pos[sc_index][0].z();	
	// double sc_pos1 = sc_pos[sc_index][1].z();

	// Start Counter geometry in hall coordinates, obtained from xml file 
	double sc_pos_soss = sc_pos[sc_index][0].z();   // Start of straight section
	double sc_pos_eoss = sc_pos[sc_index][1].z();   // End of straight section
	double sc_pos_eobs = sc_pos[sc_index][11].z();  // End of bend section
	double sc_pos_eons = sc_pos[sc_index][12].z();  // End of nose section
	// Grab the time-walk corrected start counter hit time, and the pulse integral
	double locCorrectedHitTime   = locSCHit->t;
	double locCorrectedHitEnergy = locSCHit->dE;

	//if(myz <= sc_pos1)

	// Check to see if hit occured in the straight section
	if (myz <= sc_pos_eoss)
	{
		//L=myz;
		//locCorrectedHitTime -= L/C_EFFECTIVE;

		// Apply a user-specified matching cut in the leg region
		DVector3 sc_pos_at_myz = sc_pos[sc_index][0] + (myz - sc_pos[sc_index][0].z())*sc_dir[sc_index][0];
		dphi = sc_pos_at_myz.Phi() - proj_pos.Phi(); //phi could be 0 degrees & proj_phi could be 359 degrees
		while(dphi > TMath::Pi())
			dphi -= M_TWO_PI;
		while(dphi < -1.0*TMath::Pi())
			dphi += M_TWO_PI;
		if(fabs(dphi) > sc_dphi_cut)
			return false; //no match

		// Calculate hit distance along scintillator relative to upstream end
		L = myz - sc_pos_soss;
		// Apply propagation time correction
		locCorrectedHitTime -= L*sc_pt_slope[SC_STRAIGHT][sc_index] + sc_pt_yint[SC_STRAIGHT][sc_index];
		// Apply attenuation correction
		locCorrectedHitEnergy *= 1.0/(exp(sc_attn_B[SC_STRAIGHT_ATTN][sc_index]*L));
	}
	else
	{
		unsigned int num = sc_norm[sc_index].size() - 1;
		for (unsigned int loc_i = 1; loc_i < num; ++loc_i)
		{
			locPathLength = 9.9E9;
			locFlightTime = 9.9E9;
			if(rt->GetIntersectionWithPlane(sc_pos[sc_index][loc_i], sc_norm[sc_index][loc_i], proj_pos, proj_mom, &locPathLength, &locFlightTime, &locFlightTimeVariance) != NOERROR)
				continue;

			myz = proj_pos.z();
			norm = sc_norm[sc_index][loc_i];
			if(myz > sc_pos[sc_index][loc_i + 1].z())
				continue;

			DVector3 sc_pos_at_myz = sc_pos[sc_index][loc_i] + (myz - sc_pos[sc_index][loc_i].z())*sc_dir[sc_index][loc_i];

			dphi = sc_pos_at_myz.Phi() - proj_pos.Phi();
			while(dphi > TMath::Pi())
				dphi -= M_TWO_PI;
			while(dphi < -1.0*TMath::Pi())
				dphi += M_TWO_PI;

			// Open up the phi cut in the nose region
			if(locInputDeltaPhiCut == NULL)
				sc_dphi_cut += SC_DPHI_CUT_SLOPE*(myz - sc_pos_eoss);

			if (fabs(dphi) > sc_dphi_cut)
				return false;
		}

		// Check for intersection point beyond nose
		if (myz > sc_pos[sc_index][num].z())
			return false;
	  
		// Note: in the following code, L does not include a correction for where the start counter starts in z ...
		// This is absorbed into locCorrectedHitTime, above.
		// if(myz < sc_pos1)

		// If location cannot be located in straight, bend, or nose section then
		// force it to be located at the end/start of the straight/bend section
		if(myz < sc_pos_eoss)  // Assume hit at end of straight section
	    {
			// L = sc_pos1;
			// locCorrectedHitTime -= L/C_EFFECTIVE;

			// Define the distance to be at the end of the end/start of the straight/bend section
			L = sc_pos_eoss - sc_pos_soss;
			// Apply propagation time correction
			locCorrectedHitTime -= L*sc_pt_slope[SC_STRAIGHT][sc_index] + sc_pt_yint[SC_STRAIGHT][sc_index];
			// Apply attenuation correction
			locCorrectedHitEnergy *= 1.0/(exp(sc_attn_B[SC_STRAIGHT_ATTN][sc_index]*L));
		}
		//else
		//{
		//	L = (myz - sc_pos1)*sc_angle_cor + sc_pos1;
		//	locCorrectedHitTime -= L/C_EFFECTIVE;
		//}
	  
		// Check to see if hit occured in bend section and apply corrections
		else if(myz > sc_pos_eoss && myz <= sc_pos_eobs)
		{
			//L = (myz - sc_pos_eoss)*sc_angle_corr);
			//locCorrectedHitTime -= L/C_EFFECTIVE;

			// Calculate the hit position relative to the upstream end
			L = (myz - sc_pos_eoss)*sc_angle_cor + (sc_pos_eoss - sc_pos_soss);
			// Apply propagation time correction
			locCorrectedHitTime -= L*sc_pt_slope[SC_BEND][sc_index] + sc_pt_yint[SC_BEND][sc_index];
			// Apply attenuation correction
			locCorrectedHitEnergy *= (sc_attn_A[SC_STRAIGHT_ATTN][sc_index] / ((sc_attn_A[SC_BENDNOSE_ATTN][sc_index]*
					exp(sc_attn_B[SC_BENDNOSE_ATTN][sc_index]*L)) + sc_attn_C[SC_BENDNOSE_ATTN][sc_index]));
		}
		// Check to see if hit occured in nose section and apply corrections
		else if(myz > sc_pos_eobs && myz <= sc_pos_eons)
		{
			//L = (myz - sc_pos_eoss)*sc_angle_cor;

			// Calculate the hit position relative to the upstream end
			L = (myz - sc_pos_eoss)*sc_angle_cor + (sc_pos_eoss - sc_pos_soss);
			// Apply propagation time correction
			locCorrectedHitTime -= L*sc_pt_slope[SC_NOSE][sc_index] + sc_pt_yint[SC_NOSE][sc_index];
			// Apply attenuation correction
			locCorrectedHitEnergy *= (sc_attn_A[SC_STRAIGHT_ATTN][sc_index] / ((sc_attn_A[SC_BENDNOSE_ATTN][sc_index]*
					exp(sc_attn_B[SC_BENDNOSE_ATTN][sc_index]*L)) + sc_attn_C[SC_BENDNOSE_ATTN][sc_index]));
		}
	}

	double ds = 0.3*proj_mom.Mag()/fabs(proj_mom.Dot(norm));

	// For the dEdx measurement we now need to take into account that L does not 
	// compensate for the position in z at which the start counter paddle starts
	locSCHitMatchParams.dTrack = locTrack;
	locSCHitMatchParams.dSCHit = locSCHit;
	//locSCHitMatchParams.dHitEnergy = (locSCHit->dE)*exp((L - sc_pos0)/ATTEN_LENGTH);
	locSCHitMatchParams.dHitEnergy = locCorrectedHitEnergy;
	locSCHitMatchParams.dEdx = locSCHitMatchParams.dHitEnergy/ds;
	locSCHitMatchParams.dHitTime = locCorrectedHitTime;
	locSCHitMatchParams.dHitTimeVariance = sc_paddle_resols[sc_index]*sc_paddle_resols[sc_index];
	locSCHitMatchParams.dFlightTime = locFlightTime;
	locSCHitMatchParams.dFlightTimeVariance = locFlightTimeVariance;
	locSCHitMatchParams.dPathLength = locPathLength;
	locSCHitMatchParams.dDeltaPhiToHit = dphi;

	// Optionally output intersection position	
	if (IntersectionPoint!=NULL)
		*IntersectionPoint=proj_pos;
	if (IntersectionDir!=NULL)
	{
		*IntersectionDir=proj_mom;
		IntersectionDir->SetMag(1.);
	}

	return true;
}

// Predict the start counter paddle that would match a track whose reference 
// trajectory is given by rt.
unsigned int DParticleID::PredictSCSector(const DReferenceTrajectory* rt, const double dphi_cut, DVector3* locProjPos, bool* locProjBarrelRegion, double* locMinDPhi) const
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
    if (myz<sc_pos[sc_index][0].z()+1e-4) continue;

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
      if (fabs(dphi)<fabs(min_dphi)){
	best_sc_index=sc_index;
   if(locProjBarrelRegion != NULL)
     *locProjBarrelRegion = true;
   if(locProjPos != NULL)
     *locProjPos = proj_pos;
	min_dphi=dphi;
      }
    }
    else{
      bool got_match=false;
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
	  if (fabs(dphi)<SC_DPHI_CUT_WB){
	    got_match=true;
	    break;
	  }
	}
      }	
      if (got_match==false) continue; //doesn't match well with nose region

      // Check for intersection point beyond nose
      if (myz> sc_pos[sc_index][num].z()) continue;

      	

      if (fabs(dphi)<fabs(min_dphi)){
	best_sc_index=sc_index;
   if(locProjBarrelRegion != NULL)
     *locProjBarrelRegion = false;
   if(locProjPos != NULL)
     *locProjPos = proj_pos;
	min_dphi=dphi;
      }
    }
  }

  if(locMinDPhi != NULL)
    *locMinDPhi = min_dphi;

  if (fabs(min_dphi)<dphi_cut) return best_sc_index+1;

  return 0;
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
	locDistance = rt->DistToRTwithTime(bcal_pos, &locPathLength, &locFlightTime,NULL,SYS_BCAL);
	if(!isfinite(locDistance))
		return false;

	// Check that the hit is not out of time with respect to the track
	double locDeltaT = locBCALShower->t - locFlightTime - locInputStartTime;
	if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
		return false;

	DVector3 proj_pos = rt->GetLastDOCAPoint();
	if(proj_pos.Perp() < 65.0)
		return false;  // not inside BCAL!

	locDeltaZ = bcal_pos.z() - proj_pos.z();
	locDeltaPhi = bcal_pos.Phi() - proj_pos.Phi();
	while(locDeltaPhi > M_PI)
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
	if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime,NULL,SYS_FCAL) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	double locDeltaT = locFCALShower->getTime() - locFlightTime - locInputStartTime;
	if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
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
	if(rt->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, NULL, SYS_TOF) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	double locDeltaT = locTOFPoint->t - locFlightTime - locInputStartTime;
	if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
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

	// Find intersection with a "barrel" approximation for the start counter
	DVector3 proj_pos(NaN,NaN,NaN), proj_mom(NaN,NaN,NaN);
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(rt->GetIntersectionWithRadius(sc_pos[0][1].x(), proj_pos, &locPathLength, &locFlightTime, &proj_mom) != NOERROR)
		return false;

	// Check that the hit is not out of time with respect to the track
	double locDeltaT = locSCHit->t - locFlightTime - locInputStartTime;
	if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
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

const DTOFPoint* DParticleID::Get_ClosestToTrack_TOFPoint(const DKinematicData* locTrack, vector<const DTOFPoint*>& locTOFPoints, double& locBestDeltaX, double& locBestDeltaY) const
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

const DTOFPaddleHit* DParticleID::Get_ClosestTOFPaddleHit_Horizontal(const DReferenceTrajectory* locReferenceTrajectory, const vector<const DTOFPaddleHit*>& locTOFPaddleHits, double locInputStartTime, double& locBestDeltaY) const
{
	if(locReferenceTrajectory == nullptr)
		return nullptr;

	// Evaluate matching solely by physical geometry of the paddle: NOT the distance along the paddle of the hit
	DVector3 tof_pos(0.0, 0.0, dTOFGeometry->CenterHPlane); //a point on the TOF plane
	DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(locReferenceTrajectory->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, NULL, SYS_TOF) != NOERROR)
		return nullptr;

	const DTOFPaddleHit* locClosestPaddleHit = nullptr;
	locBestDeltaY = 999.0;
	for(auto& locTOFPaddleHit : locTOFPaddleHits)
	{
		if(locTOFPaddleHit->orientation != 1)
			continue; //horizontal orientation is 1

		bool locNorthIsGoodHitFlag = (locTOFPaddleHit->E_north > TOF_E_THRESHOLD);
		bool locSouthIsGoodHitFlag = (locTOFPaddleHit->E_south > TOF_E_THRESHOLD);
		if(!locNorthIsGoodHitFlag && !locSouthIsGoodHitFlag)
			continue; //hit is junk

		// Check geometric distance, if better than before
		double locDeltaY = dTOFGeometry->bar2y(locTOFPaddleHit->bar) - proj_pos.Y();
		if(fabs(locDeltaY) > fabs(locBestDeltaY))
			continue;

		// Check that the hit is not out of time with respect to the track

		//Construct spacetime hit: averages times, or if only one end with hit, reports time at center of paddle
		DTOFPoint_factory::tof_spacetimehit_t* locSpacetimeHit = dTOFPointFactory->Build_TOFSpacetimeHit_Horizontal(locTOFPaddleHit);
		double locHitTime = locSpacetimeHit->t;

		// if single-ended paddle, or only one side has a hit: time reported at center: must propagate to track location
		if(locNorthIsGoodHitFlag != locSouthIsGoodHitFlag)
		{
			//Paddle midpoint
			double locPaddleMidPoint = 0.0; //is 0 except when is single-ended bar (22 & 23)
			if(!locSpacetimeHit->dIsDoubleEndedBar)
				locPaddleMidPoint = locNorthIsGoodHitFlag ? ONESIDED_PADDLE_MIDPOINT_MAG : -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;

			//correct the time
			double locDistanceToMidPoint = locNorthIsGoodHitFlag ? locPaddleMidPoint - proj_pos.X() : proj_pos.X() - locPaddleMidPoint;
			int id = 44 + locTOFPaddleHit->bar - 1; //for propation speed
			locHitTime -= locDistanceToMidPoint/propagation_speed[id];
		}

		//time cut
		double locDeltaT = locHitTime - locFlightTime - locInputStartTime;
		if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
			continue;

		locBestDeltaY = locDeltaY;
		locClosestPaddleHit = locTOFPaddleHit;
	}

	return locClosestPaddleHit;
}

const DTOFPaddleHit* DParticleID::Get_ClosestTOFPaddleHit_Vertical(const DReferenceTrajectory* locReferenceTrajectory, const vector<const DTOFPaddleHit*>& locTOFPaddleHits, double locInputStartTime, double& locBestDeltaX) const
{
	if(locReferenceTrajectory == nullptr)
		return nullptr;

	// Evaluate matching solely by physical geometry of the paddle: NOT the distance along the paddle of the hit
	DVector3 tof_pos(0.0, 0.0, dTOFGeometry->CenterVPlane); //a point on the TOF plane
	DVector3 norm(0.0, 0.0, 1.0); //normal vector to TOF plane
	DVector3 proj_pos, proj_mom;
	double locPathLength = 9.9E9, locFlightTime = 9.9E9;
	if(locReferenceTrajectory->GetIntersectionWithPlane(tof_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, NULL, SYS_TOF) != NOERROR)
		return nullptr;

	const DTOFPaddleHit* locClosestPaddleHit = nullptr;
	locBestDeltaX = 999.0;
	for(auto& locTOFPaddleHit : locTOFPaddleHits)
	{
		if(locTOFPaddleHit->orientation != 0)
			continue; //vertical orientation is 0

		bool locNorthIsGoodHitFlag = (locTOFPaddleHit->E_north > TOF_E_THRESHOLD);
		bool locSouthIsGoodHitFlag = (locTOFPaddleHit->E_south > TOF_E_THRESHOLD);
		if(!locNorthIsGoodHitFlag && !locSouthIsGoodHitFlag)
			continue; //hit is junk

		// Check geometric distance, if better than before
		double locDeltaX = dTOFGeometry->bar2y(locTOFPaddleHit->bar) - proj_pos.X();
		if(fabs(locDeltaX) > fabs(locBestDeltaX))
			continue;

		// Check that the hit is not out of time with respect to the track

		//Construct spacetime hit: averages times, or if only one end with hit, reports time at center of paddle
		DTOFPoint_factory::tof_spacetimehit_t* locSpacetimeHit = dTOFPointFactory->Build_TOFSpacetimeHit_Vertical(locTOFPaddleHit);
		double locHitTime = locSpacetimeHit->t;

		// if single-ended paddle, or only one side has a hit: time reported at center: must propagate to track location
		if(locNorthIsGoodHitFlag != locSouthIsGoodHitFlag)
		{
			//Paddle midpoint
			double locPaddleMidPoint = 0.0; //is 0 except when is single-ended bar (22 & 23)
			if(!locSpacetimeHit->dIsDoubleEndedBar)
				locPaddleMidPoint = locNorthIsGoodHitFlag ? ONESIDED_PADDLE_MIDPOINT_MAG : -1.0*ONESIDED_PADDLE_MIDPOINT_MAG;

			//correct the time
			double locDistanceToMidPoint = locNorthIsGoodHitFlag ? locPaddleMidPoint - proj_pos.Y() : proj_pos.Y() - locPaddleMidPoint;
			int id = locTOFPaddleHit->bar - 1; //for propation speed
			locHitTime -= locDistanceToMidPoint/propagation_speed[id];
		}

		//time cut
		double locDeltaT = locHitTime - locFlightTime - locInputStartTime;
		if(fabs(locDeltaT) > OUT_OF_TIME_CUT)
			continue;

		locBestDeltaX = locDeltaX;
		locClosestPaddleHit = locTOFPaddleHit;
	}

	return locClosestPaddleHit;
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
		if(fabs(locDeltaPhi) > fabs(locBestDeltaPhi))
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

void DParticleID::Get_CDCNumHitRingsPerSuperlayer(int locBitPattern, map<int, int>& locNumHitRingsPerSuperlayer) const
{
	set<int> locCDCRings;
	Get_CDCRings(locBitPattern, locCDCRings);
	Get_CDCNumHitRingsPerSuperlayer(locCDCRings, locNumHitRingsPerSuperlayer);
}

void DParticleID::Get_CDCNumHitRingsPerSuperlayer(const set<int>& locCDCRings, map<int, int>& locNumHitRingsPerSuperlayer) const
{
	locNumHitRingsPerSuperlayer.clear();
	for(int locCDCSuperlayer = 1; locCDCSuperlayer <= 7; ++locCDCSuperlayer)
		locNumHitRingsPerSuperlayer[locCDCSuperlayer] = 0;

	set<int>::const_iterator locIterator = locCDCRings.begin();
	for(; locIterator != locCDCRings.end(); ++locIterator)
	{
		int locCDCSuperlayer = ((*locIterator) - 1)/4 + 1;
		++locNumHitRingsPerSuperlayer[locCDCSuperlayer];
	}
}

void DParticleID::Get_FDCNumHitPlanesPerPackage(int locBitPattern, map<int, int>& locNumHitPlanesPerPackage) const
{
	set<int> locFDCPlanes;
	Get_FDCPlanes(locBitPattern, locFDCPlanes);
	Get_FDCNumHitPlanesPerPackage(locFDCPlanes, locNumHitPlanesPerPackage);
}

void DParticleID::Get_FDCNumHitPlanesPerPackage(const set<int>& locFDCPlanes, map<int, int>& locNumHitPlanesPerPackage) const
{
	locNumHitPlanesPerPackage.clear();
	for(int locFDCPackage = 1; locFDCPackage <= 4; ++locFDCPackage)
		locNumHitPlanesPerPackage[locFDCPackage] = 0;

	set<int>::const_iterator locIterator = locFDCPlanes.begin();
	for(; locIterator != locFDCPlanes.end(); ++locIterator)
	{
		int locFDCPackage = ((*locIterator) - 1)/6 + 1;
		map<int, int>::iterator locMapIterator = locNumHitPlanesPerPackage.find(locFDCPackage);
		++locNumHitPlanesPerPackage[locFDCPackage];
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

