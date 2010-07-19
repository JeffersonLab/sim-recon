//************************************************************************
// DTrackFitterRiemann.cc
//************************************************************************
// Uses drift time information to refine the results of the circle and line
// fits in the Riemann fit formalism.

#include "DTrackFitterRiemann.h"
#include "START_COUNTER/DSCHit.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include <TDecompLU.h>
#include <math.h>

#define MOLIERE_FRACTION 0.99
#define ONE_THIRD  0.33333333333333333
#define ONE_SIXTH  0.16666666666666667
#define TWO_THIRDS 0.66666666666666667
#define EPS 1e-8
#define Z_MIN 50.0
#define Z_VERTEX 65.0
#define Z_MAX 80.0

DTrackFitterRiemann::DTrackFitterRiemann(JEventLoop *loop):DTrackFitter(loop){

}

DTrackFitter::fit_status_t DTrackFitterRiemann::FitTrack(void)
{
  // Use the current best knowledge for the track parameters at the "vertex"
  // to set the seed (initial) values for the fit
  jerror_t error = SetSeed(input_params.charge(), input_params.position(), 
                           input_params.momentum(),input_params.mass());
  if (error!=NOERROR) return kFitFailed;

  // Clear the hits
  for (unsigned int i=0;i<my_line_hits.size();i++){
    if (my_line_hits[i]->fdc==NULL)
      delete my_line_hits[i];
  }
  my_line_hits.clear(); 
  for (unsigned int i=0;i<my_circle_hits.size();i++){
    delete my_circle_hits[i];
  }
  my_circle_hits.clear();

  // Clear the projection and s vectors
  projections.clear();
  s.clear();

  // Deal with cdc axial hits
  DVector2 XY(-D*sin(phi0),D*cos(phi0)); // Starting radial coords.
  double sperp=0.; // perpendicular arc length
  for (unsigned int i=0;i<cdchits.size();i++){
    // Axial wires
    if (fabs(cdchits[i]->wire->stereo)<EPS){ 
      DRiemannHit_t *hit= new DRiemannHit_t;
      // Pointers to fdc/cdc hit objects
      hit->fdc=NULL;
      hit->cdc=cdchits[i];

      GetAxialPosition(sperp,XY,hit);
      XY=hit->XY;

      my_circle_hits.push_back(hit);
    }
  }

  // Add the FDC hits to the circle hit list
  for (unsigned int i=0;i<fdchits.size();i++){
    DRiemannHit_t *hit= new DRiemannHit_t;
    // Pointers to fdc/cdc hit objects
    hit->fdc=fdchits[i];
    hit->cdc=NULL;

    hit->z=hit->fdc->wire->origin.z();
    GetFDCPosition(hit);

    // Measurement covariance
    double cosa=hit->fdc->wire->udir.y();
    double sina=hit->fdc->wire->udir.x();
    double cosa2=cosa*cosa;
    double sina2=sina*sina;
    double sigx2=0.015*0.015;
    double sigy2=0.02*0.02;
    hit->covx=sigx2*cosa2+sigy2*sina2;
    hit->covy=sigx2*sina2+sigy2*cosa2;
    hit->covxy=(sigy2-sigx2)*sina*cosa;
    
    my_circle_hits.push_back(hit);
  }

  // Check that we have enough hits on the circle to proceed
  if (my_circle_hits.size()<3) return kFitFailed;

  // Otherwise proceed with the fit...
    
  // Using the new combined list of hits, compute the covariance matrix 
  // for RPhi associated with these hits
  unsigned int ncirclehits=my_circle_hits.size();
  CRPhi.ResizeTo(ncirclehits,ncirclehits);
  ComputeCRPhi();    
  
  // Do the circle fit
  if (FitCircle()!=NOERROR) return kFitFailed;
 
  // Deal with cdc stereo hits
  // .. First compute starting radial coords
  XY.Set(-D*sin(phi0),D*cos(phi0));
  sperp=0.; // perpendicular arc length
  for (unsigned int i=0;i<cdchits.size();i++){
    // Stereo wires
    if (fabs(cdchits[i]->wire->stereo)>EPS){ 
      DRiemannHit_t *hit= new DRiemannHit_t;
      // Pointers to fdc/cdc hit objects
      hit->fdc=NULL;
      hit->cdc=cdchits[i];
      
      GetStereoPosition(sperp,XY,hit);
 
      my_line_hits.push_back(hit);
    }
  }
  
  // Add the fdc hits to the my_line_hits vector
  for (unsigned int i=0;i<my_circle_hits.size();i++){
    DRiemannHit_t *hit= my_circle_hits[i];
    if (hit->fdc!=NULL){
      GetFDCPosition(hit);

      my_line_hits.push_back(hit);
    }
  }
    
  // Check that we have enough hits for the line fit
  if (my_line_hits.size()<2) return kFitFailed;

  // If we have enough hits, proceed...
  unsigned int nhits=my_line_hits.size();
  projections.resize(nhits);
  s.resize(nhits);
    
  if (ComputeIntersections()!=NOERROR) return kFitFailed;
   
  // If there are no fdc hits, then the data needed for the line fit is coming
  // purely from the stereo wires, for which the relevant error is in z 
  if (fdchits.size()==0){
    unsigned int nstereo=my_line_hits.size();
    Cz.ResizeTo(nstereo,nstereo);
    ComputeCz();
  }
  else{
    // If we have FDC hits, then the relevant error for the line fit is in R:
    // Size the matrices according to the number of hits
    CR.ResizeTo(nhits,nhits);
    ComputeCR();
  }

  // Do the line fit
  FitLine();
 
  if (cdchits.size()>0){
    // Get the field value
    B=bfield->GetBz(-D*sin(phi0),D*cos(phi0),z_vertex);
    
    // Compute the magnitude of the momentum at this stage of the fit:
    p=0.003*fabs(B)*rc/cos(atan(tanl));
    one_over_v=sqrt(1.+mass2/(p*p))/29.98;
    
    // reset sperp to zero
    sperp=0.; 
    // Get FDC and CDC axial positions with refined circle and line parameters
    for (unsigned int i=0;i<my_circle_hits.size();i++){
      DRiemannHit_t *hit=my_circle_hits[i];
      if (hit->fdc!=NULL){
	GetFDCPosition(hit);
      }
      if (hit->cdc!=NULL){
	GetAxialPosition(sperp,XY,hit);
	XY=hit->XY;
      }
    }   
    // Compute the covariance matrix for RPhi
    ComputeCRPhi();    
    
    // Do the circle fit
    if (FitCircle()!=NOERROR) return kFitFailed;
    
    // Recompute the intersections with the circle and redo the line fit, 
    // this time with drift time information for the stereo wires
    double cosl=cos(atan(tanl));
    //double s_old=0;
    for (unsigned int i=0;i<nhits;i++){
      DRiemannHit_t *hit=my_line_hits[i];
      // Use the results of the previous line fit to predict the z-value 
      // for a certain arc length value s and use this to determine the 
      // direction of the correction to the wire position for the drift
      // time.
      if (hit->cdc!=NULL){
	DVector2 XYp=GetHelicalPosition(s[i]);
	DVector2 dXY(XYp.X()-xc,XYp.Y()-yc);
	double drw=dXY.Mod();
	DVector2 dir=(1./drw)*dXY;
	double tdrift=hit->cdc->tdrift-s[i]*one_over_v/cosl;
	// prediction for z
	double zpred=z_vertex+s[i]*tanl;
	DVector2 dXYtest=0.0055*tdrift*dir;
	// Two LR solutions 
	double zplus=GetStereoZ(dXYtest.X(),dXYtest.Y(),hit);
	double zminus=GetStereoZ(-dXYtest.X(),-dXYtest.Y(),hit);
	if (zpred>=167.3){
	  zpred=167.3;
	  if (zminus<zpred) hit->z=zminus;
	  else if (zplus<zpred) hit->z=zplus;
	  else hit->z=zpred;
	}
	else if (zpred<=17.){
	  zpred=17.0;
	  if (zminus>zpred) hit->z=zminus;
	  else if (zplus>zpred) hit->z=zplus;
	  else hit->z=zpred;
	}
	else{
	  if (fabs(zplus-zpred)<fabs(zminus-zpred)) hit->z=zplus;
	  else hit->z=zminus;
	}   
	// Crude approximation for covariances
	hit->covx=0.015*0.015*dir.X()*dir.X();
	hit->covy=0.015*0.015*dir.Y()*dir.Y();
	hit->covz*=(0.0015*0.0015)/0.2133;
      }
      else{
	GetFDCPosition(hit);
      }
    }
    if (ComputeIntersections()!=NOERROR) return kFitFailed;
    
    // New covariance matrix for z
    if (fdchits.size()==0){
      // New covariance matrix for z
      ComputeCz();
    }
    else{
    // New covariance matrix for R
      ComputeCR();
    }

    // Do the line fit
    FitLine();
  }

  double sinphi=sin(phi0);
  double cosphi=cos(phi0);
  double x0=-D*sinphi;
  double y0=D*cosphi;
  if (fit_type==kWireBased){
    B=bfield->GetBz(x0,y0,z_vertex);
  }
  else{
    B=0;
    for (unsigned int i=0;i<my_line_hits.size();i++){
      const DRiemannHit_t *hit=my_line_hits[i];
      B+=bfield->GetBz(hit->XY.X(),hit->XY.Y(),hit->z);
    }
    B/=my_line_hits.size();
  }
  
  double pt=0.003*fabs(B)*rc;
  fit_params.setPosition(DVector3(x0,y0,z_vertex));
  fit_params.setMomentum(DVector3(pt*cosphi,pt*sinphi,pt*tanl));
  fit_params.setCharge(q);
  this->chisq=ChiSq();
  
  return kFitSuccess;
}

//-----------------
// ChiSq
//-----------------
double DTrackFitterRiemann::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr, vector<pull_t> *pulls_ptr)
{
  double chisq = ChiSq();
  unsigned int ndf = this->Ndof;
  
  if(chisq_ptr)*chisq_ptr = chisq;
  if(dof_ptr)*dof_ptr = int(ndf);
  //if(pulls_ptr)*pulls_ptr = pulls;
  
  printf("reduced chi2 %f\n",chisq/double(ndf));
  return chisq/double(ndf);


}

// Compute the chi2 for the fit using both the line and circle hits
double DTrackFitterRiemann::ChiSq(){
  double chi2=0;
  this->Ndof=my_circle_hits.size()-5;
  double sperp=0.;
  double x0=-D*sin(phi0);
  double y0=D*cos(phi0);
  DVector2 XYold(x0,y0);
  // First take care of the FDC hits and the cdc axial wires
  for (unsigned int i=0;i<my_circle_hits.size();i++){
    if (my_circle_hits[i]->fdc!=NULL){
      sperp=(my_circle_hits[i]->z-z_vertex)/tanl;
    }
    else{
      double ratio=(my_circle_hits[i]->XY-XYold).Mod()/(2.*rc);
      sperp+=2.*rc*((ratio>1.)?M_PI_2:asin(ratio));
      XYold=my_circle_hits[i]->XY;
    }
    // Position on the fitted helix
    DVector2 XYp=GetHelicalPosition(sperp);

    double Phi=my_circle_hits[i]->XY.Phi();
    double cosPhi=cos(Phi);
    double sinPhi=sin(Phi);
    chi2+=(my_circle_hits[i]->XY-XYp).Mod2()
      /(cosPhi*cosPhi*my_circle_hits[i]->covx
	+sinPhi*sinPhi*my_circle_hits[i]->covy
	+2.*sinPhi*cosPhi*my_circle_hits[i]->covxy);
  }

  // Next take care of the cdc stereo wires
  XYold.Set(x0,y0);
  sperp=0.;
  for (unsigned int i=0;i<my_line_hits.size();i++){
    if (my_line_hits[i]->cdc!=NULL){
      this->Ndof++;
      //double sperp=(my_line_hits[i]->z-z_vertex)/tanl;
      double ratio=(my_line_hits[i]->XY-XYold).Mod()/(2.*rc);
      sperp+=2.*rc*((ratio>1.)?M_PI_2:asin(ratio));
      XYold=my_line_hits[i]->XY;

      // Position on the fitted helix
      DVector2 XYp=GetHelicalPosition(sperp);  
   
      double Phi=my_line_hits[i]->XY.Phi();
      double cosPhi=cos(Phi);
      double sinPhi=sin(Phi);
      chi2+=(my_line_hits[i]->XY-XYp).Mod2()
	/(cosPhi*cosPhi*my_line_hits[i]->covx
	  +sinPhi*sinPhi*my_line_hits[i]->covy
	  +2.*sinPhi*cosPhi*my_line_hits[i]->covxy);
    }
  }

  return chi2;
}

// Initialize the state vector
jerror_t DTrackFitterRiemann::SetSeed(double my_q,const DVector3 &pos,
				      const DVector3 &mom,
				      double mass){
  if (!isfinite(pos.Mag()) || !isfinite(mom.Mag())){
    _DBG_ << "Invalid seed data." <<endl;
    return UNRECOVERABLE_ERROR;
  }
  // mass squared
  mass2=mass*mass;

  // Momentum
  p=mom.Mag();
  if (p>8.){
    p=8.;
  }
  one_over_v=sqrt(1.+mass2/(p*p))/29.98;

  // Charge
  q=my_q;
  
  // Dip angle
  //double lambda=M_PI_2-mom.Theta();;
  //tanl=tan(lambda); 
  theta=mom.Theta();
  tanl=tan(M_PI_2-theta);
 
  // Phi angle
  phi0=mom.Phi();

  // "vertex" position along z
  z_vertex=pos.z();

  // Circle parameters
  B=bfield->GetBz(pos.x(),pos.y(),pos.z());
  if (fit_type==kTimeBased){
    if (fdchits.size()>0){
      for (unsigned int i=0;i<fdchits.size();i++){
	const DFDCPseudo *hit=fdchits[i];
	double u=hit->w;
	double v=hit->s; 
	double cosa=hit->wire->udir.y();
	double sina=hit->wire->udir.x();
	B+=bfield->GetBz(u*cosa+v*sina,-u*sina+v*cosa,hit->wire->origin.z());
      }
      B/=fdchits.size()+1;
    }
  }

  rc=p*sin(theta)/fabs(0.003*B);
  xc=pos.x()-q*rc*sin(phi0);
  yc=pos.y()+q*rc*cos(phi0);
  // Signed distance to origin
  D=-pos.x()/sin(phi0);

  return NOERROR;
}

// Correct the axial wire positions with the drift distance
jerror_t DTrackFitterRiemann::GetAxialPosition(double &sperp,
					       const DVector2 &XYold,
					       DRiemannHit_t *hit){
  DVector2 XY(hit->cdc->wire->origin.x(),hit->cdc->wire->origin.y());
  DVector2 dXY(XY.X()-xc,XY.Y()-yc);
  double drw=dXY.Mod();
  DVector2 dir=(1./drw)*dXY;

  double sign=(drw<rc)?1.:-1.;
  double ratio=(XY-XYold).Mod()/(2.*rc);
  double cosl=cos(atan(tanl));
  sperp+=2.*rc*((ratio>1.)?M_PI_2:asin(ratio));
  double tflight=sperp*one_over_v/cosl;
  double tdrift=hit->cdc->tdrift-tflight;
  hit->XY=XY+sign*0.0055*tdrift*dir;

  // Crude approximation for covariance matrix ignoring error in dir
  hit->covx=0.015*0.015*dir.X()*dir.X();
  hit->covy=0.015*0.015*dir.Y()*dir.Y();
  hit->covxy=0.;

  // Guess z-position from result of candidate fit
  hit->z=z_vertex+tanl*sperp;
  if (hit->z>167.3){
    hit->z=167.3;
    sperp=(hit->z-z_vertex)/tanl;
    hit->XY=XY+sign*0.0055*(hit->cdc->tdrift-sperp*one_over_v/cosl)*dir;
  }
  else if (hit->z<17.0){
    hit->z=17.0;
    sperp=(hit->z-z_vertex)/tanl;
    hit->XY=XY+sign*0.0055*(hit->cdc->tdrift-sperp*one_over_v/cosl)*dir;  
  }
  return NOERROR;
}

// Get the z-position of the intersection between the cylinder of radius rc
// and the stereo wire
double DTrackFitterRiemann::GetStereoZ(double dx,double dy,
				       DRiemannHit_t *hit){
  double uz=hit->cdc->wire->udir.z();
  double ux=hit->cdc->wire->udir.x()/uz;
  double uy=hit->cdc->wire->udir.y()/uz;
  double denom=ux*ux+uy*uy;
  double my_z=0.;
  
  // wire origin
  double xwire0=hit->cdc->wire->origin.x()+dx;
  double ywire0=hit->cdc->wire->origin.y()+dy;
  double zwire0=hit->cdc->wire->origin.z();
  my_z=zwire0+((xc-xwire0)*ux+(yc-ywire0)*uy)/denom;
  double rootz2=denom*rc*rc-(uy*(xwire0-xc)-ux*(ywire0-yc))
    *(uy*(xwire0-xc)-ux*(ywire0-yc));  
  if (rootz2>0.){
    double rootz=sqrt(rootz2);
    double z1=my_z+rootz/denom;
    double z2=my_z-rootz/denom;

    if (fabs(z2-Z_VERTEX)>fabs(z1-Z_VERTEX)){
      my_z=z1;
    }
    else{
      my_z=z2; 
    }
  }
  if (my_z>167.3) my_z=167.3;
  if (my_z<17.0) my_z=17.0;
  
  return my_z;
}

// Get the position of the intersection between the cylinder of radius rc
// and the stereo wire
jerror_t DTrackFitterRiemann::GetStereoPosition(double &sperp,
						DVector2 &XYold,
						DRiemannHit_t *hit){
  double uz=hit->cdc->wire->udir.z();
  double ux=hit->cdc->wire->udir.x()/uz;
  double uy=hit->cdc->wire->udir.y()/uz;
  double denom=ux*ux+uy*uy;
  
  // wire origin
  double xwire0=hit->cdc->wire->origin.x();
  double ywire0=hit->cdc->wire->origin.y();
  double zwire0=hit->cdc->wire->origin.z();
  hit->z=zwire0+((xc-xwire0)*ux+(yc-ywire0)*uy)/denom;
  double rootz2=denom*rc*rc-(uy*(xwire0-xc)-ux*(ywire0-yc))
    *(uy*(xwire0-xc)-ux*(ywire0-yc));  
  double dx0dz=ux/denom;
  double dy0dz=uy/denom;
  if (rootz2>0.){
    double rootz=sqrt(rootz2);
    double z1=hit->z+rootz/denom;
    double z2=hit->z-rootz/denom;

    if (fabs(z2-Z_VERTEX)>fabs(z1-Z_VERTEX)){
      hit->z=z1;
      dx0dz+=((uy*(xwire0-xc)-ux*(ywire0-yc))*uy/rootz)/denom;
      dy0dz-=((uy*(xwire0-xc)-ux*(ywire0-yc))*ux/rootz)/denom;
    }
    else{
      hit->z=z2; 
      dx0dz-=((uy*(xwire0-xc)-ux*(ywire0-yc))*uy/rootz)/denom; 
      dy0dz+=((uy*(xwire0-xc)-ux*(ywire0-yc))*ux/rootz)/denom;
    }
  }
  // Deal with endplates
  if (hit->z>167.3) hit->z=167.3;
  if (hit->z<17.0) hit->z=17.0;
  double dzwire=hit->z-zwire0;
  hit->XY.Set(xwire0+ux*dzwire,ywire0+uy*dzwire);

  DVector2 dXY(hit->XY.X()-xc,hit->XY.Y()-yc);
  double drw=dXY.Mod();
  DVector2 dir=(1./drw)*dXY;
  // Crude approximation for covariance matrix ignoring error in dir
  hit->covx=0.2133*dir.X()*dir.X();
  hit->covy=0.2133*dir.Y()*dir.Y();
  hit->covxy=0.;
  hit->covz=dx0dz*dx0dz*hit->covx+dy0dz*dy0dz*hit->covy;

  XYold=hit->XY;

  return NOERROR;
}

// Compute the position on the helical trajectory for a given arc length
DVector2 DTrackFitterRiemann::GetHelicalPosition(double sperp){
    double twokappa=q/rc;
    double twoks=twokappa*sperp;
    double sin2ks=sin(twoks);
    double cos2ks=cos(twoks);
    double cosphi=cos(phi0);
    double sinphi=sin(phi0);
    double one_minus_cos2ks=1.-cos2ks;
    double one_over_2k=1./(twokappa);
    return 
      DVector2(-D*sinphi+(cosphi*sin2ks-sinphi*one_minus_cos2ks)*one_over_2k,
	       D*cosphi+(sinphi*sin2ks+cosphi*one_minus_cos2ks)*one_over_2k);
}


// Compute the position from the FDC hit, correcting for the drift time and
// the Lorentz effect
jerror_t DTrackFitterRiemann::GetFDCPosition(DRiemannHit_t *hit){          
  // Position on the helical trajectory 
  hit->z=hit->fdc->wire->origin.z();
  double sperp=(hit->z-z_vertex)/tanl;
  //double my_phi=phi1+q*sperp/rc;
  //double x=xc+rc*cos(my_phi);
  //double y=yc+rc*sin(my_phi);
  DVector2 XYp=GetHelicalPosition(sperp);
    
  // Find difference between point on helical path and wire
  double cosa=hit->fdc->wire->udir.y();
  double sina=hit->fdc->wire->udir.x();
  //double w=x*cosa-y*sina-hit->fdc->w;
  double w=XYp.X()*cosa-XYp.Y()*sina-hit->fdc->w;
  // .. and use it to determine which sign to use for the drift time data
  double sign=(w>0?1.:-1.);
  
  // Correct the drift time for the flight path and convert to distance 
  // units
  double delta_x=sign*(hit->fdc->time-sperp*one_over_v/sin(theta))*55E-4;
    
  // Next find correction to y from table of deflections
  double delta_y=lorentz_def->GetLorentzCorrection(XYp.X(),XYp.Y(),hit->z,theta,delta_x);
   
  double u=hit->fdc->w+delta_x;
  double v=hit->fdc->s-delta_y;
  hit->XY.Set(u*cosa+v*sina,-u*sina+v*cosa);

  return NOERROR;
}

// Compute the error matrices for the RPhi coordinates
jerror_t DTrackFitterRiemann::ComputeCRPhi(){
  // Size the matrices according to the number of hits
  unsigned int nhits=my_circle_hits.size();
  DMatrix my_CR(nhits,nhits);// Needed for correction for non-normal incidence

  // Zero the CRPhi matrix out before filling them with new data
  CRPhi.Zero();

  // Loop over the hits, fill in diagonal elements
  for (unsigned int i=0;i<nhits;i++){
    double Phi=my_circle_hits[i]->XY.Phi();
    double cosPhi=cos(Phi);
    double sinPhi=sin(Phi);
    double Phi_sinPhi_plus_cosPhi=Phi*sinPhi+cosPhi;
    double Phi_cosPhi_minus_sinPhi=Phi*cosPhi-sinPhi;
    CRPhi(i,i)=Phi_cosPhi_minus_sinPhi*Phi_cosPhi_minus_sinPhi*my_circle_hits[i]->covx
      +Phi_sinPhi_plus_cosPhi*Phi_sinPhi_plus_cosPhi*my_circle_hits[i]->covy
      +2.*Phi_sinPhi_plus_cosPhi*Phi_cosPhi_minus_sinPhi*my_circle_hits[i]->covxy;  
    my_CR(i,i)=cosPhi*cosPhi*my_circle_hits[i]->covx+sinPhi*sinPhi*my_circle_hits[i]->covy
      +2.*sinPhi*cosPhi*my_circle_hits[i]->covxy;
  }
  //printf("Errors\n");
  //my_CR.Print();
  //CRPhi.Print();

   // Correct the covariance matrices for contributions due to multiple
  // scattering
  DMatrix CRPhi_ms(nhits,nhits);
  double lambda=atan(tanl);
  double cosl=cos(lambda);
  if (cosl<EPS) cosl=EPS;
  double cosl2=cosl*cosl;
  for (unsigned int m=0;m<nhits;m++){
    double Rm=my_circle_hits[m]->XY.Mod();
    for (unsigned int k=m;k<nhits;k++){
      double Rk=my_circle_hits[k]->XY.Mod();
      unsigned int imax=(k>m)?m:k;
      for (unsigned int i=0;i<imax;i++){ 
	double zi=my_circle_hits[i]->z;
        double sigma2_ms=GetProcessNoise(my_circle_hits[i]->XY,zi);
	if (isnan(sigma2_ms)){
	  sigma2_ms=0.;
	}
        double Ri=my_circle_hits[i]->XY.Mod();
        CRPhi_ms(m,k)+=sigma2_ms*(Rk-Ri)*(Rm-Ri)/cosl2;
      }
      CRPhi_ms(k,m)=CRPhi_ms(m,k);
    }
  }
  CRPhi+=CRPhi_ms;

  // Correction for non-normal incidence of track at the intersection R=R_i
  // The correction is 
  //    CRPhi'= C*CRPhi*C+S*CR*S, where S(i,i)=R_i*kappa/2
  //                                and C(i,i)=sqrt(1-S(i,i)^2)  
  DMatrix C(nhits,nhits);
  DMatrix S(nhits,nhits);
  for (unsigned int i=0;i<my_circle_hits.size();i++){
    double stemp=my_circle_hits[i]->XY.Mod()/(4.*rc);
    double ctemp=1.-stemp*stemp;
    if (ctemp>0){
      S(i,i)=stemp;
      C(i,i)=sqrt(ctemp);
    }
    else{
      S(i,i)=0.;
      C(i,i)=1;      
    }
  }
  CRPhi=C*CRPhi*C+S*my_CR*S;
  /*
    CR.Print();
    CRPhi.Print();
  */
  return NOERROR;
}

// Compute the error matrix for the R coordinates
jerror_t DTrackFitterRiemann::ComputeCR(){
  unsigned int nhits=my_line_hits.size();

  // Zero out the matrix before filling with new elements
  CR.Zero();

  // Loop over the hits, fill in diagonal elements
  for (unsigned int i=0;i<nhits;i++){
    double Phi=my_line_hits[i]->XY.Phi();
    double cosPhi=cos(Phi);
    double sinPhi=sin(Phi);
    CR(i,i)=cosPhi*cosPhi*my_line_hits[i]->covx+sinPhi*sinPhi*my_line_hits[i]->covy
      +2.*sinPhi*cosPhi*my_line_hits[i]->covxy;
  }
  //  printf("Errors\n");
  //CR.Print();
  //CRPhi.Print();

  // Correct the covariance matrix for contributions due to multiple
  // scattering
  DMatrix CR_ms(nhits,nhits);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  if (sinl<EPS) sinl=EPS;
  double sinl4=sinl*sinl*sinl*sinl;
  for (unsigned int m=0;m<nhits;m++){
    double zm=my_line_hits[m]->z;
    for (unsigned int k=m;k<nhits;k++){
      double zk=my_line_hits[k]->z;
      unsigned int imax=(k>m)?m:k;
      for (unsigned int i=0;i<imax;i++){
	double zi=my_line_hits[i]->z;
        double sigma2_ms=GetProcessNoise(my_line_hits[i]->XY,zi);
	if (isnan(sigma2_ms)){
	  sigma2_ms=0.;
	}      
        CR_ms(m,k)+=sigma2_ms*(zk-zi)*(zm-zi)/sinl4;
      }
      CR_ms(k,m)=CR_ms(m,k);
    }
  }
  CR+=CR_ms;
  /*
  CR.Print();
  CRPhi.Print();
  */
  return NOERROR;
}
 
// Compute the covariance matrix for z
jerror_t DTrackFitterRiemann::ComputeCz(){
  unsigned int n=my_line_hits.size();
  // First zero out the matrix
  Cz.Zero();
  
  // Next deal with the diagonal components
  for (unsigned int i=0;i<n;i++){
    Cz(i,i)=my_line_hits[i]->covz;
  }
  // Correct the Cz covariance matrix for contributions due to multiple
  // scattering
  DMatrix Cz_ms(n,n);
  double lambda=atan(tanl);
  double cosl=cos(lambda);
  if (cosl<EPS) cosl=EPS;
  double cosl4=cosl*cosl*cosl*cosl;
  for (unsigned int m=0;m<n;m++){
    for (unsigned int k=m;k<n;k++){
      unsigned int imax=(k>m)?m:k;
      for (unsigned int i=0;i<imax;i++){
	double zi=my_line_hits[k]->z;
	double sigma2_ms=GetProcessNoise(my_line_hits[i]->XY,zi);
	if (isnan(sigma2_ms)){
	  sigma2_ms=0.;
	}
	Cz_ms(m,k)+=sigma2_ms*(s[k]-s[i])*(s[m]-s[i])/cosl4;
	}
      Cz_ms(k,m)=Cz_ms(m,k);
    }
  }
  Cz+=Cz_ms;

  return NOERROR;
}

// Compute the variance of the projected multiple scattering angle following 
// the formalism of Lynch and Dahl
double DTrackFitterRiemann::GetProcessNoise(const DVector2 &XY,const double z){
  // Get the material properties for this position
  double Z,rho_Z_over_A,K_rho_Z_over_A,LnI;
  DVector3 pos(XY.X(),XY.Y(),z);
  if(geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI)!=NOERROR){
	return 0.;
  }
  
  double p2=p*p;
  double F=MOLIERE_FRACTION; // Fraction of Moliere distribution to be taken into account
  double alpha=7.29735e-03; // Fine structure constant
  double one_over_beta2=1.+mass2/p2;
  double my_ds=1.;
  double chi2c=0.157*(Z+1)*rho_Z_over_A*my_ds*one_over_beta2/p2;
  double chi2a=2.007e-5*pow(Z,TWO_THIRDS)
    *(1.+3.34*Z*Z*alpha*alpha*one_over_beta2)/p2;
  double nu=0.5*chi2c/(chi2a*(1.-F));
  return (2.*chi2c*1e-6/(1.+F*F)*((1.+nu)/nu*log(1.+nu)-1.));
}


//-----------------
// FitCircle
//-----------------
jerror_t DTrackFitterRiemann::FitCircle(){
/// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
/// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
/// task of fitting points in (x,y) to a circle is transormed to the task of
/// fitting planes in (x,y, w=x^2+y^2) space
///
  unsigned nhits=my_circle_hits.size();
  DMatrix X(nhits,3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  double B0,B1,B2,Q,Q1,R,sum,diff;
  double angle,lambda_min=0.;
  // Column and row vectors of ones
  DMatrix Ones(nhits,1),OnesT(1,nhits);
  DMatrix W_sum(1,1);
  DMatrix W(nhits,nhits);
 
  // The goal is to find the eigenvector corresponding to the smallest 
  // eigenvalue of the equation
  //            lambda=n^T (X^T W X - W_sum Xavg^T Xavg)n
  // where n is the normal vector to the plane slicing the cylindrical 
  // paraboloid described by the parameterization (x,y,w=x^2+y^2),
  // and W is the weight matrix, assumed for now to be diagonal.
  // In the absence of multiple scattering, W_sum is the sum of all the 
  // diagonal elements in W.

  for (unsigned int i=0;i<nhits;i++){
    X(i,0)=my_circle_hits[i]->XY.X();
    X(i,1)=my_circle_hits[i]->XY.Y();
    X(i,2)=my_circle_hits[i]->XY.Mod2();
    Ones(i,0)=OnesT(0,i)=1.;
  }

  // Check that CRPhi is invertible 
  TDecompLU lu(CRPhi);
  if (lu.Decompose()==false){
    return UNRECOVERABLE_ERROR; // error placeholder
  }
  W=DMatrix(DMatrix::kInverted,CRPhi);
  W_sum=OnesT*(W*Ones);
  Xavg=(1./W_sum(0,0))*(OnesT*(W*X));
  
  A=DMatrix(DMatrix::kTransposed,X)*(W*X)
    -W_sum(0,0)*(DMatrix(DMatrix::kTransposed,Xavg)*Xavg);
  if(!A.IsValid()){
    return UNRECOVERABLE_ERROR;
  }

  // The characteristic equation is 
  //   lambda^3+B2*lambda^2+lambda*B1+B0=0 
  //
  B2=-(A(0,0)+A(1,1)+A(2,2));
  B1=A(0,0)*A(1,1)-A(1,0)*A(0,1)+A(0,0)*A(2,2)-A(2,0)*A(0,2)+A(1,1)*A(2,2)
    -A(2,1)*A(1,2);
  B0=-A.Determinant();
  if(B0==0 || !finite(B0)){
    return UNRECOVERABLE_ERROR;
  }

  // The roots of the cubic equation are given by 
  //        lambda1= -B2/3 + S+T
  //        lambda2= -B2/3 - (S+T)/2 + i sqrt(3)/2. (S-T)
  //        lambda3= -B2/3 - (S+T)/2 - i sqrt(3)/2. (S-T)
  // where we define some temporary variables:
  //        S= (R+sqrt(Q^3+R^2))^(1/3)
  //        T= (R-sqrt(Q^3+R^2))^(1/3)
  //        Q=(3*B1-B2^2)/9
  //        R=(9*B2*B1-27*B0-2*B2^3)/54
  //        sum=S+T;
  //        diff=i*(S-T)
  // We divide Q and R by a safety factor to prevent multiplying together 
  // enormous numbers that cause unreliable results.

  Q=(3.*B1-B2*B2)/9.e4; 
  R=(9.*B2*B1-27.*B0-2.*B2*B2*B2)/54.e6;
  Q1=Q*Q*Q+R*R;
  if (Q1<EPS){
    if (fabs(Q1)<EPS) Q1=0.;
    else Q1=sqrt(-Q1);
    // DeMoivre's theorem for fractional powers of complex numbers:  
    //      (r*(cos(angle)+i sin(angle)))^(p/q)
    //                  = r^(p/q)*(cos(p*angle/q)+i sin(p*angle/q))
    //
    double temp=100.*pow(R*R+Q1*Q1,0.16666666666666666667);
    angle=atan2(Q1,R)/3.;
    sum=2.*temp*cos(angle);
    diff=-2.*temp*sin(angle);
    // Third root
    lambda_min=-B2/3.-sum/2.+sqrt(3.)/2.*diff;
  }
  else{
    Q1=sqrt(Q1);
    // first root
    lambda_min=-B2/3+pow(R+Q1,ONE_THIRD)+pow(R-Q1,ONE_THIRD);
  }

  // Calculate the (normal) eigenvector corresponding to the eigenvalue lambda
  N.SetXYZ(1.,
	   (A(1,0)*A(0,2)-(A(0,0)-lambda_min)*A(1,2))
	   /(A(0,1)*A(2,1)-(A(1,1)-lambda_min)*A(0,2)),
	   (A(2,0)*(A(1,1)-lambda_min)-A(1,0)*A(2,1))
	   /(A(1,2)*A(2,1)-(A(2,2)-lambda_min)*(A(1,1)-lambda_min)));
  
  // Normalize: n1^2+n2^2+n3^2=1
  N.SetMag(1.0);

  // Distance to origin
  c_origin=-(N.X()*Xavg(0,0)+N.Y()*Xavg(0,1)+N.Z()*Xavg(0,2));

  // Center and radius of the circle
  double one_over_2Nz=1./(2.*N.Z());
  xc=-N.X()*one_over_2Nz;
  yc=-N.Y()*one_over_2Nz;
  rc=sqrt(1.-N.Z()*N.Z()-4.*c_origin*N.Z())*fabs(one_over_2Nz);
 
  // Phi value at "vertex"
  phi0=atan2(-xc,yc);  
  if (q<0) phi0+=M_PI;
  if(phi0<0)phi0+=2.0*M_PI;
  if(phi0>=2.0*M_PI)phi0-=2.0*M_PI;
  
  D=yc/cos(phi0)-q*rc;

  // Calculate the chisq
  //ChisqCircle();
  //chisq_source = CIRCLE;

  return NOERROR;
}

// Compute the intersections of the fitted circle with the measurements 
// through the paraboloid transform, such that the problem becomes the 
// calculation of the intersection of two planes -- i.e., a straight line.
// Store the cumulative arc length from measurement to measurement.
jerror_t DTrackFitterRiemann::ComputeIntersections(){
  double x_int0,temp,y_int0;
  double denom=N.Perp();
  for (unsigned int m=0;m<my_line_hits.size();m++){
    if (my_line_hits[m]->cdc!=NULL){
      projections[m]=my_line_hits[m]->XY;
    }
    else{
      double r2=my_line_hits[m]->XY.Mod2();
      double numer=c_origin+r2*N.z();
      double ratio=numer/denom;
      x_int0=-N.x()*ratio;
      y_int0=-N.y()*ratio;

      temp=denom*r2-numer*numer;
      // Since we will be taking a square root next, check for positive value
      if (temp<0){  
	// Not sure what to do here.  Maybe use measurement itself??
	//projections[m]=DVector2(x_int0,y_int0);
	projections[m]=my_line_hits[m]->XY;

	continue;
      }
      temp=sqrt(temp)/denom;
      
      // Choose sign of square root based on proximity to actual measurements
      double deltax=N.y()*temp;
      double deltay=-N.x()*temp;
      DVector2 XY1(x_int0+deltax,y_int0+deltay);
      DVector2 XY2(x_int0-deltax,y_int0-deltay);
      if ((XY1-my_line_hits[m]->XY).Mod2() > (XY2-my_line_hits[m]->XY).Mod2()){
	projections[m]=XY2;
      }
      else{
	projections[m]=XY1;
      }
    }
  }

  // Fill in vector of arc lengths
  double my_s=0;
  DVector2 XYold(-D*sin(phi0),D*cos(phi0));
  for (unsigned int m=0;m<my_line_hits.size();m++){
    double chord_ratio=(projections[m]-XYold).Mod()/(2.*rc);
    my_s=my_s+(chord_ratio>1.?2.*rc*M_PI_2:2.*rc*asin(chord_ratio));
    s[m]=my_s;
    XYold=projections[m];     
  }

  return NOERROR;
}


//-----------------
// FitLine
//-----------------
jerror_t DTrackFitterRiemann::FitLine(){
/// Riemann Line fit: linear regression to determine the tangent of 
/// the dip angle and the z position of the closest approach to the beam line.
  // Also computes the average B value over the hits.
  unsigned int n=projections.size();
  double sumv=0.,sumx=0.,sumy=0.,sumxx=0.,sumxy=0.;
  double Delta;
  double z=0.;

  DVector2 old_projection=projections[0];
  if (fdchits.size()==0){
    for (unsigned int k=0;k<n;k++){
      z=my_line_hits[k]->z;

      double weight=1./Cz(k,k);
      sumv+=weight;
      sumy+=z*weight;
      sumx+=s[k]*weight;
      sumxx+=s[k]*s[k]*weight;
      sumxy+=s[k]*z*weight;
    }

    Delta=(sumv*sumxx-sumx*sumx);
    // Track parameters tan(lambda) and z-vertex
    //theta=atan(Delta/(sumv*sumxy-sumy*sumx));
    //tanl=tan(M_PI_2-theta);
    tanl=(sumv*sumxy-sumx*sumy)/Delta;
    theta=M_PI_2-atan(tanl);
    z_vertex=(sumxx*sumy-sumx*sumxy)/Delta; 

    // Try to constrain the particle to come from somewhere near the center 
    // of the target if the initial fit result goes out of range in z
    if (z_vertex<Z_MIN || z_vertex>Z_MAX){
      double weight=1000.;
      sumv+=weight;
      sumy+=Z_VERTEX*weight;
      Delta= sumv*sumxx-sumx*sumx;
      tanl=(sumv*sumxy-sumx*sumy)/Delta;
      theta=M_PI_2-atan(tanl);
      z_vertex=(sumxx*sumy-sumx*sumxy)/Delta;  
    }
  }
  else{ // Got FDC hits
    for (unsigned int k=0;k<n;k++){
      z=my_line_hits[k]->z;

      // Assume errors in s dominated by errors in R 
      double weight=1./CR(k,k);
      sumv+=weight;
      sumy+=s[k]*weight;
      sumx+=z*weight;
      sumxx+=z*z*weight;
      sumxy+=s[k]*z*weight;
    }

    Delta= sumv*sumxx-sumx*sumx;
    double Delta1=sumv*sumxy-sumy*sumx;
    // Track parameters tan(lambda) and z-vertex
    tanl=Delta/Delta1;
    theta=M_PI_2-atan(tanl);

    z_vertex=-(sumxx*sumy-sumx*sumxy)/Delta1;

    // Try to constrain the particle to come from somewhere near the center 
    // of the target if the initial fit result goes out of range in z
    if (z_vertex<Z_MIN || z_vertex>Z_MAX){
      double weight=1000.;
      sumv+=weight;
      sumx+=Z_VERTEX*weight;
      sumxx+=Z_VERTEX*Z_VERTEX*weight;
      Delta= sumv*sumxx-sumx*sumx;
      Delta1=sumv*sumxy-sumy*sumx;
      // Track parameters tan(lambda) and z-vertex
      tanl=Delta/Delta1;
      theta=M_PI_2-atan(tanl);
      z_vertex=-(sumxx*sumy-sumx*sumxy)/Delta1;
    }
    
    //z_vertex=z-s[n-1]*tanl;    
  }

  return NOERROR;
}



jerror_t DTrackFitterRiemann::GetCharge(){


  return NOERROR;
}
