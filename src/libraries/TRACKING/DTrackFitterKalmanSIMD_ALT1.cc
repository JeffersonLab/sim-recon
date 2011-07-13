// $Id$
//
//    File: DTrackFitterKalmanSIMD_ALT1.cc
// Created: Tue Mar 29 09:45:14 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DTrackFitterKalmanSIMD_ALT1.h"

// Kalman engine for forward tracks.  For FDC hits only the position along 
// the wire is used in the fit
jerror_t DTrackFitterKalmanSIMD_ALT1::KalmanForward(double anneal_factor, 
					       DMatrix5x1 &S, 
					       DMatrix5x5 &C,
					       double &chisq, 
					       unsigned int &numdof){
  DMatrix1x5 H;  // Track projection matrix for cdc hits
  DMatrix5x1 H_T; // Transpose of track projection matrix for cdc hits
  DMatrix5x5 J;  // State vector Jacobian matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // Kalman gain matrix for cdc hits
  DMatrix5x1 S0,S0_; //State vector
  DMatrix5x5 I; // identity matrix
  for (unsigned int i=0;i<5;i++)I(i,i)=1.;

  // Set the "used_in_fit" flags to false for all hits
  for (unsigned int i=0;i<my_fdchits.size();i++){
    my_fdchits[i]->used_in_fit=false;
  }
  for (unsigned int i=0;i<my_cdchits.size();i++){
    my_cdchits[i]->used_in_fit=false;
  }

  // Save the starting values for C and S in the deque
  forward_traj[0].Skk=S;
  forward_traj[0].Ckk=C;

  // Initialize chi squared
  chisq=0;
  pulls.clear();

  // Initialize number of degrees of freedom
  numdof=0;

  // Variables for estimating t0 from tracking
  //mInvVarT0=mT0wires=0.;

  int num_fdc_hits=my_fdchits.size();
  int num_cdc_hits=my_cdchits.size(); 
  int cdc_index=num_cdc_hits-1;
  double old_doca=1000.;

  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size();k++){
    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(forward_traj[k].S);
    J=(forward_traj[k].J);
    Q=(forward_traj[k].Q);

    // State S is perturbation about a seed S0
    //dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*(S-S0_);

    //C=J*(C*J_T)+Q;   
    C=Q.AddSym(C.SandwichMultiply(J));

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (num_fdc_hits>0){
      if (forward_traj[k].h_id>0){
	unsigned int id=forward_traj[k].h_id-1;
	      
	double cosa=my_fdchits[id]->cosa;
	double sina=my_fdchits[id]->sina;
	double u=my_fdchits[id]->uwire;
	double v=my_fdchits[id]->vstrip;
	double x=S(state_x);
	double y=S(state_y);
	double tx=S(state_tx);
	double ty=S(state_ty);
	double du=x*cosa-y*sina-u;
	double tu=tx*cosa-ty*sina;
	double one_plus_tu2=1.+tu*tu;
	double alpha=atan(tu);
	double cosalpha=cos(alpha);
	double sinalpha=sin(alpha);
	// (signed) distance of closest approach to wire
	double doca=du*cosalpha;
	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;

	// Variance in coordinate along wire
	double V=anneal_factor*fdc_y_variance(alpha,doca,my_fdchits[id]->dE);
		
	// Difference between measurement and projection
	double Mdiff=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);
	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	H(state_x)=H_T(state_x)=sina;
	H(state_y)=H_T(state_y)=cosa;
	
	// Terms that depend on the correction for the Lorentz effect
	H(state_x)=H_T(state_x)
	  =sina+cosa*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	H(state_y)=H_T(state_y)
	=cosa-sina*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	double temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				       -2.*nr*cosalpha*sinalpha);
	H(state_tx)=H_T(state_tx)=cosa*temp;
	H(state_ty)=H_T(state_ty)=-sina*temp;
    
	// Check to see if we have multiple hits in the same plane
	if (forward_traj[k].num_hits>1){ 
	  // If we do have multiple hits, then all of the hits within some
	  // validation region are included with weights determined by how
	  // close the hits are to the track projection of the state to the
	  // "hit space".
	  vector<DMatrix5x1> Klist;
	  vector<double> Mlist;
	  vector<DMatrix1x5> Hlist;
	  vector<double> Vlist;
	  vector<double>probs;

	  // Deal with the first hit:
	  double Vtemp=V+H*C*H_T;
	  double InvV=1./Vtemp;;
       
	  //probability
	  double chi2_hit=Mdiff*Mdiff*InvV;
	  double prob_hit=exp(-0.5*chi2_hit)/sqrt(2.*M_PI*Vtemp);

	  // Cut out outliers
	  if (sqrt(chi2_hit)<NUM_SIGMA){
	    probs.push_back(prob_hit);
	    Vlist.push_back(V);
	    Hlist.push_back(H);
	    Mlist.push_back(Mdiff);
	    Klist.push_back(InvV*(C*H_T)); // Kalman gain
	  }
	  
	  // loop over the remaining hits
	  for (unsigned int m=1;m<forward_traj[k].num_hits;m++){
	    unsigned int my_id=id-m;
	    u=my_fdchits[my_id]->uwire;
	    v=my_fdchits[my_id]->vstrip;
	    double du=x*cosa-y*sina-u;
	    doca=du*cosalpha;
	    
	    // variance for coordinate along the wire
	    V=anneal_factor*fdc_y_variance(alpha,doca,my_fdchits[my_id]->dE);
	    
	    // Difference between measurement and projection
	    Mdiff=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);

	    // Update the terms in H/H_T that depend on the particular hit
	    temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				    -2.*nr*cosalpha*sinalpha);
	    H(state_tx)=H_T(state_tx)=cosa*temp;
	    H(state_ty)=H_T(state_ty)=-sina*temp;
						
	    // Calculate the kalman gain for this hit 
	    Vtemp=V+H*C*H_T;
	    InvV=1./Vtemp;
	
	    // probability
	    chi2_hit=Mdiff*Mdiff*InvV;
	    prob_hit=exp(-0.5*chi2_hit)/sqrt(2.*M_PI*Vtemp);

	    // Cut out outliers
	    if(sqrt(chi2_hit)<NUM_SIGMA){	      
	      probs.push_back(prob_hit);	
	      Mlist.push_back(Mdiff);
	      Vlist.push_back(V);
	      Hlist.push_back(H);  
	      Klist.push_back(InvV*(C*H_T));
	    }
	  }
	  double prob_tot=0.;
	  for (unsigned int m=0;m<probs.size();m++){
	    prob_tot+=probs[m];
	  }

	  // Adjust the state vector and the covariance using the hit 
	  //information
	  DMatrix5x5 sum=I;
	  DMatrix5x5 sum2;
	  for (unsigned int m=0;m<Klist.size();m++){
	    double my_prob=probs[m]/prob_tot;
	    S+=my_prob*(Mlist[m]*Klist[m]);
	    sum+=my_prob*(Klist[m]*Hlist[m]);
	    sum2+=(my_prob*my_prob*Vlist[m])*MultiplyTranspose(Klist[m]);
	  }
	  C=C.SandwichMultiply(sum)+sum2;

	  // update number of degrees of freedom
	  numdof++;

	}
	else{
	   // Variance for this hit
	  double Vtemp=V+H*C*H_T;
	  double InvV=1./Vtemp;
	
	  // Check if this hit is an outlier
	  double chi2_hit=Mdiff*Mdiff*InvV;
	  if (sqrt(chi2_hit)<NUM_SIGMA){
	    // Flag that we used this hit in the fit
	    my_fdchits[id]->used_in_fit=true;
	    
	    // Compute Kalman gain matrix
	    K=InvV*(C*H_T);
	    
	    // Update the state vector 
	    S+=Mdiff*K;
	    
	    // Update state vector covariance matrix
	    //C=C-K*(H*C);    
	    C=C.SubSym(K*(H*C));
	    
	    // Filtered residual and covariance of filtered residual
	    double R=Mdiff*(1.-H*K);   
	    double RC=V-H*(C*H_T);
	    
	    // Update chi2 for this segment
	    chisq+=R*R/RC;
	    
	    if (DEBUG_LEVEL>2){
	      printf("hit %d p %5.2f dm %5.2f sig %5.3f chi2 %5.2f z %5.2f\n",
		     id,1./S(state_q_over_p),Mdiff,sqrt(RC),R*R/RC,
		     forward_traj[k].pos.z());
	    
	    }
	      // update number of degrees of freedom
	    numdof++;
	    
	    my_fdchits[id]->yres=R/sqrt(RC);
	    pulls.push_back(pull_t(R, sqrt(fabs(RC)/anneal_factor), forward_traj[k].s));
	      
	  }
	  else{
	    // Flag that we did not use this hit after all
	    my_fdchits[id]->used_in_fit=false;
	  }
	}
	num_fdc_hits--;
      }
    }
    else if (num_cdc_hits>0){
      DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
      double z0w=origin.z();
      DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
      double uz=dir.z();
      double z=forward_traj[k].pos.z();
      DVector3 wirepos=origin+((z-z0w)/uz)*dir;

      // doca variables
      double dx=S(state_x)-wirepos.x();
      double dy=S(state_y)-wirepos.y();
      double doca=sqrt(dx*dx+dy*dy);
     
      // Check if the doca is no longer decreasing
      if (doca>old_doca){
	if(my_cdchits[cdc_index]->status==0){	
	  // Get energy loss 
	  double dedx=GetdEdx(S(state_q_over_p), 
			      forward_traj[k].K_rho_Z_over_A,
			      forward_traj[k].rho_Z_over_A,
			      forward_traj[k].LnI);
	  double tx=S(state_tx);
	  double ty=S(state_ty);	
	  double tanl=1./sqrt(tx*tx+ty*ty);
	  double sinl=sin(atan(tanl));
	  
	  // Wire direction variables
	  double ux=dir.x();
	  double uy=dir.y();
	  // Variables relating wire direction and track direction
	  double my_ux=tx-ux/uz;
	  double my_uy=ty-uy/uz;
	  double denom=my_ux*my_ux+my_uy*my_uy;
	  double dz=0.;
	  
	  // if the path length increment is small relative to the radius 
	  // of curvature, use a linear approximation to find dz	
	  bool do_brent=false;
	  double step1=mStepSizeZ;
	  double step2=mStepSizeZ;
	  if (k>=2){
	    step1=-forward_traj[k].pos.z()+forward_traj[k-1].pos.z();
	    step2=-forward_traj[k-1].pos.z()+forward_traj[k-2].pos.z();
	  }
	  //printf("step1 %f step 2 %f \n",step1,step2);
	  double two_step=step1+step2;
	  if (fabs(qBr2p*S(state_q_over_p)
		   *bfield->GetBz(S(state_x),S(state_y),z)
		   *two_step/sinl)<0.01 
	      && denom>EPS){
	    double dzw=(z-z0w)/uz;
	    dz=-((S(state_x)-origin.x()-ux*dzw)*my_ux
	       +(S(state_y)-origin.y()-uy*dzw)*my_uy)
	      /(my_ux*my_ux+my_uy*my_uy);

	    if (fabs(dz)>two_step) do_brent=true;
	  }
	  else do_brent=true;
	  if (do_brent){
	    // We have bracketed the minimum doca:  use Brent's agorithm
	    /*
	      double step_size
	      =forward_traj_cdc[k].pos.z()-forward_traj_cdc[k-1].pos.z();
	      dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	    */
	    dz=BrentsAlgorithm(z,-0.5*two_step,dedx,origin,dir,S);
	  }
	  double newz=z+dz;
	  // Check for exiting the straw
	  if (newz>endplate_z){
	    newz=endplate_z;
	    dz=endplate_z-z;
	  }
	  
	  // Step current state by dz
	  Step(z,newz,dedx,S);
	  
	  // Step reference trajectory by dz
	  Step(z,newz,dedx,S0); 
	  
	  // propagate error matrix to z-position of hit
	  StepJacobian(z,newz,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);
	  
	  // Wire position at current z
	  wirepos=origin+((newz-z0w)/uz)*dir;
	  double xw=wirepos.x();
	  double yw=wirepos.y();
	  
	  // predicted doca taking into account the orientation of the wire
	  dy=S(state_y)-yw;
	  dx=S(state_x)-xw;      
	  double cosstereo=cos(my_cdchits[cdc_index]->hit->wire->stereo);
	  double d=sqrt(dx*dx+dy*dy)*cosstereo;
	  
	  // Track projection
	  double cosstereo2_over_d=cosstereo*cosstereo/d;
	  H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
	  H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;
      
	  //H.Print();
	  
	  // The next measurement
	  double dm=0.;
	  double Vc=0.2133; //1.6*1.6/12.;
	  
	  if (fit_type==kTimeBased){
	    /*
	    dm=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift-mT0
				-forward_traj[k-1].t);
	    */
	    double tdrift=my_cdchits[cdc_index]->hit->tdrift-mT0
		-forward_traj[k-1].t;
	    if (tdrift>0.) dm=0.02887*sqrt(tdrift)-1.315e-5*tdrift;
	 


	    // variance
	    //Vc=CDC_VARIANCE*10.0;
	    if (newz<endplate_z)
	      Vc=cdc_variance(dm)+CDC_DRIFT_SPEED*CDC_DRIFT_SPEED*mVarT0;
	    }
	  else if (USE_T0_FROM_WIRES && mInvVarT0>EPS){
	    dm=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift
				-mT0wires
				-forward_traj[k-1].t);
	    if (newz<endplate_z)
	      Vc=cdc_variance(d)+CDC_DRIFT_SPEED*CDC_DRIFT_SPEED/mInvVarT0;
	  }

	  // inverse variance including prediction
	  double InvV1=1./(Vc+H*(C*H_T));
	  if (InvV1<0.){
	    if (DEBUG_LEVEL>0)
	      _DBG_ << "Negative variance???" << endl;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  if (DEBUG_LEVEL>2)
	    printf("Ring %d straw %d pred %f meas %f V %f %f sig %f\n",
		   my_cdchits[cdc_index]->hit->wire->ring,
		   my_cdchits[cdc_index]->hit->wire->straw,
		   d,dm,Vc,1./InvV1,1./sqrt(InvV1));
	  // Check if this hit is an outlier
	  double chi2_hit=(dm-d)*(dm-d)*InvV1;
	  if (sqrt(chi2_hit)<NUM_SIGMA){
	    // Flag the place along the reference trajectory closest to the 
	    // doca position
	    forward_traj[k-1].h_id=1000+cdc_index;

	    // Flag that we used this hit
	    my_cdchits[cdc_index]->used_in_fit=true;
	    
	    // Compute KalmanSIMD gain matrix
	    K=InvV1*(C*H_T);
	    
	    // Update the state vector
	    double res=dm-d;
	    S+=res*K;
	      
	    // Update state vector covariance matrix
	    //C=C-K*(H*C);    
	    C=C.SubSym(K*(H*C));

	    // Residual
	    res*=1.-H*K;
	  
	  
	    // Update chi2 for this segment
	    double err2 = Vc-H*(C*H_T);
	    chisq+=anneal_factor*res*res/err2;
	 	      
	    // update number of degr  virtual jerror_t SmoothForward(DMatrix5x1 &S);   ees of freedom
	    numdof++;
	    
	    my_cdchits[cdc_index]->residual=res/sqrt(err2);
	    pulls.push_back(pull_t(res, sqrt(fabs(err2/anneal_factor)), forward_traj[k].s));
	    
	  }
	  else{
	    // Flag that we did not use this hit
	    my_cdchits[cdc_index]->used_in_fit=false;
	  }

	  // Step C back to the z-position on the reference trajectory
	  StepJacobian(newz,z,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);
	  
	  // Step S to current position on the reference trajectory
	  Step(newz,z,dedx,S);
	}

	// new wire origin and direction
	if (cdc_index>0){
	  cdc_index--;
	  origin=my_cdchits[cdc_index]->hit->wire->origin;
	  dir=my_cdchits[cdc_index]->hit->wire->udir;
	}
      
	// Update the wire position
	uz=dir.z();
	z0w=origin.z();
	wirepos=origin+((z-z0w)/uz)*dir;
	
	// new doca
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca=sqrt(dx*dx+dy*dy);
	num_cdc_hits--;
	if (cdc_index==0) num_cdc_hits=0;
      }
      old_doca=doca;
    }

    // Save the current state and covariance matrix in the deque
    forward_traj[k].Skk=S;
    forward_traj[k].Ckk=C;

  }
  
  // If chisq is still zero after the fit, something went wrong...
  if (chisq<EPS) return UNRECOVERABLE_ERROR;

  chisq*=anneal_factor;
  

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  if (DEBUG_LEVEL>0){
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
    cout << "Momentum " << 1./S(state_q_over_p) <<endl;
  }
    
  return NOERROR;
}
