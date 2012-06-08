// $Id$
//
//    File: DTrackFitterKalmanSIMD_ALT1.cc
// Created: Tue Mar 29 09:45:14 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DTrackFitterKalmanSIMD_ALT1.h"

// Kalman engine for forward tracks.  For FDC hits only the position along 
// the wire is used in the fit
kalman_error_t DTrackFitterKalmanSIMD_ALT1::KalmanForward(double anneal_factor, 
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
  DMatrix5x5 Ctest; // Covariance matrix

  // Vectors for cdc wires
  DVector3 origin,dir,wirepos;

  // Set used_in_fit flags to false for fdc and cdc hits
  for (unsigned int i=0;i<cdc_updates.size();i++) cdc_updates[i].used_in_fit=false;
  for (unsigned int i=0;i<fdc_updates.size();i++) fdc_updates[i].used_in_fit=false;
  
  // Save the starting values for C and S in the deque
  forward_traj[0].Skk=S;
  forward_traj[0].Ckk=C;

  // Initialize chi squared
  chisq=0;

  // Initialize number of degrees of freedom
  numdof=0;

  double fdc_chi2cut
    =anneal_factor*anneal_factor*NUM_FDC_SIGMA_CUT*NUM_FDC_SIGMA_CUT; 
  double cdc_chi2cut
    =anneal_factor*anneal_factor*NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;

  // Variables for estimating t0 from tracking
  //mInvVarT0=mT0wires=0.;

  unsigned int num_fdc_hits=my_fdchits.size();
  unsigned int num_cdc_hits=my_cdchits.size(); 
  unsigned int cdc_index=0;
  if (num_cdc_hits>0) cdc_index=num_cdc_hits-1;
  bool more_cdc_measurements=(num_cdc_hits>0);
  double old_doca=1000.;

  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size();k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (C(0,0)<0 || C(1,1)<0 || C(2,2)<0 || C(3,3)<0 || C(4,4)<0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
      return BROKEN_COVARIANCE_MATRIX;
    }

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(forward_traj[k].S);
    J=(forward_traj[k].J);
    Q=(forward_traj[k].Q);

    // State S is perturbation about a seed S0
    //dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*(S-S0_);

    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
	}
      return MOMENTUM_OUT_OF_RANGE;
    }


    //C=J*(C*J_T)+Q;   
    C=Q.AddSym(C.SandwichMultiply(J));

    // Save the current state and covariance matrix in the deque
    forward_traj[k].Skk=S;
    forward_traj[k].Ckk=C;
    
    // Save the current state of the reference trajectory
    S0_=S0;

    // Z position along the trajectory 
    double z=forward_traj[k].pos.z();

    // Add the hit
    if (num_fdc_hits>0){
      if (forward_traj[k].h_id>0 && forward_traj[k].h_id<1000){
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

	// t0 estimate
	if (fit_type==kWireBased && id==mMinDriftID-1){
	  mT0MinimumDriftTime=mMinDriftTime
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	} 
	
	// Variance in coordinate along wire
	double tanl=1./sqrt(tx*tx+ty*ty);
	double V=fdc_y_variance(tanl,doca,my_fdchits[id]->dE);
		
	// Difference between measurement and projection
	double Mdiff=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);
	
       	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
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
	  double prob_hit=exp(-0.5*chi2_hit)/sqrt(M_TWO_PI*Vtemp);

	  // Cut out outliers
	  if (chi2_hit<fdc_chi2cut){
	    probs.push_back(prob_hit);
	    Vlist.push_back(V);
	    Hlist.push_back(H);
	    Mlist.push_back(Mdiff);
	    Klist.push_back(InvV*(C*H_T)); // Kalman gain

	    fdc_updates[id].used_in_fit=true;
	  }
	  
	  // loop over the remaining hits
	  for (unsigned int m=1;m<forward_traj[k].num_hits;m++){
	    unsigned int my_id=id-m;
	    fdc_updates[my_id].used_in_fit=true;

	    u=my_fdchits[my_id]->uwire;
	    v=my_fdchits[my_id]->vstrip;
	    double du=x*cosa-y*sina-u;
	    doca=du*cosalpha;

	    // t0 estimate
	    if (fit_type==kWireBased && my_id==mMinDriftID-1){
	      mT0MinimumDriftTime=mMinDriftTime
		-forward_traj[k].t*TIME_UNIT_CONVERSION;
	    } 
	    
	    // variance for coordinate along the wire
	    V=fdc_y_variance(alpha,doca,my_fdchits[my_id]->dE);
	    
	    // Difference between measurement and projection
	    Mdiff=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);

	    // Update the terms in H/H_T that depend on the particular hit
	    temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				    -2.*nr*cosalpha*sinalpha);
	    H(state_tx)=H_T(state_tx)=cosa*temp;
	    H(state_ty)=H_T(state_ty)=-sina*temp;
						
	    // Calculate the kalman gain for this hit 
	    //Vtemp=V+H*C*H_T;
	    Vtemp=V+C.SandwichMultiply(H_T);
	    InvV=1./Vtemp;

	    // probability
	    chi2_hit=Mdiff*Mdiff*InvV;
	    prob_hit=exp(-0.5*chi2_hit)/sqrt(M_TWO_PI*Vtemp);

	    // Cut out outliers
	    if(chi2_hit<fdc_chi2cut){	      
	      probs.push_back(prob_hit);	
	      Mlist.push_back(Mdiff);
	      Vlist.push_back(V);
	      Hlist.push_back(H);   
	      Klist.push_back(InvV*(C*H_T));
	    }
	  }
	  double prob_tot=1e-100;
	  for (unsigned int m=0;m<probs.size();m++){
	    prob_tot+=probs[m];
	  }

	  // Adjust the state vector and the covariance using the hit 
	  //information
	  DMatrix5x5 sum=I5x5;
	  DMatrix5x5 sum2;
	  for (unsigned int m=0;m<Klist.size();m++){
	    double my_prob=probs[m]/prob_tot;
	    S+=my_prob*(Mlist[m]*Klist[m]);
	    sum+=my_prob*(Klist[m]*Hlist[m]);
	    sum2+=(my_prob*my_prob*Vlist[m])*MultiplyTranspose(Klist[m]);
	  }
	  C=C.SandwichMultiply(sum)+sum2;

	  if (fit_type==kTimeBased){
	    for (unsigned int m=0;m<forward_traj[k].num_hits;m++){
	      unsigned int my_id=id-m;
	      if (fdc_updates[my_id].used_in_fit){
		fdc_updates[my_id].S=S;
		fdc_updates[my_id].C=C; 
		fdc_updates[my_id].tflight
		  =forward_traj[k].t*TIME_UNIT_CONVERSION;  
		fdc_updates[my_id].pos=forward_traj[k].pos;
		fdc_updates[my_id].dEdx=0.;
		fdc_updates[my_id].B=forward_traj[k].B;
		fdc_updates[my_id].s=forward_traj[k].s;
	      }
	    }
	  }
	  
	  // update chi2
	  for (unsigned int m=0;m<Klist.size();m++){
	    chisq+=(probs[m]/prob_tot)*(1.-Hlist[m]*Klist[m])*Mlist[m]*Mlist[m]/Vlist[m];
	  }

	  // update number of degrees of freedom
	  numdof++;
	}
	else{
	   // Variance for this hit
	  //double Vtemp=V+H*C*H_T;
	  double Vtemp=V+C.SandwichMultiply(H_T);
	  double InvV=1./Vtemp;
	
	  // Check if this hit is an outlier
	  double chi2_hit=Mdiff*Mdiff*InvV;
	  if (chi2_hit<fdc_chi2cut){
	    if (DEBUG_HISTS && fit_type==kTimeBased && id==0){
	      fdc_yres->Fill(doca,Mdiff*my_fdchits[id]->dE/
			     (2.6795e-4*FDC_CATHODE_SIGMA));
	    }     

	    // Compute Kalman gain matrix
	    K=InvV*(C*H_T);
	    
	    // Update the state vector 
	    S+=Mdiff*K;
	    
	    // Update state vector covariance matrix
	    //C=C-K*(H*C);    
	    C=C.SubSym(K*(H*C));

	    if (fit_type==kTimeBased){
	      fdc_updates[id].S=S;
	      fdc_updates[id].C=C;
	      fdc_updates[id].tflight
		=forward_traj[k].t*TIME_UNIT_CONVERSION;  
	      fdc_updates[id].pos=forward_traj[k].pos;
	      fdc_updates[id].dEdx=0.;
	      fdc_updates[id].B=forward_traj[k].B;
	      fdc_updates[id].s=forward_traj[k].s;
	    }
	    fdc_updates[id].used_in_fit=true;
	    
	    // Update chi2 for this segment
	    chisq+=(1.-H*K)*Mdiff*Mdiff/V;
		    
	    if (DEBUG_LEVEL>2){
	      printf("hit %d p %5.2f dm %5.2f chi2 %5.2f z %5.2f\n",
		     id,1./S(state_q_over_p),Mdiff,(1.-H*K)*Mdiff*Mdiff/V,
		     forward_traj[k].pos.z());
	    
	    }
	      // update number of degrees of freedom
	    numdof++;


	    break_point_fdc_index=id;
	    break_point_step_index=k;
	  }
	}
	if (num_fdc_hits>=forward_traj[k].num_hits)
	  num_fdc_hits-=forward_traj[k].num_hits;
      }
    }
    else if (more_cdc_measurements && z<endplate_z){
      origin=my_cdchits[cdc_index]->hit->wire->origin;
      double z0w=origin.z();
      dir=my_cdchits[cdc_index]->hit->wire->udir;
      double uz=dir.z();   
      wirepos=origin+((z-z0w)/uz)*dir;
    
      // doca variables
      double dx=S(state_x)-wirepos.x();
      double dy=S(state_y)-wirepos.y();
      double doca=sqrt(dx*dx+dy*dy);
     
      // Check if the doca is no longer decreasing
      if (doca>old_doca+EPS3){
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
	    step1=-forward_traj[k].pos.z()+forward_traj[k_minus_1].pos.z();
	    step2=-forward_traj[k_minus_1].pos.z()+forward_traj[k-2].pos.z();
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

	    if (fabs(dz)>two_step || dz<0) do_brent=true;
	  }
	  else do_brent=true;
	  if (do_brent){
	    // We have bracketed the minimum doca:  use Brent's agorithm
	    /*
	      double step_size
	      =forward_traj_cdc[k].pos.z()-forward_traj_cdc[k_minus_1].pos.z();
	      dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	    */
	    //dz=BrentsAlgorithm(z,-0.5*two_step,dedx,origin,dir,S);
	    dz=BrentsAlgorithm(z,-mStepSizeZ,dedx,origin,dir,S);
	  }
	  double newz=z+dz;
	  // Check for exiting the straw
	  if (newz>endplate_z){
	    newz=endplate_z;
	    dz=endplate_z-z;
	  }
	  // Step the state and covariance through the field
	  int num_steps=0;
	  double dz3=0.;
	  double my_dz=0.;
	  if (fabs(dz)>mStepSizeZ){
	    my_dz=(dz>0?1.0:-1.)*mStepSizeZ;
	    num_steps=int(fabs(dz/my_dz));
	    dz3=dz-num_steps*my_dz;
	    
	    double my_z=z;
	    for (int m=0;m<num_steps;m++){
	      newz=my_z+my_dz;
	      
	      // Step current state by my_dz
	      Step(z,newz,dedx,S);
	      	      
	      // propagate error matrix to z-position of hit
	      StepJacobian(z,newz,S0,dedx,J);
	      //C=J*C*J.Transpose();
	      C=C.SandwichMultiply(J);
	      
	      // Step reference trajectory by my_dz
	      Step(z,newz,dedx,S0); 
	      
	      my_z=newz;
	    }
	    
	    newz=my_z+dz3;
	    
	    // Step current state by dz3
	    Step(my_z,newz,dedx,S);	  
	    
	    // propagate error matrix to z-position of hit
	    StepJacobian(my_z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz3
	    Step(my_z,newz,dedx,S0); 	    
	  }
	  else{
	    // Step current state by dz
	    Step(z,newz,dedx,S);
	    
	    // propagate error matrix to z-position of hit
	    StepJacobian(z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz
	    Step(z,newz,dedx,S0); 
	  }
	  
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
	    double tdrift=my_cdchits[cdc_index]->hit->tdrift-mT0
		-forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	    dm=cdc_drift_distance(tdrift,forward_traj[k].B);
	    
	    // variance
	    double tx=S(state_tx),ty=S(state_ty);
	    double tanl=1./sqrt(tx*tx+ty*ty);
	    Vc=cdc_forward_variance(tanl,tdrift);	  
	    double temp=1./(1131.+2.*140.7*dm);
	    Vc+=mVarT0*temp*temp;
	  }
	  // t0 estimate
	  else if (cdc_index==mMinDriftID-1000){
	    mT0MinimumDriftTime=mMinDriftTime
	      -forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	  } 

	  // Residual
	  double res=dm-d;

	  // inverse variance including prediction
	  double InvV1=1./(Vc+H*(C*H_T));
	  if (InvV1<0.){
	    if (DEBUG_LEVEL>0)
	      _DBG_ << "Negative variance???" << endl;
	    return NEGATIVE_VARIANCE;
	  }
	  	 
	  // Check if this hit is an outlier
	  double chi2_hit=res*res*InvV1;
	  if (chi2_hit<cdc_chi2cut){
	    // Compute KalmanSIMD gain matrix
	    K=InvV1*(C*H_T);
	    
        
	    // Update state vector covariance matrix
	    //C=C-K*(H*C);    
	    Ctest=C.SubSym(K*(H*C));
	    //Ctest=C.SandwichMultiply(I5x5-K*H)+Vc*MultiplyTranspose(K);	 
	    // Check that Ctest is positive definite
	    if (Ctest(0,0)>0 && Ctest(1,1)>0 && Ctest(2,2)>0 && Ctest(3,3)>0 
		&& Ctest(4,4)>0){
	      C=Ctest;

	      // Flag the place along the reference trajectory with hit id
	      forward_traj[k_minus_1].h_id=1000+cdc_index;
      	      
	      // Update the state vector
	      double res=dm-d;
	      S+=res*K;
	      
	      // Store the "improved" values of the state and covariance matrix
	      if (fit_type==kTimeBased){
		cdc_updates[cdc_index].S=S;
		cdc_updates[cdc_index].C=C;	 
		cdc_updates[cdc_index].tflight
		  =forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
		cdc_updates[cdc_index].pos=forward_traj[k_minus_1].pos;
		cdc_updates[cdc_index].dEdx=dedx;
		cdc_updates[cdc_index].B=forward_traj[k_minus_1].B;
		cdc_updates[cdc_index].s=forward_traj[k_minus_1].s;
	     	      
		// propagate error matrix to z-position of hit
		StepJacobian(newz,forward_traj[k_minus_1].pos.z(),
			     cdc_updates[cdc_index].S,dedx,J);
		cdc_updates[cdc_index].C
		  =cdc_updates[cdc_index].C.SandwichMultiply(J);	 
		
		// Step state back to previous z position
		Step(newz,forward_traj[k_minus_1].pos.z(),dedx,cdc_updates[cdc_index].S);
	      }
	      cdc_updates[cdc_index].used_in_fit=true;
 
	      // Update chi2 for this segment
	      chisq+=(1.-H*K)*res*res/Vc;

	      if (DEBUG_LEVEL>2)
		printf("Ring %d straw %d pred %f meas %f chi2 %f\n",
		       my_cdchits[cdc_index]->hit->wire->ring,
		       my_cdchits[cdc_index]->hit->wire->straw,
		       d,dm,(1.-H*K)*res*res/Vc);
	      
	      // update number of degrees of freedom
	      numdof++;

	      break_point_cdc_index=cdc_index;
	      break_point_step_index=k_minus_1;
	    }
	  }

	  if (num_steps==0){
	    // Step C back to the z-position on the reference trajectory
	    StepJacobian(newz,z,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step S to current position on the reference trajectory
	    Step(newz,z,dedx,S);
	  }
	  else{
	    double my_z=newz;
	    for (int m=0;m<num_steps;m++){
	      z=my_z-my_dz;
	      
	      // Step C along z
	      StepJacobian(my_z,z,S0,dedx,J);
	      //C=J*C*J.Transpose();
	      C=C.SandwichMultiply(J);
	    
	      // Step S along z
	      Step(my_z,z,dedx,S); 

	      // Step S0 along z
	      Step(my_z,z,dedx,S0);
	    
	      my_z=z;
	    }
	    z=my_z-dz3;
	    
	    // Step C back to the z-position on the reference trajectory
	    StepJacobian(my_z,z,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step S to current position on the reference trajectory
	    Step(my_z,z,dedx,S);
	  }
	  
	}

	// new wire origin and direction
	if (cdc_index>0){
	  cdc_index--;
	  origin=my_cdchits[cdc_index]->hit->wire->origin;
	  dir=my_cdchits[cdc_index]->hit->wire->udir;
	}
	else more_cdc_measurements=false;
      
	// Update the wire position
	uz=dir.z();
	z0w=origin.z();
	wirepos=origin+((z-z0w)/uz)*dir;
	
	// new doca
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca=sqrt(dx*dx+dy*dy);
      }
      old_doca=doca;
    }
  }
  
  // Check that there were enough hits to make this a valid fit
  if (numdof<6){
    chisq=MAX_CHI2;
    numdof=0;
   
    return INVALID_FIT;
  }

  //  chisq*=anneal_factor;
  numdof-=5;

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  if (DEBUG_LEVEL>0){
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
    cout << "Momentum " << 1./S(state_q_over_p) <<endl;
  }

  // Check if we have a kink in the track or threw away too many cdc hits
  if (cdc_updates.size()>=6){
    unsigned int num_good=0; 
    for (unsigned int j=0;j<cdc_updates.size();j++){
      if (cdc_updates[j].used_in_fit) num_good++;
    }
    if (break_point_cdc_index>4) return BREAK_POINT_FOUND;
    if (double(num_good)/double(my_cdchits.size())<MINIMUM_HIT_FRACTION) return PRUNED_TOO_MANY_HITS;
  }
    
  return FIT_SUCCEEDED;
}
