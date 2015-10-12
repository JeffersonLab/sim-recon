// $Id$
//
//    File: DTrackFitterKalmanSIMD_ALT1.cc
// Created: Tue Mar 29 09:45:14 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DTrackFitterKalmanSIMD_ALT1.h"

// Kalman engine for forward tracks.  For FDC hits only the position along 
// the wire is used in the fit
kalman_error_t DTrackFitterKalmanSIMD_ALT1::KalmanForward(double fdc_anneal_factor, double cdc_anneal_factor,
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
  //  double Vc=0.2028; // covariance for cdc wires =1.56*1.56/12.;
  double Vc=0.0507*1.15;

  // Vectors for cdc wires
  DVector2 origin,dir,wirepos;
  double z0w=0.; // origin in z for wire

  // Set used_in_fit flags to false for fdc and cdc hits
  unsigned int num_cdc=cdc_updates.size();
  unsigned int num_fdc=fdc_updates.size();
  for (unsigned int i=0;i<num_cdc;i++) cdc_updates[i].used_in_fit=false;
  for (unsigned int i=0;i<num_fdc;i++) fdc_updates[i].used_in_fit=false;
  for (unsigned int i=0;i<forward_traj.size();i++){
    if (forward_traj[i].h_id>999)
      forward_traj[i].h_id=0;
  }
  
  // Save the starting values for C and S in the deque
  forward_traj[break_point_step_index].Skk=S;
  forward_traj[break_point_step_index].Ckk=C;

  // Initialize chi squared
  chisq=0;

  // Initialize number of degrees of freedom
  numdof=0;
  
  double my_cdc_anneal=cdc_anneal_factor*cdc_anneal_factor;
  double my_fdc_anneal=fdc_anneal_factor*fdc_anneal_factor;
 
  double var_fdc_cut=NUM_FDC_SIGMA_CUT*NUM_FDC_SIGMA_CUT;
  double fdc_chi2cut=my_fdc_anneal*var_fdc_cut;

  double var_cdc_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
  double cdc_chi2cut=my_cdc_anneal*var_cdc_cut;

  // Variables for estimating t0 from tracking
  //mInvVarT0=mT0wires=0.;

  unsigned int num_fdc_hits=break_point_fdc_index+1;
  unsigned int max_num_fdc_used_in_fit=num_fdc_hits;
  unsigned int num_cdc_hits=my_cdchits.size(); 
  unsigned int cdc_index=0;
  if (num_cdc_hits>0) cdc_index=num_cdc_hits-1;
  bool more_cdc_measurements=(num_cdc_hits>0);
  double old_doca2=1e6;

  if (num_fdc_hits+num_cdc_hits<MIN_HITS_FOR_REFIT){
    cdc_chi2cut=1000.0;
    fdc_chi2cut=1000.0;
  }

  if (more_cdc_measurements){
    origin=my_cdchits[cdc_index]->origin;  
    dir=my_cdchits[cdc_index]->dir;   
    z0w=my_cdchits[cdc_index]->z0wire;
    wirepos=origin+(forward_traj[break_point_step_index].z-z0w)*dir;
  }

  S0_=(forward_traj[break_point_step_index].S);
  for (unsigned int k=break_point_step_index+1;k<forward_traj.size();k++){
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
      break_point_fdc_index=(3*num_fdc)/4;
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
    double z=forward_traj[k].z;

    // Add the hit
    if (num_fdc_hits>0){
      if (forward_traj[k].h_id>0 && forward_traj[k].h_id<1000){
	unsigned int id=forward_traj[k].h_id-1;
	
	double cosa=my_fdchits[id]->cosa;
	double sina=my_fdchits[id]->sina;
	double u=my_fdchits[id]->uwire;
	double v=my_fdchits[id]->vstrip;

	// Position and direction from state vector
	double x=S(state_x);
	double y=S(state_y);
	double tx=S(state_tx);
	double ty=S(state_ty);

	// Projected position along the wire without doca-dependent corrections
	double vpred_uncorrected=x*sina+y*cosa;

	// Projected postion in the plane of the wires transverse to the wires
	double upred=x*cosa-y*sina;

	// Direction tangent in the u-z plane
	double tu=tx*cosa-ty*sina;
	double alpha=atan(tu);
	double cosalpha=cos(alpha);
	double cosalpha2=cosalpha*cosalpha;
	double sinalpha=sin(alpha);

	// (signed) distance of closest approach to wire
	double doca=(upred-u)*cosalpha;

	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;
	
	// Variance in coordinate along wire
	double V=my_fdchits[id]->vvar;
		
	// Difference between measurement and projection
	double tv=tx*sina+ty*cosa;
	double Mdiff=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
						-tv*sinalpha
						));

	if (DEBUG_HISTS && fit_type==kTimeBased){
	  fdc_dy_vs_d->Fill(doca,v-vpred_uncorrected);
	}
	
       	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	double temp2=nz_sinalpha_plus_nr_cosalpha
	  -tv*sinalpha
	  ;
	H_T(state_x)=sina+cosa*cosalpha*temp2;	
	H_T(state_y)=cosa-sina*cosalpha*temp2;	
       
	double cos2_minus_sin2=cosalpha2-sinalpha*sinalpha;
	double fac=nz*cos2_minus_sin2-2.*nr*cosalpha*sinalpha;
	double doca_cosalpha=doca*cosalpha;
	double temp=doca_cosalpha*fac;	
	H_T(state_tx)=cosa*temp
	  -doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
	  ;
	H_T(state_ty)=-sina*temp
	  -doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
	  ;

	// Matrix transpose H_T -> H
	H(state_x)=H_T(state_x);
	H(state_y)=H_T(state_y);
	H(state_tx)=H_T(state_tx);
	H(state_ty)=H_T(state_ty);
    
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
	  vector<unsigned int>used_ids;

	  // Deal with the first hit:
	  //double Vtemp=V+H*C*H_T;
	  double Vtemp=V+C.SandwichMultiply(H_T);
	  double InvV=1./Vtemp;;
       
	  //probability
	  double chi2_hit=Mdiff*Mdiff*InvV;
	  double prob_hit=exp(-0.5*chi2_hit)/sqrt(M_TWO_PI*Vtemp);

	  // Cut out outliers
	  if (chi2_hit<fdc_chi2cut && my_fdchits[id]->status==good_hit){
	    probs.push_back(prob_hit);
	    Vlist.push_back(V);
	    Hlist.push_back(H);
	    Mlist.push_back(Mdiff);
	    Klist.push_back(InvV*(C*H_T)); // Kalman gain

	    used_ids.push_back(id);
	    fdc_updates[id].used_in_fit=true;
	  }
	  
	  // loop over the remaining hits
	  for (unsigned int m=1;m<forward_traj[k].num_hits;m++){
	    unsigned int my_id=id-m;
	    if (my_fdchits[my_id]->status==good_hit){
	      u=my_fdchits[my_id]->uwire;
	      v=my_fdchits[my_id]->vstrip;

	      // Doca to this wire
	      doca=(upred-u)*cosalpha;
	    
	      // variance for coordinate along the wire
	      V=my_fdchits[my_id]->vvar;
	      
	      // Difference between measurement and projection
	      Mdiff=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
					       -tv*sinalpha
					       ));
	      
	      // Update the terms in H/H_T that depend on the particular hit    
	      doca_cosalpha=doca*cosalpha;
	      temp=doca_cosalpha*fac;	
	      H_T(state_tx)=cosa*temp	 
		-doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
		;
	      H_T(state_ty)=-sina*temp
		-doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
		;
	      
	      // Matrix transpose H_T -> H
	      H(state_tx)=H_T(state_tx);
	      H(state_ty)=H_T(state_ty);
	      
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
		
		used_ids.push_back(my_id);
		fdc_updates[my_id].used_in_fit=true;
		
	      }
	    }
	  }
	  double prob_tot=1e-100;
	  for (unsigned int m=0;m<probs.size();m++){
	    prob_tot+=probs[m];
	  }
	  
	  // Adjust the state vector and the covariance using the hit 
	  //information
	  bool skip_plane=(my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP);
	  if (skip_plane==false){
	    DMatrix5x5 sum=I5x5;
	    DMatrix5x5 sum2;
	    for (unsigned int m=0;m<Klist.size();m++){
	      double my_prob=probs[m]/prob_tot;
	      S+=my_prob*(Mlist[m]*Klist[m]);
	      sum+=my_prob*(Klist[m]*Hlist[m]);
	      sum2+=(my_prob*my_prob*Vlist[m])*MultiplyTranspose(Klist[m]);
	    }
	    C=C.SandwichMultiply(sum)+sum2;
	  }
	  
	  for (unsigned int m=0;m<Hlist.size();m++){
	    unsigned int my_id=used_ids[m];
	    double scale=(skip_plane)?1.:(1.-Hlist[m]*Klist[m]);    
	    fdc_updates[my_id].S=S;
	    fdc_updates[my_id].C=C; 
	    fdc_updates[my_id].tdrift
	      =my_fdchits[my_id]->t-forward_traj[k].t*TIME_UNIT_CONVERSION-mT0;
	    fdc_updates[my_id].tcorr=fdc_updates[my_id].tdrift; // temporary!
	    fdc_updates[my_id].residual=scale*Mlist[m];
	    fdc_updates[my_id].variance=scale*Vlist[m];
	    fdc_updates[my_id].doca=doca;
	    
	    // update chi2
	    if (skip_plane==false){
	      chisq+=(probs[m]/prob_tot)*(1.-Hlist[m]*Klist[m])*Mlist[m]*Mlist[m]/Vlist[m];
	    }
	  }

	  // update number of degrees of freedom
	  if (skip_plane==false){
	    numdof++;
	  }
	}
	else{
	   // Variance for this hit
	  //double Vtemp=V+H*C*H_T;
	  double Vproj=C.SandwichMultiply(H_T);
	  double Vtemp=V+Vproj;
	  double InvV=1./Vtemp;
	
	  // Check if this hit is an outlier
	  double chi2_hit=Mdiff*Mdiff*InvV;
	  if (chi2_hit<fdc_chi2cut){
	    // Compute Kalman gain matrix
	    K=InvV*(C*H_T);
	    
	    bool skip_plane=(my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP);
	    if (skip_plane==false){
	      // Update the state vector 
	      S+=Mdiff*K;
	    
	      // Update state vector covariance matrix
	      //C=C-K*(H*C);    
	      C=C.SubSym(K*(H*C));
	    }

	    // Store the "improved" values for the state vector and covariance
	    double scale=(skip_plane)?1.:(1.-H*K);
	    fdc_updates[id].S=S;
	    fdc_updates[id].C=C;
	    fdc_updates[id].tdrift
	      =my_fdchits[id]->t-forward_traj[k].t*TIME_UNIT_CONVERSION-mT0;
	    fdc_updates[id].tcorr=fdc_updates[id].tdrift; // temporary!
	    fdc_updates[id].residual=scale*Mdiff;
	    fdc_updates[id].variance=scale*V;
	    fdc_updates[id].doca=doca;
	    fdc_updates[id].used_in_fit=true;
	    
	    if (skip_plane==false){
	      // Update chi2 for this segment
	      chisq+=scale*Mdiff*Mdiff/V;
	      // update number of degrees of freedom
	      numdof++;
	    }
		    
	    if (DEBUG_LEVEL>10){
	      printf("hit %d p %5.2f t %f dm %5.2f sig %f chi2 %5.2f z %5.2f\n",
		     id,1./S(state_q_over_p),fdc_updates[id].tdrift,Mdiff,sqrt(V),(1.-H*K)*Mdiff*Mdiff/V,
		     forward_traj[k].z);
	    
	    }
	      


	    break_point_fdc_index=id;
	    break_point_step_index=k;
	  }
	}
	if (num_fdc_hits>=forward_traj[k].num_hits)
	  num_fdc_hits-=forward_traj[k].num_hits;
      }
    }
    else if (more_cdc_measurements /* && z<endplate_z*/){   
      // new wire position
      wirepos=origin;
      wirepos+=(z-z0w)*dir;

      // doca variables
      double dx=S(state_x)-wirepos.X();
      double dy=S(state_y)-wirepos.Y();
      double doca2=dx*dx+dy*dy;
     
      // Check if the doca is no longer decreasing
      if (doca2>old_doca2 /* && z<endplate_z */){
	if(my_cdchits[cdc_index]->status==good_hit){
	  double newz=z;
	
	  // Get energy loss 
	  double dedx=0.;
	  if (CORRECT_FOR_ELOSS){
	    dedx=GetdEdx(S(state_q_over_p), 
			 forward_traj[k].K_rho_Z_over_A,
			 forward_traj[k].rho_Z_over_A,
			 forward_traj[k].LnI);
	  }
	  
	  // track direction variables
	  double tx=S(state_tx);
	  double ty=S(state_ty);	
	  double tanl=1./sqrt(tx*tx+ty*ty);
	  double sinl=sin(atan(tanl));
	  
	  // Wire direction variables
	  double ux=dir.X();
	  double uy=dir.Y();
	  // Variables relating wire direction and track direction
	  double my_ux=tx-ux;
	  double my_uy=ty-uy;
	  double denom=my_ux*my_ux+my_uy*my_uy;
	  double dz=0.;
	  
	  // if the path length increment is small relative to the radius 
	  // of curvature, use a linear approximation to find dz	
	  bool do_brent=false;
	  double step1=mStepSizeZ;
	  double step2=mStepSizeZ;
	  if (k>=2){
	    step1=-forward_traj[k].z+forward_traj[k_minus_1].z;
	    step2=-forward_traj[k_minus_1].z+forward_traj[k-2].z;
	  }
	  //printf("step1 %f step 2 %f \n",step1,step2);
	  double two_step=step1+step2;
	  if (fabs(qBr2p*S(state_q_over_p)
		   *bfield->GetBz(S(state_x),S(state_y),z)
		   *two_step/sinl)<0.01 
	      && denom>EPS)
	    {
	    double dzw=z-z0w;
	    dz=-((S(state_x)-origin.X()-ux*dzw)*my_ux
		 +(S(state_y)-origin.Y()-uy*dzw)*my_uy)/denom;

	    if (fabs(dz)>two_step || dz<0){
	      do_brent=true;
	    }
	    else{
	      newz=z+dz;
	      // Check for exiting the straw
	      if (newz>endplate_z){
		newz=endplate_z;
		dz=endplate_z-z;
	      }
	      // Step the state and covariance through the field
	      if (dz>mStepSizeZ){
		double my_z=z;
		int my_steps=int(dz/mStepSizeZ);
		double dz2=dz-my_steps*mStepSizeZ;		     
		for (int m=0;m<my_steps;m++){
		  newz=my_z+mStepSizeZ;
	      
		  // Bail if the momentum has dropped below some minimum
		  if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		    if (DEBUG_LEVEL>2)
		      {
			_DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		      }	  
		    break_point_fdc_index=(3*num_fdc)/4;		    
		    return MOMENTUM_OUT_OF_RANGE;
		  }

		  // Step current state by step size 
		  Step(my_z,newz,dedx,S);
					       
		  my_z=newz;
		}
		newz=my_z+dz2;
		// Bail if the momentum has dropped below some minimum
		if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		  if (DEBUG_LEVEL>2)
		    {
		      _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		    }
		  break_point_fdc_index=(3*num_fdc)/4;
		  return MOMENTUM_OUT_OF_RANGE;
		}

		Step(my_z,newz,dedx,S);
	      }
	      else{
		// Bail if the momentum has dropped below some minimum
		if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		  if (DEBUG_LEVEL>2)
		    {
		      _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		    }
		  break_point_fdc_index=(3*num_fdc)/4;
		  return MOMENTUM_OUT_OF_RANGE;
		}
		
		Step(z,newz,dedx,S);
	      }
	    }
	  }
	  else do_brent=true;
	  if (do_brent){
	    // We have bracketed the minimum doca:  use Brent's agorithm
	    if (BrentsAlgorithm(z,-mStepSizeZ,dedx,z0w,origin,dir,S,dz)!=NOERROR){
	      break_point_fdc_index=(3*num_fdc)/4;
	      return MOMENTUM_OUT_OF_RANGE;
	    }
	    newz=z+dz;

	    if (fabs(dz)>2.*mStepSizeZ-EPS3){    
	      // whoops, looks like we didn't actually bracket the minimum 
	      // after all.  Swim to make sure we pass the minimum doca.
	      double ztemp=newz;
	      
	      // new wire position
	      wirepos=origin;
	      wirepos+=(ztemp-z0w)*dir;

	      // doca
	      old_doca2=doca2;

	      dx=S(state_x)-wirepos.X();
	      dy=S(state_y)-wirepos.Y();
	      doca2=dx*dx+dy*dy;
	      
	      while(doca2<old_doca2){
		newz=ztemp+mStepSizeZ;
		old_doca2=doca2;
		
		// Step to the new z position
		Step(ztemp,newz,dedx,S);

		// find the new distance to the wire
		wirepos=origin;
		wirepos+=(newz-z0w)*dir;
	      
		dx=S(state_x)-wirepos.X();
		dy=S(state_y)-wirepos.Y();
		doca2=dx*dx+dy*dy;
		
		ztemp=newz;
	      }
	      // Find the true doca
	      double dz2=0.;
	      if (BrentsAlgorithm(newz,mStepSizeZ,dedx,z0w,origin,dir,S,dz2)!=NOERROR){
		break_point_fdc_index=(3*num_fdc)/4;
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      newz=ztemp+dz2;
	   
	      // Change in z relative to where we started for this wire
	      dz=newz-z;
	    }
	    
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
	      //Step(z,newz,dedx,S);
	      
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
	    //Step(my_z,newz,dedx,S);	  
	    
	    // propagate error matrix to z-position of hit
	    StepJacobian(my_z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz3
	    Step(my_z,newz,dedx,S0); 	    
	  }
	  else{
	    // Step current state by dz
	    //Step(z,newz,dedx,S);
	    
	    // propagate error matrix to z-position of hit
	    StepJacobian(z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz
	    Step(z,newz,dedx,S0); 
	  }
	  
	  // Wire position at current z
	  wirepos=origin;
	  wirepos+=(newz-z0w)*dir;
	  
	  double xw=wirepos.X();
	  double yw=wirepos.Y();
	  
	  // predicted doca taking into account the orientation of the wire
	  dy=S(state_y)-yw;
	  dx=S(state_x)-xw;      
	  double cosstereo=my_cdchits[cdc_index]->cosstereo;
	  double d=sqrt(dx*dx+dy*dy)*cosstereo;
	  
	  // Track projection
	  double cosstereo2_over_d=cosstereo*cosstereo/d;
	  H_T(state_x)=dx*cosstereo2_over_d; 
	  H(state_x)=H_T(state_x);
	  H_T(state_y)=dy*cosstereo2_over_d;	  
	  H(state_y)=H_T(state_y);
      
	  //H.Print();
	  
	  // The next measurement
	  double dm=0.39,tdrift=0.,tcorr=0.;
	  if (fit_type==kTimeBased || USE_PASS1_TIME_MODE){
	    // Find offset of wire with respect to the center of the
	    // straw at this z position
	    const DCDCWire *mywire=my_cdchits[cdc_index]->hit->wire;
	    int ring_index=mywire->ring-1;
	    int straw_index=mywire->straw-1;
	    double dz=newz-z0w;
	    double phi_d=atan2(dy,dx);
	    double delta
	      =max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
	      *sin(phi_d+sag_phi_offset[ring_index][straw_index]);
	    double dphi=phi_d-mywire->origin.Phi();
	    while (dphi>M_PI) dphi-=2*M_PI;
	    while (dphi<-M_PI) dphi+=2*M_PI;
	    if (mywire->origin.Y()<0) dphi*=-1.;

	    tdrift=my_cdchits[cdc_index]->tdrift-mT0
	      -forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	    double B=forward_traj[k_minus_1].B;
	    ComputeCDCDrift(dphi,delta,tdrift,B,dm,Vc,tcorr);

	    //_DBG_ << "t " << tdrift << " d " << d << " delta " << delta << " dphi " << atan2(dy,dx)-mywire->origin.Phi() << endl;
	    
	    //_DBG_ << tcorr << " " << dphi << " " << dm << endl;
	    
	  }

	  // Residual
	  double res=dm-d;

	  // inverse variance including prediction
	  //double InvV1=1./(Vc+H*(C*H_T));
	  double Vproj=C.SandwichMultiply(H_T);
	  double InvV1=1./(Vc+Vproj);
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
	      bool skip_ring
		=(my_cdchits[cdc_index]->hit->wire->ring==RING_TO_SKIP);
	      // update covariance matrix and state vector
	      if (my_cdchits[cdc_index]->hit->wire->ring!=RING_TO_SKIP){
		C=Ctest;
		S+=res*K;
	      }

	      // Flag the place along the reference trajectory with hit id
	      forward_traj[k].h_id=1000+cdc_index;

	      // Store updated info related to this hit
	      double scale=(skip_ring)?1.:(1.-H*K); 
	      cdc_updates[cdc_index].tdrift=tdrift;
	      cdc_updates[cdc_index].tcorr=tcorr;
	      cdc_updates[cdc_index].residual=res*scale;
	      cdc_updates[cdc_index].variance=Vc;
	      cdc_updates[cdc_index].doca=dm;
	      cdc_updates[cdc_index].used_in_fit=true;
	    
	      // Update chi2 and number of degrees of freedom for this hit
	      if (skip_ring==false){
		chisq+=scale*res*res/Vc;
		numdof++;
	      }

	      if (DEBUG_LEVEL>10)
		cout << "Ring " <<  my_cdchits[cdc_index]->hit->wire->ring
		     << " Straw " <<  my_cdchits[cdc_index]->hit->wire->straw
		     << " Pred " << d << " Meas " << dm
		     << " Sigma meas " << sqrt(Vc)
		     << " t " << tcorr
		     << " Chi2 " << (1.-H*K)*res*res/Vc << endl;
	      
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
	  
	  cdc_updates[cdc_index].S=S;
	  cdc_updates[cdc_index].C=C;	  
	}

	// new wire origin and direction
	if (cdc_index>0){
	  cdc_index--;
	  origin=my_cdchits[cdc_index]->origin;
	  z0w=my_cdchits[cdc_index]->z0wire;
	  dir=my_cdchits[cdc_index]->dir;
	}
	else more_cdc_measurements=false;
      
	// Update the wire position
	wirepos=origin+(z-z0w)*dir;
	
	// new doca
	dx=S(state_x)-wirepos.X();
	dy=S(state_y)-wirepos.Y();
	doca2=dx*dx+dy*dy;
      }
      old_doca2=doca2;
    }
  }
  // Save final z position
  z_=forward_traj[forward_traj.size()-1].z;

  // The following code segment addes a fake point at a well-defined z position
  // that would correspond to a thin foil target.  It should not be turned on
  // for an extended target.
  if (ADD_VERTEX_POINT){
    double dz_to_target=TARGET_Z-z_;
    double my_dz=mStepSizeZ*(dz_to_target>0?1.:-1.);
    int num_steps=int(fabs(dz_to_target/my_dz));

    for (int k=0;k<num_steps;k++){
      double newz=z_+my_dz;
      // Step C along z
      StepJacobian(z_,newz,S,0.,J);
      
      //C=J*C*J.Transpose();
      C=C.SandwichMultiply(J);
      
      // Step S along z
      Step(z_,newz,0.,S);
      
      z_=newz;
    }

    // Step C along z
    StepJacobian(z_,TARGET_Z,S,0.,J);
    
    //C=J*C*J.Transpose();
    C=C.SandwichMultiply(J);
    
    // Step S along z
    Step(z_,TARGET_Z,0.,S);
    
    z_=TARGET_Z;

    // predicted doca taking into account the orientation of the wire
    double dy=S(state_y);
    double dx=S(state_x);      
    double d=sqrt(dx*dx+dy*dy);
    
    // Track projection
    double one_over_d=1./d;
    H_T(state_x)=dx*one_over_d; 
    H(state_x)=H_T(state_x);
    H_T(state_y)=dy*one_over_d;	  
    H(state_y)=H_T(state_y);
  	  
    // Variance of target point
    // Variance is for average beam spot size assuming triangular distribution
    // out to 2.2 mm from the beam line.
    //   sigma_r = 2.2 mm/ sqrt(18)
    Vc=0.002689;

    // inverse variance including prediction
    //double InvV1=1./(Vc+H*(C*H_T));
    double InvV1=1./(Vc+C.SandwichMultiply(H_T));
    if (InvV1<0.){
      if (DEBUG_LEVEL>0)
	_DBG_ << "Negative variance???" << endl;
      return NEGATIVE_VARIANCE;
    }
    // Compute KalmanSIMD gain matrix
    K=InvV1*(C*H_T);
    
    // Update the state vector with the target point
    // "Measurement" is average of expected beam spot size
    double res=0.1466666667-d;
    S+=res*K;  
    // Update state vector covariance matrix
    //C=C-K*(H*C);    
    C=C.SubSym(K*(H*C));
    
    // Update chi2 for this segment
    chisq+=(1.-H*K)*res*res/Vc;
    numdof++;
  }
  
  // Check that there were enough hits to make this a valid fit
  if (numdof<6){
    chisq=MAX_CHI2;
    numdof=0;
   
    return INVALID_FIT;
  }

  //  chisq*=anneal_factor;
  numdof-=5;

  // Final positions in x and y for this leg
  x_=S(state_x);
  y_=S(state_y);

  if (DEBUG_LEVEL>1){
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
    cout << "Momentum " << 1./S(state_q_over_p) <<endl;
  }

  // Check if we have a kink in the track or threw away too many hits
  if (num_cdc>0 && break_point_fdc_index>0){ 
    if (break_point_fdc_index+num_cdc<MIN_HITS_FOR_REFIT){
      //_DBG_ << endl;
      unsigned int new_index=num_fdc/2;
      if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
	break_point_fdc_index=new_index;
      }
      else{
	break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
      }
    }
    return BREAK_POINT_FOUND;
  }
  if (num_cdc==0 && break_point_fdc_index>2){
    //_DBG_ << endl;
    if (break_point_fdc_index<num_fdc/2){
      break_point_fdc_index=(3*num_fdc)/4;
    }
    if (break_point_fdc_index<MIN_HITS_FOR_REFIT-1){
      break_point_fdc_index=MIN_HITS_FOR_REFIT-1;
    }
    return BREAK_POINT_FOUND;
  }
  if (num_cdc>5 && break_point_cdc_index>2){
    //_DBG_ << endl;  
    unsigned int new_index=num_fdc/2;
    if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
      break_point_fdc_index=new_index;
    }
    else{
      break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
    }
    return BREAK_POINT_FOUND;
  }
  unsigned int num_good=0; 
  unsigned int num_hits=num_cdc+max_num_fdc_used_in_fit;
  for (unsigned int j=0;j<num_cdc;j++){
    if (cdc_updates[j].used_in_fit) num_good++;
  }
  for (unsigned int j=0;j<num_fdc;j++){
    if (fdc_updates[j].used_in_fit) num_good++;
  }
  if (double(num_good)/double(num_hits)<MINIMUM_HIT_FRACTION){
    //_DBG_ <<endl;
    if (num_cdc==0){
      unsigned int new_index=(3*num_fdc)/4;
      break_point_fdc_index=(new_index>=MIN_HITS_FOR_REFIT)?new_index:(MIN_HITS_FOR_REFIT-1);
    }
    else{
      unsigned int new_index=num_fdc/2;
      if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
	break_point_fdc_index=new_index;
      }
      else{
	break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
      }
    }
    return PRUNED_TOO_MANY_HITS;
  }
    
  return FIT_SUCCEEDED;
}

// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD_ALT1::SmoothForward(void){ 
  if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;
  
  unsigned int max=forward_traj.size()-1;
  DMatrix5x1 S=(forward_traj[max].Skk);
  DMatrix5x5 C=(forward_traj[max].Ckk);
  DMatrix5x5 JT=(forward_traj[max].JT);
  DMatrix5x1 Ss=S;
  DMatrix5x5 Cs=C;
  DMatrix5x5 A,dC;

  for (unsigned int m=max-1;m>0;m--){
    if (forward_traj[m].h_id>0){
      if (forward_traj[m].h_id<1000){
	unsigned int id=forward_traj[m].h_id-1;
	A=fdc_updates[id].C*JT*C.InvertSym();
	Ss=fdc_updates[id].S+A*(Ss-S);
	
	if (!finite(Ss(state_q_over_p))){ 
	  if (DEBUG_LEVEL>5) _DBG_ << "Invalid values for smoothed parameters..." << endl;
	  return VALUE_OUT_OF_RANGE;
	}
	dC=A*(Cs-C)*A.Transpose();
	Cs=fdc_updates[id].C+dC;
	
	double cosa=my_fdchits[id]->cosa;
	double sina=my_fdchits[id]->sina;
	double u=my_fdchits[id]->uwire;
	double v=my_fdchits[id]->vstrip;
	
	// Position and direction from state vector
	double x=Ss(state_x);
	double y=Ss(state_y);
	double tx=Ss(state_tx);
	double ty=Ss(state_ty);
	
	// Projected position along the wire without doca-dependent corrections
	double vpred_uncorrected=x*sina+y*cosa;
	
	// Projected position in the plane of the wires transverse to the wires
	double upred=x*cosa-y*sina;
	
	// Direction tangent in the u-z plane
	double tu=tx*cosa-ty*sina;
	double alpha=atan(tu);
	double cosalpha=cos(alpha);
	//double cosalpha2=cosalpha*cosalpha;
	double sinalpha=sin(alpha);
	
	// (signed) distance of closest approach to wire
	double doca=(upred-u)*cosalpha;

	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;
	
	// Difference between measurement and projection
	double tv=tx*sina+ty*cosa;
	double resi=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
					       -tv*sinalpha));
	
	// Variance from filter step
	double V=fdc_updates[id].variance;
	// Compute projection matrix and find the variance for the residual
	DMatrix5x1 H_T;
	double temp2=nz_sinalpha_plus_nr_cosalpha-tv*sinalpha;
	H_T(state_x)=sina+cosa*cosalpha*temp2;	
	H_T(state_y)=cosa-sina*cosalpha*temp2;	
       
	double cos2_minus_sin2=cosalpha*cosalpha-sinalpha*sinalpha;
	double fac=nz*cos2_minus_sin2-2.*nr*cosalpha*sinalpha;
	double doca_cosalpha=doca*cosalpha;
	double temp=doca_cosalpha*fac;	
	H_T(state_tx)=cosa*temp
	  -doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
	  ;
	H_T(state_ty)=-sina*temp
	  -doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
	  ;

	if (my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP){
	  V+=Cs.SandwichMultiply(H_T);
	}
	else{
	  V-=dC.SandwichMultiply(H_T);
	}

	pulls.push_back(pull_t(resi,sqrt(V),
			       forward_traj[m].s,
			       fdc_updates[id].tdrift,
			       fdc_updates[id].doca,
			       NULL,my_fdchits[id]->hit,
			       forward_traj[m].z));
      }
      else{
	unsigned int id=forward_traj[m].h_id-1000;
	A=cdc_updates[id].C*JT*C.InvertSym();
	Ss=cdc_updates[id].S+A*(Ss-S);

	if (!finite(Ss(state_q_over_p))){
	  if (DEBUG_LEVEL>5) _DBG_ << "Invalid values for smoothed parameters..." << endl;
	  return VALUE_OUT_OF_RANGE;
	}

	Cs=cdc_updates[id].C+A*(Cs-C)*A.Transpose();
	
	// Fill in pulls information for cdc hits
	FillPullsVectorEntry(Ss,Cs,forward_traj[m],my_cdchits[id],
			     cdc_updates[id]);
      }
    }
    else{
      A=forward_traj[m].Ckk*JT*C.InvertSym();
      Ss=forward_traj[m].Skk+A*(Ss-S);
      Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    
    S=forward_traj[m].Skk;
    C=forward_traj[m].Ckk;
    JT=forward_traj[m].JT;
  }

  return NOERROR;
}
