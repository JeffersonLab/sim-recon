/*
 * hitFDC - registers hits for forward drift chambers
 *
 *        This is a part of the hits package for the
 *        HDGeant simulation program for Hall D.
 *
 *        version 1.0         -Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitForwardDC
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#include "calibDB.h"
extern s_HDDM_t* thisInputEvent;

const float Tau[] = {0,-45,0,45,15,60,105,-105,-60,-15};
//const float wire_dead_zone_radius[4]={2.3,3.2,3.9,4.6};
//const float strip_dead_zone_radius[4]={2.3,3.2,3.9,4.6};
const float wire_dead_zone_radius[4]={3.0,3.0,3.0,3.9};
const float strip_dead_zone_radius[4]={1.3,1.3,1.3,1.3};

#define CATHODE_ROT_ANGLE 1.309 // 75 degrees

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
static float DRIFT_SPEED           =.0055;
static float ACTIVE_AREA_OUTER_RADIUS =48.5;
static float ANODE_CATHODE_SPACING =0.5;
static float TWO_HIT_RESOL         =25.;
static int   WIRES_PER_PLANE       =96;
static float WIRE_SPACING          =1.0;
static float U_OF_WIRE_ZERO        =0;//(-((WIRES_PER_PLANE-1.)*WIRE_SPACING)/2)
static float STRIPS_PER_PLANE      =192;
static float STRIP_SPACING         =0.5;
static float U_OF_STRIP_ZERO	   =0;//  (-((STRIPS_PER_PLANE-1.)*STRIP_SPACING)/2)
static float STRIP_GAP             =0.1;
static int   MAX_HITS             =100;
static float K2                  =1.15;
static float STRIP_NODES          = 3;
static float THRESH_KEV           =1. ;
static float THRESH_STRIPS        =5. ;  /* pC */
static float ELECTRON_CHARGE =1.6022e-4; /* fC */
static float DIFFUSION_COEFF  =   1.1e-6; // cm^2/s --> 200 microns at 1 cm

/* The folowing are for interpreting grid of Lorentz deflection data */
#define PACKAGE_Z_POINTS 10
#define LORENTZ_X_POINTS 21
#define LORENTZ_Z_POINTS 4*PACKAGE_Z_POINTS

binTree_t* forwardDCTree = 0;
static int stripCount = 0;
static int wireCount = 0;
static int pointCount = 0;
static int initialized=0;
static int initializedx=0;

// Variables for implementing lorentz effect (deflections of avalanche position
// due to the magnetic field).
float lorentz_x[LORENTZ_X_POINTS];
float lorentz_z[LORENTZ_Z_POINTS];
float *lorentz_nx[LORENTZ_X_POINTS];
float *lorentz_nz[LORENTZ_X_POINTS];

// Array of parameters for drift time smearing
float fdc_smear_parms[9*26];

void gpoiss_(float*,int*,const int*); // avoid solaris compiler warnings
void rnorml_(float*,int*);

// Locate a position in array xx given x
void locate(float *xx,int n,float x,int *j){
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if (x>=xx[jm]==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x==xx[0]) *j=0;
  else if (x==xx[n-1]) *j=n-2;
  else *j=jl; 
}

// Polynomial interpolation on a grid.
// Adapted from Numerical Recipes in C (2nd Edition), pp. 121-122.
void polint(float *xa, float *ya,int n,float x, float *y,float *dy){
  int i,m,ns=0;
  float den,dif,dift,ho,hp,w;

  float *c=(float *)calloc(n,sizeof(float));
  float *d=(float *)calloc(n,sizeof(float));

  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++){
    if ((dift=fabs(x-xa[i]))<dif){
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];

  for (m=1;m<n;m++){
    for (i=1;i<=n-m;i++){
      ho=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i+1-1]-d[i-1];
      if ((den=ho-hp)==0.0) return;
      
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
      
    }
    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
}

// Function for modelling the variation of the drift distance (or time)
float drift_smear_function(float x,float parms[]){
  float sum=0;
  unsigned int i;
  for (i=0;i<3;i++){
    float sigma=parms[6+i];
    if (sigma>0.){
      float norm_xdiff=(x-parms[3+i])/sigma;
      sum+=parms[i]*exp(-0.5*norm_xdiff*norm_xdiff)/(sigma*sqrt(2.*M_PI));
    }
  }
  return sum;
}

// Rejection method for determining amount of smearing of the drift distance
float get_drift_smear(float parms[]){
  float smear,smear_max=2.0;
  float f,f_prime,f_max=drift_smear_function(0.,parms);
  int goodf=0;
  float rndno[2];
  int two=2;

  while (!goodf){
    grndm_(rndno,&two);
    smear=smear_max*(rndno[0]-0.5);
    f=drift_smear_function(smear,parms);
    f_prime=f_max*rndno[1];
    if (f>f_prime) goodf=1;
  }

  return smear;
}



/* register hits during tracking (from gustep) */

void hitForwardDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum,
                   int track, int stack, int history, int ipart)
{
  float x[3], t;
  float dx[3], dr;
  float dEdx;
  float xlocal[3];
  float xinlocal[3];
  float xoutlocal[3];
  float dradius;
  float alpha;
  int i,j;

  if (!initializedx){
      mystr_t strings[250];
      float values[250];
      int nvalues = 250;
      // Get the parameters for the drift time smearing from the data base
      int status = GetArrayConstants("FDC/drift_smear_parms", &nvalues, 
				     fdc_smear_parms,
				     strings); 
      // Get other parameters related to the geometry and the signals 
      status = GetConstants("FDC/fdc_parms", &nvalues,values,strings);
      if (!status) {
        int ncounter = 0;
        int i;
        for ( i=0;i<(int)nvalues;i++){
          //printf("%d %s %f\n",i,strings[i].str,values[i]);
          if (!strcmp(strings[i].str,"FDC_DRIFT_SPEED")) {
            DRIFT_SPEED  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_ACTIVE_AREA_OUTER_RADIUS")) {
            ACTIVE_AREA_OUTER_RADIUS  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_ANODE_CATHODE_SPACING")) {
            ANODE_CATHODE_SPACING  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_TWO_HIT_RESOL")) {
            TWO_HIT_RESOL  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_WIRES_PER_PLANE")) {
            WIRES_PER_PLANE  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_WIRE_SPACING")) {
            WIRE_SPACING  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_STRIPS_PER_PLANE")) {
            STRIPS_PER_PLANE  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_STRIP_SPACING")) {
            STRIP_SPACING  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_STRIP_GAP")) {
            STRIP_GAP  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_MAX_HITS")) {
            MAX_HITS  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_K2")) {
            K2  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_STRIP_NODES")) {
            STRIP_NODES  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_THRESH_KEV")) {
            THRESH_KEV  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_THRESH_STRIPS")) {
            THRESH_STRIPS  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_ELECTRON_CHARGE")) {
            ELECTRON_CHARGE  = values[i];
            ncounter++;
          }
          if (!strcmp(strings[i].str,"FDC_DIFFUSION_COEFF")) {
            DIFFUSION_COEFF  = values[i];
            ncounter++;
          }
	}
	U_OF_WIRE_ZERO     = (-((WIRES_PER_PLANE-1.)*WIRE_SPACING)/2);
	U_OF_STRIP_ZERO	   = (-((STRIPS_PER_PLANE-1.)*STRIP_SPACING)/2);
  
	if (ncounter==16){
          printf("FDC: ALL parameters loaded from Data Base\n");
        } else if (ncounter<16){
          printf("FDC: NOT ALL necessary parameters found in Data Base %d out of 16\n",ncounter);
        } else {
          printf("FDC: SOME parameters found more than once in Data Base\n");
        }       
      }
      initializedx = 1 ;
  }

  /* Get chamber information */ 
  int layer = getlayer_(); 
  if (layer==0){
    printf("hitFDC: FDC layer number evaluates to zero! THIS SHOULD NEVER HAPPEN! drop this particle.\n");
    return;
  }
  //int module = getmodule_();
  //int chamber = (module*10)+layer;
  //int PackNo = (chamber-11)/20;
  int PackNo = getpackage_()-1;
  if (PackNo==-1){
     printf("hitFDC: FDC package number evaluates to zero! THIS SHOULD NEVER HAPPEN! drop this particle.\n");
    return;
  }
  int module = 2*(PackNo)+(layer-1)/3+1;
  int chamber = (module*10)+(layer-1)%3+1;
  int wire1,wire2;
  int wire,dwire;

  // transform layer number into Richard's scheme
  layer=(layer-1)%3+1;

  // Initialize arrays of deflection data from the Lorentz effect
  if (!initialized){
		// Allocate memory for 2-D tables. By "faking" a 2-D array this way,
		// we can pass the pointers into GetLorentzDeflections() without
		// having to hardwire array sizes into the data types of the arguments.
		for (i=0;i<LORENTZ_X_POINTS;i++){
			lorentz_nx[i] = (float*)malloc(LORENTZ_Z_POINTS*sizeof(float));
			lorentz_nz[i] = (float*)malloc(LORENTZ_Z_POINTS*sizeof(float));
		}
		
		// Get tables from database
		GetLorentzDeflections(lorentz_x, lorentz_z, lorentz_nx, lorentz_nz, LORENTZ_X_POINTS, LORENTZ_Z_POINTS);
		initialized=1;
  }  

  transformCoord(xin,"global",xinlocal,"local");

  wire1 = ceil((xinlocal[0] - U_OF_WIRE_ZERO)/WIRE_SPACING +0.5);
  transformCoord(xout,"global",xoutlocal,"local");
  wire2 = ceil((xoutlocal[0] - U_OF_WIRE_ZERO)/WIRE_SPACING +0.5);
  // Check that wire numbers are not out of range
  if ((wire1>WIRES_PER_PLANE && wire2==WIRES_PER_PLANE) ||
      (wire2>WIRES_PER_PLANE && wire1==WIRES_PER_PLANE)) 
    wire1=wire2=WIRES_PER_PLANE;  
  if ((wire1==0 && wire2 == 1) || (wire1==1 && wire2== 0)){
    wire1=wire2=1;
  }
  // Make sure at least one wire number is valid
  if (wire1>WIRES_PER_PLANE&&wire2>WIRES_PER_PLANE) return;
  if (wire1==0 && wire2==0) return;
  dwire = (wire1 < wire2)? 1 : -1;
  alpha = atan2(xoutlocal[0]-xinlocal[0],xoutlocal[2]-xinlocal[2]);
  xlocal[0] = (xinlocal[0] + xoutlocal[0])/2;
  xlocal[1] = (xinlocal[1] + xoutlocal[1])/2;
  xlocal[2] = (xinlocal[2] + xoutlocal[2])/2;
  wire = ceil((xlocal[0] - U_OF_WIRE_ZERO)/WIRE_SPACING +0.5);
  x[0] = (xin[0] + xout[0])/2;
  x[1] = (xin[1] + xout[1])/2;
  x[2] = (xin[2] + xout[2])/2;
  t = (xin[3] + xout[3])/2 * 1e9;
  dx[0] = xin[0] - xout[0];
  dx[1] = xin[1] - xout[1];
  dx[2] = xin[2] - xout[2];
  dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
  if (dr > 1e-3)
  {
    dEdx = dEsum/dr;
  }
  else
  {
    dEdx = 0;
  }

  /* Make a fuzzy boundary around the forward dead region 
   * by killing any track segment whose midpoint is within the boundary */

  if (sqrt(xlocal[0]*xlocal[0]+xlocal[1]*xlocal[1])
      < wire_dead_zone_radius[PackNo])
  {
    return;
  }

  /* post the hit to the truth tree */
 
  if (history == 0)
  {
    int mark = (1<<16) + (chamber<<20) + pointCount;
    void** twig = getTwig(&forwardDCTree, mark);
    if (*twig == 0)
    {
      s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
      s_FdcChambers_t* chambers = make_s_FdcChambers(1);
      s_FdcTruthPoints_t* points = make_s_FdcTruthPoints(1);
      float xwire = U_OF_WIRE_ZERO + (wire-1)*WIRE_SPACING;
      float u[2];
      u[0] = xinlocal[2];
      u[1] = xinlocal[0]-xwire;
      dradius = fabs(u[1]*cos(alpha)-u[0]*sin(alpha));
      points->mult = 1;
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
      points->in[0].primary = (stack <= a);
      points->in[0].track = track;
      points->in[0].x = x[0];
      points->in[0].y = x[1];
      points->in[0].z = x[2];
      points->in[0].t = t;
      points->in[0].px = pin[0]*pin[4];
      points->in[0].py = pin[1]*pin[4];
      points->in[0].pz = pin[2]*pin[4];
      points->in[0].E = pin[3];
      points->in[0].dradius = dradius;
      points->in[0].dEdx = dEdx;
      points->in[0].ptype = ipart;
      chambers->mult = 1;
      chambers->in[0].module = module;
      chambers->in[0].layer = layer;
      chambers->in[0].fdcTruthPoints = points;
      fdc->fdcChambers = chambers;
      pointCount++;
    }
  }

  /* post the hit to the hits tree, mark cell as hit */

  if (dEsum > 0)
  {
    int nhit;
    s_FdcAnodeTruthHits_t* ahits;    
    s_FdcCathodeTruthHits_t* chits;    
    float tdrift,tdrift_unsmeared;    
    // Array of drift distance grid points needed for drift time smearing
    float xarray[26]={0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,
		      0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,
		      0.48,0.50};
    // array of parameters for drift time smearing
    float my_parms[9];
    int m=0,ind;


    for (wire=wire1; wire-dwire != wire2; wire+=dwire)
    {
      int valid_hit=1;
      float dE,dt;
      float u[2];
      float x0[3],x1[3];
      float avalanche_y;
      float xwire = U_OF_WIRE_ZERO + (wire-1)*WIRE_SPACING;
         
      x0[0] = xwire-0.5*dwire*WIRE_SPACING;
      x0[1] = xinlocal[1] + (x0[0]-xinlocal[0]+1e-20)*
             (xoutlocal[1]-xinlocal[1])/(xoutlocal[0]-xinlocal[0]+1e-20);
      x0[2] = xinlocal[2] + (x0[0]-xinlocal[0]+1e-20)*
             (xoutlocal[2]-xinlocal[2])/(xoutlocal[0]-xinlocal[0]+1e-20);
      if (fabs(x0[2]-xoutlocal[2]) > fabs(xinlocal[2]-xoutlocal[2]))
      {
        x0[0] = xinlocal[0];
        x0[1] = xinlocal[1];
        x0[2] = xinlocal[2];
      }
      x1[0] = xwire+0.5*dwire*WIRE_SPACING;
      x1[1] = xinlocal[1] + (x1[0]-xinlocal[0]+1e-20)*
             (xoutlocal[1]-xinlocal[1])/(xoutlocal[0]-xinlocal[0]+1e-20);
      x1[2] = xinlocal[2] + (x1[0]-xinlocal[0]+1e-20)*
             (xoutlocal[2]-xinlocal[2])/(xoutlocal[0]-xinlocal[0]+1e-20);
      if (fabs(x1[2]-xinlocal[2]) > fabs(xoutlocal[2]-xinlocal[2]))
      {
        x1[0] = xoutlocal[0];
        x1[1] = xoutlocal[1];
        x1[2] = xoutlocal[2];
      }
      u[0] = xinlocal[2];
      u[1] = xinlocal[0]-xwire;
      dradius = fabs(u[1]*cos(alpha)-u[0]*sin(alpha));
      tdrift_unsmeared = t + dradius/DRIFT_SPEED;
      dE = dEsum*(x1[2]-x0[2])/(xoutlocal[2]-xinlocal[2]);

      /* Simulate smearing of drift time due to diffusion in the gas */
      // First interpolate over the grid of smearing parameters 
      ind=(dradius<0.5?(int)floor(dradius/0.02):25);
      for (m=0;m<9;m++){
	if (ind<25){
	  float dpar;
	  float parlist[26];
	  int  num=3;
	  int p;
	  if (ind>=23){
	    ind=23;
	  }
	  for (p=0;p<num;p++){
	    parlist[p]=fdc_smear_parms[9*(ind+p)+m];
	  }
	  polint(&xarray[ind],parlist,num,dradius,&my_parms[m],&dpar);
	}
	// .. avoid doing an extrapolation...
	else my_parms[m]=fdc_smear_parms[9*ind+m];
      }
      dt=get_drift_smear(my_parms)/DRIFT_SPEED;
      tdrift=tdrift_unsmeared+dt;

      /* Compute Lorentz deflection of avalanche position */
      /* Checks to see if deflection would push avalanche out of active area */
      if (dE > 0){
	float tanr,tanz,r;
	float dist_to_wire=xlocal[0]-xwire;
	float phi=atan2(x0[1]+x1[1],x0[0]+x1[0]);
	float ytemp[LORENTZ_X_POINTS],ytemp2[LORENTZ_X_POINTS],dy;
	int imin,imax,ind,ind2;
	float rndno[2];
	int two=2;
	
	avalanche_y = (x0[1]+x1[1])/2;
	r=sqrt((x0[0]+x1[0])*(x0[0]+x1[0])+(x0[1]+x1[1])*(x0[1]+x1[1]))/2.;

	// Locate positions in x and z arrays given r and z=x[2]
	locate(lorentz_x,LORENTZ_X_POINTS,r,&ind);
	locate(lorentz_z,LORENTZ_Z_POINTS,x[2],&ind2);
	
	// First do interpolation in z direction 
	imin=PACKAGE_Z_POINTS*(ind2/PACKAGE_Z_POINTS); // Integer division...
	for (j=0;j<LORENTZ_X_POINTS;j++){
	  polint(&lorentz_z[imin],&lorentz_nx[j][imin],PACKAGE_Z_POINTS,x[2],
		 &ytemp[j],&dy);
	  polint(&lorentz_z[imin],&lorentz_nz[j][imin],PACKAGE_Z_POINTS,x[2],
		 &ytemp2[j],&dy);
	}
	// Then do final interpolation in x direction 
	if (ind>0){
	  imin=((ind+3)>LORENTZ_X_POINTS)?(LORENTZ_X_POINTS-4):(ind-1);
	}
	else imin=0;
	polint(&lorentz_x[imin],&ytemp[imin],4,r,&tanr,&dy);
	polint(&lorentz_x[imin],&ytemp2[imin],4,r,&tanz,&dy);


	// Correct avalanche position with deflection along wire	
 	avalanche_y+=-tanz*dist_to_wire*sin(alpha)*cos(phi)
	  +tanr*dist_to_wire*cos(alpha);

	// Only and always use the built-in Geant random generator,
        // otherwise debugging is a problem because sequences are not
        // reproducible from a given pair of random seeds. [rtj]

        /* rnorml_(rndno,&two); */ {
           float rho,phi;
           grndm_(rndno,&two);
           rho = sqrt(-2*log(rndno[0]));
           phi = rndno[1]*2*M_PI;
           rndno[0] = rho*cos(phi);
           rndno[1] = rho*sin(phi);
        }

	// Angular dependence from PHENIX data
	avalanche_y+=0.16*ANODE_CATHODE_SPACING*tan(alpha)*rndno[0];
	// Crude approximation for transverse diffusion
	avalanche_y+=sqrt(2.*DIFFUSION_COEFF*dradius/DRIFT_SPEED)*rndno[1];

	/* If the Lorentz effect would deflect the avalanche out of the active 
	   region, mark as invalid hit */
	r=sqrt(avalanche_y*avalanche_y+xwire*xwire);
	if (r>ACTIVE_AREA_OUTER_RADIUS || r<wire_dead_zone_radius[PackNo]){
	  valid_hit=0;
	}
      }

    /* first record the anode wire hit */

      if (dE > 0 && valid_hit)
      {
        int mark = (chamber<<20) + (2<<10) + wire;
        void** twig = getTwig(&forwardDCTree, mark);

        if (*twig == 0)
        {
          s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
          s_FdcChambers_t* chambers = make_s_FdcChambers(1);
          s_FdcAnodeWires_t* wires = make_s_FdcAnodeWires(1);
          wires->mult = 1;
          wires->in[0].wire = wire;
          wires->in[0].fdcAnodeTruthHits = ahits = make_s_FdcAnodeTruthHits(MAX_HITS);
          chambers->mult = 1;
	  chambers->in[0].module = module;
          chambers->in[0].layer = layer;
          chambers->in[0].fdcAnodeWires = wires;
          fdc->fdcChambers = chambers;
          wireCount++;          
        }
        else
        {
           s_ForwardDC_t* fdc = *twig;
           ahits = fdc->fdcChambers->in[0].fdcAnodeWires->in[0].fdcAnodeTruthHits;
        }

        for (nhit = 0; nhit < ahits->mult; nhit++)
        {
          if (fabs(ahits->in[nhit].t - tdrift) < TWO_HIT_RESOL)
          {
            break;
          }
        }
        if (nhit < ahits->mult)                 /* merge with former hit */
        {
				/* keep the earlier hit and discard the later one */
				/* Feb. 11, 2008 D. L. */
				if(ahits->in[nhit].t>tdrift){
					ahits->in[nhit].t = tdrift;
					ahits->in[nhit].t_unsmeared=tdrift_unsmeared;
					ahits->in[nhit].d = dradius;
					ahits->in[nhit].dE = dEsum;
					ahits->in[nhit].itrack = track;
					ahits->in[nhit].ptype = ipart;
				}
			
          /*ahits->in[nhit].t = 
              (ahits->in[nhit].t * ahits->in[nhit].dE + tdrift * dE)
              / (ahits->in[nhit].dE += dE);
			 */
        }
        else if (nhit < MAX_HITS)              /* create new hit */
        {
          ahits->in[nhit].t = tdrift;
	  ahits->in[nhit].t_unsmeared=tdrift_unsmeared;
          ahits->in[nhit].dE = dE;
	  ahits->in[nhit].d = dradius;
	  ahits->in[nhit].itrack = track;
	  ahits->in[nhit].ptype = ipart;
          ahits->mult++;
        }
        else
        {
          fprintf(stderr,"HDGeant error in hitForwardDC: ");
          fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
        }
      }

    /* then generate hits in the two surrounding cathode planes */

      if (dE > 0 && valid_hit)
      {
        float avalanche_x = xwire;
  
        /* Mock-up of cathode strip charge distribution */
  
        int plane, node;
        for (plane=1; plane<4; plane+=2)
        {
          float theta = (plane == 1)? -CATHODE_ROT_ANGLE : +CATHODE_ROT_ANGLE;
          float cathode_u = avalanche_x*cos(theta)+avalanche_y*sin(theta);
          int strip1 = ceil((cathode_u - U_OF_STRIP_ZERO)/STRIP_SPACING +0.5);
          float cathode_u1 = (strip1-1)*STRIP_SPACING + U_OF_STRIP_ZERO;
          float delta = cathode_u-cathode_u1;
  
          /* Variables for approximating number of ion pairs  */ 
          const float w_eff=26.0; // eV, average energy needed to produce an 
                                  // ion pair.  For simplicity use pure argon
          const float n_s_per_p=2.2; // average number of secondary ion pairs
                                     // per primary ionization
          float n_p_mean,n_s_mean;
          int n_p,n_t,n_s;
          const int one=1;
      
          /* variables for gain approximation */
          const float k=1.81;   // empirical constant related to first
                                // Townsend coefficient of argon gas
          const float N=269.*273/293;  // proportional to number of gas 
                                       // molecules/cm^3 
          const float V=1800;   // V, operation voltage 
          const float V_t=703.5; // V, threshold voltage
          const float amp_gain=46; // 2.3 mV/uA into 50 Ohms
          const float C=0.946;  // capacitance in units of episilon_0 
          const float a=0.001;  // cm, radius of sense wires 
  
          /* Sauli eq. 31: */
          float gain = exp(sqrt(2.*k*V*C*N*a/M_PI)*(sqrt(V/V_t)-1.));
          float q_anode;
  
          /* Total number of ion pairs.  On average for each primary ion 
             pair produced there are n_s secondary ion pairs produced.  The
             probability distribution is a compound poisson distribution
             that requires generating two Poisson variables.
           */
          n_p_mean = dE/w_eff/(1.+n_s_per_p)*1e9;
          gpoiss_(&n_p_mean,&n_p,&one);
          n_s_mean = ((float)n_p)*n_s_per_p;
          gpoiss_(&n_s_mean,&n_s,&one);
          n_t = n_s+n_p;
          q_anode=((float)n_t)*gain*ELECTRON_CHARGE;
  
          for (node=-STRIP_NODES; node<=STRIP_NODES; node++)
          {
          /* Induce charge on the strips according to the Mathieson 
             function tuned to results from FDC prototype
           */
            float lambda1=(((float)node-0.5)*STRIP_SPACING+STRIP_GAP/2.
                          -delta)/ANODE_CATHODE_SPACING;
            float lambda2=(((float)node+0.5)*STRIP_SPACING-STRIP_GAP/2.
                           -delta)/ANODE_CATHODE_SPACING;
            float q = q_anode*amp_gain/1000.
                      *(tanh(M_PI*K2*lambda2/4.)-tanh(M_PI*K2*lambda1/4.))/4.;
  
            int strip = strip1+node;
	    /* Throw away hits on strips falling within a certain dead-zone
	       radius */
	    float strip_outer_u=cathode_u1
	      +(STRIP_SPACING+STRIP_GAP/2.)*(int)node;
	    float cathode_v=-avalanche_x*sin(theta)+avalanche_y*cos(theta);
	    float check_radius=sqrt(strip_outer_u*strip_outer_u
				    +cathode_v*cathode_v);

            if ((strip > 0) && (check_radius>strip_dead_zone_radius[PackNo]) 
		&& (strip <= STRIPS_PER_PLANE))
            {
  	      int mark = (chamber<<20) + (plane<<10) + strip;
              void** twig = getTwig(&forwardDCTree, mark);
              if (*twig == 0)
              {
                s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
                s_FdcChambers_t* chambers = make_s_FdcChambers(1);
                s_FdcCathodeStrips_t* strips = make_s_FdcCathodeStrips(1);
                strips->mult = 1;
                strips->in[0].plane = plane;
                strips->in[0].strip = strip;
                strips->in[0].fdcCathodeTruthHits = chits
                                             = make_s_FdcCathodeTruthHits(MAX_HITS);
                chambers->mult = 1;
                chambers->in[0].module = module;
                chambers->in[0].layer = layer;
                chambers->in[0].fdcCathodeStrips = strips;
                fdc->fdcChambers = chambers;
                stripCount++;
              }
              else
              {
                s_ForwardDC_t* fdc = *twig;
                chits = fdc->fdcChambers->in[0].fdcCathodeStrips
                                        ->in[0].fdcCathodeTruthHits;
              }
          
              for (nhit = 0; nhit < chits->mult; nhit++)
              {
                if (fabs(chits->in[nhit].t - tdrift) < TWO_HIT_RESOL)
                {
                  break;
                }
              }
              if (nhit < chits->mult)          /* merge with former hit */
              {
						/* keep the earlier hit and discard the later one */
						/* Feb. 11, 2008 D. L. */
						if(chits->in[nhit].t>tdrift){
							chits->in[nhit].t = tdrift;
							chits->in[nhit].q = q;
							chits->in[nhit].itrack = track;
							chits->in[nhit].ptype = ipart;
						}
                /*chits->in[nhit].t = 
                    (chits->in[nhit].t * chits->in[nhit].q + tdrift * q)
                    / (chits->in[nhit].q += q);
					*/
              }
              else if (nhit < MAX_HITS)        /* create new hit */
              {
                chits->in[nhit].t = tdrift;
                chits->in[nhit].q = q;
		chits->in[nhit].itrack = track;
		chits->in[nhit].ptype = ipart;
                chits->mult++;
              }
              else
              {
                fprintf(stderr,"HDGeant error in hitForwardDC: ");
                fprintf(stderr,"max hit count %d exceeded, truncating!\n",
                        MAX_HITS);
              }
            }
          }
        }
      }
    }
  }
}

/* entry points from fortran */

void hitforwarddc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitForwardDC(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}


/* pick and package the hits for shipping */

s_ForwardDC_t* pickForwardDC ()
{
   s_ForwardDC_t* box;
   s_ForwardDC_t* item;

   if ((stripCount == 0) && (wireCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_ForwardDC();
   box->fdcChambers = make_s_FdcChambers(32);
   box->fdcChambers->mult = 0;
   while (item = (s_ForwardDC_t*) pickTwig(&forwardDCTree))
   {
      s_FdcChambers_t* chambers = item->fdcChambers;
      int module = chambers->in[0].module;
      int layer = chambers->in[0].layer;
      int m = box->fdcChambers->mult;

      /* compress out the hits below threshold */
      s_FdcAnodeWires_t* wires = chambers->in[0].fdcAnodeWires;
      int wire;
      s_FdcCathodeStrips_t* strips = chambers->in[0].fdcCathodeStrips;
      int strip;
      s_FdcTruthPoints_t* points = chambers->in[0].fdcTruthPoints;
      int point;
      int mok=0;
      for (wire=0; wire < wires->mult; wire++)
      {
         s_FdcAnodeTruthHits_t* ahits = wires->in[wire].fdcAnodeTruthHits;

         int i,iok;
         for (iok=i=0; i < ahits->mult; i++)
         {
            if (ahits->in[i].dE >= THRESH_KEV/1e6)
            {
               if (iok < i)
               {
                  ahits->in[iok] = ahits->in[i];
               }
               ++iok;
               ++mok;
            }
         }
         if (iok)
         {
            ahits->mult = iok;
         }
         else if (ahits != HDDM_NULL)
         {
            FREE(ahits);
         }
      }
      if ((wires != HDDM_NULL) && (mok == 0))
      {
         FREE(wires);
         wires = HDDM_NULL;
      }

      mok = 0;
      for (strip=0; strip < strips->mult; strip++)
      {
         s_FdcCathodeTruthHits_t* chits = strips->in[strip].fdcCathodeTruthHits;
         int i,iok;
         for (iok=i=0; i < chits->mult; i++)
         {
            if (chits->in[i].q >= THRESH_STRIPS)
            {
               if (iok < i)
               {
                  chits->in[iok] = chits->in[i];
               }
               ++iok;
               ++mok;
            }
         }
         if (iok)
         {
            chits->mult = iok;
         }
         else if (chits != HDDM_NULL)
         {
           FREE(chits);
         }
      }
      if ((strips != HDDM_NULL) && (mok == 0))
      {
         FREE(strips);
         strips = HDDM_NULL;
      }

      if ((wires != HDDM_NULL) || 
          (strips != HDDM_NULL) ||
          (points != HDDM_NULL))
      {
        if ((m == 0) || (module > box->fdcChambers->in[m-1].module)
                     || (layer  > box->fdcChambers->in[m-1].layer))
        {
          box->fdcChambers->in[m] = chambers->in[0];
          box->fdcChambers->in[m].fdcCathodeStrips = 
                              make_s_FdcCathodeStrips(stripCount);
          box->fdcChambers->in[m].fdcAnodeWires =
                              make_s_FdcAnodeWires(wireCount);
          box->fdcChambers->in[m].fdcTruthPoints =
                              make_s_FdcTruthPoints(pointCount);
          box->fdcChambers->mult++;
        }
        else
        {
          m--;
        }
        for (strip=0; strip < strips->mult; ++strip)
        {
           int mm = box->fdcChambers->in[m].fdcCathodeStrips->mult++;
           box->fdcChambers->in[m].fdcCathodeStrips->in[mm] = strips->in[strip];
        }
        if (strips != HDDM_NULL)
        {
           FREE(strips);
        }
        for (wire=0; wire < wires->mult; ++wire)
        {
           int mm = box->fdcChambers->in[m].fdcAnodeWires->mult++;
           box->fdcChambers->in[m].fdcAnodeWires->in[mm] = wires->in[wire];
        }
        if (wires != HDDM_NULL)
        {
           FREE(wires);
        }
        for (point=0; point < points->mult; ++point)
        {
           int mm = box->fdcChambers->in[m].fdcTruthPoints->mult++;
           box->fdcChambers->in[m].fdcTruthPoints->in[mm] = points->in[point];
        }
        if (points != HDDM_NULL)
        {
           FREE(points);
        }
     }
     FREE(chambers);
     FREE(item);
   }

   stripCount = wireCount = pointCount = 0;

   if ((box->fdcChambers != HDDM_NULL) &&
       (box->fdcChambers->mult == 0))
   {
      FREE(box->fdcChambers);
      box->fdcChambers = HDDM_NULL;
   }
   if (box->fdcChambers->mult == 0)
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
