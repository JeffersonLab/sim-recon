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
extern double asic_response(double t);
extern double Ei(double x);

typedef struct{
  int writeenohits;
  int showersincol;
  int driftclusters;
}controlparams_t;

extern controlparams_t controlparams_;


const float wire_dead_zone_radius[4]={3.0,3.0,3.9,3.9};
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
static float THRESH_ANODE         = 1.;
static float THRESH_STRIPS        =5. ;  /* pC */
static float ELECTRON_CHARGE =1.6022e-4; /* fC */
static float DIFFUSION_COEFF  =   1.1e-6; // cm^2/s --> 200 microns at 1 cm
static float FDC_TIME_WINDOW = 1000.0; //time window for accepting FDC hits, ns
static float GAS_GAIN = 8e4;


#if 0
static float wire_dx_offset[2304];
static float wire_dz_offset[2304];
#endif

binTree_t* forwardDCTree = 0;
static int stripCount = 0;
static int wireCount = 0;
static int pointCount = 0;
static int initializedx=0;

void gpoiss_(float*,int*,const int*); // avoid solaris compiler warnings
void rnorml_(float*,int*);

typedef int (*compfn)(const void*, const void*);

// Sort functions for sorting clusters
int fdc_anode_cluster_sort(const void *a,const void *b){
  const s_FdcAnodeTruthHit_t *ca=a;
  const s_FdcAnodeTruthHit_t *cb=b;
  if (ca->t<cb->t) return -1;
  else if (ca->t>cb->t) return 1;
  else return 0;
}
int fdc_cathode_cluster_sort(const void *a,const void *b){
  const s_FdcCathodeTruthHit_t *ca=a;
  const s_FdcCathodeTruthHit_t *cb=b;
  if (ca->t<cb->t) return -1;
  else if (ca->t>cb->t) return 1;
  else return 0;
}



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
      if ((den=ho-hp)==0.0) {
        free(c);
        free(d);
        return;
      }
      
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;
      
    }
    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
}

// Simulation of signal on a wire
double wire_signal(double t,s_FdcAnodeTruthHits_t* ahits){
  double t0=1.0; // ns; rough order of magnitude
  int m;
  double asic_gain=0.76; // mV/fC
  double func=0;
  for (m=0;m<ahits->mult;m++){
    if (t>ahits->in[m].t){
      double my_time=t-ahits->in[m].t;
      func+=asic_gain*ahits->in[m].dE*asic_response(my_time);
    }
  }
  return func;
}

// Simulation of signal on a cathode strip (ASIC output)
double cathode_signal(double t,s_FdcCathodeTruthHits_t* chits){
  double t0=1.0; // ns; rough order of magnitude
  int m;
  double asic_gain=2.3;
  double func=0;
  for (m=0;m<chits->mult;m++){
    if (t>chits->in[m].t){
      double my_time=t-chits->in[m].t;
      func+=asic_gain*chits->in[m].q*asic_response(my_time);
    }
  }
  return func;
}

// Generate hits in two cathode planes flanking the wire plane  
void AddFDCCathodeHits(int PackNo,float xwire,float avalanche_y,float tdrift,
		       int n_p,int track,int ipart,int chamber,int module,
		       int layer, int global_wire_number){

  s_FdcCathodeTruthHits_t* chits;    	  

  // Anode charge
  float q_anode;
  int n_t;
  // Average number of secondary ion pairs for 40/60 Ar/CO2 mixture
  float n_s_per_p=1.89; 
  if (controlparams_.driftclusters==0){    
    /* Total number of ion pairs.  On average for each primary ion 
       pair produced there are n_s secondary ion pairs produced.  The
       probability distribution is a compound poisson distribution
       that requires generating two Poisson variables.
    */
    int n_s,one=1;  
    float n_s_mean = ((float)n_p)*n_s_per_p;
    gpoiss_(&n_s_mean,&n_s,&one);
    n_t = n_s+n_p;
    q_anode=((float)n_t)*GAS_GAIN*ELECTRON_CHARGE;
  }
  else{
    // Distribute the number of secondary ionizations for this primary
    // ionization according to a Poisson distribution with mean n_s_over_p.
    // For simplicity we assume these secondary electrons and the primary
    // electron stay together as a cluster.
    int n_s;
    int one=1;
    gpoiss_(&n_s_per_p,&n_s,&one);
    // Anode charge in units of fC
    n_t=1+n_s;
    q_anode=GAS_GAIN*ELECTRON_CHARGE*((float)n_t);
  }

  /* Mock-up of cathode strip charge distribution */ 
  int plane, node;
  for (plane=1; plane<4; plane+=2){
    float theta = (plane == 1)? -CATHODE_ROT_ANGLE : +CATHODE_ROT_ANGLE;
    float cathode_u = xwire*cos(theta)+avalanche_y*sin(theta);
    int strip1 = ceil((cathode_u-U_OF_STRIP_ZERO)/STRIP_SPACING +0.5);
    float cathode_u1 = (strip1-1)*STRIP_SPACING + U_OF_STRIP_ZERO;
    float delta = cathode_u-cathode_u1;
    float half_gap=ANODE_CATHODE_SPACING;

#if 0
    half_gap+=(plane==1)?+wire_dz_offset[global_wire_number]:
	       -wire_dz_offset[global_wire_number];
#endif

    for (node=-STRIP_NODES; node<=STRIP_NODES; node++){
      /* Induce charge on the strips according to the Mathieson 
	 function tuned to results from FDC prototype
      */
      float lambda1=(((float)node-0.5)*STRIP_SPACING+STRIP_GAP/2.
		     -delta)/half_gap;
      float lambda2=(((float)node+0.5)*STRIP_SPACING-STRIP_GAP/2.
		     -delta)/half_gap;
      float factor=0.25*M_PI*K2;
      float q = 0.25*q_anode*(tanh(factor*lambda2)-tanh(factor*lambda1));
      
      int strip = strip1+node;
      /* Throw away hits on strips falling within a certain dead-zone
	 radius */
      float strip_outer_u=cathode_u1
	+(STRIP_SPACING+STRIP_GAP/2.)*(int)node;
      float cathode_v=-xwire*sin(theta)+avalanche_y*cos(theta);
      float check_radius=sqrt(strip_outer_u*strip_outer_u
			      +cathode_v*cathode_v);
      
      if ((strip > 0) 
	  && (check_radius>strip_dead_zone_radius[PackNo]) 
	  && (strip <= STRIPS_PER_PLANE)){
	int mark = (chamber<<20) + (plane<<10) + strip;
	void** cathodeTwig = getTwig(&forwardDCTree, mark);
	if (*cathodeTwig == 0){
	  s_ForwardDC_t* fdc = *cathodeTwig = make_s_ForwardDC();
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
	else{
	  s_ForwardDC_t* fdc = *cathodeTwig;
	  chits = fdc->fdcChambers->in[0].fdcCathodeStrips
	    ->in[0].fdcCathodeTruthHits;
	}
	
	int nhit;
	for (nhit = 0; nhit < chits->mult; nhit++){
	  // To cut down on the number of output clusters, combine 
	  // those that would be indistiguishable in time given the 
	  // expected timing resolution
	  if (fabs(chits->in[nhit].t - tdrift) <TWO_HIT_RESOL)
	    {
	      break;
	    }
	}
	if (nhit < chits->mult)		/* merge with former hit */
	  {
	    /* Use the time from the earlier hit but add the charge */
	    chits->in[nhit].q += q;
	    if(chits->in[nhit].t>tdrift){
	      chits->in[nhit].t = tdrift;
	      chits->in[nhit].itrack = track;
	      chits->in[nhit].ptype = ipart;
	    }
	  }
	else if (nhit < MAX_HITS){        /* create new hit */
	  chits->in[nhit].t = tdrift;
	  chits->in[nhit].q = q;
	  chits->in[nhit].itrack = track;
	  chits->in[nhit].ptype = ipart;
	  chits->mult++;
	}
	else{
	  // supress warning
	  /*
	    fprintf(stderr,"HDGeant error in hitForwardDC: ");
	    fprintf(stderr,"max hit count %d exceeded, truncating!\n",
	    MAX_HITS);
	  */
	}
	
      }
    } // loop over cathode strips
  } // loop over cathode views
}









// Add wire information
int AddFDCAnodeHit(s_FdcAnodeTruthHits_t* ahits,int layer,int ipart,int track,
		   float xwire,float xyz[3],float dE,float t,float *tdrift){
 
  // Generate 2 random numbers from a Gaussian distribution
  // 
  float rndno[2];
  int two=2;

  // Only and always use the built-in Geant random generator,
  // otherwise debugging is a problem because sequences are not
	  // reproducible from a given pair of random seeds. [rtj]
  
  /* rnorml_(rndno,&two); */ {
    float rho,phi1;
    grndm_(rndno,&two);
    rho = sqrt(-2*log(rndno[0]));
    phi1 = rndno[1]*2*M_PI;
    rndno[0] = rho*cos(phi1);
    rndno[1] = rho*sin(phi1);
  }

   // Get the magnetic field at this cluster position	    
  float x[3],B[3];
  transformCoord(xyz,"local",x,"global");
  gufld_db_(x,B);
  
  // Find the angle between the wire direction and the direction of the
  // magnetic field in the x-y plane
  float wire_dir[2];
  float wire_theta=1.0472*(float)((layer%3)-1);
  float phi=0.;;
  float Br=sqrt(B[0]*B[0]+B[1]*B[1]);
  
  wire_dir[0]=sin(wire_theta);
  wire_dir[1]=cos(wire_theta);
  if (Br>0.) phi= acos((B[0]*wire_dir[0]+B[1]*wire_dir[1])/Br);
  
  // useful combinations of dx and dz
  float dx=xyz[0]-xwire;
  float dx2=dx*dx;
  float dx4=dx2*dx2;
  float dz2=xyz[2]*xyz[2];
  float dz4=dz2*dz2;
  
  // Next compute the avalanche position along wire.  
  // Correct avalanche position with deflection along wire due to 
  // Lorentz force.
  xyz[1]+=( 0.1458*B[2]*(1.-0.048*Br) )*dx
    +( 0.1717+0.01227*B[2] )*(Br*cos(phi))*xyz[2]
    +( -0.000176 )*dx*dx2/(dz2+0.001);
  // Add transverse diffusion
  xyz[1]+=(( 0.01 )*pow(dx2+dz2,0.125)+( 0.0061 )*dx2)*rndno[1];

  // Do not use this cluster if the Lorentz force would deflect 
  // the electrons outside the active region of the detector
  if (sqrt(xyz[1]*xyz[1]+xwire*xwire)>ACTIVE_AREA_OUTER_RADIUS) 
    return 0;

  // Model the drift time and longitudinal diffusion as a function of 
  // position of the cluster within the cell 	 	  
  float tdrift_unsmeared=( 1086.0-106.7*B[2] )*dx2+( 1068.0 )*dz2
    +dx4*(( -2.675 )/(dz2+0.001)+( 2.4e4 )*dz2);	
  float dt=(( 39.44   )*dx4/(0.5-dz2)+( 56.0  )*dz4/(0.5-dx2)
      +( 0.01566 )*dx4/(dz4+0.002)/(0.251-dx2))*rndno[1];
  // Only allow deviations to longer times if near the cell border
  double dradius=sqrt(dx2+dz2);
  if (dradius-0.5>=0. && dt<0){ 
    dt=0.;
  }

  // Minimum drift time for docas near wire (very crude approximation)
  double v_max=0.08; // guess for now based on Garfield, near wire 
  double tmin=dradius/v_max;
  double tdrift_smeared=tdrift_unsmeared+dt;
  if (tdrift_smeared<tmin){
    tdrift_smeared=tmin;
  }

  // Avalanche time
  *tdrift=t+tdrift_smeared;
	  
  // Skip cluster if the time would go beyond readout window
  if (t>FDC_TIME_WINDOW) return 0;

  int nhit;

  // Record the anode hit
  for (nhit = 0; nhit < ahits->mult; nhit++)
    {
      if (fabs(ahits->in[nhit].t - *tdrift) < TWO_HIT_RESOL)
	{
	  break;
	}
    }
  if (nhit < ahits->mult)                 /* merge with former hit */
    {
      /* use the time from the earlier hit but add the energy */
      ahits->in[nhit].dE += dE;
      if(ahits->in[nhit].t>*tdrift){
	ahits->in[nhit].t = *tdrift;
	ahits->in[nhit].t_unsmeared=tdrift_unsmeared;
	ahits->in[nhit].d = sqrt(dx2+dz2);
	
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
      ahits->in[nhit].t = *tdrift;
      ahits->in[nhit].t_unsmeared=tdrift_unsmeared;
      ahits->in[nhit].dE = dE;
      ahits->in[nhit].d = sqrt(dx2+dz2);
      ahits->in[nhit].itrack = track;
      ahits->in[nhit].ptype = ipart;
      ahits->mult++;
    }
  else
    {
      fprintf(stderr,"HDGeant error in hitForwardDC: ");
      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
    }

  return 1;
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
  float alpha,sinalpha,cosalpha;
  int i,j;

  if (!initializedx){
      mystr_t strings[250];
      float values[250];
      int nvalues = 250;

      // Get parameters related to the geometry and the signals 
      int status = GetConstants("FDC/fdc_parms", &nvalues,values,strings);
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
#if 0
	{
	  int num_values=2304*2;
	  float my_values[2304*2]; 
	  mystr_t my_strings[2304*2];
	  status=GetArrayConstants("FDC/fdc_wire_offsets",&num_values,
				   my_values,my_strings);
	  if (!status){
	    int i;
	    for (i=0;i<num_values;i+=2){
	      int j=i/2;
	      wire_dx_offset[j]=my_values[i];
	      wire_dz_offset[j]=my_values[i+1];
	    }
	  }
	}
#endif
      }
      initializedx = 1 ;
  }

  /* Get chamber information */ 
  int layer = getlayer_wrapper_(); 
  if (layer==0){
    printf("hitFDC: FDC layer number evaluates to zero! THIS SHOULD NEVER HAPPEN! drop this particle.\n");
    return;
  }
  //int module = getmodule_wrapper_();
  //int chamber = (module*10)+layer;
  //int PackNo = (chamber-11)/20;
  int PackNo = getpackage_wrapper_()-1;
  if (PackNo==-1){
     printf("hitFDC: FDC package number evaluates to zero! THIS SHOULD NEVER HAPPEN! drop this particle.\n");
    return;
  }
  int glayer=6*PackNo+layer-1;
  int module = 2*(PackNo)+(layer-1)/3+1;
  int chamber = (module*10)+(layer-1)%3+1;
  int wire1,wire2;
  int wire,dwire;

  // transform layer number into Richard's scheme
  layer=(layer-1)%3+1;

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

  if (wire1>WIRES_PER_PLANE) wire1=wire2;
  else if (wire2>WIRES_PER_PLANE) wire2=wire1;
  if (wire1==0) wire1=wire2;
  else if (wire2==0) wire2=wire1;

  dwire = (wire1 < wire2)? 1 : -1;
  alpha = atan2(xoutlocal[0]-xinlocal[0],xoutlocal[2]-xinlocal[2]);
  sinalpha=sin(alpha);
  cosalpha=cos(alpha);
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
      dradius = fabs(u[1]*cosalpha-u[0]*sinalpha);
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
    float sign=1.; // for dealing with the y-position for tracks crossing two cells

    for (wire=wire1; wire-dwire != wire2; wire+=dwire)
    {
      int valid_hit=1;
      float dE,dt;
      float u[2];
      float x0[3],x1[3];
      float avalanche_y;
      float xwire = U_OF_WIRE_ZERO + (wire-1)*WIRE_SPACING;
      int global_wire_number=96*glayer+wire-1;

#if 0
      xwire+=wire_dx_offset[global_wire_number];
#endif     

      if (wire1==wire2){
	dE=dEsum; 
	x0[0] = xinlocal[0];
	x0[1] = xinlocal[1];
	x0[2] = xinlocal[2]; 
	x1[0] = xoutlocal[0];
	x1[1] = xoutlocal[1];
	x1[2] = xoutlocal[2];
      }
      else{
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
	dE = dEsum*(x1[2]-x0[2])/(xoutlocal[2]-xinlocal[2]);
      }

      if (dE > 0){
	s_FdcAnodeTruthHits_t* ahits;    

	// Create (or grab) an entry in the tree for the anode wire
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
	

	float rndno[2];
	int two=2;
      
	// Find the number of primary ion pairs:
	/* The total number of ion pairs depends on the energy deposition 
	   and the effective average energy to produce a pair, w_eff.
	   On average for each primary ion pair produced there are n_s_per_p 
	   secondary ion pairs produced. 
	*/
	int one=1;
	// Average number of secondary ion pairs for 40/60 Ar/CO2 mixture
	float n_s_per_p=1.89; 
	//Average energy needed to produce an ion pair for 50/50 mixture
	float w_eff=30.2e-9; // GeV
	// Average number of primary ion pairs
	float n_p_mean = dE/w_eff/(1.+n_s_per_p);
	int n_p; // number of primary ion pairs
	gpoiss_(&n_p_mean,&n_p,&one);
   
	// Drift time
	float tdrift=0;

	if (controlparams_.driftclusters==0){
	  float zrange=x1[2]-x0[2];
	  float tany=(x1[1]-x0[1])/zrange;
	  float tanx=(x1[0]-x0[0])/zrange;
	  float dz=ANODE_CATHODE_SPACING-dradius*sign*sinalpha;
#if 0	  
	  dz+=wire_dz_offset[global_wire_number];
#endif
	  xlocal[0]=x0[0]+tanx*dz;
	  if (fabs(xlocal[0]-xwire)>0.5){
	    xlocal[0]=x1[0];
	    xlocal[1]=x1[1];
	    xlocal[2]=x1[2];
	  }
	  else{
	    xlocal[1]=x0[1]+tany*dz;
	    xlocal[2]=x0[2]+dz;
	  }
	  
	  /* If the cluster position is within the wire-deadened region of the 
	     detector, skip this cluster 
	  */
	  if (sqrt(xlocal[0]*xlocal[0]+xlocal[1]*xlocal[1])
	      >=wire_dead_zone_radius[PackNo]){	 
	    if (AddFDCAnodeHit(ahits,layer,ipart,track,xwire,xlocal,dE,t,
			       &tdrift)){
	      AddFDCCathodeHits(PackNo,xwire,xlocal[1],tdrift,n_p,track,ipart,
				chamber,module,layer,global_wire_number);
	    }
	  
	  }
	}
	else{
	  // Loop over the number of primary ion pairs
	  int n;
	  for (n=0;n<n_p;n++){
	    // Generate a cluster at a random position along the path with cell
	    float rndno[2];
	    grndm_(rndno,&two);
	    float dzrand=(x1[2]-x0[2])*rndno[0];
	    // Position of the cluster
	    xlocal[0]=x0[0]+(x1[0]-x0[0])*rndno[0];
	    xlocal[1]=x0[1]+(x1[1]-x0[1])*rndno[0];
	    xlocal[2]=x0[2]+dzrand;
	    /* If the cluster position is within the wire-deadened region of the 
	       detector, skip this cluster 
	    */
	    if (sqrt(xlocal[0]*xlocal[0]+xlocal[1]*xlocal[1])
		>=wire_dead_zone_radius[PackNo]){	   
	      if (AddFDCAnodeHit(ahits,layer,ipart,track,xwire,xlocal,dE,t,
				 &tdrift)){
		AddFDCCathodeHits(PackNo,xwire,xlocal[1],tdrift,n_p,track,ipart,
				  chamber,module,layer,global_wire_number);
	      }
	    }

	  } // loop over primary ion pairs 
	}
      } // Check for non-zero energy

      sign*=-1; // for dealing with the y-position for tracks crossing two cells
    } // loop over wires
  } // Check that total energy deposition is not zero
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
	   
	 // Sort the clusters by time
	 qsort(ahits->in,ahits->mult,sizeof(s_FdcAnodeTruthHit_t),
	       (compfn)fdc_anode_cluster_sort);
	 
	 int i,iok=0;

	 if (controlparams_.driftclusters==0){
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
	 }
	 else{  // Simulate clusters within the cell
	
	   // printf("-------------\n");
	  
	   
	   // Temporary histogram in 1 ns bins to store waveform data
	   int num_samples=(int)FDC_TIME_WINDOW;
	   float *samples=(float *)malloc(num_samples*sizeof(float));
	   for (i=0;i<num_samples;i++){
	     samples[i]=wire_signal((float)i,ahits);
	     //printf("%f %f\n",(float)i,samples[i]);
	   }
	   
	   int returned_to_baseline=0;
	   float q=0;
	   for (i=0;i<num_samples;i++){
	     if (samples[i]>=THRESH_ANODE){
	       if (returned_to_baseline==0){
		 ahits->in[iok].itrack = ahits->in[0].itrack;
		 ahits->in[iok].ptype = ahits->in[0].ptype;
		 
		 // Do an interpolation to find the time at which the threshold
		 // was crossed.
		 float t_array[4];
		 int k;
		 float my_t,my_terr;
		 for (k=0;k<4;k++) t_array[k]=i-1+k;
		 polint(&samples[i-1],t_array,4,THRESH_ANODE,&my_t,&my_terr);
		 ahits->in[iok].t=my_t;
		 
		 returned_to_baseline=1;
		 iok++;
		 mok++;
	       }
	       q+=samples[i];
	     }
	     if (returned_to_baseline 
		 && (samples[i]<THRESH_ANODE)){
	       returned_to_baseline=0;   
	     }
	   }
	   free(samples);
	 } // Simulation of clusters within cell

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
	  
	  // Sort the clusters by time
	  qsort(chits->in,chits->mult,sizeof(s_FdcCathodeTruthHit_t),
		  (compfn)fdc_cathode_cluster_sort);
	  
	  int i,iok=0;
	  
	  if (controlparams_.driftclusters==0){
	    for (iok=i=0; i < chits->mult; i++)
	      {
		if (chits->in[i].q > 0.)
		  {
		    if (iok < i)
		      {
			chits->in[iok] = chits->in[i];
		      }
		    ++iok;
		    ++mok;
		  }
	      }	    
	    
	  }
	  else{
	   
	    // Temporary histogram in 1 ns bins to store waveform data
	    int num_samples=(int)(FDC_TIME_WINDOW);
	    float *samples=(float *)malloc(num_samples*sizeof(float));
	    for (i=0;i<num_samples;i++){
	      samples[i]=cathode_signal((float)i,chits);
		 //printf("t %f V %f\n",(float)i,samples[i]);
	    }	 
	    
	    int threshold_toggle=0;
	    int istart=0;
	    float q=0;
	    int FADC_BIN_SIZE=1;
	    for (i=0;i<num_samples;i+=FADC_BIN_SIZE){
	      if (samples[i]>=THRESH_STRIPS){
		if (threshold_toggle==0){
		  chits->in[iok].itrack = chits->in[0].itrack;
		  chits->in[iok].ptype = chits->in[0].ptype;
		  chits->in[iok].t=(float) i;
		  //chits->in[iok].q=samples[i];
		  istart=i-1;
		  threshold_toggle=1;
		  //iok++;
		  //mok++;
		}
	      }
	      if (threshold_toggle && 
		  (samples[i]<THRESH_STRIPS)){
		int j;
		// Find the first peak
		for (j=istart+1;j<i-1;j++){
		  if (samples[j]>samples[j-1] && samples[j]>samples[j+1]){
		    chits->in[iok].q=samples[j];
		    break;
		  }
		} 
		threshold_toggle=0; 
		iok++;
		mok++;
		//break;
	      }
	    }
	    i=num_samples-1;
	    if (samples[i]>=THRESH_STRIPS&&threshold_toggle){
	      int j;
	      for (j=istart+1;j<i-1;j++){
		if (samples[j]>samples[j-1] && samples[j]>samples[j+1]){
		  chits->in[iok].q=samples[j];
		  break;
		}
	      }
	    }
	    
	    free(samples);
	  }// Simulate clusters within cell
	
	  if (iok)
	    {
	      chits->mult = iok;
	      //chits->mult=1;
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
