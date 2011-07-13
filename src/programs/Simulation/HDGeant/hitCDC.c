/*
 * hitCDC - registers hits for Central Drift Chamber
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function call hitCentralDC
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#include "calibDB.h"
extern s_HDDM_t* thisInputEvent;

void gpoiss_(float*,int*,const int*); // avoid solaris compiler warnings

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
static float DRIFT_SPEED     =   0.0055;
static float TWO_HIT_RESOL   = 1.;
static int   MAX_HITS 	     = 1000;
static float THRESH_KEV	     =   1.;
static float THRESH_MV = 1.;
static float STRAW_RADIUS    =   0.8;
static float CDC_TIME_WINDOW = 1000.0; //time window for accepting CDC hits, ns
static float ELECTRON_CHARGE =1.6022e-4; /* fC */
static float GAS_GAIN = 1e5;

binTree_t* centralDCTree = 0;
static int strawCount = 0;
static int pointCount = 0;
static int stripCount = 0;
static int initialized = 0;

/* void GetDOCA(int ipart, float x[3], float p[5], float doca[3]);  disabled 6/24/2009 */

typedef int (*compfn)(const void*, const void*);

// Sort function for sorting clusters
int cdc_cluster_sort(const void *a,const void *b){
  const s_CdcStrawTruthHit_t *ca=a;
  const s_CdcStrawTruthHit_t *cb=b;
  if (ca->t<cb->t) return -1;
  else if (ca->t>cb->t) return 1;
  else return 0;
}

// Simulation of the ASIC response to a pulse due to a cluster
double asic_response(double t){
  double func=0;
  double par[11]={-0.01986,0.01802,-0.001097,10.3,11.72,-0.03701,35.84,
		  15.93,0.006141,80.95,24.77};
  if (t<par[3]){
    func=par[0]*t+par[1]*t*t+par[2]*t*t*t;
  }
  else{
    func+=(par[0]*par[3]+par[1]*par[3]*par[3]+par[2]*par[3]*par[3]*par[3])
      *exp(-(t-par[3])*(t-par[3])/(par[4]*par[4]));
    func+=par[5]*exp(-(t-par[6])*(t-par[6])/(par[7]*par[7]));
    func+=par[8]*exp(-(t-par[9])*(t-par[9])/(par[10]*par[10]));
  }
  return func;
}

// Simulation of signal on a wire
double cdc_wire_signal(double t,s_CdcStrawTruthHits_t* chits){
  double t0=1.0; // ns; rough order of magnitude
  int m;
  double asic_gain=0.5; // mV/fC
  double func=0;
  for (m=0;m<chits->mult;m++){
    if (t>chits->in[m].t){
      double my_time=t-chits->in[m].t;
      func+=asic_gain*chits->in[m].dE*asic_response(my_time);
    }
  }
  return func;
}




/* register hits during tracking (from gustep) */

void hitCentralDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum,
                   int track, int stack, int history, int ipart )
{
  float x[3], t;
  float dx[3], dr;
  float dEdx;
  float xlocal[3];
  float xinlocal[3];
  float xoutlocal[3];
  float dradius,drin,drout;
  float doca[3]; 
  float trackdir[3];
  float alpha;
  
  if (!initialized) {
    mystr_t strings[50];
    float values[50];
    int nvalues = 50;
    int status = GetConstants("CDC/cdc_parms", &nvalues, values, strings);
    
    if (!status) {
      int ncounter = 0;
      int i;
      for ( i=0;i<(int)nvalues;i++){
	//printf("%d %s \n",i,strings[i].str);
	if (!strcmp(strings[i].str,"CDC_DRIFT_SPEED")) {
	  DRIFT_SPEED  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_TWO_HIT_RESOL")) {
	  TWO_HIT_RESOL  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_MAX_HITS")) {
	  MAX_HITS  = (int)values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_THRESH_KEV")) {
	  THRESH_KEV  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_STRAW_RADIUS")) {
	  THRESH_KEV  = values[i];
	  ncounter++;
	}
      }
      if (ncounter==5){
	printf("CDC: ALL parameters loaded from Data Base\n");
      } else if (ncounter<5){
	printf("CDC: NOT ALL necessary parameters found in Data Base %d out of 5\n",ncounter);
      } else {
	printf("CDC: SOME parameters found more than once in Data Base\n");
      } 	
    }
    initialized = 1;
  }
    
  x[0] = (xin[0] + xout[0])/2;
  x[1] = (xin[1] + xout[1])/2;
  x[2] = (xin[2] + xout[2])/2;
  t    = (xin[3] + xout[3])/2 * 1e9;
  dx[0] = xin[0] - xout[0];
  dx[1] = xin[1] - xout[1];
  dx[2] = xin[2] - xout[2];
  transformCoord(xin,"global",xinlocal,"local");
  transformCoord(xout,"global",xoutlocal,"local");
  /*
  xlocal[0] = (xinlocal[0] + xoutlocal[0])/2;
  xlocal[1] = (xinlocal[1] + xoutlocal[1])/2;
  xlocal[2] = (xinlocal[2] + xoutlocal[2])/2;
  */
  
  /* For particles that range out inside the active volume, the
   * "out" time seems to be set to something enormously high.
   * This screws up the hit. Check for this case here by looking
   * at xout[3] and making sure it is less than 1 second. If it's
   * not, then just use xin[3] for "t".
   */
  if(xout[3] > 1.0) t = xin[3] * 1e9;
	 
  drin = sqrt(xinlocal[0]*xinlocal[0] + xinlocal[1]*xinlocal[1]);
  drout = sqrt(xoutlocal[0]*xoutlocal[0] + xoutlocal[1]*xoutlocal[1]);
  
  trackdir[0] =-xinlocal[0] + xoutlocal[0];
  trackdir[1] =-xinlocal[1] + xoutlocal[1];
  trackdir[2] =-xinlocal[2] + xoutlocal[2];
  alpha=-(xinlocal[0]*trackdir[0]+xinlocal[1]*trackdir[1])
    /(trackdir[0]*trackdir[0]+trackdir[1]*trackdir[1]);
  xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
  xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
  xlocal[2]=xinlocal[2]+trackdir[2]*alpha;
  
  // Deal with tracks exiting the ends of the straws
  if (fabs(xlocal[2])>=75.0){
    int ring = getring_();
    float sign=(xoutlocal[2]>0)?1.:-1.;
    if (ring<=4 || (ring>=13 && ring<=16) || ring>=25){
      alpha=(sign*75.0-xinlocal[2])/trackdir[2];
      xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
      xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
      xlocal[2]=sign*75.0;
    }
    else if (fabs(xlocal[2])>=75.415){
      alpha=(sign*75.415-xinlocal[2])/trackdir[2]; 
      xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
      xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
      xlocal[2]=sign*75.415;
    }
  } 

  /* This will get called when the particle actually passes through
   * the wire volume itself. For these cases, we should set the 
   * location of the hit to be the point on the wire itself. Do
   * determine if this is what is happening, we check drout to
   * see if it is very close to the wire and drin to see if it is
   * close to the tube.
   *
   * For the other case, when drin is close to the wire, we assume
   * it is because it is emerging from the wire volume and
   * automatically ignore those hits by returning immediately.
   */
  if(drin < 0.0050)return; /* entering straw within 50 microns of wire. ignore */
  if(drin>(STRAW_RADIUS-0.0200) && drout<0.0050){
    /* We entered within 200 microns of the straw tube and left
     * within 50 microns of the wire. Assume the track passed through
     * the wire volume.
     */
    
    x[0] = xin[0];
    x[1] = xin[1];
    x[2] = xin[2];
    t = xin[3] * 1e9;
    xlocal[0] = xinlocal[0];
    xlocal[1] = xinlocal[1];
    xlocal[2] = xinlocal[2];
    
    /* For dx, we will just assume it is twice the distance from
     * the straw to wire.
     */
    dx[0] *= 2.0;
    dx[1] *= 2.0;
    dx[2] *= 2.0;
  }
  
  /* Distance of hit from center of wire */
  dradius = sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);

  /* Calculate dE/dx */

   dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
   if (dr > 1e-3)
   {
      dEdx = dEsum/dr;
   }
   else
   {
      dEdx = 0;
   }

   /* post the hit to the truth tree */

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&centralDCTree, mark);
      if (*twig == 0)
      {
         s_CentralDC_t* cdc = *twig = make_s_CentralDC();
         s_CdcTruthPoints_t* points = make_s_CdcTruthPoints(1);
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (stack <= a);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].dradius = dradius;
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
         points->mult = 1;
         cdc->cdcTruthPoints = points;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   if (dEsum > 0)
   {
      int nhit;
      s_CdcStrawTruthHits_t* hits;
#if CATHODE_STRIPS_IN_CDC
      s_CdcCathodeStrips_t* chits;
      int cell = getcell_();
#endif
      int layer = getlayer_();
      int ring = getring_();
      int sector = getsector_();

      if (layer == 0)		/* in a straw */
      {	
	int mark = (ring<<20) + sector;
	void** twig = getTwig(&centralDCTree, mark);
	if (*twig == 0)
	{
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcStraws_t* straws = make_s_CdcStraws(1);
            straws->mult = 1;
            straws->in[0].ring = ring;
            straws->in[0].straw = sector;
            straws->in[0].cdcStrawTruthHits = hits = make_s_CdcStrawTruthHits(MAX_HITS);
            cdc->cdcStraws = straws;
            strawCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->cdcStraws->in[0].cdcStrawTruthHits;
         }

	/* Simulate number of primary ion pairs*/
	/* The total number of ion pairs depends on the energy deposition 
	   and the effective average energy to produce a pair, w_eff.
	   On average for each primary ion pair produced there are n_s_per_p 
	   secondary ion pairs produced. 
	*/
	int one=1;
	// Average number of secondary ion pairs for 50/50 Ar/CO2 mixture
	float n_s_per_p=1.94; 
	//Average energy needed to produce an ion pair for 50/50 mixture
	float w_eff=29.5e-9; // GeV
	// Average number of primary ion pairs
	float n_p_mean = dEsum/w_eff/(1.+n_s_per_p);
	int n_p; // number of primary ion pairs
	gpoiss_(&n_p_mean,&n_p,&one);
	
	// Loop over the number of primary ion pairs
	int n;
	for (n=0;n<n_p;n++){
	  // Generate a cluster at a random position along the path within the 
	  // straw
	  int two=2;
	  float rndno[2];
	  grndm_(rndno,&two);
	  xlocal[0]=xinlocal[0]+trackdir[0]*rndno[0];
	  xlocal[1]=xinlocal[1]+trackdir[1]*rndno[0];
	  xlocal[2]=xinlocal[2]+trackdir[2]*rndno[0];
	  dradius=sqrt(xlocal[0]*xlocal[0]+xlocal[1]*xlocal[1]);

	  // Find the drift time for this cluster. Drift time depends on B:
	  // (dependence derived from Garfield calculations)
	  float B[3],Bmag; 
	  transformCoord(xlocal,"local",x,"global");
	  gufld2_(x,B);
	  Bmag=sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	  float d2=dradius*dradius;
	  float d3=d2*dradius;
	  float tdrift=t+(-49.41+4.74*Bmag)*dradius+(1129.0+78.66*Bmag)*d2;

	  //Longitudinal diffusion
	  float rho = sqrt(-2*log(rndno[0]));
	  float phi = rndno[1]*2*M_PI;
	  float dt=(7.515*dradius-2.139*d2+12.63*d3)*rho*cos(phi);
	  tdrift+=dt;

	  // Skip cluster if the time would go beyond readout window
	  if (tdrift>CDC_TIME_WINDOW) continue;
	  
	  // Distribute the number of secondary ionizations for this primary
	  // ionization according to a Poisson distribution with mean n_s_over_p.
	  // For simplicity we assume these secondary electrons and the primary
	  // electron stay together as a cluster.
	  int n_s;
	  gpoiss_(&n_s_per_p,&n_s,&one);
	  // Energy deposition, equivalent to anode charge, in units of fC
	  float dE=GAS_GAIN*ELECTRON_CHARGE*(float)(1+n_s);

	  for (nhit = 0; nhit < hits->mult; nhit++)
	    {
	      if (fabs(hits->in[nhit].t - tdrift) < TWO_HIT_RESOL)
		{
		  break;
		}
	    }
	  if (nhit < hits->mult)		/* merge with former hit */
	    {
	      /* Use the time from the earlier hit but add the charge*/
	      hits->in[nhit].dE += dE;
	      if(hits->in[nhit].t>tdrift){
		hits->in[nhit].t = tdrift;
		hits->in[nhit].d = dradius;
		hits->in[nhit].itrack = track;
		hits->in[nhit].ptype = ipart;
	      }
	      
	      /* hits->in[nhit].t =
		 (hits->in[nhit].t * hits->in[nhit].dE + tdrift * dEsum)
		 / (hits->in[nhit].dE += dEsum);
	      */
	    }
	  else if (nhit < MAX_HITS)		/* create new hit */
	    {
	      hits->in[nhit].t = tdrift;
	      hits->in[nhit].dE = dE;
	      hits->in[nhit].d = dradius;
	      hits->in[nhit].itrack = track;
	      hits->in[nhit].ptype = ipart;
	      
	      hits->mult++;
	    }
	  else
	    {
	      fprintf(stderr,"HDGeant error in hitCentralDC: ");
	      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
	    }
	}
      }
#if CATHODE_STRIPS_IN_CDC
      else			/* in one of the z-readout (strip) layers */
      {
         int nchit;
         int mark = (layer << 28) + (sector << 20) + (cell << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcCathodeStrips strips = make_s_CdcCathodeStrips(1);
            strips->mult = 1;
            strips->in[0].layer = layer;
            strips->in[0].sector = sector;
            strips->in[0].strip = cell;
            strips->in[0].cdcCathodeHits = chits
                                       = make_s_CdcCathodeHits(MAX_HITS);
            cdc->cdcCathodeStrips = strips;
            stripCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            chits = cdc->cdcCathodeStrips->in[0].cdcCathodeHits;
         }

         for (nchit = 0; nchit < chits->mult; nchit++)
         {
            if (fabs(chits->in[nchit].t - tdrift) < TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nchit < chits->mult)		/* merge with former hit */
         {
            chits->in[nchit].t =
                    (chits->in[nchit].t * chits->in[nchit].dE + t * dEsum)
                  / (chits->in[nchit].dE += dEsum);
         }
         else if (nchit < MAX_HITS)		/* create new hit */
         {
            chits->in[nchit].t = t;
            chits->in[nchit].dE = dEsum;
            chits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitCentralDC: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }
      }
#endif
   }
}

/* entry points from fortran */

void hitcentraldc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitCentralDC(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


/* pick and package the hits for shipping */

s_CentralDC_t* pickCentralDC ()
{
   s_CentralDC_t* box;
   s_CentralDC_t* item;

   if ((strawCount == 0) && (stripCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_CentralDC();
   box->cdcStraws = make_s_CdcStraws(strawCount);
   box->cdcTruthPoints = make_s_CdcTruthPoints(pointCount);
#if CATHODE_STRIPS_IN_CDC
   box->cdcCathodeStrips = make_s_CdcCathodeStrips(stripCount);
#endif
   while (item = (s_CentralDC_t*) pickTwig(&centralDCTree))
   {
      s_CdcStraws_t* straws = item->cdcStraws;
      int straw;
#if CATHODE_STRIPS_IN_CDC
      s_CdcCathodetrips_t* strips = item->cdcCathodeStrips;
      int strip;
#endif
      s_CdcTruthPoints_t* points = item->cdcTruthPoints;
      int point;
      for (straw=0; straw < straws->mult; ++straw)
      {
         int m = box->cdcStraws->mult;

         s_CdcStrawTruthHits_t* hits = straws->in[straw].cdcStrawTruthHits;

	 // Sort the clusters by time
	 qsort(hits->in,hits->mult,sizeof(s_CdcStrawTruthHit_t),(compfn)cdc_cluster_sort);
	
         /* compress out the hits below threshold */
         int i,iok=0;

	 // Temporary histogram in 1 ns bins to store waveform data
	 int num_samples=(int)CDC_TIME_WINDOW;
	 float *samples=(float *)malloc(num_samples*sizeof(float));
	 for (i=0;i<num_samples;i++){
	   samples[i]=cdc_wire_signal((float)i,hits);
	   //printf("%f %f\n",(float)i,samples[i]);
	 }
	 
	 int returned_to_baseline=0;
	 float q=0.; 
	 float w_eff=0.0295; //keV	
	 int FADC_BIN_SIZE=1;
	 for (i=0;i<num_samples;i+=FADC_BIN_SIZE){
	   if (samples[i]>=THRESH_MV){
	     if (returned_to_baseline==0){
	       hits->in[iok] = hits->in[0];
	       hits->in[iok].t=(float) i;
	       returned_to_baseline=1;
	       iok++;
	     }
	     q+=(float)FADC_BIN_SIZE*samples[i];
	   }
	   if (samples[i]<THRESH_MV || i==num_samples-1){
	     returned_to_baseline=0;   
	     if (iok>0 && q>0.){
	       hits->in[iok-1].dE=q*w_eff/(GAS_GAIN*ELECTRON_CHARGE);;
	       q=0.;
	     }
	   }
	 }
	 if (q>0){
	   hits->in[iok-1].dE=q*w_eff/(GAS_GAIN*ELECTRON_CHARGE);
	 }
	 free(samples);

         if (iok)
         {
            hits->mult = iok;
            box->cdcStraws->in[m] = straws->in[straw];
            box->cdcStraws->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (straws != HDDM_NULL)
      {
         FREE(straws);
      }
#if CATHODE_STRIPS_IN_CDC
      for (strip=0; strip < strips->mult; ++strip)
      {
         int m = box->cdcCathodeStrips->mult;

         s_CdcCathodeHits* hits = strips->in[strip].cdcCathodeHits;

         /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < hits->mult; i++)
         {
            if (hits->in[i].dE >= THRESH_KEV/1e6)
            {
               if (iok < i)
               {
                  hits->in[iok] = hits->in[i];
               }
               ++iok;
            }
         }

         if (iok)
         {
            hits->mult = iok;
            box->cdcCathodeStrips->in[m] = strips->in[strip];
            box->cdcCathodeStrips->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (strips != HDDM_NULL)
      {
         FREE(strips);
      }
#endif
      for (point=0; point < points->mult; ++point)
      {
         int m = box->cdcTruthPoints->mult++;
         box->cdcTruthPoints->in[m] = points->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   strawCount = stripCount = pointCount = 0;

#if CATHODE_STRIPS_IN_CDC
   if ((box->cdcCathodeStrips != HDDM_NULL) &&
       (box->cdcCathodeStrips->mult == 0))
   {
      FREE(box->cdcCathodeStrips);
      box->cdcCathodeStrips = HDDM_NULL;
   }
#endif
   if ((box->cdcStraws != HDDM_NULL) &&
       (box->cdcStraws->mult == 0))
   {
      FREE(box->cdcStraws);
      box->cdcStraws = HDDM_NULL;
   }
   if ((box->cdcTruthPoints != HDDM_NULL) &&
       (box->cdcTruthPoints->mult == 0))
   {
      FREE(box->cdcTruthPoints);
      box->cdcTruthPoints = HDDM_NULL;
   }
   if ((box->cdcStraws->mult == 0) &&
       (box->cdcTruthPoints->mult == 0))
   {
#if CATHODE_STRIPS_IN_CDC
    if (box->cdcCathodeStrips->mult == 0) {
#endif
      FREE(box);
      box = HDDM_NULL;
#if CATHODE_STRIPS_IN_CDC
    }
#endif
   }
   return box;
}

