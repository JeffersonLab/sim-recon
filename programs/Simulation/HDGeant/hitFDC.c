/*
 * hitFDC - registers hits for forward drift chambers
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define DRIFT_SPEED	.0022
#define STRIP_SPACING	0.5
#define STRIP_GAP       0.1
#define ANODE_CATHODE_SPACING 0.5
#define TWO_HIT_RESOL	250.
#define WIRE_SPACING 	1.0
#define MAX_HITS 	100
#define K2              1.15
#define NODES           3
#define ELECTRON_CHARGE 1.6022e-4 /* fC */

binTree_t* forwardDCTree = 0;
static int stripCount = 0;
static int wireCount = 0;
static int pointCount = 0;

#define U0		100.

void rnpssn_(float*,int*,int*); // avoid solaris compiler warnings

/* register hits during tracking (from gustep) */

void hitForwardDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum,
                   int track, int stack)
{
  float x[3], t;
  float dx[3], dr;
  float dEdx;
  float xlocal[3];
  float xfcal[3];
  float zeroHat[] = {0,0,0};
  int i,j;
  

  if (dEsum == 0) return;              /* only seen if it deposits energy */

  
  {
    /* Get plane information */
    int module = getmodule_();
    int layer = getlayer_();
    int plane = getplane_();
    int player = (layer -1)*3 + plane;
    
    transformCoord(zeroHat,"local",xfcal,"global");
    
    x[0] = (xin[0] + xout[0])/2;
    x[1] = (xin[1] + xout[1])/2;
    x[2] = (xin[2] + xout[2])/2;
    t    = (xin[3] + xout[3])/2 * 1e9;
    transformCoord(x,"global",xlocal,"local");
    
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
        
    /* post the hit to the hits tree, mark cell as hit */
    {
      int nhit;
      s_Hits_t* hits;    
      const float tau[] = {0,-45,0,45,15,60,105,-105,-60,-15};
      int iwire = (xlocal[0] + U0)/WIRE_SPACING +1;
      float dradius = fabs((xlocal[0] + U0) - (iwire - 0.5)*WIRE_SPACING);
      float tdrift = t + dradius/DRIFT_SPEED;

      if (plane == 2){		/* anode drift cell plane */
	int mark = (module << 24) + (player << 20) + (iwire << 8);
	void** twig = getTwig(&forwardDCTree, mark);

	if (*twig == 0){
	  s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
	  fdc->chambers = make_s_Chambers(1);
	  fdc->chambers->mult = 1;

	  /* anode information */
	  fdc->chambers->in[0].module = module;
	  fdc->chambers->in[0].layer = layer;
	  fdc->chambers->in[0].anodePlanes = make_s_AnodePlanes(1);
	  fdc->chambers->in[0].anodePlanes->mult = 1;
	  fdc->chambers->in[0].anodePlanes->in[0].z = xfcal[2];
	  fdc->chambers->in[0].anodePlanes->in[0].tau = tau[player];
	  fdc->chambers->in[0].anodePlanes->in[0].wires = make_s_Wires(1);
	  fdc->chambers->in[0].anodePlanes->in[0].wires->mult = 1;
	  fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].u =
	    (floor(xlocal[0]/WIRE_SPACING) + 0.5)*WIRE_SPACING;
	  fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
	  wireCount++;	  
	}
	else
	{
	   s_ForwardDC_t* fdc = *twig;
	   hits = fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].hits;
	}
          
	for (nhit = 0; nhit < hits->mult; nhit++)
	{
	    if (fabs(hits->in[nhit].t - t) < TWO_HIT_RESOL)
	      {
		  break;
	      }
	  }
	if (nhit < hits->mult)                 /* merge with former hit */
	  {
	    hits->in[nhit].t = 
	      (hits->in[nhit].t * hits->in[nhit].dE + tdrift * dEsum)
	      / (hits->in[nhit].dE += dEsum);
	  }
	else if (nhit < MAX_HITS)              /* create new hit */
	  {
	    hits->in[nhit].t = tdrift;
	    hits->in[nhit].dE = dEsum;
	    hits->mult++;
	  }
	else
	  {
	    fprintf(stderr,"HDGeant error in hitForwardDC: ");
	    fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
	  }

 /* post the hit to the truth tree, once per primary track per cell */
 
         mark = (module << 24) + (player << 20) + (iwire << 8) + track;
         twig = getTwig(&forwardDCTree, mark);
         if (*twig == 0)
         {
            s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
            s_FdcPoints_t* points = make_s_FdcPoints(1);
            fdc->chambers = make_s_Chambers(1);
            fdc->chambers->mult = 1;
            fdc->chambers->in[0].module = module;
            fdc->chambers->in[0].layer = layer;
            fdc->chambers->in[0].anodePlanes = make_s_AnodePlanes(1);
            fdc->chambers->in[0].anodePlanes->mult = 1;
            fdc->chambers->in[0].anodePlanes->in[0].z = xfcal[2];
            fdc->chambers->in[0].anodePlanes->in[0].tau = tau[player];
            fdc->chambers->in[0].anodePlanes->in[0].wires = make_s_Wires(1);
            fdc->chambers->in[0].anodePlanes->in[0].wires->mult = 1;
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].u =
                            (floor(xlocal[0]/WIRE_SPACING) + 0.5)*WIRE_SPACING;
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].fdcPoints =
                                                                    points;
            points->in[0].primary = (stack == 0);
            points->in[0].track = track;
            points->in[0].x = x[0];
            points->in[0].y = x[1];
            points->in[0].z = x[2];
            points->in[0].dradius = dradius;
            points->in[0].dEdx = dEdx;
            points->mult++;
            pointCount++;

         }
      
	/* cathode information */ 
	for (j=-1;j<2;j+=2){
	  float theta=M_PI/4*(float)j;
	  int player=(layer -1)*3 + plane+j;
	  int mark = (module << 24) + (player << 20) + (iwire << 8);
	  
	  void** twig = getTwig(&forwardDCTree, mark);
	  if (*twig == 0){
	    s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
	    fdc->chambers = make_s_Chambers(1);
            fdc->chambers->mult = 1;
            fdc->chambers->in[0].module = module;
            fdc->chambers->in[0].layer = layer;
	    fdc->chambers->in[0].cathodePlanes = make_s_CathodePlanes(1);
	    fdc->chambers->in[0].cathodePlanes->mult = 1;
	    fdc->chambers->in[0].cathodePlanes->in[0].z = 
	      xfcal[2]+ANODE_CATHODE_SPACING*(float)j;
	    fdc->chambers->in[0].cathodePlanes->in[0].tau = tau[player];
	    fdc->chambers->in[0].cathodePlanes->in[0].strips 
	      = make_s_Strips(2*NODES+1);
	    fdc->chambers->in[0].cathodePlanes->in[0].strips->mult 
	      = 2*NODES+1;

	    /* Mock-up of cathode strip charge distribution */
	    {
	      float cathode_local=xlocal[1]*sin(theta)+
		(floor(xlocal[0]/WIRE_SPACING) + 0.5)*WIRE_SPACING*cos(theta);
	      float u=(floor(cathode_local/STRIP_SPACING)+0.5)*STRIP_SPACING;
	      float delta=cathode_local-u;

	      /* Variables for approximating number of ion pairs  */ 
	      float w_eff=26.0; /*eV, average energy needed to produce an 
				  ion pair.  For simplicity use pure argon */
	      float n_s_per_p=2.2; /* average number of secondary ion pairs per
				primary ionization */
	      float n_p_mean,n_s_mean;
	      int n_p,n_t,n_s,err;
      
	      /* variables for gain approximation */
	      float gain,q_anode;
	      float k=1.81;   /*empirical constant related to first Townsend
				coefficient of argon gas */
	      float N=269.*273/293;  /* proportional to number of gas 
					molecules/cm^3 */
	      float V=1800; /* V, operation voltage */
	      float V_t=703.5; /* V, threshold voltage */
	      float amp_gain=46; /* 2.3 mV/uA / 50 Ohms */
	      float C=0.946; /* capacitance in units of episilon_0 */
	      float a=0.001; /* cm, radius of sense wires */
	      /* Sauli eq. 31: */
	      gain=exp(sqrt(2.*k*V*C*N*a/M_PI)*(sqrt(V/V_t)-1.));

	      /* Total number of ion pairs.  On average for each primary ion 
		 pair produced there are n_s secondary ion pairs produced.  The
	         probability distribution is a compound poisson distribution
	         that requires calling rnpssn_ twice.*/
	      n_p_mean=dEsum/w_eff/(1.+n_s_per_p)*1e9;
	      rnpssn_(&n_p_mean,&n_p,&err);
	      n_s_mean=((float)n_p)*n_s_per_p;
	      rnpssn_(&n_s_mean,&n_s,&err);
	      n_t=n_s+n_p;
	      q_anode=((float)n_t)*gain*ELECTRON_CHARGE;

	      for (i=-NODES;i<=NODES;i++){
		float dE=0.;		

		/* Parameters needed for charge distribution */
		float lambda1=(((float)i-0.5)*STRIP_SPACING+STRIP_GAP/2.
			       -delta)/ANODE_CATHODE_SPACING;
		float lambda2=(((float)i+0.5)*STRIP_SPACING-STRIP_GAP/2.
			       -delta)/ANODE_CATHODE_SPACING;
		
		fdc->chambers->in[0].cathodePlanes->in[0].strips
		  ->in[i+NODES].u = u+(float)i*STRIP_SPACING; 
		fdc->chambers->in[0].cathodePlanes->in[0].strips
		  ->in[i+NODES].hits = hits = make_s_Hits(1);
		hits->mult=1;
		
		/* Induce charge on the strips according to the Mathieson 
		   function tuned to results from FDC prototype */
		dE=fdc->chambers->in[0].cathodePlanes->in[0].strips
		  ->in[i+NODES].hits->in[0].dE =q_anode*amp_gain/1000.
		  *(tanh(M_PI*K2*lambda2/4.)-tanh(M_PI*K2*lambda1/4.))/4.;
		 
		stripCount++;
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
                   int* track, int* stack)
{
   hitForwardDC(xin,xout,pin,pout,*dEsum,*track,*stack);
}


/* pick and package the hits for shipping */

s_ForwardDC_t* pickForwardDC ()
{
   s_ForwardDC_t* box;
   s_ForwardDC_t* item;

   if ((stripCount == 0) && (wireCount == 0) && (pointCount == 0))
   {
      return 0;
   }

   box = make_s_ForwardDC();
   box->chambers = make_s_Chambers(32);
   while (item = (s_ForwardDC_t*) pickTwig(&forwardDCTree))
   {
      int module = item->chambers->in[0].module;
      int layer = item->chambers->in[0].layer;
      int m = box->chambers->mult;

      if ((m == 0) || (module > box->chambers->in[m-1].module)
                   || (layer  > box->chambers->in[m-1].layer))
      {
         box->chambers->in[m] = item->chambers->in[0];
         box->chambers->in[m].cathodePlanes = make_s_CathodePlanes(32);
         box->chambers->in[m].anodePlanes = make_s_AnodePlanes(32);
         box->chambers->mult++;
      }
      else
      {
         m--;
      }
      if (item->chambers->in[0].cathodePlanes)
      {
         float z = item->chambers->in[0].cathodePlanes->in[0].z;
         int mm = box->chambers->in[m].cathodePlanes->mult;

         if ((mm == 0) ||
             (z > box->chambers->in[m].cathodePlanes->in[mm-1].z + 0.5))
         {
            box->chambers->in[m].cathodePlanes->in[mm] =
                           item->chambers->in[0].cathodePlanes->in[0];
            box->chambers->in[m].cathodePlanes->in[mm].strips =
                           make_s_Strips(stripCount);
            box->chambers->in[m].cathodePlanes->mult++;
         }
         else
         {
            mm--;
         }

         {
            int mmm,i=0;
	    int mult=box->chambers->in[m].cathodePlanes->in[mm].strips->mult;
	    int mult_per_track=item->chambers->in[0].cathodePlanes->in[0].strips->mult;
	    /* Loop over strips... */
	    for (mmm=0;mmm<mult_per_track;mmm++){
	      float u=item->chambers->in[0].cathodePlanes->in[0].strips->in[mmm].u;
	      float dE=item->chambers->in[0].cathodePlanes->in[0].strips->in[mmm].hits->in[0].dE;
	      int j, multi_hit=0;
	     
	      /* Look for strips that have already fired in the current plane*/
	      for (j=0;j<mult;j++){
		if (u==box->chambers->in[m].cathodePlanes->in[mm].strips->in[j].u){
		  multi_hit=1;
		  break;
		}
	      }
	      	      
	      if (!multi_hit){ /* Only one hit in the strip so far */
		box->chambers->in[m].cathodePlanes->in[mm].strips->in[mult+i] =
		  item->chambers->in[0].cathodePlanes->in[0].strips->in[mmm];
		box->chambers->in[m].cathodePlanes->in[mm].strips->mult++;
		i++;
	      }
	      else{  
		/* For the multi-hit case:  currently: just add the charge 
		   from each hit. */ 
		box->chambers->in[m].cathodePlanes->in[mm].strips->in[j].hits->in[0].dE+=dE;
	      }
	      
	    }
	 
         }
         FREE(item->chambers->in[0].cathodePlanes->in[0].strips);
         FREE(item->chambers->in[0].cathodePlanes);
      }
      else if (item->chambers->in[0].anodePlanes)
      {
         float z = item->chambers->in[0].anodePlanes->in[0].z;
         int mm = box->chambers->in[m].anodePlanes->mult;
         if ((mm == 0) ||
             (z > box->chambers->in[m].anodePlanes->in[mm-1].z + 0.5))
         {
            box->chambers->in[m].anodePlanes->in[mm] =
                           item->chambers->in[0].anodePlanes->in[0];
            box->chambers->in[m].anodePlanes->in[mm].wires =
                           make_s_Wires(wireCount);
            box->chambers->in[m].anodePlanes->mult++;  
	   
	 }
         else
         {
            mm--;
         }
         {
            int mmm = box->chambers->in[m].anodePlanes->in[mm].wires->mult;
            if ((mmm == 0) || item->chambers->in[0].anodePlanes->in[0].wires
                                                          ->in[0].hits)
            {
               box->chambers->in[m].anodePlanes->in[mm].wires->in[mmm] =
                     item->chambers->in[0].anodePlanes->in[0].wires->in[0];
               box->chambers->in[m].anodePlanes->in[mm].wires->mult++;
            }
            else if (item->chambers->in[0].anodePlanes->in[0].wires
                                                          ->in[0].fdcPoints)
            {
               int mmmm;
               if (box->chambers->in[m].anodePlanes->in[mm].wires
                                                  ->in[mmm-1].fdcPoints == 0)
               {
                  box->chambers->in[m].anodePlanes->in[mm].wires
                     ->in[mmm-1].fdcPoints = make_s_FdcPoints(pointCount);
               }
               mmmm = box->chambers->in[m].anodePlanes->in[mm].wires
                                               ->in[mmm-1].fdcPoints->mult++;
               box->chambers->in[m].anodePlanes->in[mm].wires
                   ->in[mmm-1].fdcPoints->in[mmmm] =
                     item->chambers->in[0].anodePlanes
                         ->in[0].wires->in[0].fdcPoints->in[0];
               FREE(item->chambers->in[0].anodePlanes->in[0].wires
                                  ->in[0].fdcPoints);
            }
         }
         FREE(item->chambers->in[0].anodePlanes->in[0].wires);
         FREE(item->chambers->in[0].anodePlanes);
      }
      FREE(item->chambers);
      FREE(item);
   }
   stripCount = wireCount = pointCount = 0;

   return box;
}
