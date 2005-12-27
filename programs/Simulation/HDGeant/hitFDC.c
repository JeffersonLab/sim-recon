/*
 * hitFDC - registers hits for forward drift chambers
 *
 *        This is a part of the hits package for the
 *        HDGeant simulation program for Hall D.
 *
 *        version 1.0         -Richard Jones July 16, 2001
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

const float Tau[] = {0,-45,0,45,15,60,105,-105,-60,-15};

#define DRIFT_SPEED        .0022
#define STRIP_SPACING        0.5
#define STRIP_GAP       0.1
#define ANODE_CATHODE_SPACING 0.5
#define TWO_HIT_RESOL        250.
#define WIRE_SPACING         1.0
#define MAX_HITS         100
#define K2              1.15
#define STRIP_NODES     3
#define MAX_STRIPS      240
#define THRESH_KEV        1.
#define THRESH_STRIPS     5.      /* mV */
#define ELECTRON_CHARGE 1.6022e-4 /* fC */

binTree_t* forwardDCTree = 0;
static int stripCount = 0;
static int wireCount = 0;
static int pointCount = 0;

#define U0                100.

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
  int i;
  

  if (dEsum == 0) return;              /* only seen if it deposits energy */

  
  {
    /* Get chamber information */
    int module = getmodule_();
    int layer = getlayer_();
    int chamber = (module*10)+layer;
    
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
      s_FdcAnodeHits_t* ahits;    
      s_FdcCathodeHits_t* chits;    
      int wire = (xlocal[0] + U0)/WIRE_SPACING +1;
      float dradius = fabs((xlocal[0] + U0) - (wire - 0.5)*WIRE_SPACING);
      float tdrift = t + dradius/DRIFT_SPEED;

      /* first record the anode wire hit */
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
          wires->in[0].fdcAnodeHits = ahits = make_s_FdcAnodeHits(MAX_HITS);
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
           ahits = fdc->fdcChambers->in[0].fdcAnodeWires->in[0].fdcAnodeHits;
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
          ahits->in[nhit].t = 
              (ahits->in[nhit].t * ahits->in[nhit].dE + tdrift * dEsum)
              / (ahits->in[nhit].dE += dEsum);
        }
        else if (nhit < MAX_HITS)              /* create new hit */
        {
          ahits->in[nhit].t = tdrift;
          ahits->in[nhit].dE = dEsum;
          ahits->mult++;
        }
        else
        {
          fprintf(stderr,"HDGeant error in hitForwardDC: ");
          fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
        }
      }

      /* then generate hits in the two surrounding cathode planes */

      {
        float avalanche_x = (floor(xlocal[0]/WIRE_SPACING)+0.5)*WIRE_SPACING;
        float avalanche_y = xlocal[1];

        /* Mock-up of cathode strip charge distribution */

        int plane, node;
        for (plane=1; plane<4; plane+=2)
        {
          float theta = (plane == 1)? -M_PI/4 : +M_PI/4;
          float cathode_u = avalanche_x*cos(theta)+avalanche_y*sin(theta);
          int strip1 = ceil((cathode_u + U0)/STRIP_SPACING);
          float cathode_u1 = (strip1 - 0.5)*STRIP_SPACING - U0;
          float delta = cathode_u-cathode_u1;

          /* Variables for approximating number of ion pairs  */ 
          const float w_eff=26.0; // eV, average energy needed to produce an 
                                  // ion pair.  For simplicity use pure argon
          const float n_s_per_p=2.2; // average number of secondary ion pairs
                                     // per primary ionization
          float n_p_mean,n_s_mean;
          int n_p,n_t,n_s,err;
      
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
             that requires calling rnpssn_ twice.
           */
          n_p_mean = dEsum/w_eff/(1.+n_s_per_p)*1e9;
          rnpssn_(&n_p_mean,&n_p,&err);
          n_s_mean = ((float)n_p)*n_s_per_p;
          rnpssn_(&n_s_mean,&n_s,&err);
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
            float dE = q_anode*amp_gain/1000.
                       *(tanh(M_PI*K2*lambda2/4.)-tanh(M_PI*K2*lambda1/4.))/4.;

            int strip = strip1+node;
            if ((strip > 0) && (strip <= MAX_STRIPS))
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
                strips->in[0].fdcCathodeHits = chits
                                             = make_s_FdcCathodeHits(MAX_HITS);
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
                                        ->in[0].fdcCathodeHits;
              }
          
              for (nhit = 0; nhit < ahits->mult; nhit++)
              {
                if (fabs(chits->in[nhit].t - tdrift) < TWO_HIT_RESOL)
                {
                  break;
                }
              }
              if (nhit < chits->mult)          /* merge with former hit */
              {
                chits->in[nhit].t = 
                    (chits->in[nhit].t * chits->in[nhit].dE + tdrift * dE)
                    / (chits->in[nhit].dE += dE);
              }
              else if (nhit < MAX_HITS)        /* create new hit */
              {
                chits->in[nhit].t = tdrift;
                chits->in[nhit].dE = dE;
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

 /* post the hit to the truth tree, once per primary track per anode plane */
 
    {
      int wire = (xlocal[0] + U0)/WIRE_SPACING +1;
      float dradius = fabs((xlocal[0] + U0) - (wire - 0.5)*WIRE_SPACING);
      int mark = (chamber<<20) + (track<<14);
      void** twig = getTwig(&forwardDCTree, mark);
      if (*twig == 0)
      {
        s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
        s_FdcChambers_t* chambers = make_s_FdcChambers(1);
        s_FdcTruthPoints_t* points = make_s_FdcTruthPoints(1);
        points->mult = 1;
        points->in[0].primary = (stack == 0);
        points->in[0].track = track;
        points->in[0].x = x[0];
        points->in[0].y = x[1];
        points->in[0].z = x[2];
        points->in[0].dradius = dradius;
        points->in[0].dEdx = dEdx;
        chambers->mult = 1;
        chambers->in[0].module = module;
        chambers->in[0].layer = layer;
        chambers->in[0].fdcTruthPoints = points;
        fdc->fdcChambers = chambers;
        pointCount++;
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
      s_FdcCathodeStrips_t* strips = chambers->in[0].fdcCathodeStrips;
      s_FdcTruthPoints_t* points = chambers->in[0].fdcTruthPoints;
      int iwire,istrip;
      int mok=0;
      for (iwire=0; iwire < wires->mult; iwire++)
      {
         s_FdcAnodeHits_t* ahits = wires->in[0].fdcAnodeHits;
         if (ahits != HDDM_NULL)
         {
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
            ahits->mult = iok;
            if (iok == 0)
            {
              FREE(ahits);
              FREE(wires);
              wires = HDDM_NULL;
            }
         }
      }
      for (istrip=0; istrip < strips->mult; istrip++)
      {
         s_FdcCathodeHits_t* chits = strips->in[0].fdcCathodeHits;
         if (chits != HDDM_NULL)
         {
            int i,iok;
            for (iok=i=0; i < chits->mult; i++)
            {
               if (chits->in[i].dE >= THRESH_STRIPS)
               {
                  if (iok < i)
                  {
                     chits->in[iok] = chits->in[i];
                  }
                  ++iok;
                  ++mok;
               }
            }
            chits->mult = iok;
            mok += iok;
            if (iok == 0)
            {
              FREE(chits);
              FREE(strips);
              strips = HDDM_NULL;
            }
         }
      }

      if (mok || points != HDDM_NULL)
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
        if (strips != HDDM_NULL)
        {
          int mm = box->fdcChambers->in[m].fdcCathodeStrips->mult++;
          box->fdcChambers->in[m].fdcCathodeStrips->in[mm] = strips->in[0];
          FREE(strips);
        }
        else if (wires != HDDM_NULL)
        {
          int mm = box->fdcChambers->in[m].fdcAnodeWires->mult++;
          box->fdcChambers->in[m].fdcAnodeWires->in[mm] = wires->in[0];
          FREE(wires);
       }
       else if (points != HDDM_NULL)
       {
         int mm = box->fdcChambers->in[m].fdcTruthPoints->mult++;
         box->fdcChambers->in[m].fdcTruthPoints->in[mm] = points->in[0];
         FREE(points);
       }
     }
     FREE(chambers);
     FREE(item);
   }
   stripCount = wireCount = pointCount = 0;

   return box;
}
