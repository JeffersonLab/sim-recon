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

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define U0		-60.
#define WIRE_SPACING 	1.0
#define STRIP_SPACING	1.0
#define TWO_HIT_RESOL	250.
#define MAX_HITS 	10

binTree_t* forwardDCTree = 0;
static int stripCount = 0;
static int wireCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitForwardDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   getlocalcoord_(x,xlocal,"local",5);
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
      int module = getmodule_();
      int layer = getlayer_();
      int plane = getplane_();
      int player = (layer -1)*3 + plane;
      const float tau[] = {60,0,-60,120,60,0,0,-60,-120};
      int iwire = (xlocal[0] - U0)/WIRE_SPACING +1;
      if (layer != 2)		/* cathode strip plane */
      {
         int mark = (module << 24) + (player << 20) + (iwire << 8);
         void** twig = getTwig(&forwardDCTree, mark);
         if (*twig == 0)
         {
            s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
            fdc->chambers = make_s_Chambers(1);
            fdc->chambers->mult = 1;
            fdc->chambers->in[0].module = module;
            fdc->chambers->in[0].layer = layer;
            fdc->chambers->in[0].cathodePlanes = make_s_CathodePlanes(1);
            fdc->chambers->in[0].cathodePlanes->mult = 1;
            fdc->chambers->in[0].cathodePlanes->in[0].z = x[2];
            fdc->chambers->in[0].cathodePlanes->in[0].tau = tau[player];
            fdc->chambers->in[0].cathodePlanes->in[0].strips = make_s_Strips(1);
            fdc->chambers->in[0].cathodePlanes->in[0].strips->mult = 1;
            fdc->chambers->in[0].cathodePlanes->in[0].strips->in[0].u = 
                                                                    xlocal[0];
            fdc->chambers->in[0].cathodePlanes->in[0].strips->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
            stripCount++;
         }
         else
         {
            s_ForwardDC_t* fdc = *twig;
            hits = fdc->chambers->in[0].cathodePlanes->in[0].strips->in[0].hits;
         }
      }
      else			/* anode drift cell plane */
      {
         int mark = (module << 24) + (player << 20) + (iwire << 8);
         void** twig = getTwig(&forwardDCTree, mark);
         if (*twig == 0)
         {
            s_ForwardDC_t* fdc = *twig = make_s_ForwardDC();
            fdc->chambers = make_s_Chambers(1);
            fdc->chambers->mult = 1;
            fdc->chambers->in[0].module = module;
            fdc->chambers->in[0].layer = layer;
            fdc->chambers->in[0].anodePlanes = make_s_AnodePlanes(1);
            fdc->chambers->in[0].anodePlanes->mult = 1;
            fdc->chambers->in[0].anodePlanes->in[0].z = x[2];
            fdc->chambers->in[0].anodePlanes->in[0].tau = tau[player];
            fdc->chambers->in[0].anodePlanes->in[0].wires = make_s_Wires(1);
            fdc->chambers->in[0].anodePlanes->in[0].wires->mult = 1;
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].u = xlocal[0];
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
            wireCount++;
         }
         else
         {
            s_ForwardDC_t* fdc = *twig;
            hits = fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].hits;
         }
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
         hits->in[nhit].t = (hits->in[nhit].t * hits->in[nhit].dE + t*dEsum)
                          / (hits->in[nhit].dE += dEsum);
      }
      else if (nhit < MAX_HITS)              /* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].dE = dEsum;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitForwardDC: ");
         fprintf(stderr,"max hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }

   /* post the hit to the truth tree, once per primary track per cell */
      if (layer == 2)
      {
         s_FdcPoints_t* points;
         int mark = (module << 24) + (player << 20) + (iwire << 8) + track;
         void** twig = getTwig(&forwardDCTree, mark);
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
            fdc->chambers->in[0].anodePlanes->in[0].z = x[2];
            fdc->chambers->in[0].anodePlanes->in[0].tau = tau[player];
            fdc->chambers->in[0].anodePlanes->in[0].wires = make_s_Wires(1);
            fdc->chambers->in[0].anodePlanes->in[0].wires->mult = 1;
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].u = xlocal[0];
            fdc->chambers->in[0].anodePlanes->in[0].wires->in[0].fdcPoints =
                                                                    points;
            points->in[0].x = x[0];
            points->in[0].y = x[1];
            points->in[0].z = x[2];
            points->in[0].dradius = fabs((xlocal[0] - U0)
                                       - (iwire - 0.5)*WIRE_SPACING);
            points->in[0].dEdx = dEdx;
            points->mult++;
            pointCount++;
         }
      }
   }
}

/* entry points from fortran */

void hitforwarddc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum, int* track)
{
   hitForwardDC(xin,xout,pin,pout,*dEsum,*track);
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
            int mmm = box->chambers->in[m].cathodePlanes->in[mm].strips->mult++;
            box->chambers->in[m].cathodePlanes->in[mm].strips->in[mmm] =
                     item->chambers->in[0].cathodePlanes->in[0].strips->in[0];
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
