/*
 * hitCerenkov - registers hits for Cerenkov counter
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

#define N0_FIGURE	60
#define REFR_INDEX	1.0017
#define TWO_HIT_RESOL   50.
#define MAX_HITS        100
#define THRESH_PE       2

binTree_t* cerenkovTree = 0;
static int sectionCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitCerenkov (float xin[4], float xout[4],
                  float pin[5], float pout[5], float dEsum,
                  int track, int stack)
{
   float x[3], t;
   float dx[3], dr;
   float xcere[3];
   float xHat[] = {1,0,0};
   double Eave = (pin[3] + pout[3])/2;
   double pave = (pin[4] + pout[4])/2;
   double beta = pave/Eave;
   double costheta = 1/(beta*REFR_INDEX);

   if (costheta > 1) return;		/* no light below threshold */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(xHat,"local",xcere,"CERE");
   dx[0] = xin[0] - xout[0];
   dx[1] = xin[1] - xout[1];
   dx[2] = xin[2] - xout[2];
   dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

   /* post the hit to the hits tree, mark sector as hit */
   {
      int nshot;
      s_CereHits_t* hits;
      int sector = getsector_();
      float phim = atan2(xcere[1],xcere[0]);
      double sintheta2 = 1 - costheta*costheta;
      float pe = N0_FIGURE * sintheta2 * dr;
      int mark = sector;
      void** twig = getTwig(&cerenkovTree, mark);
      if (*twig == 0)
      {
         s_Cerenkov_t* cere = *twig = make_s_Cerenkov();
         s_CereSections_t* sections = make_s_CereSections(1);
         sections->mult = 1;
         sections->in[0].sector = sector;
         sections->in[0].cereHits = hits = make_s_CereHits(MAX_HITS);
         cere->cereSections = sections;
         sectionCount++;
      }
      else
      {
         s_Cerenkov_t* cere = *twig;
         hits = cere->cereSections->in[0].cereHits;
      }

      for (nshot = 0; nshot < hits->mult; nshot++)
      {
         if (fabs(hits->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < hits->mult)            /* merge with former hit */
      {
         hits->in[nshot].t = (hits->in[nshot].t * hits->in[nshot].pe + t*pe)
                            / (hits->in[nshot].pe += pe);
      }
      else if (nshot < MAX_HITS)         /* create new shot */
      {
         hits->in[nshot].t = t;
         hits->in[nshot].pe = pe;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitCerenkov: ");
         fprintf(stderr,"max shot count %d exceeded, truncating!\n",MAX_HITS);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      int mark = (1<<30) + track;
      void** twig = getTwig(&cerenkovTree, mark);
      if (*twig == 0)
      {
         s_Cerenkov_t* cere = *twig = make_s_Cerenkov();
         s_CereTruthPoints_t* points = make_s_CereTruthPoints(1);
         cere->cereTruthPoints = points;
         float E = (pin[3]+pout[3])/2;
         float p = (pin[4]+pout[4])/2;
         points->in[0].primary = (stack == 0);
         points->in[0].track = track;
         points->in[0].x = x[0];
         points->in[0].y = x[1];
         points->in[0].z = x[2];
         points->in[0].t = t;
         points->in[0].px = p*(pin[0]+pout[0])/2;
         points->in[0].py = p*(pin[1]+pout[1])/2;
         points->in[0].pz = p*(pin[2]+pout[2])/2;
         points->in[0].E = E;
         points->mult = 1;
         pointCount++;
      }
   }
}

/* entry point from fortran */

void hitcerenkov_(float* xin, float* xout,
                  float* pin, float* pout, float* dEsum,
                  int* track, int* stack)
{
   hitCerenkov(xin,xout,pin,pout,*dEsum,*track,*stack);
}


/* pick and package the hits for shipping */

s_Cerenkov_t* pickCerenkov ()
{
   s_Cerenkov_t* box;
   s_Cerenkov_t* item;

   if ((sectionCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_Cerenkov();
   box->cereSections = make_s_CereSections(sectionCount);
   box->cereTruthPoints = make_s_CereTruthPoints(pointCount);
   while (item = pickTwig(&cerenkovTree))
   {
      s_CereSections_t* sections = item->cereSections;
      s_CereTruthPoints_t* points = item->cereTruthPoints;

      if (sections != HDDM_NULL)
      {
      /* compress out the hits below threshold */
         s_CereHits_t* hits = sections->in[0].cereHits;
         if (hits != HDDM_NULL)
         {
            int mok,i;
            for (mok=i=0; i < hits->mult; i++)
            {
              if (hits->in[i].pe >= THRESH_PE)
              {
                if (mok < i)
                {
                   hits->in[mok] = hits->in[i];
                }
                ++mok;
              }
            }
            hits->mult = mok;
            if (mok)
            {
               int m = box->cereSections->mult++;
               box->cereSections->in[m] = sections->in[0];
            }
            else
            {
               FREE(hits);
            }
         }
         FREE(sections);
      }
      else if (points != HDDM_NULL)
      {
         int m = box->cereTruthPoints->mult++;
         box->cereTruthPoints->in[m] = points->in[0];
         FREE(points);
      }
      FREE(item);
   }
   sectionCount = pointCount = 0;
   return box;
}
