/*
 * hitCerenkov - registers hits for Cerenkov counter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitCerenkov
 *
 *          Oct 12 2012, yqiang, removed changes made in revision 9720
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

extern s_HDDM_t* thisInputEvent;

#define TWO_HIT_RESOL   50.
#define MAX_HITS        100
#define THRESH_PE       2
#define OPTICAL_PHOTON  50

binTree_t* cerenkovTree = 0;
static int sectionCount = 0;
static int pointCount = 0;


/* register truth points during tracking (from gustep) */

void hitCerenkov (float xin[4], float xout[4],
                  float pin[5], float pout[5], float dEsum,
                  int track, int stack, int history, int ipart)
{
   float x[3], t;

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;

   /* post the hit to the truth tree */

   if ((history == 0) && (dEsum > 0))
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&cerenkovTree, mark);
      if (*twig == 0)
      {
         s_Cerenkov_t* cere = *twig = make_s_Cerenkov();
         s_CereTruthPoints_t* points = make_s_CereTruthPoints(1);
         cere->cereTruthPoints = points;
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (stack <= a);
         points->in[0].track = track;
         points->in[0].x = xin[0];
         points->in[0].y = xin[1];
         points->in[0].z = xin[2];
         points->in[0].t = xin[3]*1e9;
         points->in[0].px = pin[4]*pin[0];
         points->in[0].py = pin[4]*pin[1];
         points->in[0].pz = pin[4]*pin[2];
         points->in[0].E = pin[3];
         points->in[0].ptype = ipart;
         points->mult = 1;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   if (dEsum < 0)  		/* indicates a detector Cerenkov photon */
   {
      int nshot;
      s_CereHits_t* hits;
      int sector = getsector_wrapper_();
      float pe = 1;
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
}

/* entry points from fortran */

void hitcerenkov_(float* xin, float* xout,
                  float* pin, float* pout, float* dEsum,
                  int* track, int* stack, int* history, int* ipart)
{
   hitCerenkov(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
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
      int section;
      s_CereTruthPoints_t* points = item->cereTruthPoints;
      int point;

      for (section=0; section < sections->mult; ++section)
      {
         s_CereHits_t* hits = sections->in[section].cereHits;

      /* compress out the hits below threshold */
         int iok,i;
         for (iok=i=0; i < hits->mult; i++)
         {
           if (hits->in[i].pe >= THRESH_PE)
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
            int m = box->cereSections->mult++;
            box->cereSections->in[m] = sections->in[section];
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (sections != HDDM_NULL)
      {
         FREE(sections);
      }

      for (point=0; point < points->mult; ++point)
      {
         int m = box->cereTruthPoints->mult++;
         box->cereTruthPoints->in[m] = points->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   sectionCount = pointCount = 0;

   if ((box->cereSections != HDDM_NULL) &&
       (box->cereSections->mult == 0))
   {
      FREE(box->cereSections);
      box->cereSections = HDDM_NULL;
   }
   if ((box->cereTruthPoints != HDDM_NULL) &&
       (box->cereTruthPoints->mult == 0))
   {
      FREE(box->cereTruthPoints);
      box->cereTruthPoints = HDDM_NULL;
   }
   if ((box->cereSections->mult == 0) &&
       (box->cereTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
