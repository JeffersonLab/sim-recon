/*
 * hitLGD - registers hits for forward calorimeter
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

#define BLOCK_WIDTH	4.0
#define THRESH_MEV	150.
#define TWO_HIT_RESOL   75.
#define MAX_HITS        10

binTree_t* forwardEMcalTree = 0;
static int blockCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitForwardEMcal (float xin[4], float xout[4],
                      float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float xlocal[3];

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   getlocalcoord_(x,xlocal,"FCAL",4);

   /* post the hit to the hits tree, mark block as hit */
   {
      int nshot;
      s_Showers_t* shots;
      int row = getrow_();
      int column = getcolumn_();
      float xcol = (floor(xlocal[0]/BLOCK_WIDTH) + 0.5)*BLOCK_WIDTH;
      float yrow = (floor(xlocal[1]/BLOCK_WIDTH) + 0.5)*BLOCK_WIDTH;
      int mark = (row << 16) + column;
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->rows = make_s_Rows(1);
         cal->rows->mult = 1;
         cal->rows->in[0].y = yrow;
         cal->rows->in[0].columns = make_s_Columns(1);
         cal->rows->in[0].columns->mult = 1;
         cal->rows->in[0].columns->in[0].x = xcol;
         cal->rows->in[0].columns->in[0].showers =
         shots = make_s_Showers(MAX_HITS);
         blockCount++;
      }
      else
      {
         s_ForwardEMcal_t* cal = *twig;
         shots = cal->rows->in[0].columns->in[0].showers;
      }

      for (nshot = 0; nshot < shots->mult; nshot++)
      {
         if (fabs(shots->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < shots->mult)		/* merge with former hit */
      {
         shots->in[nshot].t =
                       (shots->in[nshot].t * shots->in[nshot].E + t*dEsum)
                     / (shots->in[nshot].E += dEsum);
      }
      else if (nshot < MAX_HITS)         /* create new hit */
      {
         shots->in[nshot].t =
                       (shots->in[nshot].t * shots->in[nshot].E + t*dEsum)
                     / (shots->in[nshot].E += dEsum);
         shots->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitforwardEMcal: ");
         fprintf(stderr,"maximum hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      s_ForwardShowers_t* showers;
      int mark = (track << 24);
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->forwardShowers = showers = make_s_ForwardShowers(1);
         showers->in[0].t = t;
         showers->in[0].x = x[0];
         showers->in[0].y = x[1];
         showers->in[0].E = dEsum;
         showers->mult = 1;
         showerCount++;
      }
      else
      {
         showers = ((s_ForwardEMcal_t*) *twig)->forwardShowers;
         showers->in[0].x = (showers->in[0].x * showers->in[0].E + x[0]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].y = (showers->in[0].y * showers->in[0].E + x[1]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].t = (showers->in[0].t * showers->in[0].E + t*dEsum)
                          / (showers->in[0].E += dEsum);
      }
   }
}

/* entry point from fortran */

void hitforwardemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum, int* track)
{
   hitForwardEMcal(xin,xout,pin,pout,*dEsum,*track);
}


/* pick and package the hits for shipping */

s_ForwardEMcal_t* pickForwardEMcal ()
{
   s_ForwardEMcal_t* box;
   s_ForwardEMcal_t* item;

   if ((blockCount == 0) && (showerCount == 0))
   {
      return 0;
   }

   box = make_s_ForwardEMcal();
   box->rows = make_s_Rows(blockCount);
   box->forwardShowers = make_s_ForwardShowers(showerCount);
   while (item = (s_ForwardEMcal_t*) pickTwig(&forwardEMcalTree))
   {
      if (item->rows)
      {
         float E = item->rows->in[0].columns->in[0].showers->in[0].E;
         float y = item->rows->in[0].y;
         int m = box->rows->mult;
         if (m == 0)
         {
            box->rows->in[m] = item->rows->in[0];
            box->rows->in[m].columns = make_s_Columns(blockCount);
            box->rows->mult++;
         }
         else if (y > box->rows->in[m-1].y + 0.5)
         {
            if (box->rows->in[m-1].columns->mult == 0)
            {
               FREE(box->rows->in[--m].columns);
               box->rows->mult--;
            }
            box->rows->in[m] = item->rows->in[0];
            box->rows->in[m].columns = make_s_Columns(blockCount);
            box->rows->mult++;
         }
         else
         {
            m--;
         }
         if (E > THRESH_MEV/1000)
         {
            int mm = box->rows->in[m].columns->mult++;
            box->rows->in[m].columns->in[mm] = 
                                      item->rows->in[0].columns->in[0];
         }
         else
         {
            FREE(item->rows->in[0].columns->in[0].showers);
         }
         FREE(item->rows->in[0].columns);
         FREE(item->rows);
      }
      else if (item->forwardShowers)
      {
         int m = box->forwardShowers->mult++;
         box->forwardShowers->in[m] = item->forwardShowers->in[0];
         FREE(item->forwardShowers);
      }
      FREE(item);
   }

   /* Reduce the last <row>...</row> if empty */
   {
      int m = box->rows->mult;
      if ((m > 0) && box->rows->in[m-1].columns->mult == 0)
      {
         FREE(box->rows->in[--m].columns);
         box->rows->mult--;
      }
   }

   blockCount = showerCount = 0;
   return box;
}
