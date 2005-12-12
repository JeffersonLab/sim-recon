/*
 * cdcdump - an example program for accessing the contents
 *           of events stored in a hddm file.
 *
 * Richard Jones
 * GlueX collaboration
 * January 10, 2005
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "hddm_s.h"

int process_event(s_HDDM_t *event);

int main(int argc, char **argv)
{
   s_HDDM_t *thisInputEvent = 0;
   s_iostream_t *thisInputFile = 0;
   int input;
   for (input=1; input<argc; input++) {
      if (! (thisInputFile = open_s_HDDM(argv[input]))) {
         fprintf(stderr,"Error - could not open input file %s\n",
                 argv[input]);
         exit(1);
      }
      while (thisInputEvent = read_s_HDDM(thisInputFile)) {
         process_event(thisInputEvent);
         flush_s_HDDM(thisInputEvent,0);
      }
      close_s_HDDM(thisInputFile);
   }
}

int process_event(s_HDDM_t *event)
{
   s_HitView_t *hits;
   s_Rings_t *rings;
   int ring;
   hits = event->physicsEvents->in[0].hitView;
   if (hits == 0 ||
       hits->centralDC == 0 ||
       hits->centralDC->rings == 0) {
      return 0;
   }
   printf("New event number %d,",event->physicsEvents->in[0].eventNo);
   printf(" run number %d\n",event->physicsEvents->in[0].runNo);
   rings = hits->centralDC->rings;
   for (ring=0; ring<rings->mult; ring++) {
      if (fabs(rings->in[ring].radius-19.5) < 0.5e5) {
         s_Straws_t *straws = rings->in[ring].straws;
         int straw;
         for (straw=0; straw<straws->mult; straw++) {
            s_CdcPoints_t *points = straws->in[straw].cdcPoints;
            if (points == 0) {
               printf(" orphan found!\n");
            }
            else if (points->mult != 1) {
               printf(" found %d cdcPoints!\n",points->mult);
            }
            else if (straws->in[straw].hits == 0) {
               printf(" widow found!\n");
            }
            else {
               printf("  straw hit at phi=%f,",straws->in[straw].phim);
               printf("  drift time=%f,",straws->in[straw].hits->in[0].t);
               printf("  dE/dx=%f\n",points->in[0].dEdx);
            }
         }
      }
   }
}
