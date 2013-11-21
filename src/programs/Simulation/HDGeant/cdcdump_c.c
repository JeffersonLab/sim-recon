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
#include <HDDM/hddm_s.h>

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
   s_CdcTruthPoints_t *points;
   hits = event->physicsEvents->in[0].hitView;
   if (hits == HDDM_NULL ||
       hits->centralDC == HDDM_NULL ||
       hits->centralDC->cdcTruthPoints == HDDM_NULL) {
      return 0;
   }
   printf("New event number %d,",event->physicsEvents->in[0].eventNo);
   printf(" run number %d\n",event->physicsEvents->in[0].runNo);
   points = hits->centralDC->cdcTruthPoints;
	printf(" found %d cdcTruthPoints!\n",points->mult);
	int ipoint;
   for (ipoint=0; ipoint<points->mult; ipoint++) {
		s_CdcTruthPoint_t *point = &points->in[ipoint];
      if (fabs(point->dradius-19.5) < 0.5e5) {
         printf("  dradius=%f,",point->dradius);
         printf("  phi=%f,",point->phi);
         printf("  primary=%s,",point->primary ? "true":"false");
         printf("  r=%f,",point->r);
         printf("  track=%d,",point->track);
         printf("  z=%f,",point->z);
         printf("  dE/dx=%f\n",point->dEdx * 1e6);
      }
   }
   return 1;
}
