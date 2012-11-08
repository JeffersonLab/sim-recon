/*
 * hddmOutput - functions to handle output of simulation results from HDGeant
 *		through the standard hddm i/o mechanism.
 *
 * Interface:
 *	openOutput(filename) - open output stream to file <filename>
 *      loadOutput()  - load output event from hit structures
 *      flushOutput() - flush current event structure to output stream
 *	closeOutput() - close currently open output stream
 *
 * Richard Jones
 * University of Connecticut
 * July 13, 2001
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <hddm_s.h>
#include <hddmOutput.h>

#include "memcheck.h"

extern const char* GetMD5Geom(void);

s_iostream_t* thisOutputStream = 0;
s_HDDM_t* thisOutputEvent = 0;
extern s_HDDM_t* thisInputEvent;

static unsigned int Nevents = 0;

int openOutput (char* filename)
{
   thisOutputStream = init_s_HDDM(filename);
   return (thisOutputStream == 0);
}

int flushOutput ()
{
   if (thisOutputEvent != 0)
   {
      if (flush_s_HDDM(thisOutputEvent, thisOutputStream) != 0) {
         fprintf(stderr,"Fatal error in flushOutput:");
         fprintf(stderr," write failed to hddm output file.\n");
         exit(7);
      }
      thisOutputEvent = 0;
   }
   checkpoint();
   return 0;
}

int closeOutput ()
{
   if (thisOutputStream)
   {
      close_s_HDDM(thisOutputStream);
      thisOutputStream = 0;
   }
   return 0;
}

int loadOutput ()
{
   int packages_hit=0;
   s_HitView_t *hitView;
	
	Nevents++;

   if (thisOutputEvent)
   {
      flush_s_HDDM(thisOutputEvent, 0);
   }

   thisOutputEvent = thisInputEvent;
   thisInputEvent = 0;
   if (thisOutputEvent == 0)
   {
      static int eventNo = 0;
      thisOutputEvent = make_s_HDDM();
      thisOutputEvent->physicsEvents = make_s_PhysicsEvents(1);
      thisOutputEvent->physicsEvents->mult = 1;
      thisOutputEvent->physicsEvents->in[0].eventNo = ++eventNo;
   }
	
	if (Nevents == 1) {
		if (thisOutputEvent->geometry == HDDM_NULL) {
			thisOutputEvent->geometry = make_s_Geometry();
		}
		thisOutputEvent->geometry->md5simulation = strdup(GetMD5Geom());
		thisOutputEvent->geometry->md5smear = strdup("");
		thisOutputEvent->geometry->md5reconstruction = strdup("");
	}
	
   if (thisOutputEvent->physicsEvents->in[0].hitView == HDDM_NULL)
   {
      thisOutputEvent->physicsEvents->in[0].hitView = make_s_HitView();
   }

   hitView = thisOutputEvent->physicsEvents->in[0].hitView;
   if ((hitView->centralDC = pickCentralDC()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->forwardDC = pickForwardDC()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->startCntr = pickStartCntr()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->barrelEMcal = pickBarrelEMcal()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->Cerenkov = pickCerenkov()) != HDDM_NULL) {
      ++packages_hit;
   }
   // add RICH, yqiang, Oct 11 2012
   if ((hitView->RICH = pickRich()) != HDDM_NULL) {
	   ++packages_hit;
   }
   if ((hitView->forwardTOF = pickForwardTOF()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->forwardEMcal = pickForwardEMcal()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->ComptonEMcal = pickComptonEMcal()) != HDDM_NULL) {
      ++packages_hit;
   }
#ifdef TESTING_CAL_CONTAINMENT
   if ((hitView->gapEMcal = pickGapEMcal()) != HDDM_NULL) {
      ++packages_hit;
   }
#endif
   if ((hitView->upstreamEMveto = pickUpstreamEMveto()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->tagger = pickTagger()) != HDDM_NULL) {
      ++packages_hit;
   }
   if ((hitView->mcTrajectory = pickMCTrajectory()) != HDDM_NULL) {
      ++packages_hit;
   }
   return packages_hit;
}

/* entry points from Fortran */

int openoutput_ (char* filename)
{
   int retcode = openOutput(strtok(filename," "));
   return retcode;
}

int flushoutput_ ()
{
   return flushOutput();
}

int loadoutput_ ()
{
   return loadOutput();
}

int closeoutput_ ()
{
   return closeOutput();
}
