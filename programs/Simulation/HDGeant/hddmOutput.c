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

s_iostream_t* thisOutputStream = 0;
s_HDDM_t* thisOutputEvent = 0;
extern s_HDDM_t* thisInputEvent;

int openOutput (char* filename)
{
   thisOutputStream = init_s_HDDM(filename);
   return (thisOutputStream == 0);
}

int flushOutput ()
{
   if (thisOutputEvent != 0)
   {
      flush_s_HDDM(thisOutputEvent, thisOutputStream);
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
   if (thisOutputEvent->physicsEvents->in[0].hitView == 0)
   {
      thisOutputEvent->physicsEvents->in[0].hitView = make_s_HitView();
   }
   thisOutputEvent->physicsEvents->in[0].hitView->
      					centralDC = pickCentralDC();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					forwardDC = pickForwardDC();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					startCntr = pickStartCntr();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					barrelEMcal = pickBarrelEMcal();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					Cerenkov = pickCerenkov();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					forwardTOF = pickForwardTOF();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					forwardEMcal = pickForwardEMcal();
   thisOutputEvent->physicsEvents->in[0].hitView->
      					upstreamEMveto = pickUpstreamEMveto();
   return 0;
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
