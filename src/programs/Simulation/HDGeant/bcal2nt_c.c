/*
 * bcal2nt - an simple program for accessing the bcal information
 *           of events stored in a hddm file, and storing them in
 *           a paw ntuple.
 *
 * Richard Jones
 * GlueX collaboration
 * May 15, 2005
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <HDDM/hddm_s.h>

struct bcalnt_struct {
   int event;
   int module;
   int layer;
   int sector;
   int nhit;
   float t[200];
   float E[200];
} bcalnt;

#define BCALNT_FORM "event:i,module:i,layer:i,sector:i,nhit[0,200]:i,t(nhit):r,E(nhit):r,z(nhit):r"

#define PAWC_SIZE 1000000
struct {
   float q[PAWC_SIZE];
} pawc_;

struct {
   int iq[100];
} quest_;

void hlimit(int size)
{
   void hlimit_(int *words);
   hlimit_(&size);
}
void hbset(char* name, int *value, int *istat)
{
   void hbset_(char *,int *,int *,int);
   hbset_(name, value, istat, strlen(name));
}
void hropen(int lun, char *name, char*filename, char*status, int *lrec, int *istat)
{
   void hropen_(int *,char *,char *,char *,int *,int *,int,int,int);
   hropen_(&lun, name, filename, status, lrec, istat, strlen(name),
           strlen(filename), strlen(status));
}
void hbnt(int id,char*chtitle,char*chopt)
{
   void hbnt_(int *id ,char* name, char* chmod,int ,int);
   hbnt_(&id,chtitle,chopt,strlen(chtitle),strlen(chopt));
}
void hbname(int id,char*chblok,void*variable,char*chform)
{
   void hbname_(int *id, char* chblok, void*variable, char*chform, int,int);
   hbname_(&id,chblok,variable,chform,strlen(chblok),strlen(chform));
}
void hfnt(int id)
{
   void hfnt_(int*id);
   hfnt_(&id);
}
void hrend(char*filename)
{
   void hrend_(char *,int);
   hrend_(filename, strlen(filename));
   return;
}
void hrout(int num, int *icycle, char*opt)
{
   void hrout_(int *,int *,char *,int);
   hrout_(&num, icycle, opt, strlen(opt));
}

int process_event(s_HDDM_t *event);

int main(int argc, char **argv)
{
   s_HDDM_t *thisInputEvent = 0;
   s_iostream_t *thisInputFile = 0;
   int input;
   int lrec=65536;
   int status;
   int cycle;

   hlimit(PAWC_SIZE);
   hbset("BSIZE",&lrec,&status);
   quest_.iq[9] = 256000;  // extend RZ quota to 2^32 bits
   hropen(50,"RZfile","bcal2nt.hbook","NQE",&lrec,&status);
   hbnt(1,"BCal diagnostic ntuple"," ");
   hbname(1,"bcalnt",&bcalnt,BCALNT_FORM);

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
   hrout(1,&cycle," ");
   hrend("RZfile");
}

int process_event(s_HDDM_t *event)
{
   s_HitView_t *hits;
   if ((hits = event->physicsEvents->in[0].hitView) == HDDM_NULL) {
      fprintf(stderr,"no hits information in this file, quitting!\n");
      exit(1);
   }
   
   if (hits->barrelEMcal != HDDM_NULL) {
      s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
      int cell;
      for (cell=0; cell < cells->mult; cell++) {
         s_BcalHits_t *hits = cells->in[cell].bcalHits;
         int hit;
         bcalnt.event = event->physicsEvents->in[0].eventNo;
         bcalnt.module = cells->in[cell].module;
         bcalnt.layer = cells->in[cell].layer;
         bcalnt.sector = cells->in[cell].sector;
         for (hit=0; hit < hits->mult; hit++) {
            bcalnt.t[hit] = hits->in[hit].t; 
            bcalnt.E[hit] = hits->in[hit].E; 
         }
         bcalnt.nhit = hit;
         hfnt(1);
      }
      return 1;
   }
   else {
      return 0;
   }
}
