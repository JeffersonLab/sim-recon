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

#include <fstream>
#include <stdlib.h>
#include "hddm_s.hpp"

extern "C" {
   
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
   void hbset(const char* name, int *value, int *istat)
   {
      void hbset_(const char *,int *,int *,int);
      hbset_(name, value, istat, strlen(name));
   }
   void hropen(int lun, const char *name, const char *filename, const char *status, int *lrec, int *istat)
   {
      void hropen_(int *,const char *,const char *,const char *,int *,int *,int,int,int);
      hropen_(&lun, name, filename, status, lrec, istat, strlen(name),
              strlen(filename), strlen(status));
   }
   void hbnt(int id,const char*chtitle,const char*chopt)
   {
      void hbnt_(int *id ,const char* name, const char* chmod,int ,int);
      hbnt_(&id,chtitle,chopt,strlen(chtitle),strlen(chopt));
   }
   void hbname(int id,const char*chblok,void*variable,const char*chform)
   {
      void hbname_(int *id, const char* chblok, void*variable, const char*chform, int,int);
      hbname_(&id,chblok,variable,chform,strlen(chblok),strlen(chform));
   }
   void hfnt(int id)
   {
      void hfnt_(int*id);
      hfnt_(&id);
   }
   void hrend(const char*filename)
   {
      void hrend_(const char *,int);
      hrend_(filename, strlen(filename));
      return;
   }
   void hrout(int num, int *icycle, const char *opt)
   {
      void hrout_(int *,int *,const char *,int);
      hrout_(&num, icycle, opt, strlen(opt));
   }
   
}

int process_event(hddm_s::HDDM &event);

int main(int argc, char **argv)
{
   int input;
   int lrec=65536;
   int status;
   int cycle=1; // initialize to 1 just to avoid compiler warnings

   hlimit(PAWC_SIZE);
   hbset("BSIZE",&lrec,&status);
   quest_.iq[9] = 256000;  // extend RZ quota to 2^32 bits
   hropen(50,"RZfile","bcal2nt.hbook","NQE",&lrec,&status);
   hbnt(1,"BCal diagnostic ntuple"," ");
   hbname(1,"bcalnt",&bcalnt,BCALNT_FORM);

   hddm_s::HDDM record;
   for (input=1; input<argc; input++) {
      std::ifstream ifs(argv[input]);
      if (!ifs.is_open()) {
         std::cerr << "Error - could not open input file "
                   << argv[input] << std::endl;
         exit(1);
      }
      hddm_s::istream istr(ifs);
      while (ifs.good()) {
         istr >> record;
         process_event(record);
         record.clear();
      }
   }
   hrout(1,&cycle," ");
   hrend("RZfile");
}

int process_event(hddm_s::HDDM &event)
{
   hddm_s::HitViewList views = event.getHitViews();
   if (views.size() == 0) {
      std::cerr << "no hits information in this file, quitting!"
                << std::endl;
      exit(1);
   }
   
   hddm_s::BcalCellList cells = event.getBcalCells();
   hddm_s::BcalCellList::iterator iter;
   for (iter = cells.begin(); iter != cells.end(); ++iter) {
      bcalnt.event = iter->getEventNo();
      bcalnt.module = iter->getModule();
      bcalnt.layer = iter->getLayer();
      bcalnt.sector = iter->getSector();
      hddm_s::BcalHitList hits = iter->getBcalHits();
      hddm_s::BcalHitList::iterator hiter;
      int hit=0;
      for (hiter = hits.begin(); hiter != hits.end(); ++hiter, ++hit) {
         bcalnt.t[hit] = hiter->getT(); 
         bcalnt.E[hit] = hiter->getE(); 
      }
      bcalnt.nhit = hits.size();
      hfnt(1);
   }
   return 1;
}
