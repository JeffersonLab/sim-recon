#include <iostream>
using namespace std;
#include "hddm_v.h"

int main()
{
   char hddmfile[] = "test.hddm";
   v_iostream_t *outputfp = init_v_HDDM(hddmfile);
   if (!outputfp) {
      cerr << "Test failed: error opening output file\n";
      exit(1);
   }

   v_HDDM_t *event;
   v_GenericTag_t *tag;
   event = make_v_HDDM();
   tag = make_v_GenericTag();
   tag->floatTag = make_v_FloatTag();
   tag->doubleTag = make_v_DoubleTag();
   tag->particleTag = make_v_ParticleTag();
   tag->stringTag = make_v_StringTag();
   tag->intTag = make_v_IntTag();
   tag->longTag = make_v_LongTag();
   tag->booleanTag = make_v_BooleanTag();
   tag->anyURITag = make_v_AnyURITag();
   tag->floatTag->pi = 3.141593;
   tag->doubleTag->pi = 3.141592592636;
   tag->particleTag->pi = PiPlus;

   tag->stringTag->quote = (char*)malloc(100);
   strcpy(tag->stringTag->quote,"pass the red quarks, please");
   tag->intTag->magic = 133557799;
   tag->longTag->magic = 133557799002244668LL;
   tag->booleanTag->truth = 1;
   tag->anyURITag->uri = (char*)malloc(100);
   strcpy(tag->anyURITag->uri,"http://portal.gluex.org");
   event->genericTag = tag;

   flush_v_HDDM(event,outputfp);
   close_v_HDDM(outputfp);
}
