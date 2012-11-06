#include <iostream>
#include <fstream>
using namespace std;
#include "hddm_v.h"
#include "hddm_v.hpp"

void test_hddm_c_interface(char *outfile);
void test_hddm_cpp_interface(char *outfile);

int main()
{
   // the following two output files should be identical
   test_hddm_c_interface("test.hddm");
   test_hddm_cpp_interface("test2.hddm");
   return 0;
}

void test_hddm_c_interface(char *outfile)
{
   v_iostream_t *outputfp = init_v_HDDM(outfile);
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

void test_hddm_cpp_interface(char *outfile)
{
   ofstream ofs(outfile);
   if (!ofs.is_open()) {
      cerr << "Test failed: error opening output file\n";
      exit(1);
   }

   hddm_v::HDDM record;
   hddm_v::GenericTagList generics = record.addGenericTags();
   hddm_v::FloatTagList floats = generics().addFloatTags();
   hddm_v::DoubleTagList doubles = generics().addDoubleTags();
   hddm_v::ParticleTagList particles = generics().addParticleTags();
   hddm_v::StringTagList strings = generics().addStringTags();
   hddm_v::IntTagList ints = generics().addIntTags();
   hddm_v::LongTagList longs = generics().addLongTags();
   hddm_v::BooleanTagList booleans = generics().addBooleanTags();
   hddm_v::AnyURITagList anyURIs = generics().addAnyURITags();
   floats().setPi(3.141593);
   doubles().setPi(3.141592592636);
   particles().setPi(PiPlus);
   strings().setQuote("pass the red quarks, please");
   ints().setMagic(133557799);
   longs().setMagic(133557799002244668LL);
   booleans().setTruth(1);
   anyURIs().setUri("http://portal.gluex.org");

   hddm_v::ostream ostr(ofs);
   ostr << record;
}
