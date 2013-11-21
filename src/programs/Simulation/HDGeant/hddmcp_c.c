/*
 * hddmcp - a utility program for copying from one hddm file to another,
 *          and in the process stripping out everything except a desired
 *          subset of the data fields.  The template for the output hddm
 *          file (a subset of what is found in the input hddm file) is
 *          assumed to have already been created, perhaps by grabbing the
 *          header from the input hddm file, stripping out the desired parts
 *          and dumping the resulting template into a new file called x.hddm,
 *          which then must be processed into c using the following command.
 *
 *            $ hddm-c x.hddm
 *
 * Richard Jones
 * GlueX collaboration
 * October 6, 2009
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <HDDM/hddm_s.h>

int main(int argc, char **argv)
{
   s_HDDM_t *thisInputEvent = 0;
   s_iostream_t *thisInputFile = 0;
   s_iostream_t *thisOutputFile = 0;
   int input, output;
   if (argc < 2) {
      printf("Usage: hddmcp <inputfile1.hddm> ... <outputfile.hddm>\n");
      exit(1);
   }
   else {
      output = argc-1;
   }
   thisOutputFile = init_s_HDDM(argv[output]);
   for (input=1; input<output; input++) {
      if (! (thisInputFile = open_s_HDDM(argv[input]))) {
         fprintf(stderr,"Error - could not open input file %s\n",
                 argv[input]);
         exit(1);
      }
      int count = 0;
      while (thisInputEvent = read_s_HDDM(thisInputFile)) {
         flush_s_HDDM(thisInputEvent,thisOutputFile);
      }
      close_s_HDDM(thisInputFile);
   }
}
