/*
 *  hddmcat :	tool that reads in a sequence of HDDM files
 *		and catenates them into a single HDDM stream
 *
 *  Original version - Richard Jones, February 24 2004.
 *
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE true
#endif

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/xdr.h>
#include <string.h>
#include <unistd.h>

struct stringArray;
struct stringArray
{
   char *line;
   struct stringArray *next;
};

void usage()
{
   cerr << "\nUsage:\n"
        << "    hddmcat file1.hddm [file2.hddm] ...\n\n"
        << endl;
}

int main(int argC, char* argV[])
{
   char *xFilename = 0;

   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
      {
         break;
      }
      else
      {
         usage();
         return 1;
      }
   }

   char* hddmFile;
   FILE* ifd;
   if (argInd == argC)
   {
      ifd = stdin;
      hddmFile = new char[1];
      *hddmFile = 0;
   }
   else if (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifd = fopen(hddmFile,"r");
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifd)
   {
      cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
      exit(1);
   }

   struct stringArray* head = new struct stringArray;
   struct stringArray* h = head;
   size_t lineSize = 500;
   h->line = new char [lineSize];
   if (getline(&h->line,&lineSize,ifd))
   {
      if (strstr(h->line,"<?xml") != 0)
      {
         cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if (strstr(h->line,"<HDDM") == h->line)
      {
         printf("%s",h->line);
      }
      else
      {
         cerr << "hddmcat: Input stream contains invalid hddm header"
              << endl;
         exit(1);
      }
   }
   else
   {
      cerr << "hddmcat: Error reading from input stream " << hddmFile << endl;
      exit(1);
   }
   h = h->next = new struct stringArray;
   h->line = new char [lineSize];
   while (getline(&h->line,&lineSize,ifd))
   {
      printf("%s",h->line);
      if (strstr(h->line,"</HDDM>") != 0)
      {
	 h->next = 0;
         break;
      }
      h = h->next = new struct stringArray;
      h->line = new char [lineSize];
   }

   int bufferSize = 1000000;
   int buffer[bufferSize];
   int count;
   while (count = fread(buffer,sizeof(int),bufferSize,ifd))
   {
      fwrite(buffer,sizeof(int),count,stdout);
   }

   while (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifd = fopen(hddmFile,"r");
      if (!ifd)
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
      h = head;
      char* line = new char [lineSize];
      if (getline(&line,&lineSize,ifd))
      {
         if (strstr(line,"<?xml") != 0)
         {
            cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
            cerr << "Input file appears to be an xml document!" << endl;
            exit(1);
         }
         else if (strstr(line,h->line) == line)
         {
	    h = h->next;
         }
         else
         {
            cerr << "hddmcat: Input stream contains invalid hddm header"
                 << endl;
            exit(1);
         }
      }
      else
      {
         cerr << "hddmcat: Error reading from input stream " << hddmFile << endl;
         exit(1);
      }
      while (getline(&line,&lineSize,ifd))
      {
         if (h == 0 || strstr(line,h->line) != line)
         {
            cerr << "hddmcat: Input stream contains invalid hddm header"
                 << endl;
            exit(1);
	 }
	 else if (h->next == 0 && strstr(line,"</HDDM>") == line)
         {
            break;
	 }
         h = h->next;
      }

      while (count = fread(buffer,sizeof(int),bufferSize,ifd))
      {
         fwrite(buffer,sizeof(int),count,stdout);
      }
   }
}
