/*
 *  hddmcat :	tool that reads in a sequence of HDDM files
 *		and catenates them into a single HDDM stream
 *
 *  Original version - Richard Jones, February 24 2004.
 *
 */

#include <iostream>
using namespace std;

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <rpc/xdr.h>
#include <string.h>
#include <unistd.h>

#include "fdstream.hpp"

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
   int ifd;
   boost::fdistream* ifs;
   if (argInd == argC)
   {
      hddmFile = new char[1];
      *hddmFile = 0;
      ifd = 0;
   }
   else if (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifd = open(hddmFile,O_RDONLY);
   }
   else
   {
      usage();
      return 1;
   }
   if (ifd < 0)
   {
      cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
      exit(1);
   }
   ifs = new boost::fdistream(ifd);

   struct stringArray* head = new struct stringArray;
   struct stringArray* h = head;
   size_t lineSize = 500;
   h->line = new char [lineSize];
   if (ifs->getline(h->line,lineSize))
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
   while (ifs->getline(h->line,lineSize))
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

   const int bufferSize = 1000000;
   char buffer[bufferSize];
   int count;
   while (ifs->read(buffer,bufferSize), count = ifs->gcount())
   {
      cout.write(buffer,count);
   }
   if (ifd)
      close(ifd);
   delete *ifs;

   while (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifd = open(hddmFile,O_RDONLY);
      if (ifd < 0)
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
      ifs = new boost::fdistream(ifd);
      h = head;
      char* line = new char [lineSize];
      if (ifs->getline(line,lineSize))
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
      while (ifs->getline(line,lineSize))
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

      while (ifs->read(buffer,bufferSize), count=ifs->gcount())
      {
         fwrite(buffer,sizeof(int),count,stdout);
      }
      close(ifd);
      delete *ifs;
   }
}
