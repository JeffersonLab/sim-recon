/*
 *  hddmcat :	tool that reads in a sequence of HDDM files
 *		and catenates them into a single HDDM stream
 *
 *  Original version - Richard Jones, February 24 2004.
 *
 */

#include <iostream>
using namespace std;
#include <fstream>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
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
   istream* ifs;
   if (argInd == argC)
   {
      hddmFile = new char[1];
      *hddmFile = 0;
      ifs = &cin;
   }
   else if (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifs = new ifstream(hddmFile);
      if (! ((ifstream*)ifs)->is_open())
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
   }
   else
   {
      usage();
      return 1;
   }

   struct stringArray* head = new struct stringArray;
   struct stringArray* h = head;
   h->line = new char[1000];
   if (ifs->getline(h->line,999))
   {
      if (strstr(h->line,"<?xml") != 0)
      {
         cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if (strstr(h->line,"<HDDM") == h->line)
      {
         cout << h->line << endl;
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
   h->line = new char[1000];
   while (ifs->getline(h->line,999))
   {
      cout << h->line << endl;
      if (strstr(h->line,"</HDDM>") != 0)
      {
	 h->next = 0;
         break;
      }
      h = h->next = new struct stringArray;
      h->line = new char[1000];
   }

   const int bufferSize = 1000000;
   char buffer[bufferSize];
   int count;
   while (ifs->read(buffer,bufferSize), count = ifs->gcount())
   {
      cout.write(buffer,count);
   }
   if (ifs != &cin)
   {
      delete ifs;
   }

   while (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifs = new ifstream(hddmFile);
      if (! ((ifstream*)ifs)->is_open())
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
      h = head;
      char* line = new char[1000];
      if (ifs->getline(line,999))
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
      while (ifs->getline(line,999))
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
         cout.write(buffer,count);
      }
      delete ifs;
   }
}
