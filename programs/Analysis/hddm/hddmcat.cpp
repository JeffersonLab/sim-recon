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

int getline(char* buf, int lenbuf, FILE* stream)
{
   int count = 0;
   for (count = 0; count < lenbuf-1;)
   {
      int c = fgetc(stream);
      if (c == EOF) break;
      buf[count++] = c;
      if (c == '\n') break;
   }
   buf[count] = 0;
   return count;
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
   FILE* ifp;
   if (argInd == argC)
   {
      hddmFile = new char[1];
      *hddmFile = 0;
      ifp = stdin;
   }
   else if (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifp = fopen(hddmFile,"r");
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifp)
   {
      cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
      exit(1);
   }

   struct stringArray* head = new struct stringArray;
   struct stringArray* h = head;
   h->line = new char[1000];
   if (getline(h->line,1000,ifp))
   {
      if (strstr(h->line,"<?xml") != 0)
      {
         cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if (strstr(h->line,"<HDDM") == h->line)
      {
         cout << h->line;
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
   while (getline(h->line,1000,ifp))
   {
      cout << h->line;
      if (strstr(h->line,"</HDDM>") != 0)
      {
	 h->next = 0;
         break;
      }
      h = h->next = new struct stringArray;
      h->line = new char[1000];
   }

   const int bufferSize = 65536;
   char buffer[bufferSize];
   int count;
   while (count = fread(buffer,sizeof(char),bufferSize,ifp))
   {
      cout.write(buffer,count);
   }
   if (ifp != stdin)
   {
      fclose(ifp);
   }

   while (argInd < argC)
   {
      hddmFile = argV[argInd++];
      ifp = fopen(hddmFile,"r");
      if (!ifp)
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
      h = head;
      char* line = new char[1000];
      if (getline(line,1000,ifp))
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
         cerr << "hddmcat: Error reading from input stream " << hddmFile
              << endl;
         exit(1);
      }
      while (getline(line,1000,ifp))
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

      while (count = fread(buffer,sizeof(char),bufferSize,ifp))
      {
         cout.write(buffer,count);
      }
      fclose(ifp);
   }
}
