/*
 *  hddmcat :	tool that reads in a sequence of HDDM files
 *		and catenates them into a single HDDM stream
 *
 *  Version 1.2 - Richard Jones, December 2005.
 *  - Updated code to use STL strings and vectors instead of old c-style
 *    pre-allocated arrays and strXXX functions.
 *  - Moved functions into classes grouped by function for better clarity.
 *
 *  Original version - Richard Jones, February 24 2004.
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <list>
using namespace std;

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <rpc/rpc.h>
#include <unistd.h>


void usage()
{
   cerr << "\nUsage:\n"
        << "    hddmcat file1.hddm [file2.hddm] ...\n\n"
        << endl;
}

int main(int argC, char* argV[])
{
   string xFilename;

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

   string hddmFile;
   istream* ifs;
   if (argInd == argC)
   {
      ifs = &cin;
   }
   else if (argInd < argC)
   {
      hddmFile = string(argV[argInd++]);
      ifs = new ifstream(hddmFile.c_str());
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifs->good())
   {
      cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
      exit(1);
   }

   list<std::string*> stringList;
   stringList.push_back(new string);
   list<std::string*>::iterator h;
   h = stringList.begin();
   if (std::getline(*ifs,**h))
   {
      if ((*h)->substr(0,5) == "<?xml")
      {
         cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if ((*h)->substr(0,5) == "<HDDM")
      {
         cout << **h << endl;
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
   stringList.push_back(new string);
   while (getline(*ifs,**(++h)))
   {
      cout << **h << endl;
      if (**h == "</HDDM>")
      {
         break;
      }
      stringList.push_back(new string);
   }

   const int bufferSize = 65536;
   char buffer[bufferSize];
   int count;
   while ((count = (ifs->read(buffer,bufferSize), ifs->gcount())))
   {
      cout.write(buffer,count);
   }
   if (ifs != &cin)
   {
      ((ifstream*)ifs)->close();
   }

   while (argInd < argC)
   {
      ifstream* ifs;
      hddmFile = argV[argInd++];
      ifs = new ifstream(hddmFile.c_str());
      if (!ifs->good())
      {
         cerr << "hddmcat: Error opening input stream " << hddmFile << endl;
         exit(1);
      }
      h = stringList.begin();
      string line;
      if (getline(*ifs,line))
      {
         if (line.substr(0,5) == "<?xml")
         {
            cerr << "hddmcat: Error reading input stream " << hddmFile << endl;
            cerr << "Input file appears to be an xml document!" << endl;
            exit(1);
         }
         else if (**h == line)
         {
	    ++h;
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
      while (getline(*ifs,line))
      {
         if (h == stringList.end() || **h != line)
         {
            cerr << "hddmcat: Input stream contains invalid hddm header"
                 << endl;
            exit(1);
	 }
	 else if (++h == stringList.end() && line == "</HDDM>")
         {
            break;
	 }
      }

      while ((count = (ifs->read(buffer,bufferSize), ifs->gcount())))
      {
         cout.write(buffer,count);
      }
      delete ifs;
   }
}
