/*
 *  findall :   an example utility that uses the hddsBrowser class
 *              to search the hdds (Hall D Detector Specification)
 *              geometry document for the appearance of a named volume
 *              and return the central coordinates and rotation angles
 *              of each placed copy in the geometry.
 *
 *  Original version - Richard Jones, June 3 2008.
 *
 *  Notes:
 *  ------
 * 1. As a by-product of using the DOM parser to access the xml source,
 *    findall verifies the source against the schema before translating it.
 *    Therefore it may also be used as a validator of the xml specification
 *    (see the -v option).
 */

#define APP_NAME "findall"

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

using namespace xercesc;

#include "XString.hpp"
#include "XParsers.hpp"
#include "hddsBrowser.hpp"

#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <list>
#include <map>

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

void usage()
{
    std::cerr
         << "Usage:    " << APP_NAME << " [-v] {HDDS file} {volume}"
         << std::endl <<  "Options:" << std::endl
         << "    -v   validate only" << std::endl;
}

int main(int argC, char* argV[])
{
   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      std::cerr
           << APP_NAME << " - error during initialization!"
           << std::endl << S(message) << std::endl;
      return 1;
   }

   if (argC < 3)
   {
      usage();
      return 1;
   }
   else if ((argC == 2) && (strcmp(argV[1], "-?") == 0))
   {
      usage();
      return 2;
   }

   XString xmlFile;
   XString targetVolume;
   bool dosearch = true;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
         break;

      if (strcmp(argV[argInd], "-v") == 0)
         dosearch = false;
      else
         std::cerr
              << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << std::endl;
   }

   if (argInd != argC - 2)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd++];
   targetVolume = argV[argInd++];


   if (dosearch)
   {
      hddsBrowser *browser = new hddsBrowser(xmlFile);
      std::vector<Refsys> *vlist = browser->find(targetVolume);
      for (std::vector<Refsys>::iterator it = vlist->begin();
           it < vlist->end();
           ++it) {
         std::vector<double> angles = it->getRotation();
         angles[0] *= 180/M_PI;
         angles[1] *= 180/M_PI;
         angles[2] *= 180/M_PI;
         std::cout << "found one at " << it->fOrigin[0] << ","
                   << it->fOrigin[1] << "," << it->fOrigin[2]
                   << " with rotation angles " << angles[0] << ","
                   << angles[1] << "," << angles[2]
                   << std::endl;
      }
      delete vlist;
   }

   XMLPlatformUtils::Terminate();
   return 0;
}
