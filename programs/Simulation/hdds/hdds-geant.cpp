/*
 *  hdds-geant :   an interface utility that reads in a HDDS document
 *		   (Hall D Detector Specification) and writes out a
 *		   GEANT-3 geometry description in the form of a
 *		   fortran subroutine.
 *
 *  Original version - Richard Jones, May 19 2001.
 *
 *  Notes:
 *  ------
 * 1. The HDDS specification is an xml document, as described by HDDS.dtd.
 * 2. Access by hdds-geant to the xml source is through the industry-
 *    standard DOM-1 interface.
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 * 4. Output is sent to standard out through the ordinary c++ i/o library.
 * 5. As a by-product of using the DOM parser to access the xml source,
 *    hdds-geant verifies the source against the dtd before translating it.
 *    Therefore it may also be used as a validator of the xml specification
 *    (see the -v option).
 *
 *  Implementation:
 *  ---------------
 * Most of the translation was straight-forward, but there were a couple
 * of tricky cases where decisions had to be made.  I think that these
 * choices should work out in most cases.  If not, further tuning of the
 * algorithm will be necessary.  Here are the tricky cases.
 *
 * 1. When to use divisions instead of placing multiple copies.
 * 
 *  Most of the time when a <mpos...> command appears it can be translated
 *  into a division of the mother volume in Geant.  This is good to do
 *  because it makes both more compact description and is more efficient
 *  at tracking time.  The difficulty here is that there is no easy way
 *  to check if the contents fit entirely inside the division.  This is
 *  not a problem in the HDDS geometry description because the <mpos...>
 *  command is only for positioning, and makes no statement about what
 *  slice of the mother it occupies.  But it is a problem for Geant because
 *  contents of divisions have to fit inside the division.  This can
 *  happen at any time, but it most frequently happens when the object
 *  is being rotated before placement, as in the case of stereo layers
 *  in a drift chamber.  My solution is to make a strict set of rules
 *  that are required before hdds-geant will create a division in response
 *  to a <mpos...> command, and do individual placement by default.
 *  The rules for creation of divisions are as follows:
 *      (a) the <composition> command must have a solid container, either
 *          via the envelope="..." attribute or by itself being a division.
 *      (b) the <mpos...> command must be alone inside its <composition>
 *      (c) the kind of shape of the container must be compatible with the
 *          <mpos...> type, eg. <mposPhi> would work if its container is
 *          a "tubs" (or division theroef) but not if it is a "box".
 *      (d) for <mposPhi> the impliedRot attribute must be "true".
 *      (e) the rot="..." attribute must be zeros or missing.
 *  The last condition is not logically necessary for it to work, but it
 *  avoids the problems that often occur with rotated placements failing
 *  to fit inside the division.  To bypass this limitation and place
 *  rotated volumes inside divisions, one can simply create a new volume
 *  as a <composition> into which the content is placed with rotation,
 *  and then place the new volume with the <mpos...> command without rot.
 *
 * 2. How to recognize which media contain magnetic fields.
 *
 *  There is no provision in the hdds geometry model for magnetic field
 *  information.  Ultimately that is something that will be stored as a map
 *  somewhere in a database.  Geant needs to distinguish between 4 cases:
 *       (0) no magnetic field
 *       (1) general case of inhomogenous field (Runge-Kutta)
 *       (2) quasi-homogenous field with map (helical segments)
 *       (3) uniform field (helices along local z-axis)
 *  My solution is as follows.  By default, I assume case 0.  For all
 *  contents of a composition named "fieldVolume" I assign case 2.  For
 *  all contents of a composition named "sweepMagnet" I assign case 3.
 *  For the magnitude of the field, I simply store a constant field value
 *  (kG) for each case, and rely on the GUFLD user routine in Geant3 to
 *  handle the actual field values at tracking time.
 *
 * 3. What to do about stackX/stackY/stackZ tags
 *
 *  In the case of Boolean tags (union/intersection/subtraction) the choice
 *  was easy: no support in Geant3 for these cases.  For the stacks it is
 *  possible to construct such structures in Geant3 but it is complicated
 *  by the fact that the stacks do not give information about the kind or
 *  size of volume that should be used for the mother.  Since stacks can
 *  be implemented easily using compositions anyway, I decided not to include
 *  support for them in hdds-geant.
 */

#include <util/PlatformUtils.hpp>
#include <sax/SAXException.hpp>
#include <sax/SAXParseException.hpp>
#include <parsers/DOMParser.hpp>
#include <dom/DOM_DOMException.hpp>
#include <dom/DOM_NamedNodeMap.hpp>
#include <dom/DOMString.hpp>

#include "hdds-geant.hpp"

#include <assert.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

double fieldStrength[] =
{
   0.0,	  // zero field regions
   0.0,   // inhomogenous field regions (unused)
   22.4,  // mapped field regions: solenoid (kG, approximate)
   2.0    // uniform field regions: sweep magnets (kG)
};

void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-geant [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

class Refsys
{
 public:
   DOM_Element fMother;		// current mother volume element
   int fMagField;		// flag indicating magnetic field
   double fOrigin[3];		// x,y,z coordinate of volume origin (cm)
   double fRmatrix[3][3];	// rotation matrix (daughter -> mother)
   int fRotation;		// unique Rmatrix flag

   static char* fIdentifierList; // list of identifier strings (space-sep)
   static int fVolumes;		// total number of volumes to far

   struct VolIdent
   {
      DOMString fieldS;	
      int value;
      int step;
   };
   struct VolIdent fIdentifier[999];	// identifier tag list 
   int fIdentifiers;			// length of identifier tag list

   Refsys();				// empty constructor
   Refsys(const Refsys& src);		// copy constructor
   Refsys& operator=(Refsys& src);	// copy operator
   Refsys& reset();			// reset origin, Rmatrix
   Refsys& reset(const Refsys& ref);	// reset origin, Rmatrix to ref
   Refsys& shift(const double vector[3]); // translate origin
   Refsys& shift(const Refsys& ref);	  // copy origin from ref
   Refsys& shift(const Refsys& ref,
                 const double vector[3]); // translate origin in ref frame
   Refsys& rotate(const double omega[3]); // rotate by vector omega (radians)
   Refsys& rotate(const Refsys& ref);	  // copy Rmatrix from ref
   Refsys& rotate(const Refsys& ref,
                  const double omega[3]); // rotate by omega in ref frame

   int createMaterial(DOM_Element& el);	// generate code for materials
   int createSolid(DOM_Element& el);	// generate code for solids
   int createRotation();		// generate code for rotations
   int createDivision(char* divStr,
                      int ncopy,
                      int iaxis,
                      double start,
                      double step);	// generate code for divisions
   int createVolume(DOM_Element& el);	// generate code for placement

private:
   static int fRotations;	// non-trivial rotations defined so far
};

void fortranHeader()
{
   cout << "*" 							<< endl
        << "* HDDSgeant3 - fortran geometry definition package" << endl
        << "*              for the Hall D experiment."		<< endl
        << "*"							<< endl
        << "*         WARNING: DO NOT EDIT THIS FILE"		<< endl
        << "*"							<< endl
        << "* This file was generated automatically from the"	<< endl
        << "* HDDS xml geometry definition by the hdds-geant"	<< endl
        << "* translator.  Any changes made to this file will"	<< endl
        << "* disappear as soon as it is regenerated from the"	<< endl
        << "* xml source.  To introduce Geant3 optimizations,"	<< endl
        << "* see the subroutine Goptimize() in goptimize.F."	<< endl
        << "*"							<< endl
        << "      subroutine HDDSgeant3" 			<< endl
        << "      implicit none"				<< endl
        << "      integer imate"				<< endl
        << "      character*20 chnama,namate"			<< endl
        << "      real a,z,dens,radl,absl,ubuf(99)"		<< endl
        << "      integer nwbuf"				<< endl
        << "      real amat(99),zmat(99),wmat(99)"		<< endl
        << "      integer nlmat"				<< endl
        << "      integer itmed"				<< endl
        << "      character*20 natmed"				<< endl
        << "      integer nmat,isvol,ifield"			<< endl
        << "      real fieldm,tmaxfd,stemax,deemax,epsil,stmin" << endl
        << "      character*4 chname,chshap,chmoth"		<< endl
        << "      integer nmed,npar,ivolu"			<< endl
        << "      real par(99)"					<< endl
        << "      integer irot"					<< endl
        << "      real theta1,phi1,theta2,phi2,theta3,phi3"	<< endl
        << "      integer nr,ndiv,iaxis,numed,ndvmax"		<< endl
        << "      real step,c0"					<< endl
        << "      real x,y"					<< endl
        << "      character*4 chonly"				<< endl;
}

void fortranTrailer()
{
   cout << "      end"						<< endl;
}

void fortranGetfunc(DOM_Element& el, char* ident)
{
   int* start = new int[Refsys::fVolumes + 1];
   int* table = new int[999999];
   int tableLength = 0;

   char funcName[40];
   sprintf(funcName, "get%s", ident);
   funcName[3] = toupper(funcName[3]);

   DOM_NodeList alltagsList = el.getOwnerDocument().getElementsByTagName("*");
   for (int itag = 0; itag < alltagsList.getLength(); itag++)
   {
      DOM_Node node = alltagsList.item(itag);
      DOM_Element elem = (DOM_Element&) node;
      DOMString icopyS = elem.getAttribute("Geant3icopy");
      DOMString ivoluS = elem.getAttribute("Geant3ivolu");
      if (ivoluS != 0)
      {
         char* ivoluStr = ivoluS.transcode();
         char* icopyStr = icopyS.transcode();
         int ivolu = atoi(ivoluStr);
         int icopy = atoi(icopyStr);
         delete [] icopyStr;
         delete [] ivoluStr;
         DOMString idlistS = elem.getAttribute(ident);
         if (idlistS != 0)
         {
            char* idlistStr = idlistS.transcode();
            start[ivolu] = tableLength + 1;
            char* idStr;
            char** token = &idlistStr;
            for (idStr = strsep(token, " ");
                 (idStr != 0) && (strlen(idStr) > 0);
                 idStr = strsep(token, " "))
            {
               table[tableLength++] = atoi(idStr);
               --icopy;
            }
            for (; icopy > 0; --icopy)
            {
               table[tableLength++] = 0;
            }
            delete [] idlistStr;
         }
         else
         {
            start[ivolu] = 0;
         }
      }
   }

   cout << endl
        << "      function " << funcName << "()" << endl
        << "      implicit none" << endl
        << "      integer " << funcName << endl
        << "      integer nlevel,names,number,lvolum" << endl
        << "      common /gcvolu/nlevel,names(15),number(15),lvolum(15)"
        << endl;

   if (tableLength > 0)
   {
      cout  << "      integer istart(" << Refsys::fVolumes << ")" << endl
            << "      data istart/" << endl;

      for (int i = 1; i <= Refsys::fVolumes;)
      {
         int i0 = i;
         char str[16];
         sprintf(str, "%5d", start[i]);
         cout << "     + " << str;
         for (i++; i <= Refsys::fVolumes; i++)
         {
            if (i == (i0 + 10)) break;
            sprintf(str, ",%5d", start[i]);
            cout << str;
         }
         if (i > Refsys::fVolumes)
         {
            cout << "/" << endl;
         }
         else
         {
            cout << "," << endl;
         }
      }
      cout << "      integer lookup(" << tableLength << ")" << endl
           << "      data lookup/" << endl;
      for (int i = 0; i < tableLength;)
      {
         int i0 = i;
         char str[16];
         sprintf(str, "%5d", table[i]);
         cout << "     + " << str;
         for (i++; i < tableLength; i++)
         {
            if (i == (i0 + 10)) break;
            sprintf(str, ",%5d", table[i]);
            cout << str;
         }
         if (i == tableLength)
         {
            cout << "/" << endl;
         }
         else
         {
            cout << "," << endl;
         }
      }
      cout << "      integer level,index" << endl
           << "      integer " << ident << endl
           << "      " << funcName << " = 0" << endl
           << "      do level=1,nlevel" << endl
           << "        index = istart(lvolum(level))" << endl
           << "        if (index.gt.0) then" << endl
           << "          " << ident
           << " = lookup(index + number(level) - 1)" << endl
           << "          if (" << ident << ".gt.0) then" << endl
           << "            " << funcName << " = " << ident << endl
           << "          endif" << endl
           << "        endif" << endl
           << "      enddo" << endl;
   }
   else
   {
      cout << "      " << funcName << " = 0" << endl;
   }
   cout << "      end" << endl;
   delete [] table;
   delete [] start;
}

int main(int argC, char* argV[])
{
   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      cerr << "hdds-geant: Error during initialization! :\n"
           << StrX(toCatch.getMessage()) << endl;
      return 1;
   }

   if (argC < 2)
   {
      usage();
      return 1;
   }
   else if ((argC == 2) && (strcmp(argV[1], "-?") == 0))
   {
      usage();
      return 2;
   }

   const char*  xmlFile = 0;
   bool geantOutput = true;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
         break;

      if (strcmp(argV[argInd], "-v") == 0)
         geantOutput = false;
      else
         cerr << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << endl;
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd];

   DOMParser parser;
   parser.setValidationScheme(DOMParser::Val_Auto);
   parser.setCreateEntityReferenceNodes(false);
   parser.setDoNamespaces(false);

   MyOwnErrorHandler errorHandler;
   parser.setErrorHandler(&errorHandler);

   try
   {
      parser.parse(xmlFile);
   }
   catch (const XMLException& toCatch)
   {
      cerr << "\nhdds-geant: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << StrX(toCatch.getMessage()) << "\n" << endl;
      return -1;
   }
   catch (const DOM_DOMException& toCatch)
   {
      cerr << "\nhdds-geant: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << toCatch.msg.transcode() << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 4;
   }
   catch (...)
   {
      cerr << "\nhdds-geant: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      XMLPlatformUtils::Terminate();
      return 4;
   }

   if (errorHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
   }

   DOM_Document document = parser.getDocument();
   DOM_Element rootEl = document.getElementById("everything");
   if (rootEl == 0)
   {
      cerr << "hdds-geant : Error scanning HDDS document, " << endl
           << "  no element named \"eveything\" found" << endl;
      return 1;
   }


   if (geantOutput)
   {
      Refsys mrs;
      fortranHeader();
      mrs.createVolume(rootEl);
      fortranTrailer();

      if (mrs.fIdentifierList != 0)
      {
         char* identStr;
         char** token = &mrs.fIdentifierList;
         for (identStr = strsep(token, " ");
              identStr != 0;
              identStr = strsep(token, " "))
         {
            fortranGetfunc(rootEl, identStr);
         }
      }
   }

   XMLPlatformUtils::Terminate();
   return 0;
}


MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}

// Overrides of the SAX ErrorHandler interface

void MyOwnErrorHandler::error(const SAXParseException& e)
{
   fSawErrors = true;
   cerr << "\nhdds-geant: Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   cerr << "\nhdds-geant: Fatal Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   cerr << "\nhdds-geant: Warning at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::resetErrors()
{
}

int Refsys::fVolumes = 0;
int Refsys::fRotations = 0;
char* Refsys::fIdentifierList = 0;

Refsys::Refsys()			// empty constructor
{
   fMother = 0;
   fMagField = 0;
   this->reset();
}

Refsys::Refsys(const Refsys& src)	// copy constructor
{
   reset(src);
   fMother = src.fMother;
   fRotation = src.fRotation;
   fMagField = src.fMagField;
   fIdentifiers = src.fIdentifiers;
   for (int i = 0; i < fIdentifiers; i++)
   {
      fIdentifier[i] = src.fIdentifier[i];
   }
} 

Refsys& Refsys::operator=(Refsys& src)	// copy operator (deep sematics)
{
   Refsys* dst = new Refsys(src);
   return *dst;
}

Refsys& Refsys::reset()			// reset Origin, Rmatrix to null
{
   fOrigin[0] = fOrigin[1] = fOrigin[2] = 0;
   fRmatrix[0][0] = fRmatrix[1][1] = fRmatrix[2][2] = 1;
   fRmatrix[0][1] = fRmatrix[1][0] = fRmatrix[1][2] =
   fRmatrix[0][2] = fRmatrix[2][0] = fRmatrix[2][1] = 0;
   fRotation = 0;
   return *this;
}

Refsys& Refsys::reset(const Refsys& ref) // reset Origin, Rmatrix to ref
{
   for (int i = 0; i < 3; i++)
   {
      fOrigin[i] = ref.fOrigin[i];
      fRmatrix[i][0] = ref.fRmatrix[i][0];
      fRmatrix[i][1] = ref.fRmatrix[i][1];
      fRmatrix[i][2] = ref.fRmatrix[i][2];
   }
   fRotation = ref.fRotation;
   return *this;
}

Refsys& Refsys::shift(const double vector[3])     // translate origin
{
   for (int i = 0; i < 3; i++)
   {
      fOrigin[i] += fRmatrix[i][0] * vector[0] +
                    fRmatrix[i][1] * vector[1] +
                    fRmatrix[i][2] * vector[2];
   }
   return *this;
}

Refsys& Refsys::shift(const Refsys& ref)      // copy origin from ref
{
   fOrigin[0] = ref.fOrigin[0];
   fOrigin[1] = ref.fOrigin[1];
   fOrigin[2] = ref.fOrigin[2];
   return *this;
}

Refsys& Refsys::shift(const Refsys& ref,
                      const double vector[3]) // translate origin in ref frame
{
   Refsys myRef(ref);
   myRef.shift(vector);
   return shift(myRef);
}

Refsys& Refsys::rotate(const double omega[3]) // rotate by vector omega (rad)
{
   if ( (omega[0] != 0.0) || (omega[1] != 0.0) || (omega[2] != 0.0) )
   {
      double cosx = cos(omega[0]);
      double sinx = sin(omega[0]);
      double cosy = cos(omega[1]);
      double siny = sin(omega[1]);
      double cosz = cos(omega[2]);
      double sinz = sin(omega[2]);

      for (int i = 0; i < 3; i++)
      {
         double x[3];
         double xx[3];
         x[0] = fRmatrix[i][0] * cosz + fRmatrix[i][1] * sinz;
         x[1] = fRmatrix[i][1] * cosz - fRmatrix[i][0] * sinz;
         x[2] = fRmatrix[i][2];
         xx[0] = x[0] * cosy - x[2] * siny;
         xx[1] = x[1];
         xx[2] = x[2] * cosy + x[0] * siny;
         fRmatrix[i][0] = xx[0];
         fRmatrix[i][1] = xx[1] * cosx + xx[2] * sinx;
         fRmatrix[i][2] = xx[2] * cosx - xx[1] * sinx;
      }

      fRotation = -1;
   }

   return *this;
}

Refsys& Refsys::rotate(const Refsys& ref)      // copy Rmatrix from ref
{
   for (int i = 0; i < 3; i++)
   {
      fRmatrix[i][0] = ref.fRmatrix[i][0];
      fRmatrix[i][1] = ref.fRmatrix[i][1];
      fRmatrix[i][2] = ref.fRmatrix[i][2];
   }
   fRotation = ref.fRotation;
   return *this;
}

Refsys& Refsys::rotate(const Refsys& ref,
                       const double omega[3])  // rotate by omega in ref frame
{
   Refsys myRef(ref);
   myRef.rotate(omega);
   return rotate(myRef);
}

int Refsys::createMaterial(DOM_Element& el)
{
   static int imateCount = 0;
   int imate = ++imateCount;
   char imateStr[30];
   sprintf(imateStr, "%d", imate);
   el.setAttribute("Geant3imate", imateStr);

   DOMString tagS = el.getTagName();
   DOMString matS = el.getAttribute("name");
   char* matStr = matS.transcode();
   char* paramStr;
   char* valueStr;

   paramStr = el.getAttribute("a").transcode();
   double a = atof(paramStr);
   delete [] paramStr;

   paramStr = el.getAttribute("z").transcode();
   double z = atof(paramStr);
   delete [] paramStr;

   double dens = -1;
   double radl = -1;
   double absl = -1;
   DOM_NodeList paramList = el.getElementsByTagName("real");
   for (int ip = 0; ip < paramList.getLength(); ip++)
   {
      DOM_Node node = paramList.item(ip);
      DOM_Element elem = (DOM_Element&) node;
      DOMString paramS = elem.getAttribute("name");
      DOMString valueS = elem.getAttribute("value");
      char* valueStr = valueS.transcode();
      if (paramS.equals("density"))
      {
         dens = atof(valueStr);
      }
      else if (paramS.equals("radlen"))
      {
         radl = atof(valueStr);
      }
      else if (paramS.equals("abslen"))
      {
         absl = atof(valueStr);
      }
      delete [] valueStr;
   }

   int nList = 0;
   int iList[999];
   double wList[999];
   DOM_Document document = el.getOwnerDocument();
   DOM_NodeList compList = el.getElementsByTagName("addmaterial");
   for (int ic = 0; ic < compList.getLength(); ic++)
   {
      DOM_Node node = compList.item(ic);
      DOM_Element elem = (DOM_Element&) node;
      DOMString compS = elem.getAttribute("material");
      DOM_Element compEl = document.getElementById(compS);
      DOMString cS;
      if ( (cS = compEl.getAttribute("Geant3imate")) != 0)
      {
         char* cStr = cS.transcode();
         iList[nList] = atoi(cStr);
         delete [] cStr;
      }
      else
      {
         iList[nList] = createMaterial(compEl);
      }

      DOM_NodeList atomList = elem.getElementsByTagName("natoms");
      DOM_NodeList fracList = elem.getElementsByTagName("fractionmass");
      if (atomList.getLength() == 1)
      {
         DOMString typeS = compEl.getTagName();
         if (! typeS.equals("element"))
         {
            cerr << "hdds-geant: error processing composite " << matStr << endl
                 << "natoms can only be specified for elements." << endl;
            exit(1);
         }
         else if (dens < 0)
         {
            char* tagStr = tagS.transcode();
            cerr << "hdds-geant error: " << tagStr << " " << matStr
                 << " is missing a density specification." << endl;
            exit(1);
         }
         DOM_Node node = atomList.item(0);
         DOM_Element elem = (DOM_Element&) node;
         char* nStr = elem.getAttribute("n").transcode();
         wList[nList] = atoi(nStr);
         delete [] nStr;
      }
      else if (fracList.getLength() == 1)
      {
         DOM_Node node = fracList.item(0);
         DOM_Element elem = (DOM_Element&) node;
         char* fStr = elem.getAttribute("fraction").transcode();
         wList[nList] = atof(fStr);
         delete [] fStr;

         double rho = 0;
         DOM_NodeList propList = compEl.getElementsByTagName("real");
         for (int ip = 0; ip < propList.getLength(); ip++)
         {
            DOM_Node pnode = propList.item(ip);
            DOM_Element pelem = (DOM_Element&) pnode;
            DOMString pS = pelem.getAttribute("name");
            DOMString vS = pelem.getAttribute("value");
            if (pS.equals("density"))
            {
               char* vStr = vS.transcode();
               rho = atof(vStr);
               delete [] vStr;
            }
         }
         assert (rho > 0);
         if (dens < 0)
         {
            dens -= wList[nList] / rho;
         }
      }
      else
      {
         char* tagStr = tagS.transcode();
         cerr << "hdds-geant error: " << tagStr << " " << matStr
              << " is missing some proportion data." << endl;
         exit(1);
      }
      nList++;
   }

   if (tagS.equals("element"))
   {
      if (dens < 0)
      {
         char* tagStr = tagS.transcode();
         cerr << "hdds-geant error: " << tagStr << " " << matStr
              << " is missing a density specification." << endl;
         exit(1);
      }
      if ((radl == -1) || (absl == -1))
      {
         cout << endl
              << "      imate = " << imate << endl
              << "      namate = \'" << matStr << "\'" << endl
              << "      a = " << a << endl
              << "      z = " << z << endl
              << "      dens = " << dens << endl
              << "      nlmat = 1" << endl
              << "      wmat(1) = 1" << endl
              << "      call gsmixt(imate,namate,a,z,dens,nlmat,wmat)"
              << endl;
      }
      else
      {
         cout << endl
              << "      imate = " << imate << endl
              << "      chnama = \'" << matStr << "\'" << endl
              << "      a = " << a << endl
              << "      z = " << z << endl
              << "      dens = " << dens << endl
              << "      radl = " << radl / (dens + 1e-30) << endl
              << "      absl = " << absl / (dens + 1e-30) << endl
              << "      nwbuf = 0" << endl
              << "      call gsmate(imate,chnama,a,z,dens,radl,absl,ubuf,nwbuf)"
              << endl;
      }
   }
   else
   {
      double rho = dens;
      if (rho < 0)
      {
         rho = -1/(dens + 0.999999);
      }
      cout << endl
           << "      imate = " << imate << endl
           << "      namate = \'" << matStr << "\'" << endl;
      for (int im = 0; im < nList; im++)
      {
         cout << "      wmat(" << im + 1 << ") = " << wList[im] << endl
              << "      call gfmate(" << iList[im] << ",chnama,"
              << "amat(" << im + 1 << "),zmat(" << im + 1 << "),"
              << "dens,radl,absl,ubuf,nwbuf)" << endl;
      }
      cout << "      dens = " << rho << endl
           << "      nlmat = " << ((wList[0] < 1) ? nList : -nList) << endl
           << "      call gsmixt(imate,namate,amat,zmat,dens,nlmat,wmat)"
           << endl;
   }

   
   if (dens < 0)
   {
      char densStr[30];
      sprintf(densStr, "%f", -1/(dens + 0.999999));
      DOM_Element densEl = document.createElement("real");
      densEl.setAttribute("name", "density");
      densEl.setAttribute("value", densStr);
      el.appendChild(densEl);
   }

   delete [] matStr;
   return imate;
}

void getConversions(DOM_Element& el, double& tocm, double& todeg)
{
   DOMString unitS;
   char* unitStr;

   unitS = el.getAttribute("unit_length");
   if (unitS.equals("mm"))
   {
      tocm = 0.1;
   }
   else if (unitS.equals("cm"))
   {
      tocm = 1.0;
   }
   else if (unitS.equals("m"))
   {
      tocm = 100.;
   }
   else
   {
      unitStr = unitS.transcode();
      char* tagStr = el.getTagName().transcode();
      cerr << "hdds-geant error: unknown length unit " << unitStr
           << " on tag " << tagStr << endl;
      exit(1);
   }

   unitS = el.getAttribute("unit_angle");
   if (unitS.equals("deg"))
   {
      todeg = 1.0;
   }
   else if (unitS.equals("mrad"))
   {
      todeg = 0.180/M_PI;
   }
   else
   {
      unitStr = unitS.transcode();
      char* tagStr = el.getTagName().transcode();
      cerr << "hdds-geant error: unknown angle unit " << unitStr
           << " on volume " << tagStr << endl;
      exit(1);
   }
}

int Refsys::createSolid(DOM_Element& el)
{
   DOMString shapeS = el.getTagName();
   DOMString nameS = el.getAttribute("name");
   DOMString materialS = el.getAttribute("material");
   char* nameStr = nameS.transcode();
   char* matStr = materialS.transcode();

   int imate;
   DOM_Document document = el.getOwnerDocument();
   DOM_Element matEl = document.getElementById(materialS);
   DOMString imateS = matEl.getAttribute("Geant3imate");
   if (imateS != 0)
   {
      char* imateStr = imateS.transcode();
      imate = atoi(imateStr);
      delete [] imateStr;
   }
   else
   {
      imate = createMaterial(matEl);
   }
   
   static int itmedCount = 0;
   int itmed = ++itmedCount;
   DOMString sensiS = el.getAttribute("sensitive");
   cout << endl
        << "      itmed = " << itmed << endl
        << "      natmed = \'" << nameStr << " " << matStr << "\'" << endl
        << "      nmat = " << imate << endl
        << "      isvol = " << (sensiS.equals("true") ? 1 : 0) << endl
        << "      ifield = " << fMagField << endl
        << "      fieldm = " << fieldStrength[fMagField] << endl
        << "      tmaxfd = " << ((fMagField == 0) ? 0 : 10) << endl
        << "      stemax = 0" << endl
        << "      deemax = 0" << endl
        << "      epsil = 1e-3" << endl
        << "      stmin = 0" << endl
        << "      nwbuf = 0" << endl
        << "      call gstmed(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,"
        << endl
        << "     +            stemax,deemax,epsil,stmin,ubuf,nwbuf)" << endl;

   double tocm, todeg;
   getConversions(el, tocm, todeg);

   double par[99];
   int npar = 0;
   if (shapeS.equals("box"))
   {
      shapeS = "BOX ";
      double xl, yl, zl;
      char* xyzStr = el.getAttribute("X_Y_Z").transcode();
      sscanf(xyzStr, "%lf %lf %lf", &xl, &yl, &zl);
      delete [] xyzStr;

      npar = 3;
      par[0] = xl/2 * tocm;
      par[1] = yl/2 * tocm;
      par[2] = zl/2 * tocm;
   }
   else if (shapeS.equals("tubs"))
   {
      shapeS = "TUBS";
      double ri, ro, zl, phi0, dphi;
      char* riozStr = el.getAttribute("Rio_Z").transcode();
      sscanf(riozStr, "%lf %lf %lf", &ri, &ro, &zl);
      delete [] riozStr;
      char* profStr = el.getAttribute("profile").transcode();
      sscanf(profStr, "%lf %lf", &phi0, &dphi);
      delete [] profStr;

      npar = 5;
      par[0] = ri * tocm;
      par[1] = ro * tocm;
      par[2] = zl/2 * tocm;
      par[3] = phi0 * todeg;
      par[4] = (phi0 + dphi) * todeg;
      if (dphi == 360)
      {
         shapeS = "TUBE";
         npar = 3;
      }
   }
   else if (shapeS.equals("trd"))
   {
      shapeS = "TRAP";
      double xm, ym, xp, yp, zl;
      char* xyzStr = el.getAttribute("Xmp_Ymp_Z").transcode();
      sscanf(xyzStr, "%lf %lf %lf %lf %lf", &xm, &xp, &ym, &yp, &zl);
      delete [] xyzStr;
      double alph_xz, alph_yz;
      char* incStr = el.getAttribute("inclination").transcode();
      sscanf(incStr, "%lf %lf", &alph_xz, &alph_yz);
      delete [] incStr;

      npar = 11;
      double x = tan(alph_xz * M_PI/180);
      double y = tan(alph_yz * M_PI/180);
      double r = sqrt(x*x + y*y);
      par[0] = zl/2 * tocm;
      par[1] = atan2(r,1) * 180/M_PI;
      par[2] = atan2(y,x) * 180/M_PI;
      par[3] = ym/2 * tocm;
      par[4] = xm/2 * tocm;
      par[5] = xm/2 * tocm;
      par[6] = 0;
      par[7] = yp/2 * tocm;
      par[8] = xp/2 * tocm;
      par[9] = xp/2 * tocm;
      par[10] = 0;
   }
   else if (shapeS.equals("pcon"))
   {
      shapeS = "PCON";
      double phi0, dphi;
      char* profStr = el.getAttribute("profile").transcode();
      sscanf(profStr, "%lf %lf", &phi0, &dphi);
      DOM_NodeList planeList = el.getElementsByTagName("polyplane");
      delete [] profStr;

      npar = 3;
      par[0] = phi0 * todeg;
      par[1] = dphi * todeg;
      par[2] = planeList.getLength();
      for (int p = 0; p < planeList.getLength(); p++)
      {
         double ri, ro, zl;
         DOM_Node node = planeList.item(p);
         DOM_Element elem = (DOM_Element&) node;
         char* riozStr = elem.getAttribute("Rio_Z").transcode();
         sscanf(riozStr, "%lf %lf %lf", &ri, &ro, &zl);
         delete [] riozStr;
         par[npar++] = zl * tocm;
         par[npar++] = ri * tocm;
         par[npar++] = ro * tocm;
      }
   }
   else if (shapeS.equals("cons"))
   {
      shapeS = "CONS";
      double rim, rip, rom, rop, zl;
      char* riozStr = el.getAttribute("Rio1_Rio2_Z").transcode();
      sscanf(riozStr, "%lf %lf %lf %lf %lf", &rim, &rom, &rip, &rop, &zl);
      delete [] riozStr;
      double phi0, dphi;
      char* profStr = el.getAttribute("profile").transcode();
      sscanf(profStr, "%lf %lf", &phi0, &dphi);
      delete [] profStr;

      npar = 7;
      par[0] = zl/2 * tocm;
      par[1] = rim * tocm;
      par[2] = rom * tocm;
      par[3] = rip * tocm;
      par[4] = rop * tocm;
      par[5] = phi0 * todeg;
      par[6] = (phi0 + dphi) * todeg;
      if (dphi == 360)
      {
         shapeS = "CONE";
         npar = 5;
      }
   }
   else
   {
      char* shapeStr = shapeS.transcode();
      cerr << "hdds-geant error: volume " << nameStr
           << " should be one of the valid shapes, not " << shapeStr << endl;
      exit(1);
   }

   if (strlen(nameStr) > 4)
   {
      cerr << "hdds-geant error: volume name " << nameStr
           << " should be no more than 4 characters long." << endl;
      exit(1);
   }

   char* shapeStr = shapeS.transcode();
   cout << endl
        << "      chname = \'" << nameStr << "\'" << endl
        << "      chshap = \'" << shapeStr << "\'" << endl
        << "      nmed = " << itmed << endl
        << "      npar = " << npar << endl;
   for (int ipar = 0; ipar < npar; ipar++)
   {
      cout << "      par(" << ipar + 1 << ") = " << par[ipar] << endl;
   }
   cout << "      call gsvolu(chname,chshap,nmed,par,npar,ivolu)" << endl;

   delete [] shapeStr;
   delete [] nameStr;
   delete [] matStr;

   char ivoluStr[10];
   int ivolu = ++Refsys::fVolumes;
   sprintf(ivoluStr, "%d", ivolu);
   el.setAttribute("Geant3ivolu", ivoluStr);  
   el.setAttribute("Geant3icopy", "0");  

/* consistency check #1: require Geant's volume index to match mine
 * 
 * This is required if the getX() lookup functions are going to work.
 * I count volumes in the order I define them, starting from 1.  If
 * Geant does the same thing then this error should never occur.
 */
   cout << "      if (ivolu.ne." << ivolu << ")"
        << " stop \'consistency check #1 failed\'" << endl;

   return Refsys::fVolumes;
}

int Refsys::createRotation()
{
   if (fRotation < 0)
   {
      fRotation = ++fRotations;
   }
   else
   {
      return fRotation;
   }

   cout << endl
        << "      irot = " << fRotation << endl;

   for (int i = 0; i < 3; i++)
   {
      double theta, phi;
      double r = sqrt(fRmatrix[0][i] * fRmatrix[0][i]
                    + fRmatrix[1][i] * fRmatrix[1][i]);
      theta = atan2(r, fRmatrix[2][i]) * 180/M_PI;
      phi = atan2(fRmatrix[1][i], fRmatrix[0][i]) * 180/M_PI;
      cout << "      theta" << i + 1 << " = " << theta << endl
           << "      phi" << i + 1 << " = " << phi << endl;
   }

   cout << "      "
        << "call gsrotm(irot,theta1,phi1,theta2,phi2,theta3,phi3)"
        << endl;

   return fRotation;
}

int Refsys::createDivision(char* divStr,
                           int ncopy,
                           int iaxis,
                           double start,
                           double step)
{
   divStr[0] = toupper(divStr[0]);
   divStr[1] = toupper(divStr[1]);
   divStr[2] = toupper(divStr[2]);
   divStr[3] = toupper(divStr[3]);

   assert (fMother != 0);
   char* motherStr = fMother.getAttribute("name").transcode();

   cout << endl
        << "      chname = \'" << divStr << "\'" << endl
        << "      chmoth = \'" << motherStr << "\'" << endl
        << "      ndiv = " << ncopy << endl
        << "      iaxis = " << iaxis << endl
        << "      step = " << step << endl
        << "      c0 = " << start << endl
        << "      numed = 0" << endl
        << "      ndvmax = 0" << endl
        << "      call gsdvx(chname,chmoth,ndiv,iaxis,step,c0,numed,ndvmax)"
        << endl;
   delete [] motherStr;

   char attStr[30];
   DOM_Document document = fMother.getOwnerDocument();
   DOM_Element divEl = document.createElement("Geant3division");
   divEl.setAttribute("name", divStr);
   divEl.setAttribute("volume", motherStr);
   sprintf(attStr, "%d", ++Refsys::fVolumes);
   divEl.setAttribute("Geant3ivolu", attStr);  
   sprintf(attStr, "%d", ncopy);
   divEl.setAttribute("Geant3icopy", attStr);  
   fMother.appendChild(divEl);
   fMother = divEl;

   for (int id = 0; id < fIdentifiers; id++)
   {
      DOMString fieldS = fIdentifier[id].fieldS;
      int value = fIdentifier[id].value;
      int step = fIdentifier[id].step;
      DOMString idlistS;
      for (int ic = 0; ic < ncopy; ic++)
      {
         char str[30];
         sprintf(str, "%d ", value);
         idlistS += str;
         value += step;
      }
      divEl.setAttribute(fieldS, idlistS);
   }
   fIdentifiers = 0;
   return ncopy;
}

int Refsys::createVolume(DOM_Element& el)
{
   int icopy = 0;

   Refsys myRef(*this);
   DOMString tagS = el.getTagName();
   DOMString nameS = el.getAttribute("name");
   char* nameStr = nameS.transcode();
   if (nameS.equals("fieldVolume"))
   {
      myRef.fMagField = 2;
   }
   else if (nameS.equals("sweepMagnet"))
   {
      myRef.fMagField = 3;
   }

   DOM_Element env;
   DOM_Document document = el.getOwnerDocument();
   DOMString envS = el.getAttribute("envelope");
   if (envS != 0)
   {
      env = document.getElementById(envS);
      DOMString containS = env.getAttribute("contains");
      if (containS.equals(nameS))
      {
         return myRef.createVolume(env);
      }
      else if (containS != 0)
      {
         char* envStr = envS.transcode();
         cerr << "hdds-geant error: re-use of shape " << envStr
              << " is not allowed by Geant3." << endl;
         exit(1);
      }
      env.setAttribute("contains",nameS);
      icopy = myRef.createVolume(env);
      myRef.fIdentifiers = 0;
      myRef.fMother = env;
      myRef.reset();
   }

   if (tagS.equals("intersection") ||
       tagS.equals("subtraction") ||
       tagS.equals("union"))
   {
      char* tagStr = tagS.transcode();
      cerr << "hdds-geant error: boolean " << tagStr
           << " operator is not supported by Geant3." << endl;
      exit(1);
   }
   else if (tagS.equals("composition"))
   {
      DOM_Node cont;
      int nSiblings = 0;
      for (cont = el.getFirstChild(); 
           cont != 0;
           cont = cont.getNextSibling())
      {
         if (cont.getNodeType() == ELEMENT_NODE)
         {
            ++nSiblings;
         }
      }

      for (cont = el.getFirstChild(); 
           cont != 0;
           cont = cont.getNextSibling())
      {
         if (cont.getNodeType() != ELEMENT_NODE)
         {
            continue;
         }
         DOM_Element contEl = (DOM_Element&) cont;
         DOMString comdS = contEl.getTagName();
         DOMString targS = contEl.getAttribute("volume");
         DOM_Element targEl = document.getElementById(targS);

         Refsys drs(myRef);
         double origin[3], angle[3];
         char* rotStr = contEl.getAttribute("rot").transcode();
         sscanf(rotStr, "%lf %lf %lf", &angle[0], &angle[1], &angle[2]);
         delete [] rotStr;
         double tocm, todeg;
         getConversions(contEl, tocm, todeg);
         double torad = todeg * M_PI/180;
         angle[0] *= torad;
         angle[1] *= torad;
         angle[2] *= torad;
         bool noRotation = (angle[0] == 0) &&
                           (angle[1] == 0) &&
                           (angle[2] == 0) ;

         DOM_Node ident;
         for (ident = cont.getFirstChild(); 
              ident != 0;
              ident = ident.getNextSibling())
         {
            if (ident.getNodeType() != ELEMENT_NODE)
            {
               continue;
            }
            int id = drs.fIdentifiers++;
            DOM_Element identEl = (DOM_Element&) ident;
            drs.fIdentifier[id].fieldS = identEl.getAttribute("field");
            char* fieldStr = identEl.getAttribute("field").transcode();
            char* valueStr = identEl.getAttribute("value").transcode();
            char* stepStr = identEl.getAttribute("step").transcode();
            drs.fIdentifier[id].value = atoi(valueStr);
            drs.fIdentifier[id].step = atoi(stepStr);
            if (fIdentifierList == 0)
            {
               fIdentifierList = new char[32];
               fIdentifierList[0] = 0;
               strncpy(fIdentifierList, fieldStr, 30);
            }
            else if (strstr(fIdentifierList, fieldStr) == 0)
            {
               int len = strlen(fIdentifierList);
               char* newList = new char[len + 32];
               strcpy(newList, fIdentifierList);
               strcat(newList, " ");
               strncat(newList, fieldStr, 30);
               delete [] fIdentifierList;
               fIdentifierList = newList;
            }
            delete [] fieldStr;
            delete [] valueStr;
            delete [] stepStr;
         }

         if (comdS.equals("posXYZ"))
         {
            char* xyzStr = contEl.getAttribute("X_Y_Z").transcode();
            sscanf(xyzStr, "%lf %lf %lf", &origin[0], &origin[1], &origin[2]);
            delete [] xyzStr;
            origin[0] *= tocm;
            origin[1] *= tocm;
            origin[2] *= tocm;
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS.equals("posRPhiZ"))
         {
            double r, phi, z;
            char* rphizStr = contEl.getAttribute("R_Phi_Z").transcode();
            sscanf(rphizStr, "%lf %lf %lf", &r, &phi, &z);
            delete [] rphizStr;
            double s;
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            phi *= torad;
            r *= tocm;
            z *= tocm;
            s *= tocm;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            DOMString implrotS = contEl.getAttribute("impliedRot");
            if (implrotS.equals("true") && (phi != 0))
            {
               angle[2] += phi;
            }
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS.equals("mposPhi"))
         {
            char* ncopyStr = contEl.getAttribute("ncopy").transcode();
            int ncopy = atoi(ncopyStr);
            delete ncopyStr;
            if (ncopy <= 0)
            {
               cerr << "hdds-geant error: volume " << nameStr
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double phi0, dphi;
            char* phiStr = contEl.getAttribute("Phi0").transcode();
            phi0 = atof(phiStr) * torad;
            delete [] phiStr;
            DOMString phiS = contEl.getAttribute("dPhi");
            if (phiS != 0)
            {
               char* phiStr = phiS.transcode();
               dphi = atof(phiStr) * torad;
               delete [] phiStr;
            }
            else
            {
               dphi = 2 * M_PI / ncopy;
            }

            double r, s, z;
            char* rzStr = contEl.getAttribute("R_Z").transcode();
            sscanf(rzStr, "%lf %lf", &r, &z);
            delete [] rzStr;
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            r *= tocm;
            z *= tocm;
            s *= tocm;

            DOMString containerS;
            if (env != 0)
            {
               containerS = env.getTagName();
            }
            else
            {
               containerS = el.getAttribute("divides");
            }
            DOMString implrotS = contEl.getAttribute("impliedRot");
            if (noRotation && (nSiblings == 1) &&
                (containerS.equals("pcon") ||
                 containerS.equals("cons") ||
                 containerS.equals("tubs")) &&
                implrotS.equals("true"))
            {
               static int phiDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "s%3.3x", ++phiDivisions);
               phi0 *= 180/M_PI;
               dphi *= 180/M_PI;
               DOM_Element targEnv(targEl);
               DOMString targEnvS = targEl.getAttribute("envelope");
               if (targEnvS != 0)
               {
                  targEnv = document.getElementById(targEnvS);
               }
               DOMString profS = targEnv.getAttribute("profile");
               if ((r == 0) && (profS != 0))
               {
                  double phi1, dphi1;
                  char* profStr = profS.transcode();
                  sscanf(profStr, "%lf %lf", &phi1, &dphi1);
                  delete [] profStr;
                  getConversions(targEnv, tocm, todeg);
                  double phiShift = phi1 + dphi1/2;
                  phi0 += phiShift * todeg;
                  phi1 -= phiShift;
                  profStr = new char[80];
                  sprintf(profStr, "%lf %lf", phi1, dphi1);
                  targEnv.setAttribute("profile",profStr);
                  delete [] profStr;
               }
               int iaxis = 2;
               drs.createDivision(divStr, ncopy, iaxis, phi0 - dphi/2, dphi);
               targEl.setAttribute("divides", containerS);
               origin[0] = r;
               origin[1] = s;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double phi = phi0 + inst * dphi;
                  origin[0] = r * cos(phi) - s * sin(phi);
                  origin[1] = r * sin(phi) + s * cos(phi);
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  if (implrotS.equals("true"))
                  {
                     angle[2] += ((inst == 0) ? phi0 : dphi);
                     drs.rotate(drs0, angle);
                  }
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposR"))
         {
            char* ncopyStr = contEl.getAttribute("ncopy").transcode();
            int ncopy = atoi(ncopyStr);
            delete ncopyStr;
            if (ncopy <= 0)
            {
               cerr << "hdds-geant error: volume " << nameStr
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double r0, dr;
            char* rStr = contEl.getAttribute("R0").transcode();
            r0 = atof(rStr) * tocm;
            delete [] rStr;
            rStr = contEl.getAttribute("dR").transcode();
            dr = atof(rStr) * tocm;
            delete [] rStr;

            double phi, z, s;
            char* zphiStr = contEl.getAttribute("Z_Phi").transcode();
            sscanf(zphiStr, "%lf %lf", &z, &phi);
            delete [] zphiStr;
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            phi *= torad;
            z *= tocm;
            s *= tocm;

            DOMString containerS;
            if (env != 0)
            {
               containerS = env.getTagName();
            }
            else
            {
               containerS = el.getAttribute("divides");
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS.equals("pcon") ||
                 containerS.equals("cons") ||
                 containerS.equals("tubs")))
            {
               static int rDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "r%3.3x", ++rDivisions);
               int iaxis = 1;
               drs.createDivision(divStr, ncopy, iaxis, r0 - dr/2, dr);
               targEl.setAttribute("divides", containerS);
               origin[0] = r0 * cos(phi) - s * sin(phi);
               origin[1] = r0 * sin(phi) + s * cos(phi);
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double r = r0 + inst * dr;
                  origin[0] = r * cos(phi) - s * sin(phi);
                  origin[1] = r * sin(phi) + s * cos(phi);
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposX"))
         {
            char* ncopyStr = contEl.getAttribute("ncopy").transcode();
            int ncopy = atoi(ncopyStr);
            delete ncopyStr;
            if (ncopy <= 0)
            {
               cerr << "hdds-geant error: volume " << nameStr
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double x0, dx;
            char* xStr = contEl.getAttribute("X0").transcode();
            x0 = atof(xStr) * tocm;
            delete [] xStr;
            xStr = contEl.getAttribute("dX").transcode();
            dx = atof(xStr) * tocm;
            delete [] xStr;

            double y, z, s;
            char* yzStr = contEl.getAttribute("Y_Z").transcode();
            sscanf(yzStr, "%lf %lf", &y, &z);
            delete [] yzStr;
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            y *= tocm;
            z *= tocm;
            s *= tocm;

            DOMString containerS;
            if (env != 0)
            {
               containerS = env.getTagName();
            }
            else
            {
               containerS = el.getAttribute("divides");
            }
            if (noRotation && (nSiblings == 1) && 
                containerS.equals("box"))
            {
               static int xDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "x%3.3x", ++xDivisions);
               int iaxis = 1;
               drs.createDivision(divStr, ncopy, iaxis, x0 - dx/2, dx);
               targEl.setAttribute("divides", containerS);
               origin[0] = 0;
               origin[1] = y + s;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double x = x0 + inst * dx;
                  origin[0] = x;
                  origin[1] = y + s;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposY"))
         {
            char* ncopyStr = contEl.getAttribute("ncopy").transcode();
            int ncopy = atoi(ncopyStr);
            delete ncopyStr;
            if (ncopy <= 0)
            {
               cerr << "hdds-geant error: volume " << nameStr
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double y0, dy;
            char* yStr = contEl.getAttribute("Y0").transcode();
            y0 = atof(yStr) * tocm;
            delete [] yStr;
            yStr = contEl.getAttribute("dY").transcode();
            dy = atof(yStr) * tocm;
            delete [] yStr;

            double x, z, s;
            char* xzStr = contEl.getAttribute("Z_X").transcode();
            sscanf(xzStr, "%lf %lf", &x, &z);
            delete [] xzStr;
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            x *= tocm;
            z *= tocm;
            s *= tocm;

            DOMString containerS;
            if (env != 0)
            {
               containerS = env.getTagName();
            }
            else
            {
               containerS = el.getAttribute("divides");
            }
            if (noRotation && (nSiblings == 1) && 
                containerS.equals("box"))
            {
               static int yDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "y%3.3x", ++yDivisions);
               int iaxis = 2;
               drs.createDivision(divStr, ncopy, iaxis, y0 - dy/2, dy);
               targEl.setAttribute("divides", containerS);
               origin[0] = x + s;
               origin[1] = 0;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double y = y0 + inst * dy;
                  double phi = atan2(y,x);
                  origin[0] = x;
                  origin[1] = y;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposZ"))
         {
            char* ncopyStr = contEl.getAttribute("ncopy").transcode();
            int ncopy = atoi(ncopyStr);
            delete ncopyStr;
            if (ncopy <= 0)
            {
               cerr << "hdds-geant error: volume " << nameStr
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double z0, dz;
            char* zStr = contEl.getAttribute("Z0").transcode();
            z0 = atof(zStr) * tocm;
            delete [] zStr;
            zStr = contEl.getAttribute("dZ").transcode();
            dz = atof(zStr) * tocm;
            delete [] zStr;

            double x, y, s;
            char* xyStr = contEl.getAttribute("X_Y").transcode();
            if (strlen(xyStr) > 0)
            {
               sscanf(xyStr, "%lf %lf", &x, &y);
               delete [] xyStr;
            }
            else
            {
               double r, phi;
               char* rphiStr = contEl.getAttribute("R_Phi").transcode();
               sscanf(rphiStr, "%lf %lf", &r, &phi);
               x = r * cos(phi * torad);
               y = r * sin(phi * torad);
               delete [] rphiStr;
            }
            char* sStr = contEl.getAttribute("S").transcode();
            s = atof(sStr);
            delete [] sStr;
            x *= tocm;
            y *= tocm;
            s *= tocm;

            DOMString containerS;
            if (env != 0)
            {
               containerS = env.getTagName();
            }
            else
            {
               containerS = el.getAttribute("divides");
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS != 0))
            {
               static int zDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "z%3.3x", ++zDivisions);
               int iaxis = 3;
               drs.createDivision(divStr, ncopy, iaxis, z0 - dz/2, dz);
               targEl.setAttribute("divides", containerS);
               double phi = atan2(y,x);
               origin[0] = x - s * sin(phi);
               origin[1] = y + s * cos(phi);
               origin[2] = 0;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               double phi = atan2(y,x);
               origin[0] = x - s * sin(phi);
               origin[1] = y + s * cos(phi);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double z = z0 + inst * dz;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else
         {
            char* comdStr = comdS.transcode();
            cerr << "hdds-geant error: composition of volume " << nameStr
                 << " contains unknown tag " << comdStr << endl;
            exit(1);
         }
      }
   }
   else if (tagS.equals("stackX") || 
            tagS.equals("stackY") ||
            tagS.equals("stackZ"))
   {
      cerr << "hdds-geant error: stacks are not supported by Geant3." << endl
           << "Use compositions instead." << endl;
      exit(1);
   }
   else
   {
      DOMString icopyS = el.getAttribute("Geant3icopy");
      if (icopyS != 0)
      {
         char* icopyStr = icopyS.transcode();
         icopy = atoi(icopyStr);
         delete [] icopyStr;
      }
      else
      {
         myRef.createSolid(el);
         icopy = 0;
      }

      if (myRef.fMother != 0)
      {
         char* motherStr = myRef.fMother.getAttribute("name").transcode();
         cout << endl
              << "      chname = \'" << nameStr << "\'" << endl
              << "      nr = " << ++icopy << endl
              << "      chmoth = \'" << motherStr << "\'" << endl
              << "      x = " << myRef.fOrigin[0] << endl
              << "      y = " << myRef.fOrigin[1] << endl
              << "      z = " << myRef.fOrigin[2] << endl
              << "      irot = " << myRef.createRotation() << endl
              << "      chonly = \'ONLY\'" << endl
              << "      call gspos(chname,nr,chmoth,x,y,z,irot,chonly)"
              << endl;
         delete [] motherStr;
      }

      char icopyStr[30];
      sprintf(icopyStr, "%d", icopy);
      el.setAttribute("Geant3icopy", icopyStr);
      for (int id = 0; id < myRef.fIdentifiers; id++)
      {
         DOMString fieldS = myRef.fIdentifier[id].fieldS;
         DOMString idlistS = el.getAttribute(fieldS);
         char* idlistStr = idlistS.transcode();
         char* token = idlistStr;
         int count = icopy;
         for (char* idStr = strsep(&token, " ");
              (idStr != 0) && (strlen(idStr) > 0);
              idStr = strsep(&token, " "))
         {
            count--;
         }
         delete [] idlistStr;
         for ( ; count > 1; --count)
         {
            idlistS += "0 ";
         }
         char str[30];
         sprintf(str, "%d ", myRef.fIdentifier[id].value);
         if (idlistS == 0)
         {
            idlistS = str;
         }
         else
         {
            idlistS += str;
         }
         el.setAttribute(fieldS, idlistS);
         myRef.fIdentifier[id].value += myRef.fIdentifier[id].step;
      }
   }

   delete [] nameStr;
   return icopy;
}
