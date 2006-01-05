/*
 *  hdds-geant :   an interface utility that reads in a HDDS document
 *		   (Hall D Detector Specification) and writes out a
 *		   GEANT-3 geometry description in the form of a
 *		   fortran subroutine.
 *
 *  Revision - Richard Jones, January 25, 2005.
 *   -added the sphere section as a new supported volume type
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
 *  contents of a composition named "*fieldVolume" I assign case 2.  For
 *  all contents of a composition named "*Magnet*" I assign case 3.
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

/*
 * FIX_XERCES_getElementById_BUG does a store/load cycle at parsing time
 * to fully instantiate entity references on the document tree.
 * See xerces-c++ bug 12800 at http://nagoya.apache.org
 */
#define FIX_XERCES_getElementById_BUG true


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

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()

double fieldStrength[] =
{
   0.0,	  // zero field regions
   0.0,   // inhomogenous field regions (unused)
   22.4,  // mapped field regions: solenoid (kG, approximate)
   2.0    // uniform field regions: sweep magnets (kG)
};

struct material_desc_
{
   int gindex;
   double wfact;
};

void usage()
{
    std::cerr
         << "\nUsage:\n"
            "    hdds-geant [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << std::endl;
}

class Refsys
{
 public:
   DOMElement* fMother;		// current mother volume element
   int fMagField;		// flag indicating magnetic field
   double fOrigin[3];		// x,y,z coordinate of volume origin (cm)
   double fPhiOffset;		// azimuthal angle of volume origin (deg)
   double fRmatrix[3][3];	// rotation matrix (daughter -> mother)
   int fRotation;		// unique Rmatrix flag

   static int fVolumes;		// total number of volumes to far
   static XString fIdentifierList; // list of identifier strings (space-sep)

   struct VolIdent_
   {
      XString fieldS;	
      int value;
      int step;
   };
   std::list<struct VolIdent_> fIdentifier;	// identifier tag list 

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

   int createMaterial(DOMElement* el);	// generate code for materials
   int createSolid(DOMElement* el);	// generate code for solids
   int createRotation();		// generate code for rotations
   int createDivision(XString& divStr,
                      int ncopy,
                      int iaxis,
                      double start,
                      double step);	// generate code for divisions
   int createVolume(DOMElement* el);	// generate code for placement

 private:
   static int fRotations;	// non-trivial rotations defined so far
};

class Units
{
 /* The Units class provides conversion constants for creating readable
  * units-aware code.  For example, to ensure that the quantity "velocity"
  * is provided in units of km/hr regardless of the internal representation
  * simply write:
  *	distance_in_km_per_hr = distance * unit.km/unit.hr;
  * where unit is a Units class instance that has been previously set.
  * The user of the Units class should generally treat its data members
  * as constants and use Units::getConversions() to manage the values. 
  */
 public:
   double s,ns,ms;
   double min,hr,days,weeks;
   double m,km,cm,mm,um,nm;
   double in,ft,miles,mils;
   double rad,mrad,urad;
   double deg,arcmin,arcsec;
   double eV,KeV,MeV,GeV,TeV;
   double g,kg,mg;
   double m2,cm2,mm2;
   double b,mb,ub,nb,pb;
   double l,ml,cm3;
   double G,kG,Tesla;
   double percent;

   Units();				// empty constructor
   Units(Units& u);			// copy constructor
   void getConversions(DOMElement* el);	// get conversion constants from tag

 private:
   void set_1s(double tu);
   void set_1cm(double lu);
   void set_1rad(double au);
   void set_1deg(double au);
   void set_1KeV(double eu);
   void set_1MeV(double eu);
   void set_1GeV(double eu);
   void set_1g(double mu);
   void set_1cm2(double l2u);
   void set_1cm3(double l3u);
   void set_1G(double bfu);
};


class FortranWriter
{
 public:
   FortranWriter() {};
   void header();
   void trailer();
   void getFunctions(DOMElement* el, XString& ident);
};


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
           << "hdds-geant: Error during initialization! :\n"
           << S(message) << std::endl;
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

   XString xmlFile;
   bool geantOutput = true;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
         break;

      if (strcmp(argV[argInd], "-v") == 0)
         geantOutput = false;
      else
         std::cerr
              << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << std::endl;
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd];

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(xmlFile,false);
#else
   DOMDocument* document = buildDOMDocument(xmlFile,false);
#endif
   if (document == 0)
   {
      std::cerr
           << "hdds-geant : Error parsing HDDS document, "
           << "cannot continue" << std::endl;
      return 1;
   }

   DOMNode* docEl;
   try {
      docEl = document->getDocumentElement();
   }
   catch (DOMException& e) {
      XString msgS(e.msg);
      std::cerr << "Woops " << S(msgS) << std::endl;
      return 1;
   }

   XString everythingS("everything");
   DOMElement* rootEl = document->getElementById(X(everythingS));
   if (rootEl == 0)
   {
      std::cerr
           << "hdds-geant : Error scanning HDDS document, " << std::endl
           << "  no element named \"everything\" found" << std::endl;
      return 1;
   }

   if (geantOutput)
   {
      Refsys mrs;
      FortranWriter fout;
      fout.header();
      mrs.createVolume(rootEl);
      fout.trailer();

      XString::size_type start;
      XString::size_type stop;
      for (start = 0; start < mrs.fIdentifierList.size(); start = stop+1)
      {
         stop = mrs.fIdentifierList.find(" ",start);
         stop = (stop == XString::npos)? mrs.fIdentifierList.size() : stop;
         XString identStr = mrs.fIdentifierList.substr(start,stop-start);
         fout.getFunctions(rootEl, identStr);
      }
      
   }

   XMLPlatformUtils::Terminate();
   return 0;
}


int Refsys::fVolumes = 0;
int Refsys::fRotations = 0;
XString Refsys::fIdentifierList;

Refsys::Refsys()			// empty constructor
 : fIdentifier(),
   fMother(0),
   fMagField(0),
   fPhiOffset(0)
{
   reset();
}

Refsys::Refsys(const Refsys& src)	// copy constructor
 : fIdentifier(src.fIdentifier),
   fMother(src.fMother),
   fMagField(src.fMagField),
   fPhiOffset(src.fPhiOffset)
{
   reset(src);
}

Refsys& Refsys::operator=(Refsys& src)	// copy operator (deep sematics)
{
   Refsys* dst = new Refsys(src);
   return *dst;
}

Refsys& Refsys::reset()			// reset origin, Rmatrix to null
{
   fOrigin[0] = fOrigin[1] = fOrigin[2] = 0;
   fRmatrix[0][0] = fRmatrix[1][1] = fRmatrix[2][2] = 1;
   fRmatrix[0][1] = fRmatrix[1][0] = fRmatrix[1][2] =
   fRmatrix[0][2] = fRmatrix[2][0] = fRmatrix[2][1] = 0;
   fRotation = 0;
   return *this;
}

Refsys& Refsys::reset(const Refsys& ref) // reset origin, Rmatrix to ref
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

Refsys& Refsys::shift(const double vector[3])  // translate origin
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
   if ( (omega[0] != 0) || (omega[1] != 0) || (omega[2] != 0) )
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

int Refsys::createMaterial(DOMElement* el)
{
   static int imateCount = 0;
   int imate = ++imateCount;
   std::stringstream imateStr;
   imateStr << imate;
   XString imateAttS("Geant3imate");
   XString imateS(imateStr.str());
   el->setAttribute(X(imateAttS),X(imateS));

   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString matS(el->getAttribute(X(nameAttS)));

   XString aAttS("a");
   XString aS(el->getAttribute(X(aAttS)));
   double a = atof(S(aS));

   XString zAttS("z");
   XString zS(el->getAttribute(X(zAttS)));
   double z = atof(S(zS));

   double dens = -1;
   double radl = -1;
   double absl = -1;
   XString realAttS("real");
   DOMNodeList* paramList = el->getElementsByTagName(X(realAttS));
   for (int ip = 0; ip < paramList->getLength(); ip++)
   {
      DOMNode* node = paramList->item(ip);
      DOMElement* elem = (DOMElement*) node;
      XString paramAttS("name");
      XString valueAttS("value");
      XString paramS(elem->getAttribute(X(paramAttS)));
      XString valueS(elem->getAttribute(X(valueAttS)));
      if (paramS == "density")
      {
         dens = atof(S(valueS));
      }
      else if (paramS == "radlen")
      {
         radl = atof(S(valueS));
      }
      else if (paramS == "abslen")
      {
         absl = atof(S(valueS));
      }
   }

   DOMDocument* document = el->getOwnerDocument();
   XString addmaterialTagS("addmaterial");
   DOMNodeList* compList = el->getElementsByTagName(X(addmaterialTagS));
   std::list<struct material_desc_*> matList;
   for (int ic = 0; ic < compList->getLength(); ic++)
   {
      material_desc_* desc = new material_desc_;
      matList.push_back(desc);
      DOMNode* node = compList->item(ic);
      DOMElement* elem = (DOMElement*) node;
      XString compAttS("material");
      XString compS(elem->getAttribute(X(compAttS)));
      DOMElement* compEl = document->getElementById(X(compS));
      XString imateAttS("Geant3imate");
      XString cS(compEl->getAttribute(X(imateAttS)));
      if (cS.size() != 0)
      {
         desc->gindex = atoi(S(cS));
      }
      else
      {
         desc->gindex = createMaterial(compEl);
      }

      XString natomsTagS("natoms");
      XString frmassTagS("fractionmass");
      DOMNodeList* atomList = elem->getElementsByTagName(X(natomsTagS));
      DOMNodeList* fracList = elem->getElementsByTagName(X(frmassTagS));
      if (atomList->getLength() == 1)
      {
         XString typeS(compEl->getTagName());
         if (typeS != "element")
         {
            std::cerr
                 << "hdds-geant: error processing composite " << S(matS)
                 << std::endl
                 << "natoms can only be specified for elements."
                 << std::endl;
            exit(1);
         }
         else if (dens < 0)
         {
            std::cerr
                 << "hdds-geant error: " << S(tagS) << " " << S(matS)
                 << " is missing a density specification." << std::endl;
            exit(1);
         }
         DOMNode* node = atomList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString nAttS("n");
         XString nS(elem->getAttribute(X(nAttS)));
         desc->wfact = atoi(S(nS));
      }
      else if (fracList->getLength() == 1)
      {
         DOMNode* node = fracList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString fractionAttS("fraction");
         XString fS(elem->getAttribute(X(fractionAttS)));
         desc->wfact = atof(S(fS));

         double rho = 0;
         XString realTagS("real");
         DOMNodeList* propList = compEl->getElementsByTagName(X(realTagS));
         for (int ip = 0; ip < propList->getLength(); ip++)
         {
            DOMNode* pnode = propList->item(ip);
            DOMElement* pelem = (DOMElement*) pnode;
            XString pAttS("name");
            XString vAttS("value");
            XString pS(pelem->getAttribute(X(pAttS)));
            XString vS(pelem->getAttribute(X(vAttS)));
            if (pS == "density")
            {
               rho = atof(S(vS));
            }
         }
         assert (rho > 0);
         if (dens < 0)
         {
            dens -= desc->wfact / rho;
         }
      }
      else
      {
         std::cerr
              << "hdds-geant error: " << S(tagS) << " " << S(matS)
              << " is missing some proportion data." << std::endl;
         exit(1);
      }
   }

   if (tagS == "element")
   {
      if (dens < 0)
      {
         std::cerr
              << "hdds-geant error: " << S(tagS) << " " << S(matS)
              << " is missing a density specification." << std::endl;
         exit(1);
      }
      if ((radl == -1) || (absl == -1))
      {
         std::cout
              << std::endl
              << "      imate = " << imate << std::endl
              << "      namate = \'" << S(matS) << "\'" << std::endl
              << "      a = " << a << std::endl
              << "      z = " << z << std::endl
              << "      dens = " << dens << std::endl
              << "      nlmat = 1" << std::endl
              << "      wmat(1) = 1" << std::endl
              << "      call gsmixt(imate,namate,a,z,dens,nlmat,wmat)"
              << std::endl;
      }
      else
      {
         std::cout
              << std::endl
              << "      imate = " << imate << std::endl
              << "      chnama = \'" << S(matS) << "\'" << std::endl
              << "      a = " << a << std::endl
              << "      z = " << z << std::endl
              << "      dens = " << dens << std::endl
              << "      radl = " << radl / (dens + 1e-30) << std::endl
              << "      absl = " << absl / (dens + 1e-30) << std::endl
              << "      nwbuf = 0" << std::endl
              << "      call gsmate(imate,chnama,a,z,dens,radl,absl,ubuf,nwbuf)"
              << std::endl;
      }
   }
   else
   {
      double rho = dens;
      if (rho < 0)
      {
         rho = -1/(dens + 0.999999);
      }
      std::cout
           << std::endl
           << "      imate = " << imate << std::endl
           << "      namate = \'" << S(matS) << "\'" << std::endl;
      std::list<material_desc_*>::iterator iter = matList.begin();
      for (int im = 0; im < matList.size(); im++, iter++)
      {
         std::cout
              << "      wmat(" << im + 1 << ") = " << (*iter)->wfact
              << std::endl
              << "      call gfmate(" << (*iter)->gindex << ",chnama,"
              << "amat(" << im + 1 << "),zmat(" << im + 1 << "),"
              << "dens,radl,absl,ubuf,nwbuf)" << std::endl;
      }
      iter = matList.begin();
      std::cout
           << "      dens = " << rho << std::endl
           << "      nlmat = " << (((*iter)->wfact < 1)? "" : "-")
	   << matList.size() << std::endl
           << "      call gsmixt(imate,namate,amat,zmat,dens,nlmat,wmat)"
           << std::endl;
   }
   std::list<material_desc_*>::iterator iter;
   for (iter = matList.begin(); iter != matList.end(); ++iter)
   {
      delete *iter;
   }
   
   if (dens < 0)
   {
      std::stringstream densStr;
      densStr << -1/(dens + 0.999999);
      XString realTagS("real");
      XString nameAttS("name");
      XString valueAttS("value");
      XString nameS("density");
      XString valueS(densStr.str());
      DOMElement* densEl = document->createElement(X(realTagS));
      densEl->setAttribute(X(nameAttS),X(nameS));
      densEl->setAttribute(X(valueAttS),X(valueS));
      el->appendChild(densEl);
   }

   return imate;
}

Units::Units()
{
   set_1s(1.);
   set_1cm(1.);
   set_1rad(1.);
   set_1MeV(1.);
   set_1g(1.);
   set_1cm2(1.);
   set_1cm3(1.);
   set_1G(1.);
   percent = 100;
}

Units::Units(Units& u)
{
   set_1s(1/u.s);
   set_1cm(1/u.cm);
   set_1rad(1/u.rad);
   set_1MeV(1/u.MeV);
   set_1g(1/u.g);
   set_1cm2(1/u.cm2);
   set_1cm3(1/u.cm3);
   set_1G(1/u.G);
   percent = 100;
}

void Units::set_1s(double tu)
{
   s=1/tu;
   ns=s*1e9; ms=s*1e3;
   min=s*60; hr=min*60; days=hr*24; weeks=days*7;
}

void Units::set_1cm(double lu)
{
   cm=1/lu;
   m=cm*1e-2; mm=m*1e3; um=m*1e6; nm=m*1e9; km=m*1e-3;
   in=cm/2.54; ft=in/12; miles=ft/5280; mils=in*1e3;
}

void Units::set_1rad(double au)
{
   rad=1/au;
   mrad=rad*1e3; urad=rad*1e6;
   deg=rad*180/M_PI; arcmin=deg*60; arcsec=arcmin*60;
}

void Units::set_1deg(double au)
{
   deg=1/au;
   arcmin=deg*60; arcsec=arcmin*60;
   rad=deg*M_PI/180; mrad=rad*1e3; urad=rad*1e6;
}

void Units::set_1MeV(double eu)
{
   MeV=1/eu;
   eV=MeV*1e6; KeV=eV*1e-3; GeV=eV*1e-9; TeV=eV*1e-12;
}

void Units::set_1g(double mu)
{
   g=1/mu;
   kg=g*1e-3; mg=g*1e3;
}

void Units::set_1cm2(double l2u)
{
   cm2=1/l2u;
   m2=cm2*1e-4; mm2=cm2*1e2;
   b=cm2*1e24; mb=b*1e3; ub=b*1e6; nb=b*1e9; pb=b*1e12;
}

void Units::set_1cm3(double l3u)
{
   cm3=1/l3u;
   ml=cm3; l=ml*1e-3;
}

void Units::set_1G(double bfu)
{
   G=1/bfu;
   kG=G*1e-3; Tesla=kG*1e-1;
}

void Units::getConversions(DOMElement* el)
{
   XString unitlAttS("unit_length");
   XString unitlS(el->getAttribute(X(unitlAttS)));
   if (unitlS.size() == 0)
   {
      ;
   }
   else if (unitlS == "mm")
   {
      set_1cm(10);
   }
   else if (unitlS == "cm")
   {
      set_1cm(1);
   }
   else if (unitlS == "m")
   {
      set_1cm(0.01);
   }
   else
   {
      XString tagS(el->getTagName());
      std::cerr
           << "hdds-geant error: unknown length unit " << S(unitlS)
           << " on tag " << S(tagS) << std::endl;
      exit(1);
   }

   XString unitaAttS("unit_angle");
   XString unitaS = el->getAttribute(X(unitaAttS));
   if (unitaS.size() == 0)
   {
      ;
   }
   else if (unitaS == "deg")
   {
      set_1deg(1);
   }
   else if (unitaS == "rad")
   {
      set_1rad(1);
   }
   else if (unitaS == "mrad")
   {
      set_1rad(1e3);
   }
   else
   {
      XString tagS(el->getTagName());
      std::cerr
           << "hdds-geant error: unknown angle unit " << S(unitaS)
           << " on volume " << S(tagS) << std::endl;
      exit(1);
   }

   XString unitAttS("unit");
   XString unitS = el->getAttribute(X(unitAttS));
   if (unitS.size() == 0)
   {
      ;
   }
   else if (unitS == "mm")
   {
      set_1cm(10);
   }
   else if (unitS == "cm")
   {
      set_1cm(1);
   }
   else if (unitS == "m")
   {
      set_1cm(0.01);
   }
   else if (unitaS == "deg")
   {
      set_1deg(1);
   }
   else if (unitaS == "rad")
   {
      set_1rad(1);
   }
   else if (unitaS == "mrad")
   {
      set_1rad(1e3);
   }
   else if (unitS == "eV")
   {
      set_1MeV(1e6);
   }
   else if (unitS == "KeV")
   {
      set_1MeV(1e3);
   }
   else if (unitS == "MeV")
   {
      set_1MeV(1);
   }
   else if (unitS == "GeV")
   {
      set_1MeV(1e-3);
   }
   else if (unitS == "g/cm^2")
   {
      set_1g(1);
      set_1cm2(1);
   }
   else if (unitS == "g/cm^3")
   {
      set_1g(1);
      set_1cm3(1);
   }
   else if (unitS == "MeV/g/cm^2")
   {
      set_1MeV(1);
      set_1g(1);
      set_1cm2(1);
   }
   else if (unitS == "Tesla")
   {
      set_1G(1e-4);
   }
   else if (unitS == "kG")
   {
      set_1G(1e-3);
   }
   else if (unitS == "G")
   {
      set_1G(1);
   }
   else if (unitS == "percent")
   {
      ;
   }
   else if (unitS == "none")
   {
      ;
   }
   else
   {
      XString tagS(el->getTagName());
      std::cerr
           << "hdds-geant error: unknown unit " << S(unitS)
           << " on volume " << S(tagS) << std::endl;
      exit(1);
   }
}

int Refsys::createSolid(DOMElement* el)
{
   XString nameTagS("name");
   XString matAttS("material");
   XString nameS(el->getAttribute(X(nameTagS)));
   XString matS(el->getAttribute(X(matAttS)));

   int imate;
   DOMDocument* document = el->getOwnerDocument();
   DOMElement* matEl = document->getElementById(X(matS));
   XString imateAttS("Geant3imate");
   XString imateS(matEl->getAttribute(X(imateAttS)));
   if (imateS.size() != 0)
   {
      imate = atoi(S(imateS));
   }
   else
   {
      imate = createMaterial(matEl);
   }
   
   static int itmedCount = 0;
   int itmed = ++itmedCount;
   XString sensiAttS("sensitive");
   XString sensiS(el->getAttribute(X(sensiAttS)));
   std::cout
        << std::endl
        << "      itmed = " << itmed << std::endl
        << "      natmed = \'" << S(nameS) << " " << S(matS) << "\'"
        << std::endl
        << "      nmat = " << imate << std::endl
        << "      isvol = " << (sensiS == "true" ? 1 : 0) << std::endl
        << "      ifield = " << ((fMagField == 0) ? 0 : 2) << std::endl
        << "      fieldm = " << fieldStrength[fMagField] << std::endl
        << "      tmaxfd = " << ((fMagField == 0) ? 0 : 1) << std::endl
        << "      stemax = 0" << std::endl
        << "      deemax = 0" << std::endl
        << "      epsil = 1e-3" << std::endl
        << "      stmin = 0" << std::endl
        << "      nwbuf = 0" << std::endl
        << "      call gstmed(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,"
        << std::endl
        << "     +            stemax,deemax,epsil,stmin,ubuf,nwbuf)"
        << std::endl;

   Units unit;
   unit.getConversions(el);

   double par[99];
   int npar = 0;
   XString shapeS(el->getTagName());
   if (shapeS == "box")
   {
      shapeS = "BOX ";
      double xl, yl, zl;
      XString xyzAttS("X_Y_Z");
      XString xyzS(el->getAttribute(X(xyzAttS)));
      sscanf(S(xyzS), "%lf %lf %lf", &xl, &yl, &zl);

      npar = 3;
      par[0] = xl/2 * unit.cm;
      par[1] = yl/2 * unit.cm;
      par[2] = zl/2 * unit.cm;
   }
   else if (shapeS == "tubs")
   {
      shapeS = "TUBS";
      double ri, ro, zl, phi0, dphi;
      XString riozAttS("Rio_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      sscanf(S(riozS), "%lf %lf %lf", &ri, &ro, &zl);
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);

      npar = 5;
      par[0] = ri * unit.cm;
      par[1] = ro * unit.cm;
      par[2] = zl/2 * unit.cm;
      par[3] = phi0 * unit.deg;
      par[4] = (phi0 + dphi) * unit.deg;
      if (dphi*unit.deg == 360)
      {
         shapeS = "TUBE";
         npar = 3;
      }
   }
   else if (shapeS == "trd")
   {
      shapeS = "TRAP";
      double xm, ym, xp, yp, zl;
      XString xyzAttS("Xmp_Ymp_Z");
      XString xyzS(el->getAttribute(X(xyzAttS)));
      sscanf(S(xyzS), "%lf %lf %lf %lf %lf", &xm, &xp, &ym, &yp, &zl);
      double alph_xz, alph_yz;
      XString incAttS("inclination");
      XString incS(el->getAttribute(X(incAttS)));
      sscanf(S(incS), "%lf %lf", &alph_xz, &alph_yz);

      npar = 11;
      double x = tan(alph_xz * unit.rad);
      double y = tan(alph_yz * unit.rad);
      double r = sqrt(x*x + y*y);
      par[0] = zl/2 * unit.cm;
      par[1] = atan2(r,1) * unit.deg/unit.rad;
      par[2] = atan2(y,x) * unit.deg/unit.rad;
      par[3] = ym/2 * unit.cm;
      par[4] = xm/2 * unit.cm;
      par[5] = xm/2 * unit.cm;
      par[6] = 0;
      par[7] = yp/2 * unit.cm;
      par[8] = xp/2 * unit.cm;
      par[9] = xp/2 * unit.cm;
      par[10] = 0;
   }
   else if (shapeS == "pcon")
   {
      shapeS = "PCON";
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);
      XString planeTagS("polyplane");
      DOMNodeList* planeList = el->getElementsByTagName(X(planeTagS));

      npar = 3;
      par[0] = phi0 * unit.deg;
      par[1] = dphi * unit.deg;
      par[2] = planeList->getLength();
      for (int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozAttS("Rio_Z");
         XString riozS(elem->getAttribute(X(riozAttS)));
         sscanf(S(riozS), "%lf %lf %lf", &ri, &ro, &zl);
         par[npar++] = zl * unit.cm;
         par[npar++] = ri * unit.cm;
         par[npar++] = ro * unit.cm;
      }
   }
   else if (shapeS == "pgon")
   {
      shapeS = "PGON";
      int segments;
      XString segAttS("segments");
      XString segS(el->getAttribute(X(segAttS)));
      sscanf(S(segS), "%d", &segments);
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);
      XString planeTagS("polyplane");
      DOMNodeList* planeList = el->getElementsByTagName(X(planeTagS));

      npar = 4;
      par[0] = phi0 * unit.deg;
      par[1] = dphi * unit.deg;
      par[2] = segments;
      par[3] = planeList->getLength();
      for (int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozAttS("Rio_Z");
         XString riozS(elem->getAttribute(X(riozAttS)));
         sscanf(S(riozS), "%lf %lf %lf", &ri, &ro, &zl);
         par[npar++] = zl * unit.cm;
         par[npar++] = ri * unit.cm;
         par[npar++] = ro * unit.cm;
      }
   }
   else if (shapeS == "cons")
   {
      shapeS = "CONS";
      double rim, rip, rom, rop, zl;
      XString riozAttS("Rio1_Rio2_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      sscanf(S(riozS), "%lf %lf %lf %lf %lf", &rim, &rom, &rip, &rop, &zl);
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);

      npar = 7;
      par[0] = zl/2 * unit.cm;
      par[1] = rim * unit.cm;
      par[2] = rom * unit.cm;
      par[3] = rip * unit.cm;
      par[4] = rop * unit.cm;
      par[5] = phi0 * unit.deg;
      par[6] = (phi0 + dphi) * unit.deg;
      if (dphi*unit.deg == 360)
      {
         shapeS = "CONE";
         npar = 5;
      }
   }
   else if (shapeS == "sphere")
   {
      shapeS = "SPHE";
      double ri, ro;
      XString rioAttS("Rio");
      XString rioS(el->getAttribute(X(rioAttS)));
      sscanf(S(rioS), "%lf %lf", &ri, &ro);
      double theta0, theta1;
      XString polarAttS("polar_bounds");
      XString polarS(el->getAttribute(X(polarAttS)));
      sscanf(S(polarS), "%lf %lf", &theta0, &theta1);
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);

      npar = 6;
      par[0] = ri * unit.cm;
      par[1] = ro * unit.cm;
      par[2] = theta0 * unit.deg;
      par[3] = theta1 * unit.deg;
      par[4] = phi0 * unit.deg;
      par[5] = (phi0 + dphi) * unit.deg;
   }
   else
   {
      std::cerr
           << "hdds-geant error: volume " << S(nameS)
           << " should be one of the valid shapes, not " << S(shapeS)
           << std::endl;
      exit(1);
   }
   assert(npar < 99);

   if (nameS.size() > 4)
   {
      std::cerr
           << "hdds-geant error: volume name " << S(nameS)
           << " should be no more than 4 characters long." << std::endl;
      exit(1);
   }

   std::cout
        << std::endl
        << "      chname = \'" << S(nameS) << "\'" << std::endl
        << "      chshap = \'" << S(shapeS) << "\'" << std::endl
        << "      nmed = " << itmed << std::endl
        << "      npar = " << npar << std::endl;
   for (int ipar = 0; ipar < npar; ipar++)
   {
      std::cout
           << "      par(" << ipar + 1 << ") = " << par[ipar] << std::endl;
   }
   std::cout
        << "      call gsvolu(chname,chshap,nmed,par,npar,ivolu)" << std::endl;


   std::stringstream ivoluStr;
   int ivolu = ++Refsys::fVolumes;
   ivoluStr << ivolu;
   XString ivoluAttS("Geant3ivolu");
   XString icopyAttS("Geant3icopy");
   XString ivoluS(ivoluStr.str());
   XString icopyS("0");
   el->setAttribute(X(ivoluAttS),X(ivoluS));  
   el->setAttribute(X(icopyAttS),X(icopyS));  

/* consistency check #1: require Geant's volume index to match mine
 * 
 * This is required if the getX() lookup functions are going to work.
 * I count volumes in the order I define them, starting from 1.  If
 * Geant does the same thing then this error should never occur.
 */
   std::cout
        << "      if (ivolu.ne." << ivolu << ")"
        << " stop \'consistency check #1 failed\'" << std::endl;

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

   std::cout
        << std::endl
        << "      irot = " << fRotation << std::endl;

   for (int i = 0; i < 3; i++)
   {
      double theta, phi;
      double r = sqrt(fRmatrix[0][i] * fRmatrix[0][i]
                    + fRmatrix[1][i] * fRmatrix[1][i]);
      theta = atan2(r, fRmatrix[2][i]) * 180/M_PI;
      phi = atan2(fRmatrix[1][i], fRmatrix[0][i]) * 180/M_PI;
      std::cout
           << "      theta" << i + 1 << " = " << theta << std::endl
           << "      phi" << i + 1 << " = " << phi << std::endl;
   }

   std::cout
        << "      "
        << "call gsrotm(irot,theta1,phi1,theta2,phi2,theta3,phi3)"
        << std::endl;

   return fRotation;
}

int Refsys::createDivision(XString& divStr,
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

   XString nameAttS("name");
   XString motherS(fMother->getAttribute(X(nameAttS)));

   std::cout
        << std::endl
        << "      chname = \'" << divStr << "\'" << std::endl
        << "      chmoth = \'" << S(motherS) << "\'" << std::endl
        << "      ndiv = " << ncopy << std::endl
        << "      iaxis = " << iaxis << std::endl
        << "      step = " << step << std::endl
        << "      c0 = " << start << std::endl
        << "      numed = 0" << std::endl
        << "      ndvmax = 0" << std::endl
        << "      call gsdvx(chname,chmoth,ndiv,iaxis,step,c0,numed,ndvmax)"
        << std::endl;

   std::stringstream attStr;
   DOMDocument* document = fMother->getOwnerDocument();
   XString divTagS("Geant3division");
   DOMElement* divEl = document->createElement(X(divTagS));
   XString divS(divStr);
   XString voluAttS("volume");
   divEl->setAttribute(X(nameAttS),X(divS));
   divEl->setAttribute(X(voluAttS),X(motherS));
   XString ivoluAttS("Geant3ivolu");
   attStr << ++Refsys::fVolumes;
   XString ivoluS(attStr.str());
   divEl->setAttribute(X(ivoluAttS),X(ivoluS));  
   XString icopyAttS("Geant3icopy");
   std::stringstream copyStr;
   copyStr << ncopy;
   XString icopyS(copyStr.str());
   divEl->setAttribute(X(icopyAttS),X(icopyS));  
   fMother->appendChild(divEl);
   fMother = divEl;

   std::list<VolIdent_>::iterator iter;
   for (iter = fIdentifier.begin(); iter != fIdentifier.end(); ++iter)
   {
      XString fieldS(iter->fieldS);
      int value = iter->value;
      int step = iter->step;
      std::stringstream idlistStr;
      for (int ic = 0; ic < ncopy; ic++)
      {
         idlistStr << value << " ";
         value += step;
      }
      XString idlistS(idlistStr.str());
      divEl->setAttribute(X(fieldS),X(idlistS));
   }
   fIdentifier.clear();
   return ncopy;
}

int Refsys::createVolume(DOMElement* el)
{
   int icopy = 0;

   Refsys myRef(*this);
   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString nameS(el->getAttribute(X(nameAttS)));
   if (strstr(S(nameS),"fieldVolume"))
   {
      myRef.fMagField = 2;
   }
   else if (strstr(S(nameS),"Magnet"))
   {
      myRef.fMagField = 3;
   }

   DOMElement* env = 0;
   DOMDocument* document = el->getOwnerDocument();
   XString envAttS("envelope");
   XString envS(el->getAttribute(X(envAttS)));
   if (envS.size() != 0)
   {
      env = document->getElementById(X(envS));
      XString contAttS("contains");
      XString containS(env->getAttribute(X(contAttS)));
      if (containS == nameS)
      {
         return myRef.createVolume(env);
      }
      else if (containS.size() != 0)
      {
         std::cerr
              << "hdds-geant error: re-use of shape " << S(envS)
              << " is not allowed by Geant3." << std::endl;
         exit(1);
      }
      env->setAttribute(X(contAttS),X(nameS));
      icopy = myRef.createVolume(env);
      myRef.fIdentifier.clear();
      myRef.fMother = env;
      myRef.reset();
   }

   if (tagS == "intersection" ||
       tagS == "subtraction" ||
       tagS == "union")
   {
      std::cerr
           << "hdds-geant error: boolean " << S(tagS)
           << " operator is not supported by Geant3." << std::endl;
      exit(1);
   }
   else if (tagS == "composition")
   {
      DOMNode* cont;
      int nSiblings = 0;
      for (cont = el->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
         {
            ++nSiblings;
         }
      }

      for (cont = el->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() != DOMNode::ELEMENT_NODE)
         {
            continue;
         }
         DOMElement* contEl = (DOMElement*) cont;
         XString comdS(contEl->getTagName());
         XString voluAttS("volume");
         XString targS(contEl->getAttribute(X(voluAttS)));
         DOMElement* targEl = document->getElementById(X(targS));

         Refsys drs(myRef);
         double origin[3], angle[3];
         XString rotAttS("rot");
         XString rotS(contEl->getAttribute(X(rotAttS)));
         sscanf(S(rotS), "%lf %lf %lf", &angle[0], &angle[1], &angle[2]);
         Units unit;
         unit.getConversions(contEl);
         angle[0] *= unit.rad;
         angle[1] *= unit.rad;
         angle[2] *= unit.rad;
         bool noRotation = (angle[0] == 0) &&
                           (angle[1] == 0) &&
                           (angle[2] == 0) ;

         DOMNode* ident;
         for (ident = cont->getFirstChild(); 
              ident != 0;
              ident = ident->getNextSibling())
         {
            if (ident->getNodeType() != DOMNode::ELEMENT_NODE)
            {
               continue;
            }
            DOMElement* identEl = (DOMElement*) ident;
            XString fieldAttS("field");
            XString valueAttS("value");
            XString stepAttS("step");
            XString fieldS(identEl->getAttribute(X(fieldAttS)));
            XString valueS(identEl->getAttribute(X(valueAttS)));
            XString stepS(identEl->getAttribute(X(stepAttS)));
            VolIdent_ id;
            id.fieldS = fieldS;
            id.value = atoi(S(valueS));
            id.step = atoi(S(stepS));
            drs.fIdentifier.push_back(id);
            if (fIdentifierList.find(fieldS) == XString::npos)
            {
               if (fIdentifierList.size() > 0)
               {
                  fIdentifierList += " ";
               }
               fIdentifierList += fieldS;
            }
         }

         if (comdS == "posXYZ")
         {
            XString xyzAttS("X_Y_Z");
            XString xyzS(contEl->getAttribute(X(xyzAttS)));
            sscanf(S(xyzS), "%lf %lf %lf", &origin[0], &origin[1], &origin[2]);
            origin[0] *= unit.cm;
            origin[1] *= unit.cm;
            origin[2] *= unit.cm;
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS == "posRPhiZ")
         {
            double r, phi, z;
            XString rphizAttS("R_Phi_Z");
            XString rphizS(contEl->getAttribute(X(rphizAttS)));
            sscanf(S(rphizS), "%lf %lf %lf", &r, &phi, &z);
            double s;
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            phi *= unit.rad;
            r *= unit.cm;
            z *= unit.cm;
            s *= unit.cm;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            XString irotAttS("impliedRot");
            XString implrotS(contEl->getAttribute(X(irotAttS)));
            if (implrotS == "true" && (phi != 0))
            {
               angle[2] += phi;
            }
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS == "mposPhi")
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << "hdds-geant error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double phi0, dphi;
            XString phi0AttS("Phi0");
            XString phi0S(contEl->getAttribute(X(phi0AttS)));
            phi0 = atof(S(phi0S)) * unit.rad;
            XString dphiAttS("dPhi");
            XString dphiS(contEl->getAttribute(X(dphiAttS)));
            if (dphiS.size() != 0)
            {
               dphi = atof(S(dphiS)) * unit.rad;
            }
            else
            {
               dphi = 2 * M_PI / ncopy;
            }

            double r, s, z;
            XString rzAttS("R_Z");
            XString rzS(contEl->getAttribute(X(rzAttS)));
            sscanf(S(rzS), "%lf %lf", &r, &z);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            r *= unit.cm;
            z *= unit.cm;
            s *= unit.cm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            XString irotAttS("impliedRot");
            XString implrotS(contEl->getAttribute(X(irotAttS)));
            if (noRotation && (nSiblings == 1) &&
                (containerS == "pcon" ||
                 containerS == "cons" ||
                 containerS == "tubs") &&
                 implrotS == "true")
            {
               static int phiDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "s" << std::setfill('0') << std::setw(3) << std::hex
                      << ++phiDivisions;
               phi0 *= 180/M_PI;
               dphi *= 180/M_PI;
               XString envAttS("envelope");
               XString targEnvS(targEl->getAttribute(X(envAttS)));
               DOMElement* targEnv;
               if (targEnvS.size() != 0)
               {
                  targEnv = document->getElementById(X(targEnvS));
               }
               else
               {
                  targEnv = targEl;
               }
               XString profAttS("profile");
               XString profS(targEnv->getAttribute(X(profAttS)));
               if ((r == 0) && (profS.size() != 0))
               {
                  double phi1, dphi1;
                  sscanf(S(profS), "%lf %lf", &phi1, &dphi1);
                  Units tunit(unit);
                  tunit.getConversions(targEnv);
                  double phiOffset = (phi1 + dphi1/2) * tunit.deg;
		  drs.fPhiOffset = phiOffset;
                  phi0 += phiOffset;
               }
               int iaxis = 2;
               XString divS(divStr.str());
               drs.createDivision(divS, ncopy, iaxis,
                                  phi0 - dphi/2 - myRef.fPhiOffset, dphi);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
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
                  if (implrotS == "true")
                  {
                     angle[2] += ((inst == 0) ? phi0 : dphi);
                     drs.rotate(drs0, angle);
                  }
                  drs.createVolume(targEl);
                  std::list<VolIdent_>::iterator iter;
                  for (iter = drs.fIdentifier.begin();
                       iter != drs.fIdentifier.end();
                       ++iter)
                  {
                     iter->value += iter->step;
                  }
               }
            }
         }
         else if (comdS == "mposR")
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << "hdds-geant error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double r0, dr;
            XString r0AttS("R0");
            XString r0S(contEl->getAttribute(X(r0AttS)));
            r0 = atof(S(r0S)) * unit.cm;
            XString drAttS("dR");
            XString drS(contEl->getAttribute(X(drAttS)));
            dr = atof(S(drS)) * unit.cm;

            double phi, z, s;
            XString zphiAttS("Z_Phi");
            XString zphiS(contEl->getAttribute(X(zphiAttS)));
            sscanf(S(zphiS), "%lf %lf", &z, &phi);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            phi *= unit.rad;
            z *= unit.cm;
            s *= unit.cm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS == "pcon" ||
                 containerS == "cons" ||
                 containerS == "tubs"))
            {
               static int rDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "r" << std::setfill('0') << std::setw(3) << std::hex
                      << ++rDivisions;
               int iaxis = 1;
               XString divS(divStr.str());
               drs.createDivision(divS, ncopy, iaxis, r0 - dr/2, dr);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
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
                  std::list<VolIdent_>::iterator iter;
                  for (iter = drs.fIdentifier.begin();
                       iter != drs.fIdentifier.end();
                       ++iter)
                  {
                     iter->value += iter->step;
                  }
               }
            }
         }
         else if (comdS == "mposX")
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << "hdds-geant error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double x0, dx;
            XString x0AttS("X0");
            XString x0S(contEl->getAttribute(X(x0AttS)));
            x0 = atof(S(x0S)) * unit.cm;
            XString dxAttS("dX");
            XString dxS(contEl->getAttribute(X(dxAttS)));
            dx = atof(S(dxS)) * unit.cm;

            double y, z, s;
            XString yzAttS("Y_Z");
            XString yzS(contEl->getAttribute(X(yzAttS)));
            sscanf(S(yzS), "%lf %lf", &y, &z);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            y *= unit.cm;
            z *= unit.cm;
            s *= unit.cm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS == "box")
            {
               static int xDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "x" << std::setfill('0') << std::setw(3) << std::hex
                      << ++xDivisions;
               int iaxis = 1;
               XString divS(divStr.str());
               drs.createDivision(divS, ncopy, iaxis, x0 - dx/2, dx);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
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
                  std::list<VolIdent_>::iterator iter;
                  for (iter = drs.fIdentifier.begin();
                       iter != drs.fIdentifier.end();
                       ++iter)
                  {
                     iter->value += iter->step;
                  }
               }
            }
         }
         else if (comdS == "mposY")
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << "hdds-geant error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double y0, dy;
            XString y0AttS("Y0");
            XString y0S(contEl->getAttribute(X(y0AttS)));
            y0 = atof(S(y0S)) * unit.cm;
            XString dyAttS("dY");
            XString dyS(contEl->getAttribute(X(dyAttS)));
            dy = atof(S(dyS)) * unit.cm;

            double x, z, s;
            XString zxAttS("Z_X");
            XString zxS(contEl->getAttribute(X(zxAttS)));
            sscanf(S(zxS), "%lf %lf", &z, &x);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            x *= unit.cm;
            z *= unit.cm;
            s *= unit.cm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS == "box")
            {
               static int yDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "y" << std::setfill('0') << std::setw(3) << std::hex 
                      << ++yDivisions;
               int iaxis = 2;
               XString divS(divStr.str());
               drs.createDivision(divS, ncopy, iaxis, y0 - dy/2, dy);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
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
                  std::list<VolIdent_>::iterator iter;
                  for (iter = drs.fIdentifier.begin();
                       iter != drs.fIdentifier.end();
                       ++iter)
                  {
                     iter->value += iter->step;
                  }
               }
            }
         }
         else if (comdS == "mposZ")
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << "hdds-geant error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double z0, dz;
            XString z0AttS("Z0");
            XString z0S(contEl->getAttribute(X(z0AttS)));
            z0 = atof(S(z0S)) * unit.cm;
            XString dzAttS("dZ");
            XString dzS(contEl->getAttribute(X(dzAttS)));
            dz = atof(S(dzS)) * unit.cm;

            double x, y, s;
            XString xyAttS("X_Y");
            XString xyS(contEl->getAttribute(X(xyAttS)));
            if (xyS.size() > 0)
            {
               sscanf(S(xyS), "%lf %lf", &x, &y);
            }
            else
            {
               double r, phi;
               XString rphiAttS("R_Phi");
               XString rphiS(contEl->getAttribute(X(rphiAttS)));
               sscanf(S(rphiS), "%lf %lf", &r, &phi);
	       phi *= unit.rad;
               x = r * cos(phi);
               y = r * sin(phi);
            }
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            x *= unit.cm;
            y *= unit.cm;
            s *= unit.cm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS.size() != 0))
            {
               static int zDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "z" << std::setfill('0') << std::setw(3) << std::hex
                      << ++zDivisions;
               int iaxis = 3;
               XString divS(divStr.str());
               drs.createDivision(divS, ncopy, iaxis, z0 - dz/2, dz);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
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
                  std::list<VolIdent_>::iterator iter;
                  for (iter = drs.fIdentifier.begin();
                       iter != drs.fIdentifier.end();
                       ++iter)
                  {
                     iter->value += iter->step;
                  }
               }
            }
         }
         else
         {
            std::cerr
                 << "hdds-geant error: composition of volume " << S(nameS)
                 << " contains unknown tag " << S(comdS) << std::endl;
            exit(1);
         }
      }
   }
   else if (tagS == "stackX" || 
            tagS == "stackY" ||
            tagS == "stackZ")
   {
      std::cerr
           << "hdds-geant error: stacks are not supported by Geant3."
           << std::endl
           << "Use compositions instead." << std::endl;
      exit(1);
   }
   else
   {
      XString icopyAttS("Geant3icopy");
      XString icopyS(el->getAttribute(X(icopyAttS)));
      if (icopyS.size() != 0)
      {
         icopy = atoi(S(icopyS));
      }
      else
      {
         XString profAttS("profile");
         XString profS(el->getAttribute(X(profAttS)));
         if (profS.size() != 0)
         {
            double phi0, dphi;
            sscanf(S(profS), "%lf %lf", &phi0, &dphi);
            Units punit;
            punit.getConversions(el);
            phi0 *= punit.deg;
            dphi *= punit.deg;
	    if ( (myRef.fOrigin[0] == 0) && (myRef.fOrigin[1] == 0) )
	    {
               phi0 -= myRef.fPhiOffset;
	    }
            std::stringstream pStr;
            pStr << phi0 << " " << dphi;
            XString pS(pStr.str());
            el->setAttribute(X(profAttS),X(pS));
	 }
         myRef.createSolid(el);
         icopy = 0;
      }

      if (myRef.fMother != 0)
      {
         XString nameAttS("name");
         XString motherS(myRef.fMother->getAttribute(X(nameAttS)));
         int irot = myRef.createRotation();
         std::cout
              << std::endl
              << "      chname = \'" << S(nameS) << "\'" << std::endl
              << "      nr = " << ++icopy << std::endl
              << "      chmoth = \'" << S(motherS) << "\'" << std::endl
              << "      x = " << myRef.fOrigin[0] << std::endl
              << "      y = " << myRef.fOrigin[1] << std::endl
              << "      z = " << myRef.fOrigin[2] << std::endl
              << "      irot = " << irot << std::endl
              << "      chonly = \'ONLY\'" << std::endl
              << "      call gspos(chname,nr,chmoth,x,y,z,irot,chonly)"
              << std::endl;
      }

      std::stringstream icopyStr;
      icopyStr << icopy;
      XString copyS(icopyStr.str());
      XString copyAttS("Geant3icopy");
      el->setAttribute(X(copyAttS),X(copyS));
      std::list<VolIdent_>::iterator iter;
      for (iter = myRef.fIdentifier.begin();
           iter != myRef.fIdentifier.end();
           ++iter)
      {
         XString fieldS(iter->fieldS);
         XString idlistS(el->getAttribute(X(fieldS)));
         XString spaceS(" ");
         XMLStringTokenizer picker(X(idlistS),X(spaceS));
         int count = icopy;
         XString idS;
         for (idS = picker.nextToken(); idS.size() != 0;
                                        idS = picker.nextToken())
         {
            count--;
         }
         XString mylistS(idlistS);
         for ( ; count > 1; --count)
         {
            mylistS += "0 ";
         }
         std::stringstream str;
         str << iter->value << " ";
         mylistS += str.str();
         el->setAttribute(X(fieldS),X(mylistS));
      }
   }
   return icopy;
}

void FortranWriter::header()
{
   std::cout
        << "*" 							<< std::endl
        << "* HDDSgeant3 - fortran geometry definition package" << std::endl
        << "*              for the Hall D experiment."		<< std::endl
        << "*"							<< std::endl
        << "*         WARNING: DO NOT EDIT THIS FILE"		<< std::endl
        << "*"							<< std::endl
        << "* This file was generated automatically from the"	<< std::endl
        << "* HDDS xml geometry definition by the hdds-geant"	<< std::endl
        << "* translator.  Any changes made to this file will"	<< std::endl
        << "* disappear as soon as it is regenerated from the"	<< std::endl
        << "* xml source.  To introduce Geant3 optimizations,"	<< std::endl
        << "* see the subroutine Goptimize() in goptimize.F."	<< std::endl
        << "*"							<< std::endl
        << "      subroutine HDDSgeant3" 			<< std::endl
        << "      implicit none"				<< std::endl
        << "      integer imate"				<< std::endl
        << "      character*20 chnama,namate"			<< std::endl
        << "      real a,z,dens,radl,absl,ubuf(99)"		<< std::endl
        << "      integer nwbuf"				<< std::endl
        << "      real amat(99),zmat(99),wmat(99)"		<< std::endl
        << "      integer nlmat"				<< std::endl
        << "      integer itmed"				<< std::endl
        << "      character*20 natmed"				<< std::endl
        << "      integer nmat,isvol,ifield"			<< std::endl
        << "      real fieldm,tmaxfd,stemax,deemax,epsil,stmin" << std::endl
        << "      character*4 chname,chshap,chmoth"		<< std::endl
        << "      integer nmed,npar,ivolu"			<< std::endl
        << "      real par(99)"					<< std::endl
        << "      integer irot"					<< std::endl
        << "      real theta1,phi1,theta2,phi2,theta3,phi3"	<< std::endl
        << "      integer nr,ndiv,iaxis,numed,ndvmax"		<< std::endl
        << "      real step,c0"					<< std::endl
        << "      real x,y"					<< std::endl
        << "      character*4 chonly"				<< std::endl;
}

void FortranWriter::trailer()
{
   std::cout << "      end"					<< std::endl;
}

void FortranWriter::getFunctions(DOMElement* el, XString& ident)
{
   std::vector<int> table;
   int* start = new int[Refsys::fVolumes + 1];

   XString funcNameStr;
   XString identCaps(ident);
   identCaps[0] = toupper(identCaps[0]);
   funcNameStr = "get" + identCaps;
   XString identS(ident);
   XString wildS("*");
   DOMNodeList* alltagsList = el->getOwnerDocument()
                                ->getElementsByTagName(X(wildS));
   for (int itag = 0; itag < alltagsList->getLength(); itag++)
   {
      DOMNode* node = alltagsList->item(itag);
      DOMElement* elem = (DOMElement*) node;
      XString icopyAttS("Geant3icopy");
      XString ivoluAttS("Geant3ivolu");
      XString icopyS(elem->getAttribute(X(icopyAttS)));
      XString ivoluS(elem->getAttribute(X(ivoluAttS)));
      if (ivoluS.size() != 0)
      {
         int ivolu = atoi(S(ivoluS));
         int icopy = atoi(S(icopyS));
         XString idlistS(elem->getAttribute(X(identS)));
         if (idlistS.size() != 0)
         {
            XString spaceS(" ");
            XMLStringTokenizer picker(X(idlistS),X(spaceS));
            start[ivolu] = table.size() + 1;
            XString idS;
            for (idS = picker.nextToken(); idS.size() != 0;
                 idS = picker.nextToken())
            {
               table.push_back(atoi(S(idS)));
               --icopy;
            }
            for (; icopy > 0; --icopy)
            {
               table.push_back(0);
            }
         }
         else
         {
            start[ivolu] = 0;
         }
      }
   }

   std::cout
        << std::endl
        << "      function " << funcNameStr << "()" << std::endl
        << "      implicit none" << std::endl
        << "      integer " << funcNameStr << std::endl
        << "      integer nlevel,names,number,lvolum" << std::endl
        << "      common /gcvolu/nlevel,names(15),number(15),lvolum(15)"
        << std::endl;

   if (table.size() > 0)
   {
      std::cout
           << "      integer i,istart(" << Refsys::fVolumes << ")"
           << std::endl;

      for (int i = 0; i < Refsys::fVolumes;)
      {
         std::stringstream str;
         if (i % 100 == 0)
         {
            int ilimit = i + 100;
            ilimit = (ilimit > Refsys::fVolumes)? Refsys::fVolumes : ilimit;
            std::cout << "      data (istart(i),i=" << i + 1 << "," << ilimit
                 << ") /" << std::endl;
         }
         if (i % 10 == 0)
         {
            std::cout << "     + ";
         }
         str << std::setw(5) << start[++i];
         std::cout << str.str();
         if (i == Refsys::fVolumes)
         {
            std::cout << "/" << std::endl;
         }
         else if (i % 100 == 0)
         {
            std::cout << "/" << std::endl;
         }
         else if (i % 10 == 0)
         {
            std::cout << "," << std::endl;
         }
         else
         {
            std::cout << ",";
         }
      }

      std::cout << "      integer lookup(" << table.size() << ")" << std::endl;

      for (int i = 0; i < table.size();)
      {
         std::stringstream str;
         if (i % 100 == 0)
         {
            int ilimit = i + 100;
            ilimit = (ilimit > table.size())? table.size() : ilimit;
            std::cout << "      data (lookup(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         }
         if (i % 10 == 0)
         {
            std::cout << "     + ";
         }
         str << std::setw(5) << table[i++];
         std::cout << str.str();
         if (i == table.size())
         {
            std::cout << "/" << std::endl;
         }
         else if (i % 100 == 0)
         {
            std::cout << "/" << std::endl;
         }
         else if (i % 10 == 0)
         {
            std::cout << "," << std::endl;
         }
         else
         {
            std::cout << ",";
         }
      }

      std::cout
           << "      integer level,index" << std::endl
           << "      integer " << ident << std::endl
           << "      " << funcNameStr << " = 0" << std::endl
           << "      do level=1,nlevel" << std::endl
           << "        index = istart(lvolum(level))" << std::endl
           << "        if (index.gt.0) then" << std::endl
           << "          " << ident
           << " = lookup(index + number(level) - 1)" << std::endl
           << "          if (" << ident << ".gt.0) then" << std::endl
           << "            " << funcNameStr << " = " << ident
           << std::endl
           << "          endif" << std::endl
           << "        endif" << std::endl
           << "      enddo" << std::endl;
   }
   else
   {
      std::cout << "      " << funcNameStr << " = 0" << std::endl;
   }
   std::cout << "      end" << std::endl;
   delete [] start;
}
