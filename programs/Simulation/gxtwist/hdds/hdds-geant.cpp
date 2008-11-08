/*
 *  hdds-geant :   an interface utility that reads in a HDDS document
 *                   (Hall D Detector Specification) and writes out a
 *                   GEANT-3 geometry description in the form of a
 *                   fortran subroutine.
 *
 *  Revision - Richard Jones, November 25, 2006.
 *   -added output of optical properties for materials with optical
 *    properties defined
 *
 *  Revision - Richard Jones, January 25, 2005.
 *   -added the sphere section as a new supported volume type
 *
 *  Original version - Richard Jones, May 19 2001.
 *
 *  Notes:
 *  ------
 * 1. Output is sent to standard out through the ordinary c++ i/o library.
 * 2. As a by-product of using the DOM parser to access the xml source,
 *    hdds-geant verifies the source against the schema before translating it.
 *    Therefore it may also be used as a validator of the xml specification
 *    (see the -v option).
 */

#define APP_NAME "hdds-geant"

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
#include "hddsCommon.hpp"

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
         << "Usage:    " << APP_NAME << " [-v] {HDDS file}"
         << std::endl <<  "Options:" << std::endl
         << "    -v   validate only" << std::endl;
}

class FortranWriter : public CodeWriter
{
 public:
   void createHeader();
   void createTrailer();
   int createMaterial(DOMElement* el);   // generate code for materials
   int createSolid(DOMElement* el,
                   Refsys& ref);    	 // generate code for solids
   int createRotation(Refsys& ref);      // generate code for rotations
   int createRegion(DOMElement* el,
                    Refsys& ref);        // generate code for regions
   int createVolume(DOMElement* el,
                    Refsys& ref);  	 // generate code for placement
   int createDivision(XString& divStr,
                      Refsys& ref);	 // generate code for divisions
   void createSetFunctions(DOMElement* el,
                  const XString& ident); // generate code for properties
   void createGetFunctions(DOMElement* el,
                  const XString& ident); // generate code for identifiers
   void createMapFunctions(DOMElement* el,
                  const XString& ident); // generate code for field maps
   void createUtilityFunctions(DOMElement* el,
                  const XString& ident); // generate utility functions
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
           << APP_NAME << " - error during initialization!"
           << std::endl << S(message) << std::endl;
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
           << APP_NAME << " - error parsing HDDS document, "
           << "cannot continue" << std::endl;
      return 1;
   }

   DOMNode* docEl;
   try {
      docEl = document->getDocumentElement();
   }
   catch (DOMException& e) {
      std::cerr << "Woops " << e.msg << std::endl;
      return 1;
   }

   DOMElement* rootEl = document->getElementById(X("everything"));
   if (rootEl == 0)
   {
      std::cerr
           << APP_NAME << " - error scanning HDDS document, " << std::endl
           << "  no element named \"everything\" found" << std::endl;
      return 1;
   }

   if (geantOutput)
   {
      FortranWriter fout;
      fout.translate(rootEl);
   }

   XMLPlatformUtils::Terminate();
   return 0;
}

#ifdef LINUX_CPUTIME_PROFILING
extern CPUtimer timer;
#endif

int FortranWriter::createMaterial(DOMElement* el)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createMaterial: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int imate = CodeWriter::createMaterial(el);

   double a = fSubst.getAtomicWeight();
   double z = fSubst.getAtomicNumber();
   double dens = fSubst.getDensity();
   double radl = fSubst.getRadLength();
   double absl = fSubst.getAbsLength();
   double coll = fSubst.getColLength();
   double dedx = fSubst.getMIdEdx();
   XString matS = fSubst.getName();

   if (fSubst.fBrewList.size() == 0)
   {
      if ((radl == 0) || (absl == 0))
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
      std::cout
           << std::endl
           << "      imate = " << imate << std::endl
           << "      namate = \'" << S(matS) << "\'" << std::endl;
      std::list<Substance::Brew>::iterator iter = fSubst.fBrewList.begin();
      for (int im = 0; im < fSubst.fBrewList.size(); im++, iter++)
      {
         std::cout
              << "      wmat(" << im + 1 << ") = "
              << ((iter->natoms > 0)? (double)iter->natoms : iter->wfact)
              << std::endl
              << "      call gfmate(" << iter->sub->fUniqueID << ",chnama,"
              << "amat(" << im + 1 << "),zmat(" << im + 1 << "),"
              << "dens,radl,absl,ubuf,nwbuf)" << std::endl;
      }
      iter = fSubst.fBrewList.begin();
      std::cout
           << "      dens = " << dens << std::endl
           << "      nlmat = " << ((iter->natoms == 0)? "" : "-")
           << fSubst.fBrewList.size() << std::endl
           << "      call gsmixt(imate,namate,amat,zmat,dens,nlmat,wmat)"
           << std::endl;
   }
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
   return imate;
}

int FortranWriter::createSolid(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createSolid: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int ivolu = CodeWriter::createSolid(el,ref);
   int imate = fSubst.fUniqueID;
   static int itmedCount = 0;

   std::map<std::string,double> defaultPar;
   defaultPar["ifield"] = 0;	// default values for tracking properties
   defaultPar["fieldm"] = 0;	// are overridden by values in region tag
   defaultPar["tmaxfd"] = 0;
   defaultPar["stemax"] = 1;
   defaultPar["deemax"] = 0;
   defaultPar["epsil"] = 1e-3;
   defaultPar["stmin"] = 0;
   std::map<std::string,double>::iterator iter;
   for (iter = defaultPar.begin(); iter != defaultPar.end(); ++iter)
   {
      if (ref.fPar.find(iter->first) == ref.fPar.end())
      {
         ref.fPar[iter->first] = defaultPar[iter->first];
      }
   }

   int itmed = ++itmedCount;
   XString nameS(el->getAttribute(X("name")));
   XString matS(el->getAttribute(X("material")));
   XString sensiS(el->getAttribute(X("sensitive")));
   std::cout
        << std::endl
        << "      itmed = " << itmed << std::endl
        << "      natmed = \'" << S(nameS) << " " << S(matS) << "\'"
        << std::endl
        << "      nmat = " << imate << std::endl
        << "      isvol = " << (sensiS == "true" ? 1 : 0) << std::endl
        << "      ifield = " << ref.fPar["ifield"] << std::endl
        << "      fieldm = " << ref.fPar["fieldm"] << std::endl
        << "      tmaxfd = " << ref.fPar["tmaxfd"] << std::endl
        << "      stemax = " << ref.fPar["stemax"] << std::endl
        << "      deemax = " << ref.fPar["deemax"] << std::endl
        << "      epsil = " << ref.fPar["epsil"] << std::endl
        << "      stmin = " << ref.fPar["stmin"] << std::endl
        << "      nwbuf = 0" << std::endl
        << "      call gstmed(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,"
        << std::endl
        << "     +            stemax,deemax,epsil,stmin,ubuf,nwbuf)"
        << std::endl;

   DOMElement* matEl = el->getOwnerDocument()->getElementById(X(matS));
   DOMNodeList* propList = matEl->getElementsByTagName(X("optical_properties"));
   if (propList->getLength() > 0)
   {
      std::cout << "      call setoptical" << imate << "(itmed)" << std::endl;
   }

   Units unit;
   unit.getConversions(el);

   double par[99];
   int npar = 0;
   XString shapeS(el->getTagName());
   if (shapeS == "box")
   {
      shapeS = "BOX ";
      double xl, yl, zl;
      XString xyzS(el->getAttribute(X("X_Y_Z")));
      std::stringstream listr(xyzS);
      listr >> xl >> yl >> zl;

      npar = 3;
      par[0] = xl/2 * unit.cm;
      par[1] = yl/2 * unit.cm;
      par[2] = zl/2 * unit.cm;
   }
   else if (shapeS == "tubs")
   {
      shapeS = "TUBS";
      double ri, ro, zl, phi0, dphi;
      XString riozS(el->getAttribute(X("Rio_Z")));
      std::stringstream listr(riozS);
      listr >> ri >> ro >> zl;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;

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
      XString xyzS(el->getAttribute(X("Xmp_Ymp_Z")));
      std::stringstream listr(xyzS);
      listr >> xm >> xp >> ym >> yp >> zl;
      double alph_xz, alph_yz;
      XString incS(el->getAttribute(X("inclination")));
      listr.clear(), listr.str(incS);
      listr >> alph_xz >> alph_yz;

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
      XString profS(el->getAttribute(X("profile")));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
      DOMNodeList* planeList = el->getElementsByTagName(X("polyplane"));

      npar = 3;
      par[0] = phi0 * unit.deg;
      par[1] = dphi * unit.deg;
      par[2] = planeList->getLength();
      for (int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozS(elem->getAttribute(X("Rio_Z")));
         std::stringstream listr1(riozS);
         listr1 >> ri >> ro >> zl;
         par[npar++] = zl * unit.cm;
         par[npar++] = ri * unit.cm;
         par[npar++] = ro * unit.cm;
      }
   }
   else if (shapeS == "pgon")
   {
      shapeS = "PGON";
      int segments;
      XString segS(el->getAttribute(X("segments")));
      segments = atoi(S(segS));
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
      DOMNodeList* planeList = el->getElementsByTagName(X("polyplane"));

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
         XString riozS(elem->getAttribute(X("Rio_Z")));
         std::stringstream listr1(riozS);
         listr1 >> ri >> ro >> zl;
         par[npar++] = zl * unit.cm;
         par[npar++] = ri * unit.cm;
         par[npar++] = ro * unit.cm;
      }
   }
   else if (shapeS == "cons")
   {
      shapeS = "CONS";
      double rim, rip, rom, rop, zl;
      XString riozS(el->getAttribute(X("Rio1_Rio2_Z")));
      std::stringstream listr(riozS);
      listr >> rim >> rom >> rip >> rop >> zl;
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;

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
      XString rioS(el->getAttribute(X("Rio")));
      std::stringstream listr(rioS);
      listr >> ri >> ro;
      double theta0, theta1;
      XString polarS(el->getAttribute(X("polar_bounds")));
      listr.clear(), listr.str(polarS);
      listr >> theta0 >> theta1;
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;

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
           << APP_NAME << " error: volume " << S(nameS)
           << " should be one of the valid shapes, not " << S(shapeS)
           << std::endl;
      exit(1);
   }

   if (nameS.size() > 4)
   {
      std::cerr
           << APP_NAME << " error: volume name " << S(nameS)
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


/* consistency check #1: require Geant's volume index to match mine
 * 
 * This is required if the getX() lookup functions are going to work.
 * I count volumes in the order I define them, starting from 1.  If
 * Geant does the same thing then this error should never occur.
 */
   std::cout
        << "      if (ivolu.ne." << ivolu << ")"
        << " stop \'consistency check #1 failed\'" << std::endl;

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
   return ivolu;
}

int FortranWriter::createRotation(Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createRotation: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int irot = CodeWriter::createRotation(ref);

   if (irot > 0)
   {
      std::cout
           << std::endl
           << "      irot = " << irot << std::endl;
      for (int i = 0; i < 3; i++)
      {
         double theta, phi;
         double r = sqrt(ref.fRmatrix[0][i] * ref.fRmatrix[0][i]
                       + ref.fRmatrix[1][i] * ref.fRmatrix[1][i]);
         theta = atan2(r, ref.fRmatrix[2][i]) * 180/M_PI;
         phi = atan2(ref.fRmatrix[1][i], ref.fRmatrix[0][i]) * 180/M_PI;
         std::cout
              << "      theta" << i + 1 << " = " << theta << std::endl
              << "      phi" << i + 1 << " = " << phi << std::endl;
      }

      std::cout
           << "      "
           << "call gsrotm(irot,theta1,phi1,theta2,phi2,theta3,phi3)"
           << std::endl;
   }
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
   return irot;
}

int FortranWriter::createRegion(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createRegion: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int iregion = CodeWriter::createRegion(el,ref);

   if (ref.fRegion)
   {
      DOMNodeList* noBfieldL = ref.fRegion->getElementsByTagName(X("noBfield"));
      DOMNodeList* uniBfieldL = ref.fRegion->getElementsByTagName(X("uniformBfield"));
      DOMNodeList* mapBfieldL = ref.fRegion->getElementsByTagName(X("mappedBfield"));
      DOMNodeList* swimL = ref.fRegion->getElementsByTagName(X("swim"));
      if (noBfieldL->getLength() > 0)
      {
         ref.fPar["ifield"] = 0;
         ref.fPar["fieldm"] = 0;
         ref.fPar["tmaxfd"] = 0;
      }
      else if (uniBfieldL->getLength() > 0)
      {
         Units funit;
         DOMElement* uniBfieldEl = (DOMElement*)uniBfieldL->item(0);
         funit.getConversions(uniBfieldEl);
         XString bvecS(uniBfieldEl->getAttribute(X("Bx_By_Bz")));
         std::stringstream str(S(bvecS));
         double B[3];
         str >> B[0] >> B[1] >> B[2];
         ref.fPar["fieldm"] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
         ref.fPar["fieldm"] *= funit.kG;
         ref.fPar["ifield"] = 2;
         ref.fPar["tmaxfd"] = 1;
      }
      else if (mapBfieldL->getLength() > 0)
      {
         Units funit;
         DOMElement* mapBfieldEl = (DOMElement*)mapBfieldL->item(0);
         funit.getConversions(mapBfieldEl);
         XString bmaxS(mapBfieldEl->getAttribute(X("maxBfield")));
         ref.fPar["fieldm"] = atof(S(bmaxS));
         ref.fPar["fieldm"] *= funit.kG;
         ref.fPar["ifield"] = 2;
         ref.fPar["tmaxfd"] = 1;
      }
      if (swimL->getLength() > 0)
      {
         DOMElement* swimEl = (DOMElement*)swimL->item(0);
         XString methodS(swimEl->getAttribute(X("method")));
         ref.fPar["ifield"] = (methodS == "RungeKutta")? 1 : 2;
         Units unit;
         unit.getConversions(swimEl);
         XString stepS(swimEl->getAttribute(X("maxArcStep")));
         ref.fPar["tmaxfd"] = atof(S(stepS)) * unit.deg;
      }
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
   return iregion;
}

int FortranWriter::createDivision(XString& divStr, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createDivision: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   divStr[0] = toupper(divStr[0]);
   divStr[1] = toupper(divStr[1]);
   divStr[2] = toupper(divStr[2]);
   divStr[3] = toupper(divStr[3]);

   int ndiv = CodeWriter::createDivision(divStr,ref);

   XString motherS(ref.fMother->getAttribute(X("name")));
   std::cout
        << std::endl
        << "      chname = \'" << divStr << "\'" << std::endl
        << "      chmoth = \'" << S(motherS) << "\'" << std::endl
        << "      ndiv = " << ref.fPartition.ncopy << std::endl
        << "      iaxis = " << ref.fPartition.iaxis << std::endl
        << "      step = " << ref.fPartition.step << std::endl
        << "      c0 = " << ref.fPartition.start << std::endl
        << "      numed = 0" << std::endl
        << "      ndvmax = 0" << std::endl
        << "      call gsdvx(chname,chmoth,ndiv,iaxis,step,c0,numed,ndvmax)"
        << std::endl;

   return ndiv;
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

int FortranWriter::createVolume(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createVolume: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int icopy = CodeWriter::createVolume(el,ref);

   if (fPending)
   {
      XString nameS(el->getAttribute(X("name")));
      XString motherS(fRef.fMother->getAttribute(X("name")));
      int irot = fRef.fRotation;
      std::cout
           << std::endl
           << "      chname = \'" << S(nameS) << "\'" << std::endl
           << "      nr = " << icopy << std::endl
           << "      chmoth = \'" << S(motherS) << "\'" << std::endl
           << "      x = " << fRef.fOrigin[0] << std::endl
           << "      y = " << fRef.fOrigin[1] << std::endl
           << "      z = " << fRef.fOrigin[2] << std::endl
           << "      irot = " << irot << std::endl
           << "      chonly = \'ONLY\'" << std::endl
           << "      call gspos(chname,nr,chmoth,x,y,z,irot,chonly)"
           << std::endl;
      fPending = false;
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
   return icopy;
}

void FortranWriter::createHeader()
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createHeader: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createHeader();

   std::cout
        << "*"                                                    << std::endl
        << "* HDDSgeant3 - fortran geometry definition package"   << std::endl
        << "*              for the Hall D experiment."            << std::endl
        << "*"                                                    << std::endl
        << "*         WARNING: DO NOT EDIT THIS FILE"             << std::endl
        << "*"                                                    << std::endl
        << "* This file was generated automatically from the"     << std::endl
        << "* HDDS xml geometry definition by the hdds-geant"     << std::endl
        << "* translator.  Any changes made to this file will"    << std::endl
        << "* disappear as soon as it is regenerated from the"    << std::endl
        << "* xml source.  To introduce Geant3 optimizations,"    << std::endl
        << "* see the subroutine Goptimize() in goptimize.F."     << std::endl
        << "*"                                                    << std::endl
        << "      subroutine HDDSgeant3"                          << std::endl
        << "      implicit none"                                  << std::endl
        << "      integer imate"                                  << std::endl
        << "      character*20 chnama,namate"                     << std::endl
        << "      real a,z,dens,radl,absl,ubuf(99)"               << std::endl
        << "      integer nwbuf"                                  << std::endl
        << "      real amat(99),zmat(99),wmat(99)"                << std::endl
        << "      integer nlmat"                                  << std::endl
        << "      integer itmed"                                  << std::endl
        << "      character*20 natmed"                            << std::endl
        << "      integer nmat,isvol,ifield"                      << std::endl
        << "      real fieldm,tmaxfd,stemax,deemax,epsil,stmin"   << std::endl
        << "      character*4 chname,chshap,chmoth"               << std::endl
        << "      integer nmed,npar,ivolu"                        << std::endl
        << "      real par(99)"                                   << std::endl
        << "      integer irot"                                   << std::endl
        << "      real theta1,phi1,theta2,phi2,theta3,phi3"       << std::endl
        << "      integer nr,ndiv,iaxis,numed,ndvmax"             << std::endl
        << "      real step,c0"                                   << std::endl
        << "      real x,y"                                       << std::endl
        << "      character*4 chonly"                             << std::endl;
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

void FortranWriter::createTrailer()
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createTrailer: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createTrailer();

   std::cout << "      end"                                       << std::endl;
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

void FortranWriter::createSetFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createSetFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createSetFunctions(el,ident);

   DOMNodeList* specL = el->getElementsByTagName(X("specify"));
   int len = specL->getLength();
   XString subNameStr;
   subNameStr = "setoptical" + ident;
   std::cout
        << std::endl
        << "      subroutine " << subNameStr << "(itmed)" << std::endl
        << "      implicit none" << std::endl
        << "      integer itmed" << std::endl
        << "      integer nbins" << std::endl
        << "      real Ephot(" << len << ")" << std::endl
        << "      real abslen(" << len << ")" << std::endl
        << "      real effic(" << len << ")" << std::endl
        << "      real rindex(" << len << ")" << std::endl
        << "      real smooth(" << len << ")" << std::endl
        << "      real refloss(" << len << ")" << std::endl
        << "      common /optical" << ident
        << "/nbins,Ephot,rindex,abslen,refloss,smooth,effic" << std::endl
        << "      integer i" << std::endl
        << "      data nbins/" << len << "/" << std::endl;

   std::vector<double> Ephot;
   std::vector<double> abslen;
   std::vector<double> effic;
   std::vector<double> rindex;
   std::vector<double> smooth;
   std::vector<double> reflect;
   for (int i=0; i < len; ++i)
   {
      DOMElement* specEl = (DOMElement*)specL->item(i);
      XString valS;
      valS = specEl->getAttribute(X("E"));
      Ephot.push_back(atof(S(valS)));
      valS = specEl->getAttribute(X("refindex"));
      rindex.push_back(atof(S(valS)));
      valS = specEl->getAttribute(X("abslen"));
      abslen.push_back(atof(S(valS)));
      valS = specEl->getAttribute(X("smooth"));
      smooth.push_back(atof(S(valS)));
      valS = specEl->getAttribute(X("reflect"));
      reflect.push_back(atof(S(valS)));
      valS = specEl->getAttribute(X("effic"));
      effic.push_back(atof(S(valS)));
   }
   std::stringstream Ephotstr;
   std::stringstream rindexstr;
   std::stringstream abslenstr;
   std::stringstream smoothstr;
   std::stringstream reflosstr;
   std::stringstream efficstr;
   for (int i=0; i < len; ++i)
   {
      if (i % 60 == 0)
      {
         int ilimit = i + 60;
         ilimit = (ilimit > len)? len : ilimit;
         Ephotstr  << "      data (Ephot(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         rindexstr << "      data (rindex(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         abslenstr << "      data (abslen(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         smoothstr << "      data (smooth(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         reflosstr << "      data (refloss(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
         efficstr  << "      data (effic(i),i=" << i + 1 << "," << ilimit
                   << ") /" << std::endl;
      }
      if (i % 6 == 0)
      {
         Ephotstr  << "     + ";
         rindexstr << "     + ";
         abslenstr << "     + ";
         smoothstr << "     + ";
         reflosstr << "     + ";
         efficstr  << "     + ";
      }
      Ephotstr  << Ephot[i]*1e-9;
      rindexstr << rindex[i];
      abslenstr << abslen[i];
      smoothstr << smooth[i];
      reflosstr << 1-reflect[i];
      efficstr  << effic[i];
      if ((i == len-1) || (i % 60 == 59))
      {
         Ephotstr  << "/" << std::endl;
         rindexstr << "/" << std::endl;
         abslenstr << "/" << std::endl;
         smoothstr << "/" << std::endl;
         reflosstr << "/" << std::endl;
         efficstr  << "/" << std::endl;
      }
      else if (i % 6 == 5)
      {
         Ephotstr  << "," << std::endl;
         rindexstr << "," << std::endl;
         abslenstr << "," << std::endl;
         smoothstr << "," << std::endl;
         reflosstr << "," << std::endl;
         efficstr  << "," << std::endl;
      }
      else
      {
         Ephotstr  << ",";
         rindexstr << ",";
         abslenstr << ",";
         smoothstr << ",";
         reflosstr << ",";
         efficstr  << ",";
      }
   }
   std::cout << Ephotstr.str()
             << rindexstr.str()
             << abslenstr.str()
             << smoothstr.str()
             << reflosstr.str()
             << efficstr.str()
             << "      call GSCKOV(itmed,nbins,Ephot,"
             << ((rindex[0] < 1)? "refloss" : "abslen")
             << ",effic,rindex)" << std::endl
             << "      end"
             << std::endl;

   subNameStr = "getoptical" + ident;
   std::cout
        << std::endl
        << "      subroutine " << subNameStr
        << "(E,refl,absl,rind,plsh,eff)" << std::endl
        << "      implicit none" << std::endl
        << "      real E,refl,absl,rind,plsh,eff" << std::endl
        << "      integer n" << std::endl
        << "      real x" << std::endl
        << "      integer nbins" << std::endl
        << "      real Ephot(" << len << ")" << std::endl
        << "      real abslen(" << len << ")" << std::endl
        << "      real effic(" << len << ")" << std::endl
        << "      real rindex(" << len << ")" << std::endl
        << "      real smooth(" << len << ")" << std::endl
        << "      real refloss(" << len << ")" << std::endl
        << "      common /optical" << ident
        << "/nbins,Ephot,rindex,abslen,refloss,smooth,effic" << std::endl
        << "      do n=1,nbins" << std::endl
        << "        if (E.lt.Ephot(n)) go to 10" << std::endl
        << "      enddo" << std::endl
        << "   10 continue" << std::endl
        << "      if (n.eq.1.or.n.gt.nbins) then" << std::endl
        << "        refl = 0" << std::endl
        << "        absl = 0" << std::endl
        << "        rind = 0" << std::endl
        << "        plsh = 1" << std::endl
        << "        eff = 0" << std::endl
        << "      else" << std::endl
        << "        x = (E-Ephot(n-1))/(Ephot(n)-Ephot(n-1))" << std::endl
        << "        refl = (1-x)*refloss(n-1)+x*refloss(n)" << std::endl
        << "        absl = (1-x)*abslen(n-1)+x*abslen(n)" << std::endl
        << "        rind = (1-x)*rindex(n-1)+x*rindex(n)" << std::endl
        << "        plsh = (1-x)*smooth(n-1)+x*smooth(n)" << std::endl
        << "        eff = (1-x)*effic(n-1)+x*effic(n)" << std::endl
        << "      endif" << std::endl
        << "      end" << std::endl;

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

void FortranWriter::createGetFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createGetFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createGetFunctions(el,ident);

   std::vector<int> table;
   std::vector<int> start;
   start.push_back(0);

   XString funcNameStr;
   XString identCaps(ident);
   identCaps[0] = toupper(identCaps[0]);
   funcNameStr = "get" + identCaps;
   for (int ivolu = 1; ivolu <= Refsys::fVolumes; ivolu++)
   {
      start.push_back(0);
      int ncopy = Refsys::fIdentifierTable[ivolu]["copy counter"].back();
      std::map<std::string,std::vector<int> >::iterator idlist = 
                  Refsys::fIdentifierTable[ivolu].find(ident);
      if (idlist != Refsys::fIdentifierTable[ivolu].end())
      {
         if (ncopy != idlist->second.size())
         {
            std::cerr
                  << APP_NAME << " warning: volume " << ivolu
                  << " has " << ncopy << " copies, but "
                  << idlist->second.size() << " " 
                  << ident << " identifiers!" << std::endl;
            for (int idx = 0; idx < idlist->second.size(); idx++)
            {
               std::cerr << idlist->second[idx]  << " ";
               if (idx/20*20 == idx)
                  std::cerr << std::endl;
            }
            std::cerr << std::endl;
         }
         start[ivolu] = table.size() + 1;
         for (int icopy = 0; icopy < ncopy; icopy++)
         {
            table.push_back(idlist->second[icopy]);
         }
      }
      else
      {
         start[ivolu] = 0;
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
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

void FortranWriter::createMapFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createMapFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createMapFunctions(el,ident);

   if (el == 0)
   {
      return;
   }
   DOMNodeList* regionsL = el->getOwnerDocument()->getDocumentElement()
                             ->getElementsByTagName(X("regions"));
   if (regionsL->getLength() == 0)
   {
      return;
   }

   std::cout
        << std::endl
        << "      subroutine gufld(r,B)" << std::endl
        << "      implicit none" << std::endl
        << "      real r(3),B(3)" << std::endl
        << "      real rr(3),rs(3),BB(3)" << std::endl
        << "      integer iregion" << std::endl
        << "      integer getMap" << std::endl
        << "      external getMap" << std::endl
        << std::endl;

   DOMElement* regionsEl = (DOMElement*)regionsL->item(0);
   DOMNodeList* regionL = regionsEl->getElementsByTagName(X("region"));
   for (int ireg=0; ireg < regionL->getLength(); ++ireg)
   {
      DOMElement* regionEl = (DOMElement*)regionL->item(ireg);
      DOMNodeList* noTagL = regionEl->getElementsByTagName(X("noBfield"));
      if (noTagL->getLength())
      {
         continue;
      }
      DOMNodeList* mapL = regionEl->getElementsByTagName(X("HDDSregion"));
      for (int imap=0; imap < mapL->getLength(); ++imap)
      {
         DOMElement* mapEl = (DOMElement*)mapL->item(imap);
         XString idS(mapEl->getAttribute(X("id")));
         XString origS(mapEl->getAttribute(X("origin")));
         XString rmatS(mapEl->getAttribute(X("Rmatrix")));
         int id = atoi(S(idS));
         double origin[3];
         std::stringstream listr(S(origS));
         listr >> origin[0] >> origin[1] >> origin[2];
         double Rmatrix[3][3];
         listr.clear(), listr.str(S(rmatS));
         listr >> Rmatrix[0][0] >> Rmatrix[0][1] >> Rmatrix[0][2]
               >> Rmatrix[1][0] >> Rmatrix[1][1] >> Rmatrix[1][2]
               >> Rmatrix[2][0] >> Rmatrix[2][1] >> Rmatrix[2][2];
         std::cout
           << "      real orig" << id << "(3),rmat" << id << "(3,3)"
           << std::endl
           << "      data orig" << id << "/"
           << origin[0] << "," << origin[1] << "," << origin[2] << "/"
           << std::endl
           << "      data rmat" << id << "/"
           << Rmatrix[0][0] << "," << Rmatrix[0][1] << "," << Rmatrix[0][2]
           << "," << std::endl << "     +          "
           << Rmatrix[1][0] << "," << Rmatrix[1][1] << "," << Rmatrix[1][2]
           << "," << std::endl << "     +          "
           << Rmatrix[2][0] << "," << Rmatrix[2][1] << "," << Rmatrix[2][2]
           << "/" << std::endl;
      }
   }

   std::cout
        << std::endl
        << "      iregion = getMap()" << std::endl
        << "      if (iregion.eq.0) then" << std::endl
        << "        B(1) = 0" << std::endl
        << "        B(2) = 0" << std::endl
        << "        B(3) = 0" << std::endl;

   std::list<DOMElement*> fieldMap;
   for (int ireg=0; ireg < regionL->getLength(); ++ireg)
   {
      DOMElement* regionEl = (DOMElement*)regionL->item(ireg);
      DOMNodeList* unifTagL = regionEl->getElementsByTagName(X("uniformBfield"));
      DOMNodeList* mapfTagL = regionEl->getElementsByTagName(X("mappedBfield"));
      DOMNodeList* mapL = regionEl->getElementsByTagName(X("HDDSregion"));
      for (int imap=0; imap < mapL->getLength(); ++imap)
      {
         DOMElement* mapEl = (DOMElement*)mapL->item(imap);
         XString idS(mapEl->getAttribute(X("id")));
         XString origS(mapEl->getAttribute(X("origin")));
         XString rmatS(mapEl->getAttribute(X("Rmatrix")));
         int id = atoi(S(idS));
         double origin[3];
         std::stringstream listr(S(origS));
         listr >> origin[0] >> origin[1] >> origin[2];
         double Rmatrix[3][3];
         listr.clear(), listr.str(S(rmatS));
         listr >> Rmatrix[0][0] >> Rmatrix[0][1] >> Rmatrix[0][2]
               >> Rmatrix[1][0] >> Rmatrix[1][1] >> Rmatrix[1][2]
               >> Rmatrix[2][0] >> Rmatrix[2][1] >> Rmatrix[2][2];
         if (unifTagL->getLength() > 0)
         {
            DOMElement* unifEl = (DOMElement*)unifTagL->item(0);
            XString bvecS(unifEl->getAttribute(X("Bx_By_Bz")));
            std::stringstream listr(S(bvecS));

            double b[3];
            listr >> b[0] >> b[1] >> b[2];
            Units unit;
            unit.getConversions(unifEl);
            b[0] *= unit.kG;
            b[1] *= unit.kG;
            b[2] *= unit.kG;

            std::cout
             << "      else if (iregion.eq." << id << ") then" << std::endl
             << "        B(1) = "
             << Rmatrix[0][0]*b[0] + Rmatrix[0][1]*b[1] + Rmatrix[0][2]*b[2]
             << std::endl
             << "        B(2) = "
             << Rmatrix[1][0]*b[0] + Rmatrix[1][1]*b[1] + Rmatrix[1][2]*b[2]
             << std::endl
             << "        B(3) = "
             << Rmatrix[2][0]*b[0] + Rmatrix[2][1]*b[1] + Rmatrix[2][2]*b[2]
             << std::endl;
         }
         else if (mapfTagL->getLength() > 0)
         {
            DOMElement* mapfEl = (DOMElement*)mapfTagL->item(0);
            fieldMap.push_back(mapfEl);
            int map = fieldMap.size();

            Units unit;
            unit.getConversions(mapfEl);

            std::cout
             << "      else if (iregion.eq." << id << ") then" << std::endl
             << "        rs(1) = r(1)-orig" << id << "(1)" << std::endl
             << "        rs(2) = r(2)-orig" << id << "(2)" << std::endl
             << "        rs(3) = r(3)-orig" << id << "(3)" << std::endl
             << "        rr(1) = rs(1)*rmat" << id << "(1,1)"
             <<                "+rs(2)*rmat" << id << "(1,2)"
             <<                "+rs(3)*rmat" << id << "(1,3)" << std::endl
             << "        rr(2) = rs(1)*rmat" << id << "(2,1)"
             <<                "+rs(2)*rmat" << id << "(2,2)"
             <<                "+rs(3)*rmat" << id << "(2,3)" << std::endl
             << "        rr(3) = rs(1)*rmat" << id << "(3,1)"
             <<                "+rs(2)*rmat" << id << "(3,2)"
	     <<                 "+rs(3)*rmat" << id << "(3,3)" << std::endl
             << "        call gufld" << map << "(rr,BB)"     << std::endl
             << "        B(1) = BB(1)*rmat" << id << "(1,1)"
             <<               "+BB(2)*rmat" << id << "(2,1)"
             <<               "+BB(3)*rmat" << id << "(3,1)" << std::endl
             << "        B(2) = BB(1)*rmat" << id << "(1,2)"
             <<               "+BB(2)*rmat" << id << "(2,2)"
             <<               "+BB(3)*rmat" << id << "(3,2)" << std::endl
             << "        B(3) = BB(1)*rmat" << id << "(1,3)"
             <<               "+BB(2)*rmat" << id << "(2,3)"
             <<               "+BB(3)*rmat" << id << "(3,3)" << std::endl;

             if (unit.kG != 1)
             {
                std::cout
                   << "        B(1) = B(1)*" << unit.kG << std::endl
                   << "        B(2) = B(2)*" << unit.kG << std::endl
                   << "        B(3) = B(3)*" << unit.kG << std::endl;
             }
         }
      }
   }
   std::cout
        << "      endif" << std::endl
        << "      end" << std::endl;

   int map = 1;
   std::list<DOMElement*>::iterator iter;
   for (iter = fieldMap.begin(); iter != fieldMap.end(); ++iter, ++map)
   {
      DOMElement* regionEl = (DOMElement*)(*iter)->getParentNode();
      XString nameS(regionEl->getAttribute(X("name")));

      std::cout
        << std::endl
        << "      subroutine gufld" << map << "(r,B)" << std::endl
        << "      implicit none" << std::endl
        << "      real r(3),B(3),Br(3)" << std::endl
        << "      real rho,phi,alpha" << std::endl
        << "      real u(3)" << std::endl
        << "      real twopi" << std::endl
        << "      parameter (twopi=6.28318530717959)" << std::endl
        << std::endl;

      Units unit;
      unit.getConversions(*iter);

      int axorder[] = {0,0,0,0};
      int axsamples[] = {0,0,0,0};

      XString gridtype;
      DOMNodeList* gridL = (*iter)->getElementsByTagName(X("grid"));
      int ngrid;
      for (ngrid = 0; ngrid < gridL->getLength(); ++ngrid)
      {
         int axsense[] = {1,1,1,1};
         double axlower[4], axupper[4];
         DOMElement* gridEl = (DOMElement*)gridL->item(ngrid);
         XString typeS(gridEl->getAttribute(X("type")));
         if (gridtype.size() > 0 && typeS != gridtype)
         {
            std::cerr
               << APP_NAME << " error: mappedBfield in region " << S(nameS)
               << " superimposes incompatible grid types." << std::endl;
            exit(1);
         }
         gridtype = typeS;

         DOMNodeList* samplesL = gridEl->getElementsByTagName(X("samples"));
         if (samplesL->getLength() != 3)
         {
            std::cerr
              << APP_NAME << " error: mappedBfield in region " << S(nameS)
              << " does not have samples for three axes." << std::endl;
            exit(1);
         }

         for (int iax = 1; iax <= 3; ++iax)
         {
            DOMElement* sampleEl = (DOMElement*)samplesL->item(iax-1);
            XString nS(sampleEl->getAttribute(X("n")));
            XString axisS(sampleEl->getAttribute(X("axis")));
            XString boundsS(sampleEl->getAttribute(X("bounds")));
            XString senseS(sampleEl->getAttribute(X("sense")));
            Units sunit;
            double bound[2];
            sunit.getConversions(sampleEl);
            std::stringstream listr(boundsS);
            listr >> bound[0] >> bound[1];
            int iaxis;
            if (gridtype == "cartesian")
            {
               if (axisS == "x" && 
                  (axorder[0] == 0 || axorder[0] == iax))
               {
                  iaxis = 1;
                  axorder[0] = iax;
                  bound[0] *= sunit.cm;
                  bound[1] *= sunit.cm;
               }
               else if (axisS == "y" && 
                       (axorder[1] == 0 || axorder[1] == iax))
               {
                  iaxis = 2;
                  axorder[1] = iax;
                  bound[0] *= sunit.cm;
                  bound[1] *= sunit.cm;
               }
               else if (axisS == "z" &&
                       (axorder[2] == 0 || axorder[2] == iax))
               {
                  iaxis = 3;
                  axorder[2] = iax;
                  bound[0] *= sunit.cm;
                  bound[1] *= sunit.cm;
               }
               else
               {
                  std::cerr
                  << APP_NAME << " error: grid in region " << S(nameS)
                  << " contains an incompatible set of samples." << std::endl;
                  exit(1);
               }
            }
            else if (gridtype == "cylindrical")
            {
               if (axisS == "r" && axorder[0] == 0)
               {
                  iaxis = 1;
                  axorder[0] = iax;
                  bound[0] *= sunit.cm;
                  bound[1] *= sunit.cm;
               }
               else if (axisS == "phi" && axorder[1] == 0)
               {
                  iaxis = 2;
                  axorder[1] = iax;
                  bound[0] *= sunit.rad;
                  bound[1] *= sunit.rad;
               }
               else if (axisS == "z" && axorder[2] == 0)
               {
                  iaxis = 3;
                  axorder[2] = iax;
                  bound[0] *= sunit.cm;
                  bound[1] *= sunit.cm;
               }
               else
               {
                  std::cerr
                  << APP_NAME << " error: grid in region " << S(nameS)
                  << " contains an incompatible set of samples." << std::endl;
                  exit(1);
               }
            }

            int n = atoi(S(nS));
            if (axsamples[iaxis] == 0 ||
                axsamples[iaxis] == n)
            {
               axsamples[iaxis] = n;
            }
            else
            {
               std::cerr
                 << APP_NAME << " error: mappedBfield in region " << S(nameS)
                 << " combines incompatible grid elements." << std::endl;
               exit(1);
            }

            axlower[iaxis] = bound[0];
            axupper[iaxis] = bound[1];
            if (senseS == "reverse")
            {
               axsense[iaxis] = -1;
            }
         }
         std::cout
              << "      real bound" << ngrid << "(3,2)" << std::endl
              << "      data bound" << ngrid << "/"
              << axlower[1] << "," << axlower[2] << "," << axlower[3] << ","
              << std::endl << "     +            "
              << axupper[1] << "," << axupper[2] << "," << axupper[3] << "/"
              << std::endl
              << "      integer reverse" << ngrid << "(3)" << std::endl
              << "      data reverse" << ngrid << "/"
              << axsense[1] << "," << axsense[2] << "," << axsense[3] << "/"
              << std::endl;
      }
      std::cout
           << "      real Bmap(3,"
           << axsamples[1] << "," << axsamples[2] << "," << axsamples[3]
           << ")" << std::endl
           << "      integer nsites(3)" << std::endl
           << "      data nsites/" 
           << axsamples[1] << "," << axsamples[2] << "," << axsamples[3]
           << "/" << std::endl
           << "      logical loaded" << std::endl
           << "      data loaded/.false./" << std::endl
           << "      save Bmap,nsites,loaded" << std::endl
           << "      integer i,i1,i2,i3" << std::endl
           << std::endl;

      XString mapS((*iter)->getAttribute(X("map")));
      XString encS((*iter)->getAttribute(X("encoding")));
      if (encS != "utf-8")
      {
         std::cerr
              << APP_NAME << " error: mappedBfield in region " << S(nameS)
              << " uses unsupported encoding " << encS << std::endl;
         exit(1);
      }
      else if (mapS.substr(0,7) != "file://")
      {
         std::cerr
              << APP_NAME << " error: mappedBfield in region " << S(nameS)
              << " uses unsupported map URL " << mapS << std::endl;
         exit(1);
      }
      mapS.erase(0,7);

      std::cout
           << "      if (.not.loaded) then" << std::endl
           << "        print*,'Reading B-field map from \"" << mapS << "\"'" << std::endl
           << "        open(unit=78,file='" << mapS << "',status='old',err=7)"
           << std::endl
           << "        read(unit=78,fmt=*,err=5,end=6)" << std::endl
           << "     +      ((((Bmap(i,i1,i2,i3),i=1,3)," << std::endl
           << "     +         i" << axorder[2] << "=1," 
           << axsamples[axorder[2]] << ")," << std::endl
           << "     +        i" << axorder[1] << "=1," 
           << axsamples[axorder[1]] << ")," << std::endl
           << "     +       i" << axorder[0] << "=1," 
           << axsamples[axorder[0]] << ")" << std::endl
           << "        go to 8" << std::endl
           << "    5   stop 'error reading magnetic field map, stop'"
           << std::endl
           << "    6   stop 'EOF encountered reading magnetic field map, stop'"
           << std::endl
           << "    7   stop 'error opening magnetic field map, stop'"
           << std::endl
           << "    8   loaded=.true." << std::endl
           << "      endif" << std::endl
           << std::endl;

      for (int igrid = 0; igrid < ngrid; igrid++)
      {
         if (gridtype == "cylindrical")
         {
            std::cout
              << "      rho = sqrt(r(1)**2+r(2)**2)" << std::endl
              << "      phi = atan2(r(2),r(1))" << std::endl
              << "      u(1) = (rho-bound" << igrid << "(1,1))/"
              << "(bound" << igrid << "(1,2)" << "-bound" << igrid << "(1,1))"
              << std::endl
              << "      u(2) = (phi-bound" << igrid << "(2,1))/"
              << "(bound" << igrid << "(2,2)" << "-bound" << igrid << "(2,1))"
              << std::endl
              << "      alpha = abs(twopi/(bound" << igrid << "(2,2)"
              << "-bound" << igrid << "(2,1)))" << std::endl
              << "      u(2) = u(2)-int(u(2)/alpha)*alpha" << std::endl
              << "      if (u(2).lt.0) then" <<std::endl
              << "        u(2) = u(2)+alpha" << std::endl
              << "      endif" <<std::endl
              << "      u(3) = (r(3)-bound" << igrid << "(3,1))/"
              << "(bound" << igrid << "(3,2)" << "-bound" << igrid << "(3,1))"
              << std::endl
              << "      if ((u(1).ge.0.and.u(1).le.1).and." << std::endl
              << "     +    (u(2).ge.0.and.u(2).le.1).and." << std::endl
              << "     +    (u(3).ge.0.and.u(3).le.1)) then" << std::endl
              << "        call interpol3(Bmap,nsites,u,Br)" << std::endl
              << "        Br(1)=Br(1)*reverse" << igrid << "(1)" << std::endl
              << "        Br(2)=Br(2)*reverse" << igrid << "(2)" << std::endl
              << "        B(1)=Br(1)*cos(phi)-Br(2)*sin(phi)" << std::endl
              << "        B(2)=Br(2)*cos(phi)+Br(1)*sin(phi)" << std::endl
              << "        B(3)=Br(3)*reverse" << igrid << "(3)" << std::endl;
         }
         else
         {
            std::cout
              << "      u(1) = (r(1)-bound" << igrid << "(1,1))/"
              << "(bound" << igrid << "(1,2)" << "-bound" << igrid << "(1,1))"
              << std::endl
              << "      u(2) = (r(2)-bound" << igrid << "(2,1))/"
              << "(bound" << igrid << "(2,2)" << "-bound" << igrid << "(2,1))"
              << std::endl
              << "      u(3) = (r(3)-bound" << igrid << "(3,1))/"
              << "(bound" << igrid << "(3,2)" << "-bound" << igrid << "(3,1))"
              << std::endl
              << "      if ((u(1).ge.0.and.u(1).le.1).and." << std::endl
              << "     +    (u(2).ge.0.and.u(2).le.1).and." << std::endl
              << "     +    (u(3).ge.0.and.u(3).le.1)) then" << std::endl
              << "        call interpol3(Bmap,nsites,u,B)" << std::endl
              << "        B(1)=B(1)*reverse" << igrid << "(1)" << std::endl
              << "        B(2)=B(2)*reverse" << igrid << "(2)" << std::endl
              << "        B(3)=B(3)*reverse" << igrid << "(3)" << std::endl;
         }
         std::cout
              << "        return" << std::endl
              << "      endif" << std::endl
              << std::endl;
      }
      std::cout
           << "      B(1) = 0" << std::endl
           << "      B(2) = 0" << std::endl
           << "      B(3) = 0" << std::endl
           << "      return" << std::endl
           << "      end" << std::endl
           << std::endl
           << "      subroutine interpol3(Bmap,nsites,u,B)" << std::endl
           << "      implicit none" << std::endl
           << "      integer nsites(3)" << std::endl
           << "      real Bmap(3,nsites(1),nsites(2),nsites(3))" << std::endl
           << "      real u(3),B(3)" << std::endl
           << "      integer ir(3),ir0(3),ir1(3)" << std::endl
           << "      real ur(3),dur(3),ugrad(3,3)" << std::endl
           << "      integer i" << std::endl
           << "      do i=1,3" << std::endl
           << "        ur(i)=u(i)*(nsites(i)-1)+1" << std::endl
           << "        ir(i)=nint(ur(i))" << std::endl
           << "        ir0(i)=max(ir(i)-1,1)" << std::endl
           << "        ir1(i)=min(ir(i)+1,nsites(i))" << std::endl
           << "        dur(i)=(ur(i)-ir(i))/(ir1(i)-ir0(i)+1e-20)" << std::endl
           << "      enddo" << std::endl
           << "      ugrad(1,1)=(Bmap(1,ir1(1),ir(2),ir(3))"
           <<                  "-Bmap(1,ir0(1),ir(2),ir(3)))" << std::endl
           << "      ugrad(2,1)=(Bmap(2,ir1(1),ir(2),ir(3))"
           <<                  "-Bmap(2,ir0(1),ir(2),ir(3)))" << std::endl
           << "      ugrad(3,1)=(Bmap(3,ir1(1),ir(2),ir(3))"
           <<                  "-Bmap(3,ir0(1),ir(2),ir(3)))" << std::endl
           << "      ugrad(1,2)=(Bmap(1,ir(1),ir1(2),ir(3))"
           <<                  "-Bmap(1,ir(1),ir0(2),ir(3)))" << std::endl
           << "      ugrad(2,2)=(Bmap(2,ir(1),ir1(2),ir(3))"
           <<                  "-Bmap(2,ir(1),ir0(2),ir(3)))" << std::endl
           << "      ugrad(3,2)=(Bmap(3,ir(1),ir1(2),ir(3))"
           <<                  "-Bmap(3,ir(1),ir0(2),ir(3)))" << std::endl
           << "      ugrad(1,3)=(Bmap(1,ir(1),ir(2),ir1(3))"
           <<                  "-Bmap(1,ir(1),ir(2),ir0(3)))" << std::endl
           << "      ugrad(2,3)=(Bmap(2,ir(1),ir(2),ir1(3))"
           <<                  "-Bmap(2,ir(1),ir(2),ir0(3)))" << std::endl
           << "      ugrad(3,3)=(Bmap(3,ir(1),ir(2),ir1(3))"
           <<                  "-Bmap(3,ir(1),ir(2),ir0(3)))" << std::endl
           << "      B(1)=Bmap(1,ir(1),ir(2),ir(3))" << std::endl
           << "     +       +ugrad(1,1)*dur(1)+ugrad(1,2)*dur(2)"
           <<              "+ugrad(1,3)*dur(3)" << std::endl
           << "      B(2)=Bmap(2,ir(1),ir(2),ir(3))" << std::endl
           << "     +       +ugrad(2,1)*dur(1)+ugrad(2,2)*dur(2)"
           <<              "+ugrad(2,3)*dur(3)" << std::endl
           << "      B(3)=Bmap(3,ir(1),ir(2),ir(3))" << std::endl
           << "     +       +ugrad(3,1)*dur(1)+ugrad(3,2)*dur(2)"
           <<              "+ugrad(3,3)*dur(3)" << std::endl
           << "      end" << std::endl;
   }
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}

void FortranWriter::createUtilityFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createUtilityFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createUtilityFunctions(el, ident);

   std::cout
        << std::endl
        << "      subroutine getoptical"
        << "(imat,E,refl,absl,rind,plsh,eff)" << std::endl
        << "      implicit none" << std::endl
        << "      integer imat" << std::endl
        << "      real E,refl,absl,rind,plsh,eff" << std::endl
        << "      ";

   DOMNodeList* propL = el->getOwnerDocument()
                          ->getElementsByTagName(X("optical_properties"));
   int ifclauses = 0;
   for (int iprop=0; iprop < propL->getLength(); ++iprop)
   {
      DOMElement* propEl = (DOMElement*)propL->item(iprop);
      DOMElement* matEl = (DOMElement*)propEl->getParentNode();
      XString imateS(matEl->getAttribute(X("HDDSmate")));
      if (imateS.size() > 0)
      {
         int imate = atoi(S(imateS));
         std::cout
            << "if (imat.eq." << imate << ") then" << std::endl
            << "        call getoptical" << imate
            << "(E,refl,absl,rind,plsh,eff)" << std::endl
            << "      else";
         ++ifclauses;
      }
   }
   if (ifclauses)
   {
      std::cout
           << std::endl
           << "        refl = 0" << std::endl
           << "        absl = 0" << std::endl
           << "        rind = 0" << std::endl
           << "        plsh = 1" << std::endl
           << "        eff = 0" << std::endl
           << "      endif" << std::endl
           << "      end" << std::endl;
   }
   else
   {
      std::cout
           << std::endl
           << "      refl = 0" << std::endl
           << "      absl = 0" << std::endl
           << "      rind = 0" << std::endl
           << "      plsh = 1" << std::endl
           << "      eff = 0" << std::endl
           << "      end" << std::endl;
   }

   std::cout
        << std::endl
        << "      subroutine guplsh(medi0,medi1)" << std::endl
        << "      implicit none" << std::endl
        << "      integer medi0,medi1" << std::endl
        << "      real E,refl,absl,rind,plsh,eff" << std::endl
        << "      integer nmat" << std::endl
        << "      character*20 natmed" << std::endl
        << "      integer isvol,ifield" << std::endl
        << "      real fieldm,tmaxfd,stemax,deemax,epsil,stmin" << std::endl
        << "      integer nwbuf" << std::endl
        << "      real ubuf(99)" << std::endl
        << "      call GFTMED(medi1,"
        << "natmed,nmat,isvol,ifield,fieldm," << std::endl
        << "     +  tmaxfd,stemax,deemax,epsil,stmin,ubuf,nwbuf)" << std::endl
        << "      E = 2.5" << std::endl
        << "      call getoptical(nmat,E,refl,absl,rind,plsh,eff)" << std::endl
        << "      end" << std::endl;

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   std::cerr << timestr.str() << std::endl;
#endif
}
