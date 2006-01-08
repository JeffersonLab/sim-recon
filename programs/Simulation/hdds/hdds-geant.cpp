/*
 *  hdds-geant :   an interface utility that reads in a HDDS document
 *                   (Hall D Detector Specification) and writes out a
 *                   GEANT-3 geometry description in the form of a
 *                   fortran subroutine.
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

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()

double fieldStrength[] =
{
   0.0,          // zero field regions
   0.0,   // inhomogenous field regions (unused)
   22.4,  // mapped field regions: solenoid (kG, approximate)
   2.0    // uniform field regions: sweep magnets (kG)
};

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
   int createMaterial(DOMElement* el);    // generate code for materials
   int createSolid(DOMElement* el,
                   const Refsys& ref);    // generate code for solids
   int createRotation(Refsys& ref);       // generate code for rotations
   int createVolume(DOMElement* el,
                    const Refsys& ref);   // generate code for placement
   int createDivision(XString& divStr,
                      Refsys& ref);	  // generate code for divisions
   void createGetFunctions(DOMElement* el,
                     XString& ident);     // generate code for identifiers
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
      XString msgS(e.msg);
      std::cerr << "Woops " << S(msgS) << std::endl;
      return 1;
   }

   XString everythingS("everything");
   DOMElement* rootEl = document->getElementById(X(everythingS));
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


int FortranWriter::createMaterial(DOMElement* el)
{
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
   return imate;
}

int FortranWriter::createSolid(DOMElement* el, const Refsys& ref)
{
   int ivolu = CodeWriter::createSolid(el,ref);
   int imate = fSubst.fUniqueID;

   static int itmedCount = 0;
   int itmed = ++itmedCount;
   XString nameAttS("name");
   XString nameS(el->getAttribute(X(nameAttS)));
   XString matAttS("material");
   XString matS(el->getAttribute(X(matAttS)));
   XString sensiAttS("sensitive");
   XString sensiS(el->getAttribute(X(sensiAttS)));
   std::cout
        << std::endl
        << "      itmed = " << itmed << std::endl
        << "      natmed = \'" << S(nameS) << " " << S(matS) << "\'"
        << std::endl
        << "      nmat = " << imate << std::endl
        << "      isvol = " << (sensiS == "true" ? 1 : 0) << std::endl
        << "      ifield = " << ((ref.fMagField == 0) ? 0 : 2) << std::endl
        << "      fieldm = " << fieldStrength[ref.fMagField] << std::endl
        << "      tmaxfd = " << ((ref.fMagField == 0) ? 0 : 1) << std::endl
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
      XString riozAttS("Rio_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      std::stringstream listr(riozS);
      listr >> ri >> ro >> zl;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
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
      XString xyzAttS("Xmp_Ymp_Z");
      XString xyzS(el->getAttribute(X(xyzAttS)));
      std::stringstream listr(xyzS);
      listr >> xm >> xp >> ym >> yp >> zl;
      double alph_xz, alph_yz;
      XString incAttS("inclination");
      XString incS(el->getAttribute(X(incAttS)));
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
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
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
      XString segAttS("segments");
      XString segS(el->getAttribute(X(segAttS)));
      segments = atoi(S(segS));
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
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
      XString riozAttS("Rio1_Rio2_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      std::stringstream listr(riozS);
      listr >> rim >> rom >> rip >> rop >> zl;
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
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
      XString rioAttS("Rio");
      XString rioS(el->getAttribute(X(rioAttS)));
      std::stringstream listr(rioS);
      listr >> ri >> ro;
      double theta0, theta1;
      XString polarAttS("polar_bounds");
      XString polarS(el->getAttribute(X(polarAttS)));
      listr.clear(), listr.str(polarS);
      listr >> theta0 >> theta1;
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
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

   return ivolu;
}

int FortranWriter::createRotation(Refsys& ref)
{
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
   return irot;
}

int FortranWriter::createDivision(XString& divStr, Refsys& ref)
{
   divStr[0] = toupper(divStr[0]);
   divStr[1] = toupper(divStr[1]);
   divStr[2] = toupper(divStr[2]);
   divStr[3] = toupper(divStr[3]);

   int ndiv = CodeWriter::createDivision(divStr,ref);

   XString nameAttS("name");
   XString motherS(ref.fMother->getAttribute(X(nameAttS)));
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
}

int FortranWriter::createVolume(DOMElement* el, const Refsys& ref)
{
   int icopy = CodeWriter::createVolume(el,ref);

   if (fPending)
   {
      XString nameAttS("name");
      XString nameS(el->getAttribute(X(nameAttS)));
      XString motherS(fRef.fMother->getAttribute(X(nameAttS)));
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

   return icopy;
}

void FortranWriter::createHeader()
{
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
}

void FortranWriter::createTrailer()
{
   std::cout << "      end"                                       << std::endl;
}

void FortranWriter::createGetFunctions(DOMElement* el, XString& ident)
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
