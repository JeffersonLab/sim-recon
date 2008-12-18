/*  HDDS Common Classes
 *
 *  Author: richard.t.jones@uconn.edu
 *
 *  Original version - Richard Jones, January 6, 2006.
 *
 *  Notes:
 *  ------
 * 1. The HDDS specification is an xml document in the standard W3C
 *    schema namespace http://www.gluex.org/hdds, as described by the
 *    HDDS-1_0.xsd schema document.
 * 2. Access to the xml source is through the industry-standard DOM interface.
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 * 4. The common classes were created originally as tools for the hdds-geant
 *    converter to write Fortran geometry definition code for the Geant3
 *    simulation.  Later they were reused for other simulation packages,
 *    but occasionally the original purpose is still visible in comments
 *    and some of the implementation details.
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
 *  There is was originally no provision in the AGDD geometry model for
 *  magnetic field information.  This information was introduced in HDDS
 *  through a new concept called a "region" and a new tag "apply" which
 *  associates a volume in the geometry definition with a region.  For more
 *  information about the meaning of regions in HDDS and the methods for
 *  specifying a region's magnetic field, see the documentation in the
 *  HDDS schema file.  Geant needs to distinguish between 4 cases:
 *     ifield=0 : no magnetic field
 *     ifield=1 : general case of inhomogenous field (Runge-Kutta)
 *     ifield=2 : quasi-homogenous field with map (helical segments)
 *     ifield=3 : uniform field (helices along local z-axis)
 *  For detector regions with no magnetic field, the volume is created with
 *  ifield=0.  If a field map is provided for the region, it is assigned
 *  ifield=1 and the field interpolator is written as a part of the output
 *  code.  If the field is uniform then it is assigned ifield=2 unless it
 *  is along the z-axis, in which case it is assigned ifield=3.
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

#define APP_NAME "hddsCommon"

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

/*  Refsys class:
 *	Stores persistent information about coordinate
 *	reference systems which are used to orient and place
 *	volumes in the simulation geometry.
 */

int Refsys::fVolumes = 0;
int Refsys::fRegions = 0;
int Refsys::fRotations = 0;
std::map<std::string,Refsys::VolIdent> Refsys::fIdentifiers;
std::vector<std::map<std::string,std::vector<int> > > Refsys::fIdentifierTable;

Refsys::Refsys()			// empty constructor
 : fIdentifier(),
   fMother(0),
   fRegion(0),
   fRegionID(0),
   fPhiOffset(0)
{
   fMOrigin[0] = fMOrigin[1] = fMOrigin[2] = 0;
   fMRmatrix[0][0] = fMRmatrix[1][1] = fMRmatrix[2][2] = 1;
   fMRmatrix[0][1] = fMRmatrix[1][0] = fMRmatrix[1][2] =
   fMRmatrix[0][2] = fMRmatrix[2][0] = fMRmatrix[2][1] = 0;
   reset();
}

Refsys::Refsys(const Refsys& src)	// copy constructor
 : fIdentifier(src.fIdentifier),
   fMother(src.fMother),
   fRegion(src.fRegion),
   fRegionID(src.fRegionID),
   fPhiOffset(src.fPhiOffset),
   fPartition(src.fPartition)
{
   for (int i=0; i<3; i++)
   {
      fMOrigin[i] = src.fMOrigin[i];
      fMRmatrix[i][0] = src.fMRmatrix[i][0];
      fMRmatrix[i][1] = src.fMRmatrix[i][1];
      fMRmatrix[i][2] = src.fMRmatrix[i][2];
   }
   std::map<std::string,double>::const_iterator iter;
   for (iter = src.fPar.begin(); iter != src.fPar.end(); ++iter)
   {
      fPar[iter->first] = iter->second;
   }
   reset(src);
}

Refsys& Refsys::operator=(Refsys& src)	// copy operator (deep sematics)
{
   fIdentifier = src.fIdentifier;
   fMother = src.fMother;
   fRegion = src.fRegion;
   fRegionID = src.fRegionID;
   fPhiOffset = src.fPhiOffset;
   fPartition = src.fPartition;
   for (int i=0; i<3; i++)
   {
      fMOrigin[i] = src.fMOrigin[i];
      fMRmatrix[i][0] = src.fMRmatrix[i][0];
      fMRmatrix[i][1] = src.fMRmatrix[i][1];
      fMRmatrix[i][2] = src.fMRmatrix[i][2];
   }
   std::map<std::string,double>::iterator iter;
   for (iter = src.fPar.begin(); iter != src.fPar.end(); ++iter)
   {
      fPar[iter->first] = iter->second;
   }
   reset(src);
   return *this;
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
      fMOrigin[i] += fMRmatrix[i][0] * vector[0] +
                     fMRmatrix[i][1] * vector[1] +
                     fMRmatrix[i][2] * vector[2];
   }
   return *this;
}

Refsys& Refsys::shift(const Refsys& ref)      // copy origin from ref
{
   fOrigin[0] = ref.fOrigin[0];
   fOrigin[1] = ref.fOrigin[1];
   fOrigin[2] = ref.fOrigin[2];
   fMOrigin[0] = ref.fMOrigin[0];
   fMOrigin[1] = ref.fMOrigin[1];
   fMOrigin[2] = ref.fMOrigin[2];
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

         x[0] = fMRmatrix[i][0] * cosz + fMRmatrix[i][1] * sinz;
         x[1] = fMRmatrix[i][1] * cosz - fMRmatrix[i][0] * sinz;
         x[2] = fMRmatrix[i][2];
         xx[0] = x[0] * cosy - x[2] * siny;
         xx[1] = x[1];
         xx[2] = x[2] * cosy + x[0] * siny;
         fMRmatrix[i][0] = xx[0];
         fMRmatrix[i][1] = xx[1] * cosx + xx[2] * sinx;
         fMRmatrix[i][2] = xx[2] * cosx - xx[1] * sinx;
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
      fMRmatrix[i][0] = ref.fMRmatrix[i][0];
      fMRmatrix[i][1] = ref.fMRmatrix[i][1];
      fMRmatrix[i][2] = ref.fMRmatrix[i][2];
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

void Refsys::clearIdentifiers()
{
   fIdentifier.clear();
}

void Refsys::incrementIdentifiers()
{
   std::map<std::string,Refsys::VolIdent>::iterator iter;
   for (iter = fIdentifier.begin(); iter != fIdentifier.end(); ++iter)
   {
      iter->second.value += iter->second.step;
   }
}

void Refsys::addIdentifier(XString ident, int value, int step)
{
   VolIdent id;
   id.value = value;
   id.step = step;
   fIdentifier[ident] = id;
   fIdentifiers[ident] = id;
}

int Refsys::nextRotationID()
{
   return ++fRotations;
}

int Refsys::nextVolumeID()
{
   int ivolu = ++fVolumes;
   while (fIdentifierTable.size() <= ivolu)
   {
      std::map<std::string,std::vector<int> > unmarked;
      fIdentifierTable.push_back(unmarked);
   }
   return ivolu;
}

int Refsys::nextRegionID()
{
   return ++fRegions;
}

/* Substance class:
 *	Computes and saves properties of materials that are used
 *	in the detector description, sometimes using directly the
 *	values in the hdds file and other times computing them
 *	based on the hdds information.
 */

Substance::Substance()
 : fMaterialEl(0),
   fAtomicWeight(0), fAtomicNumber(0), fDensity(-1),
   fRadLen(0), fAbsLen(0), fColLen(0), fMIdEdx(0),
   fUniqueID(0),fBrewList(0)
{}

Substance::Substance(Substance& src)
 : fMaterialEl(src.fMaterialEl),
   fAtomicWeight(src.fAtomicWeight),
   fAtomicNumber(src.fAtomicNumber),
   fDensity(src.fDensity),
   fRadLen(src.fRadLen),
   fAbsLen(src.fAbsLen),
   fColLen(src.fColLen),
   fMIdEdx(src.fMIdEdx),
   fUniqueID(src.fUniqueID)
{
   std::list<Brew>::iterator iter;
   for (iter = src.fBrewList.begin();
        iter != src.fBrewList.end();
        ++iter)
   {
      Brew formula = *iter;
      formula.sub = new Substance(*iter->sub);
      fBrewList.push_back(formula);
   }
}

Substance::Substance(DOMElement* elem)
 : fMaterialEl(elem),
   fAtomicWeight(0), fAtomicNumber(0), fDensity(-1),
   fRadLen(0), fAbsLen(0), fColLen(0), fMIdEdx(0),
   fUniqueID(0),fBrewList(0)
{
   XString aS(fMaterialEl->getAttribute(X("a")));
   fAtomicWeight = atof(S(aS));
   XString zS(fMaterialEl->getAttribute(X("z")));
   fAtomicNumber = atof(S(zS));

   double wfactSum = 0;

   DOMNode* cont;
   for (cont = fMaterialEl->getFirstChild(); cont != 0;
        cont = cont->getNextSibling())
   {
      if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
      {
         DOMElement* contEl = (DOMElement*)cont;
         XString tagS(contEl->getTagName());
         if (tagS == "real")
         {
            Units unit;
            unit.getConversions(contEl);
            XString nameS(contEl->getAttribute(X("name")));
            XString valueS(contEl->getAttribute(X("value")));
            if (nameS == "density")
            {
               fDensity = atof(S(valueS)) * unit.g/unit.cm3;
            }
            else if (nameS == "radlen")
            {
               fRadLen = atof(S(valueS)) * unit.cm;
            }
            else if (nameS == "abslen")
            {
               fAbsLen = atof(S(valueS)) * unit.cm;
            }
            else if (nameS == "collen")
            {
               fColLen = atof(S(valueS)) * unit.cm;
            }
            else if (nameS == "dedx")
            {
               fMIdEdx = atof(S(valueS)) * unit.MeV/unit.cm;
            }
         }
         else if (tagS == "addmaterial")
         {
            XString matS(contEl->getAttribute(X("material")));
            DOMDocument* document = fMaterialEl->getOwnerDocument();
            DOMElement* targEl = document->getElementById(X(matS));
            XString typeS(targEl->getTagName());
            Substance::Brew formula;
            formula.sub = new Substance(targEl);
            formula.natoms = 0;
            formula.wfact = 0;
            DOMNode* mix;
            for (mix = contEl->getFirstChild(); mix != 0;
                 mix = mix->getNextSibling())
            {
               if (mix->getNodeType() == DOMNode::ELEMENT_NODE)
               {
                  DOMElement* mixEl = (DOMElement*)mix;
                  XString mixS(mixEl->getTagName());
                  if (mixS == "natoms")
                  {
                     if (typeS != "element")
                     {
                        std::cerr << APP_NAME
                             << " - error processing composite " 
                             << S(matS) << std::endl
                             << "natoms can only be specified for elements."
                             << std::endl;
                        exit(1);
                     }
                     XString nS(mixEl->getAttribute(X("n")));
                     formula.natoms = atoi(S(nS));
                     formula.wfact = formula.natoms * formula.sub->getAtomicWeight();
                  }
                  else if (mixS == "fractionmass")
                  {
                     XString fS(mixEl->getAttribute(X("fraction")));
                     formula.wfact = atof(S(fS));
                  }
               }
            }
            fBrewList.push_back(formula);
            wfactSum += formula.wfact;
         }
      }
   }

   double aSum=0;	double aNorm=0;
   double zSum=0;	double zNorm=0;
   double densitySum=0;	double densityNorm=0;
   double radlenSum=0;	double radlenNorm=0;
   double abslenSum=0;	double abslenNorm=0;
   double collenSum=0;	double collenNorm=0;
   double dedxSum=0;	double dedxNorm=0;

   std::list<Substance::Brew>::iterator iter;
   for (iter = fBrewList.begin(); iter != fBrewList.end(); ++iter)
   {
      iter->wfact /= wfactSum;
      double weight = (iter->natoms > 0)?
                iter->natoms * iter->sub->fAtomicWeight : iter->wfact;
      aSum += weight*iter->sub->fAtomicWeight;
      aNorm += weight;
      zSum += weight*iter->sub->fAtomicNumber;
      zNorm += weight;
      densitySum += weight/iter->sub->fDensity;
      densityNorm += weight;
      weight /= iter->sub->fDensity;
      radlenSum += weight/iter->sub->fRadLen;
      radlenNorm += weight;
      abslenSum += weight/iter->sub->fAbsLen;
      abslenNorm += weight;
      collenSum += weight/iter->sub->fColLen;
      collenNorm += weight;
      dedxSum += weight*iter->sub->fMIdEdx;
      dedxNorm += weight;
   }

   if (fAtomicWeight == 0 && aNorm > 0)
   {
      fAtomicWeight = aSum/aNorm;
   }

   if (fAtomicNumber == 0 && zNorm > 0)
   {
      fAtomicNumber = zSum/zNorm;
   }

   if (fDensity <= 0 && densitySum > 0)
   {
      fDensity = densityNorm/densitySum;
   }

   if (fRadLen == 0 && radlenSum > 0)
   {
      fRadLen = radlenNorm/radlenSum;
   }

   if (fAbsLen == 0 && abslenSum > 0)
   {
      fAbsLen = abslenNorm/abslenSum;
   }

   if (fColLen == 0 && collenSum > 0)
   {
      fColLen = collenNorm/collenSum;
   }

   if (fMIdEdx == 0 && dedxNorm > 0)
   {
      fMIdEdx = dedxSum/dedxNorm;
   }

   if (fDensity < 0)
   {
      XString tagS(fMaterialEl->getTagName());
      XString nameS(fMaterialEl->getAttribute(X("name")));
      std::cerr
           << APP_NAME << " error: " << S(tagS) << " " << S(nameS)
           << ", atomic number " << fAtomicNumber
           << ", atomic weight " << fAtomicWeight
           << " is missing a density specification." << std::endl;
      exit(1);
   }
}

Substance::~Substance()
{
   std::list<Brew>::iterator iter;
   for (iter = fBrewList.begin(); iter != fBrewList.end(); ++iter)
   {
      delete iter->sub;
   }
}

Substance& Substance::operator=(const Substance& src)
{
   fMaterialEl = src.fMaterialEl;
   fAtomicWeight = src.fAtomicWeight;
   fAtomicNumber = src.fAtomicNumber;
   fDensity = src.fDensity;
   fRadLen = src.fRadLen;
   fAbsLen = src.fAbsLen;
   fColLen = src.fColLen;
   fMIdEdx = src.fMIdEdx;
   fUniqueID = src.fUniqueID;
   fBrewList = src.fBrewList;
   std::list<Brew>::iterator iter;
   for (iter = fBrewList.begin();
        iter != fBrewList.end();
        ++iter)
   {
      iter->sub = new Substance(*iter->sub);
   }
}

double Substance::getAtomicWeight()
{
   return fAtomicWeight;
}

double Substance::getAtomicNumber()
{
   return fAtomicNumber;
}

double Substance::getDensity()
{
   return fDensity;
}

double Substance::getRadLength()
{
   return fRadLen;
}

double Substance::getAbsLength()
{
   return fAbsLen;
}

double Substance::getColLength()
{
   return fColLen;
}

double Substance::getMIdEdx()
{
   return fMIdEdx;
}

XString Substance::getName()
{
   return XString(fMaterialEl->getAttribute(X("name")));
}

XString Substance::getSymbol()
{
   return XString(fMaterialEl->getAttribute(X("symbol")));
}

DOMElement* Substance::getDOMElement()
{
   return fMaterialEl;
}


/* Units class:
 *	Provides conversion constants for convenient extraction of
 *	dimensioned parameters in the hdds file.
 */

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
   XString unitlS(el->getAttribute(X("unit_length")));
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
   else if (unitlS == "km")
   {
      set_1cm(1e-5);
   }
   else if (unitlS == "um")
   {
      set_1cm(1e6);
   }
   else if (unitlS == "nm")
   {
      set_1cm(1e9);
   }
   else if (unitlS == "in")
   {
      set_1cm(1/2.54);
   }
   else if (unitlS == "ft")
   {
      set_1cm(1/(12*2.54));
   }
   else if (unitlS == "mil")
   {
      set_1cm(1000/2.54);
   }
   else
   {
      XString tagS(el->getTagName());
      std::cerr
           << APP_NAME << " error: unknown length unit " << S(unitlS)
           << " on tag " << S(tagS) << std::endl;
      exit(1);
   }

   XString unitaS = el->getAttribute(X("unit_angle"));
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
           << APP_NAME << " error: unknown angle unit " << S(unitaS)
           << " on volume " << S(tagS) << std::endl;
      exit(1);
   }

   XString unitS = el->getAttribute(X("unit"));
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
   else if (unitlS == "km")
   {
      set_1cm(1e-5);
   }
   else if (unitlS == "um")
   {
      set_1cm(1e6);
   }
   else if (unitlS == "nm")
   {
      set_1cm(1e9);
   }
   else if (unitlS == "in")
   {
      set_1cm(1/2.54);
   }
   else if (unitlS == "ft")
   {
      set_1cm(1/(12*2.54));
   }
   else if (unitlS == "mil")
   {
      set_1cm(1000/2.54);
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
   else if (unitS == "Tesla" || unitS == "T")
   {
      set_1G(1e-4);
   }
   else if (unitS == "kG" || unitS == "kGs")
   {
      set_1G(1e-3);
   }
   else if (unitS == "G" || unitS == "Gs")
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
           << APP_NAME << " error: unknown unit " << S(unitS)
           << " on volume " << S(tagS) << std::endl;
      exit(1);
   }
}

#ifdef LINUX_CPUTIME_PROFILING
CPUtimer::CPUtimer()
{
   getRusage();
   gettimeofday(&fClock0,&fTZ);
   fClockRef = fClock0;
   fRef = fLast;
}

double CPUtimer::getUserTime()
{
   if (getRusage() == 0)
   {
      return fLast.ru_utime.tv_sec + fLast.ru_utime.tv_usec/1e6;
   }
   else
   {
      return -1;
   }
}

double CPUtimer::getSystemTime()
{
   if (getRusage() == 0)
   {
      return fLast.ru_stime.tv_sec + fLast.ru_stime.tv_usec/1e6;
   }
   else
   {
      return -1;
   }
}

double CPUtimer::getRealTime()
{
   if (getRusage() == 0)
   {
      return (fClock.tv_sec - fClock0.tv_sec)
            +(fClock.tv_usec - fClock0.tv_usec)/1e6;
   }
   else
   {
      return -1;
   }
}

double CPUtimer::getUserDelta()
{
   if (getRusage() == 0)
   {
      return fLast.ru_utime.tv_sec + fLast.ru_utime.tv_usec/1e6
            -(fRef.ru_utime.tv_sec + fRef.ru_utime.tv_usec/1e6);
   }
   else
   {
      return -1;
   }
}

double CPUtimer::getSystemDelta()
{
   if (getRusage() == 0)
   {
      return fLast.ru_stime.tv_sec + fLast.ru_stime.tv_usec/1e6
            -(fRef.ru_stime.tv_sec + fRef.ru_stime.tv_usec/1e6);
   }
   else
   {
      return -1;
   }
}

double CPUtimer::getRealDelta()
{
   if (getRusage() == 0)
   {
      return (fClock.tv_sec - fClockRef.tv_sec)
            +(fClock.tv_usec - fClockRef.tv_usec)/1e6;
   }
   else
   {
      return -1;
   }
}

void CPUtimer::resetClocks()
{
   getRusage();
   fRef = fLast;
   fClockRef = fClock;
}

int CPUtimer::getRusage()
{
   return getrusage(RUSAGE_SELF,&fLast) + gettimeofday(&fClock,&fTZ);
}

CPUtimer timer;
#endif

int CodeWriter::createMaterial(DOMElement* el)
{
   static int imateCount = 0;
   int imate = ++imateCount;
   std::stringstream imateStr;
   imateStr << imate;
   el->setAttribute(X("HDDSmate"),X(imateStr.str()));

   Substance subst(el);
   std::list<Substance::Brew>::iterator iter;
   for (iter = subst.fBrewList.begin();
        iter != subst.fBrewList.end(); ++iter)
   {
      DOMElement* subEl = iter->sub->getDOMElement();
      XString subS(subEl->getAttribute(X("HDDSmate")));
      iter->sub->fUniqueID = atoi(S(subS));
      if (iter->sub->fUniqueID == 0)
      {
         iter->sub->fUniqueID = createMaterial(subEl);
      }
   }
   fSubst = subst;
   return imate;
}

int CodeWriter::createSolid(DOMElement* el, Refsys& ref)
{
   XString nameS(el->getAttribute(X("name")));
   XString matS(el->getAttribute(X("material")));

   DOMDocument* document = el->getOwnerDocument();
   DOMElement* matEl = document->getElementById(X(matS));
   XString imateS(matEl->getAttribute(X("HDDSmate")));
   if (imateS.size() != 0)
   {
      fSubst.fUniqueID = atoi(S(imateS));
   }
   else
   {
      fSubst.fUniqueID = createMaterial(matEl);
   }
   
   int ivolu = ref.nextVolumeID();
   std::stringstream ivoluStr;
   ivoluStr << ivolu;
   el->setAttribute(X("HDDSvolu"),X(ivoluStr.str()));  
   el->setAttribute(X("HDDScopy"),X("0"));  

   return ivolu;
}

int CodeWriter::createRotation(Refsys& ref)
{
   if (ref.fRotation < 0)
   {
      ref.fRotation = ref.nextRotationID();
   }
   return ref.fRotation;
}

int CodeWriter::createRegion(DOMElement* el, Refsys& ref)
{
   int iregion = ref.nextRegionID();

   XString regionS(el->getAttribute(X("region")));
   DOMDocument* document = el->getOwnerDocument();
   ref.fRegion = document->getElementById(X(regionS));
   ref.fRegionID = iregion;

   double origin[3], angle[3];
   XString rotS(el->getAttribute(X("rot")));
   std::stringstream listr(rotS);
   listr >> angle[0] >> angle[1] >> angle[2];
   Units unit;
   unit.getConversions(el);
   angle[0] *= unit.rad;
   angle[1] *= unit.rad;
   angle[2] *= unit.rad;
   XString xyzS(el->getAttribute(X("origin")));
   listr.clear(), listr.str(xyzS);
   listr >> origin[0] >> origin[1] >> origin[2];
   origin[0] *= unit.cm;
   origin[1] *= unit.cm;
   origin[2] *= unit.cm;
   ref.shift(origin);
   ref.rotate(angle);

   std::stringstream attStr;
   attStr << iregion;
   el->setAttribute(X("HDDSregion"),X(attStr.str()));
   DOMElement* regEl = document->createElement(X("HDDSregion"));
   regEl->setAttribute(X("id"),X(attStr.str()));

   attStr.clear(), attStr.str("");
   attStr << ref.fMOrigin[0] << " "
          << ref.fMOrigin[1] << " " 
          << ref.fMOrigin[2];
   regEl->setAttribute(X("origin"),X(attStr.str()));
   attStr.clear(), attStr.str("");
   attStr << ref.fMRmatrix[0][0] << " "
          << ref.fMRmatrix[0][1] << " "
          << ref.fMRmatrix[0][2] << " "
          << ref.fMRmatrix[1][0] << " "
          << ref.fMRmatrix[1][1] << " "
          << ref.fMRmatrix[1][2] << " "
          << ref.fMRmatrix[2][0] << " "
          << ref.fMRmatrix[2][1] << " "
          << ref.fMRmatrix[2][2];
   regEl->setAttribute(X("Rmatrix"),X(attStr.str()));
   ref.fRegion->appendChild(regEl);

   for (DOMNode* cont = ref.fRegion->getFirstChild();
        cont != 0;
        cont = cont->getNextSibling())
   {
      if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
      {
         XString tagS(((DOMElement*)cont)->getTagName());
         if (tagS.find("Bfield") != XString::npos)
         {
            ref.addIdentifier(XString("map"),iregion,0);
            break;
         }
      }
   }
   return iregion;
}

int CodeWriter::createDivision(XString& divStr, Refsys& ref)
{
   int ncopy = ref.fPartition.ncopy;
   int ivolu = ref.nextVolumeID();

   assert (ref.fMother != 0);

   std::stringstream attStr;
   DOMDocument* document = ref.fMother->getOwnerDocument();
   DOMElement* divEl = document->createElement(X("HDDSdivision"));
   divEl->setAttribute(X("name"),X(divStr));
   XString motherS(ref.fMother->getAttribute(X("name")));
   divEl->setAttribute(X("volume"),X(motherS));
   attStr << ivolu;
   divEl->setAttribute(X("HDDSvolu"),X(attStr.str()));  
   std::stringstream copyStr;
   copyStr << ncopy;
   divEl->setAttribute(X("HDDScopy"),X(copyStr.str()));  
   ref.fMother->appendChild(divEl);
   ref.fPartition.divEl = divEl;

   std::map<std::string,Refsys::VolIdent>::iterator iter;
   for (iter = ref.fIdentifier.begin();
        iter != ref.fIdentifier.end();
        ++iter)
   {
      int value = iter->second.value;
      int step = iter->second.step;
      XString fieldS(iter->first);
      std::vector<int>* idlist = &Refsys::fIdentifierTable[ivolu][fieldS];
      for (int ic = 0; ic < ncopy; ic++)
      {
         idlist->push_back(value);
         value += step;
      }
   }
   Refsys::fIdentifierTable[ivolu]["copy counter"].push_back(ncopy);
   ref.clearIdentifiers();
   return ncopy;
}

int CodeWriter::createVolume(DOMElement* el, Refsys& ref)
{
   fPending = false;
   int icopy = 0;

   XString tagS(el->getTagName());
   XString nameS(el->getAttribute(X("name")));

   Refsys myRef(ref);
   DOMElement* env = 0;
   DOMDocument* document = el->getOwnerDocument();
   XString envS(el->getAttribute(X("envelope")));
   if (envS.size() != 0)
   {
      env = document->getElementById(X(envS));
      XString containS(env->getAttribute(X("contains")));
      if (containS == nameS)
      {
         return createVolume(env,myRef);
      }
      else if (containS.size() != 0)
      {
         std::cerr
              << APP_NAME << " error: re-use of shape " << S(envS)
              << " is not allowed by " << APP_NAME << std::endl;
         exit(1);
      }

      DOMNode* cont;
      for (cont = env->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() != DOMNode::ELEMENT_NODE)
         {
            continue;
         }
         DOMElement* contEl = (DOMElement*) cont;
         XString comdS(contEl->getTagName());
         if (comdS == "apply")
         {
            Refsys drs(myRef);
            myRef.fRegionID = createRegion(contEl,drs);
            myRef.fIdentifier["map"] = drs.fIdentifier["map"];
            myRef.fRegion = drs.fRegion;
            myRef.fPar = drs.fPar;
         }
      }

      env->setAttribute(X("contains"),X(nameS));
      icopy = createVolume(env,myRef);
      myRef.clearIdentifiers();
      myRef.fMother = env;
      myRef.reset();
   }

   if (tagS == "intersection" ||
       tagS == "subtraction" ||
       tagS == "union")
   {
      std::cerr
           << APP_NAME << " error: boolean " << S(tagS)
           << " operator is not supported by "<< APP_NAME << std::endl;
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
         XString targS(contEl->getAttribute(X("volume")));
         DOMElement* targEl = document->getElementById(X(targS));

         Refsys drs(myRef);
         double origin[3], angle[3];
         XString rotS(contEl->getAttribute(X("rot")));
         std::stringstream listr1(rotS);
         listr1 >> angle[0] >> angle[1] >> angle[2];
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
            XString fieldS(identEl->getAttribute(X("field")));
            XString valueS(identEl->getAttribute(X("value")));
            XString stepS(identEl->getAttribute(X("step")));
            drs.addIdentifier(fieldS,atoi(S(valueS)),atoi(S(stepS)));
         }

         if (comdS == "posXYZ")
         {
            XString xyzS(contEl->getAttribute(X("X_Y_Z")));
            std::stringstream listr(xyzS);
            listr >> origin[0] >> origin[1] >> origin[2];
            origin[0] *= unit.cm;
            origin[1] *= unit.cm;
            origin[2] *= unit.cm;
            drs.shift(origin);
            drs.rotate(angle);
            createVolume(targEl,drs);
         }
         else if (comdS == "posRPhiZ")
         {
            double r, phi, z;
            XString rphizS(contEl->getAttribute(X("R_Phi_Z")));
            std::stringstream listr(rphizS);
            listr >> r >> phi >> z;
            double s;
            XString sS(contEl->getAttribute(X("S")));
            s = atof(S(sS));
            phi *= unit.rad;
            r *= unit.cm;
            z *= unit.cm;
            s *= unit.cm;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            XString implrotS(contEl->getAttribute(X("impliedRot")));
            if (implrotS == "true" && (phi != 0))
            {
               angle[2] += phi;
            }
            drs.shift(origin);
            drs.rotate(angle);
            createVolume(targEl,drs);
         }
         else if (comdS == "mposPhi")
         {
            XString ncopyS(contEl->getAttribute(X("ncopy")));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << APP_NAME << " error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double phi0, dphi;
            XString phi0S(contEl->getAttribute(X("Phi0")));
            phi0 = atof(S(phi0S)) * unit.rad;
            XString dphiS(contEl->getAttribute(X("dPhi")));
            if (dphiS.size() != 0)
            {
               dphi = atof(S(dphiS)) * unit.rad;
            }
            else
            {
               dphi = 2 * M_PI / ncopy;
            }

            double r, s, z;
            XString rzS(contEl->getAttribute(X("R_Z")));
            std::stringstream listr(rzS);
            listr >> r >> z;
            XString sS(contEl->getAttribute(X("S")));
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
               containerS = el->getAttribute(X("divides"));
            }
            XString implrotS(contEl->getAttribute(X("impliedRot")));
            XString targTagS(targEl->getTagName());
            XString targEnvS(targEl->getAttribute(X("envelope")));
            DOMElement* targEnv;
            if (targEnvS.size() != 0)
            {
               targEnv = document->getElementById(X(targEnvS));
            }
            else
            {
               targEnv = targEl;
            }
            XString targProfS(targEnv->getAttribute(X("profile")));
            if (noRotation && (nSiblings == 1) &&
                (containerS == "pcon" ||
                 containerS == "cons" ||
                 containerS == "tubs") &&
                (targProfS.size() > 0) &&
                (implrotS == "true"))
            {
               double phi1, dphi1;
               std::stringstream profstr(targProfS);
               profstr >> phi1 >> dphi1;
               Units tunit;
               tunit.getConversions(targEnv);
               phi1 *= tunit.deg;
               dphi1 *= tunit.deg;
               static int phiDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "s" << std::setfill('0') << std::setw(3) << std::hex
                      << ++phiDivisions;
               phi0 *= unit.deg/unit.rad;
               dphi *= unit.deg/unit.rad;
               drs.fPartition.ncopy = ncopy;
               drs.fPartition.iaxis = 2;
               drs.fPartition.start = phi0 + phi1 - myRef.fPhiOffset;
               drs.fPartition.step = dphi;
               XString divS(divStr.str());
               createDivision(divS, drs);
               drs.fMother = drs.fPartition.divEl;
               targEl->setAttribute(X("divides"),X(containerS));
               drs.reset();

               double phioffset = phi1 + dphi/2;
               if (fabs(phioffset) > 0.001)
               {
                  angle[0] = 0;
                  angle[1] = 0;
                  angle[2] = -phioffset*unit.rad/unit.deg;
                  drs.rotate(angle);
               }
               origin[0] = r;
               origin[1] = s;
               origin[2] = z;
               drs.shift(origin);
               createVolume(targEl,drs);
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
                  createVolume(targEl,drs);
                  drs.incrementIdentifiers();
               }
            }
         }
         else if (comdS == "mposR")
         {
            XString ncopyS(contEl->getAttribute(X("ncopy")));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << APP_NAME << " error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double r0, dr;
            XString r0S(contEl->getAttribute(X("R0")));
            r0 = atof(S(r0S)) * unit.cm;
            XString drS(contEl->getAttribute(X("dR")));
            dr = atof(S(drS)) * unit.cm;

            double phi, z, s;
            XString zphiS(contEl->getAttribute(X("Z_Phi")));
            std::stringstream listr(zphiS);
            listr >> z >> phi;
            XString sS(contEl->getAttribute(X("S")));
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
               containerS = el->getAttribute(X("divides"));
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
               drs.fPartition.ncopy = ncopy;
               drs.fPartition.iaxis = 1;
               drs.fPartition.start = r0 - dr/2;
               drs.fPartition.step = dr;
               XString divS(divStr.str());
               createDivision(divS, drs);
               drs.fMother = drs.fPartition.divEl;
               targEl->setAttribute(X("divides"),X(containerS));
               origin[0] = r0 * cos(phi) - s * sin(phi);
               origin[1] = r0 * sin(phi) + s * cos(phi);
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               createVolume(targEl,drs);
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
                  createVolume(targEl,drs);
                  drs.incrementIdentifiers();
               }
            }
         }
         else if (comdS == "mposX")
         {
            XString ncopyS(contEl->getAttribute(X("ncopy")));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << APP_NAME << " error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double x0, dx;
            XString x0S(contEl->getAttribute(X("X0")));
            x0 = atof(S(x0S)) * unit.cm;
            XString dxS(contEl->getAttribute(X("dX")));
            dx = atof(S(dxS)) * unit.cm;

            double y, z, s;
            XString yzS(contEl->getAttribute(X("Y_Z")));
            std::stringstream listr(yzS);
            listr >> y >> z;
            XString sS(contEl->getAttribute(X("S")));
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
               containerS = el->getAttribute(X("divides"));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS == "box")
            {
               static int xDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "x" << std::setfill('0') << std::setw(3) << std::hex
                      << ++xDivisions;
               drs.fPartition.ncopy = ncopy;
               drs.fPartition.iaxis = 1;
               drs.fPartition.start = x0 - dx/2;
               drs.fPartition.step = dx;
               XString divS(divStr.str());
               createDivision(divS, drs);
               drs.fMother = drs.fPartition.divEl;
               targEl->setAttribute(X("divides"),X(containerS));
               origin[0] = 0;
               origin[1] = y + s;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               createVolume(targEl,drs);
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
                  createVolume(targEl,drs);
                  drs.incrementIdentifiers();
               }
            }
         }
         else if (comdS == "mposY")
         {
            XString ncopyS(contEl->getAttribute(X("ncopy")));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << APP_NAME << " error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double y0, dy;
            XString y0S(contEl->getAttribute(X("Y0")));
            y0 = atof(S(y0S)) * unit.cm;
            XString dyS(contEl->getAttribute(X("dY")));
            dy = atof(S(dyS)) * unit.cm;

            double x, z, s;
            XString zxS(contEl->getAttribute(X("Z_X")));
            std::stringstream listr(zxS);
            listr >> z >> x;
            XString sS(contEl->getAttribute(X("S")));
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
               containerS = el->getAttribute(X("divides"));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS == "box")
            {
               static int yDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "y" << std::setfill('0') << std::setw(3) << std::hex 
                      << ++yDivisions;
               drs.fPartition.ncopy = ncopy;
               drs.fPartition.iaxis = 2;
               drs.fPartition.start = y0 - dy/2;
               drs.fPartition.step = dy;
               XString divS(divStr.str());
               createDivision(divS, drs);
               drs.fMother = drs.fPartition.divEl;
               targEl->setAttribute(X("divides"),X(containerS));
               origin[0] = x + s;
               origin[1] = 0;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               createVolume(targEl,drs);
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
                  createVolume(targEl,drs);
                  drs.incrementIdentifiers();
               }
            }
         }
         else if (comdS == "mposZ")
         {
            XString ncopyS(contEl->getAttribute(X("ncopy")));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               std::cerr
                    << APP_NAME << " error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!"
                    << std::endl;
               exit(1);
            }

            double z0, dz;
            XString z0S(contEl->getAttribute(X("Z0")));
            z0 = atof(S(z0S)) * unit.cm;
            XString dzS(contEl->getAttribute(X("dZ")));
            dz = atof(S(dzS)) * unit.cm;

            double x, y, s;
            XString xyS(contEl->getAttribute(X("X_Y")));
            if (xyS.size() > 0)
            {
               std::stringstream listr(xyS);
               listr >> x >> y;
            }
            else
            {
               double r, phi;
               XString rphiS(contEl->getAttribute(X("R_Phi")));
               std::stringstream listr(rphiS);
               listr >> r >> phi;
               phi *= unit.rad;
               x = r * cos(phi);
               y = r * sin(phi);
            }
            XString sS(contEl->getAttribute(X("S")));
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
               containerS = el->getAttribute(X("divides"));
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS.size() != 0))
            {
               static int zDivisions = 0xd00;
               std::stringstream divStr;
               divStr << "z" << std::setfill('0') << std::setw(3) << std::hex
                      << ++zDivisions;
               drs.fPartition.ncopy = ncopy;
               drs.fPartition.iaxis = 3;
               drs.fPartition.start = z0 - dz/2;
               drs.fPartition.step = dz;
               XString divS(divStr.str());
               createDivision(divS, drs);
               drs.fMother = drs.fPartition.divEl;
               targEl->setAttribute(X("divides"),X(containerS));
               double phi = atan2(y,x);
               origin[0] = x - s * sin(phi);
               origin[1] = y + s * cos(phi);
               origin[2] = 0;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               createVolume(targEl,drs);
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
                  createVolume(targEl,drs);
                  drs.incrementIdentifiers();
               }
            }
         }
         else if (comdS == "apply")
         {
            myRef.fRegionID = createRegion(contEl,drs);
            myRef.fIdentifier["map"] = drs.fIdentifier["map"];
            myRef.fRegion = drs.fRegion;
            myRef.fPar = drs.fPar;
         }
         else
         {
            std::cerr
                 << APP_NAME << " error: composition of volume " << S(nameS)
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
           << APP_NAME << " error: stacks are not supported by " << APP_NAME
           << std::endl
           << "Use compositions instead." << std::endl;
      exit(1);
   }
   else
   {
      XString icopyS(el->getAttribute(X("HDDScopy")));
      if (icopyS.size() != 0)
      {
         icopy = atoi(S(icopyS));
      }
      else
      {
         XString profS(el->getAttribute(X("profile")));
         if (profS.size() != 0)
         {
            double phi0, dphi;
            std::stringstream listr(profS);
            listr >> phi0 >> dphi;
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
            el->setAttribute(X("profile"),X(pStr.str()));
         }
         createSolid(el,myRef);
         icopy = 0;
      }

      if (myRef.fMother != 0)
      {
         createRotation(myRef);
         fPending = true;
         fRef = myRef;
         ++icopy;
      }

      std::stringstream icopyStr;
      icopyStr << icopy;
      el->setAttribute(X("HDDScopy"),X(icopyStr.str()));
      XString voluS(el->getAttribute(X("HDDSvolu")));
      int ivolu = atoi(S(voluS));
      std::map<std::string,Refsys::VolIdent>::iterator iter;
      for (iter = myRef.fIdentifier.begin();
           iter != myRef.fIdentifier.end();
           ++iter)
      {
         XString fieldS(iter->first);
         std::vector<int>* idlist = &Refsys::fIdentifierTable[ivolu][fieldS];
         while (idlist->size() < icopy-1)
         {
            idlist->push_back(0);
         }
         idlist->push_back(iter->second.value);
      }
      Refsys::fIdentifierTable[ivolu]["copy counter"].push_back(icopy);
   }
   return icopy;
}

void CodeWriter::createHeader()
{
}

void CodeWriter::createTrailer()
{
}

void CodeWriter::createSetFunctions(DOMElement* el, const XString& ident)
{
}

void CodeWriter::createGetFunctions(DOMElement* el, const XString& ident)
{
}

void CodeWriter::createMapFunctions(DOMElement* el, const XString& ident)
{
}

void CodeWriter::createUtilityFunctions(DOMElement* el, const XString& ident)
{
}

void CodeWriter::translate(DOMElement* topel)
{
   Refsys mrs;
   createHeader();
   createVolume(topel,mrs);
   createTrailer();

   DOMNodeList* propL = topel->getOwnerDocument()
                             ->getElementsByTagName(X("optical_properties"));
   for (int iprop=0; iprop < propL->getLength(); ++iprop)
   {
      DOMElement* propEl = (DOMElement*)propL->item(iprop);
      DOMElement* matEl = (DOMElement*)propEl->getParentNode();
      XString imateS(matEl->getAttribute(X("HDDSmate")));
      if (imateS.size() > 0)
      {
         createSetFunctions(propEl,imateS);
      }
   }

   std::map<std::string,Refsys::VolIdent>::iterator iter;
   for (iter = mrs.fIdentifiers.begin();
        iter != mrs.fIdentifiers.end(); ++iter)
   {
      createGetFunctions(topel, iter->first);
   }

   DOMNodeList* regionsL = topel->getOwnerDocument()
                                ->getElementsByTagName(X("regions"));
   createMapFunctions((DOMElement*)regionsL->item(0),XString("map"));

   createUtilityFunctions(topel,XString("user"));
}

void CodeWriter::dump(DOMElement* el, int level=0) // useful debug function
{
   XString tagS(el->getTagName());
   for (int i=0; i<level; i++)
   {
      std::cerr << "  ";
   }
   std::cerr << "<" << tagS;
   DOMNamedNodeMap* attribL = el->getAttributes();
   for (int i=0; i < attribL->getLength(); i++)
   {
      XString nameS(attribL->item(i)->getNodeName());
      XString valueS(attribL->item(i)->getNodeValue());
      std::cerr << " " << nameS << "=\"" << valueS << "\"";
   }
   int n=0;
   for (DOMNode* cont = el->getFirstChild(); 
        cont != 0;
        cont = cont->getNextSibling())
   {
      if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
      {
         ++n;
      }
   }
   if (n > 0)
   {
      std::cerr << ">" << std::endl;
      for (DOMNode* cont = el->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
         {
            dump((DOMElement*)cont,level+1);
         }
      }
      for (int i=0; i<level; i++)
      {
         std::cerr << "  ";
      }
      std::cerr << "</" << tagS << ">" << std::endl;
   }
   else
   {
      std::cerr << "/>" << std::endl;
   }
}
