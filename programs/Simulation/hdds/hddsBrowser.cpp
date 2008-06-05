/*  HDDS Browser Classes
 *
 *  Author: richard.t.jones@uconn.edu
 *
 *  Original version - Richard Jones, June 3, 2008.
 *
 */

#include "XString.hpp"
#include "XParsers.hpp"
#include "hddsBrowser.hpp"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <list>

#define APP_NAME "hddsBrowser"

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()


// constructor
hddsBrowser::hddsBrowser(const XString xmlFile)
{
   fGeomDoc = buildDOMDocument(xmlFile,false);
   if (fGeomDoc == 0) {
      std::cerr
           << APP_NAME << " - error parsing HDDS document, "
           << "cannot continue" << std::endl;
      return;
   }
}

// Look up a volume in the geometry and return a pointer to
// a vector of Refsys objects containing the origin and rotation
// parameters for the placement of the object in the global
// reference system.  The user is responsible for deleting the
// returned vector when he is done with it.
std::vector<Refsys>* hddsBrowser::find(const XString volume,
                                       const Refsys *ref,
                                       const DOMElement *contEl)
{
   std::vector<Refsys> *result = new std::vector<Refsys>;
   if (contEl == 0) {
      contEl = fGeomDoc->getElementById(X("everything"));
      if (contEl == 0) {
         std::cerr
           << APP_NAME << " - error scanning HDDS document, " << std::endl
           << "  no element named \"everything\" found" << std::endl;
         return result;
      }
   }
   XString nameS = contEl->getAttribute(X("name"));
   XString envelS = contEl->getAttribute(X("envelope"));
   if (nameS == volume || envelS == volume) {
      if (ref == 0) {
         Refsys ref0;
         result->push_back(ref0);
      }
      else {
         result->push_back(*ref);
      }
      return result;
   }
#if DEBUGGING_FIND
   std::cout << "exploring volume " << nameS.c_str()
             << ", envelope " << envelS.c_str() << std::endl;
#endif

   // At this point contEl is pointing to an element in the geometry document
   // that is a container of positioning tags:
   //    * composition
   //    * union | intersection | subtraction
   //    * stackX | stackY | stackZ
   // The task of find() is to go through the list of positioning tags:
   //    * posXYZ (composition or booleans)
   //    * posRPhiZ (composition or booleans)
   //    * mposX (composition only)
   //    * mposY (composition only)
   //    * mposZ (composition only)
   //    * mposR (composition only)
   //    * mposPhi (composition only)
   //    * axisPos (stacks only)
   //    * axisMPos (stacks only)
   // and update the positioning parameters in ref, and look up the daughter
   // volume in the geometry tree.  There are three possibilities:
   //    1. The daughter volume name matches the argument "volume":
   //         result: Form a vector of Refsys objects describing
   //                 the placement of the named daughter in the MRS.
   //    2. else the daughter volume is another container element:
   //         result: Update the placement information in ref and call the
   //                 find() method recursively for each positioning step.
   //    3. else the daughter volume is one of the basic shapes:
   //         result: Return an empty vector of Refsys objects.

   for (DOMNode *node = contEl->getFirstChild(); 
        node != 0;
        node = node->getNextSibling())
   {
      if (node->getNodeType() != DOMNode::ELEMENT_NODE) {
         continue;
      }
      DOMElement* el = (DOMElement*) node;
      XString elS(el->getTagName());
      XString volumeS(el->getAttribute(X("volume")));
      double origin[] = {0, 0, 0};
      double angle[] = {0, 0, 0};
      int ncopy = 1;
      Units unit;
      unit.getConversions(el);
      XString rotS(el->getAttribute(X("rot")));
      if (rotS.size() > 0) {
         std::stringstream listr1(rotS);
         listr1 >> angle[0] >> angle[1] >> angle[2];
         angle[0] *= unit.rad;
         angle[1] *= unit.rad;
         angle[2] *= unit.rad;
      }
      XString implrotS(contEl->getAttribute(X("impliedRot")));
      XString sS(el->getAttribute(X("S")));
      double s = 0;
      if (sS.size() > 0) {
         s = atof(S(sS));
         s *= unit.cm;
      }
      XString ncopyS(el->getAttribute(X("ncopy")));
      if (ncopyS.size() > 0) {
         ncopy = atoi(S(ncopyS));
      }
      Refsys *dref = (ref == 0)? new Refsys : new Refsys(*ref);
      if (elS == "posXYZ") {
         XString xyzS(el->getAttribute(X("X_Y_Z")));
         if (xyzS.size() > 0) {
            std::stringstream listr1(xyzS);
            listr1 >> origin[0] >> origin[1] >> origin[2];
            origin[0] *= unit.cm;
            origin[1] *= unit.cm;
            origin[2] *= unit.cm;
            dref->shift(origin);
            dref->rotate(angle);
            if (volumeS == volume) {
               result->push_back(*dref);
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,dref,childEl);
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
               }
               delete res;
            }
         }
      }
      else if (elS == "posRPhiZ") {
         XString rphizS(el->getAttribute(X("R_Phi_Z")));
         if (rphizS.size() > 0) {
            double r=0, phi=0;
            if (rphizS.size() > 0) {
               std::stringstream listr1(rphizS);
               listr1 >> r >> phi >> origin[2];
               r *= unit.cm;
               phi *= unit.rad;
            }
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] *= unit.cm;
            if (implrotS == "true") {
               angle[2] += phi;
            }
            dref->shift(origin);
            dref->rotate(angle);
            if (volumeS == volume) {
               result->push_back(*dref);
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,dref,childEl);
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
               }
               delete res;
            }
         }
      }
      else if (elS == "mposPhi") {
         double phi0, dphi;
         XString phi0S(contEl->getAttribute(X("Phi0")));
         phi0 = atof(S(phi0S)) * unit.rad;
         XString dphiS(contEl->getAttribute(X("dPhi")));
         if (dphiS.size() != 0) {
            dphi = atof(S(dphiS)) * unit.rad;
         }
         else {
            dphi = 2 * M_PI / ncopy;
         }
         double r=0, z=0;
         XString rzS(contEl->getAttribute(X("R_Z")));
         if (rzS.size() > 0) {
            std::stringstream listr(rzS);
            listr >> r >> z;
            r *= unit.cm;
            z *= unit.cm;
         }
         dref->rotate(angle);
         for (int inst = 0; inst < ncopy; inst++) {
            double phi = phi0 + inst * dphi;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            Refsys *mref = new Refsys(*dref);
            mref->shift(origin);
            if (implrotS == "true") {
               angle[2] += ((inst == 0) ? phi0 : dphi);
               mref->rotate(angle);
            }
            if (volumeS == volume) {
               result->push_back(*mref);
               delete mref;
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,mref,childEl);
               delete mref;
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
                  delete res;
               }
               else {
                  delete res;
                  break;
               }
            }
         }
      }
      else if (elS == "mposR") {
         double r0, dr;
         XString r0S(contEl->getAttribute(X("R0")));
         r0 = atof(S(r0S)) * unit.cm;
         XString drS(contEl->getAttribute(X("dR")));
         dr = atof(S(drS)) * unit.cm;
         double phi=0, z=0;
         XString zphiS(contEl->getAttribute(X("Z_Phi")));
         if (zphiS.size() > 0) {
            std::stringstream listr(zphiS);
            listr >> z >> phi;
            phi *= unit.rad;
            z *= unit.cm;
         }
         dref->rotate(angle);
         for (int inst = 0; inst < ncopy; inst++) {
            double r = r0 + inst * dr;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            Refsys *mref = new Refsys(*dref);
            mref->shift(origin);
            if (volumeS == volume) {
               result->push_back(*mref);
               delete mref;
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,mref,childEl);
               delete mref;
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
                  delete res;
               }
               else {
                  delete res;
                  break;
               }
            }
         }
      }
      else if (elS == "mposX") {
         double x0, dx;
         XString x0S(contEl->getAttribute(X("X0")));
         x0 = atof(S(x0S)) * unit.cm;
         XString dxS(contEl->getAttribute(X("dX")));
         dx = atof(S(dxS)) * unit.cm;
         double y=0, z=0;
         XString yzS(contEl->getAttribute(X("Y_Z")));
         if (yzS.size() > 0) {
            std::stringstream listr(yzS);
            listr >> y >> z;
            y *= unit.cm;
            z *= unit.cm;
         }
         dref->rotate(angle);
         for (int inst = 0; inst < ncopy; inst++) {
            double x = x0 + inst * dx;
            origin[0] = x;
            origin[1] = y + s;
            origin[2] = z;
            Refsys *mref = new Refsys(*dref);
            mref->shift(origin);
            if (volumeS == volume) {
               result->push_back(*mref);
               delete mref;
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,mref,childEl);
               delete mref;
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
                  delete res;
               }
               else {
                  delete res;
                  break;
               }
            }
         }
      }
      else if (elS == "mposY") {
         double y0, dy;
         XString y0S(contEl->getAttribute(X("Y0")));
         y0 = atof(S(y0S)) * unit.cm;
         XString dyS(contEl->getAttribute(X("dY")));
         dy = atof(S(dyS)) * unit.cm;
         double x=0, z=0;
         XString zxS(contEl->getAttribute(X("Z_X")));
         if (zxS.size() > 0) {
            std::stringstream listr(zxS);
            listr >> z >> x;
            x *= unit.cm;
            z *= unit.cm;
         }
         dref->rotate(angle);
         for (int inst = 0; inst < ncopy; inst++) {
            double y = y0 + inst * dy;
            double phi = atan2(y,x);
            origin[0] = x;
            origin[1] = y;
            origin[2] = z;
            Refsys *mref = new Refsys(*dref);
            mref->shift(origin);
            if (volumeS == volume) {
               result->push_back(*mref);
               delete mref;
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,mref,childEl);
               delete mref;
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
                  delete res;
               }
               else {
                  delete res;
                  break;
               }
            }
         }
      }
      else if (elS == "mposZ") {
         double z0, dz;
         XString z0S(contEl->getAttribute(X("Z0")));
         z0 = atof(S(z0S)) * unit.cm;
         XString dzS(contEl->getAttribute(X("dZ")));
         dz = atof(S(dzS)) * unit.cm;
         double x=0, y=0;
         XString xyS(contEl->getAttribute(X("X_Y")));
         XString rphiS(contEl->getAttribute(X("R_Phi")));
         if (xyS.size() > 0) {
            std::stringstream listr(xyS);
            listr >> x >> y;
         }
         else if (rphiS.size() > 0) {
            double r, phi;
            std::stringstream listr(rphiS);
            listr >> r >> phi;
            phi *= unit.rad;
            x = r * cos(phi);
            y = r * sin(phi);
         }
         x *= unit.cm;
         y *= unit.cm;
         dref->rotate(angle);
         double phi = atan2(y,x);
         origin[0] = x - s * sin(phi);
         origin[1] = y + s * cos(phi);
         for (int inst = 0; inst < ncopy; inst++) {
            double z = z0 + inst * dz;
            origin[2] = z;
            Refsys *mref = new Refsys(*dref);
            mref->shift(origin);
            if (volumeS == volume) {
               result->push_back(*mref);
               delete mref;
            }
            else {
               DOMElement* childEl = fGeomDoc->getElementById(X(volumeS));
               std::vector<Refsys> *res = find(volume,mref,childEl);
               delete mref;
               if (res->size() > 0) {
                  std::vector<Refsys>::iterator it = result->end();
                  result->insert(it,res->begin(),res->end());
                  delete res;
               }
               else {
                  delete res;
                  break;
               }
            }
         }
      }
      delete dref;
   }
   return result;
}
