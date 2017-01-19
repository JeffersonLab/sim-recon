// $Id$
//
//    File: DCDCWire.h
// Created: Wed Oct 18 11:39:35 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DCDCWire_
#define _DCDCWire_

#include <DCoordinateSystem.h>

enum CDCWireD {
   dOriginXddeltaX=0,
   dOriginXddeltaY,
   dOriginXddeltaZ,
   dOriginXddeltaPhiX,
   dOriginXddeltaPhiY,
   dOriginXddeltaPhiZ,
   dOriginXddeltaXu,
   dOriginXddeltaYu,
   dOriginXddeltaXd,
   dOriginXddeltaYd,
   dOriginYddeltaX,
   dOriginYddeltaY,
   dOriginYddeltaZ,
   dOriginYddeltaPhiX,
   dOriginYddeltaPhiY,
   dOriginYddeltaPhiZ,
   dOriginYddeltaXu,
   dOriginYddeltaYu,
   dOriginYddeltaXd,
   dOriginYddeltaYd,
   dOriginZddeltaX,
   dOriginZddeltaY,
   dOriginZddeltaZ,
   dOriginZddeltaPhiX,
   dOriginZddeltaPhiY,
   dOriginZddeltaPhiZ,
   dOriginZddeltaXu,
   dOriginZddeltaYu,
   dOriginZddeltaXd,
   dOriginZddeltaYd,
   dDirXddeltaX,
   dDirXddeltaY,
   dDirXddeltaZ,
   dDirXddeltaPhiX,
   dDirXddeltaPhiY,
   dDirXddeltaPhiZ,
   dDirXddeltaXu,
   dDirXddeltaYu,
   dDirXddeltaXd,
   dDirXddeltaYd,
   dDirYddeltaX,
   dDirYddeltaY,
   dDirYddeltaZ,
   dDirYddeltaPhiX,
   dDirYddeltaPhiY,
   dDirYddeltaPhiZ,
   dDirYddeltaXu,
   dDirYddeltaYu,
   dDirYddeltaXd,
   dDirYddeltaYd,
   dDirZddeltaX,
   dDirZddeltaY,
   dDirZddeltaZ,
   dDirZddeltaPhiX,
   dDirZddeltaPhiY,
   dDirZddeltaPhiZ,
   dDirZddeltaXu,
   dDirZddeltaYu,
   dDirZddeltaXd,
   dDirZddeltaYd
};

class DCDCWire:public DCoordinateSystem{
   public:
      float phi;		// phi angle of wire at midplane in lab coordinates

      int ring;
      int straw;
      float stereo;
      // Added for alignment derivatives
      // These are values that go into the calculation of origin and udir
      double phiX,phiY,phiZ;
      double x0, y0, z0;
      double udir_mag;
      double stereo_raw;
      double stereo_sign;
      double r0;
      double phiStraw;

      vector<double> derivatives;

      vector<double> GetCDCWireDerivatives(){
         // Only need to calculate this once
         if (derivatives.size() != 0) return derivatives;
         derivatives.resize(60);
         
         double sinPhiX=sin(phiX); double cosPhiX=cos(phiX);
         double sinPhiY=sin(phiY); double cosPhiY=cos(phiY);
         double sinPhiZ=sin(phiZ); double cosPhiZ=cos(phiZ);
         double sinPhiStraw=sin(phiStraw); double cosPhiStraw=cos(phiStraw);
         double tanStereo=tan(stereo_raw);

         double x0Shift = x0+r0*cosPhiStraw;
         double y0Shift = y0+r0*sinPhiStraw;

         derivatives[CDCWireD::dOriginXddeltaX]=cosPhiY*cosPhiZ;
         derivatives[CDCWireD::dOriginXddeltaY]=cosPhiZ*sinPhiX*sinPhiY-cosPhiX*sinPhiZ;
         derivatives[CDCWireD::dOriginXddeltaZ]=cosPhiX*cosPhiZ*sinPhiY+sinPhiX*sinPhiZ;
         derivatives[CDCWireD::dOriginXddeltaPhiX]=cosPhiX*(cosPhiZ*sinPhiY*y0Shift+z0*sinPhiZ)
            +sinPhiX*(-z0*cosPhiZ*sinPhiY+y0Shift*sinPhiZ);
         derivatives[CDCWireD::dOriginXddeltaPhiY]=cosPhiZ*(z0*cosPhiX*cosPhiY
               +cosPhiY*y0Shift*sinPhiX-x0Shift*sinPhiY);
         derivatives[CDCWireD::dOriginXddeltaPhiZ]=z0*cosPhiZ*sinPhiX-(x0Shift*cosPhiY+y0Shift*sinPhiX*sinPhiY)*sinPhiZ
            -cosPhiX*(cosPhiZ*y0Shift+z0*sinPhiY*sinPhiZ);
         derivatives[CDCWireD::dOriginXddeltaXu]=0.5*derivatives[CDCWireD::dOriginXddeltaX];
         derivatives[CDCWireD::dOriginXddeltaYu]=0.5*derivatives[CDCWireD::dOriginXddeltaY];
         derivatives[CDCWireD::dOriginXddeltaXd]=0.5*derivatives[CDCWireD::dOriginXddeltaX];
         derivatives[CDCWireD::dOriginXddeltaYd]=0.5*derivatives[CDCWireD::dOriginXddeltaY];
         derivatives[CDCWireD::dOriginYddeltaX]=cosPhiY*sinPhiZ;
         derivatives[CDCWireD::dOriginYddeltaY]=cosPhiX*cosPhiZ+sinPhiX*sinPhiY*sinPhiZ;
         derivatives[CDCWireD::dOriginYddeltaZ]=-cosPhiZ*sinPhiX+cosPhiX*sinPhiY*sinPhiZ;
         derivatives[CDCWireD::dOriginYddeltaPhiX]=-cosPhiZ*(z0*cosPhiX+y0Shift*sinPhiX)
            +sinPhiY*sinPhiZ*(y0Shift*cosPhiX-z0*sinPhiX);
         derivatives[CDCWireD::dOriginYddeltaPhiY]=(z0*cosPhiX*cosPhiY+y0Shift*cosPhiY*sinPhiX-x0Shift*sinPhiY)*sinPhiZ;
         derivatives[CDCWireD::dOriginYddeltaPhiZ]=x0Shift*cosPhiY*cosPhiZ
            +sinPhiX*(y0Shift*cosPhiZ*sinPhiY+z0*sinPhiZ)
            +cosPhiX*(z0*cosPhiZ*sinPhiY-y0Shift*sinPhiZ);
         derivatives[CDCWireD::dOriginYddeltaXu]=0.5*derivatives[CDCWireD::dOriginYddeltaX];
         derivatives[CDCWireD::dOriginYddeltaYu]=0.5*derivatives[CDCWireD::dOriginYddeltaY];
         derivatives[CDCWireD::dOriginYddeltaXd]=0.5*derivatives[CDCWireD::dOriginYddeltaX];
         derivatives[CDCWireD::dOriginYddeltaYd]=0.5*derivatives[CDCWireD::dOriginYddeltaY];
         derivatives[CDCWireD::dOriginZddeltaX]=-sinPhiY;
         derivatives[CDCWireD::dOriginZddeltaY]=cosPhiY*sinPhiX;
         derivatives[CDCWireD::dOriginZddeltaZ]=cosPhiX*cosPhiY;
         derivatives[CDCWireD::dOriginZddeltaPhiX]=cosPhiY*(y0Shift*cosPhiX-z0*sinPhiX);
         derivatives[CDCWireD::dOriginZddeltaPhiY]=-x0Shift*cosPhiY-sinPhiY*(z0*cosPhiX+y0Shift*sinPhiX);
         derivatives[CDCWireD::dOriginZddeltaPhiZ]=0.0;
         derivatives[CDCWireD::dOriginZddeltaXu]=0.5*derivatives[CDCWireD::dOriginZddeltaX];
         derivatives[CDCWireD::dOriginZddeltaYu]=0.5*derivatives[CDCWireD::dOriginZddeltaY];
         derivatives[CDCWireD::dOriginZddeltaXd]=0.5*derivatives[CDCWireD::dOriginZddeltaX];
         derivatives[CDCWireD::dOriginZddeltaYd]=0.5*derivatives[CDCWireD::dOriginZddeltaY];
         derivatives[CDCWireD::dDirXddeltaX]=0.0;
         derivatives[CDCWireD::dDirXddeltaY]=0.0;
         derivatives[CDCWireD::dDirXddeltaZ]=0.0;
         derivatives[CDCWireD::dDirXddeltaPhiX]=-L*(cosPhiZ*sinPhiY*(sinPhiX+cosPhiStraw*cosPhiX*tanStereo)
               +sinPhiZ*(cosPhiStraw*sinPhiX*tanStereo-cosPhiX));
         derivatives[CDCWireD::dDirXddeltaPhiY]=L*cosPhiZ*(cosPhiX*cosPhiY-(cosPhiStraw*cosPhiY*sinPhiX+sinPhiStraw*sinPhiY)*tanStereo);
         derivatives[CDCWireD::dDirXddeltaPhiZ]=L*(cosPhiZ*(sinPhiX+cosPhiStraw*cosPhiX*tanStereo)
               -sinPhiZ*(cosPhiX*sinPhiY
                  +(cosPhiY*sinPhiStraw-cosPhiStraw*sinPhiX*sinPhiY)*tanStereo));
         derivatives[CDCWireD::dDirXddeltaXu]=-cosPhiY*cosPhiZ;
         derivatives[CDCWireD::dDirXddeltaYu]=-cosPhiZ*sinPhiX*sinPhiY+cosPhiX*sinPhiZ;
         derivatives[CDCWireD::dDirXddeltaXd]=-derivatives[CDCWireD::dDirXddeltaXu];
         derivatives[CDCWireD::dDirXddeltaYd]=-derivatives[CDCWireD::dDirXddeltaYu];
         derivatives[CDCWireD::dDirYddeltaX]=0.0;
         derivatives[CDCWireD::dDirYddeltaY]=0.0;
         derivatives[CDCWireD::dDirYddeltaZ]=0.0;
         derivatives[CDCWireD::dDirYddeltaPhiX]=-L*(sinPhiX*(sinPhiY*sinPhiZ-cosPhiStraw*cosPhiZ*tanStereo)
               +cosPhiX*(cosPhiZ+cosPhiStraw*sinPhiY*sinPhiZ*tanStereo));
         derivatives[CDCWireD::dDirYddeltaPhiY]=L*sinPhiZ*(cosPhiX*cosPhiY
               -(cosPhiStraw*cosPhiY*sinPhiX+sinPhiStraw*sinPhiY)*tanStereo);
         derivatives[CDCWireD::dDirYddeltaPhiZ]=L*(cosPhiY*cosPhiZ*sinPhiStraw*tanStereo
               +sinPhiX*(sinPhiZ-cosPhiStraw*cosPhiZ*sinPhiY*tanStereo)
               +cosPhiX*(cosPhiZ*sinPhiY+cosPhiStraw*sinPhiZ*tanStereo));
         derivatives[CDCWireD::dDirYddeltaXu]=-cosPhiY*sinPhiZ;
         derivatives[CDCWireD::dDirYddeltaYu]=-cosPhiX*cosPhiZ-sinPhiX*sinPhiY*sinPhiZ;
         derivatives[CDCWireD::dDirYddeltaXd]=-derivatives[CDCWireD::dDirYddeltaXu];
         derivatives[CDCWireD::dDirYddeltaYd]=-derivatives[CDCWireD::dDirYddeltaYu];
         derivatives[CDCWireD::dDirZddeltaX]=0.0;
         derivatives[CDCWireD::dDirZddeltaY]=0.0;
         derivatives[CDCWireD::dDirZddeltaZ]=0.0;
         derivatives[CDCWireD::dDirZddeltaPhiX]=-L*cosPhiY*(sinPhiX+cosPhiStraw*cosPhiX*tanStereo);
         derivatives[CDCWireD::dDirZddeltaPhiY]=-L*(cosPhiX*sinPhiY
               +(cosPhiY*sinPhiStraw-cosPhiStraw*sinPhiX*sinPhiY)*tanStereo);
         derivatives[CDCWireD::dDirZddeltaPhiZ]=0.0;
         derivatives[CDCWireD::dDirZddeltaXu]=sinPhiY;
         derivatives[CDCWireD::dDirZddeltaYu]=-cosPhiY*sinPhiX;
         derivatives[CDCWireD::dDirZddeltaXd]=-derivatives[CDCWireD::dDirZddeltaXu];
         derivatives[CDCWireD::dDirZddeltaYd]=-derivatives[CDCWireD::dDirZddeltaYu];

         return derivatives;
      }
};

#endif // _DCDCWire_

