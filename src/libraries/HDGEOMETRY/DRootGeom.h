// $Id$
//
//    File: DRootGeom.h
// Created: Fri Feb 13 08:43:39 EST 2009
// Creator: zihlmann
//

#ifndef _DRootGeom_
#define _DRootGeom_

#include <JANA/jerror.h>
#include <DANA/DApplication.h>

#include <DVector3.h>

// the root part
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoPcon.h>
#include <TGeoPgon.h>
#include <TGeoMatrix.h>


struct VolMat
{
  double A ;
  double Z ;
  double Density ;
  double RadLen ;
};

class DRootGeom{

 public:

  DRootGeom();
  virtual ~DRootGeom();
  
  virtual const char* className(void){return static_className();}
  static const char* static_className(void){return "DRootGeom";}  

  TGeoNode* GetCurrentNode(){return Current_Node;};
  TGeoVolume* GetCurrentVolume(){return Current_Volume;};
  struct VolMat GetCurrentMat(){return Mat;};
  
  TGeoNode* FindNode(double *x);
  TGeoVolume* FindVolume(double *x);
  struct VolMat FindMat(double *x);

  void FindMat(DVector3 pos,double &density, double &A, double &Z,
	       double &RadLen) const;

 private:
  
  TGeoManager *DRGeom;
  TGeoNode *Current_Node ;
  TGeoVolume *Current_Volume ;
  TGeoMaterial *Current_Material ;
  Int_t Mat_Index ;
  struct VolMat Mat; // material property : A, Z, Density, RadLen   

};

#endif // _DRootGeom_
