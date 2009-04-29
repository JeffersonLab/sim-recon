// $Id$
//
//    File: DRootGeom.h
// Created: Fri Feb 13 08:43:39 EST 2009
// Creator: zihlmann
//

#include "DRootGeom.h"
#include "hddsroot.h"

using namespace std;


//---------------------------------
// DRootGeom    (Constructor)
//---------------------------------
DRootGeom::DRootGeom()
{

  DRGeom = hddsroot();
 
}

//---------------------------------
// DRootGeom    (Destructor)
//---------------------------------
DRootGeom::~DRootGeom()
{
  delete DRGeom;
}

TGeoNode* DRootGeom::FindNode(double *x)
{

  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  
  return cnode;  

}

TGeoVolume* DRootGeom::FindVolume(double *x)
{

  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  TGeoVolume *cvol = cnode->GetVolume();

  return cvol;  

}

jerror_t DRootGeom::FindMat(DVector3 pos,double &density, double &A, double &Z,
			double &RadLen) const{
  density=RadLen=A=Z=0.;

  TGeoNode *cnode = DRGeom->FindNode(pos.X(),pos.Y(),pos.Z());
  if (cnode==NULL){
    _DBG_<<"Missing cnode at position (" <<pos.X()<<","<<pos.Y()<<","
	 <<pos.Z()<<")"<<endl;
    return RESOURCE_UNAVAILABLE;
  }
  TGeoVolume *cvol = cnode->GetVolume();
  if (cvol==NULL){
    _DBG_<<"Missing cvol" <<endl;
    return RESOURCE_UNAVAILABLE;
  }
  TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

  density=cmat->GetDensity();
  RadLen=cmat->GetRadLen();
  A=cmat->GetA();
  Z=cmat->GetZ();
  if (A<1.){ // try to prevent division by zero problems
    A=1.;
  } 
  return NOERROR;
}


struct VolMat DRootGeom::FindMat(double *x)
{

  TGeoNode *cnode = DRGeom->FindNode(x[0],x[1],x[2]);
  TGeoVolume *cvol = cnode->GetVolume();
  TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

  Current_Node = cnode;
  Current_Volume = cvol;
  Current_Material = cmat;
  Mat_Index = cmat->GetIndex();

  Mat.A = cmat->GetA();
  Mat.Z = cmat->GetZ();
  Mat.Density = cmat->GetDensity();
  Mat.RadLen = cmat->GetRadLen();

  //cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
  //cout<<DRGeom->FindNode(x[0],x[1],x[2])->GetVolume()->GetName()<<endl; 

  return Mat;  

}


