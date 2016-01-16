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
  double rhoZ_overA;			// density*Z/A
  double rhoZ_overA_logI;	// density*Z/A * log(mean excitation energy)
};

class DRootGeom{

 public:

  DRootGeom(JApplication *japp, unsigned int run_number=1);
  virtual ~DRootGeom();
  
  virtual const char* className(void){return static_className();}
  static const char* static_className(void){return "DRootGeom";}  

	int ReadMap(string namepath, int32_t runnumber);
	void InitTable(void);

  TGeoNode* GetCurrentNode(){return Current_Node;};
  TGeoVolume* GetCurrentVolume(){return Current_Volume;};
  struct VolMat GetCurrentMat(){return Mat;};
  
  TGeoNode* FindNode(double *x);
  TGeoVolume* FindVolume(double *x);
  struct VolMat FindMat(double *x);

  jerror_t FindMat(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
  jerror_t FindMat(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const;

  jerror_t FindMatLL(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
  jerror_t FindMatLL(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const; 
  jerror_t FindMatLL(DVector3 pos,double &density, double &A, double &Z, double &RadLen,double &LnI) const;

  jerror_t FindMatTable(DVector3 pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
  jerror_t FindMatTable(DVector3 pos,double &density, double &A, double &Z, double &RadLen) const;

  jerror_t FindMat(const char* matname,double &rhoZ_overA,
		   double &rhoZ_overA_logI, double &RadLen) const;

 private:
  
  void InitDRGeom(void);
  
  TGeoManager *DRGeom;
  TGeoNode *Current_Node ;
  TGeoVolume *Current_Volume ;
  TGeoMaterial *Current_Material;
  JCalibration *jcalib;
  Int_t Mat_Index ;
  struct VolMat Mat; // material property : A, Z, Density, RadLen   
	pthread_mutex_t mutex;
	pthread_mutexattr_t mutex_attr;
	
	bool table_initialized;
	VolMat **MatTable; // MatTable[R][Z];
	VolMat *buff;
	int Nr, Nz;		// Number of points in R and Z
	double dr, dz; // Distance between points in R and Z
	double r0, z0;	// Location of first point in R and Z

};

#endif // _DRootGeom_
