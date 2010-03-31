// Class for dealing with a map of the material in the detector

#ifndef _DMaterialMap_
#define _DMaterialMap_

#include <JANA/jerror.h>

#include <DVector3.h>
#include <DVector2.h>

class DMagneticFieldMap;

class DMaterialMap{
	public:
		DMaterialMap(string namepath, JCalibration *jcalib);
		virtual ~DMaterialMap(){};

		class MaterialNode
		{
			public:
				double A ;
				double Z ;
				double Density ;
				double RadLen ;
				double LogI; // log(mean excitation energy
				double rhoZ_overA; // density*Z/A
				double KrhoZ_overA; // 0.1535e-3*density*Z/A
				double rhoZ_overA_logI;	// density*Z/A * log(mean excitation energy)
		};
		
		const MaterialNode* FindNode(DVector3 &pos) const;
		jerror_t FindMat(DVector3 &pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
		jerror_t FindMatALT1(DVector3 &pos,double &KrhoZ_overA,
				 double &rhoZ_overA, double &logI, 
				 double &RadLen) const;
		jerror_t FindMat(DVector3 &pos, double &density, double &A, double &Z, double &RadLen) const;
		jerror_t FindMatKalman(DVector3 &pos,double &Z,
				       double &K_rho_Z_over_A,
				       double &rho_Z_over_A,double &LogI) const;
		bool IsInMap(const DVector3 &pos) const;
		double EstimatedDistanceToBoundary(const DVector3 &pos, const DVector3 &mom, const DMagneticFieldMap *bfield);
		double DistanceToBox(DVector2 pos, DVector2 dir, double xmin, double xmax, double ymin, double ymax);

		string GetNamepath(void) const {return namepath;}
		double GetRmin(void) const {return rmin;}
		double GetRmax(void) const {return rmax;}
		double GetZmin(void) const {return zmin;}
		double GetZmax(void) const {return zmax;}
		double GetNr(void) const {return Nr;}
		double GetNz(void) const {return Nz;}
		double Getdr(void) const {return dr;}
		double Getdz(void) const {return dz;}
		
	private:
		DMaterialMap(); // Forbid default constructor

		string namepath;
		vector<vector<MaterialNode> > nodes; // nodes[ir][iz]
		int Nr, Nz;		// Number of nodes in R and Z
		double dr, dz; // Distance between nodes in R and Z
		double r0, z0;	// Location of first nodes in R and Z
		
		double rmin, rmax; // Range limits in R of this map
		double zmin, zmax; // Range limits in Z of this map
		
		int MAX_BOUNDARY_SEARCH_STEPS;

		JCalibration *jcalib;
};

#endif // _DMaterialMap_
