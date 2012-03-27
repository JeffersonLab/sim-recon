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
				double chi2c_factor;  // prefactor for the characteristic angle for multiple scattering:  chi2c^2=(chi2c_factor)*X/p^/beta^2
				double chi2a_factor; // prefactor for the screening angle 
				double chi2a_corr;
		};
		
		inline const MaterialNode* FindNode(DVector3 &pos) const;
		
		jerror_t FindMat(DVector3 &pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
		jerror_t FindMatALT1(DVector3 &pos,double &KrhoZ_overA,
				 double &rhoZ_overA, double &logI, 
				 double &RadLen) const;
		jerror_t FindMat(DVector3 &pos, double &density, double &A, double &Z, double &RadLen) const;
		jerror_t FindMatKalman(DVector3 &pos,double &Z,
				       double &K_rho_Z_over_A,
				       double &rho_Z_over_A,double &LogI,double &chi2c_factor,
				       double &chi2a_factor,double &chi2a_corr) const;
		bool IsInMap(const DVector3 &pos) const;
		double EstimatedDistanceToBoundary(const DVector3 &pos, const DVector3 &mom);
		double EstimatedDistanceToBoundarySearch(double r, double z, double p_hatR, double p_hatZ, double &s_to_boundary);
		double DistanceToBox(double &x, double &y, double &xdir, double &ydir, double xmin, double xmax, double ymin, double ymax);

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

		void FindBoundaries(void);

		string namepath;
		vector<vector<MaterialNode> > nodes; // nodes[ir][iz]
		int Nr, Nz;		// Number of nodes in R and Z
		double dr, dz; // Distance between nodes in R and Z
		double r0, z0;	// Location of first nodes in R and Z
		
		double rmin, rmax; // Range limits in R of this map
		double zmin, zmax; // Range limits in Z of this map
		
		int MAX_BOUNDARY_SEARCH_STEPS;
		bool ENABLE_BOUNDARY_CHECK;
		
		bool irregular_density_profile; // Set to true if significant density changes do not follow r or z lines
		vector<double> r_boundaries;
		vector<double> z_boundaries;

		JCalibration *jcalib;
};

//-----------------
// FindNode
//-----------------
inline const DMaterialMap::MaterialNode* DMaterialMap::FindNode(DVector3 &pos) const
{
	// For now, this just finds the bin in the material map the given position is in
	// (i.e. no interpolation )
	double pos_x = pos.X();
	double pos_y = pos.Y();
	double r = sqrt(pos_x*pos_x + pos_y*pos_y);
	double z = pos.Z();
	int ir = (int)floor((r-rmin)/dr);
	int iz = (int)floor((z-zmin)/dz);
	if(ir<0 || ir>=Nr || iz<0 || iz>=Nz)return NULL;
	
	return &nodes[ir][iz];
}


#endif // _DMaterialMap_
