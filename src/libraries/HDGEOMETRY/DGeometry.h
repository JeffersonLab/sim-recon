// $Id$
//
//    File: DGeometry.h
// Created: Thu Apr  3 08:43:06 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#ifndef _DGeometry_
#define _DGeometry_

#include <JANA/jerror.h>
#include <JANA/JGeometry.h>
using namespace jana;

#include <DANA/DApplication.h>
#include "FDC/DFDCGeometry.h"
#include "FDC/DFDCWire.h"

#include <particleType.h>
#include <DVector3.h>
#include "DMaterial.h"
#include "DMaterialMap.h"
using namespace jana;

class DApplication;
class DMagneticFieldMap;
class DLorentzDeflections;


class DGeometry{
	public:
		DGeometry(JGeometry *jgeom, DApplication *dapp, unsigned int runnumber);
		virtual ~DGeometry();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DGeometry";}

		JGeometry* GetJGeometry(void) {return jgeom;}
		DMagneticFieldMap* GetBfield(void) const;
		DLorentzDeflections *GetLorentzDeflections(void);

		// These methods just map to the same methods in JGeometry. Ideally, we'd
		// base DGeometry on JGeometry and so we'd get these automatically.
		// However, that would require a more complicated generator mechanism
		// where the geometry objects are made outside of JANA.
		bool Get(string xpath, string &sval) const {return jgeom->Get(xpath, sval);}
		bool Get(string xpath, map<string, string> &svals) const {return jgeom->Get(xpath, svals);}
		template<class T> bool Get(string xpath, T &val) const {return jgeom->Get(xpath, val);}
		template<class T> bool Get(string xpath, vector<T> &vals, string delimiter=" ") const {return jgeom->Get(xpath, vals, delimiter);}
		template<class T> bool Get(string xpath, map<string,T> &vals) const {return jgeom->Get(xpath, vals);}

		// The GNU 3.2.3 compiler has a problem resolving the ambiguity between
		// Get(string, T&val) and Get(string, vector<T> &vals, string) above.
		// This does not seem to be a problem with the 4.0 compiler. To get
		// around this, some non-templated versions are provided (eeech!).
		bool Get(string xpath, vector<double> &vals, string delimiter=" ") const {return jgeom->Get(xpath, vals, delimiter);}
		bool Get(string xpath, vector<int> &vals, string delimiter=" ") const {return jgeom->Get(xpath, vals, delimiter);}
		bool Get(string xpath, vector<float> &vals, string delimiter=" ") const {return jgeom->Get(xpath, vals, delimiter);}

		bool GetMultiple(string xpath,vector<vector<double> >&vals,
				 string delimiter=" ") const
		{return jgeom->GetMultiple(xpath,vals,delimiter);}
	

		

		typedef pair<string, map<string,string> > node_t;
		typedef vector<node_t> xpathparsed_t;
		
		void FindNodes(string xpath, vector<xpathparsed_t> &matched_xpaths) const;

		// Methods for accessing material map tables obtained from calibDB
		jerror_t FindMat(DVector3 &pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const;
		jerror_t FindMat(DVector3 &pos, double &density, double &A, double &Z, double &RadLen) const;
		jerror_t FindMatALT1(DVector3 &pos, DVector3 &mom, double &KrhoZ_overA, 
				     double &rhoZ_overA,double &LnI,
				     double &X0, double *s_to_boundary=NULL) const;
		jerror_t FindMatKalman(DVector3 &pos, DVector3 &mom,double &Z,
				       double &KrhoZ_overA,
				       double &rhoZ_overA,double &LnI,
				       double &chi2c_factor,double &chi2a_factor,
				       double &chi2a_factor2,
				       unsigned int &last_index,
				       double *s_to_boundary=NULL) const;
		jerror_t FindMatKalman(DVector3 &pos, double &Z,
				      double &KrhoZ_overA,
				       double &rhoZ_overA,double &LnI,
				       double &chi2c_factor,double &chi2a_factor,
				       double &chi2a_factor2,
				       unsigned int &last_index) const;

		const DMaterialMap::MaterialNode* FindMatNode(DVector3 &pos) const;
		const DMaterialMap* FindDMaterialMap(DVector3 &pos) const;
		
		// Convenience methods
		const DMaterial* GetDMaterial(string name);

		bool GetFDCWires(vector<vector<DFDCWire *> >&fdcwires) const;
		bool GetFDCZ(vector<double> &z_wires) const; ///< z-locations for each of the FDC wire planes in cm
		bool GetFDCStereo(vector<double> &stereo_angles) const; ///< stereo angles of each of the FDC wire layers
		bool GetFDCRmin(vector<double> &rmin_packages) const; ///< beam hole size for each FDC package in cm
		bool GetFDCRmax(double &rmax_active_fdc) const; ///< outer radius of FDC active area in cm
		
		bool GetCDCOption(string &cdc_option) const; ///< get the centralDC_option-X string
		bool GetCDCCenterZ(double &cdc_center_z) const; ///< z-location of center of CDC wires in cm
		bool GetCDCAxialLength(double &cdc_axial_length) const; ///< length of CDC axial wires in cm
		bool GetCDCStereo(vector<double> &cdc_stereo) const; ///< stereo angle for each CDC layer in degrees
		bool GetCDCRmid(vector<double> &cdc_rmid) const; ///< Distance of the center of CDC wire from beamline for each layer in cm
		bool GetCDCNwires(vector<int> &cdc_nwires) const; ///< Number of wires for each CDC layer
		bool GetCDCEndplate(double &z,double &dz,double &rmin,double &rmax) const; 

		bool GetBCALRmin(double &bcal_rmin) const; ///< minimum distance of BCAL module from beam line
		bool GetBCALNmodules(unsigned int &bcal_nmodules) const; ///< Number of BCAL modules
		bool GetBCALCenterZ(double &bcal_center_z) const; ///< z-location of center of BCAL module in cm
		bool GetBCALLength(double &bcal_length) const; ///< length of BCAL module in cm
		bool GetBCALDepth(double &bcal_depth) const; ///< depth (or height) of BCAL module in cm
		
		bool GetFCALZ(double &z_fcal) const; ///< z-location of front face of FCAL in cm
		bool GetTOFZ(vector<double> &z_tof) const; ///< z-location of front face of each of TOF in cm
		bool GetTargetZ(double &z_target) const; ///< z-location og center of target
		bool GetTargetLength(double &target_length) const; ///< z-location of center of target
		
	protected:
		DGeometry(){}
		void GetMaterials(void);
		bool GetCompositeMaterial(const string &name, double &density, double &radlen);
	
	private:
		JGeometry *jgeom;
		DApplication *dapp;
		DMagneticFieldMap *bfield;
		vector<DMaterial*> materials;			/// Older implementation to keep track of material specs without ranges
		vector<DMaterialMap*> materialmaps;	/// Material maps generated automatically(indirectly) from XML with ranges and specs
};

#endif // _DGeometry_

