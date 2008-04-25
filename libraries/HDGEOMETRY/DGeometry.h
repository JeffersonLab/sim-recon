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

class DApplication;
class DMagneticFieldMap;
class DLorentzDeflections;

class DGeometry{
	public:
		DGeometry(JGeometry *jgeom, DApplication *dapp);
		virtual ~DGeometry();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DGeometry";}
		
		JGeometry* GetJGeometry(void){return jgeom;}
		DMagneticFieldMap* GetBfield(void);
		DLorentzDeflections *GetLorentzDeflections(void);

		// These methods just map to the same methods in JGeometry. Ideally, we'd
		// base DGeometry on JGeometry and so we'd get these automatically.
		// However, that would require a more complicated generator mechanism
		// where the geometry objects are made outside of JANA.
		bool Get(string xpath, string &sval){return jgeom->Get(xpath, sval);}
		bool Get(string xpath, map<string, string> &svals){return jgeom->Get(xpath, svals);}
		template<class T> bool Get(string xpath, T &val){return jgeom->Get(xpath, val);}
		template<class T> bool Get(string xpath, vector<T> &vals, string delimiter=" "){return jgeom->Get(xpath, vals, delimiter);}
		template<class T> bool Get(string xpath, map<string,T> &vals){return jgeom->Get(xpath, vals);}

		// The GNU 3.2.3 compiler has a problem resolving the ambiguity between
		// Get(string, T&val) and Get(string, vector<T> &vals, string) above.
		// This does not seem to be a problem with the 4.0 compiler. To get
		// around this, some non-templated versions are provided (eeech!).
		bool Get(string xpath, vector<double> &vals, string delimiter=" "){return jgeom->Get(xpath, vals, delimiter);}
		bool Get(string xpath, vector<int> &vals, string delimiter=" "){return jgeom->Get(xpath, vals, delimiter);}
		bool Get(string xpath, vector<float> &vals, string delimiter=" "){return jgeom->Get(xpath, vals, delimiter);}
		
		// Convenience methods
		bool GetFDCZ(vector<double> &z_wires); ///< z-locations for each of the FDC wire planes in cm
		bool GetFDCStereo(vector<double> &stereo_angles); ///< stereo angles of each of the FDC wire layers
		bool GetFDCRmin(vector<double> &rmin_packages); ///< beam hole size for each FDC package in cm
		bool GetFDCRmax(double &rmax_active_fdc); ///< outer radius of FDC active area in cm
		
		bool GetCDCOption(string &cdc_option); ///< get the centralDC_option-X string
		bool GetCDCCenterZ(double &cdc_center_z); ///< z-location of center of CDC wires in cm
		bool GetCDCAxialLength(double &cdc_axial_length); ///< length of CDC axial wires in cm
		bool GetCDCStereo(vector<double> &cdc_stereo); ///< stereo angle for each CDC layer in degrees
		bool GetCDCRmid(vector<double> &cdc_rmid); ///< Distance of the center of CDC wire from beamline for each layer in cm
		bool GetCDCNwires(vector<int> &cdc_nwires); ///< Number of wires for each CDC layer
		
		bool GetBCALRmin(double &bcal_rmin); ///< minimum distance of BCAL module from beam line
		bool GetBCALNmodules(unsigned int &bcal_nmodules); ///< Number of BCAL modules
		bool GetBCALCenterZ(double &bcal_center_z); ///< z-location of center of BCAL module in cm
		bool GetBCALLength(double &bcal_length); ///< length of BCAL module in cm
		bool GetBCALDepth(double &bcal_depth); ///< depth (or height) of BCAL module in cm
		
		bool GetFCALZ(double &z_fcal); ///< z-location of front face of FCAL in cm
		bool GetTOFZ(vector<double> &z_tof); ///< z-location of front face of each of TOF in cm
		bool GetTargetZ(double &z_target); ///< z-location og center of target
		bool GetTargetLength(double &target_length); ///< z-location of center of target
		
	protected:
		DGeometry(){}
	
	private:
		JGeometry *jgeom;
		DApplication *dapp;
};

#endif // _DGeometry_

