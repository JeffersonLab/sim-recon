// $Id$
//
//    File: DGeometry.cc
// Created: Thu Apr  3 08:43:06 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#include "DGeometry.h"
using namespace std;

//---------------------------------
// DGeometry    (Constructor)
//---------------------------------
DGeometry::DGeometry(JGeometry *jgeom, DApplication *dapp)
{
	this->jgeom = jgeom;
	this->dapp = dapp;
}

//---------------------------------
// ~DGeometry    (Destructor)
//---------------------------------
DGeometry::~DGeometry()
{

}

//---------------------------------
// GetBfield
//---------------------------------
DMagneticFieldMap* DGeometry::GetBfield(void)
{
	return dapp->GetBfield();
}

//---------------------------------
// GetLorentzDeflections
//---------------------------------
DLorentzDeflections* DGeometry::GetLorentzDeflections(void)
{
	return dapp->GetLorentzDeflections();
}

//====================================================================
// Convenience Methods
//
// Below are defined some methods to make it easy to extract certain
// key values about the GlueX detector geometry from the XML source.
// Note that one can still use the generic Get(string xpath, ...) 
// methods. This just packages some of them up for convenience.
//
// The one real gotcha here is that these methods must be kept
// in sync with the XML structure by hand. If volumes are renamed
// or their location within the hierachy is modified, then these
// routines will need to be modified as well. That, or course, is
// also true if you are using the generic Get(...) methods.
//
// What these methods are useful for is when minor changes are made
// to the XML (such as the locations of the FDC packages) they
// are automatically reflected here.
//====================================================================


//---------------------------------
// GetFDCZ
//---------------------------------
bool DGeometry::GetFDCZ(vector<double> &z_wires)
{
	// The FDC geometry is defined as 4 packages, each containing 2
	// "module"s and each of those containing 3 "chambers". The modules
	// are placed as multiple copies in Z using mposZ, but none of the
	// others are (???).
	//
	// This method is currently hardwired to assume 4 packages and
	// 3 chambers. (The number of modules is discovered via the
	// "ncopy" attribute of mposZ.)

	vector<double> ForwardDC;
	vector<double> forwardDC;
	vector<double> forwardDC_package[4];
	map<string,double> forwardDC_module[4];
	vector<double> forwardDC_chamber[4][3];

	bool good = true;
	
	good |= Get("//section/composition/posXYZ[@volume='ForwardDC']/@X_Y_Z", ForwardDC);
	good |= Get("//composition[@name='ForwardDC']/posXYZ[@volume='forwardDC']/@X_Y_Z", forwardDC);
	good |= Get("//posXYZ[@volume='forwardDC_package_1']/@X_Y_Z", forwardDC_package[0]);
	good |= Get("//posXYZ[@volume='forwardDC_package_2']/@X_Y_Z", forwardDC_package[1]);
	good |= Get("//posXYZ[@volume='forwardDC_package_3']/@X_Y_Z", forwardDC_package[2]);
	good |= Get("//posXYZ[@volume='forwardDC_package_4']/@X_Y_Z", forwardDC_package[3]);
	good |= Get("//mposZ[@volume='forwardDC_module_1']", forwardDC_module[0]);
	good |= Get("//mposZ[@volume='forwardDC_module_2']", forwardDC_module[1]);
	good |= Get("//mposZ[@volume='forwardDC_module_3']", forwardDC_module[2]);
	good |= Get("//mposZ[@volume='forwardDC_module_4']", forwardDC_module[3]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[0][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[0][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[0][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[1][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[1][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[1][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[2][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[2][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[2][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[3][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[3][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[3][2]);
	
	if(!good){
		_DBG_<<"Unable to retrieve ForwardDC positions."<<endl;
		return good;
	}
	
	// Offset due to global FDC envelopes
	double zfdc = ForwardDC[2] + forwardDC[2];
	
	// Loop over packages
	for(int package=1; package<=4; package++){
		double z_package = forwardDC_package[package-1][2];
		
		// Loop over modules for this package
		double Z0 = forwardDC_module[package-1]["Z0"];
		double dZ = forwardDC_module[package-1]["dZ"];
		for(int module=0; module<forwardDC_module[package-1]["ncopy"]; module++){
			double z_module = Z0 + (double)module*dZ;
			
			// Loop over chambers in
			for(int chamber=0; chamber<3; chamber++){
				double z_chamber = forwardDC_chamber[package-1][chamber][2];
				
				double z = zfdc + z_package + z_module + z_chamber;				
				z_wires.push_back(z);
			}
		}
	}
	
	return good;
}

//---------------------------------
// GetFDCStereo
//---------------------------------
bool DGeometry::GetFDCStereo(vector<double> &stereo_angles)
{
	// The FDC geometry is defined as 4 packages, each containing 2
	// "module"s and each of those containing 3 "chambers". The modules
	// are placed as multiple copies in Z using mposZ, but none of the
	// others are (???).
	//
	// This method is currently hardwired to assume 4 packages and
	// 3 chambers. (The number of modules is discovered via the
	// "ncopy" attribute of mposZ.)
	//
	// Stereo angles are assumed to be rotated purely about the z-axis
	// and the units are not specified, but the XML currently uses degrees.

	map<string,double> forwardDC_module[4];
	vector<double> forwardDC_chamber[4][3];

	bool good = true;
	
	good |= Get("//mposZ[@volume='forwardDC_module_1']", forwardDC_module[0]);
	good |= Get("//mposZ[@volume='forwardDC_module_2']", forwardDC_module[1]);
	good |= Get("//mposZ[@volume='forwardDC_module_3']", forwardDC_module[2]);
	good |= Get("//mposZ[@volume='forwardDC_module_4']", forwardDC_module[3]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='1']", forwardDC_chamber[0][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='2']", forwardDC_chamber[0][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='3']", forwardDC_chamber[0][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='1']", forwardDC_chamber[1][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='2']", forwardDC_chamber[1][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='3']", forwardDC_chamber[1][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='1']", forwardDC_chamber[2][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='2']", forwardDC_chamber[2][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='3']", forwardDC_chamber[2][2]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='1']", forwardDC_chamber[3][0]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='2']", forwardDC_chamber[3][1]);
	good |= Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='3']", forwardDC_chamber[3][2]);
	
	if(!good){
		_DBG_<<"Unable to retrieve ForwardDC positions."<<endl;
		return good;
	}
	
	// Loop over packages
	for(int package=1; package<=4; package++){
		
		// Loop over modules for this package
		for(int module=0; module<forwardDC_module[package-1]["ncopy"]; module++){
			
			// Loop over chambers
			for(int chamber=0; chamber<3; chamber++){
				stereo_angles.push_back(forwardDC_chamber[package-1][chamber][2]);
			}
		}
	}

	return good;
}

//---------------------------------
// GetFDCRmin
//---------------------------------
bool DGeometry::GetFDCRmin(vector<double> &rmin_packages)
{
	vector<double> FDA[4];

	bool good = true;
	
	good |= Get("//section[@name='ForwardDC']/tubs[@name='FDA1']/@Rio_Z", FDA[0]);
	good |= Get("//section[@name='ForwardDC']/tubs[@name='FDA2']/@Rio_Z", FDA[1]);
	good |= Get("//section[@name='ForwardDC']/tubs[@name='FDA3']/@Rio_Z", FDA[2]);
	good |= Get("//section[@name='ForwardDC']/tubs[@name='FDA4']/@Rio_Z", FDA[3]);
	
	if(!good){
		_DBG_<<"Unable to retrieve FDC Rmin values."<<endl;
		return good;
	}

	rmin_packages.push_back(FDA[0][0]);
	rmin_packages.push_back(FDA[1][0]);
	rmin_packages.push_back(FDA[2][0]);
	rmin_packages.push_back(FDA[3][0]);

	return good;
}

//---------------------------------
// GetFDCRmax
//---------------------------------
bool DGeometry::GetFDCRmax(double &rmax_active_fdc)
{
	// We assume that all packages have the same outer radius of the
	// active area.
	vector<double> FDA1;

	bool good = true;
	
	good |= Get("//section[@name='ForwardDC']/tubs[@name='FDA1']/@Rio_Z", FDA1);
	
	if(!good){
		_DBG_<<"Unable to retrieve FDC Rmax values."<<endl;
		return good;
	}

	rmax_active_fdc = FDA1[1];

	return good;
}

//---------------------------------
// GetCDCOption
//---------------------------------
bool DGeometry::GetCDCOption(string &cdc_option)
{
	bool good = Get("//CentralDC_s/section/composition/posXYZ/@volume", cdc_option);
	
	if(!good){
		_DBG_<<"Unable to retrieve CDC option string."<<endl;
	}

	return good;
}

//---------------------------------
// GetCDCCenterZ
//---------------------------------
bool DGeometry::GetCDCCenterZ(double &cdc_center_z)
{

	return false;
}

//---------------------------------
// GetCDCAxialLength
//---------------------------------
bool DGeometry::GetCDCAxialLength(double &cdc_axial_length)
{

	return false;
}

//---------------------------------
// GetCDCStereo
//---------------------------------
bool DGeometry::GetCDCStereo(vector<double> &cdc_stereo)
{

	return false;
}

//---------------------------------
// GetCDCRmid
//---------------------------------
bool DGeometry::GetCDCRmid(vector<double> &cdc_rmid)
{

	return false;
}

//---------------------------------
// GetCDCNwires
//---------------------------------
bool DGeometry::GetCDCNwires(vector<int> &cdc_nwires)
{

	return false;
}

//---------------------------------
// GetBCALRmin
//---------------------------------
bool DGeometry::GetBCALRmin(double &bcal_rmin)
{

	return false;
}

//---------------------------------
// GetBCALNmodules
//---------------------------------
bool DGeometry::GetBCALNmodules(unsigned int &bcal_nmodules)
{

	return false;
}

//---------------------------------
// GetBCALCenterZ
//---------------------------------
bool DGeometry::GetBCALCenterZ(double &bcal_center_z)
{

	return false;
}

//---------------------------------
// GetBCALLength
//---------------------------------
bool DGeometry::GetBCALLength(double &bcal_length)
{

	return false;
}

//---------------------------------
// GetBCALDepth
//---------------------------------
bool DGeometry::GetBCALDepth(double &bcal_depth)
{

	return false;
}

//---------------------------------
// GetFCALZ
//---------------------------------
bool DGeometry::GetFCALZ(double &z_fcal)
{
	vector<double> ForwardEMcalpos;
	bool good = Get("//section/composition/posXYZ[@volume='ForwardEMcal']/@X_Y_Z", ForwardEMcalpos);
	
	if(!good){
		_DBG_<<"Unable to retrieve ForwardEMcal position."<<endl;
		z_fcal=0.0;
		return false;
	}else{
		z_fcal = ForwardEMcalpos[2];
		return true;
	}
}

//---------------------------------
// GetTOFZ
//---------------------------------
bool DGeometry::GetTOFZ(vector<double> &z_tof)
{
	vector<double> ForwardTOF;
	vector<double> forwardTOF[2];
	vector<double> FTOC;
	bool good = true;
	good |= Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z", ForwardTOF);
	good |= Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", forwardTOF[0]);
	good |= Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", forwardTOF[1]);
	good |= Get("//box[@name='FTOC' and sensitive='true']/@X_Y_Z", FTOC);
	
	z_tof.push_back(ForwardTOF[2] + forwardTOF[0][2] - FTOC[2]/2.0);
	z_tof.push_back(ForwardTOF[2] + forwardTOF[1][2] - FTOC[2]/2.0);

	return good;
}

//---------------------------------
// GetTargetZ
//---------------------------------
bool DGeometry::GetTargetZ(double &z_target)
{

	return false;
}

//---------------------------------
// GetTargetLength
//---------------------------------
bool DGeometry::GetTargetLength(double &target_length)
{

	return false;
}

