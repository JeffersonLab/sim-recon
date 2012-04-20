// $Id$
//
//    File: DGeometry.cc
// Created: Thu Apr  3 08:43:06 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#include <algorithm>
using namespace std;

#include "DGeometry.h"
#include <JANA/JGeometryXML.h>
#include "FDC/DFDCWire.h"
#include "FDC/DFDCGeometry.h"
#include <ansi_escape.h>

using namespace std;


//---------------------------------
// DGeometry    (Constructor)
//---------------------------------
DGeometry::DGeometry(JGeometry *jgeom, DApplication *dapp, unsigned int runnumber)
{
	this->jgeom = jgeom;
	this->dapp = dapp;
	this->bfield = dapp->GetBfield();

	JCalibration *jcalib = dapp->GetJCalibration(runnumber);
	if(!jcalib){
		_DBG_<<"ERROR:  Unable to get JCalibration object!"<<endl;
		return;
	}
	
	// Get a list of all namepaths so we can find all of the material maps.
	// We want to read in all material maps with a standard naming scheme
	// so the number and coverage of the piece-wise maps can be changed
	// without requiring recompilation.
	//vector<DMaterialMap*> material_maps;
	vector<string> namepaths;
	jcalib->GetListOfNamepaths(namepaths);
	vector<string> material_namepaths;
	for(unsigned int i=0; i<namepaths.size(); i++){
		if(namepaths[i].find("Material/material_map")==0)material_namepaths.push_back(namepaths[i]);
	}
	
	// Sort alphabetically so user controls search sequence by map name
	sort(material_namepaths.begin(), material_namepaths.end());
	
	// Inform user what's happening
	if(material_namepaths.size()==0){
		jerr<<"No material maps found in calibration DB!!"<<endl;
		return;
	}
	jout<<"Found "<<material_namepaths.size()<<" material maps in calib. DB for run "<<runnumber<<endl;
	
	if(false){ // save this to work off configuration parameter
		jout<<"Will read in the following:"<<endl;
		for(unsigned int i=0; i<material_namepaths.size(); i++){
			jout<<"  "<<material_namepaths[i]<<endl;
		}
	}

	// Actually read in the maps
	unsigned int Npoints_total=0;
	//cout<<endl; // make empty line so material map can overwrite it below
	for(unsigned int i=0; i<material_namepaths.size(); i++){
		// DMaterialMap constructor prints line so we conserve real
		// estate by having each recycle the line
		//cout<<ansi_up(1)<<string(85, ' ')<<"\r";
		DMaterialMap *mat = new DMaterialMap(material_namepaths[i], jcalib);
		materialmaps.push_back(mat);
		Npoints_total += (unsigned int)(mat->GetNr()*mat->GetNz());
	}
	//cout<<ansi_up(1)<<string(85, ' ')<<"\r";
	jout<<"Read in "<<materialmaps.size()<<" material maps containing "<<Npoints_total<<" grid points total"<<endl;
}

//---------------------------------
// ~DGeometry    (Destructor)
//---------------------------------
DGeometry::~DGeometry()
{
	for(unsigned int i=0; i<materials.size(); i++)delete materials[i];
	materials.clear();
}

//---------------------------------
// GetBfield
//---------------------------------
DMagneticFieldMap* DGeometry::GetBfield(void) const
{
	return bfield;
}

//---------------------------------
// GetLorentzDeflections
//---------------------------------
DLorentzDeflections* DGeometry::GetLorentzDeflections(void)
{
	return dapp->GetLorentzDeflections();
}

//---------------------------------
// FindNodes
//---------------------------------
void DGeometry::FindNodes(string xpath, vector<xpathparsed_t> &matched_xpaths) const
{
	/// Find all nodes that match the specified xpath and return them as
	/// fully parsed lists of the nodes and attributes.
	///
	/// The matched_xpaths variable has 4 levels of STL containers nested
	/// together! The node_t data type contains 2 of them and represents
	/// a single tag with the "first" member being the tag name
	/// and the "second" member being a map of all of the attributes
	/// of that tag along with their values.
	///
	/// The xpathparsed_t data type is a STL vector of node_t objects
	/// that comprises a complete xpath. The data type of matched_xpaths
	/// is a STL vector if xpathparsed_t objects and so comprises a
	/// list of complete xpaths.
	
	/// We do this by calling the GetXPaths() method of JGeometry to get
	/// a list of all xpaths. Then we pass all of those in to 
	/// JGeometryXML::ParseXPath() to get a parsed list for each. This
	/// is compared to the parsed values of the xpath passed to us
	/// (also parsed by JGeometryXML::ParseXPath()) to find matches.
	
	// Make sure matched_xpaths is empty
	matched_xpaths.clear();
	
	// Cast JGeometry into a JGeometryXML
	JGeometryXML *jgeomxml = dynamic_cast<JGeometryXML*>(jgeom);
	
	// Parse our target xpath
	xpathparsed_t target;
	string unused_string;
	unsigned int unused_int;
	jgeomxml->ParseXPath(xpath, target, unused_string, unused_int);
	
	// Get all xpaths for current geometry source
	vector<string> allxpaths;
	jgeom->GetXPaths(allxpaths, JGeometry::attr_level_all);
	
	// Loop over xpaths
	for(unsigned int i=0; i<allxpaths.size(); i++);
}

//---------------------------------
// FindMatALT1 - parameters for alt1 fitter.
//---------------------------------
jerror_t DGeometry::FindMatALT1(DVector3 &pos, DVector3 &mom,double &KrhoZ_overA, 
				double &rhoZ_overA,double &LnI,
				double &X0, double *s_to_boundary) const
{
	for(unsigned int i=0; i<materialmaps.size(); i++){
		jerror_t err = materialmaps[i]->FindMatALT1(pos,KrhoZ_overA, rhoZ_overA,LnI,X0);
		if(err==NOERROR){
			// We found the material map containing this point. If a non-NULL 
			// pointer was passed in for s_to_boundary, then search through all
			// material maps above and including this one to find the estimated
			// distance to the nearest boundary the swimmer should step to. Maps
			// after this one would only be considered outside of this one so
			// there is no need to check them.
			if(s_to_boundary==NULL)return NOERROR;	// User doesn't want distance to boundary
			*s_to_boundary = 1.0E6;
			for(unsigned int j=0; j<=i; j++){
				double s = materialmaps[j]->EstimatedDistanceToBoundary(pos, mom);
				if(s<*s_to_boundary)*s_to_boundary = s;
			}
			return NOERROR;
		}
	}
	return RESOURCE_UNAVAILABLE;
}




//---------------------------------
// FindMatKalman - Kalman filter needs slightly different set of parms.
//---------------------------------
jerror_t DGeometry::FindMatKalman(DVector3 &pos, DVector3 &mom,double &Z,
				  double &KrhoZ_overA, 
				  double &rhoZ_overA, 
				  double &LnI,
				  double &chi2c_factor,double &chi2a_factor,double &chi2a_corr,
				  unsigned int &last_index,
				  double *s_to_boundary) const
{
  //last_index=0;
  for(unsigned int i=last_index; i<materialmaps.size(); i++){
    jerror_t err = materialmaps[i]->FindMatKalman(pos,Z,KrhoZ_overA,
						  rhoZ_overA,LnI,chi2c_factor,chi2a_factor,chi2a_corr);
    if(err==NOERROR){
      if(i==materialmaps.size()-1) last_index=0;
      else last_index=i;
      if(s_to_boundary==NULL)return NOERROR;	// User doesn't want distance to boundary

      *s_to_boundary = 1.0E6;
      // If we are in the main mother volume, search through all the maps for
      // the nearest boundary
      if(last_index==0){
	for(unsigned int j=0; j<materialmaps.size();j++){
	  double s = materialmaps[j]->EstimatedDistanceToBoundary(pos, mom);
	  if(s<*s_to_boundary){
	    *s_to_boundary = s;
	  }
	}
      }
      else{
	// otherwise, we found the material map containing this point. 
	double s = materialmaps[last_index]->EstimatedDistanceToBoundary(pos, mom);
	if(s<*s_to_boundary)*s_to_boundary = s;
      }
      return NOERROR;
    }
  }
       
return RESOURCE_UNAVAILABLE;
}

//---------------------------------
jerror_t DGeometry::FindMatKalman(DVector3 &pos,double &Z,
				  double &KrhoZ_overA, 
				  double &rhoZ_overA, 
				  double &LnI, double &chi2c_factor, double &chi2a_factor, double &chi2a_corr,
				  unsigned int &last_index) const
{
  //last_index=0;
  for(unsigned int i=last_index; i<materialmaps.size(); i++){
    jerror_t err = materialmaps[i]->FindMatKalman(pos,Z,KrhoZ_overA,
						  rhoZ_overA,LnI,chi2c_factor,chi2a_factor,chi2a_corr);
    if(err==NOERROR){
      if(i==materialmaps.size()-1) last_index=0;
      else last_index=i;
      return err;
    }
  }
       
  return RESOURCE_UNAVAILABLE;
}

//---------------------------------
// FindMat
//---------------------------------
jerror_t DGeometry::FindMat(DVector3 &pos, double &rhoZ_overA, double &rhoZ_overA_logI, double &RadLen) const
{
	for(unsigned int i=0; i<materialmaps.size(); i++){
		jerror_t err = materialmaps[i]->FindMat(pos, rhoZ_overA, rhoZ_overA_logI, RadLen);
		if(err==NOERROR)return NOERROR;
	}
	return RESOURCE_UNAVAILABLE;
}

//---------------------------------
// FindMat
//---------------------------------
jerror_t DGeometry::FindMat(DVector3 &pos, double &density, double &A, double &Z, double &RadLen) const
{
	for(unsigned int i=0; i<materialmaps.size(); i++){
		jerror_t err = materialmaps[i]->FindMat(pos, density, A, Z, RadLen);
		if(err==NOERROR)return NOERROR;
	}
	return RESOURCE_UNAVAILABLE;
}

//---------------------------------
// FindMatNode
//---------------------------------
//const DMaterialMap::MaterialNode* DGeometry::FindMatNode(DVector3 &pos) const
//{
//
//}

//---------------------------------
// FindDMaterialMap
//---------------------------------
const DMaterialMap* DGeometry::FindDMaterialMap(DVector3 &pos) const
{
	for(unsigned int i=0; i<materialmaps.size(); i++){
		const DMaterialMap* map = materialmaps[i];
		if(map->IsInMap(pos))return map;
	}
	return NULL;
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
// GetDMaterial
//---------------------------------
const DMaterial* DGeometry::GetDMaterial(string name)
{
	/// Get a pointer to the DMaterial object with the specified name.
	// Only fill materials member when one is actually requested
	// and then, only fill it once.
	if(materials.size()==0)GetMaterials();

	for(unsigned int i=0; i<materials.size(); i++){
		if(materials[i]->GetName() == name)return materials[i];
	}
	
	//_DBG_<<"No material \""<<name<<"\" found ("<<materials.size()<<" materials defined)"<<endl;
	
	return NULL;
}

//---------------------------------
// GetMaterials
//---------------------------------
void DGeometry::GetMaterials(void)
{
	/// Read in all of the materials from the geometry source and create
	/// a DMaterial object for each one.
	
	//=========== elements ===========
	string filter = "//materials/element/real[@name=\"radlen\"]";
	
	// Get list of all xpaths
	vector<string> xpaths;
	jgeom->GetXPaths(xpaths, JGeometry::attr_level_all, filter);
	
	// Look for xpaths that have "/materials[" in them
	for(unsigned int i=0; i<xpaths.size(); i++){
		// Get start of "element" section
		string::size_type pos = xpaths[i].find("/element[");
		if(pos == string::npos)continue;
		
		// Get name attribute
		string::size_type start_pos = xpaths[i].find("@name=", pos);
		start_pos = xpaths[i].find("'", start_pos);
		string::size_type end_pos = xpaths[i].find("'", start_pos+1);
		if(end_pos==string::npos)continue;
		string name = xpaths[i].substr(start_pos+1, end_pos-(start_pos+1));

		// Get values
		char xpath[256];

		double A;
		sprintf(xpath,"//materials/element[@name='%s']/[@a]", name.c_str());
		if(!Get(xpath, A))continue;

		double Z;
		sprintf(xpath,"//materials/element[@name='%s']/[@z]", name.c_str());
		if(!Get(xpath, Z))continue;

		double density;
		sprintf(xpath,"//materials/element[@name='%s']/real[@name='density']/[@value]", name.c_str());
		if(!Get(xpath, density))continue;

		double radlen;
		sprintf(xpath,"//materials/element[@name='%s']/real[@name='radlen']/[@value]", name.c_str());
		if(!Get(xpath, radlen))continue;

		DMaterial *mat = new DMaterial(name, A, Z, density, radlen);
		materials.push_back(mat);
		
		cout<<mat->toString();
	}

	//=========== composites ===========
	filter = "//materials/composite[@name]";
	
	// Get list of all xpaths
	jgeom->GetXPaths(xpaths, JGeometry::attr_level_all, filter);

	// Look for xpaths that have "/materials[" in them
	for(unsigned int i=0; i<xpaths.size(); i++){
		// Get start of "composite" section
		string::size_type pos = xpaths[i].find("/composite[");
		if(pos == string::npos)continue;
		
		// Get name attribute
		string::size_type start_pos = xpaths[i].find("@name=", pos);
		start_pos = xpaths[i].find("'", start_pos);
		string::size_type end_pos = xpaths[i].find("'", start_pos+1);
		if(end_pos==string::npos)continue;
		string name = xpaths[i].substr(start_pos+1, end_pos-(start_pos+1));

		if(GetDMaterial(name))continue; // skip duplicates

		// Get values
		char xpath[256];

		// We should calculate an effective A and Z .... but we don't
		bool found_all=true;
		double A=0;
		double Z=0;

		double density;
		sprintf(xpath,"//materials/composite[@name='%s']/real[@name='density']/[@value]", name.c_str());
		found_all &= Get(xpath, density);

		double radlen;
		sprintf(xpath,"//materials/composite[@name='%s']/real[@name='radlen']/[@value]", name.c_str());
		found_all &= Get(xpath, radlen);
		
		// If we didn't find the info we need (radlen and density) in the 
		// composite tag itself. Try calculating them from the components
		if(!found_all)found_all = GetCompositeMaterial(name, density, radlen);
		
		// If we weren't able to get the values needed to make the DMaterial object
		// then skip this one.
		if(!found_all)continue;

		DMaterial *mat = new DMaterial(name, A, Z, density, radlen);
		materials.push_back(mat);
		
		cout<<mat->toString();
	}
}

//---------------------------------
// GetCompositeMaterial
//---------------------------------
bool DGeometry::GetCompositeMaterial(const string &name, double &density, double &radlen)
{
	// Get list of all xpaths with "addmaterial" and "fractionmass"
	char filter[512];
	sprintf(filter,"//materials/composite[@name='%s']/addmaterial/fractionmass[@fraction]", name.c_str());
	vector<string> xpaths;
	jgeom->GetXPaths(xpaths, JGeometry::attr_level_all, filter);
	
	// Loop over components of this composite
_DBG_<<"Components for compsite "<<name<<endl;
	for(unsigned int i=0; i<xpaths.size(); i++){
		// Get component material name
		string::size_type start_pos = xpaths[i].find("@material=", 0);
		start_pos = xpaths[i].find("'", start_pos);
		string::size_type end_pos = xpaths[i].find("'", start_pos+1);
		if(end_pos==string::npos)continue;
		string mat_name = xpaths[i].substr(start_pos+1, end_pos-(start_pos+1));

		// Get component mass fraction
		start_pos = xpaths[i].find("fractionmass[", 0);
		start_pos = xpaths[i].find("@fraction=", start_pos);
		start_pos = xpaths[i].find("'", start_pos);
		end_pos = xpaths[i].find("'", start_pos+1);
		if(end_pos==string::npos)continue;
		string mat_frac_str = xpaths[i].substr(start_pos+1, end_pos-(start_pos+1));
		double fractionmass = atof(mat_frac_str.c_str());

		_DBG_<<"   "<<xpaths[i]<<"  fractionmass="<<fractionmass<<endl;
	}
	
	return true;
}

//---------------------------------
// GetTraversedMaterialsZ
//---------------------------------
//void DGeometry::GetTraversedMaterialsZ(double q, const DVector3 &pos, const DVector3 &mom, double z_end, vector<DMaterialStep> &materialsteps)
//{
	/// Find the materials traversed by a particle swimming from a specific
	/// position with a specific momentum through the magnetic field until
	/// it reaches a specific z-location. Energy loss is not included in
	/// the swimming since this method itself is assumed to be one of the
	/// primary means of determining energy loss. As such, one should not
	/// pass in a value of z_end that is far from pos.
	///
	/// The vector materialsteps will be filled with DMaterialStep objects
	/// corresponding to each of the materials traversed.
	///
	/// The real work here is done by the DMaterialStepper class
	
//}


//---------------------------------
// GetFDCWires
//---------------------------------
bool DGeometry::GetFDCWires(vector<vector<DFDCWire *> >&fdcwires) const{
  // Get geometrical information from database
  vector<double>z_wires;
  vector<double>stereo_angles;

  if(!GetFDCZ(z_wires)) return false;
  if(!GetFDCStereo(stereo_angles)) return false;

  for(int layer=1; layer<=FDC_NUM_LAYERS; layer++){
    double angle=-stereo_angles[layer-1]*M_PI/180.;
    vector<DFDCWire *>temp;
    for(int wire=1; wire<=WIRES_PER_PLANE; wire++){
      
      DFDCWire *w = new DFDCWire;
      w->layer = layer;
      w->wire = wire;
      w->angle = angle;

      // find coordinates of center of wire in rotated system
      float u = U_OF_WIRE_ZERO + WIRE_SPACING*(float)(wire-1);
      
      // Rotate coordinates into lab system and set the wire's origin
      // Note that the FDC measures "angle" such that angle=0
      // corresponds to the anode wire in the vertical direction
      // (i.e. at phi=90 degrees).
      float x = u*sin(angle + M_PI/2.0);
      float y = u*cos(angle + M_PI/2.0);
      w->origin.SetXYZ(x,y,z_wires[layer-1]);

      // Length of wire is set by active radius
      w->L = 2.0*sqrt(pow(FDC_ACTIVE_RADIUS,2.0) - u*u);
			
      // Set directions of wire's coordinate system with "udir"
      // along wire.
      w->udir.SetXYZ(sin(angle),cos(angle),0.0);
			
      // "s" points in direction from beamline to midpoint of
      // wire. This happens to be the same direction as "origin"
      w->sdir = w->origin;
      w->sdir.SetMag(1.0);
      
      w->tdir = w->udir.Cross(w->sdir);
      w->tdir.SetMag(1.0); // This isn't really needed
      temp.push_back(w);
    }
    fdcwires.push_back(temp);
  }

  return true;
}

//---------------------------------
// GetFDCZ
//---------------------------------
bool DGeometry::GetFDCZ(vector<double> &z_wires) const
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
	vector<double> forwardDC_module[4];
	vector<double> forwardDC_chamber[4][6];
	
	if(!Get("//section/composition/posXYZ[@volume='ForwardDC']/@X_Y_Z", ForwardDC)) return false;
	if(!Get("//composition[@name='ForwardDC']/posXYZ[@volume='forwardDC']/@X_Y_Z", forwardDC)) return false;
	if(!Get("//posXYZ[@volume='forwardDC_package_1']/@X_Y_Z", forwardDC_package[0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_package_2']/@X_Y_Z", forwardDC_package[1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_package_3']/@X_Y_Z", forwardDC_package[2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_package_4']/@X_Y_Z", forwardDC_package[3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_1']/@X_Y_Z", forwardDC_module[0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_2']/@X_Y_Z", forwardDC_module[1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_3']/@X_Y_Z", forwardDC_module[2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_4']/@X_Y_Z", forwardDC_module[3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[0][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[0][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[0][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='4']", forwardDC_chamber[0][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='5']", forwardDC_chamber[0][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='6']", forwardDC_chamber[0][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[1][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[1][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[1][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='4']", forwardDC_chamber[1][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='5']", forwardDC_chamber[1][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@X_Y_Z/layer[@value='6']", forwardDC_chamber[1][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[2][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[2][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[2][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='4']", forwardDC_chamber[2][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='5']", forwardDC_chamber[2][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@X_Y_Z/layer[@value='6']", forwardDC_chamber[2][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='1']", forwardDC_chamber[3][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='2']", forwardDC_chamber[3][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='3']", forwardDC_chamber[3][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='4']", forwardDC_chamber[3][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='5']", forwardDC_chamber[3][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@X_Y_Z/layer[@value='6']", forwardDC_chamber[3][5])) return false;
		
	// Offset due to global FDC envelopes
	double zfdc = ForwardDC[2] + forwardDC[2];
	
	// Loop over packages
	for(int package=1; package<=4; package++){
		double z_package = forwardDC_package[package-1][2];
		
		// Each "package" has 1 "module" which is currently positioned at 0,0,0 but
		// that could in principle change so we add the module z-offset
		double z_module = forwardDC_module[package-1][2];

		// Loop over chambers in this module
		for(int chamber=1; chamber<=6; chamber++){
			double z_chamber = forwardDC_chamber[package-1][chamber-1][2];
			
			double z = zfdc + z_package + z_module + z_chamber;				
			z_wires.push_back(z);
		}
	}
	
	return true;
}

//---------------------------------
// GetFDCStereo
//---------------------------------
bool DGeometry::GetFDCStereo(vector<double> &stereo_angles) const
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

	vector<double> forwardDC_module[4];
	vector<double> forwardDC_chamber[4][6];
	
	if(!Get("//posXYZ[@volume='forwardDC_module_1']/@X_Y_Z", forwardDC_module[0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_2']/@X_Y_Z", forwardDC_module[1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_3']/@X_Y_Z", forwardDC_module[2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_module_4']/@X_Y_Z", forwardDC_module[3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='1']", forwardDC_chamber[0][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='2']", forwardDC_chamber[0][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='3']", forwardDC_chamber[0][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='4']", forwardDC_chamber[0][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='5']", forwardDC_chamber[0][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_1']/@rot/layer[@value='6']", forwardDC_chamber[0][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='1']", forwardDC_chamber[1][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='2']", forwardDC_chamber[1][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='3']", forwardDC_chamber[1][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='4']", forwardDC_chamber[1][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='5']", forwardDC_chamber[1][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_2']/@rot/layer[@value='6']", forwardDC_chamber[1][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='1']", forwardDC_chamber[2][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='2']", forwardDC_chamber[2][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='3']", forwardDC_chamber[2][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='4']", forwardDC_chamber[2][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='5']", forwardDC_chamber[2][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_3']/@rot/layer[@value='6']", forwardDC_chamber[2][5])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='1']", forwardDC_chamber[3][0])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='2']", forwardDC_chamber[3][1])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='3']", forwardDC_chamber[3][2])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='4']", forwardDC_chamber[3][3])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='5']", forwardDC_chamber[3][4])) return false;
	if(!Get("//posXYZ[@volume='forwardDC_chamber_4']/@rot/layer[@value='6']", forwardDC_chamber[3][5])) return false;
	
	// Loop over packages
	for(int package=1; package<=4; package++){
		
		// Loop over chambers
		for(int chamber=1; chamber<=6; chamber++){
			stereo_angles.push_back(forwardDC_chamber[package-1][chamber-1][2]);
		}
	}

	return true;
}

//---------------------------------
// GetFDCRmin
//---------------------------------
bool DGeometry::GetFDCRmin(vector<double> &rmin_packages) const
{
	vector<double> FDA[4];
	
	if(!Get("//section[@name='ForwardDC']/tubs[@name='FDA1']/@Rio_Z", FDA[0])) return false;
	if(!Get("//section[@name='ForwardDC']/tubs[@name='FDA2']/@Rio_Z", FDA[1])) return false;
	if(!Get("//section[@name='ForwardDC']/tubs[@name='FDA3']/@Rio_Z", FDA[2])) return false;
	if(!Get("//section[@name='ForwardDC']/tubs[@name='FDA4']/@Rio_Z", FDA[3])) return false;

	rmin_packages.push_back(FDA[0][0]);
	rmin_packages.push_back(FDA[1][0]);
	rmin_packages.push_back(FDA[2][0]);
	rmin_packages.push_back(FDA[3][0]);

	return true;
}

//---------------------------------
// GetFDCRmax
//---------------------------------
bool DGeometry::GetFDCRmax(double &rmax_active_fdc) const
{
	// We assume that all packages have the same outer radius of the
	// active area.
	vector<double> FDA1;

	bool good = Get("//section[@name='ForwardDC']/tubs[@name='FDA1']/@Rio_Z", FDA1);
	
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
bool DGeometry::GetCDCOption(string &cdc_option) const
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
bool DGeometry::GetCDCCenterZ(double &cdc_center_z) const
{

	return false;
}

//---------------------------------
// GetCDCAxialLength
//---------------------------------
bool DGeometry::GetCDCAxialLength(double &cdc_axial_length) const
{
	vector<double> Rio_Z;
	bool good = Get("//section[@name='CentralDC']/tubs[@name='STRW']/@Rio_Z", Rio_Z);
	cdc_axial_length = Rio_Z[2];

	if(!good){
		_DBG_<<"Unable to retrieve CDC axial wire length"<<endl;
	}

	return false;
}

//---------------------------------
// GetCDCStereo
//---------------------------------
bool DGeometry::GetCDCStereo(vector<double> &cdc_stereo) const
{

	return false;
}

//---------------------------------
// GetCDCRmid
//---------------------------------
bool DGeometry::GetCDCRmid(vector<double> &cdc_rmid) const
{

	return false;
}

//---------------------------------
// GetCDCNwires
//---------------------------------
bool DGeometry::GetCDCNwires(vector<int> &cdc_nwires) const
{

	return false;
}


//---------------------------------
// GetCDCEndplate
//---------------------------------
bool DGeometry::GetCDCEndplate(double &z,double &dz,double &rmin,double &rmax)
  const{

  vector<double>cdc_origin;
  vector<double>cdc_center;
  vector<double>cdc_endplate_pos;
  vector<double>cdc_endplate_dim;
  
  if(!Get("//posXYZ[@volume='CentralDC'/@X_Y_Z",cdc_origin)) return false;
  if(!Get("//posXYZ[@volume='centralDC_option-1']/@X_Y_Z",cdc_center)) return false;
  if(!Get("//posXYZ[@volume='CDPD']/@X_Y_Z",cdc_endplate_pos)) return false;
  if(!Get("//tubs[@name='CDPD']/@Rio_Z",cdc_endplate_dim)) return false;
  
	if(cdc_origin.size()<3){
  		_DBG_<<"cdc_origin.size()<3 !"<<endl;
		return false;
	}
	if(cdc_center.size()<3){
  		_DBG_<<"cdc_center.size()<3 !"<<endl;
		return false;
	}
	if(cdc_endplate_pos.size()<3){
  		_DBG_<<"cdc_endplate_pos.size()<3 !"<<endl;
		return false;
	}
	if(cdc_endplate_dim.size()<3){
  		_DBG_<<"cdc_endplate_dim.size()<3 !"<<endl;
		return false;
	}
  
  z=cdc_origin[2]+cdc_center[2]+cdc_endplate_pos[2];
  dz=cdc_endplate_dim[2];
  rmin=cdc_endplate_dim[0];
  rmax=cdc_endplate_dim[1];

  return true;
}
//---------------------------------
// GetBCALRmin
//---------------------------------
bool DGeometry::GetBCALRmin(double &bcal_rmin) const
{

	return false;
}

//---------------------------------
// GetBCALNmodules
//---------------------------------
bool DGeometry::GetBCALNmodules(unsigned int &bcal_nmodules) const
{

	return false;
}

//---------------------------------
// GetBCALCenterZ
//---------------------------------
bool DGeometry::GetBCALCenterZ(double &bcal_center_z) const
{

	return false;
}

//---------------------------------
// GetBCALLength
//---------------------------------
bool DGeometry::GetBCALLength(double &bcal_length) const
{

	return false;
}

//---------------------------------
// GetBCALDepth
//---------------------------------
bool DGeometry::GetBCALDepth(double &bcal_depth) const
{

	return false;
}

//---------------------------------
// GetFCALZ
//---------------------------------
bool DGeometry::GetFCALZ(double &z_fcal) const
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
bool DGeometry::GetTOFZ(vector<double> &z_tof) const
{
	vector<double> ForwardTOF;
	vector<double> forwardTOF[2];
	vector<double> FTOC;
       
	if(!Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z", ForwardTOF)) return false;
	if(!Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", forwardTOF[0])) return false;
	if(!Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", forwardTOF[1])) return false;
	if(!Get("//box[@name='FTOC' and sensitive='true']/@X_Y_Z", FTOC)) return false;
	
	z_tof.push_back(ForwardTOF[2] + forwardTOF[0][2] - FTOC[2]/2.0);
	z_tof.push_back(ForwardTOF[2] + forwardTOF[1][2] - FTOC[2]/2.0);

	return true;
}

//---------------------------------
// GetTargetZ
//---------------------------------
bool DGeometry::GetTargetZ(double &z_target) const
{
	vector<double> xyz_vessel;
	vector<double> xyz_target;
	vector<double> xyz_detector;
       
	z_target=0.;

	if(!Get("//composition[@name='targetVessel']/posXYZ[@volume='targetTube']/@X_Y_Z", xyz_vessel)) return false;
	if(!Get("//composition[@name='Target']/posXYZ[@volume='targetVessel']/@X_Y_Z", xyz_target)) return false;
	if(!Get("//composition[@name='GlueXdetector']/posXYZ[@volume='Target']/@X_Y_Z", xyz_detector)) return false;
	
	z_target = xyz_vessel[2] + xyz_target[2] + xyz_detector[2];

	return true;
}

//---------------------------------
// GetTargetLength
//---------------------------------
bool DGeometry::GetTargetLength(double &target_length) const
{
 	vector<double> Rio_Z;
	bool good = Get("//section[@name='Target']/tubs[@name='LIH2']/@Rio_Z", Rio_Z);

	target_length = good ? Rio_Z[2]:0.0;
 
	return good;
}

