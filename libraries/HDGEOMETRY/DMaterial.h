// $Id$
//
//    File: DMaterial.h
// Created: Fri Apr 25 15:44:12 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#ifndef _DMaterial_
#define _DMaterial_

#include <string>
#include <string.h>
using std::string;

#include <particleType.h>


/// The DMaterial class holds information on a single material type. The
/// main purpose is to hold information needed to estimate energy loss and
/// multiple scattering for particles being swum through the detector.
/// The GetdEdx() method can be used to calculate dE/dx through this
/// material for a given particle type and mometum using the Bethe
/// formula. 

class DMaterial{
	public:
		DMaterial(string &name, double A, double Z, double density, double radlen);
		virtual ~DMaterial(){}
	
		double GetRadiationLength(void) const {return kXo;}
		double GetDensity(void) const {return kdensity;}
		double GetAtomicWeight(void) const {return kA;}
		double GetAtomicNumber(void) const {return kZ;}
		double GetdEdx(Particle_t ptype, double kp);
		const string& GetName(void) const {return kname;}

		double Xo(void) const {return kXo;}
		double rho(void) const {return kdensity;}
		double A(void) const {return kA;}
		double Z(void) const {return kZ;}
		
		string toString(void);

	private:
		DMaterial(); ///< Prevent use of trivial constructor
	
		double kXo; ///< radiation length of this material
		double kdensity; ///< denisty of material in g/cm^3
		double kA; ///< Effective atomic weight
		double kZ; ///< Effective atomic number
		string kname;

		Particle_t last_ptype;
		double last_p;
		double last_dEdx;
		
		void GetMaterials(void);
};

#endif // _DMaterial_

