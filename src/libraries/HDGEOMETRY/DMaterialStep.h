// $Id$
//
//    File: DMaterialStep.h
// Created: Fri Apr 25 15:44:18 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#ifndef _DMaterialStep_
#define _DMaterialStep_

#include <HDGEOMETRY/DMaterial.h>
#include <DVector3.h>


/// The DMaterialStep class is just a container for holding information
/// about what material was traversed by a track while swimming a small
/// segment. This class does not actually swim the particle or calculate
/// anything, see DGeometry::GetTraversedMaterialsZ() for that.

class DMaterialStep{
	public:
		DMaterialStep(){}
		DMaterialStep(const DMaterial *material
			, const DVector3 &startpos, const DVector3 &endpos
			, const DVector3 &startmom, const DVector3 &endmom)
		{
				this->material = material;
				this->startpos = startpos;
				this->endpos = endpos;
				this->startmom = startmom;
				this->endmom = endmom;
		}
		virtual ~DMaterialStep(){}

		const DMaterial* GetDMaterial(void) const {return material;}
		const DVector3& GetStartPosition(void) const {return startpos;}
		const DVector3& GetEndPosition(void) const {return endpos;}
		const DVector3& GetStartMomentum(void) const {return startmom;}
		const DVector3& GetEndMomentum(void) const {return endmom;}
		double GetDistanceAlongPath(void) const {return s;}
		
		void SetDMaterial(const DMaterial *material){this->material = material;}
		void SetStartPosition(const DVector3 &startpos){this->startpos = startpos;}
		void SetEndPosition(const DVector3 &endpos){this->endpos = endpos;}
		void SetStartMomentum(const DVector3 &startmom){this->startmom = startmom;}
		void SetEndMomentum(const DVector3 &endmom){this->endmom = endmom;}
		
	private:
		const DMaterial *material;
		DVector3 startpos;
		DVector3 endpos;
		DVector3 startmom;
		DVector3 endmom;
		double s;
};

#endif // _DMaterialStep_

