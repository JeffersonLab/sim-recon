


#ifndef __DMAGNETICFIELDSTEPPER_H__
#define __DMAGNETICFIELDSTEPPER_H__

#include <math.h>

#include <DVector3.h>
#include "JANA/jerror.h"

#include "HDGEOMETRY/DMagneticFieldMap.h"


/// DMagneticFieldStepper class
///
/// This class will step a particle track through a magnetic
/// field. It has methods to find the point on the track which
/// comes closest to a specified point in space.


class DMagneticFieldStepper
{
	public:

		DMagneticFieldStepper(const DMagneticFieldMap *map, double q=1.0);
		DMagneticFieldStepper(const DMagneticFieldMap *map, double q, const DVector3 *x, const DVector3 *p);
		~DMagneticFieldStepper();
	
		jerror_t SetStartingParams(double q, const DVector3 *x, const DVector3 *p);
		jerror_t SetMagneticFieldMap(const DMagneticFieldMap *map);
		jerror_t SetStepSize(double step);
		void SetCharge(double q){this->q = q;}
		double Step(DVector3 *newpos=NULL, double stepsize=0.0);
		void GetDirs(DVector3 &xdir, DVector3 &ydir, DVector3 &zdir);
		void GetBField(DVector3 &B){B = this->B;}
		void GetMomentum(DVector3 &mom){mom = this->mom;}
		void GetPosition(DVector3 &pos){pos = this->pos;}
		double GetCharge(void){return q;}
		void GetPosMom(DVector3 &pos, DVector3 &mom){pos=this->pos; mom=this->mom;}
		bool SwimToPlane(DVector3 &pos, DVector3 &mom, const DVector3 &origin, const DVector3 &norm, double *pathlen=NULL);
		bool DistToPlane(DVector3 &pos, const DVector3 &origin, const DVector3 &norm);
		bool SwimToRadius(DVector3 &pos, DVector3 &mom, double R, double *pathlen=NULL);
		bool DistToRadius(DVector3 &pos, double R);

		inline double GetRo(void){return fabs(Ro);}
		inline double Getdz_dphi(void){return Ro*mom.Dot(zdir)/mom.Dot(ydir);}
		inline double GetStepSize(void) const{return stepsize;}
	
	private:
		const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
		double stepsize;		///< maximum distance(cm) to move particle when Step() is called
		double last_stepsize;///< stepsize (cm) used for last step
		double q;				///< electric charge in units of e
		DVector3 pos;			///< current position of particle
		DVector3 mom;			///< current location of particle
		DVector3 start_pos;	///< starting position of track
		DVector3 start_mom;	///< starting momentum of track
		DVector3 B;
		double Ro, Rp;
		double cos_theta, sin_theta;
		
		DVector3 xdir, ydir, zdir;
		
		void CalcDirs(double *Bvals=NULL);
};

#endif // __DMAGNETICFIELDSTEPPER_H__
