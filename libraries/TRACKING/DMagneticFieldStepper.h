


#ifndef __DMAGNETICFIELDSTEPPER_H__
#define __DMAGNETICFIELDSTEPPER_H__

#include <math.h>

#include <TVector3.h>
#include "JANA/jerror.h"

#include "DMagneticFieldMap.h"


/// DMagneticFieldStepper class
///
/// This class will step a particle track through a magnetic
/// field. It has methods to find the point on the track which
/// comes closest to a specified point in space.


class DMagneticFieldStepper
{
	public:

		DMagneticFieldStepper(const DMagneticFieldMap *map);
		DMagneticFieldStepper(const DMagneticFieldMap *map, double q, TVector3 *x, TVector3 *p);
		~DMagneticFieldStepper();
	
		jerror_t SetStartingParams(double q, const TVector3 *x, const TVector3 *p);
		jerror_t SetMagneticFieldMap(const DMagneticFieldMap *map);
		jerror_t SetStepSize(double step);
		double Step(TVector3 *newpos=NULL);
		const DBfieldPoint_t* GetDBfieldPoint(void);
		void GetDirs(TVector3 &xdir, TVector3 &ydir, TVector3 &zdir);
		void GetMomentum(TVector3 &mom){mom = this->mom;}
		void GetPosition(TVector3 &pos){pos = this->pos;}
		double GetCharge(void){return q;}
		void GetPosMom(TVector3 &pos, TVector3 &mom){pos=this->pos; mom=this->mom;}
		inline double GetRo(void){return fabs(Ro);}
		inline double Getdz_dphi(void){return Ro*mom.Dot(zdir)/mom.Dot(ydir);}
		
		inline double GetStepSize(void) const{return stepsize;}
	
	private:
		const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
		double stepsize;		///< distance(cm) to move particle when Step() is called
		double q;				///< electric charge in units of e
		TVector3 pos;			///< current position of particle
		TVector3 mom;			///< current location of particle
		TVector3 start_pos;	///< starting position of track
		TVector3 start_mom;	///< starting momentum of track
		double Ro;
		
		TVector3 xdir, ydir, zdir;
		
		void CalcDirs(void);
		void CalcDirs(TVector3 *B);
};

#endif // __DMAGNETICFIELDSTEPPER_H__
