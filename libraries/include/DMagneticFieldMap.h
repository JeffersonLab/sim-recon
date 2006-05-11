//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifndef __DMAGNETICFIELDMAP_H__
#define __DMAGNETICFIELDMAP_H__

/*! \file DMagneticFieldMap.h
 *  \brief A class that holds the magnetic field map for the detector.
 */

#include <fstream>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c

 
typedef struct{
	float x,y,z,Bx,By,Bz;
}DBfieldPoint_t;


class DMagneticFieldMap
{

	public:

		DMagneticFieldMap();	/*!< Constructor */
		~DMagneticFieldMap();	/*!< Destructor  */

		const DBfieldPoint_t* getQuick(double x, double y, double z) const;
		
		double Bz_avg(double x, double y, double x0, double y0, double delta_phi) const;
		double Bz_avg(double x, double y, double z, double x0, double y0, double theta, double zmax) const;
    
	private:

		DBfieldPoint_t *Bmap;	/*!< serialized array for Mag Field values  */
		int Npoints;
		double rMin;   /*!< minimum double radius coordinate  */
		double rMax;   /*!< maximum double radius coordinate  */
		double zMin;   /*!< minimum double z coordinate  */
		double zMax;   /*!< maximum double z coordinate  */
		int    rDim;   /*!< radial dimension of internal data structure  */
		int    zDim;   /*!< z dimension of internal data structure  */
		int    phiDim; /*!< phi dimension of internal data structure */
		float  constBz;/*!< Use this (in tesla) for homogeneous field if gt -100 */

		DBfieldPoint_t null_point;
};

#endif
    

                          
                                  


    






