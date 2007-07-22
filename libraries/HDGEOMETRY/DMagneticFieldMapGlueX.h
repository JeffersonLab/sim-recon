// $Id$
//
// This originally defined the DMagneticFieldMap class but was converted 
// over to be a specific implementation of a virtual class of the same name.
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifndef __DMAGNETICFIELDMAPGLUEX_H__
#define __DMAGNETICFIELDMAPGLUEX_H__

/*! \file DMagneticFieldMapGlueX.h
 *  \brief A class that holds the magnetic field map for the detector.
 */
#include <cmath>
#include <fstream>


#include "DMagneticFieldMap.h"
 


class DMagneticFieldMapGlueX:public DMagneticFieldMap
{

	public:

		typedef struct{
			float x,y,z,Bx,By,Bz;
		}DBfieldPoint_t;

		typedef struct{
			double dBrdr, dBrdz;
			double dBzdr, dBzdz;
		}DGradient_t;

		DMagneticFieldMapGlueX();	/*!< Constructor */
		~DMagneticFieldMapGlueX();	/*!< Destructor  */

		int GetIndices(double x, double y, double z, int &index_r, int &index_z) const;
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
    
	 	void GetTable(const DBfieldPoint_t* &Bmap, int &Npoints);
		
	private:


		DBfieldPoint_t *Bmap;	/*!< serialized array for Mag Field values  */
		DGradient_t *GradMap;	/*!< array of gradients (numerically calculated) */
		int Npoints;
		double rMin;   /*!< minimum double radius coordinate  */
		double rMax;   /*!< maximum double radius coordinate  */
		double zMin;   /*!< minimum double z coordinate  */
		double zMax;   /*!< maximum double z coordinate  */
		int    rDim;   /*!< radial dimension of internal data structure  */
		int    zDim;   /*!< z dimension of internal data structure  */
		int    phiDim; /*!< phi dimension of internal data structure */
		float  BZ_CONST;/*!< Use this (in tesla) for homogeneous field if gt -100 */
		double dr;		/*!< change in r between map points */
		double dz;		/*!< change in z between map points */

		DBfieldPoint_t null_point;
		float BMAP_Z_OFFSET;
		float BZ_AVG_Z;
};

#endif
    

                          
                                  


    






