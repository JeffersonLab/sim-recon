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
 
typedef struct{
	float x,y,z,Bx,By,Bz;
}DBfieldPoint_t;


class DMagneticFieldMapGlueX:public DMagneticFieldMap
{

	public:

		DMagneticFieldMapGlueX();	/*!< Constructor */
		~DMagneticFieldMapGlueX();	/*!< Destructor  */

		const DBfieldPoint_t* getQuick(double x, double y, double z) const;
		int GetIndices(double x, double y, double z, int &index_r, int &index_z) const;
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
		
		double Bz_avg(double x, double y, double x0, double y0, double delta_phi) const;
		double Bz_avg(double x, double y, double z, double x0, double y0, double theta, double zmax) const;
    
	 	void GetTable(const DBfieldPoint_t* &Bmap, int &Npoints);
		
		/// Set this map to use a constant B-field in the z- direction. Bz is in Tesla.
		void SetConstField(float Bz){BZ_CONST = Bz;}
		bool isConst(void){return fabs(BZ_CONST)<100.0;}
	 
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
		float  BZ_CONST;/*!< Use this (in tesla) for homogeneous field if gt -100 */
		double dr;		/*!< change in r between map points */
		double dz;		/*!< change in z between map points */

		DBfieldPoint_t null_point;
		float BMAP_Z_OFFSET;
		float BZ_AVG_Z;
};

#endif
    

                          
                                  


    






