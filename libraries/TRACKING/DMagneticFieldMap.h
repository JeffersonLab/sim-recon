//-------------------------------------------------------------------------
//
//      Author:     Brad Futch 
//      Date:       Wed 01 Dec 2004 02:04:52 
//
//      Copyright (C) 2004 Brad Futch 
//      
//      This program is free software; you can redistribute it and/or
//      modify it under the terms of the GNU General Public License
//      as published by the Free Software Foundation; either version 2
//      of the License, or (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
//      MA  02111-1307, USA.
//
//
//-------------------------------------------------------------------------

#ifndef __DMAGNETICFIELDMAP_H__
#define __DMAGNETICFIELDMAP_H__

/*! \file DMagneticFieldMap.h
 *  \brief A class that holds the magnetic field map for the detector.
 */

#include <fstream>

using std::ifstream;

typedef struct DVEC
{
  double x;
  double y;
  double z;
} D3Vector_t;

/*! DMagneticFiledMap class. 
 * 
 *  This class provides an interface to the map that is created for 
 *  the Solenoid used in this project.  The current implementation reads
 *  a solenoid table in the format described in the dsolenoid.table file.
 *
 *  NOTE:  There is an assumed radial symmetry, so there is actually only 
 *         an x and z entry in the table.  The x coordinate is assumed 
 *         to be equivalent to the radius.
 */

class DMagneticFieldMap
{

  public:

    DMagneticFieldMap(int rDim, int zDim);   /*!< Constructor */
    ~DMagneticFieldMap();                    /*!< Destructor  */

    /*!
     * A method to read a solenoid map from a file.
     * @param in - ifstream input
    */

    void          readMap(ifstream &in);

    /*! defgroup g1 getMagVec */
    //@{
    /*! 
     * These methods will perform a trilinear interpolation on the 
     * data in the map.  
     * */
    D3Vector_t    getMagVec(const D3Vector_t &vec);
    D3Vector_t    getMagVec(double x, double y, double z);
    D3Vector_t    getMagVec(double r, double z);
    //@}

    /*! defgroup g2 getQuick */
    //@{
    /*!
     * These methods will get the nearst index to to the requested position.
     * This can also be described as 0th order interpolation.
     */
    D3Vector_t    getQuick(const D3Vector_t &vec);
    D3Vector_t    getQuick(double x, double y, double z);
    D3Vector_t    getQuick(double r, double z);
    //@}
    

    /*! Set the Debug state.
     * @param val - state to set
     */
    void          setDebug(const bool val);

    /*!
     * A utility function to return nearest indices.
     * @param r - a double radius
     * @param z - a double z value
     * @param ind - a 2D int vector of return indices
     * @param rho - amount to interpolate in the radial direction
     * @param zeta - amount to interpolte in the z direction
     *
     *  This will be used for the trilinear interpolation in order 
     *  to fetch the closest indices from which to interpolate.
     */
    
    void          getInds(const double &r, const double &z, int ind[2],
                          double &rho, double &zeta);
    /*!
     *  This utility funciton based on the assumption of radial symmetry.
     *  @param r - integer radius index
     *  @param z - integer z index
     *
     *  This will return the base index for the value.  The actual xyz values
     *  will be:
     *
     *  x = map[ind]
     *  y = map[ind+1]
     *  z = map[ind+2]
     */
    
    int           serialize(int r, int z);


  private:

    double *map;   /*!< serialized array for Mag Field values  */
    double rMin;   /*!< minimum double radius coordinate  */
    double rMax;   /*!< maximum double radius coordinate  */
    double zMin;   /*!< minimum double z coordinate  */
    double zMax;   /*!< maximum double z coordinate  */
    int    rDim;   /*!< radial dimension of internal data structure  */
    int    zDim;   /*!< z dimension of internal data structure  */
    bool   DEBUG;  /*!< if true output debug informations  */

    

};

#endif
    

                          
                                  


    






