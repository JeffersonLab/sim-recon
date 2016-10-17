// DFCALShower member functions

#include "DFCALShower.h"

DFCALShower::DFCALShower():ExyztCovariance(5)
{
   fEnergy = 0.;
   fTime = 0.;
   fPosition.SetXYZ(0., 0., 0.) ;
}

DFCALShower::~DFCALShower()
{
}

//
 
void DFCALShower::setEnergy(const double energy) 
{

        fEnergy = energy;

}

void DFCALShower::setTime(const double time) 
{

        fTime = time;

}

void DFCALShower::setPosition(const DVector3 aPosition ) 
{

        fPosition = aPosition;

}


