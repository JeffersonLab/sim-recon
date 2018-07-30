// DFCALShower member functions

#include "DFCALShower.h"

DFCALShower::DFCALShower():ExyztCovariance(5)
{
  fEnergy = 0.;
  fTime = 0.;
  fPosition.SetXYZ(0., 0., 0.);
  fTimeTr = 1E3;
  fDocaTr = 1E3;
  fSumU = 0;
  fSumV = 0;
  fE9E25 = 0;
  fE1E9 = 0;
  iNumBlocks = 0;
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

void DFCALShower::setPosition(const DVector3& aPosition ) 
{
  fPosition = aPosition;
}

void DFCALShower::setDocaTrack( const double docaTrack ){ fDocaTr = docaTrack; }
void DFCALShower::setTimeTrack( const double tTrack ){ fTimeTr = tTrack; }
void DFCALShower::setSumU( const double sumU ){ fSumU = sumU; }
void DFCALShower::setSumV( const double sumV ){ fSumV = sumV; }
void DFCALShower::setE9E25( const double e9e25 ){ fE9E25 = e9e25; }
void DFCALShower::setE1E9( const double e1e9 ){ fE1E9 = e1e9; }
void DFCALShower::setNumBlocks( const int numBlocks ){ iNumBlocks = numBlocks; }


