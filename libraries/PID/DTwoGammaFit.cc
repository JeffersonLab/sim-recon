// DPhoton member functions

#include "DTwoGammaFit.h"


DTwoGammaFit::DTwoGammaFit()
{
    DVector3 null( 0. ,0. ,0. );
    fChildFit[0].setPosition( null );
    fChildFit[0].setMomentum( null );
    fChildFit[0].setMass(0);
    fChildFit[1].setPosition( null );
    fChildFit[1].setMomentum( null );
    fChildFit[1].setMass(0);
    fChi2 = 0;
    fProb = 0;
    fPulls[0] = 100;
    fPulls[1] = 100;
}

DTwoGammaFit::~DTwoGammaFit()
{
}

/* Set data of fitted children
void DTwoGammaFit::setChildFit(const DKinematicData& aChildFit, const int i)
{
     fChildFit[i] = aChildFit;
}

// Set pulls from DKinFit
void DTwoGammaFit::setPull(const double aPull, const int i)
{
     fPulls[i] = aPull;
}

// Set chi2 from DKinFit
void DTwoGammaFit::setChi2(const double aChi2)
{
     fChi2 = aChi2;
}

// Set confidence from DKinFit
void DTwoGammaFit::setProb(const double aProb)
{
     fProb = aProb;
}

*/
