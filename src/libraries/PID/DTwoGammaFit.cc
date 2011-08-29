#include "DTwoGammaFit.h"

DTwoGammaFit::DTwoGammaFit()
{
    for ( int i = 0; i < 2; i++ ) {
        DVector3 null( 0. ,0. ,0. );
        fChildFits[i].setPosition( null );
        fChildFits[i].setMomentum( null );
        fChildFits[i].setMass(0.);
        fIDs[i] = 0;
    }
    fChi2 = 0;
    fProb = 0;
// pulls : px,py,pz for two photons ?
    for ( int i = 0; i < 6; i++) {
        fPulls[0] = -1000;
    }
}

DTwoGammaFit::DTwoGammaFit(const oid_t id) : 
    DKinematicData(id)
{
    for ( int i = 0; i < 2; i++ ) {
        DVector3 null( 0. ,0. ,0. );
        fChildFits[i].setPosition( null );
        fChildFits[i].setMomentum( null );
        fChildFits[i].setMass(0.);
        fIDs[i] = 0;
    }
    fChi2 = 0;
    fProb = 0;
// pulls : px,py,pz for two photons ?
    for ( int i = 0; i < 6; i++) {
        fPulls[0] = -1000;
    }
}

DTwoGammaFit::~DTwoGammaFit()
{
}
