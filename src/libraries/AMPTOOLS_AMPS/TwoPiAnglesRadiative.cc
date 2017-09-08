
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAnglesRadiative.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

#include <CobremsGeneration.hh>

TwoPiAnglesRadiative::TwoPiAnglesRadiative( const vector< string >& args ) :
    UserAmplitude< TwoPiAnglesRadiative >( args )
{
    assert( args.size() == 11 );

    rho000  = AmpParameter( args[0] );
    rho100  = AmpParameter( args[1] );
    rho1m10 = AmpParameter( args[2] );

    rho111  = AmpParameter( args[3] );
    rho001  = AmpParameter( args[4] );
    rho101  = AmpParameter( args[5] );
    rho1m11 = AmpParameter( args[6] );

    rho102  = AmpParameter( args[7] );
    rho1m12 = AmpParameter( args[8] );

    polAngle = AmpParameter( args[9] );

    polFraction = atof(args[10].c_str());

    // need to register any free parameters so the framework knows about them
    registerParameter( rho000 );
    registerParameter( rho100 );
    registerParameter( rho1m10 );

    registerParameter( rho111 );
    registerParameter( rho001 );
    registerParameter( rho101 );
    registerParameter( rho1m11 );

    registerParameter( rho102 );
    registerParameter( rho1m12 );

    registerParameter( polAngle );

    // Initialize coherent brem table
    // Do this over the full range since we will be using this as a lookup
    float Emax  = 12.0;
    float Epeak = 9.0;
    float Elow  = 0.139*2;
    float Ehigh = 12.0;

    int doPolFlux=0;  // want total flux (1 for polarized flux)
    float emitmr=10.e-9; // electron beam emittance
    float radt=20.e-6; // radiator thickness in m
    float collDiam=0.0034; // meters
    float Dist = 76.0; // meters
    CobremsGeneration cobrems(Emax, Epeak);
    cobrems.setBeamEmittance(emitmr);
    cobrems.setTargetThickness(radt);
    cobrems.setCollimatorDistance(Dist);
    cobrems.setCollimatorDiameter(collDiam);
    cobrems.setPolarizedFlag(doPolFlux);

    // Create histogram
    totalFlux_vs_E = new TH1D("totalFlux_vs_E", "Total Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
    polFlux_vs_E   = new TH1D("polFlux_vs_E", "Polarized Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
    polFrac_vs_E   = new TH1D("polFrac_vs_E", "Polarization Fraction vs. E_{#gamma}", 1000, Elow, Ehigh);

    // Fill totalFlux
    for(int i=1;i<=totalFlux_vs_E->GetNbinsX(); i++){
        double x = totalFlux_vs_E->GetBinCenter(i)/Emax;
        double y = 0;
        //if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
        y = cobrems.Rate_dNtdx(x);
        totalFlux_vs_E->SetBinContent(i, y);
    }

    doPolFlux=1;
    cobrems.setPolarizedFlag(doPolFlux);
    // Fill totalFlux
    for(int i=1;i<=polFlux_vs_E->GetNbinsX(); i++){
        double x = polFlux_vs_E->GetBinCenter(i)/Emax;
        double y = 0;
        //if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
        y = cobrems.Rate_dNcdx(x);
        polFlux_vs_E->SetBinContent(i, y);
    }

    polFrac_vs_E->Divide(polFlux_vs_E, totalFlux_vs_E);
}


complex< GDouble >
TwoPiAnglesRadiative::calcAmplitude( GDouble** pKin ) const {

    TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
    TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
    TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
    TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
    TLorentzVector target (0.0, 0.0, 0.0, 0.938272);

    TLorentzVector locOmegaP4 = p1 + p2;
    //Polarization plane:
    //Beam is in the lab z-direction
    //(FYI) Circularly polarized photon beam: Polarization rotates through the plane perpendicular to the direction of the photon: The XY Plane
    //The polarization vector is perpendicular to the direction of the photon
    //Linearly polarized photon beam: Polarization is confined to a plane along the direction of the photon
    //Plane defined by z-direction & some angle phi. Thus, polarization vector defined by phi.

    //Production CM frame: The center-of-mass frame of the production step. Here: g, p -> omega, p
    //In general, the beam energy is measured more accurately than the combination of all of the final-state particles
    //So define the production CM frame using the initial state
    TLorentzVector locInitialStateP4_KinFit = beam + target;
    TVector3 locBoostVector_ProdCM = -1.0*(locInitialStateP4_KinFit.BoostVector()); //negative due to coordinate system convention

    //boost beam & target proton to production CM frame
    TLorentzVector locBeamP4_ProdCM(beam);
    locBeamP4_ProdCM.Boost(locBoostVector_ProdCM);
    TLorentzVector locProtonP4_ProdCM(recoil);
    locProtonP4_ProdCM.Boost(locBoostVector_ProdCM);

    //Production plane:
    //The production plane is the plane containing the produced particles. Here: Defined by the proton and the omega
    //However, when you boost to the production CM frame, the production plane is no longer well defined: the particles are back-to-back
    //So, by convention, define the production plane in the production CM frame by the beam and the vector meson.

    //Production CM frame axes: "HELICITY SYSTEM"
    //The z-axis is defined as the direction of the meson (omega): z = Omega
    //The y-axis is defined by the vector cross product: y = Beam X Omega
    //The x-axis is defined by the vector cross product: x = y cross z
    //However, the proton momentum is in general better known than the omega momentum, so use it instead (they are back-to-back)
    //z = -1 * Proton
    //y = -1 * (Beam X Proton)
    //x = y cross z
    //Thus the production plane in the production frame is the XZ plane, and the normal vector is the Y-axis

    //Define production CM frame helicity axes
    TVector3 locHelicityZAxis_ProdCM = -1.0*locProtonP4_ProdCM.Vect().Unit();
    TVector3 locHelicityYAxis_ProdCM = -1.0*locBeamP4_ProdCM.Vect().Cross(locProtonP4_ProdCM.Vect()).Unit();
    TVector3 locHelicityXAxis_ProdCM = locHelicityYAxis_ProdCM.Cross(locHelicityZAxis_ProdCM).Unit();

    //Since the beam is in PARA configuration (Run 3185), the polarization vector is along the lab x-axis
    //Since the boost is in the z-direction, this vector is the same in the production CM frame
    TVector3 locPolUnit(cos(polAngle), sin(polAngle), 0.0);

    //In the production CM frame, locPHI is the angle between the polarization vector and the production plane
    double locCosPHI = locPolUnit.Dot(locHelicityYAxis_ProdCM.Cross(locBeamP4_ProdCM.Vect().Unit()));
    double locPHI = acos(locCosPHI); //reports phi between 0 and pi: sign ambiguity
    //Resolve the sign ambiguity
    double locSinPHI = locPolUnit.Dot(locHelicityYAxis_ProdCM);
    if(locSinPHI < 0.0) locPHI *= -1.0;

    //Now, we need the theta, phi angles between the omega decay plane and the production plane
    //The omega decay plane is defined by decay products in the omega CM frame
    //2 particles (vectors) define a plane.
    //However, to conserve momentum, the third particle cannot be out of that plane (so must also be in it)
    //So, use the pi+ and the pi- to define the plane (pi0 measurement has less resolution)
    //By the way, for rho decays, the theta & phi angles are those of the pi+ in the rho CM frame, with respect to the helicity axes

    //boost pi+/- to omega CM frame
    TVector3 locBoostVector_OmegaCM = -1.0*(locOmegaP4.BoostVector()); //negative due to coordinate system convention
    TLorentzVector locBeamP4_OmegaCM(beam);
    locBeamP4_OmegaCM.Boost(locBoostVector_OmegaCM);
    TLorentzVector locProtonP4_OmegaCM(recoil);
    locProtonP4_OmegaCM.Boost(locBoostVector_OmegaCM);
    TLorentzVector locPi0P4_OmegaCM(p1);
    locPi0P4_OmegaCM.Boost(locBoostVector_OmegaCM);
    TLorentzVector locGammaP4_OmegaCM(p2);
    locGammaP4_OmegaCM.Boost(locBoostVector_OmegaCM);

    //Define omega CM frame helicity axes
    //These are defined the same way as before, but with the boost, the direction of the x & y axes has changed
    TVector3 locHelicityZAxis_OmegaCM = -1.0*locProtonP4_OmegaCM.Vect().Unit();
    TVector3 locHelicityYAxis_OmegaCM = -1.0*locBeamP4_OmegaCM.Vect().Cross(locProtonP4_OmegaCM.Vect()).Unit();
    TVector3 locHelicityXAxis_OmegaCM = locHelicityYAxis_OmegaCM.Cross(locHelicityZAxis_OmegaCM).Unit();

    // Need the direction of the bachekor photon
    TVector3 locOmegaNormal = locGammaP4_OmegaCM.Vect();

    //Compute the theta angle to the omega decay plane
    double locCosTheta = locOmegaNormal.Dot(locHelicityZAxis_OmegaCM)/locOmegaNormal.Mag();
    double locTheta = acos(locCosTheta);

    //Compute the phi angle to the omega decay plane
    TVector3 locZCrossOmegaNormal = locHelicityZAxis_OmegaCM.Cross(locOmegaNormal);
    double locZCrossOmegaNormalMag = locZCrossOmegaNormal.Mag();
    double locCosPhi = locHelicityYAxis_OmegaCM.Dot(locZCrossOmegaNormal)/locZCrossOmegaNormalMag;
    double locPhi = acos(locCosPhi); //reports phi between 0 and pi: sign ambiguity
    //Resolve the sign ambiguity
    double locSinPhi = -1.0*locHelicityXAxis_OmegaCM.Dot(locZCrossOmegaNormal)/locZCrossOmegaNormalMag;
    if(locSinPhi < 0.0)
       locPhi *= -1.0;

    GDouble cosTheta = locCosTheta;
    GDouble sinSqTheta = sin(locTheta)*sin(locTheta);
    GDouble sin2Theta = sin(2.*locTheta);

    GDouble phi = locPhi;
    GDouble Phi = locPHI;

    //GDouble psi = phi - Phi;
    //if(psi < -1*PI) psi += 2*PI;
    //if(psi > PI) psi -= 2*PI;

    // vector meson production from K. Schilling et. al.

    GDouble Pgamma;
    if(polFraction >= 0.) Pgamma = polFraction;
    else{
       int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
       if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
          Pgamma = 0.;
       }
       else Pgamma = polFrac_vs_E->GetBinContent(bin);
    }
/*
    GDouble W = 1.0 - 0.5*(1. - rho000)*sinSqTheta - rho000*cosTheta*cosTheta + sqrt(2.)*rho100*sin2Theta*cos(phi) + rho1m10*sinSqTheta*cos(2.*phi);

    W -= Pgamma*cos(2.*Phi) * (2.*rho111 + (rho001-rho111)*sinSqTheta + sqrt(2.)*rho101*sin2Theta*cos(phi) + rho1m11*sinSqTheta*cos(2.*phi));

    W += Pgamma*sin(2.*Phi) * (sqrt(2.)*rho102*sin2Theta*sin(phi) + rho1m12*sinSqTheta*sin(2.*phi));

    */

   double rho110 = 0.5*(1.-rho000);

   GDouble W = 1.0 - sinSqTheta * rho110 - cosTheta*cosTheta*rho000 + sinSqTheta*cos(2.*phi)*rho1m10 + sqrt(2.)*rho100*sin2Theta*cos(phi);

   W -= Pgamma*cos(2.*Phi) * (2.*rho111 + sinSqTheta*(rho001-rho111) + sinSqTheta*cos(2.*phi)*rho1m11 + sqrt(2.)*rho101*sin2Theta*cos(phi));

   W += Pgamma*sin(2.*Phi) * (rho1m12*sinSqTheta*sin(2.*phi) + sqrt(2.)*rho102*sin2Theta*sin(phi));

    W *= 3./(8.*PI);

    return complex< GDouble > ( sqrt(fabs(W)) );
}

