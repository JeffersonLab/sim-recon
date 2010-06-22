/*
 *  DKinematicData.cc
 *  Hall D
 *
 *  Created by Matthew Shepherd on 5/28/07.
 *
 */

#include "PID/DKinematicData.h"

#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "DRandom.h"

const double kLarge = 1.e20;

// B field constant to convert curvature to momentum
// const double kBfieldConstant = -0.0299792458;

//
// static data member definitions
//
DMatrixDSym* DKinematicData::nullMatrix()
{
  //static DMatrixDSym* sNullMatrix = new DMatrixDSym(0,7);
  static DMatrixDSym* sNullMatrix = new DMatrixDSym(7);
  return sNullMatrix;
}                         
DMatrixDSym* DKinematicData::null5x5Matrix(){
  static DMatrixDSym* sNullMatrix = new DMatrixDSym(5);
  return sNullMatrix;
}

//
// constructors and destructor
//
DKinematicData::DKinematicData() :
  m_hasFixedMass( !false ) ,
  m_mass( kDefaultMass ) ,
  m_charge( kDefaultCharge ) ,
  m_momentum( DVector3( 0.0 , 0.0 , 0.0 ) ) ,
  m_position( DVector3(  0.0 , 0.0 , 0.0 ) ) ,
  m_errorMatrix( nullMatrix() ),
  m_TrackingErrorMatrix( null5x5Matrix() )
{
  m_t0 = -999.0;
  m_t0_err = -999.0;
  m_t0_detector = SYS_NULL; /// No detector
  m_t1 = -999.0;
  m_t1_err = -999.0;
  m_t1_detector = SYS_NULL; /// No detector
  m_pathLength = -999.0;
  m_pathLength_err = -999.0;
  m_dedx=0.;
  m_use_forward_parameters=false;
  return ;
}

// constructor with Jana id
DKinematicData::DKinematicData(const oid_t id) :
  JObject ( id ) ,
  m_hasFixedMass( !false ) ,
  m_mass( kDefaultMass ) ,
  m_charge( kDefaultCharge ) ,
  m_momentum( DVector3( 0.0 , 0.0 , 0.0 ) ) ,
  m_position( DVector3(  0.0 , 0.0 , 0.0 ) ) ,
  m_errorMatrix( nullMatrix() ),
  m_TrackingErrorMatrix( null5x5Matrix() )
{
  return ;
}


// Copy constructor with optional argument for copying over error matrix
DKinematicData::DKinematicData( const DKinematicData& aKinematicData) :
  m_hasFixedMass( aKinematicData.m_hasFixedMass ) ,
  m_mass( aKinematicData.m_mass ) ,
  m_charge( aKinematicData.m_charge ) ,
  m_momentum( aKinematicData.m_momentum ) ,
  m_position( aKinematicData.m_position ) ,
  m_errorMatrix( nullMatrix() ),
  m_TrackingErrorMatrix( null5x5Matrix() )
{
  // Copy error matrix if it exists
  if(! aKinematicData.hasNullErrorMatrix() ){
    setErrorMatrix( *aKinematicData.m_errorMatrix ) ;
  }
  if (! aKinematicData.hasNull5x5Matrix() ) {
    setTrackingErrorMatrix(*aKinematicData.m_TrackingErrorMatrix);
  }

  m_t0 = aKinematicData.m_t0;
  m_t0_err = aKinematicData.m_t0_err;
  m_t0_detector = aKinematicData.m_t0_detector;
  m_t1 =aKinematicData.m_t1;
  m_t1_err = aKinematicData.m_t1_err;
  m_t1_detector = aKinematicData.m_t1_detector;
  m_pathLength = aKinematicData.m_pathLength;
  m_pathLength_err = aKinematicData.m_pathLength_err;
  m_dedx=aKinematicData.m_dedx;
  m_use_forward_parameters=aKinematicData.m_use_forward_parameters;
  return ;
}


// Similar to copy constructor with optional argument for copying
// over error matrix
DKinematicData::DKinematicData( const DKinematicData& aKinematicData,
    const bool aCopyErrorMatrix) :
  m_hasFixedMass( aKinematicData.m_hasFixedMass ) ,
  m_mass( aKinematicData.m_mass ) ,
  m_charge( aKinematicData.m_charge ) ,
  m_momentum( aKinematicData.m_momentum ) ,
  m_position( aKinematicData.m_position ) ,
  m_errorMatrix( nullMatrix() ),
  m_TrackingErrorMatrix(null5x5Matrix() )
{
  // Copy error matrix if is requested and if it exists
  if(aCopyErrorMatrix && ! aKinematicData.hasNullErrorMatrix() ){
    setErrorMatrix( *aKinematicData.m_errorMatrix ) ;
  }
  if(aCopyErrorMatrix && ! aKinematicData.hasNull5x5Matrix() ){
    setTrackingErrorMatrix( *aKinematicData.m_TrackingErrorMatrix ) ;
  }

  m_t0 = aKinematicData.m_t0;
  m_t0_err = aKinematicData.m_t0_err;
  m_t0_detector = aKinematicData.m_t0_detector;
  m_t1 =aKinematicData.m_t1;
  m_t1_err = aKinematicData.m_t1_err;
  m_t1_detector = aKinematicData.m_t1_detector;
  m_pathLength = aKinematicData.m_pathLength;
  m_pathLength_err = aKinematicData.m_pathLength_err;
  m_dedx=aKinematicData.m_dedx;
  m_use_forward_parameters=aKinematicData.m_use_forward_parameters;
  return ;
}

DKinematicData::DKinematicData(const DVector3& aMomentum ,
    const DVector3&  aPosition ,
    const ValueType aMass ,
    const ValueType aCharge):
  m_hasFixedMass( !false ) ,
  m_mass( aMass ) ,
  m_charge( aCharge ) ,
  m_momentum( aMomentum ) ,
  m_position( aPosition ) ,
  m_errorMatrix( nullMatrix() ),
  m_TrackingErrorMatrix(null5x5Matrix() )

{
  return ;
}

DKinematicData::DKinematicData(const DVector3& aMomentum ,
    const DVector3&  aPosition ,
    const ValueType aMass ,
    const ValueType aCharge,
    const DMatrixDSym& aErrorMatrix):
  m_hasFixedMass( !false ) ,
  m_mass( aMass ) ,
  m_charge( aCharge ) ,
  m_momentum( aMomentum ) ,
  m_position( aPosition ) ,
  m_errorMatrix( nullMatrix()),
  m_TrackingErrorMatrix(null5x5Matrix())
{
  setErrorMatrix( aErrorMatrix);
  return ;
}

DKinematicData::~DKinematicData()
{
  clearErrorMatrix();
  clearTrackingErrorMatrix();
  return ;
}

//
// assignment operators
//
  const DKinematicData& 
DKinematicData::operator=( const DKinematicData& aOtherKinematicData )
{
  if(this != &aOtherKinematicData) {
    m_hasFixedMass = aOtherKinematicData.m_hasFixedMass ;
    m_mass = aOtherKinematicData.m_mass ;
    m_charge = aOtherKinematicData.m_charge ;
    m_momentum = aOtherKinematicData.m_momentum ;
    m_position = aOtherKinematicData.m_position ;
    if ( ! aOtherKinematicData.hasNullErrorMatrix() ) {
      setErrorMatrix( *aOtherKinematicData.m_errorMatrix ) ;
    }
    else {
      clearErrorMatrix();
    } 
    if ( ! aOtherKinematicData.hasNull5x5Matrix() ) {
      setTrackingErrorMatrix( *aOtherKinematicData.m_TrackingErrorMatrix ) ;
    }
    else {
      clearTrackingErrorMatrix();
    }
    m_t0 = aOtherKinematicData.m_t0;
    m_t0_err = aOtherKinematicData.m_t0_err;
    m_t0_detector = aOtherKinematicData.m_t0_detector;
    m_t1 =aOtherKinematicData.m_t1;
    m_t1_err = aOtherKinematicData.m_t1_err;
    m_t1_detector = aOtherKinematicData.m_t1_detector;
    m_pathLength = aOtherKinematicData.m_pathLength;
    m_pathLength_err = aOtherKinematicData.m_pathLength_err;
    m_dedx=aOtherKinematicData.m_dedx;
    m_use_forward_parameters=aOtherKinematicData.m_use_forward_parameters;
  }
  return ( *this ) ;
}

bool
DKinematicData::operator==( const DKinematicData& rhs ) const
{
  if( !( m_hasFixedMass == rhs.m_hasFixedMass &&
	 m_mass == rhs.m_mass &&
	 m_charge == rhs.m_charge &&
	 m_momentum == rhs.m_momentum &&
	 m_position == rhs.m_position &&
	 hasNullErrorMatrix() == rhs.hasNullErrorMatrix() &&
	 hasNull5x5Matrix()==rhs.hasNull5x5Matrix()
	 ) )
  {
    return false ;
  }
  else if( hasNullErrorMatrix() )
  {
    return true ;
  }
  else
  {
    return *m_errorMatrix == *( rhs.m_errorMatrix ) ;
  }
}

bool
DKinematicData::operator!=( const DKinematicData& rhs ) const
{
  return !( *this == rhs ) ;
}

//
// member functions
//

void DKinematicData::setForwardParmFlag(bool aFlag){
  m_use_forward_parameters=aFlag;
}


  void 
DKinematicData::setMass( const ValueType aMass )
{
  m_mass = aMass ;
}

  void 
DKinematicData::setMomentum( const DVector3& aMomentum )
{
  m_momentum = aMomentum ;
}

  void 
DKinematicData::setPosition( const DVector3& aPosition )
{
  m_position = aPosition ;
}

  void 
DKinematicData::setCharge( const ValueType aCharge )
{
  m_charge = aCharge ;
}

void DKinematicData::setdEdx(const ValueType adedx){
  m_dedx=adedx;
}

  void 
DKinematicData::setMassFixed( void )
{
  m_hasFixedMass = !false ;
}

  void 
DKinematicData::setMassFloat( void )
{
  m_hasFixedMass = false ;
}

  void 
DKinematicData::clearErrorMatrix( void )
{
  if( !hasNullErrorMatrix() ) {
    delete m_errorMatrix ;
    m_errorMatrix = nullMatrix() ;
  }
}

void DKinematicData::clearTrackingErrorMatrix(void){
  if (!hasNull5x5Matrix()){
    delete m_TrackingErrorMatrix;
    m_TrackingErrorMatrix=null5x5Matrix();
  }
}

void DKinematicData::setTrackingErrorMatrix(const DMatrixDSym& aMatrix){
  const int kMatrixSize = 5;
  assert( kMatrixSize == aMatrix.GetNrows() );

  if (&aMatrix!= null5x5Matrix()){
    clearTrackingErrorMatrix();
    m_TrackingErrorMatrix = new DMatrixDSym(aMatrix);
  }
}


  void 
DKinematicData::setErrorMatrix( const DMatrixDSym& aMatrix )
{
  const int kMatrixSize = 7;
  assert( kMatrixSize == aMatrix.GetNrows() );

  //Check to see if new Matrix is null
  bool isNullErrorMatrix = true;
  if( &aMatrix != nullMatrix() ) {
    //for(int column = 1; column <= kMatrixSize; ++column) {
      //for(int row = column; row <= kMatrixSize; ++row) {
    for(int column = 0; column < kMatrixSize; ++column) {
      for(int row = column; row < kMatrixSize; ++row) {
        if(0 != aMatrix(row,column)) {
          isNullErrorMatrix = false;
          goto endOfLoop;
        }
      }
    }
  }

endOfLoop: clearErrorMatrix();
           if( !isNullErrorMatrix ) {
             m_errorMatrix = new DMatrixDSym(aMatrix);
           }
}


/// Setters for the timing information
void DKinematicData::setT0( const ValueType at0, const ValueType at0_err, const DetectorSystem_t at0_detector)
{
  m_t0 = at0;
  m_t0_err = at0_err;
  m_t0_detector = at0_detector;
}

void DKinematicData::setT1( const ValueType at1, const ValueType at1_err, const DetectorSystem_t at1_detector)
{
  m_t1 = at1;
  m_t1_err = at1_err;
  m_t1_detector = at1_detector;
}

void DKinematicData::setPathLength( const ValueType apathLength, const ValueType apathLength_err)
{
  m_pathLength = apathLength;
  m_pathLength_err = apathLength_err;
}

/////////////////////
  DMatrixDSym*
DKinematicData::takeOwnershipOfPointer( void )
{
  // Remove the error matrix pointer from control of the object.
  // The matrix contents are not touched
  if( hasNullErrorMatrix() ) {
    return 0;
  }
  else {
    DMatrixDSym* tempPointer = m_errorMatrix;
    m_errorMatrix = nullMatrix();
    return tempPointer;
  }
}

  void
DKinematicData::restoreOwnershipOfPointer( DMatrixDSym* const aPointer)
{
  // Restore the error matrix pointer to the object.
  // The matrix contents are not touched
  if(aPointer) {
    m_errorMatrix = aPointer;
  }
}

  void 
DKinematicData::smearMCThrownMomentum( double smearPct )
{
  DRandom rnd;
  double errMatValueX = smearPct*smearPct*m_momentum.X()*m_momentum.X();
  double errMatValueY = smearPct*smearPct*m_momentum.Y()*m_momentum.Y();
  double errMatValueZ = smearPct*smearPct*m_momentum.Z()*m_momentum.Z();

  rnd.SetSeed((unsigned int)(1000*m_momentum.X()));

  DMatrixDSym sigmas(7);
  sigmas(0,0) = errMatValueX;
  sigmas(1,1) = errMatValueY;
  sigmas(2,2) = errMatValueZ;

  setErrorMatrix( sigmas );

  double newPx = m_momentum.X() + rnd.Gaus( 0.0 , smearPct*m_momentum.X() );
  double newPy = m_momentum.Y() + rnd.Gaus( 0.0 , smearPct*m_momentum.Y() );
  double newPz = m_momentum.Z() + rnd.Gaus( 0.0 , smearPct*m_momentum.Z() );

  m_momentum = DVector3(newPx, newPy, newPz);
}

/******

Note:  The following two member functions were left here for reference.  They
where used to convert tracking parameters and covariance matrix to
KinematicData covariance matrix.  Depending on tracking implemenation
in GlueX it is likely that a similar function will be needed.

template <class T, unsigned int IRow, unsigned int IColumn>
class KTMatrixArray
{
public:
KTMatrixArray(bool iClearMemory = true ) { 
if( iClearMemory) {clearMemory();} }
//NOTE: indicies begin with 1
T& operator()(unsigned int iRow, unsigned int iColumn) {
return m_array[ (iRow-1) + IRow*(iColumn-1) ];
}
private:
void clearMemory() {
memset(m_array, 0, sizeof(m_array));
}
T m_array[IRow*IColumn];
};


void
DKinematicData::calculate7x7ErrorMatrixFrom5x5ErrorMatrix(
const KTHelix& aHelix)
{
// This code converted from the Fortran code in kwfit.

// Build the 7x7 covariance matrix, V_w, in the Kinematic representation
// using the 5x5 Helix covariance matrix, V_c.
//
// Let W be the 7 Kinematic parameters
// Let C be the 5 Helix (CLEO) parameters
// Let A be the matrix of derivatives relating the 7 Kinematic
//      parameters to the 5 Helix parameters. The elements of
//      A are given by A_ij = d(W_i) / d(C_j).
//
// The 7x7 covariance matrix is then given by
//
//     V_w = A * V_c * A(t)
//
// The equations relating the W parameters to the C parameters are shown below.
//
// W = (px, py, pz, E, x, y, z)
//
// C = (cu, phi0, d0, ct, z0)
// cu = curvature
// phi0 = phi angle at point of closest approach to reference point
// d0 = distance of closest approach to reference point
// ct = cot(theta)
// z0 = z position at point of closest approach to reference point
//
// BfieldConstant = -0.0299792458
//
// pt = BfieldConstant * Bfield * charge / cu
//
// px = pt * cos(phi0)
// py = pt * sin(phi0)
// pz = pt * ct
// E  = sqrt(px^2 + py^2 + pz^2 + mass^2)
//  x = xref - d0 * sin(phi0)
//  y = yref + d0 * cos(phi0)
//  z = zref + z0
//
// The derivatives are obtained from the differential expressions
//       dpx = -px/cu*dcu - py*dphi0             = A11*dcu + A12*dphi0
//       dpy = -py/cu*dcu + px*dphi0             = A21*dcu + A22*dphi0
//       dpz = -pz/cu*dcu + pt*dct               = A31*dcu + A34*dct
//       de  = -p^2/E/cu*dcu + ct*pt^2/E*dct     = A41*dcu + A44*dct
//       dx  = -py/pt*dd0 - y*dphi0              = A52*dphi0+ A53*dd0
//       dy  =  px/pt*dd0 + x*dphi0              = A62*dphi0+ A63*dd0
//       dz  =  dz0                              = A75*dz0 = dz0

if(aHelix.hasNullErrorMatrix()) {return;}

double cuinv = aHelix.curvature() == 0. ? kLarge : 1. / aHelix.curvature();

double ct = aHelix.cotTheta();
double px = momentum().x();
double py = momentum().y();
double pz = momentum().z();
double pt = sqrt(px*px + py*py);
double p = momentum().mag();
double Einv = 1. / energy();
double x = position().x();
double y = position().y();

// Create 7x5 transformation matrix filled with zeros
KTMatrixArray<double, 7,5> A;

// This makes the matrix more understandable
const unsigned int kCurvature = KTHelix::kCurvature;
const unsigned int kPhi0      = KTHelix::kPhi0;
const unsigned int kD0        = KTHelix::kD0;
const unsigned int kCotTheta  = KTHelix::kCotTheta;
const unsigned int kZ0        = KTHelix::kZ0;

A(kPx, kCurvature) = -px * cuinv;
A(kPy, kCurvature) = -py * cuinv;
A(kPz, kCurvature) = -pz * cuinv;
A(kEnergy, kCurvature) = -p*p * Einv * cuinv;
// A(kX, kCurvature) = 0.;
// A(kY, kCurvature) = 0.;
// A(kZ, kCurvature) = 0.;


A(kPx,kPhi0) = -py;
A(kPy,kPhi0) =  px;
// A(kPz,kPhi0) = 0.;
// A(kEnergy,kPhi0) = 0.;
A(kX,kPhi0) = -y;
A(kY,kPhi0) =  x;
// A(kZ,kPhi0) = 0.;

// A(kPx, kD0) = 0.;
// A(kPy, kD0) = 0.;
// A(kPz, kD0) = 0.;
// A(kEnergy, kD0) = 0.;
A(kX, kD0) = -py / pt;
A(kY, kD0) =  px / pt;
// A(kZ, kD0) = 0.;

// A(kPx, kCotTheta) = 0.;
// A(kPy, kCotTheta) = 0.;
A(kPz, kCotTheta) = pt;
A(kEnergy, kCotTheta) = ct * pt*pt * Einv;
// A(kX, kCotTheta) = 0.;
// A(kY, kCotTheta) = 0.;
// A(kZ, kCotTheta) = 0.;

// A(kPx, kZ0) = 0.;
// A(kPy, kZ0) = 0.;
// A(kPz, kZ0) = 0.;
// A(kEnergy, kZ0) = 0.;
// A(kX, kZ0) = 0.;
// A(kY, kZ0) = 0.;
A(kZ, kZ0) = 1.;

// Calculate V_w = A * V_c * A(t).
//   setErrorMatrix ( aHelix.errorMatrix().similarity( A ) );

// The following code is fast because it only multiplies the non-zero
// elements of A.
const DMatrixDSym& Vtkc = aHelix.errorMatrix();
KTMatrixArray<double, 7,5> AVtkc(false);
unsigned int i;
for (i=1; i < 6; ++i) {
  double Vtkc_1_i = Vtkc.fast(i,1);
  double Vtkc_2_i;
  if( i == 1) {
    Vtkc_2_i = Vtkc.fast(2,1);
  } else {
    Vtkc_2_i = Vtkc.fast(i,2);
  }
  double Vtkc_3_i = Vtkc(3,i);
  double Vtkc_4_i = Vtkc(4,i);
  AVtkc(1,i) = Vtkc_1_i*A(1,1) + Vtkc_2_i*A(1,2);
  AVtkc(2,i) = Vtkc_1_i*A(2,1) + Vtkc_2_i*A(2,2);
  AVtkc(3,i) = Vtkc_1_i*A(3,1) + Vtkc_4_i*A(3,4);
  AVtkc(4,i) = Vtkc_1_i*A(4,1) + Vtkc_4_i*A(4,4);
  AVtkc(5,i) = Vtkc_2_i*A(5,2) + Vtkc_3_i*A(5,3);
  AVtkc(6,i) = Vtkc_2_i*A(6,2) + Vtkc_3_i*A(6,3);
  AVtkc(7,i) = Vtkc.fast(5,i);
}

if(m_errorMatrix == nullMatrix()){
  m_errorMatrix = new DMatrixDSym(7,0);
}
DMatrixDSym& Vtkw = *m_errorMatrix;

//when  using the 'fast' method the row index must be >= column index
const unsigned int kMatrixSize = 7;
for (i=1; i<kMatrixSize+1; ++i) {
  Vtkw.fast(i,1) = A(1,1)*AVtkc(i,1) + A(1,2)*AVtkc(i,2);
}
for (i=2; i<kMatrixSize+1; ++i) {
  Vtkw.fast(i,2) = A(2,1)*AVtkc(i,1) + A(2,2)*AVtkc(i,2);
}
for (i=3; i<kMatrixSize+1; ++i) {
  Vtkw.fast(i,3) = A(3,1)*AVtkc(i,1) + A(3,4)*AVtkc(i,4);
}
for (i=4; i<kMatrixSize+1; ++i) {
  Vtkw.fast(i,4) = A(4,1)*AVtkc(i,1) + A(4,4)*AVtkc(i,4);
}
Vtkw.fast(5,5) = A(5,2)*AVtkc(5,2) + A(5,3)*AVtkc(5,3);
Vtkw.fast(6,5) = A(5,2)*AVtkc(6,2) + A(5,3)*AVtkc(6,3);
Vtkw.fast(7,5) = A(5,2)*AVtkc(7,2) + A(5,3)*AVtkc(7,3);
Vtkw.fast(6,6) = A(6,2)*AVtkc(6,2) + A(6,3)*AVtkc(6,3);
Vtkw.fast(7,6) = A(6,2)*AVtkc(7,2) + A(6,3)*AVtkc(7,3);
Vtkw.fast(7,7) = AVtkc(7,5);

}

********/

//
// const member functions
//

DKinematicData::ValueType 
DKinematicData::mass( void ) const
{
  return ( m_mass ) ;
}

DKinematicData::ValueType 
DKinematicData::charge( void ) const
{
  return ( m_charge ) ;
}

DKinematicData::ValueType 
DKinematicData::px( void ) const
{
  return ( momentum().x() ) ;
}

DKinematicData::ValueType 
DKinematicData::py( void ) const
{
  return ( momentum().y() ) ;
}

DKinematicData::ValueType 
DKinematicData::pz( void ) const
{
  return ( momentum().z() ) ;
}

DKinematicData::ValueType 
DKinematicData::energy( void ) const
{
  return ( sqrt( ( mass() * mass() ) +
        ( momentum().Mag2() ) ) ) ;
}

DKinematicData::ValueType
DKinematicData::x( void ) const
{
  return ( position().x() ) ;
}

DKinematicData::ValueType
DKinematicData::y( void ) const
{
  return ( position().y() ) ;
}

DKinematicData::ValueType
DKinematicData::z( void ) const
{
  return ( position().z() ) ;
}

DKinematicData::ValueType 
DKinematicData::pperp( void ) const
{
  return ( sqrt( momentum().x()*momentum().x()
        + momentum().y()*momentum().y() ) );
}

DKinematicData::ValueType 
DKinematicData::pperp2( void ) const
{
  return ( momentum().x()*momentum().x()
      + momentum().y()*momentum().y() );
}

DKinematicData::ValueType 
DKinematicData::pmag( void ) const
{
  return ( momentum().Mag() );
}

DKinematicData::ValueType 
DKinematicData::pmag2( void ) const
{
  return ( momentum().Mag2() );
}

const DVector3& 
DKinematicData::momentum( void ) const
{
  return ( m_momentum ) ;
}

const DVector3& 
DKinematicData::position( void ) const
{
  return ( m_position ) ;
}

const DLorentzVector 
DKinematicData::lorentzMomentum( void ) const
{
  return ( DLorentzVector( momentum() , energy() ) ) ;
}

bool DKinematicData::forwardParmFlag(void) const{
  return (m_use_forward_parameters);
}


bool 
DKinematicData::hasFixedMass( void ) const
{
  return ( m_hasFixedMass ) ;
}

const DMatrixDSym& 
DKinematicData::errorMatrix( void ) const
{
  return ( *m_errorMatrix ) ;
}

const DMatrixDSym& DKinematicData::TrackingErrorMatrix(void) const
{
  return ( *m_TrackingErrorMatrix ) ;
}

/// Timing info
DKinematicData::ValueType DKinematicData::t0( void ) const { return ( m_t0 ) ; }
DKinematicData::ValueType DKinematicData::t0_err( void ) const { return ( m_t0_err ) ; }
DetectorSystem_t DKinematicData::t0_detector( void ) const { return ( m_t0_detector ) ; }
DKinematicData::ValueType DKinematicData::t1( void ) const { return ( m_t1 ) ; }
DKinematicData::ValueType DKinematicData::t1_err( void ) const { return ( m_t1_err ) ; }
DetectorSystem_t DKinematicData::t1_detector( void ) const { return ( m_t1_detector ) ; }
DKinematicData::ValueType DKinematicData::pathLength( void ) const { return ( m_pathLength ) ; }
DKinematicData::ValueType DKinematicData::pathLength_err( void ) const { return ( m_pathLength_err ) ; }

DKinematicData::ValueType DKinematicData::TOF( void ) const { return ( (m_t1 - m_t0) ) ; }
DKinematicData::ValueType DKinematicData::TOF_err( void ) const { return sqrt(m_t1_err*m_t1_err + m_t0_err*m_t0_err);}


DKinematicData::ValueType DKinematicData::dEdx( void ) const
{ return ( m_dedx); } 


DKinematicData::ValueType DKinematicData::measuredBeta( void ) const { return ( (m_pathLength/(m_t1 - m_t0) ))/SPEED_OF_LIGHT ; }
DKinematicData::ValueType DKinematicData::measuredBeta_err( void ) const
{
  double a = pow(m_pathLength_err/m_pathLength, 2);
  double b = (pow(m_t1_err,2) + pow(m_t0_err,2))/pow(m_t1 - m_t0,2);

  double err = measuredBeta()*sqrt(a + b);

  return err;
}

DKinematicData::ValueType DKinematicData::deltaBeta( void ) const { return ( lorentzMomentum().Beta() - measuredBeta() ) ; }
DKinematicData::ValueType DKinematicData::deltaInvBeta( void ) const { return ( 1.0/lorentzMomentum().Beta() - 1.0/measuredBeta() ) ; }
DKinematicData::ValueType DKinematicData::measuredInvBeta_err( void ) const
{
  double a = (m_t1_err*m_t1_err + m_t0_err*m_t0_err)/pow(m_t1 - m_t0,2);
  double b = pow(m_pathLength_err/m_pathLength,2);

  double err = (1.0/measuredBeta())*sqrt(a+b);

  return err;
}


//
// static member functions
//

