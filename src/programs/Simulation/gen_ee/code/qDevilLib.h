#ifndef _QDEVILLIB_
#define _QDEVILLIB_

using std::string;
using std::cout;
using std::endl;

class PhasePT
{
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector q1;
  TLorentzVector q2;
  TLorentzVector recoil;
  TLorentzVector q12;
  TLorentzVector q23;
  double m12Cut;
  double rCut;
  double event_weight;
 public:
  TLorentzVector GetBeam(){return beam;}
  TLorentzVector GetTarget(){return target;}
  TLorentzVector GetQ1(){return q1;}
  TLorentzVector GetQ2(){return q2;}
  TLorentzVector GetQ12(){return q12;}
  TLorentzVector GetQ23(){return q23;}
  TLorentzVector GetRecoil(){return recoil;}
  double GetWeight(){return event_weight;};
  void SetBeam(TLorentzVector beamIn){beam = beamIn;}
  void SetTarget(TLorentzVector targetIn){target = targetIn;}
  void SetM12Cut(double m12CutIn){m12Cut = m12CutIn;}
  void SetRCut(double rCutIn){rCut = rCutIn;}
  int Gen(TRandom3 *random);
};

class Complx
{
  double real,
    imag;
 public:
  Complx( double r = 0., double i = 0.)  { real = r; imag = i; } // constructor
  Complx operator+(const Complx&) const;       // operator+()
  Complx operator-(const Complx&) const;       // operator-()
  Complx operator*(const Complx&) const;       // operator*()
  Complx operator=(const Complx&);       // operator=()
  void Show();
  void Set(double rVal,double iVal);
  double r(){return real;}
  double i(){return imag;}
  Complx Star(){
    Complx starVal(real,-imag);
    return starVal;
  }
  double Abs(){return sqrt(pow(real,2)+pow(imag,2));}
};
/*
// define constructor
Complx::Complx( double r, double i )
{
  real = r; imag = i;
}
*/
/////////////////////////////
class QedElement
{
  Complx scalar,scalar0,scalarX,scalarY,scalarZ,scalar5;
  Complx matrix[4][4],
    matrix0[4][4],
    matrixX[4][4],
    matrixY[4][4],
    matrixZ[4][4],
    matrix5[4][4],
    gamma0[4][4],
    gammaX[4][4],
    gammaY[4][4],
    gammaZ[4][4],
    gamma5[4][4];
  string QedType;
  string lIndexName;
  int lIndexPosition;
 public:
  QedElement()  { // constructor
  QedType = "none";
  lIndexName = "none";
  lIndexPosition = 0;
  scalar.Set(0,0);
  scalar0.Set(0,0);
  scalarX.Set(0,0);
  scalarY.Set(0,0);
  scalarZ.Set(0,0);
  scalar5.Set(0,0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix[i][j].Set(0,0); 
      matrix0[i][j].Set(0,0); 
      matrixX[i][j].Set(0,0); 
      matrixY[i][j].Set(0,0); 
      matrixZ[i][j].Set(0,0); 
      matrix5[i][j].Set(0,0); 
      gamma0[i][j].Set(0,0); 
      gammaX[i][j].Set(0,0); 
      gammaY[i][j].Set(0,0); 
      gammaZ[i][j].Set(0,0); 
      gamma5[i][j].Set(0,0); 
    }
  }
  //Gamma0                                                                      
  gamma0[0][0].Set(1,0);
  gamma0[1][1].Set(1,0);
  gamma0[2][2].Set(-1,0);
  gamma0[3][3].Set(-1,0);
  //GammaX                                                                      
  gammaX[0][3].Set(1,0);
  gammaX[1][2].Set(1,0);
  gammaX[2][1].Set(-1,0);
  gammaX[3][0].Set(-1,0);
  //GammaY                                                                      
  gammaY[0][3].Set(0,-1);
  gammaY[1][2].Set(0,1);
  gammaY[2][1].Set(0,1);
  gammaY[3][0].Set(0,-1);
  //GammaZ                                                                      
  gammaZ[0][2].Set(1,0);
  gammaZ[1][3].Set(-1,0);
  gammaZ[2][0].Set(-1,0);
  gammaZ[3][1].Set(1,0);
  //Gamma5                                                                      
  gamma5[0][2].Set(1,0);
  gamma5[1][3].Set(1,0);
  gamma5[2][0].Set(1,0);
  gamma5[3][1].Set(1,0);
  }
  QedElement operator=(const QedElement&);             // operator=()
  QedElement operator*(const QedElement&) const;       // operator*()
  QedElement operator+(const QedElement&) const;       // operator+()
  QedElement operator-(const QedElement&) const;       // operator-()
  void SetGamma(string lIndexNameVal);
  void SetZeroNone();
  void Transpose();
  void Star();
  void SetScalar(Complx sVal);
  void SetU(double px,double py,
	    double pz,double mass,int spin);
  void SetV(double px,double py,
	    double pz,double mass,int spin);
  void SetUBar(double px,double py,
	       double pz,double mass,int spin);
  void SetVBar(double px,double py,
	       double pz,double mass,int spin);
  void SetMomentumSlash(double px,double py,
		 double pz,double mass);
  void SetEpsilonSlash(double cx,double cy,
		       double cz);
  Complx GetScalar(){return scalar;};
  void ShowAll();
  void ShowMatrix(int matNumber);
};
/*
// define constructor
QedElement::QedElement()
{
  QedType = "none";
  lIndexName = "none";
  lIndexPosition = 0;
  scalar.Set(0,0);
  scalar0.Set(0,0);
  scalarX.Set(0,0);
  scalarY.Set(0,0);
  scalarZ.Set(0,0);
  scalar5.Set(0,0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix[i][j].Set(0,0); 
      matrix0[i][j].Set(0,0); 
      matrixX[i][j].Set(0,0); 
      matrixY[i][j].Set(0,0); 
      matrixZ[i][j].Set(0,0); 
      matrix5[i][j].Set(0,0); 
      gamma0[i][j].Set(0,0); 
      gammaX[i][j].Set(0,0); 
      gammaY[i][j].Set(0,0); 
      gammaZ[i][j].Set(0,0); 
      gamma5[i][j].Set(0,0); 
    }
  }
  //Gamma0                                                                      
  gamma0[0][0].Set(1,0);
  gamma0[1][1].Set(1,0);
  gamma0[2][2].Set(-1,0);
  gamma0[3][3].Set(-1,0);
  //GammaX                                                                      
  gammaX[0][3].Set(1,0);
  gammaX[1][2].Set(1,0);
  gammaX[2][1].Set(-1,0);
  gammaX[3][0].Set(-1,0);
  //GammaY                                                                      
  gammaY[0][3].Set(0,-1);
  gammaY[1][2].Set(0,1);
  gammaY[2][1].Set(0,1);
  gammaY[3][0].Set(0,-1);
  //GammaZ                                                                      
  gammaZ[0][2].Set(1,0);
  gammaZ[1][3].Set(-1,0);
  gammaZ[2][0].Set(-1,0);
  gammaZ[3][1].Set(1,0);
  //Gamma5                                                                      
  gamma5[0][2].Set(1,0);
  gamma5[1][3].Set(1,0);
  gamma5[2][0].Set(1,0);
  gamma5[3][1].Set(1,0);
}
*/

#endif // _QDEVILLIB_
