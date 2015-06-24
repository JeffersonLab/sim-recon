#if !defined(NBODYPHASESPACEFACTORY)
#define NBODYPHASESPACEFACTORY

#include <vector>

#include "TLorentzVector.h"

using namespace std;

class NBodyPhaseSpaceFactory
{
	
 public:
	
  NBodyPhaseSpaceFactory( double parentMass, const vector<double>& childMass );


  /**
   * Generates N-body phase space decays using the Raubold Lynch
   * method (F.James CERN 68-15 (1968).
   * Borrows liberally from ROOT class TGenPhaseSpace as well as AcquRoot
   * class TMCGenerator (J.R.M.Annand -- http://nuclear.gla.ac.uk/~acqusys/doc)
   *
   * \param[in] uniformWeights - boolean value selecting whether to perform
   *  accept/reject on generated events, resulting in uniform event weights,
   *  or to return the first generated event and set its corresponding weight.
  */	
  vector<TLorentzVector> generateDecay(bool uniformWeights = true );

  /**
   * Returns the weight of the last-generated event. This is relevant
   * in case the last event was generated with uniformWeights=false
   */
  double getLastGeneratedWeight() const {return m_lastWt;};

 private:
        
  static const double kPi;
	
  double pdk( double a, double b, double c ) const;
  double random( double low, double hi ) const;
	
  double m_parentMass;
  vector<double> m_childMass;    // vector of daughter masses
  int m_Nd;                      // number of decay products
  double m_lastWt;

};

#endif
