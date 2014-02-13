#if !defined(BARRIERFACTOR)
#define BARRIERFACTOR

// mass0 = mass of parent
// spin  = angular momentum of the decay
// mass1 = mass of first daughter
// mass2 = mass of second daughter

double barrierFactor( double mass0, int spin, double mass1, double mass2 );

// q     = breakup momentum
// spin  = angular momentum of the decay

double barrierFactor( double q, int spin );

#endif
