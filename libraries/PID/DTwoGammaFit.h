//
//    File: DTwoGammaFit.h
// Created: Tue Avg 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_
#define _DTwoGammaFit_

///#include <TLorentzVector.h>
#include "DKinematicData.h"
#include "DPhoton.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTwoGammaFit:public DKinematicData {
	public:
		HDCLASSDEF(DTwoGammaFit);
                
                DTwoGammaFit();
		~DTwoGammaFit();

               double getChi2() const;
               double getProb() const ;
               inline double getPull(const int i) const { return fPulls[i]; }
;
               DKinematicData* getChildFit(const int i) ;

               void setChi2(double const aChi2);  
               void setProb(double const aProb);  
               void setPulls(double const aPull, const int i);  
               void setChildFit(const DKinematicData& aChildFit, const int i);  

	private:

               double fProb; // 
               double fChi2; // 
               double fPulls[9];
               DKinematicData fChildFit[2];

};


// Getters
// return confidence of the fit
inline double DTwoGammaFit::getProb() const
{
      return fProb;
}

// return chi2 of the fit
inline double DTwoGammaFit::getChi2() const
{
      return fChi2;
}

inline DKinematicData* DTwoGammaFit::getChildFit(const int i) 
{
      return &fChildFit[i];
}

/* return pull of fitted child
inline double DTwoGammaFit::getPull(const int i) const
{
      return fPulls[i];
}*/


// Setters
// Set data of fitted children
inline void DTwoGammaFit::setChildFit(const DKinematicData& aChildFit, const int i)
{
     fChildFit[i] = aChildFit;
}

// Set pulls from DKinFit
inline void DTwoGammaFit::setPulls(const double aPull, const int i)
{
     fPulls[i] = aPull;
}

// Set chi2 from DKinFit
inline void DTwoGammaFit::setChi2(const double aChi2)
{
     fChi2 = aChi2;
}

// Set confidence from DKinFit
inline void DTwoGammaFit::setProb(const double aProb)
{
     fProb = aProb;
}

/*
// return origin of pi0  (Fcal=0, Bcal=1)
inline unsigned int DPi0::getChildrenTag(int child) const
{
      return fTags[child];
}

// return child ID's
inline oid_t DPi0::getChildrenID(int child) const
{
      return fIDs[child];
}
*/

#endif // _DTwoGammaFit_

