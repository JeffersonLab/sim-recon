//
//    File: DTwoGammaFit.h
// Created: Tue Avg 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_
#define _DTwoGammaFit_

#include "DKinematicData.h"
#include "DPhoton.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTwoGammaFit:public DKinematicData {
	public:
		HDCLASSDEF(DTwoGammaFit);
                
                DTwoGammaFit();
		~DTwoGammaFit();

               inline double getChi2() const { return fChi2; }
               inline double getProb() const { return fProb; }
               inline double getPull(const int i) const { return fPulls[i]; }
// the origin of the  photon (FCAL, BCAL, charged)
               inline int getChildTag(int child) const { return fTags[child]; }
               inline oid_t getChildID(int child)  const { return fIDs[child]; }

               DKinematicData* getChildFit(const int i) ;

               void setChi2(double const aChi2);  
               void setProb(double const aProb);  
               void setPulls(double const aPull, const int i);  
               void setChildTag(const int aTag, const int i ); 
               void setChildID(const oid_t aID, const int i ); 
               void setChildFit(const DKinematicData& aChildFit, const int i);  

	private:

               oid_t fIDs[2];  
               int fTags[2]; // tag children origin (FCAL/BCAL/charged)
               double fProb;  
               double fChi2;  
               double fMass;  
               double fPulls[6]; // needs specification
               DKinematicData fChildFits[2];

};

// Getters

inline DKinematicData* DTwoGammaFit::getChildFit(const int i) 
{
      return &fChildFits[i];
}

// Setters
// Set data of fitted children
inline void DTwoGammaFit::setChildFit(const DKinematicData& aChildFit, const int i)
{
     fChildFits[i] = aChildFit;
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

// set Pi0 bits with respect to the photon detection 
inline void DTwoGammaFit::setChildTag(const int aTag, const int i )
{
   fTags[i] = aTag;
}

inline void DTwoGammaFit::setChildID(const oid_t aID, const int i )
{
   fIDs[i] = aID;
}


#endif // _DTwoGammaFit_

