//
//    File: DTwoGammaFit.h
// Created: Tue Avg 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_
#define _DTwoGammaFit_

#include "DKinematicData.h"

#include <JANA/JObject.h>
using namespace jana;

class DTwoGammaFit:public DKinematicData {
	public:
		JOBJECT_PUBLIC(DTwoGammaFit);
                
                DTwoGammaFit();
                DTwoGammaFit(const JObject::oid_t id);
		~DTwoGammaFit();

// Getters: 
               inline double getChi2() const { return fChi2; }

               inline double getNdf() const { return fNdf; }

               inline double getProb() const { return fProb; }

// Puls of px1, ... ,pz2
               inline double getPull(const int i) const { return fPulls[i]; } 

// mass before the fit
               inline double getUMass() const { return fUMass; }  

// look-up children
               inline oid_t getChildID(int child)  const { return fIDs[child]; }

// photon after the fit
               const DKinematicData* getChildFit(const int i) const ;

// Get momenta before the fit
               const DLorentzVector* getChildMom(const int i) const ;

// Setters:
               void setChi2(double const aChi2);  
	       void setNdf(int const aNdf);  
               void setProb(double const aProb);  
               void setPulls(double const aPull, const int i);  
               void setUMass(double const uMass);
               void setChildID(const JObject::oid_t aID, const int i ); 
               void setChildFit(const DKinematicData& aChildFit, const int i);
               void setChildMom(const DLorentzVector& aChildMom, const int i);


			void toStrings(vector<pair<string,string> > &items)const{
				DKinematicData::toStrings(items);
				AddString(items, "prob" , "%4.3lf", fProb);
				AddString(items, "chisq", "%4.3lf", fChi2);
				AddString(items, "Ndof" , "%d"   , fNdf);
				AddString(items, "mass" , "%4.3lf", fMass);
				AddString(items, "umass", "%4.3lf", fUMass);
			}
	private:

               JObject::oid_t fIDs[2];
               double fProb;  
               double fChi2;  
	       int fNdf;
               double fMass;  
               double fUMass;  // unconstrained mass
               double fPulls[6]; 
               DKinematicData fChildFits[2];
               DLorentzVector fChildMoms[2];
};

// Getters

inline const DKinematicData* DTwoGammaFit::getChildFit(const int i) const
{
      return &fChildFits[i];
}


inline const DLorentzVector* DTwoGammaFit::getChildMom(const int i) const
{
      return &fChildMoms[i];
}

// Setters
// Set data of fitted children
inline void DTwoGammaFit::setChildFit(const DKinematicData& aChildFit, const int i)
{
     fChildFits[i] = aChildFit;
}
// Set data of fitted children
inline void DTwoGammaFit::setChildMom(const DLorentzVector& aChildMom, const int i)
{
     fChildMoms[i] = aChildMom;
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

// Set Ndf from DKinFit
inline void DTwoGammaFit::setNdf(const int aNdf)
{
     fNdf = aNdf;
}

// Set confidence from DKinFit
inline void DTwoGammaFit::setProb(const double aProb)
{
     fProb = aProb;
}

// Set unconstrained mass 
inline void DTwoGammaFit::setUMass(const double uMass)
{
     fUMass = uMass;
}

inline void DTwoGammaFit::setChildID(const oid_t aID, const int i )
{
   fIDs[i] = aID;
}


#endif // _DTwoGammaFit_
