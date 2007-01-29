// $Id: DEventProcessor_invariant_mass_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_invariant_mass_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_invariant_mass_hists_
#define _DEventProcessor_invariant_mass_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

class DTrack;
class DMCThrown;

class DEventProcessor_invariant_mass_hists:public JEventProcessor{

	public:
		DEventProcessor_invariant_mass_hists();
		~DEventProcessor_invariant_mass_hists();
		
		TH1F *mass_1part, *mass_2part, *mass_3part;
		TH1F *mass2_1part, *mass2_2part, *mass2_3part;

		TH1F *thrown_mass_1part, *thrown_mass_2part, *thrown_mass_3part;
		TH1F *thrown_mass2_1part, *thrown_mass2_2part, *thrown_mass2_3part;
		
		TH1F *missing_mass, *thrown_missing_mass;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void MakeTLorentz(const DTrack *track, TLorentzVector &v);
		void MakeTLorentz(const DMCThrown *thrown, TLorentzVector &v);
		double SampleGaussian(double sigma);

};

#endif // _DEventProcessor_invariant_mass_hists_

