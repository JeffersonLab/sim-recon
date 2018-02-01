// $Id$
//
//    File: DEventWriterROOT_CLcomp.h
// Created: Thu Feb  1 13:48:42 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#ifndef _DEventWriterROOT_CLcomp_
#define _DEventWriterROOT_CLcomp_

#include <map>
#include <string>

#include <ANALYSIS/DEventWriterROOT.h>
#include <ANALYSIS/DAnalysisUtilities.h>
#include <KINFITTER/DKinFitter.h>
#include <ANALYSIS/DKinFitUtils_GlueX.h>

using namespace std;
using namespace jana;

class DEventWriterROOT_CLcomp : public DEventWriterROOT
{
	public:
		virtual ~DEventWriterROOT_CLcomp(void){};

	protected:

		//CUSTOM FUNCTIONS: //Inherit from this class and write custom code in these functions

		virtual void Create_CustomBranches_ThrownTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop) const;
		virtual void Fill_CustomBranches_ThrownTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const;

		virtual void Create_CustomBranches_DataTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction, bool locIsMCDataFlag) const;
		virtual void Fill_CustomBranches_DataTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DReaction* locReaction, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
				const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
				const vector<const DNeutralParticleHypothesis*>& locNeutralHypos, const deque<const DParticleCombo*>& locParticleCombos) const;
};

#endif //_DEventWriterROOT_CLcomp_
