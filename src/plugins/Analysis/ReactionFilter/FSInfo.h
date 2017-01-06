#ifndef FSINFO_H
#define FSINFO_H

using namespace std;

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <ctype.h>

#include "JANA/JStreamLog.h"
#include "particleType.h"

// ********************************************************************
// ********************************************************************
// ********************************************************************
//
//   The FSInfo class:  this class holds information about a
//                      given final state (FS)
//
// ********************************************************************
// ********************************************************************
// ********************************************************************

class FSInfo
{
	public:

		FSInfo(string FSName);

		// the name of the final state
		//  (this encodes a list of the final state particles
		//   and directions about whether to reconstruct them
		//   inclusively or exclusively)
		string FSName(void) { return m_FSName; }

		// compact string name of the reaction
		string ReactionName(void) { return m_ReactionName; }

		// is this an inclusive or exclusive final state?
		bool exclusive(void) const{ return (m_FSName.find("EXC") != string::npos); }
		bool inclusive(void) const{ return (m_FSName.find("INC") != string::npos); }
		bool missingN(void)  const{ return (m_FSName.find("MIN") != string::npos); }

		// use intermediate mass fits (true, unless NIMF is in the FS name)
		bool intermediateMassFits(void) const{ return !(m_FSName.find("NIMF") != string::npos); }

		// a list of particles associated with this final state
		vector<Particle_t>& PIDs(void) { return m_PIDs; }

		// total electric charge of all final state particles
		int totalCharge(void) const{ return m_totalCharge; };

	private:

		// private member data
		string m_FSName;
		vector<Particle_t> m_PIDs;
		int m_totalCharge;
		string m_ReactionName;

		// get short name (for building reaction/tree names)
		string Get_ShortName(Particle_t locPID) const;

		// functions to unpack the final state name and parse strings
		vector<Particle_t> getPIDsFromFSName( const string& FSName );
		vector<Particle_t> getPIDsFromFSName( const string& FSName, string& locReactionName);
		int getTotalChargeFromPIDs( const vector<Particle_t>& locPIDs ) const;
};

#endif

inline int FSInfo::getTotalChargeFromPIDs( const vector<Particle_t>& locPIDs ) const
{
	int totalCharge = 0;
	for(auto locPID : locPIDs)
		totalCharge += ParticleCharge(locPID);
	return totalCharge;
}

inline vector<Particle_t> FSInfo::getPIDsFromFSName( const string& FSName)
{
	string locReactionName;
	return getPIDsFromFSName(FSName, locReactionName);
}

inline string FSInfo::Get_ShortName(Particle_t locPID) const
{
	switch (locPID) {
	case Gamma:
		return "g";
	case Positron:
		return "ep";
	case Electron:
		return "em";
	case MuonPlus:
		return "mup";
	case MuonMinus:
		return "mum";
	case Pi0:
		return "pi0";
	case PiPlus:
		return "pip";
	case PiMinus:
		return "pim";
	case KPlus:
		return "kp";
	case KMinus:
		return "km";
	case Proton:
		return ""; //is understood
	case AntiProton:
		return "pbar";
	case KShort:
		return "ks";
	case Eta:
		return "eta";
	case Lambda:
		return "lambda";
	case AntiLambda:
		return "lbar";
	default:
		return "";
	}
}

