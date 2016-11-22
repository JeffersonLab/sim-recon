#include "FSInfo.h"

// **********************************************
//
//   FSInfo constructor
//
// **********************************************
FSInfo::FSInfo(string FSName) : m_FSName(FSName)
{
	// print out final state name
	// jout << "FSInfo:  Initializing Final State " << FSName << endl;

	// set up member data
	m_PIDs         = getPIDsFromFSName(FSName, m_ReactionName);
	m_totalCharge  = getTotalChargeFromPIDs(m_PIDs);

	// print out particle list
	for (size_t i = 0; i < m_PIDs.size(); i++)
		jout << "FSInfo:      " << ParticleType(m_PIDs[i]) << endl;

	// do some checks
	if (missingN() && totalCharge() != 0 && totalCharge() != 1)
	{
		jerr << "FSInfo ERROR 5: wrong charge for final state name: " << FSName << endl;
		exit(1);
	}
	if (exclusive() && totalCharge() != 1)
	{
		jerr << "FSInfo ERROR 6: wrong charge for final state name: " << FSName << endl;
		exit(1);
	}
}


// **********************************************
//
//   FSInfo: MODE NUMBERING UTILITIES, ETC.
//
// **********************************************
vector<Particle_t> FSInfo::getPIDsFromFSName( const string& FSName, string& locReactionName)
{
	// some quick checks
	if ((FSName.size() == 0) || (FSName.find("_") == string::npos))
	{
		jerr << "FSInfo ERROR 1: error in final state name: " << FSName << endl;
		exit(1);
	}

	// a list of allowed particle names
	map<int, Particle_t> locIndexToPIDMap;
	locIndexToPIDMap[0] = Pi0;
	locIndexToPIDMap[1] = PiMinus;
	locIndexToPIDMap[2] = PiPlus;
	locIndexToPIDMap[3] = KShort;
	locIndexToPIDMap[4] = KMinus;
	locIndexToPIDMap[5] = KPlus;
	locIndexToPIDMap[6] = Gamma;
	locIndexToPIDMap[7] = Eta;
	locIndexToPIDMap[8] = AntiProton;
	locIndexToPIDMap[9] = Proton;
	locIndexToPIDMap[10] = MuonMinus;
	locIndexToPIDMap[11] = MuonPlus;
	locIndexToPIDMap[12] = Electron;
	locIndexToPIDMap[13] = Positron;
	locIndexToPIDMap[14] = AntiLambda;
	locIndexToPIDMap[15] = Lambda;

	// parse FSName digit by digit, starting at the end
	int index = 0;
	bool lcode2 = false;
	vector<Particle_t> locPIDs;
	ostringstream locReactionStream;
	for (int i = FSName.size()-1; i >= 0 && index < 16; i--)
	{
		string digit = FSName.substr(i,1);
		if(isdigit(digit[0]))
		{
			Particle_t locPID = locIndexToPIDMap[index];
			int num = atoi(digit.c_str());
			for (int j = 0; j < num; j++)
				locPIDs.push_back(locPID);

			if(num > 1)
				locReactionStream << num;
			locReactionStream << Get_ShortName(locPID);

			index++;
		}
		else if (digit == "_" && !lcode2)
		{
			if (index > 7)
			{
				jerr << "FSInfo ERROR 2: error in final state name: " << FSName << endl;
				exit(1);
			}
			lcode2 = true;
			index = 7;
		}
		else
		{
			if (!lcode2)
			{
				jerr << "FSInfo ERROR 3: error in final state name: " << FSName << endl;
				exit(1);
			}
			break;
		}
	}

	// make sure we have particles
	if (locPIDs.empty())
	{
		jerr << "FSInfo ERROR 4: error in final state name: " << FSName << endl;
		exit(1);
	}

	locReactionName = locReactionStream.str();
	return locPIDs;
}


