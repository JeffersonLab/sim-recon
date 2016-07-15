
#include "mcsmear_config.h"
#include "RCDB/ConfigParser.h"

#include <iostream>
#include <fstream>

//-----------
// mcsmear_config_t (Constructor)
//-----------
mcsmear_config_t::mcsmear_config_t() 
{
	// default values
	DROP_TRUTH_HITS = false;
	ADD_NOISE      = false;
	SMEAR_HITS     = true;
	//SMEAR_BCAL     = true;
	IGNORE_SEEDS   = false;
    DUMP_RCDB_CONFIG = false;
	APPLY_EFFICIENCY_CORRECTIONS = true;
		
	TRIGGER_LOOKBACK_TIME = -100; // ns
		
#ifdef HAVE_RCDB
	// RCDB configuration
	// first determine which database to connect to
	if( getenv("RCDB_CONNECTION")!= NULL )
		RCDB_CONNECTION = getenv("RCDB_CONNECTION");
	else
		RCDB_CONNECTION = "mysql://rcdb@hallddb/rcdb";   // default to outward-facing MySQL DB

    // load RCDB later, so that the DApplication interface is initialized
    rcdb_connection = NULL;
#endif //HAVE_RCDB
}
	
//-----------
// mcsmear_config_t (Destructor)
//-----------
mcsmear_config_t::~mcsmear_config_t() {
	delete rcdb_connection;
}

//-----------
// SetSeeds
//-----------
void mcsmear_config_t::SetSeeds(const char *vals) 
{
  	/// This is called from the command line parser to
   	/// set the initial seeds based on user input from
   	/// the command line.
   	//
   	stringstream ss(vals);
   	Int_t seed1, seed2, seed3;
   	ss >> seed1 >> seed2 >> seed3;
   	UInt_t *useed1 = reinterpret_cast<UInt_t*>(&seed1);
   	UInt_t *useed2 = reinterpret_cast<UInt_t*>(&seed2);
   	UInt_t *useed3 = reinterpret_cast<UInt_t*>(&seed3);
   	gDRandom.SetSeeds(*useed1, *useed2, *useed3);

   	cout << "Seeds set from command line. Any random number" << endl;
   	cout << "seeds found in the input file will be ignored!" << endl;
   	IGNORE_SEEDS = true;
}


#ifdef HAVE_RCDB

//-----------
// LoadRCDBConnection
//-----------
void  mcsmear_config_t::LoadRCDBConnection() 
{
    // We want to connect to RCDB as late as possible for two reasons:
    //  1) No need to connect to the database unless we are actually using this information,
    //     which should speed things up
    //  2) To use the JANA command line parameter interface, we need to make sure that the
    //     DApplication is initialized, and it might not be when the mcsmear_config_t
    //     constructor is called

    // if we're already connected, then stop now
    if(rcdb_connection != NULL)
        return;

	gPARMS->SetDefaultParameter("RCDB_CONNECTION", RCDB_CONNECTION, "URL used to access RCDB.");
    
	// load connection to RCDB
	rcdb_connection = new rcdb::Connection(RCDB_CONNECTION);
}

//-----------
// ParseRCDBConfigFile
//-----------
bool mcsmear_config_t::ParseRCDBConfigFile(int runNumber)
{
	// This is just for testing (right now)
	// To get a lot of the configuration parameters we need, we need to parse the CODA configuration files
	// which are saved in RCDB in the JSON format

    // Lazily connect to RCDB
    LoadRCDBConnection();
    
    // The "rtvs" condition contains the file name of the CODA configuration file
    // What else does this contain??
    auto rtvsCondition = rcdb_connection->GetCondition(runNumber, "rtvs");    // Get condition by run and name
    if(!rtvsCondition) {
        jerr << "RCDB: 'rtvs' condition is not set for run " << runNumber << endl;
        return false;
    }

    auto json = rtvsCondition->ToJsonDocument();               // The CODA rtvs is serialized as JSon dictionary.
    string fileName(json["%(config)"].GetString());            // The file name is stored in '%(config)'


    // Get file out of RCDB (indexed by run number and name)
    auto file = rcdb_connection->GetFile(runNumber, fileName);
    if(!file) {                                                      // If there is no such file, null is returned
        jerr << "File with name \"" << fileName
             << "\" doesn't exist for run "<< runNumber << endl;
        return false;
    }

    // Parse CODA config file 
    vector<string> SectionNames = {"TRIGGER", "GLOBAL", "FCAL", "BCAL", "TOF", "ST", "TAGH",
                                         "TAGM", "PS", "PSC", "TPOL", "CDC", "FDC"};
    string fileContent = file->GetContent();                               // Get file content
    auto result = rcdb::ConfigParser::Parse(fileContent, SectionNames);    // Parse it!

    // EXAMPLE
    double CDC_FADC125_DAC = stod(result.Sections["CDC"].NameValues["FADC125_DAC"]);
    double CDC_FADC125_THR = stod(result.Sections["CDC"].NameValues["FADC125_THR"]);

    if(DUMP_RCDB_CONFIG) {
        //// DEBUG ////
        //cout << "CODA config file contents:" << endl;
        ofstream coda_ofile("rcdb_coda.config");
        coda_ofile << "Full file: " << endl;
        coda_ofile << fileContent << endl;


        coda_ofile << endl << "Parsed:" << endl;
        for(auto sectionName: SectionNames) {
            coda_ofile << "Section: " << sectionName << endl;
            auto sectionData = result.Sections[sectionName];
            
            for(auto value : sectionData.NameValues)
                coda_ofile << value.first << " = " << value.second << endl;
            for(auto data_vec : sectionData.NameVectors) {
                coda_ofile << data_vec.first << " = ";
                for(auto value: data_vec.second) 
                    coda_ofile << value << " ";
                coda_ofile << endl;
            } 

            coda_ofile << endl;
        }
        coda_ofile.close();
    }

	// then we get stuff out of it

	return true;
}

#endif // HAVE_RCDB

