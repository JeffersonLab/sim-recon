// $Id$
//
//    File: DEventProcessor_run_summary.cc
// Created: Tue Nov 18 15:44:17 EST 2014
// Creator: sdobbs (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "DEventProcessor_run_summary.h"
#include "DEPICSstore.h"

#include <DAQ/DEPICSvalue.h>

#include <TDirectory.h>

static int VERBOSE = 0;

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_run_summary()); //register this plugin
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_run_summary::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... create historgrams or trees ...
	// japp->RootUnLock();
	//

	current_run_number = -1;
	conditions_tree = NULL;
	epics_info = NULL;
	
	// deal with command line parameters
	gPARMS->SetDefaultParameter("SUMMARY:VERBOSE", VERBOSE,
				    "Verbosity level for creating summary values.  0=no messages, 10=all messages");
	if(VERBOSE)
		cout << "Run summary verbosity level = " << VERBOSE << endl;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_run_summary::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes
	// keep track of current run number for use in our end run function
	current_run_number = locRunNumber;

	// make tree if it doesn't exist
	japp->RootWriteLock();

	//Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this plugin
	//If another thread has already created the folder, it just changes to it. 
	TDirectory *main_dir = gDirectory;
	gDirectory->cd("/");
	
	TDirectory* plugin_dir = static_cast<TDirectory*>(gDirectory->GetDirectory("conditions"));   // TDirectoryFile
	if(plugin_dir == NULL)
		plugin_dir = gDirectory->mkdir("conditions");
	plugin_dir->cd();
	
	if(gDirectory->Get("conditions") == NULL) //check to see if already created by another thread
	        conditions_tree = new TTree("conditions","");
	else //already created by another thread
		conditions_tree = static_cast<TTree*>(gDirectory->Get("conditions"));

	main_dir->cd();

	japp->RootUnLock();

	// reset EPICS summary info each run
	if(epics_info)
		delete epics_info;
	epics_info = new DEPICSstore;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_run_summary::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// read in whatever epics values are in this event
	vector<const DEPICSvalue*> epicsvalues;
	locEventLoop->Get(epicsvalues);	

	// save their values
	for(vector<const DEPICSvalue*>::const_iterator val_itr = epicsvalues.begin();
	    val_itr != epicsvalues.end(); val_itr++) {
		const DEPICSvalue* epics_val = *val_itr;
		if(VERBOSE)
		cout << "EPICS:  " << epics_val->name << " = " << epics_val->sval << endl;
		
		if(epics_info) 
			epics_info->AddValue(epics_val);
		else 
			jerr << "EPICS store object not loaded!" << endl;
	}
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_run_summary::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	
	// if we couldn't create the tree earlier, then throw it away
	if(conditions_tree == NULL)
		return NOERROR;

	// Although we are only filling objects local to this plugin, TTree::Fill() periodically writes to file: Global ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK

	// make a branch for the run number
	TBranch *run_branch = conditions_tree->FindBranch("run_number");
	if(run_branch == NULL)
		conditions_tree->Branch("run_number", &current_run_number, "run_number/D");
	else
		conditions_tree->SetBranchAddress("run_number", &current_run_number);
	
	// make branches for each of the values
	const map<string, DEPICSvalue_data_t> &epics_store = epics_info->GetStore();
	for(map<string, DEPICSvalue_data_t>::const_iterator epics_val_itr = epics_store.begin();
	    epics_val_itr != epics_store.end(); epics_val_itr++) {
		string branch_name = epics_val_itr->first;
		// ROOT branches assume that colons define the different leaves in a branch
		// so replace them in the name with underscores
		std::replace( branch_name.begin(), branch_name.end(), ':', '_');
		// also clean numerical expressions
		std::replace( branch_name.begin(), branch_name.end(), '-', '_');  
		std::replace( branch_name.begin(), branch_name.end(), '+', '_');  
		std::replace( branch_name.begin(), branch_name.end(), '*', '_');  
		std::replace( branch_name.begin(), branch_name.end(), '/', '_');  
		TBranch *the_branch = conditions_tree->FindBranch(branch_name.c_str());
		if(the_branch == NULL) {
			string branch_def = branch_name + "/D";
			conditions_tree->Branch(branch_name.c_str(), &(epics_val_itr->second.value->fval), branch_def.c_str());
		} else {
			conditions_tree->SetBranchAddress(branch_name.c_str(), &(epics_val_itr->second.value->fval));
		}
	}

	// save the values for this run
	conditions_tree->Fill();

	japp->RootUnLock(); //RELEASE ROOT LOCK

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_run_summary::fini(void)
{
	// Called before program exit after event processing is finished.
	//cout << "=================================================" << endl;
	//cout << "Summary of processed runs:" << endl;
	//cout << "=================================================" << endl;

	return NOERROR;
}

