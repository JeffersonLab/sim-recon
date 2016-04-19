#!/usr/bin/env perl

$ActionName = "NULL";
$ActionType = 0; # 0 for reaction-independent, else reaction-dependent

if($#ARGV < 1)
{
	&Usage();
	exit;
}
else
{
	$ActionName = $ARGV[0];
	$ActionType = $ARGV[1];
}

print "\n";
print "Generating files for DAnalysisAction \"DCustomAction_$ActionName\" \n";

# Create C++ Header File for Factory Class
$dhfile = $fname = "DCustomAction_${ActionName}.h";
open(FILE, ">$dhfile");
&PrintFileHeader();
&PrintActionClass();
close(FILE);
print " - $dhfile\n";

# Create C++ Implementation file for Factory Class
$ccfile = $fname = "DCustomAction_${ActionName}.cc";
open(FILE, ">$ccfile");
&PrintFileHeader();
&PrintActionMethods();
close(FILE);
print " - $ccfile\n";

###############
# PrintFileHeader
###############
sub PrintFileHeader()
{
	# print a few lines at the very top of the file
	$uname = `uname -nprs`;
	chomp($uname);
	print FILE "// \$Id\$\n";
	print FILE "//\n";
	print FILE "//    File: $fname\n";
	print FILE "// Created: ".`date`;
	print FILE "// Creator: ".$ENV{"USER"}." (on $uname)\n";
	print FILE "//\n";
}

###############
# Usage
###############
sub Usage()
{
	print "\n";
	print "Usage:\n\t MakeAnalysisAction action_name action_type\n";
	print "\n";
	print "Generate the C++ source and header files to implement a custom DAnalysisAction\n";
	print "for use with the JANA framework.\n";
	print "\n";
	print "\"action_type\" should be 0 (zero) to create a reaction-independent action.\n";
	print "\"action_type\" should be 1 to create a reaction-dependent action.\n";
	print "\n";
	print "Add these files to your plugin, and then add the analysis action to your DReaction.\n";
	print "Be sure to \#include the header file in your DReaction factory.\n";
	print "\n";
}

###############
# PrintActionClass
###############
sub PrintActionClass()
{
	$content = "
\#ifndef _DCustomAction_${ActionName}_
\#define _DCustomAction_${ActionName}_

\#include <string>
\#include <iostream>

\#include \"JANA/JEventLoop.h\"
\#include \"JANA/JApplication.h\"

\#include \"ANALYSIS/DAnalysisAction.h\"
\#include \"ANALYSIS/DReaction.h\"
\#include \"ANALYSIS/DParticleCombo.h\"
\#include \"ANALYSIS/DAnalysisUtilities.h\"

using namespace std;
using namespace jana;

class DCustomAction_${ActionName} : public DAnalysisAction
{
	public:
";

	if(${ActionType} == 0) # reaction-independent
	{
		$content .=	"
		//user can call any of these three constructors
		DCustomAction_${ActionName}(const DReaction* locReaction, string locActionUniqueString = \"\") : 
		DAnalysisAction(locReaction, \"Custom_${ActionName}\", false, locActionUniqueString) {}

		DCustomAction_${ActionName}(string locActionUniqueString) : 
		DAnalysisAction(NULL, \"Custom_${ActionName}\", false, locActionUniqueString) {}

		DCustomAction_${ActionName}(void) : 
		DAnalysisAction(NULL, \"Custom_${ActionName}\", false, \"\") {}";
	}
	else # reaction-dependent
	{
		$content .=	"
		DCustomAction_${ActionName}(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = \"\") : 
		DAnalysisAction(locReaction, \"Custom_${ActionName}\", locUseKinFitResultsFlag, locActionUniqueString) {}";
	}
	$content .=	"

		void Initialize(JEventLoop* locEventLoop);

	private:
";
	if(${ActionType} == 0) # reaction-independent
	{
		$content .=	"
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);";
	}
	else # reaction-dependent
	{
		$content .=	"
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);";
	}
	$content .=	"

		//Store any histograms as member variables here
};

\#endif // _DCustomAction_${ActionName}_

";
	print FILE $content;
}

###############
# PrintActionMethods
###############
sub PrintActionMethods()
{
	$content = "
\#include \"DCustomAction_${ActionName}.h\"

void DCustomAction_${ActionName}::Initialize(JEventLoop* locEventLoop)
{
	/*
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock
";
	if(${ActionType} == 0) # reaction-independent
	{
		$content .=	"
	//When creating a reaction-independent action, only modify member variables within a ROOT lock. 
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously. 
";
	}

		$content .=	"
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		// Optional: Create a ROOT subfolder.
			//If another thread has already created the folder, it just changes to it. 
		// CreateAndChangeTo_Directory(\"MyDirName\", \"MyDirTitle\");
			//make sub-directory content here
		// gDirectory->cd(\"..\"); //return to the action directory

		//	(Optional) Example: Create a histogram. 
			// This function will return the histogram if already created by another thread. If not pre-existing, it will create and return it. 
			// Function arguments are identical to those used for the histogram constructors
		// dMyHist = GetOrCreate_Histogram<TH1I>(\"MyHistName\", \"MyHistTitle\", 100, 0.0, 1.0);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
	*/
}

bool DCustomAction_${ActionName}::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

";
	if(${ActionType} == 0) # reaction-independent
	{
		$content .=	"
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.
";
	}

	$content .=	"

	// Optional: Useful utility functions.
	// const DAnalysisUtilities* locAnalysisUtilities;
	// locEventLoop->GetSingle(locAnalysisUtilities);

	//Optional: check whether the user wanted to use the kinematic fit results when performing this action
	// bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();

";
	if(${ActionType} == 0) # reaction-independent
	{
		$content .=	"
	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	// if(Get_NumPreviousParticleCombos() != 0)
	//		return true;
";
	}

	$content .=	"
	/*
	//Optional: Fill histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill any histograms here
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
	*/

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
";
	print FILE $content;
}

