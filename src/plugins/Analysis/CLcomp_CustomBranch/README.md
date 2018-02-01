# CL comparison via CustomBranch template
This is proof-of-concept code and is not yet ready for production.

This plugin is a template for performing multiple kinematic fits and  comparing the confidence levels through a DCustomBranch.
The motivation for this was the large background of rho events in the KK spectrum.

The documentation for this plugin is contained within this README. It provides a basic overview of the various things that need to be done in order to create your own plugin, however, the code will be required to see the finer details.

By using a CustomBranch, you will need to run over REST data once to save information to a reduced TTree. You can then make any cuts on the CL through your DSelector. If you would rather cut upon the CL immediately in the DReaction, use the CustomAction method instead (separate example plugin, CLcomp_CustomAction).

# Setting up the reaction plugin with a CustomBranch
The following describes the process to create your own plugin with a CustomBranch to compare
the confidence levels of your choice of hypotheses. The examples will be specific to KK and pi pi.

## Making a custom plugin
- Run the command *MakeReactionPlugin.pl* without any arguments to see what's required.
- For the example, *MakeReactionPlugin.pl CLcomp_CustomBranch p2k*
- cd into this new plugin directory

## Making a custom EventWriter
- Run the command *MakeEventWriterROOT.pl* without any arguments to see what's required.
- For the example, *MakeEventWriterROOT.pl CLcomp*
- The script prints out what to do to implement the DEventWriterROOT. These changes will be detailed in the next section.

# Setting up the plugin

## DFactoryGenerator.h
- Include the following header
  - #include "DEventWriterROOT_factory_CLcomp.h"
- Add the following line to GenerateFactories()
  - locEventLoop->AddFactory(new DEventWriterROOT_factory_CLcomp());

## DEventProcessor.h
- Include the following header
  - #include "DEventWriterROOT_CLcomp.h"

## DEventProcessor.cc
- In the *evnt* method
  - Uncomment the DEventWriterROOT section
  - Replace the 2 EventWriterROOT lines with your following
    - const DEventWriterROOT_CLcomp* locEventWriterROOT = NULL;
    - locEventLoop->GetSingle(locEventWriterROOT, "CLcomp");

## DReaction.cc
- This is the file that contains the reaction information.
- Follow the instructions in the commented out code to add your reaction.
- Make sure that the kinfit is uncommented and set to the desired type.
- Uncomment the line that allows TTree output.
- In the Analysis Actions section, add/enable whatever basic cuts and histograms you would like for your reaction.

## DEventWriterROOT_CLcomp.h
- Add the following headers
  - #include "KINFITTER/KinFitter.h"
  - #include "ANALYSIS/DKinFitUtils_Gluex.h"
  - #include "PID/DKinematicData.h"

## DEventWriterROOT_CLcomp.cc
- In *Create_CustomBranches_DataTree()*, add the following lines to create branches.
  - unsigned int locInitArraySize = 10;
  - locBranchRegister.Register_Single<UInt_t>( "ComboArraySize" );
  - locBranchRegister.Register_FundamentalArray<double>("kpkm_CL__double", "ComboArraySize", locInitArraySize );
  - locBranchRegister.Register_FundamentalArray<double>("p2pi_CL__double", "ComboArraySize", locInitArraySize );
  - locBranchRegister.Register_FundamentalArray<double>("IM__double", "ComboArraySize", locInitArraySize );
- In *Fill_CustomBranches_DataTree()*
  - This is where the bulk of the code is located. See the actual code for details.
  - Create the following objects:
    - DAnalysisUtilities
    - DKinFitUtils_Gluex
    - DKinFitter
  - For uniqueness tracking (this will likely be changed to the DSelector), create a map<pair<oid_t, oid_t>, double> for:
    - KK pair and associated CL
    - PiPi pair and associated CL
    - KK pair and associated invariant mass
  - Loop over the particle combos
  - Get the confidence level from the DReaction
  - dKinFitter->Reset_NewEvent()
  - Get the initial and final state particles
    - Use the combostep index for unchanging particles
    - Get the FinalState_SourceObjects to get the DChargedTrack information and change hypotheses
  - dKinFitter->Reset_NewFit()
  - Set the P4 constraints
  - Set the vertex constraints
  - Run the kinfit
  - Check for uniqueness
  - Fill branches
