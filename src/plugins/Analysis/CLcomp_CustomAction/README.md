# CL comparison via CustomAction template
This plugin is a template for performing multiple kinematic fits and  comparing the confidence levels through a DCustomAction.
The motivation for this was the large background of rho events in the KK spectrum.

The documentation for this plugin is contained within this README. It provides a basic overview of the various things that need to be done in order to create your own plugin, however, the code will be required to see the finer details.

By using a CustomAction, you will need to run over REST data which can take more time. Cuts on the confidence level comparison will need to be performed within the CustomAction since these are not capable of writing to the output TTree. If you would rather decide upon the CL cut later in the DSelector, use the CustomBranch method instead (separate example plugin, CLcomp_CustomBranch).

# Setting up the reaction plugin with a CustomAction
The following describes the process to create your own plugin with a CustomAction to compare
the confidence levels of your choice of hypotheses. The examples will be specific to KK and pi pi.

## Making a custom plugin
- Run the command *MakeReactionPlugin.pl* without any arguments to see what's required.
- For the example, *MakeReactionPlugin.pl CLcomp_CustomAction p2k*
- cd into this new plugin directory

## Making a custom action
- Run the command *MakeAnalysisAction.pl* without any arguments to see what's required.
- For the example, *MakeAnalysisAction.pl CLcomp 1*

# Setting up the DReaction and output files

## Enable TTree output
- In order to make use of DSelectors, TTree output needs to be enabled
- Edit DEventProcessor.cc and find the *evnt* method
- Scan down to the commented region labled *REQUIRED*
- Uncomment the first chunk of code that involves DEventWriterROOT and Fill_DataTrees

## DReaction.cc
- This is the file that contains the reaction information.
- Add *#include "DCustomAction_CLcomp.h"* to the top of this file
- Follow the instructions in the commented out code to add your reaction.
- Make sure that the kinfit is uncommented and set to the desired type.
- Uncomment the line that allows TTree output.
- In the Analysis Actions section, add/enable whatever basic cuts and histograms you would like for your reaction.
- Apply the HistogramAction_KinFitResults and CutAction_KinFitFOM(locReaction, 0.0)
- Apply the CustomAction at this point.
- It's useful to histogram the invariant mass after each action that cuts combos.

## DCustomAction.h
- Add the following headers
  - #include "KINFITTER/KinFitter.h"
  - #include "ANALYSIS/DKinFitUtils_Gluex.h"
  - #include "PID/DKinematicData.h"
  - Any ROOT related headers (such as TH1.h or TH2.h)
- In the *private* section
  - Declare dAnalysisUtils, dKinFitUtils_Gluex, and dKinFitter
  - Declare any histograms, etc.

## DCustomAction.cc
- In *Initialize*
  - Define histograms in *Initialize*
  - Set up:
    - dAnalysisUtils
    - dKinFitUtils
    - dKinFitter
- In *Perform_Action*
  - Verify that the kinfit was enabled
  - Get the CL from the reaction
  - Reset_NewEvent()
  - Get the original set of particle hypotheses
  - Get the FinalParticle_SourceObjects() which are in the order of the reaction
  - Loop over the source objects and assign the appropriate hypothesis to the base track.
  - Reset_NewFit()
  - Construct P4 and/or vertex constraints
  - Set the KinFit debug level and perform the fit
  - Get the new hypotheses' CL
  - Fill any desired histograms
  - Return ( original_CL / new_CL > 1 ) (true keeps the combo, false cuts it)
