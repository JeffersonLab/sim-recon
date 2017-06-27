// $Id$
//
//    File: JEventProcessor_pi0calib.cc
// Created: Wed May 24 13:46:12 EDT 2017
// Creator: mashephe (on Linux stanley.physics.indiana.edu 2.6.32-642.6.2.el6.x86_64 unknown)
//

#include "JEventProcessor_pi0calib.h"
#include "DFactoryGenerator_pi0calib.h"

using namespace jana;

#include "evio_writer/DEventWriterEVIO.h"

#include <ANALYSIS/DAnalysisResults.h>
#include <ANALYSIS/DEventWriterROOT.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_pi0calib());
    app->AddFactoryGenerator( new DFactoryGenerator_pi0calib() );
  }
} // "C"


//------------------
// JEventProcessor_pi0calib (Constructor)
//------------------
JEventProcessor_pi0calib::JEventProcessor_pi0calib()
{

  WRITE_EVIO_FILE = 1;
  gPARMS->SetDefaultParameter( "WRITE_EVIO_FILE", WRITE_EVIO_FILE );

  WRITE_ROOT_TREE = 0;
  gPARMS->SetDefaultParameter( "WRITE_ROOT_TREE", WRITE_ROOT_TREE );
}

//------------------
// ~JEventProcessor_pi0calib (Destructor)
//------------------
JEventProcessor_pi0calib::~JEventProcessor_pi0calib()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_pi0calib::init(void)
{
  // This is called once at program startup. 

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_pi0calib::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_pi0calib::evnt(JEventLoop *loop, uint64_t eventnumber)
{

  vector<const DAnalysisResults*> analysisResultsVector;
  loop->Get( analysisResultsVector );
  
  if( WRITE_EVIO_FILE ){
    
    const DEventWriterEVIO* eventWriterEVIO = NULL;
    loop->GetSingle(eventWriterEVIO);

    // write out BOR events
    if(loop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT)) {

      eventWriterEVIO->Write_EVIOEvent(loop, "exclusivepi0");
      return NOERROR;
    }

    if( analysisResultsVector[0]->Get_NumPassedParticleCombos() != 0 ){

      eventWriterEVIO->Write_EVIOEvent(loop, "exclusivepi0");
    }
  }
  
  if( WRITE_ROOT_TREE ){
    
    //Recommended: Write surviving particle combinations (if any) to output ROOT TTree
    //If no cuts are performed by the analysis actions added to a DReaction, then this saves all of its particle combinations.
    //The event writer gets the DAnalysisResults objects from JANA, performing the analysis.
    // string is DReaction factory tag: will fill trees for all DReactions that are defined in the specified factory

    const DEventWriterROOT* eventWriterROOT = NULL;
    loop->GetSingle(eventWriterROOT);
    eventWriterROOT->Fill_DataTrees(loop, "pi0calib");
  }
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_pi0calib::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_pi0calib::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

