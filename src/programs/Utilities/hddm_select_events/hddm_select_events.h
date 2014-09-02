// $Id$
//
// Created Oct 09, 2013  Kei Moriya

/*************************************************************
 *
 * 2013/10/09 Kei Moriya
 *
 * This program is based on hddm_cull_events and will
 * select events for output based on a user-specified
 * criteria.
 *
 * The purpose is for example events containing Lambda
 * events, where the particle must go through GEANT
 * to have decay products, but some selection criteria
 * on the daughters is wanted.
 *
 * I have rewritten some of the original code
 * so that the options are easier to use (in my mind)
 * but this deviates from the behavior of other
 * hddm_xxx programs.
 *
 * This program will not work on REST format.
 * To do this, a new program will have to be written
 * based on the hddm_r file.
 *
 * Usage:
 * hddm_select_events [-i Inputfile] [-o Outputfile] [-r input is REST] \
 *                    [-a save remainder events] [-s selection type] \
 *                    [-M maximum number of events] [-d debug]
 *
 * selection types:
 * 1. select Lambda -> p pi-
 * 2. select Lambda -> p pi-, eta -> gamma gamma
 * 3. select Lambda -> p pi-, decay length less than 2 cm
 *
 *************************************************************/

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include "unistd.h" // to use optarg

using namespace std;

#include <signal.h>
#include <time.h>
#include <stdlib.h>

#include <HDDM/hddm_s.hpp>
#include <HDDM/hddm_r.hpp>

#include "TRandom2.h"
#include "TLorentzVector.h"

extern bool HDDM_USE_COMPRESSION;
extern bool HDDM_USE_INTEGRITY_CHECKS;

bool selectEvent_s(int select_type, hddm_s::HDDM &record, int nevents, bool debug);
bool selectEvent_r(int select_type, hddm_r::HDDM &record, int nevents, bool debug);

// Lambda decay constant
const double alpha = 0.642;
