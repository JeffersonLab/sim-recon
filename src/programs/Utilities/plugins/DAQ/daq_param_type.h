// $Id$
// $HeadURL$
//
//    File: daq_param_type.h
// Created: Fri Sep 26 08:06:38 EDT 2014
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0)
//
//
// These specify the daq configuration parameter codes used
// to identify the the parameter types in the DAQ config.
// parameter EVIO bank. These are written into the event by
// the DAQ system and parsed in the method:
//
// JEventSource_EVIO::ParseModuleConfiguration
//
// The rawevent plugin also uses these when encoding this 
// bank. See more details here:
//
//  https://halldweb1.jlab.org/wiki/images/2/21/20140910_config_in_datastream.pdf
//

enum daq_param_type{

	kPARAM250_NSA             = 0x0501,
	kPARAM250_NSB             = 0x0502,
	kPARAM250_NSA_NSB         = 0x0503,  // NSA+NSB
	kPARAM250_NPED            = 0x0504,

	kPARAM125_NSA             = 0x0F01,
	kPARAM125_NSB             = 0x0F02,
	kPARAM125_NSA_NSB         = 0x0F03,  // NSA+NSB
	kPARAM125_NPED            = 0x0F04,
	kPARAM125_WINWIDTH        = 0x0F05,

	kPARAMF1_REFCNT           = 0x0601,
	kPARAMF1_TRIGWIN          = 0x0602,
	kPARAMF1_TRIGLAT          = 0x0603,
	kPARAMF1_HSDIV            = 0x0604,
	kPARAMF1_BINSIZE          = 0x0605,  // in picoseconds
	
	kPARAMCAEN1290_WINWIDTH   = 0x1001,
	kPARAMCAEN1290_WINOFFSET  = 0x1002,

	kPARAM_NONE               = 0x0000
};
