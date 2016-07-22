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

#ifndef _daq_param_type_
#define _daq_param_type_

enum daq_param_type{

	kPARAM250_NSA             = 0x0501,
	kPARAM250_NSB             = 0x0502,
	kPARAM250_NSA_NSB         = 0x0503,  // NSA+NSB
	kPARAM250_NPED            = 0x0504,

	kPARAM125_NSA             = 0x0F01,
	kPARAM125_NSB             = 0x0F02,
	kPARAM125_NSA_NSB         = 0x0F03,  // NSA+NSB
	kPARAM125_NPED            = 0x0F04,  // (1 + 2^P2)/(2^PBIT)
	kPARAM125_WINWIDTH        = 0x0F05,
	kPARAM125_PL              = 0x0F06,
	kPARAM125_NW              = 0x0F07,
	kPARAM125_NPK             = 0x0F08,
	kPARAM125_P1              = 0x0F09,
	kPARAM125_P2              = 0x0F0A,
	kPARAM125_PG              = 0x0F0B,
	kPARAM125_IE              = 0x0F0C,
	kPARAM125_H               = 0x0F0D,
	kPARAM125_TH              = 0x0F0E,
	kPARAM125_TL              = 0x0F0F,
	kPARAM125_IBIT            = 0x0F10,
	kPARAM125_ABIT            = 0x0F11,
	kPARAM125_PBIT            = 0x0F12,

	kPARAMF1_REFCNT           = 0x0601,
	kPARAMF1_TRIGWIN          = 0x0602,
	kPARAMF1_TRIGLAT          = 0x0603,
	kPARAMF1_HSDIV            = 0x0604,
	kPARAMF1_BINSIZE          = 0x0605,  // in picoseconds
	kPARAMF1_REFCLKDIV        = 0x0606,
	
	kPARAMCAEN1290_WINWIDTH   = 0x1001,
	kPARAMCAEN1290_WINOFFSET  = 0x1002,

	kPARAM_NONE               = 0x0000
};

#endif // _daq_param_type_