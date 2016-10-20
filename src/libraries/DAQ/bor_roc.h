//
// Data structures for holding module configuration information
// read back from the modules themselves for packaging into the
// Beginning Of Run (BOR) record.
//
// IMPORTANT: THIS FILE IS MAINTAINED IN 2 REPOSITORIES!  
//
// Before making changes to this file, keep these points in mind:
//
// 1. This file is maintained in separate repositories since the
//    code that generates the BOR during data taking resides
//    in the online repository and the code that must interpret
//    it resides in the offline repository. The locations of the
//    file in the two repositories are:
//
//     https://halldsvn.jlab.org/repos/trunk/online/daq/vme/src/rcm/monitor/bor_roc.h
//     https://github.com/JeffersonLab/sim-recon/tree/master/src/libraries/DAQ/bor_roc.h
//
//    Any changes MUST be made in both places!
//
// 2. The BOR data structures are stored directly in the raw data 
//    files. Any changes made will not be available in existing
//    data files. Any modifications here should be added to the 
//    end of the data structures so there is some (small) hope
//    of being able to handle both old and new raw data files.
//
// 3. Any items added should either be uint32_t or an array of 
//    uint16_t with and even number of elements. This is so
//    the data will stay close packed on 64bit platforms where
//    the compiler will tend to add padding bytes so that 64bit
//    types end up on 64bit address boundaries.
//
//
// The CODA driver libraries for these modules define structs
// that map the entire register space of the module (e.g.
// "fadc_struct" in fadcLib.h). We are only interested in a
// subset of those values so we define different structures here.
// Where possible, these use the same names for data members
// that the CODA drivers used.
//
// An instance of the ModulesConfigBOR structure is included
// in the ROC shared memory section defined in shmem_roc.h .
// The values are filled from the ROL (e.g. vme_list.c) using
// calls to the bor_utils.cc library in the rol_1 source directory.
//
// Access to the ModulesConfigBOR record from outside of the ROC
// is made via the shared memory server. The hdBOR program is 
// run at the start of data taking to gather this information,
// format it into EVIO and write it to the data stream. Flags
// are set to tell the Event Recorder to write it to the
// beginning of every file it creates.
//
// The BOR is uppacked and put into JANA objects in the DAQ
// library when the data is read in. They show up in classes
// like "Df125BORConfig", and "DF1TDCBORConfig" that inherit 
// from these structures.
//
// Any questions regarding this system can be referred to:
//
// davidl@jlab.org
//

#ifndef _BOR_ROC_H_
#define _BOR_ROC_H_

#include <stdint.h>

typedef struct{
	
	uint32_t rocid;
	uint32_t slot;
	
	// This based on information in "Programming_the_FADCV2_11_25_13.pdf"
		
	uint32_t version;          // 0x0000
	uint32_t ctrl1;            // 0x0008
	uint32_t ctrl2;            // 0x000C
	uint32_t blk_level;        // 0x0010
	uint32_t delay;            // 0x0024
	uint32_t itrig_cfg;        // 0x0028
	uint16_t dac[16];          // 0x0050 - 0x006C
	uint32_t trig21_delay;     // 0x0088
	uint32_t serial_number[3]; // 0x00E4 - 0x00EC
	
	uint32_t adc_status[3];    // 0x0100 - 0x0108
	uint32_t adc_config[4];    // 0x010C - 0x0118
	uint32_t adc_ptw;          // 0x011C
	uint32_t adc_pl;           // 0x0120
	uint32_t adc_nsb;          // 0x0124
	uint32_t adc_nsa;          // 0x0128
	uint16_t adc_thres[16];    // 0x012C - 0x0148
	uint32_t adc_pedestal[16]; // 0x0158 - 0x0194
	uint32_t config6;          // 0x014C
	uint32_t config7;          // 0x0150
	uint32_t config3;          // 0x0198
}f250config;

typedef struct {

	uint32_t version;             // 0xN000
	uint32_t nw;                  // 0xN058
	uint32_t pl;                  // 0xN05C
	uint32_t threshold[6];        // 0xN070
	uint32_t config1;             // 0xN088
	uint32_t config2;             // 0xN090
	uint32_t ped_sf;              // 0xN0A0
	uint32_t timing_thres_lo[3];  // 0xN0A0 - 0xN0AC
	uint32_t ie;                  // 0xN0B0
	uint32_t timing_thres_hi[2];  // 0xN0B4 - 0xN0B8
}f125config_fe;

typedef struct {
	
	uint32_t rocid;
	uint32_t slot;
	
	// This based on information in "fadc125_register_map_9.pdf"
	// and info found in fa125Lib.h
		
	uint32_t board_id;        // 0x0000
	uint32_t version;         // 0x0008
	uint32_t clock;           // 0x000C
	uint32_t serial[4];       // 0x0020 - 0x002C
	uint32_t temperature[2];  // 0x0030 - 0x0034
	uint32_t ctrl1;           // 0x0044
	
	f125config_fe fe[12];
	
	uint32_t proc_version;    // 0xD000
	uint32_t proc_csr;        // 0xD004
	uint32_t proc_trigsrc;    // 0xD008
	uint32_t proc_ctrl2;      // 0xD00C
	uint32_t proc_blocklevel; // 0xD014
}f125config;

typedef struct {

	uint32_t rocid;
	uint32_t slot;

	uint32_t version;         // 0x0000
	uint32_t ctrl;            // 0x0008
	uint32_t blocklevel;      // 0x0010
	
	uint32_t nchips;
	uint32_t f1registers[8][16]; 
}F1TDCconfig;


typedef struct {

	uint32_t rocid;
	uint32_t slot;

	uint32_t edge_resolution;
	uint32_t double_hit_resolution;
	uint32_t max_hits_per_event;
	uint16_t almostFullLevel; // 0x1022
	uint16_t bltEventNumber;  // 0x1024
	uint16_t firmwareRev;     // 0x1026
	
	
}caen1190config;

// Container struct that can hold all module
// configurations for this crate in a statically
// sized memory space appropriate for using with
// a shared memory section.
typedef struct {
	uint32_t       Nf250;
	uint32_t       Nf125;
	uint32_t       NF1TDC;
	uint32_t       Ncaen1190;

	f250config     f250[21];
	f125config     f125[21];
	F1TDCconfig    F1TDC[21];
	caen1190config caen1190[21];
}ModulesConfigBOR;

#endif //_BOR_ROC_H_
