# 

import os,sys
import pprint
from optparse import OptionParser

from ROOT import TFile,TH1F
import ccdb
from ccdb import Directory, TypeTable, Assignment, ConstantSet

CCDB_TABLE_NAME = { "SC":   "/START_COUNTER/base_time_offset",
                    "TOF":  "/TOF/base_time_offset",
                    "FCAL": "/FCAL/base_time_offset",
                    "BCAL": "/BCAL/base_time_offset",
                    "FDC":  "/FDC/base_time_offset",
                    "CDC":  "/CDC/base_time_offset",
                    "PS": "/PHOTON_BEAM/pair_spectrometer/base_time_offset",
                    "TAGH": "/PHOTON_BEAM/hodoscope/base_time_offset",
                    "TAGM": "/PHOTON_BEAM/microscope/base_time_offset"  }


def LoadCCDB():
    sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
    provider = ccdb.AlchemyProvider()                        # this class has all CCDB manipulation functions
    provider.connect(sqlite_connect_str)                     # use usual connection string to connect to database
    provider.authentication.current_user_name = "sdobbs"  # to have a name in logs

    return provider

def main():
    pp = pprint.PrettyPrinter(indent=4)
    
    #BEGINRUN = 40000
    ENDRUN = ccdb.INFINITE_RUN
    VARIATION = "calib"
    DRY_RUN = False
    SHIFT_VAL = 2.

    # Define command line options
    parser = OptionParser(usage = "check_2ns_shift.py [run] [root_file]")
    parser.add_option("-F","--run_file", dest="run_file", 
                      help="File of runs to look at")
    parser.add_option("-V","--variation", dest="variation", 
                      help="CCDB variation to use")
    parser.add_option("-e","--end_run", dest="end_run",
                      help="Ending run when scanning RCDB")
    parser.add_option("-Y","--dry_run", dest="dry_run", action="store_true",
                      help="Don't actually make any changes")
    
    (options, args) = parser.parse_args(sys.argv)

    if(len(args) < 3):
        parser.print_help()
        sys.exit(0)

    BEGINRUN = int(args[1])
    ROOT_FILE = args[2]

    if options.variation:
        VARIATION = options.variation
    if options.end_run:
        ENDRUN = int(options.end_run)
    if options.dry_run:
        DRY_RUN = options.dry_run

    # check 2ns shift:  assume that if both the SC and TOF are shifted, then everything is
    # we'll assume the shift is -2ns for now
    f = TFile(ROOT_FILE)
    SC_RF_T_Hist = f.Get("HLDetectorTiming/TRACKING/SC - RF Time")
    TOF_RF_T_Hist = f.Get("HLDetectorTiming/TRACKING/TOF - RF Time")

    # fit to Start Counter times
    try:
        maximum = SC_RF_T_Hist.GetBinCenter(SC_RF_T_Hist.GetMaximumBin())
        fr = SC_RF_T_Hist.Fit("gaus", "SQ", "", maximum - 0.3, maximum + 0.3)
        sc_mean = fr.Parameter(1)
        #print "shift = %6.3f"%mean
    except:
        print "bad fit to SC, quitting..."
        sys.exit(0)

    # fit to TOF times
    try:
        maximum = TOF_RF_T_Hist.GetBinCenter(TOF_RF_T_Hist.GetMaximumBin())
        fr = TOF_RF_T_Hist.Fit("gaus", "SQ", "", maximum - 0.3, maximum + 0.3)
        tof_mean = fr.Parameter(1)
        #print "shift = %6.3f"%mean
    except:
        print "bad fit to TOF, quitting..."
        sys.exit(0)


    # now see if both shifts are in the 2ns range
    if abs(sc_mean)<1.7 or abs(sc_mean)>2.3 or abs(tof_mean)<1.7 or abs(tof_mean)>2.3:
        print "No 2ns shift!"
        sys.exit(0)

    print "Applying overall 2ns shift ..."
    if sc_mean<0.:
        SHIFT_VAL = -SHIFT_VAL
    
    if DRY_RUN:
        sys.exit(0)

    # apply the shifts
    ccdb_conn = LoadCCDB()
    DETECTOR_SYSTEMS = [ "SC", "TOF", "FCAL", "BCAL", "FDC", "CDC", "TAGH", "TAGM", "PS" ]
    for detector in DETECTOR_SYSTEMS:
        off_assignment = ccdb_conn.get_assignment(ccdbtable, run, VARIATION)
        base_offsets = off_assignment.constant_set.data_table
        for x in xrange(len(base_offsets[0])):
            base_offsets[0][x] = "%7.3f"%(float(base_offsets[0][x]) - SHIFT_VAL)
        ccdb_conn.create_assignment(
                data=base_offsets,
                path=CCDB_TABLE_NAME[detector],
                variation_name=VARIATION,
                min_run=run,
                max_run=ccdb.INFINITE_RUN,
                comment="Fixed calibrations due to 2ns shift")



## main function 
if __name__ == "__main__":
    main()
