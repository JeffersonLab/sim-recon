#!/usr/bin/python

# This script calculates the various timing offsets for the CCDB.
# Created by Alex Barnes, October 18, 2016

# Usage:
# python timing.py -b <root file> <run number> <calibration type>
# The calibration type can be: self or rf,
# where self refers to TDC-ADC
# Specifying the calibration type determines which histogram
# the script uses for the calibration.
##################################################################

import os,sys
from ROOT import *

# Set up CCDB connection for reading tables
import ccdb
from ccdb import Directory, TypeTable, Assignment, ConstantSet

# CCDB_HOME and CCDB_CONNECTION environment variables must be set!
ccdb_home = os.environ["CCDB_HOME"]
sqlite_connect_str = os.environ["CCDB_CONNECTION"]

# Create CCDB API class
## This class has all CCDB manipulation functions
provider = ccdb.AlchemyProvider()
## Use usual connection string to connect to database
provider.connect(sqlite_connect_str)
print getattr(provider,"is_connected")
## Provide a username for CCDB updates
provider.authentication.current_user_name = "anonymous"	

# Define TAGM constants
NCOLUMNS = 102
NROWS = 5

def main():
	if (len(sys.argv) != 5):
		print "\nUsage: python timing.py -b <root file> <run number> <calibration type>"
		print "The calibration type can be: self or rf,"
		print "where self refers to TDC-ADC\n"
		return

	rootfile = TFile.Open(str(sys.argv[2]))
	run = sys.argv[3]
	calib_type = str(sys.argv[4])
	adcfile = open('adc_offsets-' + str(run) + '.txt','w')
	tdcfile = open('tdc_offsets-' + str(run) + '.txt','w')

	baseDir = 'TAGM_TW/'

	beamPeriod_assignment = provider.get_assignment("/PHOTON_BEAM/RF/beam_period",run,"default")
	beamPeriod = float(beamPeriod_assignment.constant_set.data_table[0][0])
	print "The beam period is: " + str(beamPeriod)
	fadc_assignment = provider.get_assignment("/PHOTON_BEAM/microscope/fadc_time_offsets",run,"default")
	tdc_assignment = provider.get_assignment("/PHOTON_BEAM/microscope/tdc_time_offsets",run,"default")

	ind_cols = [9, 27, 81, 99]
	channel = 0
	for i in range(1,NCOLUMNS+1):
		# Get summed channels
		if (calib_type == 'self'):
			print 'Doing TDC-ADC timing'
			hist = rootfile.Get(baseDir+"t_adc_all").Clone().ProjectionX("_"+str(i),i,i)
			new_offset = GetSelfTiming(hist)
		elif (calib_type == 'rf'):
			print 'Doing RF timing'
			hist = rootfile.Get(baseDir+"adc_rf_all").Clone().ProjectionX("_"+str(i),i,i)
			new_offset = GetRFTiming(hist,beamPeriod)
		else:
			print "\nIncorrect calibration type. Please use self, ps, or rf.\n"
			return
		print 'new_offset: ' + str(new_offset)

		# Write fADC offset
		if (calib_type == 'self'):
			offset = float(tdc_assignment.constant_set.data_table[channel][2]) + new_offset
			tdcfile.write(' 0\t' + str(i) + '\t' + str(offset) + '\n')
		else:
			offset = float(fadc_assignment.constant_set.data_table[channel][2]) + new_offset
			print 'final offset: ' + str(offset)
			adcfile.write(' 0\t' + str(i) + '\t' + str(offset) + '\n')

		channel += 1

		# Get individual channels
		if i in ind_cols:
			for j in range(5):
				id = ind_cols.index(i)*5+j+1
				if (calib_type == 'self'):
					h_name = baseDir+"t_adc_all_ind"
					hist = rootfile.Get(h_name).Clone().ProjectionX("_"+str(id),id,id)
					new_offset = GetSelfTiming(hist)
				elif (calib_type == 'rf'):
					h_name = baseDir+"adc_rf_all_ind"
					hist = rootfile.Get(h_name).Clone().ProjectionX("_"+str(id),id,id)
					new_offset = GetRFTiming(hist,beamPeriod)
				else:
					print "\nIncorrect calibration type. Please use self, ps, or rf.\n"
					return

				# Write fADC offset
				if (calib_type == 'self'):
					offset = float(tdc_assignment.constant_set.data_table[channel][2]) + new_offset
					tdcfile.write(' ' + str(j+1) + '\t' + str(i) + '\t' + str(offset) + '\n')
				else:
					offset = float(fadc_assignment.constant_set.data_table[channel][2]) + new_offset
					adcfile.write(' ' + str(j+1) + '\t' + str(i) + '\t' + str(offset) + '\n')

				# For RF, adjust and write TDC offsets
				if (calib_type == 'rf'):
					offset = float(tdc_assignment.constant_set.data_table[channel][2]) - new_offset
					tdcfile.write(' ' + str(j+1) + '\t' + str(i) + '\t' + str(offset) + '\n')

				channel += 1


	adcfile.close()
	tdcfile.close()
	rootfile.Close()

def GetRFTiming(hist,beamPeriod):
	# This method will take in a histogram and return the histogram's offset from 0
	histname = hist.GetName()
	print histname
	mean = hist.GetMean()
	print 'Histogram mean: ' + str(mean)
	maxBin = hist.GetMaximumBin()
	maxContent = hist.GetBinContent(maxBin)
	maximum = hist.GetBinCenter(maxBin)
	print 'Histogram maximum: ' + str(maximum)
	try:
		if (maximum < -0.5):
			fitfunc = TF1("fitfunc","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))+[0]*TMath::Exp(-0.5*pow((x-[1]-[3])/[2],2))",-beamPeriod/2,0)
			fitMin = -beamPeriod/2
			fitMax = maximum + 0.5
			print 'fitMin: ' + str(fitMin)
			print 'fitMax: ' + str(fitMax)
			print 'beamPeriod: ' + str(beamPeriod)
		elif (maximum > 0.5):
			fitfunc = TF1("fitfunc","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))+[0]*TMath::Exp(-0.5*pow((x-[1]+[3])/[2],2))",0,beamPeriod/2)
			fitMin = maximum - 0.5
			fitMax = beamPeriod/2
			print 'fitMin: ' + str(fitMin)
			print 'fitMax: ' + str(fitMax)
		else:
			fitfunc = TF1("fitfunc","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))",-beamPeriod/2,beamPeriod/2)
			fitMin = maximum - 0.4
			fitMax = maximum + 0.4
			print 'fitMin: ' + str(fitMin)
			print 'fitMax: ' + str(fitMax)
		fitfunc.SetParameter(1,maximum)
		fitfunc.SetParLimits(1,-beamPeriod/2,beamPeriod/2)
		fitfunc.SetParameter(2,0.4)
		if (abs(maximum) > 0.5):
			fitfunc.FixParameter(3,beamPeriod)
		FitResult = hist.Fit(fitfunc,"sRq","",fitMin,fitMax)
			
		offset = FitResult.Parameters()[1]
		sigma = FitResult.Parameters()[2]
		print 'mean: ' + str(offset)
		print 'sigma: ' + str(sigma)
		if (abs(maximum) > 0.5):
			print '[3]: ' + str(FitResult.Parameters()[3])
	except:
		print "fit failed for histogram " + histname
		offset = 0
	
	return offset

def GetSelfTiming(hist):
	# This method will take in a histogram and return the histogram's offset from 0
	histname = hist.GetName()
	print histname
	mean = hist.GetMean()
	print 'Histogram mean: ' + str(mean)
	maxBin = hist.GetMaximumBin()
	maxContent = hist.GetBinContent(maxBin)
	maximum = hist.GetBinCenter(maxBin)
	print 'Histogram maximum: ' + str(maximum)
	try:
		FitResult = hist.Fit("gaus","sRWq","",-20,20)
		offset = FitResult.Parameters()[1]
		sigma = FitResult.Parameters()[2]
		print 'mean: ' + str(offset)
		print 'sigma: ' + str(sigma)
	except:
		print "fit failed for histogram " + histname
		offset = 0
	
	return offset

if __name__ == "__main__":
	main()
