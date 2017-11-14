#!/usr/bin/python

# This script finds bad constants for the TAGM
# Created by Alex Barnes, November 13, 2017

##################################################################

import os,sys

# Set up CCDB connection for reading tables
import ccdb, rcdb
from ccdb import Directory, TypeTable, Assignment, ConstantSet

# CCDB_HOME and CCDB_CONNECTION environment variables must be set!
ccdb_home = os.environ["CCDB_HOME"]
sqlite_connect_str = os.environ["CCDB_CONNECTION"]

# Create CCDB API class
## This class has all CCDB manipulation functions
provider = ccdb.AlchemyProvider()
## Use usual connection string to connect to database
provider.connect(sqlite_connect_str)
print 'CCDB is connected? ' + str(getattr(provider,"is_connected"))
## Provide a username for CCDB updates
provider.authentication.current_user_name = "anonymous"	

def main():
	db = rcdb.RCDBProvider("mysql://rcdb@hallddb.jlab.org/rcdb")
	rcdbQuery = '@is_production and @status_approved'
	rcdbRunList = db.select_runs(rcdbQuery, 30274, 31057)

	dirs = ['bad-adcs', 'bad-tdcs']
	for i in dirs:
		if not os.path.exists(i):
			os.makedirs(i)

	colfile = open('bad-cols.txt','w')
	colfile.write('Problem runs\n')

	variation = 'default'

	for entry in rcdbRunList:
		run = entry.number
		adcfile = open('bad-adcs/bad-adc-' + str(run) + '.txt', 'w')
		adcfile.write('Row\tCol\n')

		tdcfile = open('bad-tdcs/bad-tdc-' + str(run) + '.txt', 'w')
		tdcfile.write('Row\tCol\n')


		# check for wrong 101 and 102 columns
		fadc_assignment = provider.get_assignment("/PHOTON_BEAM/microscope/fadc_time_offsets",run,variation)
		adc_101 = float(fadc_assignment.constant_set.data_table[120][1])
		adc_102 = float(fadc_assignment.constant_set.data_table[121][1])

		tdc_assignment = provider.get_assignment("/PHOTON_BEAM/microscope/tdc_time_offsets",run,variation)
		tdc_101 = float(tdc_assignment.constant_set.data_table[120][1])
		tdc_102 = float(tdc_assignment.constant_set.data_table[121][1])
		if (adc_101 == 120 or adc_102 == 121 or tdc_101 == 120 or tdc_102 == 121):
			colfile.write(str(run) + '\n')

		# check for large adc and tdc offsets
		ind_cols = [9, 27, 81, 99]
		channel = 0
		for i in range(1,103):
			adc = float(fadc_assignment.constant_set.data_table[channel][2])
			if (abs(adc) >= 50.0):
				adcfile.write( '0\t' + str(i) + '\n')
			tdc = float(tdc_assignment.constant_set.data_table[channel][2])
			if (abs(tdc) >= 50.0):
				tdcfile.write( '0\t' + str(i) + '\n')
			channel += 1
			
			if i in ind_cols:
				for j in range(5):
					adc = float(fadc_assignment.constant_set.data_table[channel][2])
					if (abs(adc) >= 50.0):
						adcfile.write( str(j+1) + '\t' + str(i) + '\n')
					tdc = float(tdc_assignment.constant_set.data_table[channel][2])
					if (abs(tdc) >= 50.0):
						tdcfile.write( str(j+1) + '\t' + str(i) + '\n')
					channel += 1
		adcfile.close()
		tdcfile.close()
				

if __name__ == "__main__":
	main()
