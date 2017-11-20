# This script generates default tables for the TAGM to remove large offsets.
# The tables generated are fadc_time_offsets and tdc_time_offsets.
# A copy of the current table is used and only problem channels are replaced.

import os,sys
from os import listdir
from os.path import isfile, join
import argparse

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
	parser = argparse.ArgumentParser( description = 'Generates default values for bad timing constants' )
	#parser.add_argument( '-a', '--adc', metavar = 'adc', type = str, nargs = None,
	#			help = 'file containing row and column information for bad adc channels')
	#parser.add_argument( '-t', '--tdc', metavar = 'tdc', type = str, nargs = None,
	#			help = 'file containing row and column information for bad tdc channels')
	#parser.add_argument( '-r', '--run', metavar = 'run', type = str, nargs = None,
	#			help = 'runNumber to be modified' )
	parser.add_argument( '-v', '--var', metavar = 'variation', type = str, nargs = None,
				help = 'variation of CCDB table' )
	args = parser.parse_args()
	#adcfile = args.adc
	#tdcfile = args.tdc
	#runNumber = args.run
	variation = args.var
	
	baseDir = os.getcwd() + '/bad-adcs/'
	adcfiles = [f for f in listdir( baseDir ) if isfile( join( baseDir, f ) ) ]
	baseDir = os.getcwd() + '/bad-tdcs/'
	tdcfiles = [f for f in listdir( baseDir ) if isfile( join( baseDir, f ) ) ]
	colsfile = open('bad-cols.txt', 'r')

	if not os.path.exists( os.getcwd() + '/adc_offsets'):
		os.makedirs( os.getcwd() + '/adc_offsets')
	if not os.path.exists( os.getcwd() + '/tdc_offsets'):
		os.makedirs( os.getcwd() + '/tdc_offsets')

	runList = []
	print('Making adc files')
	for f in adcfiles:
		runNumber = f.split('-')[2].split('.')[0]
		if runNumber not in runList:
			runList.append( runNumber )
		MakeFile( 'bad-adcs/' + f, runNumber, 'adc', variation)
	print('Making tdc files')
	for f in tdcfiles:
		runNumber = f.split('-')[2].split('.')[0]
		if runNumber not in runList:
			runList.append( runNumber )
		MakeFile( 'bad-tdcs/' + f, runNumber, 'tdc', variation)
	print('Making col files')
	for line in colsfile:
		runNumber = line.split()[0]
		try:
			int(runNumber)
		except ValueError:
			continue
		if runNumber not in runList:
			MakeFile('', runNumber, 'adc', variation)
			MakeFile('', runNumber, 'tdc', variation)

	colsfile.close()

def MakeFile(filename, runNumber, filetype, variation):
	# Empty filename is for adjusting mis-labeled columns 101, 102
	if filename:
		badfile = open( filename, 'r')

	chans = []
	if filename:
		for line in badfile:
			row = line.split()[0]
			col = line.split()[1]
			try:
				int(row)
			except ValueError:
				continue
			chans.append( [row, col] )

	table = open( filetype + '_offsets/' + filetype + '_offsets-' + runNumber + '.txt', 'w')
	if filetype == 'adc':
		assignment = provider.get_assignment("/PHOTON_BEAM/microscope/fadc_time_offsets", int(runNumber), variation)
	else:
		assignment = provider.get_assignment("/PHOTON_BEAM/microscope/tdc_time_offsets", int(runNumber), variation)

	ind_cols = [9, 27, 81, 99]
	channel = 0
	for i in range(1, 103):
		offset = str(assignment.constant_set.data_table[channel][2])

		# Replace bad constants
		for j in range(len(chans)):
			if (int(chans[j][0]) == 0) and (int(chans[j][1]) == i):
				offset = '10'

		table.write( ' 0\t' + str(i) + '\t' + offset + '\n' )

		channel += 1

		if i in ind_cols:
			for j in range(1, 6):
				offset = str(assignment.constant_set.data_table[channel][2])

				# Replace bad constants
				for k in range(len(chans)):
					if (int(chans[k][0]) == j) and (int(chans[k][1]) == i):
						offset = '10'
			
				table.write( ' ' + str(j) + '\t' + str(i) + '\t' + offset + '\n' )
				channel += 1

	table.close()
	if filename:
		badfile.close()

if __name__ == "__main__":
	main()
