#! /usr/bin/python

# This script will take the fits-thresh.out file and select the proper
# column based on command line input.

# Created 4/29/16 aeb

import sys
args = sys.argv

def main():
	infile = open('fits-thresh.out','r')
	outfile = open('tagm-thresh.txt','w')

	if (len(args) != 2):
		print "Usage: python MakeCCDB.py <# of sigmas from 1 to 5>"
		return
	if (int(args[1]) < 1 or int(args[1]) > 5):
		print "Valid multiples of sigma are 1, 2, 3, 4, or 5"
		return
	print "Sigma " + str(args[1]) + " selected"

	for line in infile:
		if (line.split()[0].isdigit()):
			row = int(line.split()[0])
			col = int(line.split()[1])
			thresh = float(line.split()[2+int(args[1])])
			thresh = 10**thresh
			outfile.write(str(row) + "   " + str(col) + "   " + str(thresh) + "   0\n")
					
	outfile.close()

if __name__ == "__main__":
	main()
