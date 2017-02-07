#! /usr/bin/python

# Python script to fit the timewalk curve for the TAGM
# Usage: python tw.py -b <filename> <run number>
# Created 5/3/16 barnes

from ROOT import *
import os,sys

def main():
	# Check for proper usage
	if (len(sys.argv) != 4):
		print 'Usage: python tw.py -b <filename> <run number>'
		return

	# Use command line arguments
	filename = sys.argv[2]
	run = int(sys.argv[3])

	# Check if the run used the old or new bias voltage scheme
	if (run < 11572):
		newV = False
	else:
		newV = True

	# Open input and output ROOT files
	rootfile = TFile.Open(str(filename))
	outfile = TFile.Open("results.root","recreate")
	outfile.cd()

	# If the histogram is empty use the summed output hist instead
	indCol = [9,27,81,99]
	base = "TAGM_TW/tdc-rf/h_dt_vs_pp_tdc_"
	for i in range(1,103):
		# Summed outputs
		h = rootfile.Get(base+str(i))
		h.Write()
		p = tw_corr(h,0,i,newV)
		p.Write()

		# Individual fiber readouts
		if i in indCol:
			for j in range(5):
				h = rootfile.Get(base+"ind_"+str(j+1)+"_"+str(indCol.index(i)+1))
				h.Write()
				if not (h.GetEntries() > 300):
					h = rootfile.Get(base + str(i)).Clone()
					h.SetName(base+"ind_"+str(j+1)+"_"+str(indCol.index(i)+1))
				p = tw_corr(h,j+1,i,newV)
                                p.Write()
			
	# Include defaults for columns 101 and 102
	#file1 = open('tw-corr.txt','a')
	#for i in range(2):
	#	file1.write('0   ' + str(101 + i) + '   ' + '1   ' + '-1   ' +
        #            	    '0   ' + '8   ' + '0\n')
	outfile.Close()

def tw_corr(h,row,col,newV):
	# Open files for writing constants
	if (row == 0 and col == 1):
		file1 = open('tw-corr.txt','w')
	else:
		file1 = open('tw-corr.txt','a')

	# Find the reference time difference
	py = h.ProjectionY()
	fit = py.Fit("gaus","sq")
	dtmean = fit.Parameters()[1]
	
	# Make timewalk fit function and apply to hist
	# New voltage scheme has larger pulse height, adjust the range if needed
	if (newV):
		f1 = TF1("f1","[0]+[1]*(1/(x+[3]) )**[2]",400,2000)
	else:
		f1 = TF1("f1","[0]+[1]*(1/(x+[3]) )**[2]",125,2000)
	f1.SetParameter(0,-1)
	f1.SetParameter(1,100)
	f1.SetParameter(2,0.7)
	f1.SetParameter(3,-90)
	f1.SetParName(0,"c0")
	f1.SetParName(1,"c1")
	f1.SetParName(2,"c2")
	f1.SetParName(3,"c3")

	p = h.ProfileX()
	fitResult = p.Fit("f1","sRWq")

	c0 = fitResult.Parameters()[0]
	c1 = fitResult.Parameters()[1]
	c2 = fitResult.Parameters()[2]
	c3 = fitResult.Parameters()[3]

	# Write constants to file
	file1.write(str(row) + '   ' + str(col) + '   ' + str(c0) + '   ' + str(c1) + '   ' +
                    str(c2) + '   ' + str(c3) + '   ' + str(dtmean) + '\n')
	file1.close()

	return p


if __name__ == "__main__":
	main()
