#! /usr/bin/python

# Python script to fit the timewalk curve for the PSC
# Usage: python tw.py -b <filename> <run number>
# Created 5/3/16 barnes

from ROOT import *
import os,sys

def main():
	# Check for proper usage
	if (len(sys.argv) != 3):
		print 'Usage: python tw.py -b <filename>'
		return

	# Use command line arguments
	filename = sys.argv[2]

	# Open input and output ROOT files
	rootfile = TFile.Open(str(filename))
	outfile = TFile.Open("results.root","recreate")
	outfile.cd()

	# If the histogram is empty use the summed output hist instead
	base = "PSC_TW/tdc-rf/h_dt_vs_pp_tdc_"
	for i in range(1,9):
		# Left modules
		h = rootfile.Get(base+'l_'+str(i))
		h.Write()
		p = tw_corr(h,0,i)
		p.Write()
		# Right modules
		h = rootfile.Get(base+'r_'+str(i))
		h.Write()
		p = tw_corr(h,1,i)
		p.Write()

	outfile.Close()

def tw_corr(h,side,module):
	# Open files for writing constants
	if (side == 0 and module == 1):
		file1 = open('tw-corr.txt','w')
	else:
		file1 = open('tw-corr.txt','a')

	# Get MPV
	try:
		p = h.ProjectionX()
		fitResult = p.Fit("landau","sq")
		MPV = fitResult.Parameters()[1]
	except:
		MPV = 0

	# Use TH2 to create a cleaner TH2 for timewalk corrections
	# by fitting pulse peak slices to get the average time and
	# associated error.
	hist = h.ProjectionX().Clone()
	hist.Reset()
	try:
		nbinsX = h.GetNbinsX()
		for slice in range(100):
			peak = (2*10*slice + 10)
			p = h.ProjectionY("_py",10*slice,10*slice+10)
			maxbin = p.GetMaximumBin()
			if maxbin == 1 or p.GetEntries() < 20:
				continue
			maxval = p.GetBinCenter(maxbin)
			fitResult = p.Fit("gaus", "sq", "", maxval-0.5, maxval+0.5)
			mean = fitResult.Parameters()[1]
			sig = fitResult.Parameters()[2]
			hist.Fill(peak, mean)
			bin = hist.FindBin(peak)
			hist.SetBinError(bin,sig)
	except:
		if side == 0:
			print 'Histogram empty for left module ' + str(module)
		else:
			print 'Histogram empty for right module ' + str(module)

	# Make timewalk fit function and apply to hist
	# New voltage scheme has larger pulse height, adjust the range if needed
	try:
		f1 = TF1("f1","[0]+[1]*(x/16)**[2]",100,2000)
		f1.SetParameter(0,-1)
		f1.SetParameter(1,5)
		f1.SetParameter(2,-0.5)
		f1.SetParName(0,"c0")
		f1.SetParName(1,"c1")
		f1.SetParName(2,"c2")

		fitResult = hist.Fit("f1","sRWq")

		c0 = fitResult.Parameters()[0]
		c1 = fitResult.Parameters()[1]
		c2 = fitResult.Parameters()[2]

	except:
		c0 = 1
		c1 = -1
		c2 = 0

	# Write constants to file
	file1.write(str(side) + '   ' + str(module) + '   ' + str(c0) + '   ' + str(c1) + '   ' +
                    str(c2) + '   ' + str(16) + '   ' + str(MPV) + '\n')
	file1.close()

	return hist


if __name__ == "__main__":
	main()
