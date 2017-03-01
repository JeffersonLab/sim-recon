#!/usr/bin/python

# Usage: python thresholds.py -b rootfile

from ROOT import *
import sys
args = sys.argv

fitsAll = open('fits-thresh.out','w')
#fitsGaus.write("col   mean   1 sig   FWHM/2   2 sig\n")
fitsAll.write("row   col   mean   1 sig   2 sig   3 sig   4 sig   5 sig\n")

def main():
	#outfile = MakeOutputFile("results-thresh.root")
	outfile = TFile.Open("results-thresh.root","recreate")

	#rootfile = 'hd_rawdata_010390.root'
	rootfile = str(args[2])

	final_gaus = [ [] for i in range(3) ]

	#GetPeaks(rootfile,outfile)
	FitAll(rootfile,outfile,final_gaus)

	#fitsGaus.close()
	fitsAll.close()

def GetPeaks(rootfile,outfile):
	infile = TFile.Open(rootfile)
	outfile.cd()
	
	c = TCanvas("c","c",0,20,600,500)

	for i in range(1,101):
		h = infile.Get("h_pint_"+str(i))

		if (h.GetEntries() < 1):
			continue
		else:
			fitResult = h.Fit("gaus","sq","",3.0,3.5)

		if not (fitResult.IsEmpty()):
			const = fitResult.Parameters()[0]
			mean = fitResult.Parameters()[1]
			sig = fitResult.Parameters()[2]
		else: 
			const = 0
			mean = 0
			sig = 0

		h.Write()

		fitsGaus.write("   " + str(i) + "   " + str(round(mean,3)) + "   " + 
                               str(round(mean+sig,3)) + "   " + str(round(mean+1.1775*sig,3)) +
                               "   " + str(round(mean+2*sig,3)) + '\n')

def GetLandau(rootfile,row,outfile,init_land,init_gaus):
	infile = TFile.Open(rootfile)
	outfile.cd("row" + str(row))
	
	c = TCanvas("c","c",0,20,600,500)

	for i in range(1,101):
		element = (int(row)-1)*100 + (i-1)
		h = infile.Get("h_pint_"+str(i))
		if (h.GetEntries() > 1):
			#fitResult = h.Fit("landau","sq","",0,maxRange)
			#fitResult = h.Fit("landau","sq","",0,400)
			#fitResult = h.Fit("landaun","sq","",100,500)
			fitResult = h.Fit("landau","sq","",2.5,3.3)
			#fitResult = h.Fit("landau","sq","",0,40)
		else:
			const = 0
			mpv = 0
			sig = 0
			init_land[0].append(const)
			init_land[1].append(mpv)
			init_land[2].append(sig)
			h.Write()
			fitsLand.write("   " + str(row) + "   " + str(i) + "   " + str(round(mpv,3)) + '\n')
			continue
		
		if not (fitResult.IsEmpty()):
			const = fitResult.Parameters()[0]
			mpv = fitResult.Parameters()[1]
			sig = fitResult.Parameters()[2]
		elif (fitResult.IsEmpty()) or (fitResult.Parameters()[1] < 0):
			const = 0
			mpv = 0
			sig = 0
		fitsLand.write("   " + str(row) + "   " + str(i) + "   " + str(round(mpv,3)) + '\n')

		init_land[0].append(const)
		init_land[1].append(mpv)
		init_land[2].append(sig)

		h.GetXaxis().SetRangeUser(2.5,5)
		h.Write()

def GetExpo(rootfile,row,outfile,init_land,init_gaus):
	infile = TFile.Open(rootfile)
	outfile.cd("row" + str(row))
	
	c = TCanvas("c","c",0,20,600,500)

	for i in range(1,101):
		element = (int(row)-1)*100 + (i-1)
		h = infile.Get("h_pint_"+str(i))
		bin4 = h.FindBin(3.8)

		if (h.GetEntries() > 1):
			if (h.GetBinContent(bin4) < 100):
				fitResult = h.Fit("expo","sq","",3.0,3.25)
			else:
				fitResult = h.Fit("expo","sq","",3.0,3.3)
		else:
			const = 0
			mpv = 0
			init_land[0].append(const)
			init_land[1].append(mpv)
			h.Write()
			fitsLand.write("   " + str(row) + "   " + str(i) + "   " + str(round(mpv,3)) + '\n')
			continue
		
		if not (fitResult.IsEmpty()):
			const = fitResult.Parameters()[0]
			mpv = fitResult.Parameters()[1]
		elif (fitResult.IsEmpty()) or (fitResult.Parameters()[1] < 0):
			const = 0
			mpv = 0
		fitsLand.write("   " + str(row) + "   " + str(i) + "   " + str(round(mpv,3)) + '\n')

		init_land[0].append(const)
		init_land[1].append(mpv)

		h.GetXaxis().SetRangeUser(2.5,5)
		h.Write()

def FitAll(rootfile,outfile,final_gaus):
	infile = TFile.Open(rootfile)
	outfile.cd()
	
	c = TCanvas("c","c",0,20,600,500)

	function = TF1("function","gaus(0)+expo(3)",2.5,4.1)

	baseName = 'TAGM_thresh/integrals/'
	for i in range(1,103):
		#h = infile.Get("h_pint_"+str(i))
		h = infile.Get(baseName+"h_int_"+str(i))
		bins = h.GetNbinsX()
		if (i < 43):
			h.GetXaxis().SetRangeUser(2.2,2.9)
		else:
			h.GetXaxis().SetRangeUser(2.5,3.05)
		fitmin = h.GetBinCenter( h.GetMinimumBin() )
		h.GetXaxis().SetRangeUser(fitmin,4)
		fitmean = h.GetBinCenter( h.GetMaximumBin() )
		h.GetXaxis().SetRangeUser(2,5)
		try:
			function.SetParameter(0,40)
			function.SetParLimits(0,0,1000)
			function.SetParameter(1,fitmean)
			function.SetParLimits(1,2.7,4.0)
			function.SetParameter(2,0.1)
			function.SetParLimits(2,0.0,1.0)
			function.SetParameter(3,10)
			function.SetParameter(4,-2.5)
			function.SetParLimits(4,-100,0)
			fitResult = h.Fit(function,"sqRW","",fitmin,fitmean+0.2)
			constGaus = fitResult.Parameters()[0]
			meanGaus = fitResult.Parameters()[1]
			sigGaus = fitResult.Parameters()[2]
			constLand = fitResult.Parameters()[3]
			mpvLand = fitResult.Parameters()[4]
			final_gaus[0].append(constGaus)
			final_gaus[1].append(meanGaus)
			final_gaus[2].append(sigGaus)
			fitsAll.write("0   " + str(i) + "   " + str(round(meanGaus,3)) + "   " + str(round(meanGaus-sigGaus,3)) + 
                	              "   " + str(round(meanGaus-2*sigGaus,3)) + "   " + str(round(meanGaus-3*sigGaus,3)) +
                	              #"   " + str(round(meanGaus-4*sigGaus,3)) + "   " + str(round(meanGaus-5*sigGaus,3)) + '\n')
                	              "   " + str(round(meanGaus-4*sigGaus,3)) + "   " + str(round(meanGaus-7*sigGaus,3)) + '\n')

			h.Write()
		#else:
		except:
			constGaus = 0
			meanGaus = 0
			sigGaus = 0
			constLand = 0
			mpvLand = 0
			final_gaus[0].append(constGaus)
			final_gaus[1].append(meanGaus)
			final_gaus[2].append(sigGaus)
			fitsAll.write("0   " + str(i) + "   " + str(round(meanGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + '\n')
			h.Write()
			continue

		for j in range(1,6):
			#if i == 9: col = 1
			#elif i == 27: col = 2
			#elif i == 81: col = 3
			#elif i == 99: col = 4
			if i == 9: j = j
			elif i == 27: j = j+5
			elif i == 81: col = j+10
			elif i == 99: col = j+15
			else: continue
			#h = infile.Get(baseName+"h_int_ind_"+str(j)+"_"+str(col))
			h = infile.Get(baseName+"h_int_ind_"+str(j))
			bins = h.GetNbinsX()
			if (i < 43):
				h.GetXaxis().SetRangeUser(2.2,2.9)
			else:
				h.GetXaxis().SetRangeUser(2.5,3.5)
			fitmin = h.GetBinCenter( h.GetMinimumBin() )
			h.GetXaxis().SetRangeUser(fitmin,4)
			fitmean = h.GetBinCenter( h.GetMaximumBin() )
			h.GetXaxis().SetRangeUser(2,5)
			try:
				function.SetParameter(0,40)
				function.SetParLimits(0,0,1000)
				function.SetParameter(1,fitmean)
				function.SetParLimits(1,2.7,4.0)
				function.SetParameter(2,0.1)
				function.SetParLimits(2,0.0,1.0)
				function.SetParameter(3,10)
				function.SetParameter(4,-2.5)
				fitResult = h.Fit(function,"sqRW","",fitmin,fitmean+0.2)
				constGaus = fitResult.Parameters()[0]
				meanGaus = fitResult.Parameters()[1]
				sigGaus = fitResult.Parameters()[2]
				constLand = fitResult.Parameters()[3]
				mpvLand = fitResult.Parameters()[4]
				final_gaus[0].append(constGaus)
				final_gaus[1].append(meanGaus)
				final_gaus[2].append(sigGaus)
				fitsAll.write(str(j) + "   " + str(i) + "   " + str(round(meanGaus,3)) + "   " + str(round(meanGaus-sigGaus,3)) + 
                		              "   " + str(round(meanGaus-2*sigGaus,3)) + "   " + str(round(meanGaus-3*sigGaus,3)) +
                		              #"   " + str(round(meanGaus-4*sigGaus,3)) + "   " + str(round(meanGaus-5*sigGaus,3)) + '\n')
                		              "   " + str(round(meanGaus-4*sigGaus,3)) + "   " + str(round(meanGaus-7*sigGaus,3)) + '\n')

				h.Write()
			#else:
			except:
				constGaus = 0
				meanGaus = 0
				sigGaus = 0
				constLand = 0
				mpvLand = 0
				final_gaus[0].append(constGaus)
				final_gaus[1].append(meanGaus)
				final_gaus[2].append(sigGaus)
				fitsAll.write(str(j) + "   " + str(i) + "   " + str(round(meanGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + "   " + str(round(sigGaus,3)) + '\n')
				h.Write()
				continue
		
def MakeOutputFile(rootfile):
	outfile = TFile.Open(rootfile,"RECREATE")
	outfile.mkdir("row1")
	outfile.mkdir("row2")
	outfile.mkdir("row3")
	outfile.mkdir("row4")
	outfile.mkdir("row5")

	return outfile

def MakeFinalGaus(rootfile,row,outfile,final_gaus,final_land):
	infile = TFile.Open(rootfile)
	outfile.cd("row" + str(row))


	for i in range(1,101):
		element = (int(row)-1)*100 + (i-1)

		c = TCanvas("c"+str(i),"c"+str(i)+"_"+str(row),0,20,600,500)
		gPad.SetLogy()

		h = infile.Get("h_spectra_"+str(i))
		h.Rebin(16)
		h.Draw()

		f1 = TF1("f1","gaus",0,2000)
		f1.SetParameter(0,final_gaus[0][element])
		f1.SetParameter(1,final_gaus[1][element])
		f1.SetParameter(2,final_gaus[2][element])
		f1.Draw("same")

		#f2 = TF1("f2","expo",0,2000)
		#f2 = TF1("f2","landaun",0,2000)
		f2 = TF1("f2","landau",0,2000)
		f2.SetParameter(0,final_land[0][element])
		f2.SetParameter(1,final_land[1][element])
		f2.SetParameter(2,final_land[2][element])
		f2.Draw("same")

		c.Write()

if __name__ == "__main__":
	main()
