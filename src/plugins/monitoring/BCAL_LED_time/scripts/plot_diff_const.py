#!/usr/bin/python
"""
plot_diff_const.py
read in two files of time const offsets and compare
"""
import os
import numpy
from matplotlib import pyplot as plt
from datetime import datetime, date, time, timedelta
from matplotlib.backends.backend_pdf import PdfPages
#from pylab import savefig
filelist= []

# Take means from previous runs
mean0 = 10.29
meanM4 = 6.20
runlist = []
runmax = 100

figs = []
nfig = 0


# open file for output
fout = open("dat/plot_diff_const.err","w")

filelist = sorted(os.listdir("dat"))
# filelist = ["Run030660.txt"]
print filelist
ndx = 1
for run in filelist:
    #  Initialize arrays
    Chid = []
    Const1 = []
    Const2 = []
    Const0 = []         # create arrays to assign time offsets
    ConstM4 = []


    runno = 0
    runlist.append(runno)
    filename = "dat/"+run
    # print "filename=",filename
    if ".dat" in filename:
		continue
    if ".err" in filename:
        continue
    if ".py" in filename:
        continue
    if ".txte" in filename:
        continue
    runno = int(run[4:9])
    file = open (filename,"r")
    nfig += 1
    # print "filename=",filename," runno=",runno, " nfig=",nfig
    if nfig > runmax or runno > 30789:
        continue
    ndx = -1
    for line in file:
        linew = line.split()
        # print "Linew=",int(linew[0]),float(linew[1])
    	ndx += 1
    	if (ndx%2 == 0):
            continue
    	Chid.append(int(linew[0]))
    	if (float(linew[1]) > 8):
            Const0.append(float(linew[1]))
            Const1.append(float(linew[1])-mean0+0)
    	else:
            Const1.append(float(linew[1])-meanM4+4)
            ConstM4.append(float(linew[1]))

    file.close()
    for (ndx,id) in enumerate(Chid):
        value = Const1[ndx]
        # print id, value

    const0 = numpy.array(Const0)
    constM4 = numpy.array(ConstM4)
    if abs(mean0-numpy.mean(const0)) > 0.1 or abs(meanM4-numpy.mean(constM4)) > 0.1:
        print " *** mean0=", mean0, numpy.mean(const0), " meanM4=",meanM4, numpy.mean(constM4)

    file = "ccdb/"+str(runno)+".txt"
    file = open(file,"r")
    ndx = -1
    for line in file:
        linew = line.split()
    	ndx += 1
    	if ndx%2 == 0:
            continue
    	Const2.append(float(linew[0]))
    file.close()
    for (ndx,id) in enumerate(Chid):
        value = Const2[ndx]
    	# print id, value

    const1 = numpy.array(Const1)
    const2 = numpy.array(Const2)
    constM4 = numpy.array(ConstM4)
    constdiff = numpy.zeros(len(const1))
	
    for (ndx,id) in enumerate (Chid):
        module = Chid[ndx]/32 + 1
        layer = (Chid[ndx] - (module-1)*32)/8 + 1
        sector = (Chid[ndx] - (module-1)*32 - (layer-1)*8)/2 +1
        end = Chid[ndx] - (module-1)*32 - (layer-1)*8 - (sector-1)*2 + 1
        constdiff[ndx] = const1[ndx] - const2[ndx]
        if constdiff[ndx] < -1.5 or constdiff[ndx] > 1.5:
            print '{0} {1} {2:5} {3:5} {4:3} {5:3} {6:5.2f} {7:5.2f} {8:5.2f} {9}'.format("*** Run", runno, Chid[ndx],module, layer, sector,  const1[ndx], const2[ndx], constdiff[ndx],"\n")
            fout.write('{0} {1} {2:5} {3:5} {4:3} {5:3} {6:5.2f} {7:5.2f} {8:5.2f} {9}'.format("Run", runno, Chid[ndx],module, layer, sector,  const1[ndx], const2[ndx], constdiff[ndx],"\n"))
        # else:
        #    print '{0} {1} {2:5} {3:5} {4:3} {5:3} {6:5.2f} {7:5.2f} {8:5.2f} {9}'.format("Run", runno, Chid[ndx],module, layer, sector,  const1[ndx], const2[ndx], constdiff[ndx],"\n")

    fig = plt.figure(nfig)
    figs.append(fig)
    plt.subplot(2,2,1)
    plt.title("LED "+"Run"+str(runno))
    nbins = 100
    n, bins, patches = plt.hist(const1,nbins,log=True)
    plt.xlabel("Time Offset (ns)")
    plt.ylabel("Entries")
    plt.axis([-15,25,0.1,2000])
    # plt.autoscale(enable=True,axis='y')

    plt.subplot(2,2,2)
    plt.title("Data")
    nbins = 100
    n, bins, patches = plt.hist(const2,nbins,log=True)
    plt.xlabel("Time Offset (ns)")
    plt.ylabel("Entries")
    plt.axis([-15,25,0.1,2000])
    # plt.autoscale(enable=True,axis='y')

    plt.subplot(2,2,3)
    # plt.title("Difference")
    nbins = 100
    n, bins, patches = plt.hist(constdiff,nbins,log=True)
    plt.xlabel("Difference (ns)")
    # plt.ylabel("Entries")
    plt.axis([-15,25,0.1,1000])
    # plt.autoscale(enable=True,axis='y')

    plt.subplot(2,2,4)
    # plt.title("Correlations")
    plt.plot(const2,const1,'r.')
    plt.xlabel("Data")
    plt.ylabel("LED")
    plt.axis([-1,5,-5,25])
    # plt.autoscale(enable=True)

fout.close()


# plt.show()
with PdfPages("pdf/"+"plot_diff_summary.pdf") as pdf:
    for fig in figs:
        print "Page in pdf=",fig
        pdf.savefig(fig)

