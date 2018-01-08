"""
python script to loop over root files in specified folder and process Read_bcal_hadronic_eff.C for each run and layers 1-4
"""


import subprocess
import numpy
from matplotlib import pyplot as plt
from datetime import datetime, date, time, timedelta
from matplotlib.backends.backend_pdf import PdfPages
import os
from os.path import isfile, join

filenames = []
runs = []

minrunno = 30000
maxrunno = 31000
maxruns = 10
minfiles = 1
maxfiles = 100


filelist = sorted(os.listdir("/Users/elton/scratch/Bcal_hadronic_eff/"))

for file in filelist:
    if ".root" in file:
        for layer in range(1,5):
            filerun = file.replace("tree_bcal_hadronic_eff_","").replace(".root","")
            option = '{0} {1:d}'.format(filerun,layer)
            command = "root -b -q "+file+" 'call_bcal_hadronic_eff.C(\"Read_bcal_hadronic_eff2.C\",\""+option+"\")'"
            print command
            os.system(command)








