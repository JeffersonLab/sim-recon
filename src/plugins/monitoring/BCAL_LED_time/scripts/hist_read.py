"""
python script to loop over root files in specified folder and process hist_read.C for each run
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


filelist = sorted(os.listdir(indir))

ndx = 0
for file in filelist:
    print file
    command = "root -b -q 'hist_read.C(\""+file[:9]+"\")'"
    print command
    if ndx == 0:
        os.system(command)
    ndx += 1







