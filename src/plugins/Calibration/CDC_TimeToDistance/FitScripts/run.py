import os
import subprocess

dirs = ['Before', 'After', 'Combined', 'ccdb', 'Monitoring', 'Proj', 'ResVsT']
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

filelist = subprocess.check_output(["ls", "hists_merged"]).splitlines()
for file in filelist:
    subprocess.call(["root", "-l", "-b", "FitTimeToDistance.C(\"hists_merged/" + file + "\")"])
