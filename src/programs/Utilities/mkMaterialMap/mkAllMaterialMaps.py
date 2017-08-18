#!/usr/bin/env python
#
# This version of making all material maps will run all
# of mkMaterialMap processes in parallel. There are 
# about 17 maps produced, though some are done very fast.
# You should probably only use this on a computer with
# at least 2.5GB of FREE RAM.
#
import subprocess
import time
import os

procs = {}


def AddProc(mmap, cmd):
	wd  = 'dir_%s' % mmap
	os.mkdir(wd)
	p = subprocess.Popen(args=cmd.split(), cwd=wd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	procs[mmap] = p



# Target
AddProc('material_map00_target', 'mkMaterialMap -Nr 2 -Nz 2 -rmin 0 -rmax 1.8 -zmin 64.9 -zmax 65.1 -n_r 5 -n_z 5 -n_phi 10')

# Target wall
AddProc('material_map01_target_wall', 'mkMaterialMap -Nr 20 -Nz 2 -rmin 1.82 -rmax 1.9 -zmin 23.0 -zmax 65.1 -n_r 100 -n_z 5 -n_phi 10')

# Scattering chamber
AddProc('material_map02_scattering_chamber', 'mkMaterialMap -Nr 100 -Nz 400 -rmin 0.0 -rmax 4.7 -zmin 17.0 -zmax 86.1 -n_r 100 -n_z 5 -n_phi 10')

# Start Counter barrel
AddProc('material_map03_startcounter_barrel', 'mkMaterialMap -Nr 100 -Nz 100 -rmin 4.7 -rmax 9.0 -zmin 17.0 -zmax 86.1 -n_r 30 -n_z 30 -n_phi 100')

# Start Counter nose
AddProc('material_map04_startcounter_nose', 'mkMaterialMap -Nr 100 -Nz 50 -rmin 2.0 -rmax 9.0 -zmin 86.1 -zmax 97.5  -n_r 30 -n_z 30 -n_phi 20')

# CDC endplate
AddProc('material_map10_CDC_endplate', 'mkMaterialMap -Nr 60 -Nz 60 -rmin 9.0 -rmax 60.0 -zmin 167 -zmax 171 -n_r 10 -n_z 10 -n_phi 60')

# CDC inner shell
AddProc('material_map11_CDC_inner_shell', 'mkMaterialMap -Nr 50 -Nz 5 -rmin 9.0 -rmax 9.75 -zmin 17 -zmax 167 -n_r 1000 -n_z 50 -n_phi 10')

# CDC
AddProc('material_map12_CDC', 'mkMaterialMap -Nr 10 -Nz 10 -rmin 9.75 -rmax 56.0 -zmin 17 -zmax 167 -n_r 100 -n_z 100 -n_phi 200')

# CDC outer shell
AddProc('material_map13_CDC_outer_shell', 'mkMaterialMap -Nr 50 -Nz 5 -rmin 56.0 -rmax 65.0 -zmin 17 -zmax 167 -n_r 1000 -n_z 2 -n_phi 10')

# FDC package 1
AddProc('material_map21_FDC1', 'mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 174 -zmax 189 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC1-2 interpackage spacer
AddProc('material_map41_FDC_inter1', 'mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 189 -zmax 232.5 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC package 2
AddProc('material_map22_FDC2', 'mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 232.5 -zmax 247.5 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC2-3 interpackage spacer
AddProc('material_map42_FDC_inter2', 'mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 247.5 -zmax 291 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC package 3
AddProc('material_map23_FDC3', 'mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 291 -zmax 306 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC3-4 interpackage spacer
AddProc('material_map43_FDC_inter3', 'mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 306 -zmax 329.5 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC package 4
AddProc('material_map23_FDC4', 'mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 329.5 -zmax 344.5 -n_r 1000 -n_z 1000 -n_phi 10')

# FDC cables
AddProc('material_map31_cables', 'mkMaterialMap -Nr 30 -Nz 10 -rmin 61.0 -rmax 65.0 -zmin 167 -zmax 365 -n_r 1000 -n_z 1000 -n_phi 10')


print 'Waiting for all processes to complete ...'
finished = []
while True:
	Ndone = 0
	Ntotal = 0
	for lab in procs:
		Ntotal += 1
		if procs[lab].poll() is not None :
			Ndone += 1
			if lab not in finished:
				finished.append(lab)
				print ' finished: %s' % lab
				subprocess.call(['mv', 'dir_%s/material_map' % lab, '%s' % lab])
	if Ndone >= Ntotal: break
	print '%d/%d processes complete' % (Ndone,Ntotal)
	time.sleep(2)

print 'All processes complete'


# Make simple map with air everywhere as default
f = open('material_map99_course_default', 'w')
f.write("#    r       z       A      Z   density  radlen   rhoZ_overA  rhoZ_overA_logI  chi2c_factor  chi2a_factor  chi2a_corr\n")
f.write("#%  00      01      02     03        04      05           06               07            08            09          10\n")
f.write("   0.5    -49.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126\n")
f.write(" 119.5    -49.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126\n")
f.write("   0.5    649.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126\n")
f.write(" 119.5    649.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126\n")
f.close()


