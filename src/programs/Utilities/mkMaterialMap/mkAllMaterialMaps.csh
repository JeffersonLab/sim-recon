#!/bin/tcsh -f

# Target
mkMaterialMap -Nr 2 -Nz 2 -rmin 0 -rmax 1.5 -zmin 50 -zmax 80 -n_r 5 -n_z 5 -n_phi 10
mv material_map material_map00_target

# Target wall
mkMaterialMap -Nr 20 -Nz 2 -rmin 1.50 -rmax 1.55 -zmin 50 -zmax 80 -n_r 100 -n_z 5 -n_phi 10
mv material_map material_map01_target_wall

# Scattering chamber
mkMaterialMap -Nr 100 -Nz 400 -rmin 0.0 -rmax 4.7 -zmin 43.9 -zmax 86.1 -n_r 100 -n_z 5 -n_phi 10
mv material_map material_map02_scattering_chamber

# Start Counter barrel
mkMaterialMap -Nr 100 -Nz 100 -rmin 4.7 -rmax 9.0 -zmin 17.0 -zmax 86.1 -n_r 30 -n_z 30 -n_phi 100
mv material_map material_map03_startcounter_barrel

# Start Counter nose
mkMaterialMap -Nr 100 -Nz 50 -rmin 2.0 -rmax 9.0 -zmin 86.1 -zmax 97.5  -n_r 30 -n_z 30 -n_phi 20
mv material_map material_map04_startcounter_nose

# CDC endplate
mkMaterialMap -Nr 60 -Nz 60 -rmin 9.0 -rmax 60.0 -zmin 167 -zmax 171 -n_r 10 -n_z 10 -n_phi 60
mv material_map material_map10_CDC_endplate

# CDC inner shell
mkMaterialMap -Nr 50 -Nz 5 -rmin 9.0 -rmax 9.75 -zmin 17 -zmax 167 -n_r 1000 -n_z 50 -n_phi 10
mv material_map material_map11_CDC_inner_shell

# CDC
mkMaterialMap -Nr 10 -Nz 10 -rmin 9.75 -rmax 56.0 -zmin 17 -zmax 167 -n_r 100 -n_z 100 -n_phi 200
mv material_map material_map12_CDC

# CDC outer shell
mkMaterialMap -Nr 50 -Nz 5 -rmin 56.0 -rmax 65.0 -zmin 17 -zmax 167 -n_r 1000 -n_z 2 -n_phi 10
mv material_map material_map13_CDC_outer_shell

# FDC package 1
mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 174 -zmax 189 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map21_FDC1

# FDC1-2 interpackage spacer
mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 189 -zmax 232.5 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map41_FDC_inter1

# FDC package 2
mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 232.5 -zmax 247.5 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map22_FDC2

# FDC2-3 interpackage spacer
mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 247.5 -zmax 291 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map42_FDC_inter2

# FDC package 3
mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 291 -zmax 306 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map23_FDC3

# FDC3-4 interpackage spacer
mkMaterialMap -Nr 10 -Nz 10 -rmin 50.0 -rmax 55.0 -zmin 306 -zmax 329.5 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map43_FDC_inter3

# FDC package 4
mkMaterialMap -Nr 30 -Nz 20 -rmin 0.0 -rmax 60.5 -zmin 329.5 -zmax 344.5 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map23_FDC4

# FDC cables
mkMaterialMap -Nr 30 -Nz 10 -rmin 61.0 -rmax 65.0 -zmin 167 -zmax 365 -n_r 1000 -n_z 1000 -n_phi 10
mv material_map material_map31_cables

# Everything else
#mkMaterialMap -Nr 120 -Nz 700 -rmin 0.0 -rmax 120.0 -zmin -50 -zmax 650 -n_r 3 -n_z 3 -n_phi 100
#mv material_map material_map99_course_default

# Make simple map with air everywhere as default
touch material_map
echo "#    r       z       A      Z   density  radlen   rhoZ_overA  rhoZ_overA_logI  chi2c_factor  chi2a_factor  chi2a_corr" >> material_map
echo "#%  00      01      02     03        04      05           06               07            08            09          10" >> material_map
echo "   0.5    -49.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126" >> material_map
echo " 119.5    -49.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126" >> material_map
echo "   0.5    649.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126" >> material_map
echo " 119.5    649.5  14.803  7.374  0.001214  30035  0.000604743      -0.00977523   0.000795067   7.60355e-05  0.00967126" >> material_map
mv material_map material_map99_course_default
