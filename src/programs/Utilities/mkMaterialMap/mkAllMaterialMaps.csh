#!/bin/tcsh -f

# Target
mkMaterialMap -Nr 2 -Nz 2 -rmin 0 -rmax 1.5 -zmin 50 -zmax 80 -n_r 5 -n_z 5 -n_phi 10
mv material_map material_map00_target

# Target wall
mkMaterialMap -Nr 20 -Nz 2 -rmin 1.50 -rmax 1.55 -zmin 50 -zmax 80 -n_r 100 -n_z 5 -n_phi 10
mv material_map material_map01_target_wall

# Scattering chamber
mkMaterialMap -Nr 100 -Nz 400 -rmin 0.0 -rmax 4.0 -zmin 35 -zmax 90 -n_r 100 -n_z 5 -n_phi 10
mv material_map material_map02_scattering_chamber

# Start Counter barrel
mkMaterialMap -Nr 30 -Nz 2 -rmin 7.0 -rmax 9.0 -zmin 0 -zmax 90 -n_r 3 -n_z 3 -n_phi 100
mv material_map material_map03_startcounter_barrel

# Start Counter nose
mkMaterialMap>mkMaterialMap -Nr 100 -Nz 50 -rmin 1.0 -rmax 9.0 -zmin 90 -zmax 100 -n_r 30 -n_z 30 -n_phi 20
mv material_map material_map04_startcounter_nose

# CDC endplate
mkMaterialMap -Nr 60 -Nz 60 -rmin 9.0 -rmax 60.0 -zmin 165 -zmax 171 -n_r 3 -n_z 3 -n_phi 60
mv material_map material_map10_CDC_endplate

# CDC inner shell
mkMaterialMap -Nr 50 -Nz 5 -rmin 9.0 -rmax 9.75 -zmin 17 -zmax 165 -n_r 1000 -n_z 50 -n_phi 10
mv material_map material_map11_CDC_inner_shell

# CDC
mkMaterialMap -Nr 50 -Nz 5 -rmin 9.75 -rmax 56.0 -zmin 17 -zmax 165 -n_r 100 -n_z 50 -n_phi 200
mv material_map material_map12_CDC

# CDC outer shell
mkMaterialMap -Nr 50 -Nz 5 -rmin 56.0 -rmax 65.0 -zmin 17 -zmax 165 -n_r 1000 -n_z 2 -n_phi 10
mv material_map material_map13_CDC_outer_shell

# FDC package 1
mkMaterialMap -Nr 6 -Nz 2 -rmin 2.0 -rmax 60.0 -zmin 175 -zmax 189 -n_r 50 -n_z 20000 -n_phi 10
mv material_map material_map21_FDC1

# FDC package 2
mkMaterialMap -Nr 6 -Nz 2 -rmin 2.0 -rmax 60.0 -zmin 232 -zmax 246 -n_r 50 -n_z 20000 -n_phi 10
mv material_map material_map22_FDC2

# FDC package 3
mkMaterialMap -Nr 6 -Nz 2 -rmin 2.0 -rmax 60.0 -zmin 290 -zmax 304 -n_r 50 -n_z 20000 -n_phi 10
mv material_map material_map23_FDC3

# FDC package 4
mkMaterialMap -Nr 6 -Nz 2 -rmin 2.0 -rmax 60.0 -zmin 347 -zmax 361 -n_r 50 -n_z 20000 -n_phi 10
mv material_map material_map23_FDC4

# Everything else
#mkMaterialMap -Nr 120 -Nz 700 -rmin 0.0 -rmax 120.0 -zmin -50 -zmax 650 -n_r 3 -n_z 3 -n_phi 100
#mv material_map material_map99_course_default

# Make simple map with air everywhere as default
touch material_map
echo "#  r	z	A	Z	density radlen rhoZ_overA rhoZ_overA_logI" >> material_map
echo "0.5	-49.5	14.803	7.374	0.001214	30035	0.000604743	-0.00977523" >> material_map
echo "119.5	-49.5	14.803	7.374	0.001214	30035	0.000604743	-0.00977523" >> material_map
echo "0.5	649.5	14.803	7.374	0.001214	30035	0.000604743	-0.00977523" >> material_map
echo "119.5	649.5	14.803	7.374	0.001214	30035	0.000604743	-0.00977523" >> material_map
mv material_map material_map99_course_default
