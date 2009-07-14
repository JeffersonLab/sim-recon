#!/bin/tcsh -f

# Target, start-counter
mkMaterialMap -Nr 500 -Nz 1500 -rmin 0 -rmax 9.75 -zmin 15 -zmax 100 -n_r 3 -n_z 3 -n_phi 60
mv material_map material_map_sec1

# CDC
mkMaterialMap -Nr 50 -Nz 150 -rmin 9.75 -rmax 65.0 -zmin 15 -zmax 165 -n_r 5 -n_z 5 -n_phi 100
mv material_map material_map_sec2

# CDC endplate and 1st FDC package
mkMaterialMap -Nr 200 -Nz 1000 -rmin 0.0 -rmax 65.0 -zmin 165 -zmax 190 -n_r 3 -n_z 3 -n_phi 60
mv material_map material_map_sec3

# Downstream of start-counter near beamline (just air)
mkMaterialMap -Nr 2 -Nz 2 -rmin 0.0 -rmax 9.75 -zmin 100 -zmax 165 -n_r 1 -n_z 1 -n_phi 360
mv material_map material_map_sec4

# FDC, 2nd to 4th packages
mkMaterialMap -Nr 65 -Nz 2000 -rmin 0.0 -rmax 65 -zmin 190 -zmax 380 -n_r 3 -n_z 3 -n_phi 60
mv material_map material_map_sec5
