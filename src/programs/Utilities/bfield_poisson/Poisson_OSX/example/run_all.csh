#!/bin/tcsh


# Run the various poisson/superfish programs to generate a magnetic
# field map compatible with the GlueX code base.
#
# Inputs:
#
# gluex_sol_example1.am    Geometery configuration file
#
# gluex_sol_example1.in7   Output grid configuration file
#
#
# NOTE: It is assumed that "wine" is in your PATH and that the
# LANL environment variable is set prior to running this script

if (! $?LANL) then
	echo " "
	echo "You must set your LANL environment variable to point"
	echo "to your Poisson/Superfish installation directory!"
	echo " "
	exit(-1);
endif


# Run AUTOMESH.EXE
echo
echo "}}}}}}}}}}}}}}  Running AUTOMESH.EXE"
wine $LANL/AUTOMESH.EXE ./gluex_sol_example1.am

# Run POISSON.EXE
echo
echo "}}}}}}}}}}}}}}  Running POISSON.EXE"
wine $LANL/POISSON.EXE ./GLUEX_SOL_EXAMPLE1.T35

# Run SF7.EXE
echo
echo "}}}}}}}}}}}}}}  Running SF7"
wine $LANL/SF7 ./gluex_sol_example1.in7

# Transform map into GlueX ASCII format
echo
echo "}}}}}}}}}}}}}}  Running poisson2calibDB.pl"
./poisson2calibDB.pl OUTSF7.TXT

echo
echo
echo "To view a graphical representation of the map using WSFPLOT.EXE"
echo "do the following:"
echo
echo "wine $LANL/WSFPLOT.EXE ./GLUEX_SOL_EXAMPLE1.T35"
echo
echo
echo
echo "To generate a ROOT file of the new map, execute something"
echo "like the following:"
echo
echo "cp solenoid_1500_poisson ~/HallD/calib/Magnets/Solenoid/solenoid_test"
echo "bfield2root -PBFIELD_MAP=Magnets/Solenoid/solenoid_test"
echo
