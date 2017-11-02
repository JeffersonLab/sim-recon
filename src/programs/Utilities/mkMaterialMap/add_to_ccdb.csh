#!/bin/tcsh -f

# This script can be used to add all of the material maps
# to the CCDB indicated by the CCDB_CONNECTION envrionment
# variable.
#
# Use these variables to set the run range.
set RUN_MIN='30000'
set RUN_MAX='inf'


foreach f ( material_map??_* )
	echo "Writing Material/$f to $CCDB_CONNECTION"
	ccdb add Material/$f -r ${RUN_MIN}-${RUN_MAX} $f 
end

