#!/bin/tcsh -f

# Get full path to current script
set savedir=$PWD
set scriptdir=`dirname $0`
set scriptdir=`cd $scriptdir && pwd`
cd $savedir

# Try and set BMS_OSNAME if it's not already set
if ( ! $?BMS_OSNAME ) then

	# Look for osrelease.pl relative to script dir
	set osrelease='none'
	if ( -x $scriptdir/src/BMS/osrelease.pl ) then
		set osrelease=$scriptdir/src/BMS/osrelease.pl
		
	# Look relative to current directory
	else if ( -x src/BMS/osrelease.pl ) then
		set osrelease=src/BMS/osrelease.pl
	else if ( -x BMS/osrelease.pl ) then
		set osrelease=BMS/osrelease.pl
	else if ( -x ./osrelease.pl ) then
		set osrelease=./osrelease.pl
	endif
	
	# If osrelease.pl script is not found, then
	# we can't do anything
	if ( $osrelease == 'none' ) then
		echo "BMS_OSNAME not set and no osrelease.pl script found!"
	else
		set BMS_OSNAME=`$osrelease`
	endif

endif

# BMS_OSNAME may be set either upon entry to script
# or by block above. If neither set it, then we already
# printed an error message. If it is set, the use it 
# to source the platform-specific setenv.csh
if ( $?BMS_OSNAME ) then

	if ( -e $scriptdir/$BMS_OSNAME/setenv.csh ) then
		source $scriptdir/$BMS_OSNAME/setenv.csh
	else
		echo "Cannot find $scriptdir/$BMS_OSNAME/setenv.csh"
	endif

endif

