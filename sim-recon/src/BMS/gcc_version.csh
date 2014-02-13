#!/bin/tcsh
#
# gcc_version.csh
# May 2, 2008 David Lawrence
#

# Return the gcc version number of the current system.
# This will extract the gcc version by running gcc -v
# and parsing the contents of its output. Only the 
# version number is printed followed by a newline.
#
# By default, the entire version number is printed
# which typically includes 3 numbers, seperated by
# periods. If, however, the first argument of this
# script is the string "majoronly" (no quotes) then
# only the first character of the version number is
# printed which should be the major revision number.

if ($1 =~ 'majoronly') then
	gcc -v | & awk '/gcc version/ {printf("%c\n",$3)}'
else
	gcc -v | & awk '/gcc version/ {print $3}'
endif

