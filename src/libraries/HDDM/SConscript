#
# Nov. 20, 2013  David Lawrence
#
# This is pretty ugly (they way we handle this in scons, not HDDM itself).
# The way HDDM is set up in the repository is the following:
#
# utilities used for generating files like hddm_s.c and hddm_s.h
# from the input event.xml have their source in:
# src/programs/Utilities/hddm
#
# The HDDM library is made up of code maintained in the
# src/libraries/HDDM directory as well as code generated
# from XML files in that same directory using the above
# utilities.
#
# What this means is that we need to tell scons that a file like "hddm_s.c"
# depends on the utility "hddm-c" being made first. Because those live
# in very separate directories we have to somehow get the Node objects for
# both the utility and the source in one place so we can explicitly
# declare this dependency using the scons "Requires()" routine. The way
# this is done is by adding the list of dependent nodes to a special
# variable in the top-level build environment. (In this case two lists
# because the C and C++ interfaces must be handled independently.)
# These variables, "HDDMC_SRC" and "HDDMCPP_SRC" are lists of nodes that
# can be used in the src/programs/Utilities/hddm/SConscript file to
# declare this dependency. It is done there since the top-level SConstruct
# file sources this file (the one you're reading) first (via the
# libraries/SConscript file) and the one there second (via the
# programs/SConscript file). Thus, it is the only place where both the
# source Nodes and utilitiy Nodes are both available.
#
# To actually generate the C/C++ source from the XML we create two builders
# (one for C and the other for C++) here. The builders run the utilities
# hddm-c and hddm-cpp from the variant dir tree. This path, relative to
# the top-level directory (the one containing SConstruct) is hardwired
# here (e.g. "#.%s/programs/Utilities/hddm/hddm-c"). To further complicate
# matters, the output of hddm-c is two files (e.g. hddm_s.c and hddm_s.h)
# which do not depend on the name of the input file ("event.xml"). We specify
# two products of the input XML, but the value that $TARGET gets set to
# in the builder is just the first of these. We choose to always use the header
# files as the first target so the external shell "basename" command can
# be used and hardwired to strip off the ".h" (or ".hpp") suffix to get the
# base name of the target. This is what is passed into the "-o" option
# of hddm-c as it will automatically add ".c" and ".h" to this.
#
# Another complication is that by default, the C++ files produced will
# have a ".cpp" suffix, but the same base name. This would normally lead
# to object files for hddm_s.c and hddm_s.cpp both being hddm_s.o.
# We follow the historical solution to this by renaming the hddm_s.cpp
# to hddm_s++.cpp.
#
# It should be possible to use a python function instead of an
# SCons.Script.Action possibly making the code here a little easier to
# understand and almost certainly making the output when running
# with SHOWBUILD=1 more readable. This seems to be working for the
# moment and is able to print a simplified command line (i.e. when
# SHOWBUILD is not set) that is consistent with other commands.


import SCons
import sbms

# get env object and clone it
Import('*')

# Locations of the built (not "installed") utilities for
# generating hddm C/C++ source from XML input files.
hddmc   = str(env.File("#.%s/programs/Utilities/hddm/hddm-c" % osname))
hddmcpp = str(env.File("#.%s/programs/Utilities/hddm/hddm-cpp" % osname))

# These are used multiple times in the builder below. They
# have to be scripts that are run when the builder issues
# the command so that the file names can be derived from
# $TARGET which isn't set until then.
basenamec   = "`dirname $TARGET`/`basename $TARGET .h`"
basenamecpp = "`dirname $TARGET`/`basename $TARGET .hpp`"

# Create builders to generated C/C++ HDDM (de)serializer
# code using the XML input.
if env['SHOWBUILD']==0:
	hddmcaction   = SCons.Script.Action("%s -o %s $SOURCE" % (hddmc, basenamec), 'HDDM-C     [$SOURCE]')
	hddmcppaction = SCons.Script.Action("%s -o %s $SOURCE ; mv %s.cpp %s++.cpp" % (hddmcpp, basenamecpp, basenamecpp, basenamecpp), 'HDDM-CPP   [$SOURCE]')
else:
	hddmcaction   = SCons.Script.Action("%s -o %s $SOURCE" % (hddmc, basenamec))
	hddmcppaction = SCons.Script.Action("%s -o %s $SOURCE ; mv %s.cpp %s++.cpp" % (hddmcpp, basenamecpp, basenamecpp, basenamecpp))
bldc   = SCons.Script.Builder(action = hddmcaction)
bldcpp = SCons.Script.Builder(action = hddmcppaction)
env.Append(BUILDERS = {'HDDMC'   : bldc})
env.Append(BUILDERS = {'HDDMCPP' : bldcpp})

# Add the C/C++ HDDM (de)serializer source dependcies
# by hand. The output file names depend on the class
# tag defined inside the XML file. So, specifying the
# inputs and outputs explicitly here is just easier.
env.AppendUnique(HDDMC_SRC   = env.HDDMC(['hddm_s.h', 'hddm_s.c'], 'event.xml'))
env.AppendUnique(HDDMC_SRC   = env.HDDMC(['hddm_r.h', 'hddm_r.c'], 'rest.xml'))
env.AppendUnique(HDDMCPP_SRC = env.HDDMCPP(['hddm_s.hpp', 'hddm_s++.cpp'], 'event.xml'))
env.AppendUnique(HDDMCPP_SRC = env.HDDMCPP(['hddm_r.hpp', 'hddm_r++.cpp'], 'rest.xml'))

# Finally, clone the build environment and make a library
# out of all source. This should include the generated
# HDDM (de)serializer routines.
env = env.Clone()

sbms.AddDANA(env)
sbms.library(env)

