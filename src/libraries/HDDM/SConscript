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
# source Nodes and utility Nodes are both available.
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
# files as the first target so we can easily strip off the ".h" (or ".hpp")
# suffix to get the base name of the target. This is what is passed into
# the "-o" option of hddm-c as it will automatically add ".c" and ".h" to this.
#
# Another complication is that by default, the C++ files produced will
# have a ".cpp" suffix, but the same base name. This would normally lead
# to object files for hddm_s.c and hddm_s.cpp both being hddm_s.o.
# We follow the historical solution used in BMS of renaming the hddm_s.cpp
# to hddm_s++.cpp.


import re
import subprocess
import SCons
import sbms

# get env object and clone it
Import('*')

#========================================================================
# Python functions used by the hddm-c and hddm-cpp builders

#---------------
# HDDM_C
#---------------
def HDDM_C(target, source, env):

	# Get full path to tool
	hddmc   = str(env.File("#.%s/programs/Utilities/hddm/hddm-c" % osname))

	# Get basename with full path for target.
	# The first target should always be the header
	# file name to be generated so we just drop the
	# suffix from that.
	target_base = re.sub('\.h$', '', str(target[0]))

	# Form command to be executed and execute it
	cmd = [hddmc, '-o', target_base, str(source[0])]
	if( int(env['SHOWBUILD']) > 0): print ' '.join(cmd)
	cmdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

#---------------
# HDDM_CPP
#---------------
def HDDM_CPP(target, source, env):

	# Get full path to tool
	hddmcpp = str(env.File("#.%s/programs/Utilities/hddm/hddm-cpp" % osname))

	# Get basename with full path for target.
	# The first target should always be the header
	# file name to be generated so we just drop the
	# suffix from that.
	target_base = re.sub('\.hpp$', '', str(target[0]))

	# Form command to be executed and execute it
	cmd = [hddmcpp, '-o', target_base, str(source[0])]
	if( int(env['SHOWBUILD']) > 0): print ' '.join(cmd)
	cmdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

	# C++ interface must have source file names modified
	# so as not to generate same .o filenames as their C counteparts
	cmd = ['mv', '%s.cpp' % target_base, '%s++.cpp' % target_base]
	if( int(env['SHOWBUILD']) > 0): print ' '.join(cmd)
	cmdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()


# end of Python functions. Global-level scons continues below
#========================================================================

# To get succinct output, define the builder actions based
# on the value of SHOWBUILD
if env['SHOWBUILD']==0:
	hddmcaction   = SCons.Script.Action(HDDM_C  , 'HDDM-C     [$SOURCE]')
	hddmcppaction = SCons.Script.Action(HDDM_CPP, 'HDDM-CPP   [$SOURCE]')
else:
	hddmcaction   = SCons.Script.Action(HDDM_C)
	hddmcppaction = SCons.Script.Action(HDDM_CPP)

# Create the actual scons builders from the actions defined above
bldc   = SCons.Script.Builder(action = hddmcaction)
bldcpp = SCons.Script.Builder(action = hddmcppaction)
env.Append(BUILDERS = {'HDDMC'   : bldc})
env.Append(BUILDERS = {'HDDMCPP' : bldcpp})


# Add the C/C++ HDDM (de)serializer source dependencies
# by hand. The output file names depend on the class
# tag defined inside the XML file. So, specifying the
# inputs and outputs explicitly here is just easier.
# The HDDMC_SRC and HDDMCPP_SRC variables are used by
# the SConscript in src/programs/Utilities/hddm.
env.AppendUnique(HDDMC_SRC   = env.HDDMC(['hddm_s.h', 'hddm_s.c'], 'event.xml'))
env.AppendUnique(HDDMC_SRC   = env.HDDMC(['hddm_mc_s.h', 'hddm_mc_s.c'], 'mc.xml'))
env.AppendUnique(HDDMC_SRC   = env.HDDMC(['hddm_r.h', 'hddm_r.c'], 'rest.xml'))
env.AppendUnique(HDDMCPP_SRC = env.HDDMCPP(['hddm_s.hpp', 'hddm_s++.cpp'], 'event.xml'))
env.AppendUnique(HDDMCPP_SRC = env.HDDMCPP(['hddm_mc_s.hpp', 'hddm_mc_s++.cpp'], 'mc.xml'))
env.AppendUnique(HDDMCPP_SRC = env.HDDMCPP(['hddm_r.hpp', 'hddm_r++.cpp'], 'rest.xml'))

# Finally, clone the build environment and make a library
# out of all source. This should include the generated
# HDDM (de)serializer routines.
env = env.Clone()

sbms.AddDANA(env)
sbms.library(env)

#========================================================================

# now we try to build wrapper libraries - these are only built if the swig
# executable exists and that building these executables is enabled
# these are needed for other systems that work with HDDM files, e.g. EventStore
sbms.AddSWIG(env)
swig_env = env.Clone()
swig_env.AppendUnique(SWIGFLAGS = ["-c++","-python"])
swig_env.AppendUnique(LIBS = ["z","bz2"])
sbms.swig_library(swig_env, "pyhddm_r", ["pyhddm_r.i"])

