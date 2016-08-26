
import os
import sys
import subprocess
import glob

# Add SBMS directory to PYTHONPATH
sbmsdir = "%s/SBMS" % (os.getcwd())
sys.path.append(sbmsdir)

import sbms

# Get command-line options
SHOWBUILD = ARGUMENTS.get('SHOWBUILD', 0)
OPTIMIZATION = ARGUMENTS.get('OPTIMIZATION', 2)
DEBUG = ARGUMENTS.get('DEBUG', 1)
PROFILE = ARGUMENTS.get('PROFILE', 0)
BUILDSWIG = ARGUMENTS.get('BUILDSWIG', 0)

# Get platform-specific name
osname = os.getenv('BMS_OSNAME', 'build')

# Get architecture name
arch = ROOT_CFLAGS = subprocess.Popen(["uname"], stdout=subprocess.PIPE).communicate()[0].strip()

# Setup initial environment
installdir = "#../%s" %(osname)
include = "%s/include" % (installdir)
bin = "%s/bin" % (installdir)
lib = "%s/lib" % (installdir)
plugins = "%s/plugins" % (installdir)
python2 = "%s/python2" % (installdir)
python3 = "%s/python3" % (installdir)
env = Environment(        ENV = os.environ,  # Bring in full environement, including PATH
                      CPPPATH = [include],
                      LIBPATH = [lib],
                  variant_dir = ".%s" % (osname))

# These are SBMS-specific variables (i.e. not default scons ones)
env.Replace(INSTALLDIR    = installdir,
				OSNAME        = osname,
				INCDIR        = include,
				BINDIR        = bin,
				LIBDIR        = lib,
				PLUGINSDIR    = plugins,
				PYTHON2DIR    = python2,
				PYTHON3DIR    = python3,
				ALL_SOURCES   = [],        # used so we can add generated sources
				MISC_OBJECTS  = [],        # used so we can add custom built objects
				SHOWBUILD     = SHOWBUILD,
				OPTIMIZATION  = OPTIMIZATION,
				DEBUG         = DEBUG,
				BUILDSWIG     = BUILDSWIG,
		  		COMMAND_LINE_TARGETS = COMMAND_LINE_TARGETS)

# Use terse output unless otherwise specified
if SHOWBUILD==0:
	env.Replace(CCCOMSTR       = "Compiling  [$SOURCE]",
				  CXXCOMSTR       = "Compiling  [$SOURCE]",
				  FORTRANPPCOMSTR = "Compiling  [$SOURCE]",
				  FORTRANCOMSTR   = "Compiling  [$SOURCE]",
				  SHCCCOMSTR      = "Compiling  [$SOURCE]",
				  SHCXXCOMSTR     = "Compiling  [$SOURCE]",
				  LINKCOMSTR      = "Linking    [$TARGET]",
				  SHLINKCOMSTR    = "Linking    [$TARGET]",
				  INSTALLSTR      = "Installing [$TARGET]",
				  ARCOMSTR        = "Archiving  [$TARGET]",
				  RANLIBCOMSTR    = "Ranlib     [$TARGET]")


# Get compiler from environment variables (if set)
env.Replace( CXX = os.getenv('CXX', 'c++'),
             CC  = os.getenv('CC' , 'cc'),
             FC  = os.getenv('FC' , 'gfortran') )

# Get compiler name
compiler = 'unknown'
compiler_string = subprocess.Popen([env['CC'],"-v"], stderr=subprocess.PIPE).communicate()[1]
if 'clang' in compiler_string:
	compiler = 'clang'
if 'gcc' in compiler_string and 'clang' not in compiler_string:
	compiler = 'gcc'
env.Replace(COMPILER = compiler)

# Add libraries and libraries/include to include search path
env.PrependUnique(CPPPATH = ['#', '#libraries', '#libraries/include'])

# Use C++11
env.PrependUnique(    CXXFLAGS = ['-std=c++11'])

# Standard flags (optimization level and warnings)
env.PrependUnique(      CFLAGS = ['-O%s' % OPTIMIZATION, '-fPIC', '-Wall'])
env.PrependUnique(    CXXFLAGS = ['-O%s' % OPTIMIZATION, '-fPIC', '-Wall'])
env.PrependUnique(FORTRANFLAGS = ['-O%s' % OPTIMIZATION, '-fPIC', '-Wall'])

# Turn on debug symbols unless user told us not to
if not DEBUG=='0':
	env.PrependUnique(      CFLAGS = ['-g'])
	env.PrependUnique(    CXXFLAGS = ['-g'])
	env.PrependUnique(FORTRANFLAGS = ['-g'])

# Turn on profiling if user asked for it
if PROFILE=='1':
	env.PrependUnique(      CFLAGS = ['-pg'])
	env.PrependUnique(    CXXFLAGS = ['-pg'])
	env.PrependUnique(FORTRANFLAGS = ['-pg'])
	env.PrependUnique(   LINKFLAGS = ['-pg'])

# Apply any platform/architecture specific settings
sbms.ApplyPlatformSpecificSettings(env, arch)
sbms.ApplyPlatformSpecificSettings(env, osname)

# "external" packages (aka xstream)
SConscript('external/SConscript', variant_dir=".%s/external" % (osname), exports='env osname', duplicate=0)

# build libraries
SConscript('libraries/SConscript', variant_dir=".%s/libraries" % (osname), exports='env osname', duplicate=0)

# build programs
program_subdirs = ['programs']
SConscript(dirs=program_subdirs, variant_dir=".%s/programs" % (osname), exports='env osname', duplicate=0)

# build plugins
program_subdirs = ['plugins']
SConscript(dirs=program_subdirs, variant_dir=".%s/plugins" % (osname), exports='env osname', duplicate=0)

# Make install target
env.Alias('install', installdir)

# Create setenv if user explicitly specified "install" target
build_targets = map(str,BUILD_TARGETS)
if len(build_targets)>0:
	if 'install' in build_targets:
		import sbms_setenv
		sbms_setenv.mk_setenv_csh(env)
		sbms_setenv.mk_setenv_bash(env)

