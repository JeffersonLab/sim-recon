
import os
import subprocess
import SCons

##################################
# library
##################################
def library(env, installdir):

	# Library name comes from directory name
	libname = os.path.split(os.getcwd())[1]

	env.PrependUnique(CPPPATH = ['.'])

	# Build static library from all source
	mylib = env.Library(target = libname, source = env.Glob('*.c*'))

	# Installation directories for library and headers
	includedir = "%s/%s" %(env.subst('$INCDIR'), libname)
	libdir = env.subst('$LIBDIR')

	# Install targets 
	env.Install(libdir, mylib)
	env.Install(includedir, env.Glob('*.h*'))
	env.Alias('install', installdir)


##################################
# executable
##################################
def executable(env, installdir):

	# Executable name comes from directory name
	exename = os.path.split(os.getcwd())[1]

	env.PrependUnique(CPPPATH = ['.'])

	# Build program from all source
	myexe = env.Program(target = exename, source = env.Glob('*.c*'))

	# Installation directories for library and headers
	includedir = env.subst('$INCDIR')
	bindir = env.subst('$BINDIR')

	# Install targets 
	env.Install(bindir, myexe)
	env.Install(includedir, env.Glob('*.h*'))
	env.Alias('install', installdir)


##################################
# JANA
##################################
def AddJANA(env):
	JANA_CFLAGS = subprocess.Popen(["jana-config", "--cflags"], stdout=subprocess.PIPE).communicate()[0]
	JANA_LINKFLAGS = subprocess.Popen(["jana-config", "--libs"], stdout=subprocess.PIPE).communicate()[0]
	env.PrependUnique(CCFLAGS = JANA_CFLAGS.split())
	env.PrependUnique(LINKFLAGS = JANA_LINKFLAGS.split())


##################################
# HDDS
##################################
def AddHDDS(env):
	hdds_home = os.getenv('HDDS_HOME', 'hdds')
	env.PrependUnique(CPPPATH = ["%s/src" % hdds_home])


##################################
# DANA
##################################
def AddDANA(env):
	AddJANA(env)
	AddCCDB(env)
	AddHDDS(env)
	AddXERCES(env)
	Add_xstream(env)
	DANA_LIBS  = "DANA ANALYSIS PID TAGGER TRACKING START_COUNTER"
	DANA_LIBS += " CERE RICH CDC TRIGGER"
	DANA_LIBS += " FDC TOF BCAL FCAL CCAL HDGEOMETRY HDDM JANA"
	env.PrependUnique(LIBS = DANA_LIBS.split())

##################################
# xstream
##################################
def Add_xstream(env):
	env.PrependUnique(CPPPATH = ['#external/xstream/include'])
	env.AppendUnique(LIBS=['xstream', 'bz2'])


##################################
# CCDB
##################################
def AddCCDB(env):
	ccdb_home = os.getenv('CCDB_HOME', 'ccdb')
	CCDB_CPPPATH = "%s/include" % (ccdb_home)
	CCDB_LINKFLAGS = "-L%s/lib -lccdb" % (ccdb_home)
	env.PrependUnique(CPPPATH = CCDB_CPPPATH.split())
	env.PrependUnique(LINKFLAGS = CCDB_LINKFLAGS.split())


##################################
# Xerces
##################################
def AddXERCES(env):
	xercescroot = os.getenv('XERCESCROOT', 'xerces')
	XERCES_CPPPATH = "%s/include" % (xercescroot)
	XERCES_LINKFLAGS = "-L%s/lib -lxerces-c" % (xercescroot)
	env.PrependUnique(CPPPATH = XERCES_CPPPATH.split())
	env.PrependUnique(LINKFLAGS = XERCES_LINKFLAGS.split())


##################################
# ROOT
##################################
def AddROOT(env):
	ROOT_CFLAGS = subprocess.Popen(["root-config", "--cflags"], stdout=subprocess.PIPE).communicate()[0]
	ROOT_LINKFLAGS = subprocess.Popen(["root-config", "--glibs"], stdout=subprocess.PIPE).communicate()[0]
	env.PrependUnique(CCFLAGS = ROOT_CFLAGS.split())
	env.PrependUnique(LINKFLAGS = ROOT_LINKFLAGS.split())
	env.PrependUnique(LINKFLAGS = "-lGeom")

