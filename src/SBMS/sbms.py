
import os
import subprocess
import SCons
import glob

#===========================================================
# The first 3 sections provide routines for building a
# library, program, or plugin from all the source in the
# current directory. The routines that follow add support
# for various packages
#===========================================================


##################################
# library
##################################
def library(env):

	# Library name comes from directory name
	libname = os.path.split(os.getcwd())[1]

	env.PrependUnique(CPPPATH = ['.'])

	# Add C/C++ targets
	env.AppendUnique(ALL_SOURCES = env.Glob('*.c*'))

	sources = env['ALL_SOURCES']

	# Build static library from all source
	myobjs = env.Object(sources)
	mylib = env.Library(target = libname, source = myobjs)

	# Cleaning and installation are restricted to the directory
	# scons was launched from or its descendents
	CurrentDir = env.Dir('.').srcnode().abspath
	if not CurrentDir.startswith(env.GetLaunchDir()):
		# Not in launch directory. Tell scons no to clean these targets
		env.NoClean([myobjs, mylib])
	else:
		# We're in launch directory (or descendent) schedule installation

		# Installation directories for library and headers
		installdir = env.subst('$INSTALLDIR')
		includedir = "%s/%s" %(env.subst('$INCDIR'), libname)
		libdir = env.subst('$LIBDIR')

		# Install targets 
		env.Install(libdir, mylib)
		env.Install(includedir, env.Glob('*.h*'))


##################################
# executable
##################################
def executable(env):

	# Executable name comes from directory name
	exename = os.path.split(os.getcwd())[1]

	env.PrependUnique(CPPPATH = ['.'])

	# Add C/C++ targets
	env.AppendUnique(ALL_SOURCES = env.Glob('*.c*'))

	sources = env['ALL_SOURCES']

	# Build program from all source
	myobjs = env.Object(sources)
	myexe = env.Program(target = exename, source = myobjs)

	# Cleaning and installation are restricted to the directory
	# scons was launched from or its descendents
	CurrentDir = env.Dir('.').srcnode().abspath
	if not CurrentDir.startswith(env.GetLaunchDir()):
		# Not in launch directory. Tell scons no to clean these targets
		env.NoClean([myobjs, myexe])
	else:
		# We're in launch directory (or descendent) schedule installation

		# Installation directories for executable and headers
		installdir = env.subst('$INSTALLDIR')
		includedir = env.subst('$INCDIR')
		bindir = env.subst('$BINDIR')

		# Install targets 
		env.Install(bindir, myexe)


##################################
# plugin
##################################
def plugin(env):

	# Library name comes from directory name
	pluginname = os.path.split(os.getcwd())[1]

	env.PrependUnique(CPPPATH = ['.'])

	# Add C/C++ targets
	env.AppendUnique(ALL_SOURCES = env.Glob('*.c*'))

	sources = env['ALL_SOURCES']

	# Build static library from all source
	myobjs = env.SharedObject(sources)
	myplugin = env.SharedLibrary(target = pluginname, source = myobjs, SHLIBPREFIX='', SHLIBSUFFIX='.so')

	# Cleaning and installation are restricted to the directory
	# scons was launched from or its descendents
	CurrentDir = env.Dir('.').srcnode().abspath
	if not CurrentDir.startswith(env.GetLaunchDir()):
		# Not in launch directory. Tell scons no to clean these targets
		env.NoClean([myobjs, myplugin])
	else:
		# We're in launch directory (or descendent) schedule installation

		# Installation directories for plugin and headers
		installdir = env.subst('$INSTALLDIR')
		includedir = "%s/%s" %(env.subst('$INCDIR'), pluginname)
		pluginsdir = env.subst('$PLUGINSDIR')

		# Install targets 
		installed = env.Install(pluginsdir, myplugin)
		env.Install(includedir, env.Glob('*.h*'))




#===========================================================
# Package support follows
#===========================================================


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
	env.PrependUnique(CCFLAGS = ['-fPIC'])
	env.AppendUnique(LIBS=['xstream', 'bz2', 'z'])


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

	# Create Builder that can convert .h file into _Dict.cc file
	rootsys = os.getenv('ROOTSYS', '/usr/local/root/PRO')
	env.AppendENVPath('LD_LIBRARY_PATH', '%s/lib' % rootsys )
	if env['SHOWBUILD']==0:
		rootcintaction = SCons.Script.Action("%s/bin/rootcint -f $TARGET -c $SOURCE" % (rootsys), 'ROOTCINT   [$SOURCE]')
	else:
		rootcintaction = SCons.Script.Action("%s/bin/rootcint -f $TARGET -c $SOURCE" % (rootsys))
	bld = SCons.Script.Builder(action = rootcintaction, suffix='_Dict.cc', src_suffix='.h')
	env.Append(BUILDERS = {'ROOTDict' : bld})
	env.Append(LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH'])

	# Generate ROOT dictionary file targets for each header
	# containing "ClassDef"
	#
	# n.b. It seems if scons is run when the build directory doesn't exist,
	# then the cwd is set to the source directory. Otherwise, it is the
	# build directory. Since the headers will only exist in the source
	# directory, we must temporarily cd into that to look for headers that
	# we wish to generate dictionaries for. (This took a long time to figure
	# out!)
	curpath = os.getcwd()
	srcpath = env.Dir('.').srcnode().abspath
	if(env['SHOWBUILD']!=0):
		print "---- Scanning for headers to generate ROOT dictionaries in: %s" % srcpath
	os.chdir(srcpath)
	for f in glob.glob('*.h*'):
		if 'ClassDef' in open(f).read():
			env.AppendUnique(ALL_SOURCES = env.ROOTDict(f))
			if(env['SHOWBUILD']!=0):
				print "       ROOT dictionary for %s" % f
	os.chdir(curpath)


