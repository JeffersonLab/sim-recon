import os
import glob

Import('*')
env = env.Clone()

env.PrependUnique(CPPPATH = ['include'])

# Build static library from all source
myobjs = env.Object(Glob('src/*.c*'))
mylib = env.Library(target = "xstream", source = myobjs)

# Cleaning and installation are restricted to the directory
# scons was launched from or its descendents
CurrentDir = env.Dir('.').srcnode().abspath
if not CurrentDir.startswith(env.GetLaunchDir()):
	# Not in launch directory. Tell scons no to clean these targets
	env.NoClean([myobjs, mylib])
else:
	# We're in launch directory (or descendent) schedule installation

	# Installation directories for library and headers
	incdir = env.subst('$INCDIR')
	libdir = env.subst('$LIBDIR')

	# Install targets 
	env.Install(libdir, mylib)
	env.Install(incdir, Glob('include/*.h*'))
	env.Install("%s/xstream" %(incdir), Glob('include/xstream/*.h*'))

