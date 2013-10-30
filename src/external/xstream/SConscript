import os
import glob

Import('env osname installdir')

env.PrependUnique(CPPPATH = ['include'])

# Build static library from all source
myobjs = env.Object(Glob('src/*.c*'))
mylib = env.Library(target = "xstream", source = myobjs)

# Installation directories for library and headers
incdir = env.subst('$INCDIR')
libdir = env.subst('$LIBDIR')

# Install targets 
installed = env.Install(libdir, mylib)
env.Install(incdir, Glob('include/*.h*'))
env.Install("%s/xstream" %(incdir), Glob('include/xstream/*.h*'))
env.Alias('install', installdir)

# Only clean these sources when scons -c is invoked in
# this directory or in a direct ancestor
CurrentDir = env.Dir('.').srcnode().abspath
if not CurrentDir.startswith(env.GetLaunchDir()):
	env.NoClean([myobjs, mylib, installed])

