import os
import glob

Import('env osname installdir')

env.PrependUnique(CPPPATH = ['include'])

# Build static library from all source
mylib = env.Library(target = "xstream", source = Glob('src/*.c*'))

# Installation directories for library and headers
incdir = env.subst('$INCDIR')
libdir = env.subst('$LIBDIR')

# Install targets 
env.Install(libdir, mylib)
env.Install(incdir, Glob('include/*.h*'))
env.Install("%s/xstream" %(incdir), Glob('include/xstream/*.h*'))
env.Alias('install', installdir)
