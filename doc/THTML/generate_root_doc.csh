#!/bin/tcsh

rm -rf source
mkdir -p source
cp -r $JANA_HOME/include source
cd source/include/JANA

rm *_Dict.*
rm *LinkDef.h
rm -r *JIL*

foreach f (*.h)
	echo Generating dictionary for $f ...
	rootcint `basename $f .h`_Dict.C -c $f
end

echo " "
echo " "
echo "=========== Compiling exectuable ========="
g++ `root-config --cflags --libs` -lHtml -lThread -o gendoc *_Dict.C ../../../gendoc.cc -L$JANA_HOME/lib/$OSNAME -lJANA
./gendoc

rm -rf ../../../htmldoc
mv htmldoc ../../../
