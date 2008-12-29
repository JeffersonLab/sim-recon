#!/bin/tcsh

rm -rf source
mkdir -p source
#cp -r $JANA_HOME/include source
#cd source/include/JANA
cp -r $HALLD_HOME/include source
cd source/include/TOF

rm *_Dict.*
rm *LinkDef.h
rm -r *JIL*

foreach f (`ls *.h | grep -v jerror.h | grep -v JFactory.h`)
	echo Generating dictionary for $f ...
	rootcint `basename $f .h`_Dict.C -c -I.. -I$JANA_HOME/include $f
end
#cat JEventLoop_Dict.C | sed -e s/vector\<error_call_stack_t/vector\<JEventLoop::error_XXX_t/g > tmp.C
#mv tmp.C JEventLoop_Dict.C
#cat JEventLoop_Dict.C | sed -e s/vector\<call_stack_t/vector\<JEventLoop::call_stack_t/g > tmp.C
#mv tmp.C JEventLoop_Dict.C
#cat JEventLoop_Dict.C | sed -e s/vector\<JEventLoop::error_XXX_t/vector\<JEventLoop::error_call_stack_t/g > tmp.C
#mv tmp.C JEventLoop_Dict.C

echo " "
echo " "
echo "=========== Compiling exectuable ========="
g++ `root-config --cflags --libs` -lHtml -lThread -o gendoc *_Dict.C ../../../gendoc.cc -L$JANA_HOME/lib/$OSNAME -lJANA
./gendoc

rm -rf ../../../htmldoc
mv htmldoc ../../../
