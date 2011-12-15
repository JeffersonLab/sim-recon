#!/bin/tcsh -f

set savedir=$PWD

rm -rf source
mkdir -p source
cp -r $JANA_HOME/include source
cd source/include/JANA

rm -f *_Dict.*
rm -f *LinkDef.h
rm -rf *JIL*

# Tweak files a bit in preparing for rootcint
foreach f (`ls *.h`)
	# Comment out namespace jana lines since rootcint doesn't handle it well
	#cat $f | sed -e 's/namespace jana/\/\/namespace jana/'  | sed -e 's/JObject:://g' | sed -e 's/jana:://g' > tmp.h
	#cat $f | sed -e 's/namespace jana/namespace std/'  | sed -e 's/JObject:://g' | sed -e 's/jana:://g' > tmp.h
	cat $f | grep -v "namespace jana" | grep -v "JANA namespace" | sed -e 's/JObject::oid_t/oid_t/g' | sed -e 's/jana:://g' > tmp.h
	mv tmp.h $f
end


foreach f (`ls J*.h | grep -v JFactory.h`)
#foreach f (JGeometryXML.h)
	echo Generating dictionary for $f ...
	#echo rootcint `basename $f .h`_Dict.C -c -I.. -I$JANA_HOME/include $f
	rootcint `basename $f .h`_Dict.C -c -I.. -I$JANA_HOME/include -I$XERCESCROOT/include $f
end


echo " "
echo " "
echo "=========== Compiling exectuable ========="
g++ `root-config --cflags --libs` -lHtml -lThread -I.. -o gendoc *_Dict.C ../../../gendoc.cc -L$JANA_HOME/lib/$BMS_OSNAME -lJANA
./gendoc

rm -rf ../../../htmldoc
mv htmldoc ../../../

cd $savedir
