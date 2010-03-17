@ECHO OFF
rem launching the ded
ECHO Launching ded for 12.0GeV

rem print out the java version to make sure we have version 1.6
echo requires java 1.6 checking...
java -version
echo ===========
echo

rem jevio jar 
set jeviocp=.\bin\jevio\jevio.jar

rem image jar
set imagecp=.\bin\bCNU\bcnuimages.jar

rem bCNU jar
set bcnucp=.\bin\bCNU\bCNU.jar

rem ded12GeV jar
set dedcp=.\bin\ded\ded.jar

rem combine to one classpath
set cp=%jeviocp%;%imagecp%;%bcnucp%;%dedcp%
echo classpath: %cp%
echo ===========
echo

ECHO ON
java -Xmx256M -Xss512k -cp %cp% cnuphys.ded.frame.Ded