#!/usr/bin/env python


#  TextFile2String.py

#  reads file textFileName.*, creates header file textFileName.h that defines static string textFileName_String

#  E. Wolin, JLab, 1-Apr-2010



import sys


#  enough arguments?
if (len(sys.argv)<2):
    print '\n?file2string...Missing argsument'
    print '\nUsage:\n\n    python file2string fileName\n\n'
    sys.exit(1)



#  get names
textFileName   = sys.argv[1]
baseName       = sys.argv[1].split('.')[0]
headerFileName = baseName + '.h'
stringName     = baseName + '_String'


# open header file and insert preamble
header=open(headerFileName,'w')
header.write('//   ' + headerFileName + '\n')
header.write('\n')
header.write('\n')
header.write( '// *** DO NOT EDIT ***\n' )
header.write('\n')
header.write( '//  This file is automaticaly generated via a python script\n' )
header.write( '//      E. Wolin, JLab, 1-Apr-2010' )
header.write('\n')
header.write('\n')
header.write('\n')
header.write( '#ifndef _' + baseName + '_\n')
header.write( '#define _' + baseName + '_\n')
header.write('\n')
header.write('\n')
header.write('\n')
header.write( 'static string ' + stringName + ' =\n\n')


#  open text file, read/modify/write lines
textFile = open(textFileName,'r')
for line in textFile:
    header.write('"' + (line.replace('"','\\"')).rstrip('\n') + '\\n"\n')
textFile.close()


#  write final header file stuff and close
header.write( ';\n' )
header.write('\n')
header.write( '#endif // _' + baseName + '_\n' )
header.write('\n')

header.close()


#  done
sys.exit(0)
