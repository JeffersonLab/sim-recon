#!/usr/bin/env python


#  makeTranslationTable.py

#   makes a fake translation table in XML format


# Conventions:
#    start with crate 1 slot 1
#    start new crate with each new detector
#    start new crate with each new module type
#    reserved slots:
#        1 is cpu
#       10 is empty
#       11 is switch slot
#       12 is switch slot
#       21 is TI
#    for Hall D max of 16 payload slots per VXS crate, 8 on each side of switch slots, slot 10 unused


#  still to do:
#    should include cpu, switch slots and TI


#  E. Wolin, JLab, 8-Jul-2011



import sys


# for testing:  1 is on, 0 is off
detectorOn = {
    'CDC':     1,
    'FDC':     0,
    'BCAL':    0,
    'FCAL':    0,
    'SC':      0,
    'TOF':     0,
    'TAGGER':  0
    }


# channel count by module type
channelCount = {
    'FADC250':  16,
    'FADC125':  72,
    'F1TDC32':  32,
    'F1TDC48':  48
    }    


# vme module number by count (-999 is dummy since Python arrays start at zero)
vmeModuleNumber = [-999,2,3,4,5,6,7,8,9,13,14,15,16,17,18,19,20]



#  cdc straw counts by ring
cdcStrawCount = [42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158,
                 158, 170, 170, 182, 182, 197, 197, 209, 209]




#  get output file name
if (len(sys.argv)>1):
    fileName   = sys.argv[1]
else:
    fileName   = "fakeTranslationTable.xml"
    


# open file and insert preamble
file=open(fileName,'w')
file.write('<!--   ' + fileName + ' -->\n')
file.write('\n')
file.write( '<!--  *** This file contains FAKE translation table data *** -->\n' )
file.write('\n')
file.write( '<!--            E. Wolin, JLab, 8-Jul-2011 -->' )
file.write('\n')
file.write('\n')
file.write('\n<translation_table version="0.1">')
file.write('\n')
file.write('\n')



# loop over detectors, some have both ADC and TDC channels
max_payload = 16
crate       = 1
vmeModule   = 0
file.write('  <crate number="%i"  type="VXS">\n\n' % crate)



# CDC:  FADC125, 72 channels/slot
if (detectorOn['CDC']==1):
    type = 'FADC125'
    vmeModule=1
    channel=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for ring in range(len(cdcStrawCount)):
        for straw in range(1,cdcStrawCount[ring]+1):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                vmeModule = vmeModule+1
                channel = 1
                if (vmeModule>max_payload):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    vmeModule = 1
                file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="CDC" ring="%i" straw="%i" />\n' % (channel,ring+1,straw) )
            channel = channel+1
    file.write('    </slot>\n\n')
                



# SC:  FADC250, 16 channels/slot
if (detectorOn['SC']==1):
    type = 'FADC250'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for sector in range(1,41):
        if (channel>channelCount[type]):
            file.write('    </slot>\n\n')
            vmeModule = vmeModule+1
            channel = 1
            if (vmeModule>max_payload):
                file.write('  </crate>\n\n\n')
                crate = crate+1
                file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                vmeModule = 1
            file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
        file.write('      <channel number="%i" detector="SC" sector="%i"  />\n' % (channel,sector) )
        channel = channel+1
    file.write('    </slot>\n\n')

            

# SC: F1TDC, 32 channels/slot
if (detectorOn['SC']==1):
    type = 'F1TDC32'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for sector in range(1,41):
        if (channel>channelCount[type]):
            file.write('    </slot>\n\n')
            vmeModule = vmeModule+1
            channel = 1
            if (vmeModule>max_payload):
                file.write('  </crate>\n\n\n')
                crate = crate+1
                file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                vmeModule = 1
            file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
        file.write('      <channel number="%i" detector="SC" sector="%i"  />\n' % (channel,sector) )
        channel = channel+1
    file.write('    </slot>\n\n')
            


# FCAL: FADC250, 16 channels/slot
if (detectorOn['FCAL']==1):
    type = 'FADC250'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for row in range(59):
        for column in range(59):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                vmeModule = vmeModule+1
                channel = 1
                if (vmeModule>max_payload):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    vmeModule = 1
                file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="FCAL" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            


# BCAL: FADC250, 16 channels/slot
if (detectorOn['BCAL']==1):
    type = 'FADC250'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for end in range(2):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,11):
                    if (channel>channelCount[type]):
                        file.write('    </slot>\n\n')
                        vmeModule = vmeModule+1
                        channel = 1
                        if (vmeModule>max_payload):
                            file.write('  </crate>\n\n\n')
                            crate = crate+1
                            file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                            vmeModule = 1
                        file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%i" />\n'
                               % (channel,module,sector,layer,end) )
                    channel = channel+1
    file.write('    </slot>\n\n')




# BCAL: F1TDC, 32 channels/slot
if (detectorOn['BCAL']==1):
    type = 'F1TDC32'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for end in range(2):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,11):
                    if (channel>channelCount[type]):
                        file.write('    </slot>\n\n')
                        vmeModule = vmeModule+1
                        channel = 1
                        if (vmeModule>max_payload):
                            file.write('  </crate>\n\n\n')
                            crate = crate+1
                            file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                            vmeModule = 1
                        file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%i" />\n'
                               % (channel,module,sector,layer,end) )
                    channel = channel+1
    file.write('    </slot>\n\n')
            





# FDC: FADC125, 72 channels/slot, 4 packages, 6 triplets(cathode,anode,cathode) per package, 216 cathode strips/layer
if (detectorOn['FDC']==1):
    type = 'FADC125'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for package in range(4):
        for triplet in range(6):
            for plane in [1,3]:
                for element in range(1,217):
                    if (channel>channelCount[type]):
                        file.write('    </slot>\n\n')
                        vmeModule = vmeModule+1
                        channel = 1
                        if (vmeModule>max_payload):
                            file.write('  </crate>\n\n\n')
                            crate = crate+1
                            file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                            vmeModule = 1
                        file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                    gplane = package*18 + triplet*3 + plane
                    file.write('      <channel number="%i" detector="FDCCathode" gPlane ="%i" element="%i" />\n'
                               % (channel,gplane,element) )
                    channel = channel+1
    file.write('    </slot>\n\n')




            
# FDC: F1TDC48, 48 channels/slot, 4 packages, 6 triplets(cathode,anode,cathode) per package, 96 anodes/layer
if (detectorOn['FDC']==1):
    type = 'F1TDC48'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for package in range(4):
        for triplet in range(6):
            for element in range(1,97):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    vmeModule = vmeModule+1
                    channel = 1
                    if (vmeModule>max_payload):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        vmeModule = 1
                    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                gplane = package*18 + triplet*3 + 2
                file.write('      <channel number="%i" detector="FDCAnode" gPlane="%i" element="%i" />\n'
                           % (channel,gplane,element) )
                channel = channel+1
    file.write('    </slot>\n\n')
            



# TOF: FADC250, 16 channels/slot
if (detectorOn['TOF']==1):
    type = 'FADC250'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for end in range(2):
        for plane in range(2):
            for bar in range(45):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    vmeModule = vmeModule+1
                    channel = 1
                    if (vmeModule>max_payload):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        vmeModule = 1
                    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%i" />\n'
                           % (channel,plane,bar,end) )
                channel = channel+1
    file.write('    </slot>\n\n')
            


# TOF: F1TDC, 32 channels/slot
if (detectorOn['TOF']==1):
    type = 'F1TDC32'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for end in range(2):
        for plane in range(2):
            for bar in range(45):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    vmeModule = vmeModule+1
                    channel = 1
                    if (vmeModule>max_payload):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        vmeModule = 1
                    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%i" />\n'
                           % (channel,plane,bar,end) )
                channel = channel+1
    file.write('    </slot>\n\n')
            


# TAGGER: FADC250, 16 channels/slot
if (detectorOn['TAGGER']==1):
    type = 'FADC250'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for row in range(9):
        for column in range(101):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                vmeModule = vmeModule+1
                channel = 1
                if (vmeModule>max_payload):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    vmeModule = 1
                file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            


# TAGGER:  F1TDC, 32 channels/slot
if (detectorOn['TAGGER']==1):
    type = 'F1TDC32'
    vmeModule = vmeModule+1
    channel = 1
    if (vmeModule>max_payload):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        vmeModule=1
    file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')

    for row in range(9):
        for column in range(101):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                vmeModule = vmeModule+1
                channel = 1
                if (vmeModule>max_payload):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    vmeModule = 1
                file.write(('    <slot number="%i"'  % vmeModuleNumber[vmeModule]) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            



#  done
file.write('  </crate>\n\n\n')
file.write('\n')
file.write('</translation_table>\n')
file.close()
sys.exit(0)
