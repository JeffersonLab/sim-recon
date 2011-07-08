#!/usr/bin/env python


#  makeTranslationTable.py

#   makes a fake translation table in XML format

#  E. Wolin, JLab, 8-Jul-2011



import sys


# constants
channelCount = {
    'FADC250':  16,
    'FADC125':  72,
    'F1TDC32':  32,
    'F1TDC64':  64
    }    


# 1 is on, 0 is off
detectorOn = {
    'CDC':     1,
    'FDC':     1,
    'BCAL':    1,
    'FCAL':    1,
    'SC':      1,
    'TOF':     1,
    'TAGGER':  1
    }


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


# for line in textFile:
#     header.write('"' + (line.replace('"','\\"')).rstrip('\n') + '\\n"\n')
# textFile.close()



# loop over detectors, some have both ADC and TDC channels
# start with crate 1 slot 1
# start new slot with each new detector
max_slot  = 16
crate     = 1
slot      = 1
channel   = 1
file.write('  <crate number="%i"  type="VXS">\n\n' % crate)



# CDC:  FADC125, 72 channels/slot
if (detectorOn['CDC']==1):
    type = 'FADC125'
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for ring in range(1,11):
        for straw in range(1,21):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                slot = slot+1
                channel = 1
                if (slot>max_slot):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    slot = 1
                file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="CDC" ring="%i" straw="%i" />\n' % (channel,ring,straw) )
            channel = channel+1
    file.write('    </slot>\n\n')
                



# SC:  FADC250, 16 channels/slot
if (detectorOn['SC']==1):
    type = 'FADC250'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for sector in range(1,41):
        if (channel>channelCount[type]):
            file.write('    </slot>\n\n')
            slot = slot+1
            channel = 1
            if (slot>max_slot):
                file.write('  </crate>\n\n\n')
                crate = crate+1
                file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                slot = 1
            file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
        file.write('      <channel number="%i" detector="SC" sector="%i"  />\n' % (channel,sector) )
        channel = channel+1
    file.write('    </slot>\n\n')

            

# SC: F1TDC, 32 channels/slot
if (detectorOn['SC']==1):
    type = 'F1TDC32'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for sector in range(1,41):
        if (channel>channelCount[type]):
            file.write('    </slot>\n\n')
            slot = slot+1
            channel = 1
            if (slot>max_slot):
                file.write('  </crate>\n\n\n')
                crate = crate+1
                file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                slot = 1
            file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
        file.write('      <channel number="%i" detector="SC" sector="%i"  />\n' % (channel,sector) )
        channel = channel+1
    file.write('    </slot>\n\n')
            


# FCAL: FADC250, 16 channels/slot
if (detectorOn['FCAL']==1):
    type = 'FADC250'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for row in range(1,21):
        for column in range(1,21):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                slot = slot+1
                channel = 1
                if (slot>max_slot):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    slot = 1
                file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="FCAL" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            


# BCAL: FADC250, 16 channels/slot
if (detectorOn['BCAL']==1):
    type = 'FADC250'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for end in ("upstream","downstream"):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,7):
                    if (channel>channelCount[type]):
                        file.write('    </slot>\n\n')
                        slot = slot+1
                        channel = 1
                        if (slot>max_slot):
                            file.write('  </crate>\n\n\n')
                            crate = crate+1
                            file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                            slot = 1
                        file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%s" />\n'
                               % (channel,module,sector,layer,end) )
                    channel = channel+1
    file.write('    </slot>\n\n')




# BCAL: F1TDC, 32 channels/slot
if (detectorOn['BCAL']==1):
    type = 'F1TDC32'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for end in ("upstream","downstream"):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,7):
                    if (channel>channelCount[type]):
                        file.write('    </slot>\n\n')
                        slot = slot+1
                        channel = 1
                        if (slot>max_slot):
                            file.write('  </crate>\n\n\n')
                            crate = crate+1
                            file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                            slot = 1
                        file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%s" />\n'
                               % (channel,module,sector,layer,end) )
                    channel = channel+1
    file.write('    </slot>\n\n')
            


# FDC: FADC125, 72 channels/slot
if (detectorOn['FDC']==1):
    type = 'FADC125'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for module in range(1,5):
        for layer in range(1,16):
            for element in range(1,31):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    slot = slot+1
                    channel = 1
                    if (slot>max_slot):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        slot = 1
                    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="FDCCathode" module="%i" layer="%i" element="%i" />\n'
                           % (channel,module,layer,element) )
                channel = channel+1
    file.write('    </slot>\n\n')
            

# FDC: F1TDC64, 64 channels/slot
if (detectorOn['FDC']==1):
    type = 'F1TDC64'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for module in range(1,5):
        for layer in range(1,16):
            for element in range(1,31):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    slot = slot+1
                    channel = 1
                    if (slot>max_slot):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        slot = 1
                    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="FDCAnode" module="%i" layer="%i" element="%i" />\n'
                           % (channel,module,layer,element) )
                channel = channel+1
    file.write('    </slot>\n\n')
            



# TOF: FADC250, 16 channels/slot
if (detectorOn['TOF']==1):
    type = 'FADC250'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for end in ("top","bottom","left","right"):
        for plane in range(1,2):
            for bar in range(1,21):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    slot = slot+1
                    channel = 1
                    if (slot>max_slot):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        slot = 1
                    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%s" />\n'
                           % (channel,plane,bar,end) )
                channel = channel+1
    file.write('    </slot>\n\n')
            


# TOF: F1TDC, 32 channels/slot
if (detectorOn['TOF']==1):
    type = 'F1TDC32'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for end in ("top","bottom","left","right"):
        for plane in range(1,2):
            for bar in range(1,21):
                if (channel>channelCount[type]):
                    file.write('    </slot>\n\n')
                    slot = slot+1
                    channel = 1
                    if (slot>max_slot):
                        file.write('  </crate>\n\n\n')
                        crate = crate+1
                        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                        slot = 1
                    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%s" />\n'
                           % (channel,plane,bar,end) )
                channel = channel+1
    file.write('    </slot>\n\n')
            


# TAGGER: FADC250, 16 channels/slot
if (detectorOn['TAGGER']==1):
    type = 'FADC250'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for row in range(1,50):
        for column in range(1,8):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                slot = slot+1
                channel = 1
                if (slot>max_slot):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    slot = 1
                file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            


# TAGGER:  F1TDC, 32 channels/slot
if (detectorOn['TAGGER']==1):
    type = 'F1TDC32'
    slot = slot+1
    channel = 1
    if (slot>max_slot):
        file.write('  </crate>\n\n\n')
        crate = crate+1
        file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
        slot=1
    file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')

    for row in range(1,50):
        for column in range(1,8):
            if (channel>channelCount[type]):
                file.write('    </slot>\n\n')
                slot = slot+1
                channel = 1
                if (slot>max_slot):
                    file.write('  </crate>\n\n\n')
                    crate = crate+1
                    file.write('  <crate number="%i"  type="VXS">\n\n' % crate)
                    slot = 1
                file.write(('    <slot number="%i"'  % slot) + '  type="' + type + '">\n')
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel,row,column) )
            channel = channel+1
    file.write('    </slot>\n\n')
            



#  done
file.write('  </crate>\n\n\n')
file.write('\n')
file.write('</translation_table>\n')
file.close()
sys.exit(0)
