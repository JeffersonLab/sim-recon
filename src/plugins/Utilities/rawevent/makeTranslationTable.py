#!/usr/bin/env python

#  makeTranslationTable.py

#   makes a fake translation table in XML format


# Conventions:
#    start with crate 1 slot 1
#    start new crate with each new detector
#    start new crate with each new module type
#    reserved slots:
#        1 is vme cpu
#       10 is empty
#       11 is switch slot
#       12 is switch slot
#       21 is TID
#    for Hall D max of 16 payload slots per VXS crate, 8 on each side of switch slots, slot 10 unused
#    channel number starts with 0 in table


#  still to do:
#    get rid of hard-coded channel counts
#    finish tagger microscope
#    complete detector list:  PS, TAC, POLAR, others?


#  E. Wolin, JLab, 26-Jul-2012



import sys
import datetime


# for testing:  1 is on, 0 is off
detectorOn = {
    'CDC':     1,
    'FDC':     1,
    'BCAL':    1,
    'FCAL':    1,
    'ST':      1,
    'TOF':     1,
    'TAGGER':  1,
    'PSPEC':   0,
    'TAC':     0,
    'POLAR':   0
    }


# channel count by module type
channelCount = {
    'FADC250':  16,
    'FADC125':  72,
    'F1TDC32':  32,
    'F1TDC48':  48,
    'CAENTDC':  32
    }    


# max 16 daq boards per crate, 8 per side
# array holds vme slot number given daq (ADC,TDC) module number
# missing slots hold CPU, CTP, SD and TID
max_payload   = 16
vmeSlotNumber = [-9999,2,3,4,5,6,7,8,9,13,14,15,16,17,18,19,20]



#  cdc straw counts by ring
cdcStrawCount = [42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158,
                 158, 170, 170, 182, 182, 197, 197, 209, 209]


#  misc functions


def startCrate(acrate,adaqModule,achannel):
    acrate     = acrate+1
    adaqModule = 1
    achannel   = 0
    file.write('  <crate number="%i"  type="VXS">\n\n' % acrate)
    file.write('    <slot number="1"  type="VMECPU"/>\n\n')
    return(acrate,adaqModule,achannel)

def startSlot(adaqModule,atype):
    file.write('    <slot number="%i"'  % vmeSlotNumber[adaqModule] + '  type="' + atype + '">\n')

def incrChannel(achannel,adaqModule,acrate,atype):
    achannel=achannel+1
    if (achannel>channelCount[atype]):
        endSlot()
        adaqModule = adaqModule+1
        if (adaqModule>max_payload):
            endCrate()
            (acrate,adaqModule,achannel) = startCrate(acrate,adaqModule,achannel)
        elif (vmeSlotNumber[adaqModule]==13):
            addSwitchSlots()
        startSlot(adaqModule,atype)
        achannel = 1
    return(achannel,adaqModule,acrate)


def addSwitchSlots():
    file.write('    <slot number="11"  type="CTP"/>\n\n')
    file.write('    <slot number="12"  type="SD"/>\n\n')

def endSlot():
    file.write('    </slot>\n\n')
    
def endCrate():
    file.write('    <slot number="21"  type="TID"/>\n\n')
    file.write('  </crate>\n\n\n')
    
def closeCrate(adaqModule):
    endSlot()
    if(vmeSlotNumber[adaqModule]<13):
        addSwitchSlots()
    endCrate()




#  get output file name
if (len(sys.argv)>1):
    fileName   = sys.argv[1]
else:
    fileName   = "fakeTranslationTable.xml"
    


# open file and insert preamble
file=open(fileName,'w')
file.write('<!--   ' + fileName + ' -->\n')
file.write('\n')
file.write('<!-- (need XML boilerplate here) -->\n')
file.write('\n')
file.write( '<!--  *** This file contains FAKE translation table data *** -->\n' )
file.write('\n')
file.write('<!-- E. Wolin, JLab, ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + ' -->')
file.write('\n')
file.write('\n')
file.write('\n')
file.write('\n<halld_online_translation_table version="0.1">')
file.write('\n')
file.write('\n')
file.write('\n')
file.write('\n')



# loop over detectors, some have both ADC and TDC channels
crate     = 0
daqModule = 0
channel   = 0


# CDC:  FADC125, 72 channels/slot
if (detectorOn['CDC']>0):
    type = 'FADC125'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for ring in range(len(cdcStrawCount)):
        for straw in range(1,cdcStrawCount[ring]+1):
            (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
            file.write('      <channel number="%i" detector="CDC" ring="%i" straw="%i" />\n' % (channel-1,ring+1,straw) )
    closeCrate(daqModule)



# ST:  FADC250, 16 channels/slot
if (detectorOn['ST']>0):
    type = 'FADC250'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for sector in range(1,41):
        (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
        file.write('      <channel number="%i" detector="ST" sector="%i"  />\n' % (channel-1,sector) )
    closeCrate(daqModule)

            

# ST: F1TDC, 32 channels/slot
if (detectorOn['ST']>0):
    type = 'F1TDC32'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for sector in range(1,41):
        (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
        file.write('      <channel number="%i" detector="ST" sector="%i"  />\n' % (channel-1,sector) )
    closeCrate(daqModule)
            


# FCAL: FADC250, 16 channels/slot
if (detectorOn['FCAL']>0):
    type = 'FADC250'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for row in range(59):
        for column in range(59):
            (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
            file.write('      <channel number="%i" detector="FCAL" row="%i" column="%i" />\n' % (channel-1,row,column) )
    closeCrate(daqModule)
            


# BCAL: FADC250, 16 channels/slot
if (detectorOn['BCAL']>0):
    type = 'FADC250'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for end in range(2):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,11):
                    (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%i" />\n'
                               % (channel-1,module,sector,layer,end) )
    closeCrate(daqModule)




# BCAL: F1TDC, 32 channels/slot
if (detectorOn['BCAL']>0):
    type = 'F1TDC32'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for end in range(2):
        for module in range(1,49):
            for sector in range(1,5):
                for layer in range(1,11):
                    (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                    file.write('      <channel number="%i" detector="BCAL" module="%i" sector="%i" layer="%i" end="%i" />\n'
                               % (channel-1,module,sector,layer,end) )
    closeCrate(daqModule)





# FDC: FADC125, 72 channels/slot, 4 packages, 6 triplets(cathode,anode,cathode) per package, 216 cathode strips/layer
if (detectorOn['FDC']>0):
    type = 'FADC125'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for package in range(4):
        for triplet in range(6):
            for plane in [1,3]:
                for element in range(1,217):
                    (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                    gplane = package*18 + triplet*3 + plane
                    file.write('      <channel number="%i" detector="FDCCathode" gPlane ="%i" element="%i" />\n'
                               % (channel-1,gplane,element) )
    closeCrate(daqModule)




            
# FDC: F1TDC48, 48 channels/slot, 4 packages, 6 triplets(cathode,anode,cathode) per package, 96 anodes/layer
if (detectorOn['FDC']>0):
    type = 'F1TDC48'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for package in range(4):
        for triplet in range(6):
            for element in range(1,97):
                (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                gplane = package*18 + triplet*3 + 2
                file.write('      <channel number="%i" detector="FDCAnode" gPlane="%i" element="%i" />\n'
                           % (channel-1,gplane,element) )
    closeCrate(daqModule)
            



# TOF: FADC250, 16 channels/slot
if (detectorOn['TOF']>0):
    type = 'FADC250'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for end in range(2):
        for plane in range(2):
            for bar in range(45):
                (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%i" />\n'
                           % (channel-1,plane,bar,end) )
    closeCrate(daqModule)
            


#### TOF: CAENTDC, 32 channels/slot
# TOF: F1TDCC, 32 channels/slot
if (detectorOn['TOF']>0):
    type = 'F1TDC32'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for end in range(2):
        for plane in range(2):
            for bar in range(45):
                (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
                file.write('      <channel number="%i" detector="TOF" plane="%i" bar="%i" end="%i" />\n'
                           % (channel-1,plane,bar,end) )
    closeCrate(daqModule)
            


# TAGGER: FADC250, 16 channels/slot
if (detectorOn['TAGGER']>0):
    type = 'FADC250'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
###    for row in range(9):
    for row in range(1):
        for column in range(128):
            (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel-1,row,column) )
    closeCrate(daqModule)
            


# TAGGER:  F1TDC, 32 channels/slot
if (detectorOn['TAGGER']>0):
    type = 'F1TDC32'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for row in range(1):
###    for row in range(1):
        for column in range(128):
            (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
            file.write('      <channel number="%i" detector="TAGGER" row="%i" column="%i" />\n' % (channel-1,row,column) )
    closeCrate(daqModule)
            


# PSPEC: ???
if (detectorOn['PSPEC']>0):
    type = 'F1TDC48'
    (crate,daqModule,channel) = startCrate(crate,daqModule,channel)
    startSlot(daqModule,type)
    for counter in range(128):
        (channel,daqModule,crate) = incrChannel(channel,daqModule,crate,type)
        file.write('      <channel number="%i" detector="PSPEC" counter="%i" />\n' % (channel-1,counter) )
    closeCrate(daqModule)
            



#  done
file.write('\n')
file.write('</halld_online_translation_table>\n')
file.close()
sys.exit(0)
