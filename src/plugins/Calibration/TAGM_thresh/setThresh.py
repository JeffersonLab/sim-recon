#! /usr/bin/python

# This script will generate a new .cnf file containing
# new individual channel thresholds based off of the 
# original roctagm1_test_v1.cnf file. Run this script
# with a single argument that will act as a multiple

#infile = open('roctagm1_test_v1.cnf','r')
infile = open('roctagm1_spring_2017_v3.cnf','r')
threshfile = open('tagm-thresh.txt','r')
outfile = open('thresh.cnf','w')
thresh = [] #thresholds based on pulse integral

for line in threshfile:
	col = int(line.split()[1])
	thr = int(float(line.split()[2])/4.3)
	#print str(col) + ' ' + str(thr)
	thresh.append(thr)

threshTT = list(thresh) #thresholds adjusted for translation table swaps
for i in range(6): # fix columns 1-6
	new_index = 3*(i/3)+(2-i%3)
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i+1) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(6,9): # fix the set around the ind. channels
	new_index = 3*(i/3)+(2-i%3)+5
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i+1) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(9,14): # get the ind. channels
	threshTT[i-3] = thresh[i]
	print 'ind col: 9 adc: ' + str(i-3+1) + ' thr: ' + str(thresh[i])
for i in range(14,29): # fix columns 10-24
	new_index = 3*((i+1)/3)+(2-(i+1)%3)-1
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i+1-5) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(29,32): # fix the set around the 2nd set of ind. channels
	new_index = 3*((i+1)/3)+(2-(i+1)%3)+5-1
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i+1-5) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(32,37): # fix second set ind channels
	threshTT[i-3] = thresh[i]
	print 'ind col: 27 adc: ' + str(i-3+1) + ' thr: ' + str(thresh[i])
for i in range(37,88):
	new_index = 3*((i-1)/3)+(2-(i-1)%3)+1
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i-9) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(88,91): # fix the set around the 3rd set of ind. channels
	new_index = 3*((i-1)/3)+(2-(i-1)%3)+1+5
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i-9) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(91,96): # fix third set ind channels
	threshTT[i-3] = thresh[i]
	print 'ind col: 81 adc: ' + str(i-3+1) + ' thr: ' + str(thresh[i])
for i in range(96,111):
	new_index = 3*((i)/3)+(2-(i)%3)
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i-14) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(111,114): # fix the set around the 4th set of ind. channels
	new_index = 3*((i)/3)+(2-(i)%3)+5
	threshTT[new_index] = thresh[i]
	print 'col: ' + str(i-14) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])
for i in range(114,119): # fix third set ind channels
	threshTT[i-3] = thresh[i]
	print 'ind col: 99 adc: ' + str(i-3+1) + ' thr: ' + str(thresh[i])
for i in range(119,122):
	new_index = i
	print 'col: ' + str(i-19) + ' adc: ' + str(new_index+1) + ' thr: ' + str(thresh[i])

counter = -1
print str(len(thresh))
for line in infile:
	if (line.startswith('FADC250_ALLCH_THR')):
		counter += 1
		print counter
		outfile.write('FADC250_ALLCH_THR   ')
		for i in range(16):
			# Set unused thresholds to 1000 adc
			if (counter == 7 and i > 9 and i < 15):
				outfile.write(str(1000) + '   ')
			elif (counter == 7 and i == 15):
				outfile.write(str(1000) + '\n')
			elif (i == 15):
				outfile.write(str(100+threshTT[counter*16+i]) + '\n')
			else:
				outfile.write(str(100+threshTT[counter*16+i]) + '   ')
	elif (line.startswith('FADC250_ALLCH_DAC')):
		continue
	else:
		outfile.write(line)
