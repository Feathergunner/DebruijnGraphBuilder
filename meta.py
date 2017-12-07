#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import matplotlib.pyplot as plt
import math

import data_io as dio

'''
def print_progress(part, total):
	print ("Progress: "+str("%.2f" % ((float(part)/(float(total)/100)))) + "%")
'''
	
def get_adaptive_k(readlength):
	'''
	if readlength < 100:
		return 25
	elif readlength < 200:
		return 30
	elif readlength < 500:
		return 35
	elif readlength < 1000:
		return 40
	elif readlength < 3000:
		return 45
	else:
		return 50
	'''
	return int(math.log(readlength, 2)*4)

def print_progress(part, total, front_string="Progress:", end_string=""):
	if not total == 0:
		print front_string+" "+str("%6.2f" % ((float(part)/(float(total)/100)))) + "% "+end_string+"\r",
		if part >= total:
			print 
		sys.stdout.flush()

def get_inverse_sequence(sequence, alphabet={"A":"T", "C":"G", "G":"C", "T":"A"}):
	n = len(sequence)
	inv_sequence = [""]*n
	for char_position in range(len(sequence)):
		current_char = sequence[char_position]
		if current_char in alphabet:
			inv_sequence[n-char_position-1] = alphabet[current_char]
		else:
			print (sequence)
			print (current_char)
			print ("Error! Incorrect Alphabet!")
			break
	return ''.join(inv_sequence)

def compute_insert_distance(sequence_1, sequence_2, maxdist = -1):
	# algorithm may not work properly for arbitrary large insert-distances,
	# but is correct if local insert-distace is <= 1
	# returns insert-distance if insert_distance <= maxdist,
	# otherwise returns maxdist + x for a x >= 1
	if not len(sequence_1) == len(sequence_2):
		return -1
	index_1 = 0
	index_2 = 0
	
	insert_distance = 0
	n = len(sequence_1)
	if maxdist < 0:
		maxdist = n
	while (insert_distance < (maxdist+1) and index_1 < n and index_2 < n):
		t1 = 0
		t2 = 0
		
		while (index_1 + t1 < n and (not sequence_1[index_1+t1] == sequence_2[index_2])):
			t1 += 1
		while (index_2 + t2 < n and (not sequence_1[index_1] == sequence_2[index_2+t2])):
			t2 += 1
			
		if t1 <= t2:
			index_1 += t1
			insert_distance += t1
		elif t2 < t1:
			index_2 += t2
			insert_distance += t2
			
		index_1 += 1
		index_2 += 1
	return insert_distance
	
def get_readlength_distribution_from_fastq_file(filename, bucketsize=1000):
	readlengths = sorted([x[0] for x in dio.get_readlengths(filename)])
	minlength = min(readlengths)
	maxlength = max(readlengths)
	
	x = []
	y = []
	current_index = 0
	for b in range(minlength, maxlength, bucketsize):
		x.append(b)
		y_val = 0
		while (current_index < len(readlengths) and readlengths[current_index] < b+bucketsize):
			y_val += 1
			current_index += 1
		#print (str(b)+" : "+str(y_val))
		y.append(y_val)
			
	#plt.plot(x, [math.log10(i+1) for i in y])
	plt.plot(x, y)
	plt.show()
	return x, y

def get_readlength_distribution(reads, bucketsize=1000):
	n = len(reads)
	readlengths = [0]*n
	for i in range(n):
		readlengths[i] = [len(reads[i]),i]
	sorted_readlengths = sorted([x[0] for x in readlengths])
	minlength = min(sorted_readlengths)
	maxlength = max(sorted_readlengths)
	avglength = sum(sorted_readlengths)/n
	
	x = []
	y = []
	current_index = 0
	for b in range(minlength, maxlength, bucketsize):
		x.append(b)
		y_val = 0
		while (current_index < n and sorted_readlengths[current_index] < b+bucketsize):
			y_val += 1
			current_index += 1
		#print (str(b)+" : "+str(y_val))
		y.append(y_val)
	
	print ("minlength: "+str(minlength))
	print ("maxlength: "+str(maxlength))
	print ("avglength: "+str(avglength))

	plt.plot(x, [math.log10(i+1) for i in y])
	plt.show()
	plt.plot(x, y)
	plt.show()
	return x, y