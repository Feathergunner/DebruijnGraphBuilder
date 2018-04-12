#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
from random import randint
import re
from numpy import random as nura
import math

import meta

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
	
def get_random_mutation(nucleotide, alphabet=['A', 'C', 'G', 'T']):
	possible_mutations = [c for c in alphabet if not c == nucleotide]
	return random.choice(possible_mutations)

def generate_dna(length=100, alphabet=["A","C","G","T"]):
	print ("Generating dna sequence ...")
	dna = []
	for i in range(length):
		dna.append(random.choice(alphabet))
	return dna
	
def generate_set_of_related_dnas(length=10000, mean_partsize=100, variation_partsize=0, number_of_dnas=2, alphabet=["A","C","G","T"], verbose=False):
	dnas = []
	dnas.append(generate_dna(length=length, alphabet=alphabet))
	for i in range(1, number_of_dnas):
		meta.print_progress(i, number_of_dnas-1)
		# choose an already created dna as base:
		base_dna = random.choice(dnas)
		new_dna = ['a']*length
		current_pos = 0
		while current_pos < length:
			#meta.print_progress(current_pos, length-1, front_string="\tModify dna: ")
			# get random part size by gaussian distribution
			next_part_size = min(max(10, int(random.gauss(mu=mean_partsize, sigma=variation_partsize))), length-current_pos)
			if verbose:
				print ("Size of next part: "+str(next_part_size))
			# randomly determine whether part is original or variation:
			if random.random() > 0.5:
				# case original:
				if verbose:
					print ("use original sequence")
				for i in range(current_pos, current_pos+next_part_size):
					new_dna[i] = base_dna[i]
			else:
				# case variation:
				if verbose:
					print ("use variation")
				for i in range(current_pos, current_pos+next_part_size):
					new_dna[i] = random.choice(alphabet)
			current_pos += next_part_size
		dnas.append(new_dna)
	return dnas		
	
def write_dna_to_file(filename, dna):
	with open(filename, 'w') as outf:
		outf.write(''.join(dna) + '\n')

def genereate_reads(dna, coverage=50, avg_read_length=50, mutation_pct=0.1, mutation_alphabet=["A","C","G","T"], remove_pct=5, variation=0, readlength_distribution='gaussian', paired_end=False, both_directions=False, verbose=False):
	print ("Generating reads ...")
	#print ("mutation_pct = "+str(mutation_pct) + ", remove_pct = "+str(remove_pct))
	dna_length = len(dna)
	reads = []
	alignment = []
	for i in range(coverage):
		if verbose:
			print ("iteration: " + str(i))
		cur_index = 0
		while cur_index < dna_length:
			#next_index = min(cur_index+avg_read_length+random.choice(range(-10,10)), dna_length)
			if variation < 0:
				variation = 1+avg_read_length/5
			if readlength_distribution == 'gaussian':
				read_length_variation = int(random.gauss(mu=0, sigma=variation))
			elif readlength_distribution == 'exponential':
				read_length_variation = int(nura.exponential(scale=variation))
			if verbose:
				print "read_length_variation: "+str(read_length_variation)
			next_index = min(cur_index+avg_read_length+read_length_variation, dna_length)
			new_read = dna[cur_index:next_index]
			for s_i in range(len(new_read)):
				if random.choice(range(10000)) < mutation_pct*100:
					if verbose:
						print ("mutation at read "+str(s_i))
					new_read[s_i] = random.choice(mutation_alphabet)
			if random.choice(range(100)) > remove_pct:
				if both_directions and random.choice(range(100)) > 50:
					reads.append([get_inverse_sequence(new_read)])
				else:
					reads.append([''.join(new_read)])
				alignment.append([cur_index, next_index-1])
			cur_index = next_index

	### TO DO: generate paired-end information

	return reads, alignment
	
def samplereads(dna, 
				number_of_reads,
				replace_error_percentage=0.0,
				indel_error_percentage=0.0,
				mutation_alphabet=["A","C","G","T"],
				read_length_mean=50,
				read_length_stddev=0,
				readlength_distribution='gaussian',
				inverted_reads=False,
				uniform_coveragedepth=False,
				verbose=False):
				
	print ("Sample reads from dna ...")
	
	print "read_length_mean: "+str(read_length_mean)
	reads = []
	# construct reads:
	for n in range(number_of_reads):
		meta.print_progress(n, number_of_reads-1)
		
		if readlength_distribution == 'gaussian':
			readlength = int(random.gauss(mu=read_length_mean, sigma=read_length_stddev))
		elif readlength_distribution == 'exponential':
			readlength = int(0.1*read_length_mean + 0.9*nura.exponential(scale=read_length_mean))
			if readlength > read_length_mean:
				scalefactor = 1+(1.0*readlength)/len(dna)
				readlength = int(readlength*scalefactor)
		else:
			readlength = read_length_mean
		
		# ensure that reads are not longer than dna:
		readlength = min(readlength, len(dna))
		# choose a random start position of read:
		if uniform_coveragedepth:
			read_start_index_tmp = randint(-readlength, len(dna))
			read_start_index = max(0, read_start_index_tmp)
			readlength = max(len(dna)-read_start_index, readlength)
			
		read_start_index = randint(0, len(dna)-readlength)
		
		sampleread_base = dna[read_start_index:read_start_index+readlength]
		sampleread_errors = sampleread_base
		
		j = 0
		for i in range(len(sampleread_base)):
			# replacement error:
			if randint(0,10000) < int(100*replace_error_percentage):
				sampleread_errors[j] = get_random_mutation(sampleread_base[i], mutation_alphabet)
			# indel error:
			if randint(0,10000) < int(100*indel_error_percentage):
				# in/del 50/50:
				if randint(0,100) < 50:
					# insert:
					insertbase = random.choice(mutation_alphabet)
					sampleread_errors = sampleread_errors[:j] + [insertbase] + sampleread_errors[j:]
					j += 1
				else:
					sampleread_errors = sampleread_errors[:j] + sampleread_errors[(j+1):]
					j -= 1
			j += 1

		sampleread = ''.join(sampleread_errors)
		
		if (inverted_reads and randint(0,100) < 50):
			reads.append(get_inverse_sequence(sampleread.upper()))
		else:
			reads.append(sampleread.upper())
			
	return reads
	
if __name__ == "__main__":
	# test:
	dna = generate_dna()
	reads = samplereads(dna, 1000, 1.0, 1.0)